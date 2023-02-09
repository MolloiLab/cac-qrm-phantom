### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 345ed6c8-06a1-4ece-8867-6b5a433fe83c
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-qrm-phantom"
	using PlutoUI, Statistics, CSV, DataFrames, CairoMakie, Colors, GLM, MLJBase
	using StatsBase: quantile!, rmsd
end

# ╔═╡ 600bab60-d643-49a6-94c3-77e0eccc727f
include(srcdir("plot_utils.jl"));

# ╔═╡ 7f6c5f2c-58ee-4c53-ba60-02ed0a75cd6f
TableOfContents()

# ╔═╡ 04013b46-69b8-476f-b93c-187622e7883a
root_path = datadir("output")

# ╔═╡ 2f9accf4-0b6d-4f32-9b64-df499fb17d15
medphys_theme = Theme(
    Axis = (
        backgroundcolor = :white,
		xgridcolor = :gray,
		xgridwidth = 0.1,
		xlabelsize = 20,
		xticklabelsize = 20,
		ygridcolor = :gray,
		ygridwidth = 0.1,
		ylabelsize = 20,
		yticklabelsize = 20,
		bottomsplinecolor = :black,
		leftspinecolor = :black,
		titlesize = 30
	)
);

# ╔═╡ 7e2f297c-c129-4c99-8ae3-63ea6271f3b0
md"""
## Accuracy (Motion)
"""

# ╔═╡ 41b322e8-f0b1-404d-94de-6910b0f481a0
md"""
### Integrated
"""

# ╔═╡ f65471f2-8f81-46f0-85fe-7aa5b07e1bf6
df_i1 = CSV.read(joinpath(root_path, "integrated", "motion.csv"), DataFrame); # 14.35

# ╔═╡ da76517e-6032-4617-bd86-9f500200fe9d
df_i1

# ╔═╡ 2b2adef3-c59c-4a2e-b758-2d2fa7cc8194
let
	df = df_i1
	gt_array = vec(df[!, :ground_truth_mass])
	calc_array = vec(df[!, :calculated_mass])
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_i1
	model_i1 = lm(@formula(Y ~ X), data)
	global r2_1
	r2_1 = GLM.r2(model_i1)
	global rms_values1
	rms_values1 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_i1))
	]
end

# ╔═╡ dcefe122-007f-4f6d-80d9-1a0c3223e81e
begin
	newX1 = DataFrame(X=collect(1:1000));
	pred_i1 = GLM.predict(model_i1, newX1)
end

# ╔═╡ bf8867ee-ec1c-4aea-b0b1-d8a847b4251b
co1 = coef(model_i1)

# ╔═╡ 360a09c7-a926-4fdf-be8e-ed429dbcde67
md"""
### Agatston
"""

# ╔═╡ e37d1952-c084-46da-afed-00a4bf6541db
df_a1 = CSV.read(joinpath(root_path, "agatston", "motion.csv"), DataFrame); # 14.35

# ╔═╡ 7b1bf319-bca0-487e-b4a8-f7e56ee2d217
let
	df = df_a1
	gt_array = vec(df[!, :ground_truth_mass])
	calc_array = vec(df[!, :calculated_mass])
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_a1
	model_a1 = lm(@formula(Y ~ X), data)
	global r2_2
	r2_2 = GLM.r2(model_a1)
	global rms_values2
	rms_values2 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_a1))
	]
end

# ╔═╡ d6adf4e3-1118-4358-ac59-bd714fd9b244
begin
	newX2 = DataFrame(X=collect(1:1000));
	pred_a1 = GLM.predict(model_a1, newX2)
end

# ╔═╡ 3302d354-a86a-492c-872c-1b1aa0957ba3
co2 = coef(model_a1)

# ╔═╡ c1fce9ed-9cdb-4547-b502-e7d3a23187bd
function lin_reg()
    f = Figure()

    ##-- A --##
	ax = Axis(
		f[1, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Integrated Calcium Scoring",
	)

    df = df_i1
    sc1 = scatter!(df[!, :ground_truth_mass], df[!, :calculated_mass])
    errorbars!(df[!, :ground_truth_mass], df[!, :calculated_mass], rms(df[!, :ground_truth_mass], df[!, :calculated_mass]))
    ln1 = lines!([-1000, 1000], [-1000, 1000])
    ln2 = lines!(collect(1:1000), pred_i1; linestyle=:dashdot)
	create_textbox(f[1, 1], co1, r2_1, rms_values1)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)

    ##-- B --##
	ax = Axis(
		f[2, 1],
		xticks = [0, 50, 100, 150, 200],
		yticks = [0, 50, 100, 150, 200],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Agatston Scoring",
	)

    df3 = df_a1
    sc1 = scatter!(df3[!, :ground_truth_mass], df3[!, :calculated_mass])
    errorbars!(
        df3[!, :ground_truth_mass],
        df3[!, :calculated_mass],
        rms(df3[!, :ground_truth_mass], df3[!, :calculated_mass]),
    )
    ln1 = lines!([-1000, 1000], [-1000, 1000])
    ln2 = lines!(collect(1:1000), pred_a1; linestyle=:dashdot)
	create_textbox(f[2, 1], co2, r2_2, rms_values2)

    xlims!(ax; low=0, high=200)
    ylims!(ax; low=0, high=200)

    ##-- LABELS --##
    f[1:2, 2] = Legend(
        f,
        [sc1, ln1, ln2],
        ["Inserts", "Unity", "Fitted Line"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B"], [f[1, 1], f[2, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            fontsize=25,
            padding=(0, 60, 25, 0),
            halign=:right,
        )
    end

	save(plotsdir("linear_reg_motion.png"), f)
    return f
end

# ╔═╡ 15a2c529-adf3-4d7a-89b7-a80418629507
with_theme(medphys_theme) do
    lin_reg()
end

# ╔═╡ Cell order:
# ╠═345ed6c8-06a1-4ece-8867-6b5a433fe83c
# ╠═600bab60-d643-49a6-94c3-77e0eccc727f
# ╠═7f6c5f2c-58ee-4c53-ba60-02ed0a75cd6f
# ╠═04013b46-69b8-476f-b93c-187622e7883a
# ╠═2f9accf4-0b6d-4f32-9b64-df499fb17d15
# ╟─7e2f297c-c129-4c99-8ae3-63ea6271f3b0
# ╟─c1fce9ed-9cdb-4547-b502-e7d3a23187bd
# ╠═15a2c529-adf3-4d7a-89b7-a80418629507
# ╟─41b322e8-f0b1-404d-94de-6910b0f481a0
# ╠═f65471f2-8f81-46f0-85fe-7aa5b07e1bf6
# ╠═da76517e-6032-4617-bd86-9f500200fe9d
# ╠═2b2adef3-c59c-4a2e-b758-2d2fa7cc8194
# ╠═dcefe122-007f-4f6d-80d9-1a0c3223e81e
# ╠═bf8867ee-ec1c-4aea-b0b1-d8a847b4251b
# ╟─360a09c7-a926-4fdf-be8e-ed429dbcde67
# ╠═e37d1952-c084-46da-afed-00a4bf6541db
# ╠═7b1bf319-bca0-487e-b4a8-f7e56ee2d217
# ╠═d6adf4e3-1118-4358-ac59-bd714fd9b244
# ╠═3302d354-a86a-492c-872c-1b1aa0957ba3
