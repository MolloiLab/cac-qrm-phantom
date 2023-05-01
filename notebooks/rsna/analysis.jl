### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 504edef5-02f1-4ddb-8118-4e92a3753791
using DrWatson

# ╔═╡ dc732e3a-53fa-4e89-9f25-8132f521079d
# ╠═╡ show_logs = false
@quickactivate "cac-qrm-phantom"

# ╔═╡ f194debc-4c3a-496d-9934-24a301f28c3d
# ╠═╡ show_logs = false
begin
	using CairoMakie, PlutoUI, Statistics, CSV, DataFrames, Colors, GLM, MLJBase, Printf
	using StatsBase: quantile!, rmsd
end

# ╔═╡ 2b113e97-2324-45ad-bc7e-892e9e659e20
include(srcdir("plot_utils.jl"));

# ╔═╡ d3107e37-c7c7-4013-a367-e433bc73bc78
TableOfContents()

# ╔═╡ c42e5a93-4bd7-42e1-b25f-10cc26b6f89f
root_path = datadir("rsna")

# ╔═╡ 57a14f9b-4d44-4294-9551-872d61b00da5
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

# ╔═╡ 7ee4e678-a93b-434f-aa66-addc508e345f
md"""
# Accuracy
"""

# ╔═╡ 8b9feada-bdcc-4c36-bb2d-62dc39153e54
df_vf = CSV.read(joinpath(root_path, "volume_fraction", "physical.csv"), DataFrame);

# ╔═╡ 9dd9fc84-c9ee-48d1-b444-0cf3ae5d435d
let
	df = df_vf
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
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

# ╔═╡ 5361e63c-a01f-407c-8f40-73184dd711f0
begin
	newX1 = DataFrame(X=collect(1:1000));
	pred_i1 = GLM.predict(model_i1, newX1)
end

# ╔═╡ f302b2b3-da98-4db8-ba4e-5befbf8bd6cd
co1 = coef(model_i1)

# ╔═╡ 2041102a-08f0-46b5-a649-d52570cc0da6
df_a = CSV.read(joinpath(root_path, "agatston", "physical.csv"), DataFrame);

# ╔═╡ a7b86ad5-45ca-4ced-88f4-c23d4dd0df7c
let
	df = df_a
	gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
	calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_a
	model_a = lm(@formula(Y ~ X), data)
	global r2a
	r2a = GLM.r2(model_a)
	global rms_valuesa
	rms_valuesa = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_a))
	]
end

# ╔═╡ 4841200d-64a5-414a-893f-6240436ff01c
begin
	newX3 = DataFrame(X=collect(1:1000));
	pred_a = GLM.predict(model_a, newX3)
end

# ╔═╡ 047deb85-e8f7-4637-9f64-3bac2a44e8b5
co3 = coef(model_a)

# ╔═╡ 8bcaad17-4a5e-4c73-bb4d-04f4372e4143
function accuracy()
	f = Figure()

	##-- A --##
	ax = Axis(
		f[1, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Volume Fraction Calcium Mass",
	)
	
	df = df_vf
	sc1=scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	sc2=scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	sc3=scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	ln1=lines!([-1000, 1000], [-1000, 1000])
	ln2=lines!(collect(1:1000), pred_i1, linestyle=:dashdot)
	create_textbox(f[1, 1], co1, r2_1, rms_values1)
	
	xlims!(ax, low=0, high=125)
	ylims!(ax, low=0, high=125)

	ax = Axis(
		f[2, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Agatson Scoring",
	)

	##-- B --##
	df = df_a
	scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	lines!([-1000, 1000], [-1000, 1000],)
	lines!(collect(1:1000), pred_a, linestyle=:dashdot)
	create_textbox(f[2, 1], co3, r2a, rms_valuesa)
	xlims!(ax, low=0, high=125)
	ylims!(ax, low=0, high=125)

	#-- LABELS --##
	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end
	
	save(plotsdir("rsna", "linear_regression.png"), f)
	f
end

# ╔═╡ 95a57c8f-7d55-4619-bf83-55a1802c3a4e
with_theme(medphys_theme) do
    accuracy()
end

# ╔═╡ 3951d189-d3bb-457e-912f-5ca80dc26bd4
md"""
# Sensitivity and Specificity
"""

# ╔═╡ 71ccd8bf-ac90-4278-9f22-9ac536f33f48
md"""
False-positive expectations
- 0.0 std => ~50% false-positive
- 0.5 std => ~30% false-positive
- 1.0 std => ~16% false-positive
- 1.5 std => ~7% false-positive
- 2.0 std => ~2% false-positive
"""

# ╔═╡ 3311fd03-36f2-41b5-b261-4df42c448a36
std_level = 1.5

# ╔═╡ d186f370-48ab-408f-a1b7-e6a8a61f6b31
md"""
## False Negative
"""

# ╔═╡ f142b1b1-e861-48f1-bc32-6e8cf00ac716
md"""
#### Agatston
"""

# ╔═╡ 1d072000-ed66-4a96-9ee8-9bdadafc92cc
array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small]);

# ╔═╡ 11d03fd2-a1ff-4c3f-95c5-08027013bcf2
total_cac = length(array_a)

# ╔═╡ a20dbe6b-1c1b-4d03-ac10-a8ee9c935064
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ 2d65ee51-2771-42ef-83eb-a5ab740ddbf1
length(findall(x -> x <= 0, df_a[!, :calculated_mass_large])), length(findall(x -> x <= 0, df_a[!, :calculated_mass_medium])), length(findall(x -> x <= 0, df_a[!, :calculated_mass_small]))

# ╔═╡ d35cfb16-3938-4334-833a-e5e8f41741b6
40/360, 38/360, 25/360

# ╔═╡ c06ea1d6-83ef-41b5-a6eb-34aff7b89287
begin
	df_a_ld, df_a_md, df_a_hd = groupby(df_a, :inserts)
	
	length(findall(x -> x <= 0, hcat(df_a_ld[!, :calculated_mass_large], df_a_ld[!, :calculated_mass_medium], df_a_ld[!, :calculated_mass_small]))), length(findall(x -> x <= 0, hcat(df_a_md[!, :calculated_mass_large], df_a_md[!, :calculated_mass_medium], df_a_md[!, :calculated_mass_small]))), length(findall(x -> x <= 0, hcat(df_a_hd[!, :calculated_mass_large], df_a_hd[!, :calculated_mass_medium], df_a_hd[!, :calculated_mass_small])))
end

# ╔═╡ 0a781663-cf52-4131-9f86-7054f9ee8284
md"""
#### Volume Fraction
"""

# ╔═╡ 8312248f-8367-4a11-a9c8-5431695532b3
begin
	false_negative_i = []
	for i in 1:3:nrow(df_vf)-2
		mean_i, std_i = mean(df_vf[i:i+2, :mass_bkg]), std(df_vf[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_vf[i:i+2, :calculated_mass_large], df_vf[i:i+2, :calculated_mass_medium], df_vf[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i, neg)
	end
end

# ╔═╡ 3b47b01a-5f6d-4e62-bb6d-64a8b8b68aad
total_zero_i = sum(false_negative_i)

# ╔═╡ f4eb27f6-578e-46f2-b32e-314ab27e62b8
md"""
#### SWCS
"""

# ╔═╡ 089b38be-b693-4e7b-87f2-35999827c1c5
df_s = CSV.read(joinpath(root_path, "swcs", "physical.csv"), DataFrame);

# ╔═╡ 17aecabb-752b-4db5-b01e-18da6b2d4f01
begin
	false_negative_s = []
	for i in 1:3:nrow(df_s)-2
		mean_i, std_i = mean(df_s[i:i+2, :swcs_bkg]), std(df_s[i:i+2, :swcs_bkg])*std_level 
		array_i = hcat(df_s[i:i+2, :calculated_swcs_large], df_s[i:i+2, :calculated_swcs_medium], df_s[i:i+2, :calculated_swcs_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_s, neg)
	end
end

# ╔═╡ 4ed53ec4-aa12-470b-b6d9-8bc7acdf0ae6
total_zero_s = sum(false_negative_s)

# ╔═╡ e75a9cc0-03ea-4de9-9c5e-b419bb30b68c
total_zero_i, total_zero_s, num_zero_a

# ╔═╡ 7cc4e287-1240-42d6-9c7b-b5f3601b3c12
md"""
## False Positive
"""

# ╔═╡ 072feb69-adb1-408e-a920-23a409205703
md"""
#### Agatston
"""

# ╔═╡ 57b2d301-d7c8-4b87-9b1f-9f751984d626
array_a_pos = df_a[!, :mass_bkg]

# ╔═╡ a459f8db-a2cb-4b27-9f47-9096504695c6
total_cac_pos = length(array_a_pos)

# ╔═╡ 36480397-b605-4cea-b63f-34bd86861a98
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ d33126d5-0c4f-4497-ade9-82701e5c6ef8
md"""
#### Volume Fraction
"""

# ╔═╡ aa6214d2-4e17-4667-a88f-c05d11b10149
begin
	false_positive_i = []
	for i in 1:3:nrow(df_vf)-2
		mean_i, std_i = mean(df_vf[i:i+2, :mass_bkg]), std(df_vf[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_vf[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i, pos)
	end
end

# ╔═╡ be13830f-abe0-448e-86eb-f47f0954bc75
total_zero_i_pos = sum(false_positive_i)

# ╔═╡ 93e7669a-d446-4877-ac5c-23717edaddc3
md"""
#### SWCS
"""

# ╔═╡ 4a15a604-6799-4185-8b69-ea140dfbafbc
begin
	false_positive_s = []
	for i in 1:3:nrow(df_s)-2
		mean_i, std_i = mean(df_s[i:i+2, :swcs_bkg]), std(df_s[i:i+2, :swcs_bkg])*std_level
		array_i_pos = df_s[i:i+2, :swcs_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_s, pos)
	end
end

# ╔═╡ 80c575d4-6ea4-46e9-bd63-e78dc2c5dcd7
total_zero_s_pos = sum(false_positive_s)

# ╔═╡ 018331c8-ee76-475f-80d8-1795c9f4a857
function sensitivity_specificity()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:3, ["Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_i / total_cac) * 100
	h2 = (total_zero_s / total_cac) * 100
	h3 = (num_zero_a / total_cac) * 100 
	heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[1:3], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	##-- B --##
	ax = Axis(
		f[2, 1]; 
		xticks = (1:3, ["Volume Fraction", "Spatially Weighted", "Agatston"]),
		title = "False-Positive (CAC>0)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_i_pos / total_cac_pos) * 100
	h2 = (total_zero_s_pos / total_cac_pos) * 100
	h3 = (total_zero_a_pos / total_cac_pos) * 100
    heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[1:3], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save(plotsdir("rsna", "sensitivity_specificity.png"), f)

    return f
end

# ╔═╡ 6c476510-3168-4979-ab42-eb9cbff2f971
with_theme(medphys_theme) do
    sensitivity_specificity()
end

# ╔═╡ e8a0a64f-bef3-414d-b7c0-0c3d5122f8d1
total_zero_i_pos, total_zero_s_pos, total_zero_a_pos

# ╔═╡ Cell order:
# ╠═504edef5-02f1-4ddb-8118-4e92a3753791
# ╠═dc732e3a-53fa-4e89-9f25-8132f521079d
# ╠═f194debc-4c3a-496d-9934-24a301f28c3d
# ╠═2b113e97-2324-45ad-bc7e-892e9e659e20
# ╠═d3107e37-c7c7-4013-a367-e433bc73bc78
# ╠═c42e5a93-4bd7-42e1-b25f-10cc26b6f89f
# ╠═57a14f9b-4d44-4294-9551-872d61b00da5
# ╟─7ee4e678-a93b-434f-aa66-addc508e345f
# ╠═8b9feada-bdcc-4c36-bb2d-62dc39153e54
# ╠═9dd9fc84-c9ee-48d1-b444-0cf3ae5d435d
# ╠═5361e63c-a01f-407c-8f40-73184dd711f0
# ╠═f302b2b3-da98-4db8-ba4e-5befbf8bd6cd
# ╠═2041102a-08f0-46b5-a649-d52570cc0da6
# ╠═a7b86ad5-45ca-4ced-88f4-c23d4dd0df7c
# ╠═4841200d-64a5-414a-893f-6240436ff01c
# ╠═047deb85-e8f7-4637-9f64-3bac2a44e8b5
# ╟─8bcaad17-4a5e-4c73-bb4d-04f4372e4143
# ╟─95a57c8f-7d55-4619-bf83-55a1802c3a4e
# ╟─3951d189-d3bb-457e-912f-5ca80dc26bd4
# ╟─71ccd8bf-ac90-4278-9f22-9ac536f33f48
# ╠═3311fd03-36f2-41b5-b261-4df42c448a36
# ╟─d186f370-48ab-408f-a1b7-e6a8a61f6b31
# ╟─f142b1b1-e861-48f1-bc32-6e8cf00ac716
# ╠═1d072000-ed66-4a96-9ee8-9bdadafc92cc
# ╠═11d03fd2-a1ff-4c3f-95c5-08027013bcf2
# ╠═a20dbe6b-1c1b-4d03-ac10-a8ee9c935064
# ╠═2d65ee51-2771-42ef-83eb-a5ab740ddbf1
# ╠═d35cfb16-3938-4334-833a-e5e8f41741b6
# ╠═c06ea1d6-83ef-41b5-a6eb-34aff7b89287
# ╟─0a781663-cf52-4131-9f86-7054f9ee8284
# ╠═8312248f-8367-4a11-a9c8-5431695532b3
# ╠═3b47b01a-5f6d-4e62-bb6d-64a8b8b68aad
# ╠═e75a9cc0-03ea-4de9-9c5e-b419bb30b68c
# ╟─f4eb27f6-578e-46f2-b32e-314ab27e62b8
# ╠═17aecabb-752b-4db5-b01e-18da6b2d4f01
# ╠═4ed53ec4-aa12-470b-b6d9-8bc7acdf0ae6
# ╠═089b38be-b693-4e7b-87f2-35999827c1c5
# ╟─7cc4e287-1240-42d6-9c7b-b5f3601b3c12
# ╟─072feb69-adb1-408e-a920-23a409205703
# ╠═57b2d301-d7c8-4b87-9b1f-9f751984d626
# ╠═a459f8db-a2cb-4b27-9f47-9096504695c6
# ╠═36480397-b605-4cea-b63f-34bd86861a98
# ╟─d33126d5-0c4f-4497-ade9-82701e5c6ef8
# ╠═aa6214d2-4e17-4667-a88f-c05d11b10149
# ╠═be13830f-abe0-448e-86eb-f47f0954bc75
# ╟─93e7669a-d446-4877-ac5c-23717edaddc3
# ╠═4a15a604-6799-4185-8b69-ea140dfbafbc
# ╠═80c575d4-6ea4-46e9-bd63-e78dc2c5dcd7
# ╟─018331c8-ee76-475f-80d8-1795c9f4a857
# ╟─6c476510-3168-4979-ab42-eb9cbff2f971
# ╠═e8a0a64f-bef3-414d-b7c0-0c3d5122f8d1
