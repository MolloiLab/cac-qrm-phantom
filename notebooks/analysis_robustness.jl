### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 665791ce-e039-49ed-9eaa-2a064cee02ab
# ╠═╡ show_logs = false
using DrWatson; @quickactivate "cac-qrm-phantom"

# ╔═╡ 0d6339fb-dd53-4834-9720-4d906a6c4720
begin
	using PlutoUI, Statistics, CSV, DataFrames, CairoMakie, Colors, GLM, MLJBase
	using StatsBase: quantile!, rmsd
end

# ╔═╡ 7785d06f-33b7-48db-ab44-df69069be5b9
include(srcdir("plot_utils.jl"));

# ╔═╡ df7b3ae7-86cf-499c-a1f9-e467df65366d
TableOfContents()

# ╔═╡ fe46466b-6fcb-400c-854b-e7048e1c9e2f
root_path = datadir("output")

# ╔═╡ b05002fd-f68c-4f15-b5c1-f2d3f27b40e1
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
		titlesize = 22
	)
);

# ╔═╡ e18acb15-754a-4d2a-9b73-b26b62435c16
md"""
# Robustness
"""

# ╔═╡ 6717dda0-af8f-4a0f-b03c-3213d8ff6fd6
df_i_spec = CSV.read(joinpath(root_path, "integrated", "robustness.csv"), DataFrame);

# ╔═╡ 184b20c0-de97-46a0-ba66-46c2c258b8ae
begin
    df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12 = groupby(df_i_spec, :scan)

    df_thic = df1, df3
    df_kv = df2, df6
    df_ma = df4, df5
    df_fov = df7, df8
    df_ir = df9, df10
    df_qr = df11, df12
end;

# ╔═╡ 7aee5d7c-1aab-4393-936a-5d7e78d78aa4
df_a_spec = CSV.read(joinpath(root_path, "agatston", "robustness.csv"), DataFrame);

# ╔═╡ ae3cfba3-73e3-48fc-8504-760fbbbb09e0
begin
    df1a, df2a, df3a, df4a, df5a, df6a, df7a, df8a, df9a, df10a, df11a, df12a = groupby(df_a_spec, :scan)

    df_thica = df1a, df3a
    df_kva = df2a, df6a
    df_maa = df4a, df5a
    df_fova = df7a, df8a
    df_ira = df9a, df10a
    df_qra = df11a, df12a
end;

# ╔═╡ cc672236-3508-475c-b205-1c4db4b177df
md"""
### Slice thickness
"""

# ╔═╡ f1302e02-eef7-4b18-a6bd-0aeebb33f67d
md"""
#### Integrated
"""

# ╔═╡ 003a1493-3b0b-4239-bbac-55ee258ecc7d
let
    df1 = df_thic[1]
    df2 = df_thic[1]
    gt_array = vec(hcat(df1[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df1[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_thic1i
    model_thic1i = lm(@formula(Y ~ X), data)
    global r2thic1i
    r2thic1i = GLM.r2(model_thic1i)
    global rms_valuesthic1i
    rms_valuesthic1i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_thic1i))
    ]
end

# ╔═╡ 53ac74e7-f808-4e1f-94bf-f1bac854d6a3
begin
    newthic1i = DataFrame(X=collect(1:1000))
    pred_thic1i = GLM.predict(model_thic1i, newthic1i)
end

# ╔═╡ 315a6741-bd98-4247-9759-704d7accef0d
cothic1i = coef(model_thic1i)

# ╔═╡ 4c74f025-0b6a-4990-b482-b348ddba1908
let
    df1 = df_thic[2]
    df2 = df_thic[2]
    gt_array = vec(hcat(df1[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df1[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_thic2i
    model_thic2i = lm(@formula(Y ~ X), data)
    global r2thic2i
    r2thic2i = GLM.r2(model_thic2i)
    global rms_valuesthic2i
    rms_valuesthic2i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_thic2i))
    ]
end

# ╔═╡ 9540ac00-157a-4d1f-b732-be4c790d2cf8
begin
    newthic2i = DataFrame(X=collect(1:1000))
    pred_thic2i = GLM.predict(model_thic2i, newthic2i)
end

# ╔═╡ 69a1a12d-8512-4057-ad34-c3d17af6e411
cothic2i = coef(model_thic2i)

# ╔═╡ 34a52dc9-4802-464b-8955-9481b2323a9a
md"""
#### Agatston
"""

# ╔═╡ 4471308e-edb9-41bf-a7e5-3da80582e6d2
let
    df1 = df_thica[1]
    df2 = df_thica[1]
    gt_array = vec(hcat(df1[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df1[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_thic1a
    model_thic1a = lm(@formula(Y ~ X), data)
    global r2thic1a
    r2thic1a = GLM.r2(model_thic1a)
    global rms_valuesthic1a
    rms_valuesthic1a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_thic1a))
    ]
end

# ╔═╡ 17417e03-e55b-45ff-ae97-609ddf4ee781
begin
    newthic1a = DataFrame(X=collect(1:1000))
    pred_thic1a = GLM.predict(model_thic1a, newthic1a)
end

# ╔═╡ d15a4554-95d4-4fa0-9260-e489403c709d
cothic1a = coef(model_thic1a)

# ╔═╡ 4f9b37c5-200c-4c53-82b0-6f1b03e4d225
let
    df1 = df_thica[2]
    df2 = df_thica[2]
    gt_array = vec(hcat(df1[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df1[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df2[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df2[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_thic2a
    model_thic2a = lm(@formula(Y ~ X), data)
    global r2thic2a
    r2thic2a = GLM.r2(model_thic2a)
    global rms_valuesthic2a
    rms_valuesthic2a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_thic2a))
    ]
end

# ╔═╡ 5d7a1b2d-4a73-4d8a-8c84-06f5ab459b68
begin
    newthic2a = DataFrame(X=collect(1:1000))
    pred_thic2a = GLM.predict(model_thic2a, newthic2a)
end

# ╔═╡ 1c425226-97b7-476b-b5bb-41c785e5d0ca
cothic2a = coef(model_thic2a)

# ╔═╡ 4c90d8df-29a9-416a-9b17-2b5ad6c5434e
function thickness()
    f = Figure()

    ##-- A --##
    axtop = Axis(f[1, 1])

    df = df_thic[1]
    scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtop, collect(1:1000), pred_thic1i, linestyle=:dashdot)
	create_textbox(f[1, 1], cothic1i, r2thic1i, rms_valuesthic1i)

    xlims!(axtop, low=0, high=125)
    ylims!(axtop, low=0, high=125)
    axtop.xticks = [0, 25, 50, 75, 100, 125]
    axtop.yticks = [0, 25, 50, 75, 100, 125]
    axtop.title = "Integrated (1 mm)"

    ##-- B --##
    axtopright = Axis(f[2, 1])

    df3 = df_thic[2]
    scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtopright, collect(1:1000), pred_thic2i, linestyle=:dashdot)
	create_textbox(f[2, 1], cothic2i, r2thic2i, rms_valuesthic2i)

    xlims!(axtopright, low=0, high=125)
    ylims!(axtopright, low=0, high=125)
    axtopright.xticks = [0, 25, 50, 75, 100, 125]
    axtopright.yticks = [0, 25, 50, 75, 100, 125]
    axtopright.title = "Integrated (2 mm)"

    ##-- C --##
    ax3 = Axis(f[1, 2])

    df = df_thica[1]
    scatter!(ax3, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(ax3, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(ax3, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_thic1a, linestyle=:dashdot)
	create_textbox(f[1, 2], cothic1a, r2thic1a, rms_valuesthic1a)

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.title = "Agatston (1 mm)"

    ##-- B --##
    ax4 = Axis(f[2, 2])

    df3 = df_thica[2]
    scatter!(ax4, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_thic2a, linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], cothic2a, r2thic2a, rms_valuesthic2a)

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.title = "Agatston (2 mm)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Slice Thickness",
        fontsize=40)


    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 20, 40, 0),
            halign=:right)
    end

	Label(f[3, 1:2], "Known Mass (mg)", fontsize = 25)
	Label(f[1:2, 0], "Calculated Mass (mg)", fontsize = 25, rotation = pi/2)

	save(plotsdir("linear_reg_slice_thickness.png"), f)
    f
end

# ╔═╡ 56b05673-f0b0-4966-a12b-43a7c878faf0
with_theme(medphys_theme) do
    thickness()
end

# ╔═╡ 1c124a92-6f9b-4295-a53d-39edf1a71cc5
md"""
### Tube Voltage
"""

# ╔═╡ 08020bbd-9dd0-47a0-a92d-88f52d60b7fb
md"""
#### Integrated
"""

# ╔═╡ 2b774e63-84b0-47f0-98bf-1867702955ba
let
    df = df_kv[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_kv1i
    model_kv1i = lm(@formula(Y ~ X), data)
    global r2kv1i
    r2kv1i = GLM.r2(model_kv1i)
    global rms_valueskv1i
    rms_valueskv1i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_kv1i))
    ]
end

# ╔═╡ 9f56e834-056d-4c25-8261-8fd1b95f7c86
begin
    newkv1i = DataFrame(X=collect(1:1000))
    pred_kv1i = GLM.predict(model_kv1i, newkv1i)
end

# ╔═╡ 7c72602b-3a73-44a7-927f-ffdd033778c8
cokv1i = coef(model_kv1i)

# ╔═╡ 9d568aee-9c4b-44a2-9719-4b0bc281625c
let
    df = df_kv[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_kv2i
    model_kv2i = lm(@formula(Y ~ X), data)
    global r2kv2i
    r2kv2i = GLM.r2(model_kv2i)
    global rms_valueskv2i
    rms_valueskv2i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_kv2i))
    ]
end

# ╔═╡ 75240c83-ae5d-4de2-a2c1-3fab86097b95
begin
    newkv2i = DataFrame(X=collect(1:1000))
    pred_kv2i = GLM.predict(model_kv2i, newkv2i)
end

# ╔═╡ b43d05da-4bfb-42e0-91e9-1c873124182c
cokv2i = coef(model_kv2i)

# ╔═╡ df8a19d0-a151-44b6-a73a-8a778a1feb3d
md"""
#### Agatston
"""

# ╔═╡ 8b70e0ee-cc97-4933-ad0a-32b28e044206
let
    df = df_kva[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_kv1a
    model_kv1a = lm(@formula(Y ~ X), data)
    global r2kv1a
    r2kv1a = GLM.r2(model_kv1a)
    global rms_valueskv1a
    rms_valueskv1a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_kv1a))
    ]
end

# ╔═╡ cbc90a45-a4db-4c35-8c47-580c143cb7e9
begin
    newkv1a = DataFrame(X=collect(1:1000))
    pred_kv1a = GLM.predict(model_kv1a, newkv1a)
end

# ╔═╡ 6f1cf07a-5ea8-4c38-a355-c5783a07ef1e
cokv1a = coef(model_kv1a)

# ╔═╡ 9aac0919-bc9d-4059-bef2-691851c346b6
let
    df = df_kva[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_kv2a
    model_kv2a = lm(@formula(Y ~ X), data)
    global r2kv2a
    r2kv2a = GLM.r2(model_kv2a)
    global rms_valueskv2a
    rms_valueskv2a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_kv2a))
    ]
end

# ╔═╡ a3165a8a-a7ba-4c10-9382-4b0ae1fb89f5
begin
    newkv2a = DataFrame(X=collect(1:1000))
    pred_kv2a = GLM.predict(model_kv2a, newkv2a)
end

# ╔═╡ 6e000c8b-df45-4cae-b036-75ec6539566b
cokv2a = coef(model_kv2a)

# ╔═╡ d5485f9d-2941-4d13-8c71-d9efcb8b03c3
function kv()
    f = Figure()

    ##-- A --##
    axtop = Axis(f[1, 1])

    df = df_kv[1]
    scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtop, collect(1:1000), pred_kv1i, linestyle=:dashdot)
	create_textbox(f[1, 1], cokv1i, r2kv1i, rms_valueskv1i)

    xlims!(axtop, low=0, high=125)
    ylims!(axtop, low=0, high=125)
    axtop.xticks = [0, 25, 50, 75, 100, 125]
    axtop.yticks = [0, 25, 50, 75, 100, 125]
    axtop.title = "Integrated (100 kV)"

    ##-- B --##
    axtopright = Axis(f[2, 1])

    df3 = df_kv[2]
    scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtopright, collect(1:1000), pred_kv2i, linestyle=:dashdot)
	create_textbox(f[2, 1], cokv2i, r2kv2i, rms_valueskv2i)

    xlims!(axtopright, low=0, high=125)
    ylims!(axtopright, low=0, high=125)
    axtopright.xticks = [0, 25, 50, 75, 100, 125]
    axtopright.yticks = [0, 25, 50, 75, 100, 125]
    axtopright.title = "Integrated (80 kV)"

    ##-- B --##
    ax3 = Axis(f[1, 2])

    df = df_kva[1]
    scatter!(ax3, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(ax3, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(ax3, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_kv1a, linestyle=:dashdot)
	create_textbox(f[1, 2], cokv1a, r2kv1a, rms_valueskv1a)
	
    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.title = "Agatston (100 kV)"

    ##-- D --##
    ax4 = Axis(f[2, 2])

    df3 = df_kva[2]
    scatter!(ax4, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_kv2a, linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], cokv2a, r2kv2a, rms_valueskv2a)

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.title = "Agatston (80 kV)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Tube Voltage",
        fontsize=40)


    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 20, 40, 0),
            halign=:right)
    end
	Label(f[3, 1:2], "Known Mass (mg)", fontsize = 25)
	Label(f[1:2, 0], "Calculated Mass (mg)", fontsize = 25, rotation = pi/2)

	save(plotsdir("linear_reg_kv.png"), f)
    f
end

# ╔═╡ 876a2c14-0e9c-43a7-be4a-dac5e50fb4d0
with_theme(medphys_theme) do
    kv()
end

# ╔═╡ 573083c6-abe1-4dc3-aa7c-0793c39c1a6a
md"""
### Tube Current Time Product
"""

# ╔═╡ 1827d726-632e-4f29-b448-46a7c4486875
md"""
#### Integrated
"""

# ╔═╡ 0bd89223-4560-47d6-9bc3-4e586184f21e
let
    df = df_ma[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ma1i
    model_ma1i = lm(@formula(Y ~ X), data)
    global r2ma1i
    r2ma1i = GLM.r2(model_ma1i)
    global rms_valuesma1i
    rms_valuesma1i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ma1i))
    ]
end

# ╔═╡ dedc3574-39d6-4152-8de2-973178ef9d59
begin
    newma1i = DataFrame(X=collect(1:1000))
    pred_ma1i = GLM.predict(model_ma1i, newma1i)
end

# ╔═╡ 6be97623-732e-403c-9a16-b0f98de06464
coma1i = coef(model_ma1i)

# ╔═╡ 5a6c043f-c8d2-4dd5-85ae-15c39f67a706
let
    df = df_ma[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ma2i
    model_ma2i = lm(@formula(Y ~ X), data)
    global r2ma2i
    r2ma2i = GLM.r2(model_ma2i)
    global rms_valuesma2i
    rms_valuesma2i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ma2i))
    ]
end

# ╔═╡ f8e5b801-651f-4d0f-9fc1-ed8b86f1a656
begin
    newma2i = DataFrame(X=collect(1:1000))
    pred_ma2i = GLM.predict(model_ma2i, newma2i)
end

# ╔═╡ 54cc9b78-bd7f-4fe0-9492-d5f612515eff
coma2i = coef(model_ma2i)

# ╔═╡ c6b8ec67-ff1a-4615-be0f-012cf7fd6d44
md"""
#### Agatston
"""

# ╔═╡ ba5c4245-c904-4725-8a8f-5adb3d42fc1b
let
    df = df_maa[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ma1a
    model_ma1a = lm(@formula(Y ~ X), data)
    global r2ma1a
    r2ma1a = GLM.r2(model_ma1a)
    global rms_valuesma1a
    rms_valuesma1a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ma1a))
    ]
end

# ╔═╡ 66056723-3c4f-4a45-bd07-b4da248088be
begin
    newma1a = DataFrame(X=collect(1:1000))
    pred_ma1a = GLM.predict(model_ma1a, newma1a)
end

# ╔═╡ 960488a7-29d9-43f4-a0bb-72080a22c782
coma1a = coef(model_ma1a)

# ╔═╡ 15ecfa8b-c199-4cb7-ac40-4afca1927269
let
    df = df_maa[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df1[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df2[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ma2a
    model_ma2a = lm(@formula(Y ~ X), data)
    global r2ma2a
    r2ma2a = GLM.r2(model_ma2a)
    global rms_valuesma2a
    rms_valuesma2a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ma2a))
    ]
end

# ╔═╡ e007c9e4-3339-4df6-ab29-99db3e499175
begin
    newma2a = DataFrame(X=collect(1:1000))
    pred_ma2a = GLM.predict(model_ma2a, newma2a)
end

# ╔═╡ 82e7ad37-da7f-4da1-afec-f99da385d5c7
coma2a = coef(model_ma2a)

# ╔═╡ 0999a0dc-43da-42ae-8b9e-95b2d36a52e5
function ma()
    f = Figure()

    ##-- A --##
    axtop = Axis(f[1, 1])

    df = df_ma[1]
    scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtop, collect(1:1000), pred_ma1i, linestyle=:dashdot)
    Textbox(
        f[1, 1],
        placeholder="y = $(trunc(coma1i[2]; digits=3))x + $(trunc(coma1i[1]; digits=3)) \nr = $(trunc(r2ma1i; digits=3)) \nRMSE: $(trunc(rms_valuesma1i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesma1i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtop, low=0, high=125)
    ylims!(axtop, low=0, high=125)
    axtop.xticks = [0, 25, 50, 75, 100, 125]
    axtop.yticks = [0, 25, 50, 75, 100, 125]
    axtop.xlabel = "Known Mass (mg)"
    axtop.ylabel = "Calculated Mass (mg)"
    axtop.title = "Integrated (22 mAs)"
    # hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)

    ##-- B --##
    axtopright = Axis(f[2, 1])

    df3 = df_ma[2]
    scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtopright, collect(1:1000), pred_ma2i, linestyle=:dashdot)
    Textbox(
        f[2, 1],
        placeholder="y = $(trunc(coma2i[2]; digits=3))x + $(trunc(coma2i[1]; digits=3)) \nr = $(trunc(r2ma2i; digits=3)) \nRMSE: $(trunc(rms_valuesma2i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesma2i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtopright, low=0, high=125)
    ylims!(axtopright, low=0, high=125)
    axtopright.xticks = [0, 25, 50, 75, 100, 125]
    axtopright.yticks = [0, 25, 50, 75, 100, 125]
    axtopright.xlabel = "Known Mass (mg)"
    axtopright.ylabel = "Calculated Mass (mg)"
    axtopright.title = "Integrated (34 mAs)"
    # hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)

    ##-- C --##
    ax3 = Axis(f[1, 2])

    df = df_maa[1]
    scatter!(ax3, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(ax3, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(ax3, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_ma1a, linestyle=:dashdot)
    Textbox(
        f[1, 2],
        placeholder="y = $(trunc(coma1a[2]; digits=3))x + $(trunc(coma1a[1]; digits=3)) \nr = $(trunc(r2ma1a; digits=3)) \nRMSE: $(trunc(rms_valuesma1a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesma1a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Known Mass (mg)"
    ax3.ylabel = "Calculated Mass (mg)"
    ax3.title = "Agatston (22 mAs)"
    # hidedecorations!(ax3, ticklabels=false, ticks=false, label=false)

    ##-- D --##
    ax4 = Axis(f[2, 2])

    df3 = df_maa[2]
    scatter!(ax4, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_ma2a, linestyle=:dashdot, label="Fitted Line")
    Textbox(
        f[2, 2],
        placeholder="y = $(trunc(coma2a[2]; digits=3))x + $(trunc(coma2a[1]; digits=3)) \nr = $(trunc(r2ma2a; digits=3)) \nRMSE: $(trunc(rms_valuesma2a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesma2a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Known Mass (mg)"
    ax4.ylabel = "Calculated Mass (mg)"
    ax4.title = "Agatston (34 mAs)"
    # hidedecorations!(ax4, ticklabels=false, ticks=false, label=false)

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Tube Current Time Product",
        fontsize=40)


    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, -10, 5, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/linear_reg_ma.png", f)
    f
end

# ╔═╡ cf9477ab-db0b-4de6-8980-06b9aba4846a
with_theme(medphys_theme) do
    ma()
end

# ╔═╡ c20aec32-8a49-465d-bf3d-93f8622033f4
md"""
### Field-of-View
"""

# ╔═╡ 38bd180e-f76e-4d31-b43c-c81472aea052
let
    df = df_fov[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_fov1i
    model_fov1i = lm(@formula(Y ~ X), data)
    global r2fov1i
    r2fov1i = GLM.r2(model_fov1i)
    global rms_valuesfov1i
    rms_valuesfov1i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_fov1i))
    ]
end

# ╔═╡ d98011b7-00e7-4c26-9796-8eb9f29b11e2
begin
    newfov1i = DataFrame(X=collect(1:1000))
    pred_fov1i = GLM.predict(model_fov1i, newfov1i)
end

# ╔═╡ 5e6e1f88-648a-46c6-826a-3f21ab5210e0
cofov1i = coef(model_fov1i)

# ╔═╡ d0a8997d-03d3-46e0-bec4-9b8bc8afdc0e
let
    df = df_fov[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_fov2i
    model_fov2i = lm(@formula(Y ~ X), data)
    global r2fov2i
    r2fov2i = GLM.r2(model_fov2i)
    global rms_valuesfov2i
    rms_valuesfov2i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_fov2i))
    ]
end

# ╔═╡ 9b2e099a-236e-479c-aae9-ca44531bdf7c
begin
    newfov2i = DataFrame(X=collect(1:1000))
    pred_fov2i = GLM.predict(model_fov2i, newfov2i)
end

# ╔═╡ e5d7ed17-5dd4-4a24-9d05-198422944445
cofov2i = coef(model_fov2i)

# ╔═╡ 9b5e7744-ecee-42a3-93e2-5ce8db4240df
md"""
#### Agatston
"""

# ╔═╡ cae81b52-7973-4af5-ac4b-72ace8b8d7d1
let
    df = df_fova[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_fov1a
    model_fov1a = lm(@formula(Y ~ X), data)
    global r2fov1a
    r2fov1a = GLM.r2(model_ma1a)
    global rms_valuesfov1a
    rms_valuesfov1a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_fov1a))
    ]
end

# ╔═╡ da2eb942-acb4-4906-bc4b-eb949b6c72e7
begin
    newfov1a = DataFrame(X=collect(1:1000))
    pred_fov1a = GLM.predict(model_fov1a, newfov1a)
end

# ╔═╡ f250bd0b-5c6f-479e-b558-28d73565f6ae
cofov1a = coef(model_fov1a)

# ╔═╡ 90088f65-8c81-4a97-970b-5532744b8982
let
    df = df_fova[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_fov2a
    model_fov2a = lm(@formula(Y ~ X), data)
    global r2fov2a
    r2fov2a = GLM.r2(model_ma2a)
    global rms_valuesfov2a
    rms_valuesfov2a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_fov2a))
    ]
end

# ╔═╡ ee42ed78-0822-4349-a983-9654f57383b0
begin
    newfov2a = DataFrame(X=collect(1:1000))
    pred_fov2a = GLM.predict(model_fov2a, newfov2a)
end

# ╔═╡ 6b1707e5-a19a-4fd6-b042-73c6443184ef
cofov2a = coef(model_fov2a)

# ╔═╡ 5080cf77-a567-4214-9dd5-fdae43b64759
function fov()
    f = Figure()

    ##-- A --##
    axtop = Axis(f[1, 1])

    df = df_fov[1]
    scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtop, collect(1:1000), pred_fov1i, linestyle=:dashdot)
    Textbox(
        f[1, 1],
        placeholder="y = $(trunc(cofov1i[2]; digits=3))x + $(trunc(cofov1i[1]; digits=3)) \nr = $(trunc(r2fov1i; digits=3)) \nRMSE: $(trunc(rms_valuesfov1i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesfov1i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtop, low=0, high=125)
    ylims!(axtop, low=0, high=125)
    axtop.xticks = [0, 25, 50, 75, 100, 125]
    axtop.yticks = [0, 25, 50, 75, 100, 125]
    axtop.xlabel = "Known Mass (mg)"
    axtop.ylabel = "Calculated Mass (mg)"
    axtop.title = "Integrated (FOV 200)"
    # hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)

    ##-- B --##
    axtopright = Axis(f[2, 1])

    df3 = df_fov[2]
    scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtopright, collect(1:1000), pred_fov2i, linestyle=:dashdot)
    Textbox(
        f[2, 1],
        placeholder="y = $(trunc(cofov2i[2]; digits=3))x + $(trunc(cofov2i[1]; digits=3)) \nr = $(trunc(r2fov2i; digits=3)) \nRMSE: $(trunc(rms_valuesfov2i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesfov2i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtopright, low=0, high=125)
    ylims!(axtopright, low=0, high=125)
    axtopright.xticks = [0, 25, 50, 75, 100, 125]
    axtopright.yticks = [0, 25, 50, 75, 100, 125]
    axtopright.xlabel = "Known Mass (mg)"
    axtopright.ylabel = "Calculated Mass (mg)"
    axtopright.title = "Integrated (FOV 320)"
    # hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)

    ##-- C --##
    ax3 = Axis(f[1, 2])

    df = df_fova[1]
    scatter!(ax3, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(ax3, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(ax3, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_fov1a, linestyle=:dashdot)
    Textbox(
        f[1, 2],
        placeholder="y = $(trunc(cofov1a[2]; digits=3))x + $(trunc(cofov1a[1]; digits=3)) \nr = $(trunc(r2fov1a; digits=3)) \nRMSE: $(trunc(rms_valuesfov1a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesfov1a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Known Mass (mg)"
    ax3.ylabel = "Calculated Mass (mg)"
    ax3.title = "Agatston (FOV 200)"
    # hidedecorations!(ax3, ticklabels=false, ticks=false, label=false)

    ##-- D --##
    ax4 = Axis(f[2, 2])

    df3 = df_fov[2]
    scatter!(ax4, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_fov2i, linestyle=:dashdot, label="Fitted Line")
    Textbox(
        f[2, 2],
        placeholder="y = $(trunc(cofov2a[2]; digits=3))x + $(trunc(cofov2a[1]; digits=3)) \nr = $(trunc(r2fov2a; digits=3)) \nRMSE: $(trunc(rms_valuesfov2a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesfov2a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Known Mass (mg)"
    ax4.ylabel = "Calculated Mass (mg)"
    ax4.title = "Agatston (FOV 320)"
    # hidedecorations!(ax4, ticklabels=false, ticks=false, label=false)

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Field-of-View",
        fontsize=40)


    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, -10, 5, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/linear_reg_fov.png", f)
    f
end

# ╔═╡ 3cab42ea-b6f8-493a-abda-08b8e42c5b87
with_theme(medphys_theme) do
    fov()
end

# ╔═╡ f96d6ce0-e583-4994-bd20-1d33fd90a9d4
md"""
### IR
"""

# ╔═╡ 6f7684bb-a22b-46b4-ab0b-2a7df37c07dc
md"""
#### Integrated
"""

# ╔═╡ 47bcbd9b-c9c6-4086-9135-f1d36ed8016a
let
    df = df_ir[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ir1i
    model_ir1i = lm(@formula(Y ~ X), data)
    global r2ir1i
    r2ir1i = GLM.r2(model_ma1i)
    global rms_valuesir1i
    rms_valuesir1i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ir1i))
    ]
end

# ╔═╡ 81b30bc8-ae51-4db6-a998-89293c549762
begin
    newir1i = DataFrame(X=collect(1:1000))
    pred_ir1i = GLM.predict(model_ir1i, newir1i)
end

# ╔═╡ 986467a8-6467-4503-9db0-5788df133b8f
coir1i = coef(model_ir1i)

# ╔═╡ b194d8db-f97b-474f-982d-b1a21decf1ac
let
    df = df_ir[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ir2i
    model_ir2i = lm(@formula(Y ~ X), data)
    global r2ir2i
    r2ir2i = GLM.r2(model_ir2i)
    global rms_valuesir2i
    rms_valuesir2i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ir2i))
    ]
end

# ╔═╡ bdc3bc27-efb3-466a-9d39-ffdb987329d2
begin
    newir2i = DataFrame(X=collect(1:1000))
    pred_ir2i = GLM.predict(model_ir2i, newir2i)
end

# ╔═╡ 44193cc0-4dac-4c95-bdaf-a1d6a9b8c607
coir2i = coef(model_ir2i)

# ╔═╡ af634462-9d22-40ca-82aa-24bdc342e5a8
md"""
#### Agatston
"""

# ╔═╡ e9f45f2f-1229-4227-bfbc-5b2ef066c2c5
let
    df = df_ira[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ir1a
    model_ir1a = lm(@formula(Y ~ X), data)
    global r2ir1a
    r2ir1a = GLM.r2(model_ir1a)
    global rms_valuesir1a
    rms_valuesir1a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ir1a))
    ]
end

# ╔═╡ 45baf37d-3559-4706-a9b2-89b6e044aeef
begin
    newir1a = DataFrame(X=collect(1:1000))
    pred_ir1a = GLM.predict(model_ir1a, newir1a)
end

# ╔═╡ 4d8ee1fb-eb62-4970-a300-edf136dc1d29
coir1a = coef(model_ir1a)

# ╔═╡ 97bc0d26-e753-492c-b81a-649b1a28c4d1
let
    df = df_ira[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_ir2a
    model_ir2a = lm(@formula(Y ~ X), data)
    global r2ir2a
    r2ir2a = GLM.r2(model_ir2a)
    global rms_valuesir2a
    rms_valuesir2a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_ir2a))
    ]
end

# ╔═╡ 1ea8006a-f869-4ae4-8acf-618a6d2c00de
begin
    newir2a = DataFrame(X=collect(1:1000))
    pred_ir2a = GLM.predict(model_ir2a, newir2a)
end

# ╔═╡ f67f8253-28ed-4c74-aa25-b72a04870623
coir2a = coef(model_ir2a)

# ╔═╡ dc4620ac-6e82-4d53-85f8-7c8a99f4ec70
function ir()
    f = Figure()

    ##-- A --##
    axtop = Axis(f[1, 1])

    df = df_ir[1]
    scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtop, collect(1:1000), pred_ir1i, linestyle=:dashdot)
    Textbox(
        f[1, 1],
        placeholder="y = $(trunc(coir1i[2]; digits=3))x + $(trunc(coir1i[1]; digits=3)) \nr = $(trunc(r2ir1i; digits=3)) \nRMSE: $(trunc(rms_valuesir1i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesir1i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtop, low=0, high=125)
    ylims!(axtop, low=0, high=125)
    axtop.xticks = [0, 25, 50, 75, 100, 125]
    axtop.yticks = [0, 25, 50, 75, 100, 125]
    axtop.xlabel = "Known Mass (mg)"
    axtop.ylabel = "Calculated Mass (mg)"
    axtop.title = "Integrated (IR 2)"
    # hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)

    ##-- B --##
    axtopright = Axis(f[2, 1])

    df3 = df_ir[2]
    scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtopright, collect(1:1000), pred_ir2i, linestyle=:dashdot)
    Textbox(
        f[2, 1],
        placeholder="y = $(trunc(coir2i[2]; digits=3))x + $(trunc(coir2i[1]; digits=3)) \nr = $(trunc(r2ir2i; digits=3)) \nRMSE: $(trunc(rms_valuesir2i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesir2i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtopright, low=0, high=125)
    ylims!(axtopright, low=0, high=125)
    axtopright.xticks = [0, 25, 50, 75, 100, 125]
    axtopright.yticks = [0, 25, 50, 75, 100, 125]
    axtopright.xlabel = "Known Mass (mg)"
    axtopright.ylabel = "Calculated Mass (mg)"
    axtopright.title = "Integrated (IR 4)"
    # hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)

    ##-- C --##
    ax3 = Axis(f[1, 2])

    df = df_ira[1]
    scatter!(ax3, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(ax3, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(ax3, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_ir1a, linestyle=:dashdot)
    Textbox(
        f[1, 2],
        placeholder="y = $(trunc(coir1a[2]; digits=3))x + $(trunc(coir1a[1]; digits=3)) \nr = $(trunc(r2ir1a; digits=3)) \nRMSE: $(trunc(rms_valuesir1a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesir1a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Known Mass (mg)"
    ax3.ylabel = "Calculated Mass (mg)"
    ax3.title = "Agatston (IR 2)"
    # hidedecorations!(ax3, ticklabels=false, ticks=false, label=false)

    ##-- D --##
    ax4 = Axis(f[2, 2])

    df3 = df_ira[2]
    scatter!(ax4, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_ir2a, linestyle=:dashdot, label="Fitted Line")
    Textbox(
        f[2, 2],
        placeholder="y = $(trunc(coir2a[2]; digits=3))x + $(trunc(coir2a[1]; digits=3)) \nr = $(trunc(r2ir2a; digits=3)) \nRMSE: $(trunc(rms_valuesir2a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesir2a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Known Mass (mg)"
    ax4.ylabel = "Calculated Mass (mg)"
    ax4.title = "Agatston (IR 4)"
    # hidedecorations!(ax4, ticklabels=false, ticks=false, label=false)

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Iterative Reconstruction Level",
        fontsize=40)


    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, -10, 5, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/linear_reg_ir.png", f)
    f
end

# ╔═╡ c81387d6-0b9d-41ea-9782-90f877fdcb94
with_theme(medphys_theme) do
    ir()
end

# ╔═╡ ce61ccb7-4cfa-4905-b366-65985707a3ad
md"""
### QR32 vs QR44 (Convolutional Kernel)
"""

# ╔═╡ 301db02a-5647-4851-869c-ec48781b01f0
md"""
#### Integrated
"""

# ╔═╡ 0981b1db-db7e-478b-88c0-b63dbab65f58
let
    df = df_qr[1]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_qr1i
    model_qr1i = lm(@formula(Y ~ X), data)
    global r2qr1i
    r2qr1i = GLM.r2(model_ma1i)
    global rms_valuesqr1i
    rms_valuesqr1i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_qr1i))
    ]
end

# ╔═╡ d7115dbe-ff3d-45c7-b63a-a5186f59f0c7
begin
    newqr1i = DataFrame(X=collect(1:1000))
    pred_qr1i = GLM.predict(model_qr1i, newqr1i)
end

# ╔═╡ 19774429-208a-42f3-8e66-645a263a5ce4
coqr1i = coef(model_qr1i)

# ╔═╡ b4e03893-dff5-4180-b16d-061faf6aaa15
let
    df = df_qr[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_qr2i
    model_qr2i = lm(@formula(Y ~ X), data)
    global r2qr2i
    r2qr2i = GLM.r2(model_qr2i)
    global rms_valuesqr2i
    rms_valuesqr2i = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_qr2i))
    ]
end

# ╔═╡ a6b5638b-061e-4bfb-bf60-cbe53773eb30
begin
    newqr2i = DataFrame(X=collect(1:1000))
    pred_qr2i = GLM.predict(model_qr2i, newqr2i)
end

# ╔═╡ d9ab660e-27c7-470f-bbd9-5bfa6f7930c7
coqr2i = coef(model_qr2i)

# ╔═╡ 19d6b2d3-d38d-492e-a37a-90625ead90e6
md"""
#### Agatston
"""

# ╔═╡ 68419afe-ee49-4bbb-8c64-d0f4776bc8b1
let
    df = df_qra[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_qr1a
    model_qr1a = lm(@formula(Y ~ X), data)
    global r2qr1a
    r2qr1a = GLM.r2(model_qr1a)
    global rms_valuesqr1a
    rms_valuesqr1a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_qr1a))
    ]
end

# ╔═╡ 32001554-d378-4c0a-9391-eff3340f363c
begin
    newqr1a = DataFrame(X=collect(1:1000))
    pred_qr1a = GLM.predict(model_qr1a, newqr1a)
end

# ╔═╡ 68c121a1-bff9-44a0-8789-0d762b36af73
coqr1a = coef(model_qr1a)

# ╔═╡ b0920a56-5921-4126-813c-0f4eae8a678b
let
    df = df_qra[2]
    gt_array = vec(hcat(df[!, :ground_truth_mass_large], df[!, :ground_truth_mass_medium], df[!, :ground_truth_mass_small]))
    calc_array = vec(hcat(df[!, :calculated_mass_large], df[!, :calculated_mass_medium], df[!, :calculated_mass_small]))
    data = DataFrame(
        X=gt_array,
        Y=calc_array
    )
    global model_qr2a
    model_qr2a = lm(@formula(Y ~ X), data)
    global r2qr2a
    r2qr2a = GLM.r2(model_qr2a)
    global rms_valuesqr2a
    rms_valuesqr2a = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_qr2a))
    ]
end

# ╔═╡ b1318b91-8175-49fc-90b1-7f17c690c1ac
begin
    newqr2a = DataFrame(X=collect(1:1000))
    pred_qr2a = GLM.predict(model_qr2a, newqr2a)
end

# ╔═╡ 3b6f1284-806a-4f88-bec1-6e00766fe8d7
coqr2a = coef(model_qr2a)

# ╔═╡ c22f8a48-fe76-4773-94e2-36241d06a545
function qr()
    f = Figure()

    ##-- A --##
    axtop = Axis(f[1, 1])

    df = df_qr[1]
    scatter!(axtop, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(axtop, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(axtop, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(axtop, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtop, collect(1:1000), pred_qr1i, linestyle=:dashdot)
    Textbox(
        f[1, 1],
        placeholder="y = $(trunc(coqr1i[2]; digits=3))x + $(trunc(coqr1i[1]; digits=3)) \nr = $(trunc(r2qr1i; digits=3)) \nRMSE: $(trunc(rms_valuesqr1i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesqr1i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtop, low=0, high=125)
    ylims!(axtop, low=0, high=125)
    axtop.xticks = [0, 25, 50, 75, 100, 125]
    axtop.yticks = [0, 25, 50, 75, 100, 125]
    axtop.xlabel = "Known Mass (mg)"
    axtop.ylabel = "Calculated Mass (mg)"
    axtop.title = "Integrated (QR 32)"
    # hidedecorations!(axtop, ticklabels=false, ticks=false, label=false)

    ##-- B --##
    axtopright = Axis(f[2, 1])

    df3 = df_qr[2]
    scatter!(axtopright, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(axtopright, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(axtopright, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(axtopright, collect(1:1000), pred_qr2i, linestyle=:dashdot)
    Textbox(
        f[2, 1],
        placeholder="y = $(trunc(coqr2i[2]; digits=3))x + $(trunc(coqr2i[1]; digits=3)) \nr = $(trunc(r2qr2i; digits=3)) \nRMSE: $(trunc(rms_valuesqr2i[1]; digits=3)) \nRMSD: $(trunc(rms_valuesqr2i[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(axtopright, low=0, high=125)
    ylims!(axtopright, low=0, high=125)
    axtopright.xticks = [0, 25, 50, 75, 100, 125]
    axtopright.yticks = [0, 25, 50, 75, 100, 125]
    axtopright.xlabel = "Known Mass (mg)"
    axtopright.ylabel = "Calculated Mass (mg)"
    axtopright.title = "Integrated (QR 44)"
    # hidedecorations!(axtopright, ticklabels=false, ticks=false, label=false)

    ##-- C --##
    ax3 = Axis(f[1, 2])

    df = df_qra[1]
    scatter!(ax3, df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
    scatter!(ax3, df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
    scatter!(ax3, df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_qr1a, linestyle=:dashdot)
    Textbox(
        f[1, 2],
        placeholder="y = $(trunc(coqr1a[2]; digits=3))x + $(trunc(coqr1a[1]; digits=3)) \nr = $(trunc(r2qr1a; digits=3)) \nRMSE: $(trunc(rms_valuesqr1a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesqr1a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Known Mass (mg)"
    ax3.ylabel = "Calculated Mass (mg)"
    ax3.title = "Agatston (QR 32)"
    # hidedecorations!(ax3, ticklabels=false, ticks=false, label=false)

    ##-- D --##
    ax4 = Axis(f[2, 2])

    df3 = df_qra[2]
    scatter!(ax4, df3[!, :ground_truth_mass_large], df3[!, :calculated_mass_large], label="Large Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_medium], df3[!, :calculated_mass_medium], label="Medium Inserts")
    scatter!(ax4, df3[!, :ground_truth_mass_small], df3[!, :calculated_mass_small], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_qr2a, linestyle=:dashdot, label="Fitted Line")
    Textbox(
        f[2, 2],
        placeholder="y = $(trunc(coqr2a[2]; digits=3))x + $(trunc(coqr2a[1]; digits=3)) \nr = $(trunc(r2qr2a; digits=3)) \nRMSE: $(trunc(rms_valuesqr2a[1]; digits=3)) \nRMSD: $(trunc(rms_valuesqr2a[2]; digits=3))",
        tellheight=false,
        tellwidth=false,
        boxcolor=:white,
        halign=:left,
        valign=:top,
        fontsize=12
    )

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Known Mass (mg)"
    ax4.ylabel = "Calculated Mass (mg)"
    ax4.title = "Agatston (IR 4)"
    # hidedecorations!(ax4, ticklabels=false, ticks=false, label=false)

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Convolution Kernel",
        fontsize=40)


    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, -10, 5, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/linear_reg_qr.png", f)
    f
end

# ╔═╡ 93d64a98-8759-40d5-ad39-9d7cba23f537
with_theme(medphys_theme) do
    qr()
end

# ╔═╡ Cell order:
# ╠═665791ce-e039-49ed-9eaa-2a064cee02ab
# ╠═0d6339fb-dd53-4834-9720-4d906a6c4720
# ╠═7785d06f-33b7-48db-ab44-df69069be5b9
# ╠═df7b3ae7-86cf-499c-a1f9-e467df65366d
# ╠═fe46466b-6fcb-400c-854b-e7048e1c9e2f
# ╠═b05002fd-f68c-4f15-b5c1-f2d3f27b40e1
# ╟─e18acb15-754a-4d2a-9b73-b26b62435c16
# ╠═6717dda0-af8f-4a0f-b03c-3213d8ff6fd6
# ╠═184b20c0-de97-46a0-ba66-46c2c258b8ae
# ╠═7aee5d7c-1aab-4393-936a-5d7e78d78aa4
# ╠═ae3cfba3-73e3-48fc-8504-760fbbbb09e0
# ╟─cc672236-3508-475c-b205-1c4db4b177df
# ╟─f1302e02-eef7-4b18-a6bd-0aeebb33f67d
# ╠═003a1493-3b0b-4239-bbac-55ee258ecc7d
# ╠═53ac74e7-f808-4e1f-94bf-f1bac854d6a3
# ╠═315a6741-bd98-4247-9759-704d7accef0d
# ╠═4c74f025-0b6a-4990-b482-b348ddba1908
# ╠═9540ac00-157a-4d1f-b732-be4c790d2cf8
# ╠═69a1a12d-8512-4057-ad34-c3d17af6e411
# ╟─34a52dc9-4802-464b-8955-9481b2323a9a
# ╠═4471308e-edb9-41bf-a7e5-3da80582e6d2
# ╠═17417e03-e55b-45ff-ae97-609ddf4ee781
# ╠═d15a4554-95d4-4fa0-9260-e489403c709d
# ╠═4f9b37c5-200c-4c53-82b0-6f1b03e4d225
# ╠═5d7a1b2d-4a73-4d8a-8c84-06f5ab459b68
# ╠═1c425226-97b7-476b-b5bb-41c785e5d0ca
# ╟─4c90d8df-29a9-416a-9b17-2b5ad6c5434e
# ╟─56b05673-f0b0-4966-a12b-43a7c878faf0
# ╟─1c124a92-6f9b-4295-a53d-39edf1a71cc5
# ╟─08020bbd-9dd0-47a0-a92d-88f52d60b7fb
# ╠═2b774e63-84b0-47f0-98bf-1867702955ba
# ╠═9f56e834-056d-4c25-8261-8fd1b95f7c86
# ╠═7c72602b-3a73-44a7-927f-ffdd033778c8
# ╠═9d568aee-9c4b-44a2-9719-4b0bc281625c
# ╠═75240c83-ae5d-4de2-a2c1-3fab86097b95
# ╠═b43d05da-4bfb-42e0-91e9-1c873124182c
# ╟─df8a19d0-a151-44b6-a73a-8a778a1feb3d
# ╠═8b70e0ee-cc97-4933-ad0a-32b28e044206
# ╠═cbc90a45-a4db-4c35-8c47-580c143cb7e9
# ╠═6f1cf07a-5ea8-4c38-a355-c5783a07ef1e
# ╠═9aac0919-bc9d-4059-bef2-691851c346b6
# ╠═a3165a8a-a7ba-4c10-9382-4b0ae1fb89f5
# ╠═6e000c8b-df45-4cae-b036-75ec6539566b
# ╟─d5485f9d-2941-4d13-8c71-d9efcb8b03c3
# ╟─876a2c14-0e9c-43a7-be4a-dac5e50fb4d0
# ╟─573083c6-abe1-4dc3-aa7c-0793c39c1a6a
# ╟─1827d726-632e-4f29-b448-46a7c4486875
# ╠═0bd89223-4560-47d6-9bc3-4e586184f21e
# ╠═dedc3574-39d6-4152-8de2-973178ef9d59
# ╠═6be97623-732e-403c-9a16-b0f98de06464
# ╠═5a6c043f-c8d2-4dd5-85ae-15c39f67a706
# ╠═f8e5b801-651f-4d0f-9fc1-ed8b86f1a656
# ╠═54cc9b78-bd7f-4fe0-9492-d5f612515eff
# ╟─c6b8ec67-ff1a-4615-be0f-012cf7fd6d44
# ╠═ba5c4245-c904-4725-8a8f-5adb3d42fc1b
# ╠═66056723-3c4f-4a45-bd07-b4da248088be
# ╠═960488a7-29d9-43f4-a0bb-72080a22c782
# ╠═15ecfa8b-c199-4cb7-ac40-4afca1927269
# ╠═e007c9e4-3339-4df6-ab29-99db3e499175
# ╠═82e7ad37-da7f-4da1-afec-f99da385d5c7
# ╟─0999a0dc-43da-42ae-8b9e-95b2d36a52e5
# ╟─cf9477ab-db0b-4de6-8980-06b9aba4846a
# ╟─c20aec32-8a49-465d-bf3d-93f8622033f4
# ╠═38bd180e-f76e-4d31-b43c-c81472aea052
# ╠═d98011b7-00e7-4c26-9796-8eb9f29b11e2
# ╠═5e6e1f88-648a-46c6-826a-3f21ab5210e0
# ╠═d0a8997d-03d3-46e0-bec4-9b8bc8afdc0e
# ╠═9b2e099a-236e-479c-aae9-ca44531bdf7c
# ╠═e5d7ed17-5dd4-4a24-9d05-198422944445
# ╟─9b5e7744-ecee-42a3-93e2-5ce8db4240df
# ╠═cae81b52-7973-4af5-ac4b-72ace8b8d7d1
# ╠═da2eb942-acb4-4906-bc4b-eb949b6c72e7
# ╠═f250bd0b-5c6f-479e-b558-28d73565f6ae
# ╠═90088f65-8c81-4a97-970b-5532744b8982
# ╠═ee42ed78-0822-4349-a983-9654f57383b0
# ╠═6b1707e5-a19a-4fd6-b042-73c6443184ef
# ╟─5080cf77-a567-4214-9dd5-fdae43b64759
# ╟─3cab42ea-b6f8-493a-abda-08b8e42c5b87
# ╟─f96d6ce0-e583-4994-bd20-1d33fd90a9d4
# ╟─6f7684bb-a22b-46b4-ab0b-2a7df37c07dc
# ╠═47bcbd9b-c9c6-4086-9135-f1d36ed8016a
# ╠═81b30bc8-ae51-4db6-a998-89293c549762
# ╠═986467a8-6467-4503-9db0-5788df133b8f
# ╠═b194d8db-f97b-474f-982d-b1a21decf1ac
# ╠═bdc3bc27-efb3-466a-9d39-ffdb987329d2
# ╠═44193cc0-4dac-4c95-bdaf-a1d6a9b8c607
# ╟─af634462-9d22-40ca-82aa-24bdc342e5a8
# ╠═e9f45f2f-1229-4227-bfbc-5b2ef066c2c5
# ╠═45baf37d-3559-4706-a9b2-89b6e044aeef
# ╠═4d8ee1fb-eb62-4970-a300-edf136dc1d29
# ╠═97bc0d26-e753-492c-b81a-649b1a28c4d1
# ╠═1ea8006a-f869-4ae4-8acf-618a6d2c00de
# ╠═f67f8253-28ed-4c74-aa25-b72a04870623
# ╟─dc4620ac-6e82-4d53-85f8-7c8a99f4ec70
# ╟─c81387d6-0b9d-41ea-9782-90f877fdcb94
# ╟─ce61ccb7-4cfa-4905-b366-65985707a3ad
# ╟─301db02a-5647-4851-869c-ec48781b01f0
# ╠═0981b1db-db7e-478b-88c0-b63dbab65f58
# ╠═d7115dbe-ff3d-45c7-b63a-a5186f59f0c7
# ╠═19774429-208a-42f3-8e66-645a263a5ce4
# ╠═b4e03893-dff5-4180-b16d-061faf6aaa15
# ╠═a6b5638b-061e-4bfb-bf60-cbe53773eb30
# ╠═d9ab660e-27c7-470f-bbd9-5bfa6f7930c7
# ╟─19d6b2d3-d38d-492e-a37a-90625ead90e6
# ╠═68419afe-ee49-4bbb-8c64-d0f4776bc8b1
# ╠═32001554-d378-4c0a-9391-eff3340f363c
# ╠═68c121a1-bff9-44a0-8789-0d762b36af73
# ╠═b0920a56-5921-4126-813c-0f4eae8a678b
# ╠═b1318b91-8175-49fc-90b1-7f17c690c1ac
# ╠═3b6f1284-806a-4f88-bec1-6e00766fe8d7
# ╟─c22f8a48-fe76-4773-94e2-36241d06a545
# ╟─93d64a98-8759-40d5-ad39-9d7cba23f537
