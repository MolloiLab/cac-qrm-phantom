### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 40bf5e39-0b70-47c1-82d8-3033cd1ed60b
# ╠═╡ show_logs = false
using DrWatson; @quickactivate "cac-qrm-phantom"

# ╔═╡ 47b98f40-4522-4bc4-a5de-c2e99e443143
begin
	using PlutoUI, Statistics, CSV, DataFrames, CairoMakie, Colors, GLM, MLJBase
	using StatsBase: quantile!, rmsd
end

# ╔═╡ 3d7dcd97-4202-4d45-8791-cd22099d1b8c
include(srcdir("plot_utils.jl"));

# ╔═╡ 860c3023-55be-4cd5-924a-1a677d6b9281
TableOfContents()

# ╔═╡ 86b24b18-b24c-48b1-b89c-2bde53744fe5
root_path = datadir("output")

# ╔═╡ 2838451c-be95-41f7-813b-7323d4a621b4
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

# ╔═╡ 7bf364d8-cb82-4f79-8c53-0a426eb09d12
scan_names = ["large1", "large2", "large3", "large4", "large5", "small1", "small2", "small3", "small4", "small5"]

# ╔═╡ 0520c367-4ecf-4c31-a396-6b70e1a630a5
std_level = 1.5

# ╔═╡ 960ad1f4-83a9-4e78-9d2e-1517a75cb127
md"""
## True FN removal
"""

# ╔═╡ 097f3636-6fe6-41d7-ae49-dd761d3cce86
begin
    df_a = CSV.read(joinpath(root_path, "agatston", "physical.csv"), DataFrame)

    # Rename :scan
    for i in 1:120
        j = Int(round(i / 3, RoundUp))
        if i > 30 && i ≤ 60
            j = Int(round((i - 30) / 3, RoundUp))
        end
        if i > 60 && i ≤ 90
            j = Int(round((i - 60) / 3, RoundUp))
        end
        if i > 90 && i ≤ 120
            j = Int(round((i - 90) / 3, RoundUp))
        end
        new_name = scan_names[j]
        df_a[i, :scan] = new_name
    end

    # Remove False Negatives
    array_a = hcat(df_a[!, :calculated_mass_large], df_a[!, :calculated_mass_medium], df_a[!, :calculated_mass_small])

    idxs_a_large = Tuple(findall(x -> x == 0, array_a[:, 1]))
    idxs_a_med = Tuple(findall(x -> x == 0, array_a[:, 2]))
    idxs_a_small = Tuple(findall(x -> x == 0, array_a[:, 3]))

    for i in idxs_a_large
        df_a[i, :calculated_mass_large] = NaN
    end
    for i in idxs_a_med
        df_a[i, :calculated_mass_medium] = NaN
    end
    for i in idxs_a_small
        df_a[i, :calculated_mass_small] = NaN
    end

    df_a_large = hcat(df_a[:, 1:3], df_a[:, :calculated_mass_large])
    df_a_medium = hcat(df_a[:, 1:3], df_a[:, :calculated_mass_medium])
    df_a_small = hcat(df_a[:, 1:3], df_a[:, :calculated_mass_small])
end

# ╔═╡ e848a449-2983-4ead-8bad-f5a8fa735b69
begin
    df_i = CSV.read(joinpath(root_path, "integrated", "physical.csv"), DataFrame)

    # Rename :scan
    for i in 1:120
        j = Int(round(i / 3, RoundUp))
        if i > 30 && i ≤ 60
            j = Int(round((i - 30) / 3, RoundUp))
        end
        if i > 60 && i ≤ 90
            j = Int(round((i - 60) / 3, RoundUp))
        end
        if i > 90 && i ≤ 120
            j = Int(round((i - 90) / 3, RoundUp))
        end
        new_name = scan_names[j]
        df_i[i, :scan] = new_name
    end

    # Remove False Negatives
    array_i = hcat(df_i[!, :calculated_mass_large], df_i[!, :calculated_mass_medium], df_i[!, :calculated_mass_small])

    idxs_i_large = Tuple(findall(x -> x < 0, array_i[:, 1]))
    idxs_i_med = Tuple(findall(x -> x < 0, array_i[:, 2]))
    idxs_i_small = Tuple(findall(x -> x < 0, array_i[:, 3]))

    for i in idxs_a_large
        df_i[i, :calculated_mass_large] = NaN
    end
    for i in idxs_a_med
        df_i[i, :calculated_mass_medium] = NaN
    end
    for i in idxs_a_small
        df_i[i, :calculated_mass_small] = NaN
    end

    df_i_large = hcat(df_i[:, 1:3], df_i[:, :calculated_mass_large])
    df_i_medium = hcat(df_i[:, 1:3], df_i[:, :calculated_mass_medium])
    df_i_small = hcat(df_i[:, 1:3], df_i[:, :calculated_mass_small])
end

# ╔═╡ 44576c72-f38a-48de-93bc-e4b79ba1a9a3
idxs_i_small

# ╔═╡ 808d5303-7401-43b0-bc1f-42578770e668
begin
	false_negative_i = []
	for i in 1:3:nrow(df_i)-2
		mean_i, std_i = mean(df_i[i:i+2, :mass_bkg]), std(df_i[i:i+2, :mass_bkg])*std_level
		arr1 = hcat(df_i[i, :calculated_mass_large], df_i[i, :calculated_mass_medium], df_i[i, :calculated_mass_small])
		arr2 = hcat(df_i[i+1, :calculated_mass_large], df_i[i+1, :calculated_mass_medium], df_i[i+1, :calculated_mass_small])
		arr3 = hcat(df_i[i+2, :calculated_mass_large], df_i[i+2, :calculated_mass_medium], df_i[i+2, :calculated_mass_small])
		
		neg1 = findall(x -> x <= mean_i + std_i, arr1)
		neg2 = findall(x -> x <= mean_i + std_i, arr2)
		neg3 = findall(x -> x <= mean_i + std_i, arr3)

		if length(neg1) == 0
			push!(false_negative_i, 0)
		else
			push!(false_negative_i, neg1...)
		end

		if length(neg2) == 0
			push!(false_negative_i, 0)
		else
			push!(false_negative_i, neg2...)
		end

		if length(neg3) == 0
			push!(false_negative_i, 0)
		else
			push!(false_negative_i, neg3...)
		end
	end
end

# ╔═╡ 094f4999-30ab-4f7c-b7bc-accaf6c07673
begin
	negs = []
	for i in 1:length(false_negative_i)
		if false_negative_i[i] != 0
			ci = i, false_negative_i[i][2]
			push!(negs, CartesianIndex(ci...))
		end
	end
end

# ╔═╡ 16285410-d2ad-41de-b383-531f0ba1ee96
begin
	idx_i_large = []
	idx_i_med = []
	idx_i_small = []
	for i in negs
		if i[2] == 1
			push!(idx_i_large, i[1])
		elseif i[2] == 2
			push!(idx_i_med, i[1])
		elseif i[2] == 3
			push!(idx_i_small, i[1])
		end
	end
end

# ╔═╡ 9260e6cc-20e3-477c-851c-aae4a0cd5677
idx_i_small

# ╔═╡ bcbca806-78b9-4147-b068-b815201ff53d
md"""
## Group Inserts by Scan
"""

# ╔═╡ f7d7d226-1e40-4921-9af0-c065e0a6af30
begin
    # Large
    df_larges = groupby(df_i_large, :scan)
    clean_dfs = []
    for df in df_larges
        push!(clean_dfs, filter(:x1 => x -> !(isnan(x)), df))
    end
    df1_large, df2_large, df3_large, df4_large, df5_large, df6_large, df7_large, df8_large, df9_large, df10_large = clean_dfs

    # Medium
    df_mediums = groupby(df_i_medium, :scan)
    clean_dfs = []
    for df in df_mediums
        push!(clean_dfs, filter(:x1 => x -> !(isnan(x)), df))
    end
    df1_medium, df2_medium, df3_medium, df4_medium, df5_medium, df6_medium, df7_medium, df8_medium, df9_medium, df10_medium = clean_dfs

    #Small
    df_smalls = groupby(df_i_small, :scan)
    clean_dfs = []
    for df in df_smalls
        push!(clean_dfs, filter(:x1 => x -> !(isnan(x)), df))
    end
    df1_small, df2_small, df3_small, df4_small, df5_small, df6_small, df7_small, df8_small, df9_small, df10_small = clean_dfs
end;

# ╔═╡ 8745cff1-1755-4e2d-9801-21f8af8e0c36
begin
    df1_2_large_clean = innerjoin(df1_large, df2_large; on=[:vendor, :inserts], makeunique=true)
    df1_2_medium_clean = innerjoin(df1_medium, df2_medium; on=[:vendor, :inserts], makeunique=true)
    df1_2_small_clean = innerjoin(df1_small, df2_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 1afc780f-a4e0-4526-bf22-49f24eeb90ed
begin
    # Large
    df_larges_a = groupby(df_a_large, :scan)
    clean_dfs_a = []
    for df in df_larges_a
        push!(clean_dfs_a, filter(:x1 => x -> !(isnan(x)), df))
    end
    df1_large_a, df2_large_a, df3_large_a, df4_large_a, df5_large_a, df6_large_a, df7_large_a, df8_large_a, df9_large_a, df10_large_a = clean_dfs_a

    # Medium
    df_mediums_a = groupby(df_a_medium, :scan)
    clean_dfs_a = []
    for df in df_mediums_a
        push!(clean_dfs_a, filter(:x1 => x -> !(isnan(x)), df))
    end
    df1_medium_a, df2_medium_a, df3_medium_a, df4_medium_a, df5_medium_a, df6_medium_a, df7_medium_a, df8_medium_a, df9_medium_a, df10_medium_a = clean_dfs_a

    #Small
    df_smalls_a = groupby(df_a_small, :scan)
    clean_dfs_a = []
    for df in df_smalls_a
        push!(clean_dfs_a, filter(:x1 => x -> !(isnan(x)), df))
    end
    df1_small_a, df2_small_a, df3_small_a, df4_small_a, df5_small_a, df6_small_a, df7_small_a, df8_small_a, df9_small_a, df10_small_a = clean_dfs_a
end;

# ╔═╡ 7a2f23bb-c548-4410-96c0-6fa7e1630784
md"""
# Large Phantom - Integrated
"""

# ╔═╡ 8f4151a3-3baf-4126-97e7-fafb6a92b2d9
md"""
## 1 vs 2
"""

# ╔═╡ 2efd3f1c-2d88-4689-ba02-0cf66941357b
begin
    df1_2_large_clean_canon, df1_2_large_clean_geee, df1_2_large_clean_philips, df1_2_large_clean_siemens = groupby(df1_2_large_clean, :vendor)
    df1_2_medium_clean_canon, df1_2_medium_clean_geee, df1_2_medium_clean_philips, df1_2_medium_clean_siemens = groupby(df1_2_medium_clean, :vendor)
    df1_2_small_clean_canon = df1_2_small_clean
end;

# ╔═╡ 701f24a7-c333-40e7-a2a9-271fc83aa54b
md"""
### Canon
"""

# ╔═╡ a79fcdb3-a30b-4bda-9bea-f012550ec14d
let
    r_array = [Tuple(df1_2_large_clean_canon[!, :x1])..., Tuple(df1_2_medium_clean_canon[!, :x1])..., Tuple(df1_2_small_clean_canon[!, :x1])...]
    calc_array = [Tuple(df1_2_large_clean_canon[!, :x1_1])..., Tuple(df1_2_medium_clean_canon[!, :x1_1])..., Tuple(df1_2_small_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_1_canon
    model_i_r_1_canon = lm(@formula(Y ~ X), data)
    global r_i_r_1_canon
    r_i_r_1_canon = GLM.r2(model_i_r_1_canon)
    global rms_values_i_r_1_canon
    rms_values_i_r_1_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_1_canon))
    ]
end

# ╔═╡ e5478856-fae3-481c-8565-66cc1bbd47ed
begin
    X_i_r_1_canon = DataFrame(X=collect(1:1000))
    pred_i_r_1_canon = GLM.predict(model_i_r_1_canon, X_i_r_1_canon)
end

# ╔═╡ 39a7ac34-4fa1-4cd7-b736-aa6b9c05fb6e
co_i_r_1_canon = coef(model_i_r_1_canon)

# ╔═╡ 3d40a8a2-c8a9-4ba5-a9a7-9bdbba7914b6
md"""
### GE
"""

# ╔═╡ 118ad8fd-4423-4f90-a6ec-5a6a89a40637
let
    r_array = [Tuple(df1_2_large_clean_geee[!, :x1])..., Tuple(df1_2_medium_clean_geee[!, :x1])...]
    calc_array = [Tuple(df1_2_large_clean_geee[!, :x1_1])..., Tuple(df1_2_medium_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_1_geee
    model_i_r_1_geee = lm(@formula(Y ~ X), data)
    global r_i_r_1_geee
    r_i_r_1_geee = GLM.r2(model_i_r_1_geee)
    global rms_values_i_r_1_geee
    rms_values_i_r_1_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_1_geee))
    ]
end

# ╔═╡ ec158a7b-2cd3-488c-b192-74f69e34d9b7
begin
    X_i_r_1_geee = DataFrame(X=collect(1:1000))
    pred_i_r_1_geee = GLM.predict(model_i_r_1_geee, X_i_r_1_geee)
end

# ╔═╡ fac3bf90-967d-4ac8-b223-8910c3e0a423
co_i_r_1_geee = coef(model_i_r_1_geee)

# ╔═╡ 0be9d780-1097-4d23-a745-ffe0b30e1173
md"""
### Philips
"""

# ╔═╡ c0e6bdd0-f780-4af3-aa18-9878c05ab8fd
let
    r_array = [Tuple(df1_2_large_clean_philips[!, :x1])..., Tuple(df1_2_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df1_2_large_clean_philips[!, :x1_1])..., Tuple(df1_2_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_1_philips
    model_i_r_1_philips = lm(@formula(Y ~ X), data)
    global r_i_r_1_philips
    r_i_r_1_philips = GLM.r2(model_i_r_1_philips)
    global rms_values_i_r_1_philips
    rms_values_i_r_1_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_1_philips))
    ]
end

# ╔═╡ 3713e966-9bf1-4577-8b79-0ecd6ed925d1
begin
    X_i_r_1_philips = DataFrame(X=collect(1:1000))
    pred_i_r_1_philips = GLM.predict(model_i_r_1_philips, X_i_r_1_philips)
end

# ╔═╡ e307f779-b1da-43b6-a505-36e817984cb1
co_i_r_1_philips = coef(model_i_r_1_philips)

# ╔═╡ 0be916a3-5fe5-488e-91d2-81f9267b7b21
md"""
### Siemens
"""

# ╔═╡ b4a86e37-8ba7-4f11-b48c-8ee3cf43f981
let
    r_array = [Tuple(df1_2_large_clean_siemens[!, :x1])..., Tuple(df1_2_medium_clean_siemens[!, :x1])...,]
    calc_array = [Tuple(df1_2_large_clean_siemens[!, :x1_1])..., Tuple(df1_2_medium_clean_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_1_siemens
    model_i_r_1_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_1_siemens
    r_i_r_1_siemens = GLM.r2(model_i_r_1_siemens)
    global rms_values_i_r_1_siemens
    rms_values_i_r_1_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_1_siemens))
    ]
end

# ╔═╡ e596f2ed-7826-44b7-8cad-4625c9e40365
begin
    X_i_r_1_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_1_siemens = GLM.predict(model_i_r_1_siemens, X_i_r_1_siemens)
end

# ╔═╡ 4bd66f83-427f-46d3-b5e0-7d1c2d2a9cf8
co_i_r_1_siemens = coef(model_i_r_1_siemens)

# ╔═╡ f2d3a05e-17d9-4dbc-a399-8ac69dfb37c8
md"""
### Visualize
"""

# ╔═╡ d6c75439-e2bf-404e-8641-6820f88d09b1
function reprod_i_1_2()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df1_2_large_clean_canon[!, :x1], df1_2_large_clean_canon[!, :x1_1])
    scatter!(ax1, df1_2_medium_clean_canon[!, :x1], df1_2_medium_clean_canon[!, :x1_1])
    scatter!(ax1, df1_2_small_clean_canon[!, :x1], df1_2_small_clean_canon[!, :x1_1], color=:red)
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_1_canon, linestyle=:dashdot)
	create_textbox(f[1, 1], co_i_r_1_canon, r_i_r_1_canon, rms_values_i_r_1_canon)

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "Scanner 1"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_2_large_clean_geee[!, :x1], df1_2_large_clean_geee[!, :x1_1])
    scatter!(ax2, df1_2_medium_clean_geee[!, :x1], df1_2_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df1_2_small_clean_geee[!, :x1], df1_2_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_1_geee, linestyle=:dashdot)
	create_textbox(f[2, 1], co_i_r_1_geee, r_i_r_1_geee, rms_values_i_r_1_geee)

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "Scanner 2"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_2_large_clean_philips[!, :x1], df1_2_large_clean_philips[!, :x1_1])
    scatter!(ax3, df1_2_medium_clean_philips[!, :x1], df1_2_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df1_2_small_clean_philips[!, :x1], df1_2_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_1_philips, linestyle=:dashdot)
	create_textbox(f[1, 2], co_i_r_1_philips, r_i_r_1_philips, rms_values_i_r_1_philips)

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "Scanner 3"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_2_large_clean_siemens[!, :x1], df1_2_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_2_medium_clean_siemens[!, :x1], df1_2_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    # scatter!(ax4, df1_2_small_clean_siemens[!, :x1], df1_2_small_clean_siemens[!, :x1_1], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_1_siemens, linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], co_i_r_1_siemens, r_i_r_1_siemens, rms_values_i_r_1_siemens)

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "Scanner 4"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (1 vs 2)",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end
	
	save(plotsdir("repro_i_1_2.png"), f)
    f
end

# ╔═╡ 343e5452-3da7-4ed1-b1af-d281ef17192e
with_theme(medphys_theme) do
    reprod_i_1_2()
end

# ╔═╡ abc07dfe-477c-431f-92ad-ef51ac76b074
md"""
## 1 vs 3
"""

# ╔═╡ e7ad474a-9356-4ea3-a8b6-9d2397d9c4a3
begin
    df1_3_large_clean = innerjoin(df1_large, df3_large; on=[:vendor, :inserts], makeunique=true)
    df1_3_medium_clean = innerjoin(df1_medium, df3_medium; on=[:vendor, :inserts], makeunique=true)
    df1_3_small_clean = innerjoin(df1_small, df3_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 1be99b59-792b-40cc-ac8e-7cb2b19a45aa
md"""
### Canon
"""

# ╔═╡ 953444d5-ef00-4922-bbb9-9b2c64ab0e92
begin
    df1_3_large_clean_canon, df1_3_large_clean_geee, df1_3_large_clean_philips, df1_3_large_clean_siemens = groupby(df1_3_large_clean, :vendor)
    df1_3_medium_clean_canon, df1_3_medium_clean_geee, df1_3_medium_clean_philips, df1_3_medium_clean_siemens = groupby(df1_3_medium_clean, :vendor)
    df1_3_small_clean_canon = df1_3_small_clean
end;

# ╔═╡ 71e0ab27-245f-48dc-9af5-74f922b71d16
let
    r_array = [Tuple(df1_3_large_clean_canon[!, :x1])..., Tuple(df1_3_medium_clean_canon[!, :x1])..., Tuple(df1_3_small_clean_canon[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_canon[!, :x1_1])..., Tuple(df1_3_medium_clean_canon[!, :x1_1])..., Tuple(df1_3_small_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_2_canon
    model_i_r_2_canon = lm(@formula(Y ~ X), data)
    global r_i_r_2_canon
    r_i_r_2_canon = GLM.r2(model_i_r_2_canon)
    global rms_values_i_r_2_canon
    rms_values_i_r_2_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_2_canon))
    ]
end

# ╔═╡ 1a63c8c3-5489-471b-bcd4-4d0c888f12e3
begin
    X_i_r_2_canon = DataFrame(X=collect(1:1000))
    pred_i_r_2_canon = GLM.predict(model_i_r_2_canon, X_i_r_2_canon)
end

# ╔═╡ 8b9cc340-3232-44bc-a2d9-f90279cf829b
co_i_r_2_canon = coef(model_i_r_2_canon)

# ╔═╡ d5a01db2-c0d5-4488-b7f7-dae3ce95e534
md"""
### GE
"""

# ╔═╡ 82b31f78-bdad-4ac0-ac17-5de1039e24c6
let
    r_array = [Tuple(df1_3_large_clean_geee[!, :x1])..., Tuple(df1_3_medium_clean_geee[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_geee[!, :x1_1])..., Tuple(df1_3_medium_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_2_geee
    model_i_r_2_geee = lm(@formula(Y ~ X), data)
    global r_i_r_2_geee
    r_i_r_2_geee = GLM.r2(model_i_r_2_geee)
    global rms_values_i_r_2_geee
    rms_values_i_r_2_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_2_geee))
    ]
end

# ╔═╡ afe37b2a-6f33-4813-ae29-2dd6fd5f72d3
begin
    X_i_r_2_geee = DataFrame(X=collect(1:1000))
    pred_i_r_2_geee = GLM.predict(model_i_r_2_geee, X_i_r_2_geee)
end

# ╔═╡ d31758fb-5df5-4526-9c2d-aca4c9db4048
co_i_r_2_geee = coef(model_i_r_2_geee)

# ╔═╡ ca3e0281-a373-4f64-a38a-3421c94668ac
md"""
### Philips
"""

# ╔═╡ 13c80808-6270-4165-a9d7-289745282c0f
let
    r_array = [Tuple(df1_3_large_clean_philips[!, :x1])..., Tuple(df1_3_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_philips[!, :x1_1])..., Tuple(df1_3_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_2_philips
    model_i_r_2_philips = lm(@formula(Y ~ X), data)
    global r_i_r_2_philips
    r_i_r_2_philips = GLM.r2(model_i_r_2_philips)
    global rms_values_i_r_2_philips
    rms_values_i_r_2_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_2_philips))
    ]
end

# ╔═╡ c40ab351-ff02-497e-9eda-4d6cb4a58636
begin
    X_i_r_2_philips = DataFrame(X=collect(1:1000))
    pred_i_r_2_philips = GLM.predict(model_i_r_2_philips, X_i_r_2_philips)
end

# ╔═╡ 3648622c-68a9-4de9-8ad2-309b9c08afd7
co_i_r_2_philips = coef(model_i_r_2_philips)

# ╔═╡ 06bd09a5-3d97-488a-8186-d1bd920e9acf
md"""
### Siemens
"""

# ╔═╡ 2b7fc654-4d1b-465c-a2a8-699f6d2ae373
let
    r_array = [Tuple(df1_3_large_clean_siemens[!, :x1])..., Tuple(df1_3_medium_clean_siemens[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_siemens[!, :x1_1])..., Tuple(df1_3_medium_clean_siemens[!, :x1_1])...,]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_2_siemens
    model_i_r_2_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_2_siemens
    r_i_r_2_siemens = GLM.r2(model_i_r_2_siemens)
    global rms_values_i_r_2_siemens
    rms_values_i_r_2_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_2_siemens))
    ]
end

# ╔═╡ 53165e35-6836-4129-b47b-92f9a4caa422
begin
    X_i_r_2_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_2_siemens = GLM.predict(model_i_r_2_siemens, X_i_r_2_siemens)
end

# ╔═╡ 560ee5a6-9782-466d-97dd-64196b6fd3d6
co_i_r_2_siemens = coef(model_i_r_2_siemens)

# ╔═╡ b6373747-4f01-4a85-a816-8b9a417f27aa
md"""
### Visualize
"""

# ╔═╡ a8877540-3144-4068-83c1-f183eb9c5c6e
function reprod_i_1_3()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df1_3_large_clean_canon[!, :x1], df1_3_large_clean_canon[!, :x1_1])
    scatter!(ax1, df1_3_medium_clean_canon[!, :x1], df1_3_medium_clean_canon[!, :x1_1])
    scatter!(ax1, df1_3_small_clean_canon[!, :x1], df1_3_small_clean_canon[!, :x1_1], color=:red)
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_2_canon, linestyle=:dashdot)
	create_textbox(f[1, 1], co_i_r_2_canon, r_i_r_2_canon, rms_values_i_r_2_canon)

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_3_large_clean_geee[!, :x1], df1_3_large_clean_geee[!, :x1_1])
    scatter!(ax2, df1_3_medium_clean_geee[!, :x1], df1_3_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df1_3_small_clean_geee[!, :x1], df1_3_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_2_geee, linestyle=:dashdot)
	create_textbox(f[2, 1], co_i_r_2_geee, r_i_r_2_geee, rms_values_i_r_2_geee)

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_3_large_clean_philips[!, :x1], df1_3_large_clean_philips[!, :x1_1])
    scatter!(ax3, df1_3_medium_clean_philips[!, :x1], df1_3_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df1_3_small_clean_philips[!, :x1], df1_3_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_2_philips, linestyle=:dashdot)
	create_textbox(f[1, 2], co_i_r_2_philips, r_i_r_2_philips, rms_values_i_r_2_philips)

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_3_large_clean_siemens[!, :x1], df1_3_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_3_medium_clean_siemens[!, :x1], df1_3_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_2_siemens, linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], co_i_r_2_siemens, r_i_r_2_siemens, rms_values_i_r_2_siemens)

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (1 vs 3)",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ 9d677e1c-5754-44ce-beb1-1f73e9c0fffe
with_theme(medphys_theme) do
    reprod_i_1_3()
end

# ╔═╡ 9a4dd89d-5fd4-4c25-98d5-2bb1d3d85f36
md"""
## 1 vs 4
"""

# ╔═╡ f9374f75-36a5-42d6-9607-7d5622b46f4e
begin
    df1_4_large_clean = innerjoin(df1_large, df4_large; on=[:vendor, :inserts], makeunique=true)
    df1_4_medium_clean = innerjoin(df1_medium, df4_medium; on=[:vendor, :inserts], makeunique=true)
    df1_4_small_clean = innerjoin(df1_small, df4_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 31b61ad7-ac39-4610-b902-2a2943ebace2
md"""
### Canon
"""

# ╔═╡ 52f55f1b-f904-4c83-9992-b5ccb8099ec7
df1_4_small_clean

# ╔═╡ 61fccccd-dfdf-4705-9356-767d0ab2d72b
begin
    df1_4_large_clean_canon, df1_4_large_clean_geee, df1_4_large_clean_philips, df1_4_large_clean_siemens = groupby(df1_4_large_clean, :vendor)
    df1_4_medium_clean_canon, df1_4_medium_clean_geee, df1_4_medium_clean_philips, df1_4_medium_clean_siemens = groupby(df1_4_medium_clean, :vendor)
    df1_4_small_clean_canon = df1_4_small_clean
end;

# ╔═╡ 80b1e347-0550-4002-846b-424c96fbc97e
let
    r_array = [Tuple(df1_4_large_clean_canon[!, :x1])..., Tuple(df1_4_medium_clean_canon[!, :x1])..., Tuple(df1_4_small_clean_canon[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_canon[!, :x1_1])..., Tuple(df1_4_medium_clean_canon[!, :x1_1])..., Tuple(df1_4_small_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_3_canon
    model_i_r_3_canon = lm(@formula(Y ~ X), data)
    global r_i_r_3_canon
    r_i_r_3_canon = GLM.r2(model_i_r_3_canon)
    global rms_values_i_r_3_canon
    rms_values_i_r_3_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_3_canon))
    ]
end

# ╔═╡ 0b649084-457c-4339-bb19-bf3a7200262b
begin
    X_i_r_3_canon = DataFrame(X=collect(1:1000))
    pred_i_r_3_canon = GLM.predict(model_i_r_3_canon, X_i_r_3_canon)
end

# ╔═╡ 12494832-c0fa-45c1-bbd0-05b82543e138
co_i_r_3_canon = coef(model_i_r_3_canon)

# ╔═╡ 16c82c3f-7c99-496d-861e-6da382d1a75f
md"""
### GE
"""

# ╔═╡ 8dc2ac38-8a45-4799-b5d0-a849aa1f7c26
let
    r_array = [Tuple(df1_4_large_clean_geee[!, :x1])..., Tuple(df1_4_medium_clean_geee[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_geee[!, :x1_1])..., Tuple(df1_4_medium_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_3_geee
    model_i_r_3_geee = lm(@formula(Y ~ X), data)
    global r_i_r_3_geee
    r_i_r_3_geee = GLM.r2(model_i_r_3_geee)
    global rms_values_i_r_3_geee
    rms_values_i_r_3_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_3_geee))
    ]
end

# ╔═╡ 5b360dc3-213b-4c01-ab89-b4ddc888eb87
begin
    X_i_r_3_geee = DataFrame(X=collect(1:1000))
    pred_i_r_3_geee = GLM.predict(model_i_r_3_geee, X_i_r_3_geee)
end

# ╔═╡ 6529d5ad-1539-42b2-aa9a-abad76ddf96d
co_i_r_3_geee = coef(model_i_r_3_geee)

# ╔═╡ 356398c4-66e0-4b49-b1ae-2b4815d5e54a
md"""
### Philips
"""

# ╔═╡ fc6eb521-6ddc-4be2-a2d1-685db850e6ee
let
    r_array = [Tuple(df1_4_large_clean_philips[!, :x1])..., Tuple(df1_4_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_philips[!, :x1_1])..., Tuple(df1_4_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_3_philips
    model_i_r_3_philips = lm(@formula(Y ~ X), data)
    global r_i_r_3_philips
    r_i_r_3_philips = GLM.r2(model_i_r_3_philips)
    global rms_values_i_r_3_philips
    rms_values_i_r_3_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_3_philips))
    ]
end

# ╔═╡ 37620e01-de7f-4479-9014-38c0ce5407f6
begin
    X_i_r_3_philips = DataFrame(X=collect(1:1000))
    pred_i_r_3_philips = GLM.predict(model_i_r_3_philips, X_i_r_3_philips)
end

# ╔═╡ 6db477d4-197c-41bc-bf9d-78a1b83c4749
co_i_r_3_philips = coef(model_i_r_3_philips)

# ╔═╡ 6caac89d-80ed-4ecf-8d99-37181f500780
md"""
### Siemens
"""

# ╔═╡ 1384e914-d7ac-4459-903d-dc8e8f360711
let
    r_array = [Tuple(df1_4_large_clean_siemens[!, :x1])..., Tuple(df1_4_medium_clean_siemens[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_siemens[!, :x1_1])..., Tuple(df1_4_medium_clean_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_3_siemens
    model_i_r_3_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_3_siemens
    r_i_r_3_siemens = GLM.r2(model_i_r_3_siemens)
    global rms_values_i_r_3_siemens
    rms_values_i_r_3_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_3_siemens))
    ]
end

# ╔═╡ 39b386be-d0b5-4b9f-8015-450988a5255f
begin
    X_i_r_3_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_3_siemens = GLM.predict(model_i_r_3_siemens, X_i_r_3_siemens)
end

# ╔═╡ b1963053-588f-4137-9ade-a308eedd496a
co_i_r_3_siemens = coef(model_i_r_3_siemens)

# ╔═╡ 825d85a5-3898-4fb2-8fe3-30c8ba49dd0a
md"""
### Visualize
"""

# ╔═╡ f1e9e113-9531-4fb4-95aa-5dfaa146c47e
function reprod_i_1_4()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df1_4_large_clean_canon[!, :x1], df1_4_large_clean_canon[!, :x1_1])
    scatter!(ax1, df1_4_medium_clean_canon[!, :x1], df1_4_medium_clean_canon[!, :x1_1])
    scatter!(ax1, df1_4_small_clean_canon[!, :x1], df1_4_small_clean_canon[!, :x1_1], color=:red)
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_3_canon, linestyle=:dashdot)
    if co_i_r_3_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_3_canon[2]; digits=3))x+$(trunc(co_i_r_3_canon[1]; digits=3)) \nr = $(trunc(r_i_r_3_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_3_canon[2]; digits=3))x$(trunc(co_i_r_3_canon[1]; digits=3)) \nr = $(trunc(r_i_r_3_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_4_large_clean_geee[!, :x1], df1_4_large_clean_geee[!, :x1_1])
    scatter!(ax2, df1_4_medium_clean_geee[!, :x1], df1_4_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df1_4_small_clean_geee[!, :x1], df1_4_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_3_geee, linestyle=:dashdot)
    if co_i_r_3_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_3_geee[2]; digits=3))x+$(trunc(co_i_r_3_geee[1]; digits=3)) \nr = $(trunc(r_i_r_3_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_3_geee[2]; digits=3))x$(trunc(co_i_r_3_geee[1]; digits=3)) \nr = $(trunc(r_i_r_3_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_4_large_clean_philips[!, :x1], df1_4_large_clean_philips[!, :x1_1])
    scatter!(ax3, df1_4_medium_clean_philips[!, :x1], df1_4_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df1_4_small_clean_philips[!, :x1], df1_4_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_3_philips, linestyle=:dashdot)
    if co_i_r_3_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_3_philips[2]; digits=3))x+$(trunc(co_i_r_3_philips[1]; digits=3)) \nr = $(trunc(r_i_r_3_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_3_philips[2]; digits=3))x$(trunc(co_i_r_3_philips[1]; digits=3)) \nr = $(trunc(r_i_r_3_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_4_large_clean_siemens[!, :x1], df1_4_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_4_medium_clean_siemens[!, :x1], df1_4_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_3_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_i_r_3_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_3_siemens[2]; digits=3))x+$(trunc(co_i_r_3_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_3_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_3_siemens[2]; digits=3))x$(trunc(co_i_r_3_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_3_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_3_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_3_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (1 vs 4)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ ce1d1375-fe0c-436a-930e-e6e5cc12c172
with_theme(medphys_theme) do
    reprod_i_1_4()
end

# ╔═╡ 25c01496-5893-40b3-9c55-a7473180ff37
md"""
## 1 vs 5
"""

# ╔═╡ 9b36dd13-1cc6-480e-a8a3-1e36413e900a
begin
    df1_5_large_clean = innerjoin(df1_large, df5_large; on=[:vendor, :inserts], makeunique=true)
    df1_5_medium_clean = innerjoin(df1_medium, df5_medium; on=[:vendor, :inserts], makeunique=true)
    df1_5_small_clean = innerjoin(df1_small, df5_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ ca50b8dd-fe0b-408a-9989-948d07a09ffc
md"""
### Canon
"""

# ╔═╡ 430bfe52-0ebd-4de2-9a02-cc09fb4edcc0
df1_5_small_clean

# ╔═╡ 881771b5-373e-4ff3-9b01-405ae56677e5
begin
    df1_5_large_clean_canon, df1_5_large_clean_geee, df1_5_large_clean_philips, df1_5_large_clean_siemens = groupby(df1_5_large_clean, :vendor)
    df1_5_medium_clean_canon, df1_5_medium_clean_geee, df1_5_medium_clean_philips, df1_5_medium_clean_siemens = groupby(df1_5_medium_clean, :vendor)
    df1_5_small_clean_canon = df1_5_small_clean
end;

# ╔═╡ 31f3cbfb-e9c2-46dc-ad0c-a637783b421c
let
    r_array = [Tuple(df1_5_large_clean_canon[!, :x1])..., Tuple(df1_5_medium_clean_canon[!, :x1])..., Tuple(df1_5_small_clean_canon[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_canon[!, :x1_1])..., Tuple(df1_5_medium_clean_canon[!, :x1_1])..., Tuple(df1_5_small_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_4_canon
    model_i_r_4_canon = lm(@formula(Y ~ X), data)
    global r_i_r_4_canon
    r_i_r_4_canon = GLM.r2(model_i_r_4_canon)
    global rms_values_i_r_4_canon
    rms_values_i_r_4_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_4_canon))
    ]
end

# ╔═╡ 480296d4-890f-4484-b0e7-8828dd51df26
begin
    X_i_r_4_canon = DataFrame(X=collect(1:1000))
    pred_i_r_4_canon = GLM.predict(model_i_r_4_canon, X_i_r_4_canon)
end

# ╔═╡ 7bea8ddf-0cb5-4bd7-86c1-24b6e91befac
co_i_r_4_canon = coef(model_i_r_4_canon)

# ╔═╡ b5e12cdf-3ecb-4d9b-92cd-b76f9f63a5ee
md"""
### GE
"""

# ╔═╡ bc10bdfc-fa49-4b36-9020-bacf6ddd7940
let
    r_array = [Tuple(df1_5_large_clean_geee[!, :x1])..., Tuple(df1_5_medium_clean_geee[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_geee[!, :x1_1])..., Tuple(df1_5_medium_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_4_geee
    model_i_r_4_geee = lm(@formula(Y ~ X), data)
    global r_i_r_4_geee
    r_i_r_4_geee = GLM.r2(model_i_r_4_geee)
    global rms_values_i_r_4_geee
    rms_values_i_r_4_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_4_geee))
    ]
end

# ╔═╡ daa19304-8d3c-4e6b-9d60-0b5f064e5643
begin
    X_i_r_4_geee = DataFrame(X=collect(1:1000))
    pred_i_r_4_geee = GLM.predict(model_i_r_4_geee, X_i_r_4_geee)
end

# ╔═╡ 39dc0b0d-834a-43ea-a3b5-45d603bc4413
co_i_r_4_geee = coef(model_i_r_4_geee)

# ╔═╡ 025f813a-1a26-475a-994c-e7c7afcd9d15
md"""
### Philips
"""

# ╔═╡ acf5f08c-edb9-4da4-8c08-672f780cc41b
let
    r_array = [Tuple(df1_5_large_clean_philips[!, :x1])..., Tuple(df1_5_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_philips[!, :x1_1])..., Tuple(df1_5_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_4_philips
    model_i_r_4_philips = lm(@formula(Y ~ X), data)
    global r_i_r_4_philips
    r_i_r_4_philips = GLM.r2(model_i_r_4_philips)
    global rms_values_i_r_4_philips
    rms_values_i_r_4_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_4_philips))
    ]
end

# ╔═╡ bceb632d-c1b1-4c19-90fd-ea1d97c1db9c
begin
    X_i_r_4_philips = DataFrame(X=collect(1:1000))
    pred_i_r_4_philips = GLM.predict(model_i_r_4_philips, X_i_r_4_philips)
end

# ╔═╡ d85f6f5c-3160-492d-9aaa-80d19d7826fa
co_i_r_4_philips = coef(model_i_r_4_philips)

# ╔═╡ 917e1ee3-2f55-466e-aad6-16eb0241a2ff
md"""
### Siemens
"""

# ╔═╡ 4957ff4a-4457-4bd7-9b54-e86561165788
let
    r_array = [Tuple(df1_5_large_clean_siemens[!, :x1])..., Tuple(df1_5_medium_clean_siemens[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_siemens[!, :x1_1])..., Tuple(df1_5_medium_clean_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_4_siemens
    model_i_r_4_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_4_siemens
    r_i_r_4_siemens = GLM.r2(model_i_r_4_siemens)
    global rms_values_i_r_4_siemens
    rms_values_i_r_4_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_4_siemens))
    ]
end

# ╔═╡ dc57ba01-1731-4136-964a-4f4da385a967
begin
    X_i_r_4_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_4_siemens = GLM.predict(model_i_r_4_siemens, X_i_r_4_siemens)
end

# ╔═╡ c2199225-23a1-48da-b58c-fe03994dcd38
co_i_r_4_siemens = coef(model_i_r_4_siemens)

# ╔═╡ 86adbd4f-b55d-4056-a9f7-aafaebca9b3c
md"""
### Visualize
"""

# ╔═╡ 60602785-9402-49fc-a898-f8538714feee
function reprod_i_1_5()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df1_5_large_clean_canon[!, :x1], df1_5_large_clean_canon[!, :x1_1])
    scatter!(ax1, df1_5_medium_clean_canon[!, :x1], df1_5_medium_clean_canon[!, :x1_1])
    scatter!(ax1, df1_5_small_clean_canon[!, :x1], df1_5_small_clean_canon[!, :x1_1], color=:red)
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_4_canon, linestyle=:dashdot)
    if co_i_r_4_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_4_canon[2]; digits=3))x+$(trunc(co_i_r_4_canon[1]; digits=3)) \nr = $(trunc(r_i_r_4_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_4_canon[2]; digits=3))x$(trunc(co_i_r_4_canon[1]; digits=3)) \nr = $(trunc(r_i_r_4_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_5_large_clean_geee[!, :x1], df1_5_large_clean_geee[!, :x1_1])
    scatter!(ax2, df1_5_medium_clean_geee[!, :x1], df1_5_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df1_5_small_clean_geee[!, :x1], df1_5_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_4_geee, linestyle=:dashdot)
    if co_i_r_4_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_4_geee[2]; digits=3))x+$(trunc(co_i_r_4_geee[1]; digits=3)) \nr = $(trunc(r_i_r_4_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_4_geee[2]; digits=3))x$(trunc(co_i_r_4_geee[1]; digits=3)) \nr = $(trunc(r_i_r_4_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_5_large_clean_philips[!, :x1], df1_5_large_clean_philips[!, :x1_1])
    scatter!(ax3, df1_5_medium_clean_philips[!, :x1], df1_5_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df1_5_small_clean_philips[!, :x1], df1_5_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_4_philips, linestyle=:dashdot)
    if co_i_r_4_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_4_philips[2]; digits=3))x+$(trunc(co_i_r_4_philips[1]; digits=3)) \nr = $(trunc(r_i_r_4_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_4_philips[2]; digits=3))x$(trunc(co_i_r_4_philips[1]; digits=3)) \nr = $(trunc(r_i_r_4_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_5_large_clean_siemens[!, :x1], df1_5_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_5_medium_clean_siemens[!, :x1], df1_5_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_4_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_i_r_4_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_4_siemens[2]; digits=3))x+$(trunc(co_i_r_4_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_4_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_4_siemens[2]; digits=3))x$(trunc(co_i_r_4_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_4_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_4_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_4_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (1 vs 5)",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ c0c47a1b-7f84-48bd-b5eb-b80302d81514
with_theme(medphys_theme) do
    reprod_i_1_5()
end

# ╔═╡ a7363057-b7bd-410a-95fa-2b70c3a52f16
md"""
# Large Phantom - Agatston
"""

# ╔═╡ a1f8e8c6-211d-42a3-82aa-32a8560752f3
md"""
## 1 vs 2
"""

# ╔═╡ c33449f5-faba-4e86-942d-460773849ecc
begin
    df1_2_large_clean_a = innerjoin(df1_large_a, df2_large_a; on=[:vendor, :inserts], makeunique=true)
    df1_2_medium_clean_a = innerjoin(df1_medium_a, df2_medium_a; on=[:vendor, :inserts], makeunique=true)
    df1_2_small_clean_a = innerjoin(df1_small_a, df2_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ d78573c3-2ea8-49a8-81d4-077343d3a274
md"""
### Canon
"""

# ╔═╡ e60e3a6d-8815-433c-af18-1b8938e0afc6
begin
    df1_2_large_clean_a_canon, df1_2_large_clean_a_geee, df1_2_large_clean_a_philips, df1_2_large_clean_a_siemens = groupby(df1_2_large_clean_a, :vendor)
    df1_2_medium_clean_a_canon, df1_2_medium_clean_a_geee, df1_2_medium_clean_a_philips, df1_2_medium_clean_a_siemens = groupby(df1_2_medium_clean_a, :vendor)
    df1_2_small_clean_a_canon = df1_2_small_clean_a
end;

# ╔═╡ 250cd449-8bba-4854-8197-8f1d3268ee41
let
    r_array = [Tuple(df1_2_large_clean_a_canon[!, :x1])..., Tuple(df1_2_medium_clean_a_canon[!, :x1])..., Tuple(df1_2_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df1_2_large_clean_a_canon[!, :x1_1])..., Tuple(df1_2_medium_clean_a_canon[!, :x1_1])..., Tuple(df1_2_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_1_canon
    model_a_r_1_canon = lm(@formula(Y ~ X), data)
    global r_a_r_1_canon
    r_a_r_1_canon = GLM.r2(model_a_r_1_canon)
    global rms_values_a_r_1_canon
    rms_values_a_r_1_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_1_canon))
    ]
end

# ╔═╡ a9440691-ac39-4c64-b65d-83d746ebc0da
begin
    X_a_r_1_canon = DataFrame(X=collect(1:1000))
    pred_a_r_1_canon = GLM.predict(model_a_r_1_canon, X_a_r_1_canon)
end

# ╔═╡ 2bdb6b1c-e4a8-49f2-bb97-022f553373d7
co_a_r_1_canon = coef(model_a_r_1_canon)

# ╔═╡ 7adf35a0-f75b-4075-8472-425b400998c7
md"""
### GE
"""

# ╔═╡ 8a8bbef8-6806-47a2-8d8e-5e634fb45e0d
let
    r_array = [Tuple(df1_2_large_clean_a_geee[!, :x1])..., Tuple(df1_2_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df1_2_large_clean_a_geee[!, :x1_1])..., Tuple(df1_2_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_1_geee
    model_a_r_1_geee = lm(@formula(Y ~ X), data)
    global r_a_r_1_geee
    r_a_r_1_geee = GLM.r2(model_a_r_1_geee)
    global rms_values_a_r_1_geee
    rms_values_a_r_1_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_1_geee))
    ]
end

# ╔═╡ b568c1e3-e332-4713-8295-15e25b0fee83
begin
    X_a_r_1_geee = DataFrame(X=collect(1:1000))
    pred_a_r_1_geee = GLM.predict(model_a_r_1_geee, X_a_r_1_geee)
end

# ╔═╡ d8a1f7eb-6d0d-4e18-9674-d9aa169f09e2
co_a_r_1_geee = coef(model_a_r_1_geee)

# ╔═╡ 7a3ca058-d58f-4a56-901e-7b6f7547836a
md"""
### Philips
"""

# ╔═╡ 1f3c18b9-f001-4657-ab87-1fef4782f752
let
    r_array = [Tuple(df1_2_large_clean_a_philips[!, :x1])..., Tuple(df1_2_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df1_2_large_clean_a_philips[!, :x1_1])..., Tuple(df1_2_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_1_philips
    model_a_r_1_philips = lm(@formula(Y ~ X), data)
    global r_a_r_1_philips
    r_a_r_1_philips = GLM.r2(model_a_r_1_philips)
    global rms_values_a_r_1_philips
    rms_values_a_r_1_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_1_philips))
    ]
end

# ╔═╡ 60b064f0-0341-4c9d-a0e5-50900b5e502b
begin
    X_a_r_1_philips = DataFrame(X=collect(1:1000))
    pred_a_r_1_philips = GLM.predict(model_a_r_1_philips, X_a_r_1_philips)
end

# ╔═╡ e209caf6-f439-4ada-8b60-51881c572436
co_a_r_1_philips = coef(model_a_r_1_philips)

# ╔═╡ c50609b3-ac5e-43ac-93bc-cd058eb05d8c
md"""
### Siemens
"""

# ╔═╡ d4c78b5b-7501-4513-81a4-06037e190292
let
    r_array = [Tuple(df1_2_large_clean_a_siemens[!, :x1])..., Tuple(df1_2_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df1_2_large_clean_a_siemens[!, :x1_1])..., Tuple(df1_2_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_1_siemens
    model_a_r_1_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_1_siemens
    r_a_r_1_siemens = GLM.r2(model_a_r_1_siemens)
    global rms_values_a_r_1_siemens
    rms_values_a_r_1_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_1_siemens))
    ]
end

# ╔═╡ 6abaedc2-e1e9-4e4b-8c96-8f56d15a46ae
begin
    X_a_r_1_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_1_siemens = GLM.predict(model_a_r_1_siemens, X_a_r_1_siemens)
end

# ╔═╡ f972ea2c-2bdc-4caa-80a1-8aff4456b3c1
co_a_r_1_siemens = coef(model_a_r_1_siemens)

# ╔═╡ 0a42ad4a-5ad9-4bef-a5b9-3153b773957e
md"""
### Visualize
"""

# ╔═╡ 95f7ddba-c5e7-4725-aff4-9b770a2b0fc0
function reprod_a_1_2()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df1_2_large_clean_a_canon[!, :x1], df1_2_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df1_2_medium_clean_a_canon[!, :x1], df1_2_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df1_2_small_clean_a_canon[!, :x1], df1_2_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_1_canon, linestyle=:dashdot)
	create_textbox(f[1, 1], co_a_r_1_canon, r_a_r_1_canon, rms_values_a_r_1_canon)

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "Scanner 1"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_2_large_clean_a_geee[!, :x1], df1_2_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df1_2_medium_clean_a_geee[!, :x1], df1_2_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_1_geee, linestyle=:dashdot)
	create_textbox(f[2, 1], co_a_r_1_geee, r_a_r_1_geee, rms_values_a_r_1_geee)

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "Scanner 2"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_2_large_clean_a_philips[!, :x1], df1_2_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df1_2_medium_clean_a_philips[!, :x1], df1_2_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_1_philips, linestyle=:dashdot)
	create_textbox(f[1, 2], co_a_r_1_philips, r_a_r_1_philips, rms_values_a_r_1_philips)

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "Scanner 3"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_2_large_clean_a_siemens[!, :x1], df1_2_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_2_medium_clean_a_siemens[!, :x1], df1_2_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_1_siemens, linestyle=:dashdot, label="Fitted Line")
	create_textbox(f[2, 2], co_a_r_1_siemens, r_a_r_1_siemens, rms_values_a_r_1_siemens)

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "Scanner 4"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Scoring (1 vs 2)",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

	save(plotsdir("repro_a_1_2.png"), f)
    f
end

# ╔═╡ 1c68825f-118d-4c38-bbd3-1f6b8688e8a2
with_theme(medphys_theme) do
    reprod_a_1_2()
end

# ╔═╡ e830479f-7f74-459a-a531-5fc7dc6122a1
md"""
## 1 vs 3
"""

# ╔═╡ ff76ef96-57e4-49f2-9776-e2aee3f86c91
begin
    df1_3_large_clean_a = innerjoin(df1_large_a, df3_large_a; on=[:vendor, :inserts], makeunique=true)
    df1_3_medium_clean_a = innerjoin(df1_medium_a, df3_medium_a; on=[:vendor, :inserts], makeunique=true)
    df1_3_small_clean_a = innerjoin(df1_small_a, df3_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ ef3eb4ca-3ef8-4700-9efd-6e13e6069faf
md"""
### Canon
"""

# ╔═╡ 091fd9b6-ce40-437b-aa80-cf82fc3ed676
begin
    df1_3_large_clean_a_canon, df1_3_large_clean_a_geee, df1_3_large_clean_a_philips, df1_3_large_clean_a_siemens = groupby(df1_3_large_clean_a, :vendor)
    df1_3_medium_clean_a_canon, df1_3_medium_clean_a_geee, df1_3_medium_clean_a_philips, df1_3_medium_clean_a_siemens = groupby(df1_3_medium_clean_a, :vendor)
    df1_3_small_clean_a_canon = df1_3_small_clean_a
end;

# ╔═╡ ce214147-4f83-48d0-b69e-69142214b3d2
let
    r_array = [Tuple(df1_3_large_clean_a_canon[!, :x1])..., Tuple(df1_3_medium_clean_a_canon[!, :x1])..., Tuple(df1_3_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_a_canon[!, :x1_1])..., Tuple(df1_3_medium_clean_a_canon[!, :x1_1])..., Tuple(df1_3_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_2_canon
    model_a_r_2_canon = lm(@formula(Y ~ X), data)
    global r_a_r_2_canon
    r_a_r_2_canon = GLM.r2(model_a_r_2_canon)
    global rms_values_a_r_2_canon
    rms_values_a_r_2_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_2_canon))
    ]
end

# ╔═╡ 186994e9-6624-48e3-b383-f67e659ebfcc
begin
    X_a_r_2_canon = DataFrame(X=collect(1:1000))
    pred_a_r_2_canon = GLM.predict(model_a_r_2_canon, X_a_r_2_canon)
end

# ╔═╡ 9fbb0570-3d7e-4dc8-9919-0f70349980ef
co_a_r_2_canon = coef(model_a_r_2_canon)

# ╔═╡ 7f2edfc7-9c6c-4edd-bd0a-28b947ff9c91
md"""
### GE
"""

# ╔═╡ 73c43529-d5c3-4725-a7e5-5c4dea85ebe5
let
    r_array = [Tuple(df1_3_large_clean_a_geee[!, :x1])..., Tuple(df1_3_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_a_geee[!, :x1_1])..., Tuple(df1_3_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_2_geee
    model_a_r_2_geee = lm(@formula(Y ~ X), data)
    global r_a_r_2_geee
    r_a_r_2_geee = GLM.r2(model_a_r_2_geee)
    global rms_values_a_r_2_geee
    rms_values_a_r_2_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_2_geee))
    ]
end

# ╔═╡ 7a82f997-dbf3-40d9-b508-eab38c96d0c7
begin
    X_a_r_2_geee = DataFrame(X=collect(1:1000))
    pred_a_r_2_geee = GLM.predict(model_a_r_2_geee, X_a_r_2_geee)
end

# ╔═╡ 44519f33-ce8c-4ced-aba5-b7f8f52f3681
co_a_r_2_geee = coef(model_a_r_2_geee)

# ╔═╡ 2d3eec4e-7cbb-4a0e-ab16-85597af04b02
md"""
### Philips
"""

# ╔═╡ b62a9ebe-7235-4235-bef7-52f259820f75
let
    r_array = [Tuple(df1_3_large_clean_a_philips[!, :x1])..., Tuple(df1_3_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_a_philips[!, :x1_1])..., Tuple(df1_3_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_2_philips
    model_a_r_2_philips = lm(@formula(Y ~ X), data)
    global r_a_r_2_philips
    r_a_r_2_philips = GLM.r2(model_a_r_2_philips)
    global rms_values_a_r_2_philips
    rms_values_a_r_2_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_2_philips))
    ]
end

# ╔═╡ 7911378c-f43a-4b1b-a93c-e805cea9fcc6
begin
    X_a_r_2_philips = DataFrame(X=collect(1:1000))
    pred_a_r_2_philips = GLM.predict(model_a_r_2_philips, X_a_r_2_philips)
end

# ╔═╡ 735fe265-976c-45b3-b132-b4a9dc8c5e7c
co_a_r_2_philips = coef(model_a_r_2_philips)

# ╔═╡ 15db30dd-d042-4433-9baf-1d262257b09b
md"""
### Siemens
"""

# ╔═╡ 881ece53-f462-49e3-a597-21f1154e73ac
let
    r_array = [Tuple(df1_3_large_clean_a_siemens[!, :x1])..., Tuple(df1_3_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df1_3_large_clean_a_siemens[!, :x1_1])..., Tuple(df1_3_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_2_siemens
    model_a_r_2_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_2_siemens
    r_a_r_2_siemens = GLM.r2(model_a_r_2_siemens)
    global rms_values_a_r_2_siemens
    rms_values_a_r_2_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_2_siemens))
    ]
end

# ╔═╡ a9144e10-227f-44df-9a73-1bff8597e6ff
begin
    X_a_r_2_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_2_siemens = GLM.predict(model_a_r_2_siemens, X_a_r_2_siemens)
end

# ╔═╡ 709e0630-1d70-4e8b-a550-b760dad49788
co_a_r_2_siemens = coef(model_a_r_2_siemens)

# ╔═╡ ae1562ca-a019-4de8-9d17-dbee895a610d
md"""
### Visualize
"""

# ╔═╡ a7eae4c1-8fd9-42c3-8653-7572060bae2f
function reprod_a_1_3()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df1_3_large_clean_a_canon[!, :x1], df1_3_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df1_3_medium_clean_a_canon[!, :x1], df1_3_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df1_3_small_clean_a_canon[!, :x1], df1_3_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_2_canon, linestyle=:dashdot)
    if co_a_r_2_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_2_canon[2]; digits=3))x+$(trunc(co_a_r_2_canon[1]; digits=3)) \nr = $(trunc(r_a_r_2_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_2_canon[2]; digits=3))x$(trunc(co_a_r_2_canon[1]; digits=3)) \nr = $(trunc(r_a_r_2_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_3_large_clean_a_geee[!, :x1], df1_3_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df1_3_medium_clean_a_geee[!, :x1], df1_3_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_2_geee, linestyle=:dashdot)
    if co_a_r_2_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_2_geee[2]; digits=3))x+$(trunc(co_a_r_2_geee[1]; digits=3)) \nr = $(trunc(r_a_r_2_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_2_geee[2]; digits=3))x$(trunc(co_a_r_2_geee[1]; digits=3)) \nr = $(trunc(r_a_r_2_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_3_large_clean_a_philips[!, :x1], df1_3_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df1_3_medium_clean_a_philips[!, :x1], df1_3_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_2_philips, linestyle=:dashdot)
    if co_a_r_2_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_2_philips[2]; digits=3))x+$(trunc(co_a_r_2_philips[1]; digits=3)) \nr = $(trunc(r_a_r_2_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_2_philips[2]; digits=3))x$(trunc(co_a_r_2_philips[1]; digits=3)) \nr = $(trunc(r_a_r_2_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_3_large_clean_a_siemens[!, :x1], df1_3_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_3_medium_clean_a_siemens[!, :x1], df1_3_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_2_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_a_r_2_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_2_siemens[2]; digits=3))x+$(trunc(co_a_r_2_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_2_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_2_siemens[2]; digits=3))x$(trunc(co_a_r_2_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_2_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_2_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_2_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Calcium Scoring (1 vs 3)",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ a2119a75-fa6b-4b36-b70a-942856a1fe6a
with_theme(medphys_theme) do
    reprod_a_1_3()
end

# ╔═╡ 1f07992d-747c-430c-a53c-e44df17f6934
md"""
## 1 vs 4
"""

# ╔═╡ 639c2606-c2b7-4130-a005-ffa314f5096d
begin
    df1_4_large_clean_a = innerjoin(df1_large_a, df4_large_a; on=[:vendor, :inserts], makeunique=true)
    df1_4_medium_clean_a = innerjoin(df1_medium_a, df4_medium_a; on=[:vendor, :inserts], makeunique=true)
    df1_4_small_clean_a = innerjoin(df1_small_a, df4_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 8ca49384-3739-4b0c-beba-10697655a9a7
md"""
### Canon
"""

# ╔═╡ 0b30bdbd-a520-43d2-9c24-121f8912809d
begin
    df1_4_large_clean_a_canon, df1_4_large_clean_a_geee, df1_4_large_clean_a_philips, df1_4_large_clean_a_siemens = groupby(df1_4_large_clean_a, :vendor)
    df1_4_medium_clean_a_canon, df1_4_medium_clean_a_geee, df1_4_medium_clean_a_philips, df1_4_medium_clean_a_siemens = groupby(df1_4_medium_clean_a, :vendor)
    df1_4_small_clean_a_canon = df1_4_small_clean_a
end;

# ╔═╡ d3e16b7d-a15f-471e-9080-e8039bd2760f
let
    r_array = [Tuple(df1_4_large_clean_a_canon[!, :x1])..., Tuple(df1_4_medium_clean_a_canon[!, :x1])..., Tuple(df1_4_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_a_canon[!, :x1_1])..., Tuple(df1_4_medium_clean_a_canon[!, :x1_1])..., Tuple(df1_4_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_3_canon
    model_a_r_3_canon = lm(@formula(Y ~ X), data)
    global r_a_r_3_canon
    r_a_r_3_canon = GLM.r2(model_a_r_3_canon)
    global rms_values_a_r_3_canon
    rms_values_a_r_3_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_3_canon))
    ]
end

# ╔═╡ 5a42bf7b-4949-41fa-99b0-1ee8ba304f74
begin
    X_a_r_3_canon = DataFrame(X=collect(1:1000))
    pred_a_r_3_canon = GLM.predict(model_a_r_3_canon, X_a_r_3_canon)
end

# ╔═╡ 157b5e5f-decb-4c90-9c3d-e044890da547
co_a_r_3_canon = coef(model_a_r_3_canon)

# ╔═╡ 86fbe863-7404-488d-afb1-d485d97c131f
md"""
### GE
"""

# ╔═╡ a5dc5b63-8042-4fed-a4fd-199a337cc120
let
    r_array = [Tuple(df1_4_large_clean_a_geee[!, :x1])..., Tuple(df1_4_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_a_geee[!, :x1_1])..., Tuple(df1_4_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_3_geee
    model_a_r_3_geee = lm(@formula(Y ~ X), data)
    global r_a_r_3_geee
    r_a_r_3_geee = GLM.r2(model_a_r_3_geee)
    global rms_values_a_r_3_geee
    rms_values_a_r_3_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_3_geee))
    ]
end

# ╔═╡ a8fb5fea-58ce-459a-8433-822a9f2fcc50
begin
    X_a_r_3_geee = DataFrame(X=collect(1:1000))
    pred_a_r_3_geee = GLM.predict(model_a_r_3_geee, X_a_r_3_geee)
end

# ╔═╡ 51ec9457-a8a2-4d5b-96d5-d8856af3ce48
co_a_r_3_geee = coef(model_a_r_3_geee)

# ╔═╡ efdb711e-1bcd-4566-b49b-e302e1536e2a
md"""
### Philips
"""

# ╔═╡ a071e559-ed7c-4001-8c94-a73b3e8f5df9
let
    r_array = [Tuple(df1_4_large_clean_a_philips[!, :x1])..., Tuple(df1_4_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_a_philips[!, :x1_1])..., Tuple(df1_4_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_3_philips
    model_a_r_3_philips = lm(@formula(Y ~ X), data)
    global r_a_r_3_philips
    r_a_r_3_philips = GLM.r2(model_a_r_3_philips)
    global rms_values_a_r_3_philips
    rms_values_a_r_3_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_3_philips))
    ]
end

# ╔═╡ 0ae71358-9fe6-4484-bb52-fd4df9e1f4e3
begin
    X_a_r_3_philips = DataFrame(X=collect(1:1000))
    pred_a_r_3_philips = GLM.predict(model_a_r_3_philips, X_a_r_3_philips)
end

# ╔═╡ 89fdec4f-d819-4199-b6d0-3753c16e4f47
co_a_r_3_philips = coef(model_a_r_3_philips)

# ╔═╡ 616a4888-2e30-4fb8-9bda-fa9cb6dead52
md"""
### Siemens
"""

# ╔═╡ e4d6734a-7b3d-4932-bef3-e7e7bd3b1d11
let
    r_array = [Tuple(df1_4_large_clean_a_siemens[!, :x1])..., Tuple(df1_4_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df1_4_large_clean_a_siemens[!, :x1_1])..., Tuple(df1_4_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_3_siemens
    model_a_r_3_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_3_siemens
    r_a_r_3_siemens = GLM.r2(model_a_r_3_siemens)
    global rms_values_a_r_3_siemens
    rms_values_a_r_3_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_3_siemens))
    ]
end

# ╔═╡ 8a819052-9c66-450c-b4a9-31946883c95b
begin
    X_a_r_3_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_3_siemens = GLM.predict(model_a_r_3_siemens, X_a_r_3_siemens)
end

# ╔═╡ aa54f20e-a59d-4194-8a1c-46c3d1e53f21
co_a_r_3_siemens = coef(model_a_r_3_siemens)

# ╔═╡ ef07df25-3ed3-4d7e-8b61-e88a4ad520b6
md"""
### Visualize
"""

# ╔═╡ 0c5cd43c-118a-45c9-b4e3-e5aae61db595
function reprod_a_1_4()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df1_4_large_clean_a_canon[!, :x1], df1_4_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df1_4_medium_clean_a_canon[!, :x1], df1_4_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df1_4_small_clean_a_canon[!, :x1], df1_4_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_3_canon, linestyle=:dashdot)
    if co_a_r_3_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_3_canon[2]; digits=3))x+$(trunc(co_a_r_3_canon[1]; digits=3)) \nr = $(trunc(r_a_r_3_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_3_canon[2]; digits=3))x$(trunc(co_a_r_3_canon[1]; digits=3)) \nr = $(trunc(r_a_r_3_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_4_large_clean_a_geee[!, :x1], df1_4_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df1_4_medium_clean_a_geee[!, :x1], df1_4_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_3_geee, linestyle=:dashdot)
    if co_a_r_3_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_3_geee[2]; digits=3))x+$(trunc(co_a_r_3_geee[1]; digits=3)) \nr = $(trunc(r_a_r_3_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_3_geee[2]; digits=3))x$(trunc(co_a_r_3_geee[1]; digits=3)) \nr = $(trunc(r_a_r_3_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_4_large_clean_a_philips[!, :x1], df1_4_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df1_4_medium_clean_a_philips[!, :x1], df1_4_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_3_philips, linestyle=:dashdot)
    if co_a_r_3_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_3_philips[2]; digits=3))x+$(trunc(co_a_r_3_philips[1]; digits=3)) \nr = $(trunc(r_a_r_3_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_3_philips[2]; digits=3))x$(trunc(co_a_r_3_philips[1]; digits=3)) \nr = $(trunc(r_a_r_3_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_4_large_clean_a_siemens[!, :x1], df1_4_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_4_medium_clean_a_siemens[!, :x1], df1_4_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_3_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_a_r_3_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_3_siemens[2]; digits=3))x+$(trunc(co_a_r_3_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_3_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_3_siemens[2]; digits=3))x$(trunc(co_a_r_3_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_3_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_3_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_3_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Calcium Scoring (1 vs 4)",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ e0d8d102-d56b-4ce0-be63-343c2fa5202a
with_theme(medphys_theme) do
    reprod_a_1_4()
end

# ╔═╡ aeb5b3a5-45ce-42e3-b5a5-7c1ebd570324
md"""
## 1 vs 5
"""

# ╔═╡ ce698ae4-597b-4fe1-9364-f913ccbbd353
begin
    df1_5_large_clean_a = innerjoin(df1_large_a, df5_large_a; on=[:vendor, :inserts], makeunique=true)
    df1_5_medium_clean_a = innerjoin(df1_medium_a, df5_medium_a; on=[:vendor, :inserts], makeunique=true)
    df1_5_small_clean_a = innerjoin(df1_small_a, df5_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 07b17ba4-f139-4c00-b024-96b75a349b09
md"""
### Canon
"""

# ╔═╡ 184d730e-1197-4466-b95d-8a2ec1a5c4ea
begin
    df1_5_large_clean_a_canon, df1_5_large_clean_a_geee, df1_5_large_clean_a_philips, df1_5_large_clean_a_siemens = groupby(df1_5_large_clean_a, :vendor)
    df1_5_medium_clean_a_canon, df1_5_medium_clean_a_geee, df1_5_medium_clean_a_philips, df1_5_medium_clean_a_siemens = groupby(df1_5_medium_clean_a, :vendor)
    df1_5_small_clean_a_canon = df1_5_small_clean_a
end;

# ╔═╡ 38d0b5d1-974a-45a7-9780-70fa45c4dafe
let
    r_array = [Tuple(df1_5_large_clean_a_canon[!, :x1])..., Tuple(df1_5_medium_clean_a_canon[!, :x1])..., Tuple(df1_5_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_a_canon[!, :x1_1])..., Tuple(df1_5_medium_clean_a_canon[!, :x1_1])..., Tuple(df1_5_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_4_canon
    model_a_r_4_canon = lm(@formula(Y ~ X), data)
    global r_a_r_4_canon
    r_a_r_4_canon = GLM.r2(model_a_r_4_canon)
    global rms_values_a_r_4_canon
    rms_values_a_r_4_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_4_canon))
    ]
end

# ╔═╡ 5d12a5f4-ac7f-4242-b7ec-eb386c54eebf
begin
    X_a_r_4_canon = DataFrame(X=collect(1:1000))
    pred_a_r_4_canon = GLM.predict(model_a_r_4_canon, X_a_r_4_canon)
end

# ╔═╡ 3247670e-ae2f-45b1-889e-0fe968ea3fda
co_a_r_4_canon = coef(model_a_r_4_canon)

# ╔═╡ d8992eef-e000-48fa-8cc3-28d7b30c1527
md"""
### GE
"""

# ╔═╡ 1e170158-1bdc-4bae-890e-e4d1d80d0ad7
let
    r_array = [Tuple(df1_5_large_clean_a_geee[!, :x1])..., Tuple(df1_5_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_a_geee[!, :x1_1])..., Tuple(df1_5_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_4_geee
    model_a_r_4_geee = lm(@formula(Y ~ X), data)
    global r_a_r_4_geee
    r_a_r_4_geee = GLM.r2(model_a_r_4_geee)
    global rms_values_a_r_4_geee
    rms_values_a_r_4_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_4_geee))
    ]
end

# ╔═╡ 34b78437-0868-4634-a261-40bc2f120edc
begin
    X_a_r_4_geee = DataFrame(X=collect(1:1000))
    pred_a_r_4_geee = GLM.predict(model_a_r_4_geee, X_a_r_4_geee)
end

# ╔═╡ 6f252fc5-6af7-4a11-b472-e2d7211a3837
co_a_r_4_geee = coef(model_a_r_4_geee)

# ╔═╡ 6d21d3d0-62e1-4d7a-915f-c56372e63ef7
md"""
### Philips
"""

# ╔═╡ a79c64ec-5fe3-4db8-91a2-58264d73bb9a
let
    r_array = [Tuple(df1_5_large_clean_a_philips[!, :x1])..., Tuple(df1_5_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_a_philips[!, :x1_1])..., Tuple(df1_5_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_4_philips
    model_a_r_4_philips = lm(@formula(Y ~ X), data)
    global r_a_r_4_philips
    r_a_r_4_philips = GLM.r2(model_a_r_4_philips)
    global rms_values_a_r_4_philips
    rms_values_a_r_4_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_4_philips))
    ]
end

# ╔═╡ 466644b0-7c54-457a-98ce-15945d6cfacb
begin
    X_a_r_4_philips = DataFrame(X=collect(1:1000))
    pred_a_r_4_philips = GLM.predict(model_a_r_4_philips, X_a_r_4_philips)
end

# ╔═╡ b7b03903-53a0-44ec-a522-ad3683bfb0a1
co_a_r_4_philips = coef(model_a_r_4_philips)

# ╔═╡ 3398bb56-7147-4c16-831c-647a917a690b
md"""
### Siemens
"""

# ╔═╡ 457e7344-2bdb-4b86-965f-adfcf07b28b1
let
    r_array = [Tuple(df1_5_large_clean_a_siemens[!, :x1])..., Tuple(df1_5_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df1_5_large_clean_a_siemens[!, :x1_1])..., Tuple(df1_5_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_4_siemens
    model_a_r_4_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_4_siemens
    r_a_r_4_siemens = GLM.r2(model_a_r_4_siemens)
    global rms_values_a_r_4_siemens
    rms_values_a_r_4_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_4_siemens))
    ]
end

# ╔═╡ e0b819ed-7240-4a77-9a49-46bac20f884f
begin
    X_a_r_4_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_4_siemens = GLM.predict(model_a_r_4_siemens, X_a_r_4_siemens)
end

# ╔═╡ c290b1b4-6e02-42ec-b880-1619ab746391
co_a_r_4_siemens = coef(model_a_r_4_siemens)

# ╔═╡ 0c5c047e-2c4a-44a1-b0e5-cd6b97a58e56
md"""
### Visualize
"""

# ╔═╡ ad895815-ca2c-4167-bb73-ed8789ec7e6c
function reprod_a_1_5()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df1_5_large_clean_a_canon[!, :x1], df1_5_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df1_5_medium_clean_a_canon[!, :x1], df1_5_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df1_5_small_clean_a_canon[!, :x1], df1_5_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_4_canon, linestyle=:dashdot)
    if co_a_r_4_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_4_canon[2]; digits=3))x+$(trunc(co_a_r_4_canon[1]; digits=3)) \nr = $(trunc(r_a_r_4_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_4_canon[2]; digits=3))x$(trunc(co_a_r_4_canon[1]; digits=3)) \nr = $(trunc(r_a_r_4_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df1_5_large_clean_a_geee[!, :x1], df1_5_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df1_5_medium_clean_a_geee[!, :x1], df1_5_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_4_geee, linestyle=:dashdot)
    if co_a_r_4_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_4_geee[2]; digits=3))x+$(trunc(co_a_r_4_geee[1]; digits=3)) \nr = $(trunc(r_a_r_4_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_4_geee[2]; digits=3))x$(trunc(co_a_r_4_geee[1]; digits=3)) \nr = $(trunc(r_a_r_4_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df1_5_large_clean_a_philips[!, :x1], df1_5_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df1_5_medium_clean_a_philips[!, :x1], df1_5_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_4_philips, linestyle=:dashdot)
    if co_a_r_4_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_4_philips[2]; digits=3))x+$(trunc(co_a_r_4_philips[1]; digits=3)) \nr = $(trunc(r_a_r_4_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_4_philips[2]; digits=3))x$(trunc(co_a_r_4_philips[1]; digits=3)) \nr = $(trunc(r_a_r_4_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df1_5_large_clean_a_siemens[!, :x1], df1_5_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df1_5_medium_clean_a_siemens[!, :x1], df1_5_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_4_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_a_r_4_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_4_siemens[2]; digits=3))x+$(trunc(co_a_r_4_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_4_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_4_siemens[2]; digits=3))x$(trunc(co_a_r_4_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_4_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_4_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_4_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Calcium Scoring (1 vs 5)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ 30844636-edd0-4b70-b0c8-7fb3b1d1e71f
with_theme(medphys_theme) do
    reprod_a_1_5()
end

# ╔═╡ ffa9a017-24b2-4f7d-b9b0-792ed0378090
md"""
# Small Phantom - Integrated
"""

# ╔═╡ 033f9b89-4f88-42bc-9b14-92f7906808ec
md"""
## 6 vs 7
"""

# ╔═╡ e6385cf2-2c54-450c-85b2-f444cde4af6f
begin
    df6_7_large_clean = innerjoin(df6_large, df7_large; on=[:vendor, :inserts], makeunique=true)
    df6_7_medium_clean = innerjoin(df6_medium, df7_medium; on=[:vendor, :inserts], makeunique=true)
    df6_7_small_clean = innerjoin(df6_small, df7_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 64e7bd22-bc6e-4920-8c80-46ab542ed38d
md"""
### Canon
"""

# ╔═╡ 2993e116-34c2-42cf-8bf4-3aa9919e55b3
begin
    df6_7_large_clean_canon, df6_7_large_clean_geee, df6_7_large_clean_philips, df6_7_large_clean_siemens = groupby(df6_7_large_clean, :vendor)
    df6_7_medium_clean_canon, df6_7_medium_clean_geee, df6_7_medium_clean_philips, df6_7_medium_clean_siemens = groupby(df6_7_medium_clean, :vendor)
    df6_7_small_clean_canon = df6_7_small_clean
end;

# ╔═╡ 7b1292b6-5f17-4f9b-91c1-22442087834d
let
    r_array = [Tuple(df6_7_large_clean_canon[!, :x1])..., Tuple(df6_7_medium_clean_canon[!, :x1])..., Tuple(df6_7_small_clean_canon[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_canon[!, :x1_1])..., Tuple(df6_7_medium_clean_canon[!, :x1_1])..., Tuple(df6_7_small_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_5_canon
    model_i_r_5_canon = lm(@formula(Y ~ X), data)
    global r_i_r_5_canon
    r_i_r_5_canon = GLM.r2(model_i_r_5_canon)
    global rms_values_i_r_5_canon
    rms_values_i_r_5_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_5_canon))
    ]
end

# ╔═╡ 00cdfa1b-8e1f-4b43-8630-3f7718ff7ec8
begin
    X_i_r_5_canon = DataFrame(X=collect(1:1000))
    pred_i_r_5_canon = GLM.predict(model_i_r_5_canon, X_i_r_5_canon)
end

# ╔═╡ 699fbbdf-0729-4372-8617-d909f2c663eb
co_i_r_5_canon = coef(model_i_r_5_canon)

# ╔═╡ 4577afea-26c0-4be0-98ea-e0b5095b5baa
md"""
### GE
"""

# ╔═╡ 3a116456-4c27-4ce1-b870-7c4d06f34451
let
    r_array = [Tuple(df6_7_large_clean_geee[!, :x1])..., Tuple(df6_7_medium_clean_geee[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_geee[!, :x1_1])..., Tuple(df6_7_medium_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_5_geee
    model_i_r_5_geee = lm(@formula(Y ~ X), data)
    global r_i_r_5_geee
    r_i_r_5_geee = GLM.r2(model_i_r_5_geee)
    global rms_values_i_r_5_geee
    rms_values_i_r_5_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_5_geee))
    ]
end

# ╔═╡ 6a067a65-61a8-483f-9686-aa91645014eb
begin
    X_i_r_5_geee = DataFrame(X=collect(1:1000))
    pred_i_r_5_geee = GLM.predict(model_i_r_5_geee, X_i_r_5_geee)
end

# ╔═╡ fff52d0b-2956-46d6-b3a4-fa8045379b65
co_i_r_5_geee = coef(model_i_r_5_geee)

# ╔═╡ 061e0e6e-730f-45b4-bf00-5f8b91f0c507
md"""
### Philips
"""

# ╔═╡ 49cdf7dd-d6ac-499a-937b-7ec2e05b04bb
let
    r_array = [Tuple(df6_7_large_clean_philips[!, :x1])..., Tuple(df6_7_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_philips[!, :x1_1])..., Tuple(df6_7_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_5_philips
    model_i_r_5_philips = lm(@formula(Y ~ X), data)
    global r_i_r_5_philips
    r_i_r_5_philips = GLM.r2(model_i_r_5_philips)
    global rms_values_i_r_5_philips
    rms_values_i_r_5_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_5_philips))
    ]
end

# ╔═╡ 3bccc16c-d42f-4784-9e63-47a3d3e42d58
begin
    X_i_r_5_philips = DataFrame(X=collect(1:1000))
    pred_i_r_5_philips = GLM.predict(model_i_r_5_philips, X_i_r_5_philips)
end

# ╔═╡ 081bc4e5-ad09-4749-b99f-3be537efe158
co_i_r_5_philips = coef(model_i_r_5_philips)

# ╔═╡ 1d895d11-02ed-4838-abe5-8bca36fd4560
md"""
### Siemens
"""

# ╔═╡ 1357670e-a30b-45ff-8d29-9bd07000ed8a
let
    r_array = [Tuple(df6_7_large_clean_siemens[!, :x1])..., Tuple(df6_7_medium_clean_siemens[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_siemens[!, :x1_1])..., Tuple(df6_7_medium_clean_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_5_siemens
    model_i_r_5_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_5_siemens
    r_i_r_5_siemens = GLM.r2(model_i_r_5_siemens)
    global rms_values_i_r_5_siemens
    rms_values_i_r_5_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_5_siemens))
    ]
end

# ╔═╡ 65f70586-4d48-4c59-b3ec-abf7b96d6031
begin
    X_i_r_5_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_5_siemens = GLM.predict(model_i_r_5_siemens, X_i_r_5_siemens)
end

# ╔═╡ 3d1a536b-218e-4d83-847c-ce088b7bb356
co_i_r_5_siemens = coef(model_i_r_5_siemens)

# ╔═╡ 7992f0a2-9012-4801-821c-e917df0e06db
md"""
### vendors
"""

# ╔═╡ 4bd87cfe-f24a-4f14-8d29-094fab2ee6e8
function reprod_i_6_7()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df6_7_large_clean_canon[!, :x1], df6_7_large_clean_canon[!, :x1_1])
    scatter!(ax1, df6_7_medium_clean_canon[!, :x1], df6_7_medium_clean_canon[!, :x1_1])
    scatter!(ax1, df6_7_small_clean_canon[!, :x1], df6_7_small_clean_canon[!, :x1_1], color=:red)
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_5_canon, linestyle=:dashdot)
    if co_i_r_5_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_5_canon[2]; digits=3))x+$(trunc(co_i_r_5_canon[1]; digits=3)) \nr = $(trunc(r_i_r_5_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_5_canon[2]; digits=3))x$(trunc(co_i_r_5_canon[1]; digits=3)) \nr = $(trunc(r_i_r_5_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_7_large_clean_geee[!, :x1], df6_7_large_clean_geee[!, :x1_1])
    scatter!(ax2, df6_7_medium_clean_geee[!, :x1], df6_7_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df6_7_small_clean_geee[!, :x1], df6_7_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_5_geee, linestyle=:dashdot)
    if co_i_r_5_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_5_geee[2]; digits=3))x+$(trunc(co_i_r_5_geee[1]; digits=3)) \nr = $(trunc(r_i_r_5_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_5_geee[2]; digits=3))x$(trunc(co_i_r_5_geee[1]; digits=3)) \nr = $(trunc(r_i_r_5_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_7_large_clean_philips[!, :x1], df6_7_large_clean_philips[!, :x1_1])
    scatter!(ax3, df6_7_medium_clean_philips[!, :x1], df6_7_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df6_7_small_clean_philips[!, :x1], df6_7_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_5_philips, linestyle=:dashdot)
    if co_i_r_5_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_5_philips[2]; digits=3))x+$(trunc(co_i_r_5_philips[1]; digits=3)) \nr = $(trunc(r_i_r_5_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_5_philips[2]; digits=3))x$(trunc(co_i_r_5_philips[1]; digits=3)) \nr = $(trunc(r_i_r_5_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_7_large_clean_siemens[!, :x1], df6_7_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_7_medium_clean_siemens[!, :x1], df6_7_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    # scatter!(ax4, df6_7_small_clean_siemens[!, :x1], df6_7_small_clean_siemens[!, :x1_1], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_5_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_i_r_5_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_5_siemens[2]; digits=3))x+$(trunc(co_i_r_5_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_5_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_5_siemens[2]; digits=3))x$(trunc(co_i_r_5_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_5_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_5_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_5_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (6 vs 7)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_6_7.png", f)
    f
end

# ╔═╡ 97489c7d-24cb-4685-bf41-76c9457f4dcc
with_theme(medphys_theme) do
    reprod_i_6_7()
end

# ╔═╡ 4c7c637a-15ab-4103-9f43-5183ca8f8d0f
md"""
## 6 vs 8
"""

# ╔═╡ db0370ca-dde6-46e3-aa47-6669b4abe3eb
begin
    df6_8_large_clean = innerjoin(df6_large, df8_large; on=[:vendor, :inserts], makeunique=true)
    df6_8_medium_clean = innerjoin(df6_medium, df8_medium; on=[:vendor, :inserts], makeunique=true)
    df6_8_small_clean = innerjoin(df6_small, df8_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ f420b178-6daf-447c-be6e-72b243c15d0f
md"""
### Canon
"""

# ╔═╡ c74f97ce-182d-47c1-91bb-b4b56554c5b1
begin
    df6_8_large_clean_canon, df6_8_large_clean_geee, df6_8_large_clean_philips, df6_8_large_clean_siemens = groupby(df6_8_large_clean, :vendor)
    df6_8_medium_clean_canon, df6_8_medium_clean_geee, df6_8_medium_clean_philips, df6_8_medium_clean_siemens = groupby(df6_8_medium_clean, :vendor)
    df6_8_small_clean_canon = df6_8_small_clean
end;

# ╔═╡ 84345749-fbac-4682-a9bf-cb9118686938
let
    r_array = [Tuple(df6_8_large_clean_canon[!, :x1])..., Tuple(df6_8_medium_clean_canon[!, :x1])..., Tuple(df6_8_small_clean_canon[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_canon[!, :x1_1])..., Tuple(df6_8_medium_clean_canon[!, :x1_1])..., Tuple(df6_8_small_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_6_canon
    model_i_r_6_canon = lm(@formula(Y ~ X), data)
    global r_i_r_6_canon
    r_i_r_6_canon = GLM.r2(model_i_r_6_canon)
    global rms_values_i_r_6_canon
    rms_values_i_r_6_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_6_canon))
    ]
end

# ╔═╡ 81da47ef-4d90-444c-bd50-defdbd6e3440
begin
    X_i_r_6_canon = DataFrame(X=collect(1:1000))
    pred_i_r_6_canon = GLM.predict(model_i_r_6_canon, X_i_r_6_canon)
end

# ╔═╡ 30d688ea-bf2f-42b2-aa73-0fd23f7936c7
co_i_r_6_canon = coef(model_i_r_6_canon)

# ╔═╡ 0f011318-bbfa-4ae4-8998-be028c481941
md"""
### GE
"""

# ╔═╡ 005830de-16c5-4aac-bf0b-60b2cc59d735
let
    r_array = [Tuple(df6_8_large_clean_geee[!, :x1])..., Tuple(df6_8_medium_clean_geee[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_geee[!, :x1_1])..., Tuple(df6_8_medium_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_6_geee
    model_i_r_6_geee = lm(@formula(Y ~ X), data)
    global r_i_r_6_geee
    r_i_r_6_geee = GLM.r2(model_i_r_6_geee)
    global rms_values_i_r_6_geee
    rms_values_i_r_6_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_6_geee))
    ]
end

# ╔═╡ a4eb6ff4-182f-4dd3-bdba-0c46dc032a28
begin
    X_i_r_6_geee = DataFrame(X=collect(1:1000))
    pred_i_r_6_geee = GLM.predict(model_i_r_6_geee, X_i_r_6_geee)
end

# ╔═╡ 673cb1d6-7a77-4c75-851b-ac134270d5d3
co_i_r_6_geee = coef(model_i_r_6_geee)

# ╔═╡ 7c65b05a-eaa4-4c14-937b-696b395a7548
md"""
### Philips
"""

# ╔═╡ 781e21a4-bfed-4b0c-b236-6b51b29ab12d
let
    r_array = [Tuple(df6_8_large_clean_philips[!, :x1])..., Tuple(df6_8_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_philips[!, :x1_1])..., Tuple(df6_8_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_6_philips
    model_i_r_6_philips = lm(@formula(Y ~ X), data)
    global r_i_r_6_philips
    r_i_r_6_philips = GLM.r2(model_i_r_6_philips)
    global rms_values_i_r_6_philips
    rms_values_i_r_6_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_6_philips))
    ]
end

# ╔═╡ 42e07123-1010-46c8-999c-525646539a07
begin
    X_i_r_6_philips = DataFrame(X=collect(1:1000))
    pred_i_r_6_philips = GLM.predict(model_i_r_6_philips, X_i_r_6_philips)
end

# ╔═╡ 492d4cc2-0239-49d9-8ac5-4307c0ebb1ee
co_i_r_6_philips = coef(model_i_r_6_philips)

# ╔═╡ 246eadad-ec40-45f8-9808-f2413184be6e
md"""
### Siemens
"""

# ╔═╡ 75452f19-00ae-4c20-98ba-86708ab301e7
let
    r_array = [Tuple(df6_8_large_clean_siemens[!, :x1])..., Tuple(df6_8_medium_clean_siemens[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_siemens[!, :x1_1])..., Tuple(df6_8_medium_clean_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_6_siemens
    model_i_r_6_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_6_siemens
    r_i_r_6_siemens = GLM.r2(model_i_r_6_siemens)
    global rms_values_i_r_6_siemens
    rms_values_i_r_6_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_6_siemens))
    ]
end

# ╔═╡ 1f56b032-e5e5-47e4-a491-87f1149b236a
begin
    X_i_r_6_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_6_siemens = GLM.predict(model_i_r_6_siemens, X_i_r_6_siemens)
end

# ╔═╡ 0c7f45c5-0b30-461a-a382-3f9fe7f5ab16
co_i_r_6_siemens = coef(model_i_r_6_siemens)

# ╔═╡ b22547d3-3c10-4e45-b108-4a1a11119da1
md"""
### Visualize
"""

# ╔═╡ c6bc875b-6849-4286-b40b-6c23e0d78721
function reprod_i_6_8()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df6_8_large_clean_canon[!, :x1], df6_8_large_clean_canon[!, :x1_1])
    scatter!(ax1, df6_8_medium_clean_canon[!, :x1], df6_8_medium_clean_canon[!, :x1_1])
    scatter!(ax1, df6_8_small_clean_canon[!, :x1], df6_8_small_clean_canon[!, :x1_1], color=:red)
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_6_canon, linestyle=:dashdot)
    if co_i_r_6_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_6_canon[2]; digits=3))x+$(trunc(co_i_r_6_canon[1]; digits=3)) \nr = $(trunc(r_i_r_6_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_6_canon[2]; digits=3))x$(trunc(co_i_r_6_canon[1]; digits=3)) \nr = $(trunc(r_i_r_6_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_8_large_clean_geee[!, :x1], df6_8_large_clean_geee[!, :x1_1])
    scatter!(ax2, df6_8_medium_clean_geee[!, :x1], df6_8_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df6_8_small_clean_geee[!, :x1], df6_8_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_6_geee, linestyle=:dashdot)
    if co_i_r_6_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_6_geee[2]; digits=3))x+$(trunc(co_i_r_6_geee[1]; digits=3)) \nr = $(trunc(r_i_r_6_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_6_geee[2]; digits=3))x$(trunc(co_i_r_6_geee[1]; digits=3)) \nr = $(trunc(r_i_r_6_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_8_large_clean_philips[!, :x1], df6_8_large_clean_philips[!, :x1_1])
    scatter!(ax3, df6_8_medium_clean_philips[!, :x1], df6_8_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df6_8_small_clean_philips[!, :x1], df6_8_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_6_philips, linestyle=:dashdot)
    if co_i_r_6_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_6_philips[2]; digits=3))x+$(trunc(co_i_r_6_philips[1]; digits=3)) \nr = $(trunc(r_i_r_6_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_6_philips[2]; digits=3))x$(trunc(co_i_r_6_philips[1]; digits=3)) \nr = $(trunc(r_i_r_6_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_8_large_clean_siemens[!, :x1], df6_8_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_8_medium_clean_siemens[!, :x1], df6_8_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_6_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_i_r_6_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_6_siemens[2]; digits=3))x+$(trunc(co_i_r_6_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_6_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_6_siemens[2]; digits=3))x$(trunc(co_i_r_6_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_6_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_6_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_6_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (6 vs 8)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ 882745b5-96eb-4cbd-b356-6fd7a2983697
with_theme(medphys_theme) do
    reprod_i_6_8()
end

# ╔═╡ fd28b420-caa3-4638-91ff-2e98b5c7fca5
md"""
## 6 vs 9
"""

# ╔═╡ 73236865-d3b2-4ec4-b747-131ef0aec338
begin
    df6_9_large_clean = innerjoin(df6_large, df9_large; on=[:vendor, :inserts], makeunique=true)
    df6_9_medium_clean = innerjoin(df6_medium, df9_medium; on=[:vendor, :inserts], makeunique=true)
    df6_9_small_clean = innerjoin(df6_small, df9_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 36c19466-f340-4549-a531-503334fae7aa
md"""
### Canon
"""

# ╔═╡ 58c01c12-d1ec-45cd-a6ad-b3653b4ed3b0
df6_9_small_clean

# ╔═╡ b75d526e-f282-49cc-85be-fac4744638eb
begin
    df6_9_large_clean_canon, df6_9_large_clean_geee, df6_9_large_clean_philips, df6_9_large_clean_siemens = groupby(df6_9_large_clean, :vendor)
    df6_9_medium_clean_canon, df6_9_medium_clean_geee, df6_9_medium_clean_philips, df6_9_medium_clean_siemens = groupby(df6_9_medium_clean, :vendor)
    df6_9_small_clean_canon = df6_9_small_clean
end;

# ╔═╡ 5b3e9035-9cfb-4175-a781-ad7db7260d40
let
    r_array = [Tuple(df6_9_large_clean_canon[!, :x1])..., Tuple(df6_9_medium_clean_canon[!, :x1])..., Tuple(df6_9_small_clean_canon[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_canon[!, :x1_1])..., Tuple(df6_9_medium_clean_canon[!, :x1_1])..., Tuple(df6_9_small_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_7_canon
    model_i_r_7_canon = lm(@formula(Y ~ X), data)
    global r_i_r_7_canon
    r_i_r_7_canon = GLM.r2(model_i_r_7_canon)
    global rms_values_i_r_7_canon
    rms_values_i_r_7_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_7_canon))
    ]
end

# ╔═╡ 7ee2f5f8-01bb-462f-90d2-960352a53a36
begin
    X_i_r_7_canon = DataFrame(X=collect(1:1000))
    pred_i_r_7_canon = GLM.predict(model_i_r_7_canon, X_i_r_7_canon)
end

# ╔═╡ 98676111-76eb-48df-bc25-c9149a3dd833
co_i_r_7_canon = coef(model_i_r_7_canon)

# ╔═╡ 525d755c-e1fa-4be9-a220-412431900416
md"""
### GE
"""

# ╔═╡ 0fb8547e-b396-4c7a-8809-f7cc7441900a
let
    r_array = [Tuple(df6_9_large_clean_geee[!, :x1])..., Tuple(df6_9_medium_clean_geee[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_geee[!, :x1_1])..., Tuple(df6_9_medium_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_7_geee
    model_i_r_7_geee = lm(@formula(Y ~ X), data)
    global r_i_r_7_geee
    r_i_r_7_geee = GLM.r2(model_i_r_7_geee)
    global rms_values_i_r_7_geee
    rms_values_i_r_7_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_7_geee))
    ]
end

# ╔═╡ 46dd427a-345d-4487-b24f-31b6917a00d2
begin
    X_i_r_7_geee = DataFrame(X=collect(1:1000))
    pred_i_r_7_geee = GLM.predict(model_i_r_7_geee, X_i_r_7_geee)
end

# ╔═╡ 447725c6-ab92-4cc3-9a78-2d3022f1258a
co_i_r_7_geee = coef(model_i_r_7_geee)

# ╔═╡ da5ccfbd-2ec8-4f55-924a-df687711e305
md"""
### Philips
"""

# ╔═╡ a473c7f9-9885-42e6-a62e-f533b7ad4f24
let
    r_array = [Tuple(df6_9_large_clean_philips[!, :x1])..., Tuple(df6_9_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_philips[!, :x1_1])..., Tuple(df6_9_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_7_philips
    model_i_r_7_philips = lm(@formula(Y ~ X), data)
    global r_i_r_7_philips
    r_i_r_7_philips = GLM.r2(model_i_r_7_philips)
    global rms_values_i_r_7_philips
    rms_values_i_r_7_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_7_philips))
    ]
end

# ╔═╡ 6d973a09-533b-4a40-bfa1-e72e65237082
begin
    X_i_r_7_philips = DataFrame(X=collect(1:1000))
    pred_i_r_7_philips = GLM.predict(model_i_r_7_philips, X_i_r_7_philips)
end

# ╔═╡ 22196e29-f3af-4727-a729-ecac3267af90
co_i_r_7_philips = coef(model_i_r_7_philips)

# ╔═╡ 9be70624-cbff-4083-bfe5-1936a1d08c03
md"""
### Siemens
"""

# ╔═╡ e3c253ec-d7eb-4c3b-8f67-73620228f994
let
    r_array = [Tuple(df6_9_large_clean_siemens[!, :x1])..., Tuple(df6_9_medium_clean_siemens[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_siemens[!, :x1_1])..., Tuple(df6_9_medium_clean_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_7_siemens
    model_i_r_7_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_7_siemens
    r_i_r_7_siemens = GLM.r2(model_i_r_7_siemens)
    global rms_values_i_r_7_siemens
    rms_values_i_r_7_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_7_siemens))
    ]
end

# ╔═╡ d409efbf-2a66-4385-9d21-8111722832a0
begin
    X_i_r_7_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_7_siemens = GLM.predict(model_i_r_7_siemens, X_i_r_7_siemens)
end

# ╔═╡ a28122c7-df21-48ab-993b-dd1c7bb966a3
co_i_r_7_siemens = coef(model_i_r_7_siemens)

# ╔═╡ afd1a63a-4e61-4451-8d5a-031d833d3ff2
md"""
### Visualize
"""

# ╔═╡ 0513263c-720f-4e7b-842a-ba65f99b9900
function reprod_i_6_9()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df6_9_large_clean_canon[!, :x1], df6_9_large_clean_canon[!, :x1_1])
    scatter!(ax1, df6_9_medium_clean_canon[!, :x1], df6_9_medium_clean_canon[!, :x1_1])
    scatter!(ax1, df6_9_small_clean_canon[!, :x1], df6_9_small_clean_canon[!, :x1_1], color=:red)
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_7_canon, linestyle=:dashdot)
    if co_i_r_7_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_7_canon[2]; digits=3))x+$(trunc(co_i_r_7_canon[1]; digits=3)) \nr = $(trunc(r_i_r_7_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_7_canon[2]; digits=3))x$(trunc(co_i_r_7_canon[1]; digits=3)) \nr = $(trunc(r_i_r_7_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_9_large_clean_geee[!, :x1], df6_9_large_clean_geee[!, :x1_1])
    scatter!(ax2, df6_9_medium_clean_geee[!, :x1], df6_9_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df6_9_small_clean_geee[!, :x1], df6_9_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_7_geee, linestyle=:dashdot)
    if co_i_r_7_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_7_geee[2]; digits=3))x+$(trunc(co_i_r_7_geee[1]; digits=3)) \nr = $(trunc(r_i_r_7_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_7_geee[2]; digits=3))x$(trunc(co_i_r_7_geee[1]; digits=3)) \nr = $(trunc(r_i_r_7_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_9_large_clean_philips[!, :x1], df6_9_large_clean_philips[!, :x1_1])
    scatter!(ax3, df6_9_medium_clean_philips[!, :x1], df6_9_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df6_9_small_clean_philips[!, :x1], df6_9_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_7_philips, linestyle=:dashdot)
    if co_i_r_7_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_7_philips[2]; digits=3))x+$(trunc(co_i_r_7_philips[1]; digits=3)) \nr = $(trunc(r_i_r_7_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_7_philips[2]; digits=3))x$(trunc(co_i_r_7_philips[1]; digits=3)) \nr = $(trunc(r_i_r_7_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_9_large_clean_siemens[!, :x1], df6_9_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_9_medium_clean_siemens[!, :x1], df6_9_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_7_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_i_r_7_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_7_siemens[2]; digits=3))x+$(trunc(co_i_r_7_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_7_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_7_siemens[2]; digits=3))x$(trunc(co_i_r_7_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_7_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_7_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_7_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (6 vs 9)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ 5be7b2db-1b1d-4633-8c2f-52549e083fa9
with_theme(medphys_theme) do
    reprod_i_6_9()
end

# ╔═╡ 78f81079-f267-42a9-b59f-5f75eed991e8
md"""
## 6 vs 10
"""

# ╔═╡ 23ba4b57-d81c-4f7e-b7fc-7f3af1f5e121
begin
    df6_10_large_clean = innerjoin(df6_large, df10_large; on=[:vendor, :inserts], makeunique=true)
    df6_10_medium_clean = innerjoin(df6_medium, df10_medium; on=[:vendor, :inserts], makeunique=true)
    df6_10_small_clean = innerjoin(df6_small, df10_small; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ dda5c63e-42d3-4192-89fc-8730ac168842
let
    r_array = [Tuple(df6_10_large_clean[!, :x1])..., Tuple(df6_10_medium_clean[!, :x1])..., Tuple(df6_10_small_clean[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean[!, :x1_1])..., Tuple(df6_10_medium_clean[!, :x1_1])..., Tuple(df6_10_small_clean[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_8
    model_i_r_8 = lm(@formula(Y ~ X), data)
    global r_i_r_8
    r_i_r_8 = GLM.r2(model_i_r_8)
    global rms_values_i_r_8
    rms_values_i_r_8 = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_8))
    ]
end

# ╔═╡ 114bf88c-4daa-4c28-a0e1-4296cb1257dc
begin
    X_i_r_8 = DataFrame(X=collect(1:1000))
    pred_i_r_8 = GLM.predict(model_i_r_8, X_i_r_8)
end

# ╔═╡ 646ff2dc-d490-4cfd-8618-b35982fec465
co_i_r_8 = coef(model_i_r_8)

# ╔═╡ 9461e6d1-0ce9-4bd9-a68c-ece5fcebc6a5
md"""
### Canon
"""

# ╔═╡ 7d1f6a58-da38-4ff4-b0f1-822d62b71f10
df6_10_small_clean

# ╔═╡ dc9784f5-88ca-4c23-ad35-648979f35814
begin
    df6_10_large_clean_canon, df6_10_large_clean_geee, df6_10_large_clean_philips, df6_10_large_clean_siemens = groupby(df6_10_large_clean, :vendor)
    df6_10_medium_clean_canon, df6_10_medium_clean_geee, df6_10_medium_clean_philips, df6_10_medium_clean_siemens = groupby(df6_10_medium_clean, :vendor)
    df6_10_small_clean_geee = df6_10_small_clean
end;

# ╔═╡ 7ab94d89-9cee-4328-a35c-b519304ad5b8
let
    r_array = [Tuple(df6_10_large_clean_canon[!, :x1])..., Tuple(df6_10_medium_clean_canon[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_canon[!, :x1_1])..., Tuple(df6_10_medium_clean_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_8_canon
    model_i_r_8_canon = lm(@formula(Y ~ X), data)
    global r_i_r_8_canon
    r_i_r_8_canon = GLM.r2(model_i_r_8_canon)
    global rms_values_i_r_8_canon
    rms_values_i_r_8_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_8_canon))
    ]
end

# ╔═╡ 7b8ea451-02a2-4369-a2bf-37c98221d400
begin
    X_i_r_8_canon = DataFrame(X=collect(1:1000))
    pred_i_r_8_canon = GLM.predict(model_i_r_8_canon, X_i_r_8_canon)
end

# ╔═╡ 28d8585e-6f13-4ab8-8e5b-3fcee608c7fd
co_i_r_8_canon = coef(model_i_r_8_canon)

# ╔═╡ d66bac7a-4952-4c7e-ad1f-053c597a3cb1
md"""
### GE
"""

# ╔═╡ abab28f2-bb87-46ac-9470-986ca8de668f
let
    r_array = [Tuple(df6_10_large_clean_geee[!, :x1])..., Tuple(df6_10_medium_clean_geee[!, :x1])..., Tuple(df6_10_small_clean_geee[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_geee[!, :x1_1])..., Tuple(df6_10_medium_clean_geee[!, :x1_1])..., Tuple(df6_10_small_clean_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_8_geee
    model_i_r_8_geee = lm(@formula(Y ~ X), data)
    global r_i_r_8_geee
    r_i_r_8_geee = GLM.r2(model_i_r_8_geee)
    global rms_values_i_r_8_geee
    rms_values_i_r_8_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_8_geee))
    ]
end

# ╔═╡ 61429df4-ac03-40ef-91d8-8658a9d080bf
begin
    X_i_r_8_geee = DataFrame(X=collect(1:1000))
    pred_i_r_8_geee = GLM.predict(model_i_r_8_geee, X_i_r_8_geee)
end

# ╔═╡ d14099a9-7559-4940-b275-29c82a341068
co_i_r_8_geee = coef(model_i_r_8_geee)

# ╔═╡ 624d2696-eb16-434c-a8fa-93e4ec0f680c
md"""
### Philips
"""

# ╔═╡ 4fd3ce19-c1ee-452f-b098-901ba46affe5
let
    r_array = [Tuple(df6_10_large_clean_philips[!, :x1])..., Tuple(df6_10_medium_clean_philips[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_philips[!, :x1_1])..., Tuple(df6_10_medium_clean_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_8_philips
    model_i_r_8_philips = lm(@formula(Y ~ X), data)
    global r_i_r_8_philips
    r_i_r_8_philips = GLM.r2(model_i_r_8_philips)
    global rms_values_i_r_8_philips
    rms_values_i_r_8_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_8_philips))
    ]
end

# ╔═╡ f8cd4e74-72ab-42e5-ad63-0f947fcb9c72
begin
    X_i_r_8_philips = DataFrame(X=collect(1:1000))
    pred_i_r_8_philips = GLM.predict(model_i_r_8_philips, X_i_r_8_philips)
end

# ╔═╡ 67eda3a2-5b52-4927-ac9a-9c867e1b0c75
co_i_r_8_philips = coef(model_i_r_8_philips)

# ╔═╡ c0e61fb4-6e71-4b17-8444-99eb39f0db20
md"""
### Siemens
"""

# ╔═╡ 5878377f-b0b2-4e9d-97b4-379b141875ff
let
    r_array = [Tuple(df6_10_large_clean_siemens[!, :x1])..., Tuple(df6_10_medium_clean_siemens[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_siemens[!, :x1_1])..., Tuple(df6_10_medium_clean_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_i_r_8_siemens
    model_i_r_8_siemens = lm(@formula(Y ~ X), data)
    global r_i_r_8_siemens
    r_i_r_8_siemens = GLM.r2(model_i_r_8_siemens)
    global rms_values_i_r_8_siemens
    rms_values_i_r_8_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_i_r_8_siemens))
    ]
end

# ╔═╡ 337eb752-7651-4aa4-b883-80ecaae095aa
begin
    X_i_r_8_siemens = DataFrame(X=collect(1:1000))
    pred_i_r_8_siemens = GLM.predict(model_i_r_8_siemens, X_i_r_8_siemens)
end

# ╔═╡ d8fc55fe-b4c5-4662-9b66-29cca9b9f5c6
co_i_r_8_siemens = coef(model_i_r_8_siemens)

# ╔═╡ c495aecb-b875-40d8-9962-ae618ec0042a
md"""
### Visualize
"""

# ╔═╡ 418f0615-3f3b-4188-a2f9-98a9b16bf81c
function reprod_i_6_10()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    scatter!(ax1, df6_10_large_clean_canon[!, :x1], df6_10_large_clean_canon[!, :x1_1])
    scatter!(ax1, df6_10_medium_clean_canon[!, :x1], df6_10_medium_clean_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_i_r_8_canon, linestyle=:dashdot)
    if co_i_r_8_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_8_canon[2]; digits=3))x+$(trunc(co_i_r_8_canon[1]; digits=3)) \nr = $(trunc(r_i_r_8_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_i_r_8_canon[2]; digits=3))x$(trunc(co_i_r_8_canon[1]; digits=3)) \nr = $(trunc(r_i_r_8_canon; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_10_large_clean_geee[!, :x1], df6_10_large_clean_geee[!, :x1_1])
    scatter!(ax2, df6_10_medium_clean_geee[!, :x1], df6_10_medium_clean_geee[!, :x1_1])
    # scatter!(ax2, df6_10_small_clean_geee[!, :x1], df6_10_small_clean_geee[!, :x1_1], color=:red)
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_i_r_8_geee, linestyle=:dashdot)
    if co_i_r_8_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_8_geee[2]; digits=3))x+$(trunc(co_i_r_8_geee[1]; digits=3)) \nr = $(trunc(r_i_r_8_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_i_r_8_geee[2]; digits=3))x$(trunc(co_i_r_8_geee[1]; digits=3)) \nr = $(trunc(r_i_r_8_geee; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_10_large_clean_philips[!, :x1], df6_10_large_clean_philips[!, :x1_1])
    scatter!(ax3, df6_10_medium_clean_philips[!, :x1], df6_10_medium_clean_philips[!, :x1_1])
    # scatter!(ax3, df6_10_small_clean_philips[!, :x1], df6_10_small_clean_philips[!, :x1_1], color=:red)
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_i_r_8_philips, linestyle=:dashdot)
    if co_i_r_8_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_8_philips[2]; digits=3))x+$(trunc(co_i_r_8_philips[1]; digits=3)) \nr = $(trunc(r_i_r_8_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_i_r_8_philips[2]; digits=3))x$(trunc(co_i_r_8_philips[1]; digits=3)) \nr = $(trunc(r_i_r_8_philips; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_10_large_clean_siemens[!, :x1], df6_10_large_clean_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_10_medium_clean_siemens[!, :x1], df6_10_medium_clean_siemens[!, :x1_1], label="Medium Inserts")
    # scatter!(ax4, df6_10_small_clean_siemens[!, :x1], df6_10_small_clean_siemens[!, :x1_1], label="Small Inserts", color=:red)
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_i_r_8_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_i_r_8_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_8_siemens[2]; digits=3))x+$(trunc(co_i_r_8_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_8_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_i_r_8_siemens[2]; digits=3))x$(trunc(co_i_r_8_siemens[1]; digits=3)) \nr = $(trunc(r_i_r_8_siemens; digits=3)) \nRMSE: $(trunc(rms_values_i_r_8_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_i_r_8_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Integrated Calcium Mass (6 vs 10)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ b52e1a6d-853d-4b31-bb9a-127226f7946e
with_theme(medphys_theme) do
    reprod_i_6_10()
end

# ╔═╡ b244bd34-b91a-47f9-bfe4-958b4e3ad635
md"""
# Small Phantom - Agatston
"""

# ╔═╡ 6dd95339-fe3d-4efa-8f8f-c59edd4e435a
md"""
## 6 vs 7
"""

# ╔═╡ cc3fb3ee-283d-4b5b-b102-587bd0c40acf
begin
    df6_7_large_clean_a = innerjoin(df6_large_a, df7_large_a; on=[:vendor, :inserts], makeunique=true)
    df6_7_medium_clean_a = innerjoin(df6_medium_a, df7_medium_a; on=[:vendor, :inserts], makeunique=true)
    df6_7_small_clean_a = innerjoin(df6_small_a, df7_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 556e18bb-16f3-4571-af38-ca56ba90d601
md"""
### Canon
"""

# ╔═╡ b4c5ca8e-d1ab-41f5-aec2-9d6bc3138557
begin
    df6_7_large_clean_a_canon, df6_7_large_clean_a_geee, df6_7_large_clean_a_philips, df6_7_large_clean_a_siemens = groupby(df6_7_large_clean_a, :vendor)
    df6_7_medium_clean_a_canon, df6_7_medium_clean_a_geee, df6_7_medium_clean_a_philips, df6_7_medium_clean_a_siemens = groupby(df6_7_medium_clean_a, :vendor)
    df6_7_small_clean_a_canon = df6_7_small_clean_a
end;

# ╔═╡ 69182167-1cc9-4414-9402-4586fa694f4b
let
    r_array = [Tuple(df6_7_large_clean_a_canon[!, :x1])..., Tuple(df6_7_medium_clean_a_canon[!, :x1])..., Tuple(df6_7_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_a_canon[!, :x1_1])..., Tuple(df6_7_medium_clean_a_canon[!, :x1_1])..., Tuple(df6_7_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_5_canon
    model_a_r_5_canon = lm(@formula(Y ~ X), data)
    global r_a_r_5_canon
    r_a_r_5_canon = GLM.r2(model_a_r_5_canon)
    global rms_values_a_r_5_canon
    rms_values_a_r_5_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_5_canon))
    ]
end

# ╔═╡ e036ae91-541f-45d4-9b7b-9d79f72de86f
begin
    X_a_r_5_canon = DataFrame(X=collect(1:1000))
    pred_a_r_5_canon = GLM.predict(model_a_r_5_canon, X_a_r_5_canon)
end

# ╔═╡ d557365b-3700-4358-9a26-8ea131a07469
co_a_r_5_canon = coef(model_a_r_5_canon)

# ╔═╡ 5ef9cbe2-2852-483e-bdd4-62b1bdda6193
md"""
### GE
"""

# ╔═╡ 03d45881-14cf-4110-9512-12c743c160db
let
    r_array = [Tuple(df6_7_large_clean_a_geee[!, :x1])..., Tuple(df6_7_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_a_geee[!, :x1_1])..., Tuple(df6_7_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_5_geee
    model_a_r_5_geee = lm(@formula(Y ~ X), data)
    global r_a_r_5_geee
    r_a_r_5_geee = GLM.r2(model_a_r_5_geee)
    global rms_values_a_r_5_geee
    rms_values_a_r_5_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_5_geee))
    ]
end

# ╔═╡ 216c540f-2d9a-499e-b29f-5c84d3a7e22a
begin
    X_a_r_5_geee = DataFrame(X=collect(1:1000))
    pred_a_r_5_geee = GLM.predict(model_a_r_5_geee, X_a_r_5_geee)
end

# ╔═╡ 91545459-fa68-4b86-b99d-b80e2db7ea74
co_a_r_5_geee = coef(model_a_r_5_geee)

# ╔═╡ a4546417-20b2-4e1d-ab33-671d89610fdb
md"""
### Philips
"""

# ╔═╡ cebe4f0b-5fc0-49fd-9ce8-4fe024097142
let
    r_array = [Tuple(df6_7_large_clean_a_philips[!, :x1])..., Tuple(df6_7_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_a_philips[!, :x1_1])..., Tuple(df6_7_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_5_philips
    model_a_r_5_philips = lm(@formula(Y ~ X), data)
    global r_a_r_5_philips
    r_a_r_5_philips = GLM.r2(model_a_r_5_philips)
    global rms_values_a_r_5_philips
    rms_values_a_r_5_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_5_philips))
    ]
end

# ╔═╡ 1668c2c6-b278-419e-8eae-bf6ee8da7b6a
begin
    X_a_r_5_philips = DataFrame(X=collect(1:1000))
    pred_a_r_5_philips = GLM.predict(model_a_r_5_philips, X_a_r_5_philips)
end

# ╔═╡ b10ade70-d733-4338-8fe1-090128ec172f
co_a_r_5_philips = coef(model_a_r_5_philips)

# ╔═╡ b135ec0b-dec9-442b-9ccf-63122eca6cfa
md"""
### Siemens
"""

# ╔═╡ f43b75b1-c6c2-4530-a73f-02b1151bd0b6
let
    r_array = [Tuple(df6_7_large_clean_a_siemens[!, :x1])..., Tuple(df6_7_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df6_7_large_clean_a_siemens[!, :x1_1])..., Tuple(df6_7_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_5_siemens
    model_a_r_5_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_5_siemens
    r_a_r_5_siemens = GLM.r2(model_a_r_5_siemens)
    global rms_values_a_r_5_siemens
    rms_values_a_r_5_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_5_siemens))
    ]
end

# ╔═╡ a212e5a5-631f-49fa-88b1-5a3af605d1d2
begin
    X_a_r_5_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_5_siemens = GLM.predict(model_a_r_5_siemens, X_a_r_5_siemens)
end

# ╔═╡ 8822f9e1-9829-44c9-b26e-496f1144f998
co_a_r_5_siemens = coef(model_a_r_5_siemens)

# ╔═╡ f19f69c3-e484-427a-9010-2d94446d9d2f
md"""
### Visualize
"""

# ╔═╡ 6b0bad7a-bc9d-4ae6-b5f4-46ef4cea9df9
function reprod_a_6_7()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df6_7_large_clean_a_canon[!, :x1], df6_7_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df6_7_medium_clean_a_canon[!, :x1], df6_7_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df6_7_small_clean_a_canon[!, :x1], df6_7_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_5_canon, linestyle=:dashdot)
    if co_a_r_5_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_5_canon[2]; digits=3))x+$(trunc(co_a_r_5_canon[1]; digits=3)) \nr = $(trunc(r_a_r_5_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_5_canon[2]; digits=3))x$(trunc(co_a_r_5_canon[1]; digits=3)) \nr = $(trunc(r_a_r_5_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_7_large_clean_a_geee[!, :x1], df6_7_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df6_7_medium_clean_a_geee[!, :x1], df6_7_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_5_geee, linestyle=:dashdot)
    if co_a_r_5_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_5_geee[2]; digits=3))x+$(trunc(co_a_r_5_geee[1]; digits=3)) \nr = $(trunc(r_a_r_5_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_5_geee[2]; digits=3))x$(trunc(co_a_r_5_geee[1]; digits=3)) \nr = $(trunc(r_a_r_5_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_7_large_clean_a_philips[!, :x1], df6_7_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df6_7_medium_clean_a_philips[!, :x1], df6_7_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_5_philips, linestyle=:dashdot)
    if co_a_r_5_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_5_philips[2]; digits=3))x+$(trunc(co_a_r_5_philips[1]; digits=3)) \nr = $(trunc(r_a_r_5_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_5_philips[2]; digits=3))x$(trunc(co_a_r_5_philips[1]; digits=3)) \nr = $(trunc(r_a_r_5_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_7_large_clean_a_siemens[!, :x1], df6_7_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_7_medium_clean_a_siemens[!, :x1], df6_7_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_5_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_a_r_5_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_5_siemens[2]; digits=3))x+$(trunc(co_a_r_5_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_5_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_5_siemens[2]; digits=3))x$(trunc(co_a_r_5_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_5_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_5_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_5_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Calcium Scoring (6 vs 7)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_a_6_7.png", f)
    f
end

# ╔═╡ a17f7011-315b-4b79-bda7-38f8438688bd
with_theme(medphys_theme) do
    reprod_a_6_7()
end

# ╔═╡ edd93a5d-6003-4dbf-beeb-4a7c352c69e1
md"""
## 6 vs 8
"""

# ╔═╡ 0b3dd8f7-a65b-447d-bb48-dddc7bf929ea
begin
    df6_8_large_clean_a = innerjoin(df6_large_a, df8_large_a; on=[:vendor, :inserts], makeunique=true)
    df6_8_medium_clean_a = innerjoin(df6_medium_a, df8_medium_a; on=[:vendor, :inserts], makeunique=true)
    df6_8_small_clean_a = innerjoin(df6_small_a, df8_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 1beeb7b4-fde9-4e5f-b5ed-740d712edf69
md"""
### Canon
"""

# ╔═╡ 3b16345c-080c-459f-af6c-6fe657b6d142
begin
    df6_8_large_clean_a_canon, df6_8_large_clean_a_geee, df6_8_large_clean_a_philips, df6_8_large_clean_a_siemens = groupby(df6_8_large_clean_a, :vendor)
    df6_8_medium_clean_a_canon, df6_8_medium_clean_a_geee, df6_8_medium_clean_a_philips, df6_8_medium_clean_a_siemens = groupby(df6_8_medium_clean_a, :vendor)
    df6_8_small_clean_a_canon = df6_8_small_clean_a
end;

# ╔═╡ 7df8ff9b-773e-4b58-8f1b-6f8d6461ad41
let
    r_array = [Tuple(df6_8_large_clean_a_canon[!, :x1])..., Tuple(df6_8_medium_clean_a_canon[!, :x1])..., Tuple(df6_8_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_a_canon[!, :x1_1])..., Tuple(df6_8_medium_clean_a_canon[!, :x1_1])..., Tuple(df6_8_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_6_canon
    model_a_r_6_canon = lm(@formula(Y ~ X), data)
    global r_a_r_6_canon
    r_a_r_6_canon = GLM.r2(model_a_r_6_canon)
    global rms_values_a_r_6_canon
    rms_values_a_r_6_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_6_canon))
    ]
end

# ╔═╡ 1602e2fb-801e-4e68-8bd6-24ad5b5e9671
begin
    X_a_r_6_canon = DataFrame(X=collect(1:1000))
    pred_a_r_6_canon = GLM.predict(model_a_r_6_canon, X_a_r_6_canon)
end

# ╔═╡ 2b5eaa63-d745-4875-9359-69c330d5b70d
co_a_r_6_canon = coef(model_a_r_6_canon)

# ╔═╡ 58b6786b-6ca7-4401-ac06-53520782d5e9
md"""
### GE
"""

# ╔═╡ 8e25cad9-4527-4059-8904-9ac13818e6f2
let
    r_array = [Tuple(df6_8_large_clean_a_geee[!, :x1])..., Tuple(df6_8_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_a_geee[!, :x1_1])..., Tuple(df6_8_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_6_geee
    model_a_r_6_geee = lm(@formula(Y ~ X), data)
    global r_a_r_6_geee
    r_a_r_6_geee = GLM.r2(model_a_r_6_geee)
    global rms_values_a_r_6_geee
    rms_values_a_r_6_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_6_geee))
    ]
end

# ╔═╡ da382dd0-0664-48dd-8e69-d83181354b62
begin
    X_a_r_6_geee = DataFrame(X=collect(1:1000))
    pred_a_r_6_geee = GLM.predict(model_a_r_6_geee, X_a_r_6_geee)
end

# ╔═╡ fd5b0002-ffca-42f1-8b17-541aae43b346
co_a_r_6_geee = coef(model_a_r_6_geee)

# ╔═╡ 141b6d46-8d29-4518-b55d-dec1549a8cbb
md"""
### Philips
"""

# ╔═╡ 13aed017-8de3-4b39-8b84-22f3b08fe41c
let
    r_array = [Tuple(df6_8_large_clean_a_philips[!, :x1])..., Tuple(df6_8_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_a_philips[!, :x1_1])..., Tuple(df6_8_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_6_philips
    model_a_r_6_philips = lm(@formula(Y ~ X), data)
    global r_a_r_6_philips
    r_a_r_6_philips = GLM.r2(model_a_r_6_philips)
    global rms_values_a_r_6_philips
    rms_values_a_r_6_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_6_philips))
    ]
end

# ╔═╡ 77531f28-b09e-4aa6-9f83-b2d16df123b7
begin
    X_a_r_6_philips = DataFrame(X=collect(1:1000))
    pred_a_r_6_philips = GLM.predict(model_a_r_6_philips, X_a_r_6_philips)
end

# ╔═╡ accb7d7c-1a8f-47ff-b222-99a6fcc9f162
co_a_r_6_philips = coef(model_a_r_6_philips)

# ╔═╡ 31f1c405-8460-4a76-a097-7d378c978c71
md"""
### Siemens
"""

# ╔═╡ 8771f9f8-c53f-4f9f-a51f-6bd1e6da15d4
let
    r_array = [Tuple(df6_8_large_clean_a_siemens[!, :x1])..., Tuple(df6_8_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df6_8_large_clean_a_siemens[!, :x1_1])..., Tuple(df6_8_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_6_siemens
    model_a_r_6_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_6_siemens
    r_a_r_6_siemens = GLM.r2(model_a_r_6_siemens)
    global rms_values_a_r_6_siemens
    rms_values_a_r_6_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_6_siemens))
    ]
end

# ╔═╡ 1c466b4d-2ea3-4041-b85e-cf5b1ea16d80
begin
    X_a_r_6_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_6_siemens = GLM.predict(model_a_r_6_siemens, X_a_r_6_siemens)
end

# ╔═╡ 46e84370-18e1-498c-a2b7-afdec69edc4a
co_a_r_6_siemens = coef(model_a_r_6_siemens)

# ╔═╡ d65807e2-12a4-4f71-a7c6-3a57495e7007
md"""
### Visualize
"""

# ╔═╡ 920587cd-bd4d-40ac-87c3-c9e28ffb5dea
function reprod_a_6_8()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df6_8_large_clean_a_canon[!, :x1], df6_8_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df6_8_medium_clean_a_canon[!, :x1], df6_8_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df6_8_small_clean_a_canon[!, :x1], df6_8_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_6_canon, linestyle=:dashdot)
    if co_a_r_6_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_6_canon[2]; digits=3))x+$(trunc(co_a_r_6_canon[1]; digits=3)) \nr = $(trunc(r_a_r_6_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_6_canon[2]; digits=3))x$(trunc(co_a_r_6_canon[1]; digits=3)) \nr = $(trunc(r_a_r_6_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_8_large_clean_a_geee[!, :x1], df6_8_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df6_8_medium_clean_a_geee[!, :x1], df6_8_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_6_geee, linestyle=:dashdot)
    if co_a_r_6_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_6_geee[2]; digits=3))x+$(trunc(co_a_r_6_geee[1]; digits=3)) \nr = $(trunc(r_a_r_6_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_6_geee[2]; digits=3))x$(trunc(co_a_r_6_geee[1]; digits=3)) \nr = $(trunc(r_a_r_6_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_8_large_clean_a_philips[!, :x1], df6_8_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df6_8_medium_clean_a_philips[!, :x1], df6_8_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_6_philips, linestyle=:dashdot)
    if co_a_r_6_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_6_philips[2]; digits=3))x+$(trunc(co_a_r_6_philips[1]; digits=3)) \nr = $(trunc(r_a_r_6_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_6_philips[2]; digits=3))x$(trunc(co_a_r_6_philips[1]; digits=3)) \nr = $(trunc(r_a_r_6_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_8_large_clean_a_siemens[!, :x1], df6_8_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_8_medium_clean_a_siemens[!, :x1], df6_8_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_6_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_a_r_6_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_6_siemens[2]; digits=3))x+$(trunc(co_a_r_6_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_6_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_6_siemens[2]; digits=3))x$(trunc(co_a_r_6_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_6_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_6_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_6_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Calcium Scoring (6 vs 8)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ 8e4aacb1-ee4d-4b96-b633-a009ca88ab0e
with_theme(medphys_theme) do
    reprod_a_6_8()
end

# ╔═╡ 57f2a900-aaa7-4bf8-90e9-9abcc83b5e1b
md"""
## 6 vs 9
"""

# ╔═╡ 64cc8f00-b29c-4b6f-b5fd-211e3f7dc438
begin
    df6_9_large_clean_a = innerjoin(df6_large_a, df9_large_a; on=[:vendor, :inserts], makeunique=true)
    df6_9_medium_clean_a = innerjoin(df6_medium_a, df9_medium_a; on=[:vendor, :inserts], makeunique=true)
    df6_9_small_clean_a = innerjoin(df6_small_a, df9_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ dd29fc2d-814a-4fc0-90cb-9447b47acfc2
md"""
### Canon
"""

# ╔═╡ 9b9777aa-c5de-4753-92cd-170c971aa547
begin
    df6_9_large_clean_a_canon, df6_9_large_clean_a_geee, df6_9_large_clean_a_philips, df6_9_large_clean_a_siemens = groupby(df6_9_large_clean_a, :vendor)
    df6_9_medium_clean_a_canon, df6_9_medium_clean_a_geee, df6_9_medium_clean_a_philips, df6_9_medium_clean_a_siemens = groupby(df6_9_medium_clean_a, :vendor)
    df6_9_small_clean_a_canon = df6_9_small_clean_a
end;

# ╔═╡ d890fc82-5ba8-4de6-aeb0-6821e505cd7e
let
    r_array = [Tuple(df6_9_large_clean_a_canon[!, :x1])..., Tuple(df6_9_medium_clean_a_canon[!, :x1])..., Tuple(df6_9_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_a_canon[!, :x1_1])..., Tuple(df6_9_medium_clean_a_canon[!, :x1_1])..., Tuple(df6_9_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_7_canon
    model_a_r_7_canon = lm(@formula(Y ~ X), data)
    global r_a_r_7_canon
    r_a_r_7_canon = GLM.r2(model_a_r_7_canon)
    global rms_values_a_r_7_canon
    rms_values_a_r_7_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_7_canon))
    ]
end

# ╔═╡ 70a98bb9-889c-4732-b006-4c98e1e8a89b
begin
    X_a_r_7_canon = DataFrame(X=collect(1:1000))
    pred_a_r_7_canon = GLM.predict(model_a_r_7_canon, X_a_r_7_canon)
end

# ╔═╡ 935d4fb6-26c3-44d5-9155-0f066af02824
co_a_r_7_canon = coef(model_a_r_7_canon)

# ╔═╡ 2ec88f23-5876-4760-9b3c-2dab9b232955
md"""
### GE
"""

# ╔═╡ b3165e34-fdf7-4253-ab05-bcfcddf9266a
let
    r_array = [Tuple(df6_9_large_clean_a_geee[!, :x1])..., Tuple(df6_9_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_a_geee[!, :x1_1])..., Tuple(df6_9_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_7_geee
    model_a_r_7_geee = lm(@formula(Y ~ X), data)
    global r_a_r_7_geee
    r_a_r_7_geee = GLM.r2(model_a_r_7_geee)
    global rms_values_a_r_7_geee
    rms_values_a_r_7_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_7_geee))
    ]
end

# ╔═╡ 87d2466e-6c06-4261-aced-bf879adeb2ce
begin
    X_a_r_7_geee = DataFrame(X=collect(1:1000))
    pred_a_r_7_geee = GLM.predict(model_a_r_7_geee, X_a_r_7_geee)
end

# ╔═╡ c5f98e9f-88f5-47be-a90d-afcaea1fdede
co_a_r_7_geee = coef(model_a_r_7_geee)

# ╔═╡ 81322f55-2bcf-4b53-9b69-c7a690c304e9
md"""
### Philips
"""

# ╔═╡ 1fa9b627-120e-4169-8b1d-ae1f259f3b19
let
    r_array = [Tuple(df6_9_large_clean_a_philips[!, :x1])..., Tuple(df6_9_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_a_philips[!, :x1_1])..., Tuple(df6_9_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_7_philips
    model_a_r_7_philips = lm(@formula(Y ~ X), data)
    global r_a_r_7_philips
    r_a_r_7_philips = GLM.r2(model_a_r_7_philips)
    global rms_values_a_r_7_philips
    rms_values_a_r_7_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_7_philips))
    ]
end

# ╔═╡ f09115ca-10d2-45f7-8d33-808233f12664
begin
    X_a_r_7_philips = DataFrame(X=collect(1:1000))
    pred_a_r_7_philips = GLM.predict(model_a_r_7_philips, X_a_r_7_philips)
end

# ╔═╡ b4f7e3d3-adbd-437d-b2e3-4eadf4c1eb02
co_a_r_7_philips = coef(model_a_r_7_philips)

# ╔═╡ b1753d95-3c8c-40e5-a18e-836eb7958ae4
md"""
### Siemens
"""

# ╔═╡ 0e63a78e-e2b7-4343-95fc-2762084dd23f
let
    r_array = [Tuple(df6_9_large_clean_a_siemens[!, :x1])..., Tuple(df6_9_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df6_9_large_clean_a_siemens[!, :x1_1])..., Tuple(df6_9_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_7_siemens
    model_a_r_7_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_7_siemens
    r_a_r_7_siemens = GLM.r2(model_a_r_7_siemens)
    global rms_values_a_r_7_siemens
    rms_values_a_r_7_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_7_siemens))
    ]
end

# ╔═╡ dbf7c73f-853c-415e-9c93-708227eaa45d
begin
    X_a_r_7_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_7_siemens = GLM.predict(model_a_r_7_siemens, X_a_r_7_siemens)
end

# ╔═╡ d419f6a1-2823-4496-a791-38d71f8f5e69
co_a_r_7_siemens = coef(model_a_r_7_siemens)

# ╔═╡ 2af1562d-cecd-4ab4-b4eb-b1aa624872b4
md"""
### Visualize
"""

# ╔═╡ a063450f-57f7-4ea1-b2af-8a96b21ae0b7
function reprod_a_6_9()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df6_9_large_clean_a_canon[!, :x1], df6_9_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df6_9_medium_clean_a_canon[!, :x1], df6_9_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df6_9_small_clean_a_canon[!, :x1], df6_9_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_7_canon, linestyle=:dashdot)
    if co_a_r_7_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_7_canon[2]; digits=3))x+$(trunc(co_a_r_7_canon[1]; digits=3)) \nr = $(trunc(r_a_r_7_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_7_canon[2]; digits=3))x$(trunc(co_a_r_7_canon[1]; digits=3)) \nr = $(trunc(r_a_r_7_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_9_large_clean_a_geee[!, :x1], df6_9_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df6_9_medium_clean_a_geee[!, :x1], df6_9_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_7_geee, linestyle=:dashdot)
    if co_a_r_7_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_7_geee[2]; digits=3))x+$(trunc(co_a_r_7_geee[1]; digits=3)) \nr = $(trunc(r_a_r_7_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_7_geee[2]; digits=3))x$(trunc(co_a_r_7_geee[1]; digits=3)) \nr = $(trunc(r_a_r_7_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_9_large_clean_a_philips[!, :x1], df6_9_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df6_9_medium_clean_a_philips[!, :x1], df6_9_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_7_philips, linestyle=:dashdot)
    if co_a_r_7_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_7_philips[2]; digits=3))x+$(trunc(co_a_r_7_philips[1]; digits=3)) \nr = $(trunc(r_a_r_7_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_7_philips[2]; digits=3))x$(trunc(co_a_r_7_philips[1]; digits=3)) \nr = $(trunc(r_a_r_7_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_9_large_clean_a_siemens[!, :x1], df6_9_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_9_medium_clean_a_siemens[!, :x1], df6_9_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_7_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_a_r_7_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_7_siemens[2]; digits=3))x+$(trunc(co_a_r_7_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_7_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_7_siemens[2]; digits=3))x$(trunc(co_a_r_7_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_7_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_7_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_7_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Calcium Scoring (6 vs 9)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end

    # save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-phantom/figures/repro_i_large_phantom.png", f)
    f
end

# ╔═╡ f619ec41-7354-478c-bce1-bf61d738e6f9
with_theme(medphys_theme) do
    reprod_a_6_9()
end

# ╔═╡ db9b5e37-def5-4506-82a7-b71250348b47
md"""
## 6 vs 10
"""

# ╔═╡ 4bcf7744-af84-4584-aca7-6e3a96b7c1ed
begin
    df6_10_large_clean_a = innerjoin(df6_large_a, df10_large_a; on=[:vendor, :inserts], makeunique=true)
    df6_10_medium_clean_a = innerjoin(df6_medium_a, df10_medium_a; on=[:vendor, :inserts], makeunique=true)
    df6_10_small_clean_a = innerjoin(df6_small_a, df10_small_a; on=[:vendor, :inserts], makeunique=true)
end;

# ╔═╡ 8a57c7ae-bb89-43df-b870-f44c759bedd7
md"""
### Canon
"""

# ╔═╡ d7e7ba7d-0bfc-4ee7-afb2-e73afb0e27dc
begin
    df6_10_large_clean_a_canon, df6_10_large_clean_a_geee, df6_10_large_clean_a_philips, df6_10_large_clean_a_siemens = groupby(df6_10_large_clean_a, :vendor)
    df6_10_medium_clean_a_canon, df6_10_medium_clean_a_geee, df6_10_medium_clean_a_philips, df6_10_medium_clean_a_siemens = groupby(df6_10_medium_clean_a, :vendor)
    df6_10_small_clean_a_canon = df6_10_small_clean_a
end;

# ╔═╡ 2a5f0d5f-661f-4076-a33c-8d4143d47b00
let
    r_array = [Tuple(df6_10_large_clean_a_canon[!, :x1])..., Tuple(df6_10_medium_clean_a_canon[!, :x1])..., Tuple(df6_10_small_clean_a_canon[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_a_canon[!, :x1_1])..., Tuple(df6_10_medium_clean_a_canon[!, :x1_1])..., Tuple(df6_10_small_clean_a_canon[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_8_canon
    model_a_r_8_canon = lm(@formula(Y ~ X), data)
    global r_a_r_8_canon
    r_a_r_8_canon = GLM.r2(model_a_r_8_canon)
    global rms_values_a_r_8_canon
    rms_values_a_r_8_canon = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_8_canon))
    ]
end

# ╔═╡ 20e5830b-f10d-4a91-ba28-44504865ff87
begin
    X_a_r_8_canon = DataFrame(X=collect(1:1000))
    pred_a_r_8_canon = GLM.predict(model_a_r_8_canon, X_a_r_8_canon)
end

# ╔═╡ 13944fee-0eb2-416d-9e20-0a1701469df5
co_a_r_8_canon = coef(model_a_r_8_canon)

# ╔═╡ 912ba6aa-9838-46b5-a0c9-308b507f9491
md"""
### GE
"""

# ╔═╡ 63944386-35b1-41b2-8d09-4353121d0ebb
let
    r_array = [Tuple(df6_10_large_clean_a_geee[!, :x1])..., Tuple(df6_10_medium_clean_a_geee[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_a_geee[!, :x1_1])..., Tuple(df6_10_medium_clean_a_geee[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_8_geee
    model_a_r_8_geee = lm(@formula(Y ~ X), data)
    global r_a_r_8_geee
    r_a_r_8_geee = GLM.r2(model_a_r_8_geee)
    global rms_values_a_r_8_geee
    rms_values_a_r_8_geee = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_8_geee))
    ]
end

# ╔═╡ b83cfe12-765e-40a0-9cea-7eabaa859dc5
begin
    X_a_r_8_geee = DataFrame(X=collect(1:1000))
    pred_a_r_8_geee = GLM.predict(model_a_r_8_geee, X_a_r_8_geee)
end

# ╔═╡ a8fb3427-5dc5-472e-aa6d-f6e76033e119
co_a_r_8_geee = coef(model_a_r_8_geee)

# ╔═╡ bf01f558-c98e-490c-a00f-9f039ec8fb86
md"""
### Philips
"""

# ╔═╡ 26416d72-cd1c-4ba0-9544-b7336373472e
let
    r_array = [Tuple(df6_10_large_clean_a_philips[!, :x1])..., Tuple(df6_10_medium_clean_a_philips[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_a_philips[!, :x1_1])..., Tuple(df6_10_medium_clean_a_philips[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_8_philips
    model_a_r_8_philips = lm(@formula(Y ~ X), data)
    global r_a_r_8_philips
    r_a_r_8_philips = GLM.r2(model_a_r_8_philips)
    global rms_values_a_r_8_philips
    rms_values_a_r_8_philips = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_8_philips))
    ]
end

# ╔═╡ 3ee0665c-9ff4-439a-8643-658c9a4fcf91
begin
    X_a_r_8_philips = DataFrame(X=collect(1:1000))
    pred_a_r_8_philips = GLM.predict(model_a_r_8_philips, X_a_r_8_philips)
end

# ╔═╡ 5c8c33a4-62a0-4a45-b2f1-ffd0a0398e55
co_a_r_8_philips = coef(model_a_r_8_philips)

# ╔═╡ c13d23e3-b261-4abc-ad09-51af74c8bff7
md"""
### Siemens
"""

# ╔═╡ 6f1965b2-c5cc-4e31-b6cb-644574c0841a
let
    r_array = [Tuple(df6_10_large_clean_a_siemens[!, :x1])..., Tuple(df6_10_medium_clean_a_siemens[!, :x1])...]
    calc_array = [Tuple(df6_10_large_clean_a_siemens[!, :x1_1])..., Tuple(df6_10_medium_clean_a_siemens[!, :x1_1])...]
    data = DataFrame(
        X=r_array,
        Y=calc_array
    )
    global model_a_r_8_siemens
    model_a_r_8_siemens = lm(@formula(Y ~ X), data)
    global r_a_r_8_siemens
    r_a_r_8_siemens = GLM.r2(model_a_r_8_siemens)
    global rms_values_a_r_8_siemens
    rms_values_a_r_8_siemens = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model_a_r_8_siemens))
    ]
end

# ╔═╡ 426b0ca3-a308-49ab-98c4-b5940848d62f
begin
    X_a_r_8_siemens = DataFrame(X=collect(1:1000))
    pred_a_r_8_siemens = GLM.predict(model_a_r_8_siemens, X_a_r_8_siemens)
end

# ╔═╡ 3c6e1f92-0bee-4318-b212-235ac2f6bda0
co_a_r_8_siemens = coef(model_a_r_8_siemens)

# ╔═╡ 0fd5a508-bbfb-4887-b8e9-1a609122a8cf
md"""
### Visualize
"""

# ╔═╡ 4b021012-2f4c-4b60-986a-137ebc8984cb
function reprod_a_6_10()
    f = Figure()

    ##-- A --##
    ax1 = Axis(f[1, 1])
    sc1 = scatter!(ax1, df6_10_large_clean_a_canon[!, :x1], df6_10_large_clean_a_canon[!, :x1_1])
    sc2 = scatter!(ax1, df6_10_medium_clean_a_canon[!, :x1], df6_10_medium_clean_a_canon[!, :x1_1])
    sc3 = scatter!(ax1, df6_10_small_clean_a_canon[!, :x1], df6_10_small_clean_a_canon[!, :x1_1])
    lines!(ax1, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax1, collect(1:1000), pred_a_r_8_canon, linestyle=:dashdot)
    if co_a_r_8_canon[1] > 0
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_8_canon[2]; digits=3))x+$(trunc(co_a_r_8_canon[1]; digits=3)) \nr = $(trunc(r_a_r_8_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 1],
            placeholder="y = $(trunc(co_a_r_8_canon[2]; digits=3))x$(trunc(co_a_r_8_canon[1]; digits=3)) \nr = $(trunc(r_a_r_8_canon; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_canon[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_canon[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax1, low=0, high=125)
    ylims!(ax1, low=0, high=125)
    ax1.xticks = [0, 25, 50, 75, 100, 125]
    ax1.yticks = [0, 25, 50, 75, 100, 125]
    ax1.xlabel = "Mass 1 (mg)"
    ax1.ylabel = "Mass 2 (mg)"
    ax1.title = "(Canon)"

    ##-- B --##
    ax2 = Axis(f[2, 1])
    scatter!(ax2, df6_10_large_clean_a_geee[!, :x1], df6_10_large_clean_a_geee[!, :x1_1])
    scatter!(ax2, df6_10_medium_clean_a_geee[!, :x1], df6_10_medium_clean_a_geee[!, :x1_1])
    lines!(ax2, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax2, collect(1:1000), pred_a_r_8_geee, linestyle=:dashdot)
    if co_a_r_8_geee[1] > 0
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_8_geee[2]; digits=3))x+$(trunc(co_a_r_8_geee[1]; digits=3)) \nr = $(trunc(r_a_r_8_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 1],
            placeholder="y = $(trunc(co_a_r_8_geee[2]; digits=3))x$(trunc(co_a_r_8_geee[1]; digits=3)) \nr = $(trunc(r_a_r_8_geee; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_geee[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_geee[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax2, low=0, high=125)
    ylims!(ax2, low=0, high=125)
    ax2.xticks = [0, 25, 50, 75, 100, 125]
    ax2.yticks = [0, 25, 50, 75, 100, 125]
    ax2.xlabel = "Mass 1 (mg)"
    ax2.ylabel = "Mass 2 (mg)"
    ax2.title = "(GE)"

    ##-- C --##
    ax3 = Axis(f[1, 2])
    scatter!(ax3, df6_10_large_clean_a_philips[!, :x1], df6_10_large_clean_a_philips[!, :x1_1])
    scatter!(ax3, df6_10_medium_clean_a_philips[!, :x1], df6_10_medium_clean_a_philips[!, :x1_1])
    lines!(ax3, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax3, collect(1:1000), pred_a_r_8_philips, linestyle=:dashdot)
    if co_a_r_8_philips[1] > 0
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_8_philips[2]; digits=3))x+$(trunc(co_a_r_8_philips[1]; digits=3)) \nr = $(trunc(r_a_r_8_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[1, 2],
            placeholder="y = $(trunc(co_a_r_8_philips[2]; digits=3))x$(trunc(co_a_r_8_philips[1]; digits=3)) \nr = $(trunc(r_a_r_8_philips; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_philips[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_philips[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax3, low=0, high=125)
    ylims!(ax3, low=0, high=125)
    ax3.xticks = [0, 25, 50, 75, 100, 125]
    ax3.yticks = [0, 25, 50, 75, 100, 125]
    ax3.xlabel = "Mass 1 (mg)"
    ax3.ylabel = "Mass 2 (mg)"
    ax3.title = "(Philips)"

    ##-- D --##
    ax4 = Axis(f[2, 2])
    scatter!(ax4, df6_10_large_clean_a_siemens[!, :x1], df6_10_large_clean_a_siemens[!, :x1_1], label="Large Inserts")
    scatter!(ax4, df6_10_medium_clean_a_siemens[!, :x1], df6_10_medium_clean_a_siemens[!, :x1_1], label="Medium Inserts")
    lines!(ax4, [-1000, 1000], [-1000, 1000], label="Unity")
    lines!(ax4, collect(1:1000), pred_a_r_8_siemens, linestyle=:dashdot, label="Fitted Line")
    if co_a_r_8_siemens[1] > 0
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_8_siemens[2]; digits=3))x+$(trunc(co_a_r_8_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_8_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        Textbox(
            f[2, 2],
            placeholder="y = $(trunc(co_a_r_8_siemens[2]; digits=3))x$(trunc(co_a_r_8_siemens[1]; digits=3)) \nr = $(trunc(r_a_r_8_siemens; digits=3)) \nRMSE: $(trunc(rms_values_a_r_8_siemens[1]; digits=3)) \nRMSD: $(trunc(rms_values_a_r_8_siemens[2]; digits=3))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end

    xlims!(ax4, low=0, high=125)
    ylims!(ax4, low=0, high=125)
    ax4.xticks = [0, 25, 50, 75, 100, 125]
    ax4.yticks = [0, 25, 50, 75, 100, 125]
    ax4.xlabel = "Mass 1 (mg)"
    ax4.ylabel = "Mass 2 (mg)"
    ax4.title = "(Siemens)"

    ##-- LABELS --##
    f[1:2, 3] = Legend(f, ax4, framevisible=false)
    Label(f[0, 1:2], text="Agatston Calcium Scoring (6 vs 10)", font="Helvetica",
        fontsize=40)

    for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[2, 1], f[1, 2], f[2, 2]])
        Label(layout[1, 1, TopLeft()], label,
            fontsize=25,
            padding=(0, 30, 30, 0),
            halign=:right)
    end
    f
end

# ╔═╡ c5086973-ad46-4920-aad9-0afa70d4ebe9
with_theme(medphys_theme) do
    reprod_a_6_10()
end

# ╔═╡ 61335eab-4677-4ff9-bfe7-379ee6c40008
md"""
# DataFrames
"""

# ╔═╡ e53e8e8a-bb01-4695-997d-ad1eaa1deeae
md"""
## Vendor Specific
"""

# ╔═╡ 21b1b025-3e54-40de-b840-9e9087c4098d
scans = ["1 vs 2", "1 vs 3", "1 vs 4", "1 vs 5", "6 vs 7", "6 vs 8", "6 vs 9", "6 vs 10"]

# ╔═╡ 7529b04a-09d5-4a59-971b-8a501b73e248
begin
    integrated_rmse_canon = [rms_values_i_r_1_canon[1], rms_values_i_r_2_canon[1], rms_values_i_r_3_canon[1], rms_values_i_r_4_canon[1], rms_values_i_r_5_canon[1], rms_values_i_r_6_canon[1], rms_values_i_r_7_canon[1], rms_values_i_r_8_canon[1]]

    agatston_rmse_canon = [rms_values_a_r_1_canon[1], rms_values_a_r_2_canon[1], rms_values_a_r_3_canon[1], rms_values_a_r_4_canon[1], rms_values_a_r_5_canon[1], rms_values_a_r_6_canon[1], rms_values_a_r_7_canon[1], rms_values_a_r_8_canon[1]]

    integrated_rmse_geee = [rms_values_i_r_1_geee[1], rms_values_i_r_2_geee[1], rms_values_i_r_3_geee[1], rms_values_i_r_4_geee[1], rms_values_i_r_5_geee[1], rms_values_i_r_6_geee[1], rms_values_i_r_7_geee[1], rms_values_i_r_8_geee[1]]

    agatston_rmse_geee = [rms_values_a_r_1_geee[1], rms_values_a_r_2_geee[1], rms_values_a_r_3_geee[1], rms_values_a_r_4_geee[1], rms_values_a_r_5_geee[1], rms_values_a_r_6_geee[1], rms_values_a_r_7_geee[1], rms_values_a_r_8_geee[1]]

    integrated_rmse_philips = [rms_values_i_r_1_philips[1], rms_values_i_r_2_philips[1], rms_values_i_r_3_philips[1], rms_values_i_r_4_philips[1], rms_values_i_r_5_philips[1], rms_values_i_r_6_philips[1], rms_values_i_r_7_philips[1], rms_values_i_r_8_philips[1]]

    agatston_rmse_philips = [rms_values_a_r_1_philips[1], rms_values_a_r_2_philips[1], rms_values_a_r_3_philips[1], rms_values_a_r_4_philips[1], rms_values_a_r_5_philips[1], rms_values_a_r_6_philips[1], rms_values_a_r_7_philips[1], rms_values_a_r_8_philips[1]]

    integrated_rmse_siemens = [rms_values_i_r_1_siemens[1], rms_values_i_r_2_siemens[1], rms_values_i_r_3_siemens[1], rms_values_i_r_4_siemens[1], rms_values_i_r_5_siemens[1], rms_values_i_r_6_siemens[1], rms_values_i_r_7_siemens[1], rms_values_i_r_8_siemens[1]]

    agatston_rmse_siemens = [rms_values_a_r_1_siemens[1], rms_values_a_r_2_siemens[1], rms_values_a_r_3_siemens[1], rms_values_a_r_4_siemens[1], rms_values_a_r_5_siemens[1], rms_values_a_r_6_siemens[1], rms_values_a_r_7_siemens[1], rms_values_a_r_8_siemens[1]]
end

# ╔═╡ f0e39fec-2624-4199-bff5-abac5caa4494
begin
    integrated_rmsd_canon = [rms_values_i_r_1_canon[2], rms_values_i_r_2_canon[2], rms_values_i_r_3_canon[2], rms_values_i_r_4_canon[2], rms_values_i_r_5_canon[2], rms_values_i_r_6_canon[2], rms_values_i_r_7_canon[2], rms_values_i_r_8_canon[2]]

    agatston_rmsd_canon = [rms_values_a_r_1_canon[2], rms_values_a_r_2_canon[2], rms_values_a_r_3_canon[2], rms_values_a_r_4_canon[2], rms_values_a_r_5_canon[2], rms_values_a_r_6_canon[2], rms_values_a_r_7_canon[2], rms_values_a_r_8_canon[2]]

    integrated_rmsd_geee = [rms_values_i_r_1_geee[2], rms_values_i_r_2_geee[2], rms_values_i_r_3_geee[2], rms_values_i_r_4_geee[2], rms_values_i_r_5_geee[2], rms_values_i_r_6_geee[2], rms_values_i_r_7_geee[2], rms_values_i_r_8_geee[2]]

    agatston_rmsd_geee = [rms_values_a_r_1_geee[2], rms_values_a_r_2_geee[2], rms_values_a_r_3_geee[2], rms_values_a_r_4_geee[2], rms_values_a_r_5_geee[2], rms_values_a_r_6_geee[2], rms_values_a_r_7_geee[2], rms_values_a_r_8_geee[2]]

    integrated_rmsd_philips = [rms_values_i_r_1_philips[2], rms_values_i_r_2_philips[2], rms_values_i_r_3_philips[2], rms_values_i_r_4_philips[2], rms_values_i_r_5_philips[2], rms_values_i_r_6_philips[2], rms_values_i_r_7_philips[2], rms_values_i_r_8_philips[2]]

    agatston_rmsd_philips = [rms_values_a_r_1_philips[2], rms_values_a_r_2_philips[2], rms_values_a_r_3_philips[2], rms_values_a_r_4_philips[2], rms_values_a_r_5_philips[2], rms_values_a_r_6_philips[2], rms_values_a_r_7_philips[2], rms_values_a_r_8_philips[2]]

    integrated_rmsd_siemens = [rms_values_i_r_1_siemens[2], rms_values_i_r_2_siemens[2], rms_values_i_r_3_siemens[2], rms_values_i_r_4_siemens[2], rms_values_i_r_5_siemens[2], rms_values_i_r_6_siemens[2], rms_values_i_r_7_siemens[2], rms_values_i_r_8_siemens[2]]

    agatston_rmsd_siemens = [rms_values_a_r_1_siemens[2], rms_values_a_r_2_siemens[2], rms_values_a_r_3_siemens[2], rms_values_a_r_4_siemens[2], rms_values_a_r_5_siemens[2], rms_values_a_r_6_siemens[2], rms_values_a_r_7_siemens[2], rms_values_a_r_8_siemens[2]]
end

# ╔═╡ 86e51b7f-5f25-4223-9e39-051f0ef316f0
begin
    integrated_r_corr_canon = [r_i_r_1_canon, r_i_r_2_canon, r_i_r_3_canon, r_i_r_4_canon, r_i_r_5_canon, r_i_r_6_canon, r_i_r_7_canon, r_i_r_8_canon]

    agatston_r_corr_canon = [r_a_r_1_canon, r_a_r_2_canon, r_a_r_3_canon, r_a_r_4_canon, r_a_r_5_canon, r_a_r_6_canon, r_a_r_7_canon, r_a_r_8_canon]

    integrated_r_corr_geee = [r_i_r_1_geee, r_i_r_2_geee, r_i_r_3_geee, r_i_r_4_geee, r_i_r_5_geee, r_i_r_6_geee, r_i_r_7_geee, r_i_r_8_geee]

    agatston_r_corr_geee = [r_a_r_1_geee, r_a_r_2_geee, r_a_r_3_geee, r_a_r_4_geee, r_a_r_5_geee, r_a_r_6_geee, r_a_r_7_geee, r_a_r_8_geee]

    integrated_r_corr_philips = [r_i_r_1_philips, r_i_r_2_philips, r_i_r_3_philips, r_i_r_4_philips, r_i_r_5_philips, r_i_r_6_philips, r_i_r_7_philips, r_i_r_8_philips]

    agatston_r_corr_philips = [r_a_r_1_philips, r_a_r_2_philips, r_a_r_3_philips, r_a_r_4_philips, r_a_r_5_philips, r_a_r_6_philips, r_a_r_7_philips, r_a_r_8_philips]

    integrated_r_corr_siemens = [r_i_r_1_siemens, r_i_r_2_siemens, r_i_r_3_siemens, r_i_r_4_siemens, r_i_r_5_siemens, r_i_r_6_siemens, r_i_r_7_siemens, r_i_r_8_siemens]

    agatston_r_corr_siemens = [r_a_r_1_siemens, r_a_r_2_siemens, r_a_r_3_siemens, r_a_r_4_siemens, r_a_r_5_siemens, r_a_r_6_siemens, r_a_r_7_siemens, r_a_r_8_siemens]
end

# ╔═╡ eb4be2cc-9458-4102-a2e9-6df30ba22045
df_reprod = DataFrame(
    scan=scans,
    integrated_rmse_canon=integrated_rmse_canon,
    agatston_rmse_canon=agatston_rmse_canon,
    integrated_rmse_geee=integrated_rmse_geee,
    agatston_rmse_geee=agatston_rmse_geee,
    integrated_rmse_philips=integrated_rmse_philips,
    agatston_rmse_philips=agatston_rmse_philips,
    integrated_rmse_siemens=integrated_rmse_siemens,
    agatston_rmse_siemens=agatston_rmse_siemens,
    integrated_rmsd_canon=integrated_rmsd_canon,
    agatston_rmsd_canon=agatston_rmsd_canon,
    integrated_rmsd_geee=integrated_rmsd_geee,
    agatston_rmsd_geee=agatston_rmsd_geee,
    integrated_rmsd_philips=integrated_rmsd_philips,
    agatston_rmsd_philips=agatston_rmsd_philips,
    integrated_rmsd_siemens=integrated_rmsd_siemens,
    agatston_rmsd_siemens=agatston_rmsd_siemens,
    integrated_r_corr_canon=integrated_r_corr_canon,
    agatston_r_corr_canon=agatston_r_corr_canon,
    integrated_r_corr_geee=integrated_r_corr_geee,
    agatston_r_corr_geee=agatston_r_corr_geee,
    integrated_r_corr_philips=integrated_r_corr_philips,
    agatston_r_corr_philips=agatston_r_corr_philips,
    integrated_r_corr_siemens=integrated_r_corr_siemens,
    agatston_r_corr_siemens=agatston_r_corr_siemens
)

# ╔═╡ 4e9b8f08-3875-4f55-a9f2-5c0c9edd16ed
size(df_reprod)

# ╔═╡ 73180533-9668-4f83-bc8c-4da7b492b019
md"""
## Average
"""

# ╔═╡ b59fb57e-cec5-4f3d-8808-999f274f0679
begin
    mean_rmse_i_canon = mean(df_reprod[!, :integrated_rmse_canon])
    mean_rmse_i_geee = mean(df_reprod[!, :integrated_rmse_geee])
    mean_rmse_i_philips = mean(df_reprod[!, :integrated_rmse_philips])
    mean_rmse_i_siemens = mean(df_reprod[!, :integrated_rmse_siemens])

    mean_rmse_a_canon = mean(df_reprod[!, :agatston_rmse_canon])
    mean_rmse_a_geee = mean(df_reprod[!, :agatston_rmse_geee])
    mean_rmse_a_philips = mean(df_reprod[!, :agatston_rmse_philips])
    mean_rmse_a_siemens = mean(df_reprod[!, :agatston_rmse_siemens])
end

# ╔═╡ 6eb2a227-ee65-456f-bfbc-0c727dfde3c9
begin
    mean_rmsd_i_canon = mean(df_reprod[!, :integrated_rmsd_canon])
    mean_rmsd_i_geee = mean(df_reprod[!, :integrated_rmsd_geee])
    mean_rmsd_i_philips = mean(df_reprod[!, :integrated_rmsd_philips])
    mean_rmsd_i_siemens = mean(df_reprod[!, :integrated_rmsd_siemens])

    mean_rmsd_a_canon = mean(df_reprod[!, :agatston_rmsd_canon])
    mean_rmsd_a_geee = mean(df_reprod[!, :agatston_rmsd_geee])
    mean_rmsd_a_philips = mean(df_reprod[!, :agatston_rmsd_philips])
    mean_rmsd_a_siemens = mean(df_reprod[!, :agatston_rmsd_siemens])
end

# ╔═╡ 0df13a11-c319-4590-b684-ca0d409a7a29
mean_rmse_integrated = [mean_rmse_i_canon, mean_rmse_i_geee, mean_rmse_i_philips, mean_rmse_i_siemens]

# ╔═╡ 57ba1a2b-0f50-4c71-b105-6f4a0f61410f
mean_rmse_agatston = [mean_rmse_a_canon, mean_rmse_a_geee, mean_rmse_a_philips, mean_rmse_a_siemens]

# ╔═╡ b2224158-0916-4b6e-a343-7887e21aca63
mean_rmsd_integrated = [mean_rmsd_i_canon, mean_rmsd_i_geee, mean_rmsd_i_philips, mean_rmsd_i_siemens]

# ╔═╡ dc337e00-3a48-4be5-bf1d-f5af674c8edc
mean_rmsd_agatston = [mean_rmsd_a_canon, mean_rmsd_a_geee, mean_rmsd_a_philips, mean_rmsd_a_siemens]

# ╔═╡ f41b02a2-0494-4ab1-9bfd-1f9aab535648
df_reprod_avg = DataFrame(
    vendor=["canon", "ge", "philips", "siemens"],
    mean_rmse_integrated=mean_rmse_integrated,
    mean_rmse_agatston=mean_rmse_agatston,
    mean_rmsd_integrated=mean_rmsd_integrated,
    mean_rmsd_agatston=mean_rmsd_agatston
)

# ╔═╡ 95e5116e-bf8e-407a-b356-d0f431bf2dc4
size(df_reprod_avg)

# ╔═╡ Cell order:
# ╠═40bf5e39-0b70-47c1-82d8-3033cd1ed60b
# ╠═47b98f40-4522-4bc4-a5de-c2e99e443143
# ╠═3d7dcd97-4202-4d45-8791-cd22099d1b8c
# ╠═860c3023-55be-4cd5-924a-1a677d6b9281
# ╠═86b24b18-b24c-48b1-b89c-2bde53744fe5
# ╠═2838451c-be95-41f7-813b-7323d4a621b4
# ╠═7bf364d8-cb82-4f79-8c53-0a426eb09d12
# ╠═0520c367-4ecf-4c31-a396-6b70e1a630a5
# ╟─960ad1f4-83a9-4e78-9d2e-1517a75cb127
# ╠═097f3636-6fe6-41d7-ae49-dd761d3cce86
# ╠═e848a449-2983-4ead-8bad-f5a8fa735b69
# ╠═44576c72-f38a-48de-93bc-e4b79ba1a9a3
# ╠═808d5303-7401-43b0-bc1f-42578770e668
# ╠═094f4999-30ab-4f7c-b7bc-accaf6c07673
# ╠═16285410-d2ad-41de-b383-531f0ba1ee96
# ╠═9260e6cc-20e3-477c-851c-aae4a0cd5677
# ╟─bcbca806-78b9-4147-b068-b815201ff53d
# ╠═f7d7d226-1e40-4921-9af0-c065e0a6af30
# ╠═8745cff1-1755-4e2d-9801-21f8af8e0c36
# ╠═1afc780f-a4e0-4526-bf22-49f24eeb90ed
# ╟─7a2f23bb-c548-4410-96c0-6fa7e1630784
# ╟─8f4151a3-3baf-4126-97e7-fafb6a92b2d9
# ╠═2efd3f1c-2d88-4689-ba02-0cf66941357b
# ╟─701f24a7-c333-40e7-a2a9-271fc83aa54b
# ╠═a79fcdb3-a30b-4bda-9bea-f012550ec14d
# ╠═e5478856-fae3-481c-8565-66cc1bbd47ed
# ╠═39a7ac34-4fa1-4cd7-b736-aa6b9c05fb6e
# ╟─3d40a8a2-c8a9-4ba5-a9a7-9bdbba7914b6
# ╠═118ad8fd-4423-4f90-a6ec-5a6a89a40637
# ╠═ec158a7b-2cd3-488c-b192-74f69e34d9b7
# ╠═fac3bf90-967d-4ac8-b223-8910c3e0a423
# ╟─0be9d780-1097-4d23-a745-ffe0b30e1173
# ╠═c0e6bdd0-f780-4af3-aa18-9878c05ab8fd
# ╠═3713e966-9bf1-4577-8b79-0ecd6ed925d1
# ╠═e307f779-b1da-43b6-a505-36e817984cb1
# ╟─0be916a3-5fe5-488e-91d2-81f9267b7b21
# ╠═b4a86e37-8ba7-4f11-b48c-8ee3cf43f981
# ╠═e596f2ed-7826-44b7-8cad-4625c9e40365
# ╠═4bd66f83-427f-46d3-b5e0-7d1c2d2a9cf8
# ╟─f2d3a05e-17d9-4dbc-a399-8ac69dfb37c8
# ╟─d6c75439-e2bf-404e-8641-6820f88d09b1
# ╟─343e5452-3da7-4ed1-b1af-d281ef17192e
# ╟─abc07dfe-477c-431f-92ad-ef51ac76b074
# ╠═e7ad474a-9356-4ea3-a8b6-9d2397d9c4a3
# ╟─1be99b59-792b-40cc-ac8e-7cb2b19a45aa
# ╠═953444d5-ef00-4922-bbb9-9b2c64ab0e92
# ╠═71e0ab27-245f-48dc-9af5-74f922b71d16
# ╠═1a63c8c3-5489-471b-bcd4-4d0c888f12e3
# ╠═8b9cc340-3232-44bc-a2d9-f90279cf829b
# ╟─d5a01db2-c0d5-4488-b7f7-dae3ce95e534
# ╠═82b31f78-bdad-4ac0-ac17-5de1039e24c6
# ╠═afe37b2a-6f33-4813-ae29-2dd6fd5f72d3
# ╠═d31758fb-5df5-4526-9c2d-aca4c9db4048
# ╟─ca3e0281-a373-4f64-a38a-3421c94668ac
# ╠═13c80808-6270-4165-a9d7-289745282c0f
# ╠═c40ab351-ff02-497e-9eda-4d6cb4a58636
# ╠═3648622c-68a9-4de9-8ad2-309b9c08afd7
# ╟─06bd09a5-3d97-488a-8186-d1bd920e9acf
# ╠═2b7fc654-4d1b-465c-a2a8-699f6d2ae373
# ╠═53165e35-6836-4129-b47b-92f9a4caa422
# ╠═560ee5a6-9782-466d-97dd-64196b6fd3d6
# ╠═b6373747-4f01-4a85-a816-8b9a417f27aa
# ╠═a8877540-3144-4068-83c1-f183eb9c5c6e
# ╟─9d677e1c-5754-44ce-beb1-1f73e9c0fffe
# ╠═9a4dd89d-5fd4-4c25-98d5-2bb1d3d85f36
# ╠═f9374f75-36a5-42d6-9607-7d5622b46f4e
# ╠═31b61ad7-ac39-4610-b902-2a2943ebace2
# ╠═52f55f1b-f904-4c83-9992-b5ccb8099ec7
# ╠═61fccccd-dfdf-4705-9356-767d0ab2d72b
# ╠═80b1e347-0550-4002-846b-424c96fbc97e
# ╠═0b649084-457c-4339-bb19-bf3a7200262b
# ╠═12494832-c0fa-45c1-bbd0-05b82543e138
# ╠═16c82c3f-7c99-496d-861e-6da382d1a75f
# ╠═8dc2ac38-8a45-4799-b5d0-a849aa1f7c26
# ╠═5b360dc3-213b-4c01-ab89-b4ddc888eb87
# ╠═6529d5ad-1539-42b2-aa9a-abad76ddf96d
# ╠═356398c4-66e0-4b49-b1ae-2b4815d5e54a
# ╠═fc6eb521-6ddc-4be2-a2d1-685db850e6ee
# ╠═37620e01-de7f-4479-9014-38c0ce5407f6
# ╠═6db477d4-197c-41bc-bf9d-78a1b83c4749
# ╠═6caac89d-80ed-4ecf-8d99-37181f500780
# ╠═1384e914-d7ac-4459-903d-dc8e8f360711
# ╠═39b386be-d0b5-4b9f-8015-450988a5255f
# ╠═b1963053-588f-4137-9ade-a308eedd496a
# ╠═825d85a5-3898-4fb2-8fe3-30c8ba49dd0a
# ╠═f1e9e113-9531-4fb4-95aa-5dfaa146c47e
# ╠═ce1d1375-fe0c-436a-930e-e6e5cc12c172
# ╠═25c01496-5893-40b3-9c55-a7473180ff37
# ╠═9b36dd13-1cc6-480e-a8a3-1e36413e900a
# ╠═ca50b8dd-fe0b-408a-9989-948d07a09ffc
# ╠═430bfe52-0ebd-4de2-9a02-cc09fb4edcc0
# ╠═881771b5-373e-4ff3-9b01-405ae56677e5
# ╠═31f3cbfb-e9c2-46dc-ad0c-a637783b421c
# ╠═480296d4-890f-4484-b0e7-8828dd51df26
# ╠═7bea8ddf-0cb5-4bd7-86c1-24b6e91befac
# ╠═b5e12cdf-3ecb-4d9b-92cd-b76f9f63a5ee
# ╠═bc10bdfc-fa49-4b36-9020-bacf6ddd7940
# ╠═daa19304-8d3c-4e6b-9d60-0b5f064e5643
# ╠═39dc0b0d-834a-43ea-a3b5-45d603bc4413
# ╠═025f813a-1a26-475a-994c-e7c7afcd9d15
# ╠═acf5f08c-edb9-4da4-8c08-672f780cc41b
# ╠═bceb632d-c1b1-4c19-90fd-ea1d97c1db9c
# ╠═d85f6f5c-3160-492d-9aaa-80d19d7826fa
# ╠═917e1ee3-2f55-466e-aad6-16eb0241a2ff
# ╠═4957ff4a-4457-4bd7-9b54-e86561165788
# ╠═dc57ba01-1731-4136-964a-4f4da385a967
# ╠═c2199225-23a1-48da-b58c-fe03994dcd38
# ╠═86adbd4f-b55d-4056-a9f7-aafaebca9b3c
# ╠═60602785-9402-49fc-a898-f8538714feee
# ╠═c0c47a1b-7f84-48bd-b5eb-b80302d81514
# ╟─a7363057-b7bd-410a-95fa-2b70c3a52f16
# ╟─a1f8e8c6-211d-42a3-82aa-32a8560752f3
# ╠═c33449f5-faba-4e86-942d-460773849ecc
# ╟─d78573c3-2ea8-49a8-81d4-077343d3a274
# ╠═e60e3a6d-8815-433c-af18-1b8938e0afc6
# ╠═250cd449-8bba-4854-8197-8f1d3268ee41
# ╠═a9440691-ac39-4c64-b65d-83d746ebc0da
# ╠═2bdb6b1c-e4a8-49f2-bb97-022f553373d7
# ╠═7adf35a0-f75b-4075-8472-425b400998c7
# ╠═8a8bbef8-6806-47a2-8d8e-5e634fb45e0d
# ╠═b568c1e3-e332-4713-8295-15e25b0fee83
# ╠═d8a1f7eb-6d0d-4e18-9674-d9aa169f09e2
# ╠═7a3ca058-d58f-4a56-901e-7b6f7547836a
# ╠═1f3c18b9-f001-4657-ab87-1fef4782f752
# ╠═60b064f0-0341-4c9d-a0e5-50900b5e502b
# ╠═e209caf6-f439-4ada-8b60-51881c572436
# ╠═c50609b3-ac5e-43ac-93bc-cd058eb05d8c
# ╠═d4c78b5b-7501-4513-81a4-06037e190292
# ╠═6abaedc2-e1e9-4e4b-8c96-8f56d15a46ae
# ╠═f972ea2c-2bdc-4caa-80a1-8aff4456b3c1
# ╟─0a42ad4a-5ad9-4bef-a5b9-3153b773957e
# ╠═95f7ddba-c5e7-4725-aff4-9b770a2b0fc0
# ╟─1c68825f-118d-4c38-bbd3-1f6b8688e8a2
# ╠═e830479f-7f74-459a-a531-5fc7dc6122a1
# ╠═ff76ef96-57e4-49f2-9776-e2aee3f86c91
# ╠═ef3eb4ca-3ef8-4700-9efd-6e13e6069faf
# ╠═091fd9b6-ce40-437b-aa80-cf82fc3ed676
# ╠═ce214147-4f83-48d0-b69e-69142214b3d2
# ╠═186994e9-6624-48e3-b383-f67e659ebfcc
# ╠═9fbb0570-3d7e-4dc8-9919-0f70349980ef
# ╠═7f2edfc7-9c6c-4edd-bd0a-28b947ff9c91
# ╠═73c43529-d5c3-4725-a7e5-5c4dea85ebe5
# ╠═7a82f997-dbf3-40d9-b508-eab38c96d0c7
# ╠═44519f33-ce8c-4ced-aba5-b7f8f52f3681
# ╠═2d3eec4e-7cbb-4a0e-ab16-85597af04b02
# ╠═b62a9ebe-7235-4235-bef7-52f259820f75
# ╠═7911378c-f43a-4b1b-a93c-e805cea9fcc6
# ╠═735fe265-976c-45b3-b132-b4a9dc8c5e7c
# ╠═15db30dd-d042-4433-9baf-1d262257b09b
# ╠═881ece53-f462-49e3-a597-21f1154e73ac
# ╠═a9144e10-227f-44df-9a73-1bff8597e6ff
# ╠═709e0630-1d70-4e8b-a550-b760dad49788
# ╠═ae1562ca-a019-4de8-9d17-dbee895a610d
# ╠═a7eae4c1-8fd9-42c3-8653-7572060bae2f
# ╠═a2119a75-fa6b-4b36-b70a-942856a1fe6a
# ╠═1f07992d-747c-430c-a53c-e44df17f6934
# ╠═639c2606-c2b7-4130-a005-ffa314f5096d
# ╠═8ca49384-3739-4b0c-beba-10697655a9a7
# ╠═0b30bdbd-a520-43d2-9c24-121f8912809d
# ╠═d3e16b7d-a15f-471e-9080-e8039bd2760f
# ╠═5a42bf7b-4949-41fa-99b0-1ee8ba304f74
# ╠═157b5e5f-decb-4c90-9c3d-e044890da547
# ╠═86fbe863-7404-488d-afb1-d485d97c131f
# ╠═a5dc5b63-8042-4fed-a4fd-199a337cc120
# ╠═a8fb5fea-58ce-459a-8433-822a9f2fcc50
# ╠═51ec9457-a8a2-4d5b-96d5-d8856af3ce48
# ╠═efdb711e-1bcd-4566-b49b-e302e1536e2a
# ╠═a071e559-ed7c-4001-8c94-a73b3e8f5df9
# ╠═0ae71358-9fe6-4484-bb52-fd4df9e1f4e3
# ╠═89fdec4f-d819-4199-b6d0-3753c16e4f47
# ╠═616a4888-2e30-4fb8-9bda-fa9cb6dead52
# ╠═e4d6734a-7b3d-4932-bef3-e7e7bd3b1d11
# ╠═8a819052-9c66-450c-b4a9-31946883c95b
# ╠═aa54f20e-a59d-4194-8a1c-46c3d1e53f21
# ╠═ef07df25-3ed3-4d7e-8b61-e88a4ad520b6
# ╠═0c5cd43c-118a-45c9-b4e3-e5aae61db595
# ╠═e0d8d102-d56b-4ce0-be63-343c2fa5202a
# ╠═aeb5b3a5-45ce-42e3-b5a5-7c1ebd570324
# ╠═ce698ae4-597b-4fe1-9364-f913ccbbd353
# ╠═07b17ba4-f139-4c00-b024-96b75a349b09
# ╠═184d730e-1197-4466-b95d-8a2ec1a5c4ea
# ╠═38d0b5d1-974a-45a7-9780-70fa45c4dafe
# ╠═5d12a5f4-ac7f-4242-b7ec-eb386c54eebf
# ╠═3247670e-ae2f-45b1-889e-0fe968ea3fda
# ╠═d8992eef-e000-48fa-8cc3-28d7b30c1527
# ╠═1e170158-1bdc-4bae-890e-e4d1d80d0ad7
# ╠═34b78437-0868-4634-a261-40bc2f120edc
# ╠═6f252fc5-6af7-4a11-b472-e2d7211a3837
# ╠═6d21d3d0-62e1-4d7a-915f-c56372e63ef7
# ╠═a79c64ec-5fe3-4db8-91a2-58264d73bb9a
# ╠═466644b0-7c54-457a-98ce-15945d6cfacb
# ╠═b7b03903-53a0-44ec-a522-ad3683bfb0a1
# ╠═3398bb56-7147-4c16-831c-647a917a690b
# ╠═457e7344-2bdb-4b86-965f-adfcf07b28b1
# ╠═e0b819ed-7240-4a77-9a49-46bac20f884f
# ╠═c290b1b4-6e02-42ec-b880-1619ab746391
# ╠═0c5c047e-2c4a-44a1-b0e5-cd6b97a58e56
# ╠═ad895815-ca2c-4167-bb73-ed8789ec7e6c
# ╠═30844636-edd0-4b70-b0c8-7fb3b1d1e71f
# ╠═ffa9a017-24b2-4f7d-b9b0-792ed0378090
# ╠═033f9b89-4f88-42bc-9b14-92f7906808ec
# ╠═e6385cf2-2c54-450c-85b2-f444cde4af6f
# ╠═64e7bd22-bc6e-4920-8c80-46ab542ed38d
# ╠═2993e116-34c2-42cf-8bf4-3aa9919e55b3
# ╠═7b1292b6-5f17-4f9b-91c1-22442087834d
# ╠═00cdfa1b-8e1f-4b43-8630-3f7718ff7ec8
# ╠═699fbbdf-0729-4372-8617-d909f2c663eb
# ╠═4577afea-26c0-4be0-98ea-e0b5095b5baa
# ╠═3a116456-4c27-4ce1-b870-7c4d06f34451
# ╠═6a067a65-61a8-483f-9686-aa91645014eb
# ╠═fff52d0b-2956-46d6-b3a4-fa8045379b65
# ╠═061e0e6e-730f-45b4-bf00-5f8b91f0c507
# ╠═49cdf7dd-d6ac-499a-937b-7ec2e05b04bb
# ╠═3bccc16c-d42f-4784-9e63-47a3d3e42d58
# ╠═081bc4e5-ad09-4749-b99f-3be537efe158
# ╠═1d895d11-02ed-4838-abe5-8bca36fd4560
# ╠═1357670e-a30b-45ff-8d29-9bd07000ed8a
# ╠═65f70586-4d48-4c59-b3ec-abf7b96d6031
# ╠═3d1a536b-218e-4d83-847c-ce088b7bb356
# ╠═7992f0a2-9012-4801-821c-e917df0e06db
# ╠═4bd87cfe-f24a-4f14-8d29-094fab2ee6e8
# ╠═97489c7d-24cb-4685-bf41-76c9457f4dcc
# ╠═4c7c637a-15ab-4103-9f43-5183ca8f8d0f
# ╠═db0370ca-dde6-46e3-aa47-6669b4abe3eb
# ╠═f420b178-6daf-447c-be6e-72b243c15d0f
# ╠═c74f97ce-182d-47c1-91bb-b4b56554c5b1
# ╠═84345749-fbac-4682-a9bf-cb9118686938
# ╠═81da47ef-4d90-444c-bd50-defdbd6e3440
# ╠═30d688ea-bf2f-42b2-aa73-0fd23f7936c7
# ╠═0f011318-bbfa-4ae4-8998-be028c481941
# ╠═005830de-16c5-4aac-bf0b-60b2cc59d735
# ╠═a4eb6ff4-182f-4dd3-bdba-0c46dc032a28
# ╠═673cb1d6-7a77-4c75-851b-ac134270d5d3
# ╠═7c65b05a-eaa4-4c14-937b-696b395a7548
# ╠═781e21a4-bfed-4b0c-b236-6b51b29ab12d
# ╠═42e07123-1010-46c8-999c-525646539a07
# ╠═492d4cc2-0239-49d9-8ac5-4307c0ebb1ee
# ╠═246eadad-ec40-45f8-9808-f2413184be6e
# ╠═75452f19-00ae-4c20-98ba-86708ab301e7
# ╠═1f56b032-e5e5-47e4-a491-87f1149b236a
# ╠═0c7f45c5-0b30-461a-a382-3f9fe7f5ab16
# ╠═b22547d3-3c10-4e45-b108-4a1a11119da1
# ╠═c6bc875b-6849-4286-b40b-6c23e0d78721
# ╠═882745b5-96eb-4cbd-b356-6fd7a2983697
# ╠═fd28b420-caa3-4638-91ff-2e98b5c7fca5
# ╠═73236865-d3b2-4ec4-b747-131ef0aec338
# ╠═36c19466-f340-4549-a531-503334fae7aa
# ╠═58c01c12-d1ec-45cd-a6ad-b3653b4ed3b0
# ╠═b75d526e-f282-49cc-85be-fac4744638eb
# ╠═5b3e9035-9cfb-4175-a781-ad7db7260d40
# ╠═7ee2f5f8-01bb-462f-90d2-960352a53a36
# ╠═98676111-76eb-48df-bc25-c9149a3dd833
# ╠═525d755c-e1fa-4be9-a220-412431900416
# ╠═0fb8547e-b396-4c7a-8809-f7cc7441900a
# ╠═46dd427a-345d-4487-b24f-31b6917a00d2
# ╠═447725c6-ab92-4cc3-9a78-2d3022f1258a
# ╠═da5ccfbd-2ec8-4f55-924a-df687711e305
# ╠═a473c7f9-9885-42e6-a62e-f533b7ad4f24
# ╠═6d973a09-533b-4a40-bfa1-e72e65237082
# ╠═22196e29-f3af-4727-a729-ecac3267af90
# ╠═9be70624-cbff-4083-bfe5-1936a1d08c03
# ╠═e3c253ec-d7eb-4c3b-8f67-73620228f994
# ╠═d409efbf-2a66-4385-9d21-8111722832a0
# ╠═a28122c7-df21-48ab-993b-dd1c7bb966a3
# ╠═afd1a63a-4e61-4451-8d5a-031d833d3ff2
# ╠═0513263c-720f-4e7b-842a-ba65f99b9900
# ╠═5be7b2db-1b1d-4633-8c2f-52549e083fa9
# ╠═78f81079-f267-42a9-b59f-5f75eed991e8
# ╠═23ba4b57-d81c-4f7e-b7fc-7f3af1f5e121
# ╠═dda5c63e-42d3-4192-89fc-8730ac168842
# ╠═114bf88c-4daa-4c28-a0e1-4296cb1257dc
# ╠═646ff2dc-d490-4cfd-8618-b35982fec465
# ╠═9461e6d1-0ce9-4bd9-a68c-ece5fcebc6a5
# ╠═7d1f6a58-da38-4ff4-b0f1-822d62b71f10
# ╠═dc9784f5-88ca-4c23-ad35-648979f35814
# ╠═7ab94d89-9cee-4328-a35c-b519304ad5b8
# ╠═7b8ea451-02a2-4369-a2bf-37c98221d400
# ╠═28d8585e-6f13-4ab8-8e5b-3fcee608c7fd
# ╠═d66bac7a-4952-4c7e-ad1f-053c597a3cb1
# ╠═abab28f2-bb87-46ac-9470-986ca8de668f
# ╠═61429df4-ac03-40ef-91d8-8658a9d080bf
# ╠═d14099a9-7559-4940-b275-29c82a341068
# ╠═624d2696-eb16-434c-a8fa-93e4ec0f680c
# ╠═4fd3ce19-c1ee-452f-b098-901ba46affe5
# ╠═f8cd4e74-72ab-42e5-ad63-0f947fcb9c72
# ╠═67eda3a2-5b52-4927-ac9a-9c867e1b0c75
# ╠═c0e61fb4-6e71-4b17-8444-99eb39f0db20
# ╠═5878377f-b0b2-4e9d-97b4-379b141875ff
# ╠═337eb752-7651-4aa4-b883-80ecaae095aa
# ╠═d8fc55fe-b4c5-4662-9b66-29cca9b9f5c6
# ╠═c495aecb-b875-40d8-9962-ae618ec0042a
# ╠═418f0615-3f3b-4188-a2f9-98a9b16bf81c
# ╠═b52e1a6d-853d-4b31-bb9a-127226f7946e
# ╠═b244bd34-b91a-47f9-bfe4-958b4e3ad635
# ╠═6dd95339-fe3d-4efa-8f8f-c59edd4e435a
# ╠═cc3fb3ee-283d-4b5b-b102-587bd0c40acf
# ╠═556e18bb-16f3-4571-af38-ca56ba90d601
# ╠═b4c5ca8e-d1ab-41f5-aec2-9d6bc3138557
# ╠═69182167-1cc9-4414-9402-4586fa694f4b
# ╠═e036ae91-541f-45d4-9b7b-9d79f72de86f
# ╠═d557365b-3700-4358-9a26-8ea131a07469
# ╠═5ef9cbe2-2852-483e-bdd4-62b1bdda6193
# ╠═03d45881-14cf-4110-9512-12c743c160db
# ╠═216c540f-2d9a-499e-b29f-5c84d3a7e22a
# ╠═91545459-fa68-4b86-b99d-b80e2db7ea74
# ╠═a4546417-20b2-4e1d-ab33-671d89610fdb
# ╠═cebe4f0b-5fc0-49fd-9ce8-4fe024097142
# ╠═1668c2c6-b278-419e-8eae-bf6ee8da7b6a
# ╠═b10ade70-d733-4338-8fe1-090128ec172f
# ╠═b135ec0b-dec9-442b-9ccf-63122eca6cfa
# ╠═f43b75b1-c6c2-4530-a73f-02b1151bd0b6
# ╠═a212e5a5-631f-49fa-88b1-5a3af605d1d2
# ╠═8822f9e1-9829-44c9-b26e-496f1144f998
# ╠═f19f69c3-e484-427a-9010-2d94446d9d2f
# ╠═6b0bad7a-bc9d-4ae6-b5f4-46ef4cea9df9
# ╠═a17f7011-315b-4b79-bda7-38f8438688bd
# ╠═edd93a5d-6003-4dbf-beeb-4a7c352c69e1
# ╠═0b3dd8f7-a65b-447d-bb48-dddc7bf929ea
# ╠═1beeb7b4-fde9-4e5f-b5ed-740d712edf69
# ╠═3b16345c-080c-459f-af6c-6fe657b6d142
# ╠═7df8ff9b-773e-4b58-8f1b-6f8d6461ad41
# ╠═1602e2fb-801e-4e68-8bd6-24ad5b5e9671
# ╠═2b5eaa63-d745-4875-9359-69c330d5b70d
# ╠═58b6786b-6ca7-4401-ac06-53520782d5e9
# ╠═8e25cad9-4527-4059-8904-9ac13818e6f2
# ╠═da382dd0-0664-48dd-8e69-d83181354b62
# ╠═fd5b0002-ffca-42f1-8b17-541aae43b346
# ╠═141b6d46-8d29-4518-b55d-dec1549a8cbb
# ╠═13aed017-8de3-4b39-8b84-22f3b08fe41c
# ╠═77531f28-b09e-4aa6-9f83-b2d16df123b7
# ╠═accb7d7c-1a8f-47ff-b222-99a6fcc9f162
# ╠═31f1c405-8460-4a76-a097-7d378c978c71
# ╠═8771f9f8-c53f-4f9f-a51f-6bd1e6da15d4
# ╠═1c466b4d-2ea3-4041-b85e-cf5b1ea16d80
# ╠═46e84370-18e1-498c-a2b7-afdec69edc4a
# ╠═d65807e2-12a4-4f71-a7c6-3a57495e7007
# ╠═920587cd-bd4d-40ac-87c3-c9e28ffb5dea
# ╠═8e4aacb1-ee4d-4b96-b633-a009ca88ab0e
# ╠═57f2a900-aaa7-4bf8-90e9-9abcc83b5e1b
# ╠═64cc8f00-b29c-4b6f-b5fd-211e3f7dc438
# ╠═dd29fc2d-814a-4fc0-90cb-9447b47acfc2
# ╠═9b9777aa-c5de-4753-92cd-170c971aa547
# ╠═d890fc82-5ba8-4de6-aeb0-6821e505cd7e
# ╠═70a98bb9-889c-4732-b006-4c98e1e8a89b
# ╠═935d4fb6-26c3-44d5-9155-0f066af02824
# ╠═2ec88f23-5876-4760-9b3c-2dab9b232955
# ╠═b3165e34-fdf7-4253-ab05-bcfcddf9266a
# ╠═87d2466e-6c06-4261-aced-bf879adeb2ce
# ╠═c5f98e9f-88f5-47be-a90d-afcaea1fdede
# ╠═81322f55-2bcf-4b53-9b69-c7a690c304e9
# ╠═1fa9b627-120e-4169-8b1d-ae1f259f3b19
# ╠═f09115ca-10d2-45f7-8d33-808233f12664
# ╠═b4f7e3d3-adbd-437d-b2e3-4eadf4c1eb02
# ╠═b1753d95-3c8c-40e5-a18e-836eb7958ae4
# ╠═0e63a78e-e2b7-4343-95fc-2762084dd23f
# ╠═dbf7c73f-853c-415e-9c93-708227eaa45d
# ╠═d419f6a1-2823-4496-a791-38d71f8f5e69
# ╠═2af1562d-cecd-4ab4-b4eb-b1aa624872b4
# ╠═a063450f-57f7-4ea1-b2af-8a96b21ae0b7
# ╠═f619ec41-7354-478c-bce1-bf61d738e6f9
# ╠═db9b5e37-def5-4506-82a7-b71250348b47
# ╠═4bcf7744-af84-4584-aca7-6e3a96b7c1ed
# ╠═8a57c7ae-bb89-43df-b870-f44c759bedd7
# ╠═d7e7ba7d-0bfc-4ee7-afb2-e73afb0e27dc
# ╠═2a5f0d5f-661f-4076-a33c-8d4143d47b00
# ╠═20e5830b-f10d-4a91-ba28-44504865ff87
# ╠═13944fee-0eb2-416d-9e20-0a1701469df5
# ╠═912ba6aa-9838-46b5-a0c9-308b507f9491
# ╠═63944386-35b1-41b2-8d09-4353121d0ebb
# ╠═b83cfe12-765e-40a0-9cea-7eabaa859dc5
# ╠═a8fb3427-5dc5-472e-aa6d-f6e76033e119
# ╠═bf01f558-c98e-490c-a00f-9f039ec8fb86
# ╠═26416d72-cd1c-4ba0-9544-b7336373472e
# ╠═3ee0665c-9ff4-439a-8643-658c9a4fcf91
# ╠═5c8c33a4-62a0-4a45-b2f1-ffd0a0398e55
# ╠═c13d23e3-b261-4abc-ad09-51af74c8bff7
# ╠═6f1965b2-c5cc-4e31-b6cb-644574c0841a
# ╠═426b0ca3-a308-49ab-98c4-b5940848d62f
# ╠═3c6e1f92-0bee-4318-b212-235ac2f6bda0
# ╟─0fd5a508-bbfb-4887-b8e9-1a609122a8cf
# ╟─4b021012-2f4c-4b60-986a-137ebc8984cb
# ╟─c5086973-ad46-4920-aad9-0afa70d4ebe9
# ╟─61335eab-4677-4ff9-bfe7-379ee6c40008
# ╟─e53e8e8a-bb01-4695-997d-ad1eaa1deeae
# ╠═21b1b025-3e54-40de-b840-9e9087c4098d
# ╠═7529b04a-09d5-4a59-971b-8a501b73e248
# ╠═f0e39fec-2624-4199-bff5-abac5caa4494
# ╠═86e51b7f-5f25-4223-9e39-051f0ef316f0
# ╠═eb4be2cc-9458-4102-a2e9-6df30ba22045
# ╠═4e9b8f08-3875-4f55-a9f2-5c0c9edd16ed
# ╟─73180533-9668-4f83-bc8c-4da7b492b019
# ╠═b59fb57e-cec5-4f3d-8808-999f274f0679
# ╠═6eb2a227-ee65-456f-bfbc-0c727dfde3c9
# ╠═0df13a11-c319-4590-b684-ca0d409a7a29
# ╠═57ba1a2b-0f50-4c71-b105-6f4a0f61410f
# ╠═b2224158-0916-4b6e-a343-7887e21aca63
# ╠═dc337e00-3a48-4be5-bf1d-f5af674c8edc
# ╠═f41b02a2-0494-4ab1-9bfd-1f9aab535648
# ╠═95e5116e-bf8e-407a-b356-d0f431bf2dc4
