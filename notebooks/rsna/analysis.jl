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
	using PlutoUI, Statistics, CSV, DataFrames, CairoMakie, Colors, GLM, MLJBase, Printf
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
df_i = CSV.read(joinpath(root_path, "volume_fraction", "physical.csv"), DataFrame);

# ╔═╡ 9dd9fc84-c9ee-48d1-b444-0cf3ae5d435d
let
	df = df_i
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
	
	df = df_i
	sc1=scatter!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large])
	# errorbars!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
	sc2=scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	# errorbars!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
	sc3=scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	# errorbars!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
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
	# errorbars!(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large], rms(df[!, :ground_truth_mass_large], df[!, :calculated_mass_large]))
	scatter!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium])
	# errorbars!(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium], rms(df[!, :ground_truth_mass_medium], df[!, :calculated_mass_medium]))
	scatter!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], color=:red)
	# errorbars!(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small], rms(df[!, :ground_truth_mass_small], df[!, :calculated_mass_small]))
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
	
	save(plotsdir("linear_reg.eps"), f)
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

# ╔═╡ 6164a556-e07a-46bd-8793-eff9d30bdbde
df_a

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
#### Integrated
"""

# ╔═╡ 8312248f-8367-4a11-a9c8-5431695532b3
begin
	false_negative_i = []
	for i in 1:3:nrow(df_i)-2
		mean_i, std_i = mean(df_i[i:i+2, :mass_bkg]), std(df_i[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_i[i:i+2, :calculated_mass_large], df_i[i:i+2, :calculated_mass_medium], df_i[i:i+2, :calculated_mass_small]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i, neg)
	end
end

# ╔═╡ 3b47b01a-5f6d-4e62-bb6d-64a8b8b68aad
total_zero_i = sum(false_negative_i)

# ╔═╡ e75a9cc0-03ea-4de9-9c5e-b419bb30b68c
total_zero_i, num_zero_a

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
#### Integrated
"""

# ╔═╡ aa6214d2-4e17-4667-a88f-c05d11b10149
begin
	false_positive_i = []
	for i in 1:3:nrow(df_i)-2
		mean_i, std_i = mean(df_i[i:i+2, :mass_bkg]), std(df_i[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i, pos)
	end
end

# ╔═╡ be13830f-abe0-448e-86eb-f47f0954bc75
total_zero_i_pos = sum(false_positive_i)

# ╔═╡ 018331c8-ee76-475f-80d8-1795c9f4a857
function sensitivity_specificity()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:2, ["Integrated", "Agatston"]),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2]
	h1 = (total_zero_i / total_cac) * 100
	h2 = (num_zero_a / total_cac) * 100 
	heights1 = [h1, h2]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
    barplot!(table, heights1; color=colors[1:2], bar_labels=[l1, l2])

    ylims!(ax; low=0, high=100)

	##-- B --##
	ax = Axis(
		f[2, 1]; 
		xticks = (1:2, ["Integrated", "Agatston"]),
		title = "False-Positive (CAC>0)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2]
	h1 = (total_zero_i_pos / total_cac_pos) * 100
	h2 = (total_zero_a_pos / total_cac_pos) * 100
    heights1 = [h1, h2]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
    barplot!(table, heights1; color=colors[1:2], bar_labels=[l1, l2])

    ylims!(ax; low=0, high=100)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save(plotsdir("sensitivity_specificity.eps"), f)

    return f
end

# ╔═╡ 6c476510-3168-4979-ab42-eb9cbff2f971
with_theme(medphys_theme) do
    sensitivity_specificity()
end

# ╔═╡ e8a0a64f-bef3-414d-b7c0-0c3d5122f8d1
total_zero_i_pos, total_zero_a_pos

# Cell order:
# ╠═f36e510f-ac3a-4986-bf56-d3b80a068820
# ╠═d750f7ce-a801-11ed-04b2-1575dfb3a26e
# ╠═f1ae7024-bc5e-43c0-ad67-d58787f67004
# ╠═6bd3b773-d4cd-430b-b63a-91cc04ff1c95
# ╠═c3ee18d6-847b-4eb8-ac80-888a2e69c12f
# ╠═4f43a4f5-4dd2-4366-ab5b-b8d0089b1dd9
# ╟─9cf7599f-667f-45cf-a111-95149a7b2b0a
# ╟─2d2bdc53-e225-4369-a82b-f56ac153263e
# ╟─aa9185b1-e0d6-46db-9867-d6e44bcf70f8
# ╠═f4579ec7-9013-4ed6-aa32-5fb6a7de82b9
# ╠═07b01b11-a479-4a1d-88e2-6bb3ce1952c4
# ╠═22e2507a-3153-4cde-b95c-ca50108dd78a
# ╠═38d3637a-3407-46d9-a91c-5526198f6ce9
# ╠═12c47081-bb76-4ba0-ba5d-af15b21991bd
# ╠═22994bc2-8b11-48a2-bf56-68e592b0d702
# ╠═6e9e50e7-66b8-460a-8ebf-d5cf83a2443d
# ╠═5fe1ae64-d676-479f-8d2d-463740c48ae5
# ╟─45d93a55-8f32-4a44-bdd7-abac89bb41cc
# ╟─452d14b2-cc8f-4814-9b0b-679492f8426a
# ╟─2915b2a1-416e-4365-b0b4-fd2b07c8c1a7
# ╠═2abe98b8-91a9-42ed-b18a-192e7922efd2
# ╠═124809eb-61aa-4151-8608-3907d235560f
# ╟─9dba29d0-0402-48da-94ff-34ef7b64f6e7
# ╠═781a4de9-6b5b-4791-87f1-25207398fefc
# ╟─0db59538-fc2f-438c-9d2e-858714ace83b
# ╟─dc47bcc1-66a1-4137-908a-8301021a2f72
# ╠═b3f2e899-36dc-45d7-a61d-de14f5a8836c
# ╠═9fbd3b35-6266-438e-97bd-579592f65727
# ╠═3ce554d6-28f6-42b0-8926-ad2f7543cf4f
# ╠═4d55c88e-9506-4f23-8d50-6251db07841d
# ╠═070879f7-2b29-44a2-9644-49e02f48c51e
# ╠═a43bdcbb-8abc-4e8c-83b4-29547f41ad39
# ╠═3424c91c-c592-4d1d-b226-b4d1ae4b3563
# ╟─ce8a2deb-76f5-4b58-a51e-7e02c77782b2
# ╠═20533ef0-7222-41af-a8ad-6744cf159b78
# ╠═273414ff-5e2a-4338-8527-58f0952d52ba
# ╟─ec1d2254-040f-40de-83b9-6306bc869a89
# ╟─5a3a85b4-545f-4514-88d2-5df807a25aff
# ╠═02d80c8b-0ccf-4e79-8cb0-27d97a8aa390
# ╠═a483459c-da31-4ae5-bf92-d585824c4c6a
# ╠═0544b3db-432e-406c-bc88-c3ebcc4c34a3
# ╟─2b5c1c6d-46cc-4474-8927-8309d769d617
# ╠═a497ba6b-61ee-48c6-b7af-e8454dbb613e
# ╠═233f365d-6ed2-4a4d-9533-ffde89972401

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
# ╠═8bcaad17-4a5e-4c73-bb4d-04f4372e4143
# ╠═95a57c8f-7d55-4619-bf83-55a1802c3a4e
# ╠═3951d189-d3bb-457e-912f-5ca80dc26bd4
# ╠═71ccd8bf-ac90-4278-9f22-9ac536f33f48
# ╠═3311fd03-36f2-41b5-b261-4df42c448a36
# ╠═d186f370-48ab-408f-a1b7-e6a8a61f6b31
# ╠═f142b1b1-e861-48f1-bc32-6e8cf00ac716
# ╠═6164a556-e07a-46bd-8793-eff9d30bdbde
# ╠═1d072000-ed66-4a96-9ee8-9bdadafc92cc
# ╠═11d03fd2-a1ff-4c3f-95c5-08027013bcf2
# ╠═a20dbe6b-1c1b-4d03-ac10-a8ee9c935064
# ╠═2d65ee51-2771-42ef-83eb-a5ab740ddbf1
# ╠═d35cfb16-3938-4334-833a-e5e8f41741b6
# ╠═c06ea1d6-83ef-41b5-a6eb-34aff7b89287
# ╠═0a781663-cf52-4131-9f86-7054f9ee8284
# ╠═8312248f-8367-4a11-a9c8-5431695532b3
# ╠═3b47b01a-5f6d-4e62-bb6d-64a8b8b68aad
# ╠═e75a9cc0-03ea-4de9-9c5e-b419bb30b68c
# ╠═7cc4e287-1240-42d6-9c7b-b5f3601b3c12
# ╠═072feb69-adb1-408e-a920-23a409205703
# ╠═57b2d301-d7c8-4b87-9b1f-9f751984d626
# ╠═a459f8db-a2cb-4b27-9f47-9096504695c6
# ╠═36480397-b605-4cea-b63f-34bd86861a98
# ╠═d33126d5-0c4f-4497-ade9-82701e5c6ef8
# ╠═aa6214d2-4e17-4667-a88f-c05d11b10149
# ╠═be13830f-abe0-448e-86eb-f47f0954bc75
# ╠═018331c8-ee76-475f-80d8-1795c9f4a857
# ╠═6c476510-3168-4979-ab42-eb9cbff2f971
# ╠═e8a0a64f-bef3-414d-b7c0-0c3d5122f8d1
