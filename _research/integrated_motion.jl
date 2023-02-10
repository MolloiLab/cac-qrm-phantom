### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 66611812-518f-4e1d-9123-8e33c5f8c16f
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, CairoMakie
    using StatsBase: quantile!

end

# ╔═╡ c8357f16-7ea1-4871-a7d5-f823eb5e3c18
TableOfContents()

# ╔═╡ 6bcbc55d-607b-4d06-95f8-20b5451e8e25
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ c0d9d335-b943-4b5d-8941-8fbdfcb9ad5a
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 6a193644-fdb2-40f0-972f-ac0e053c261c
BASE_PATH = joinpath(dirname(dirname(dirname(dirname(pwd())))), "Datasets", "CAC Data")

# ╔═╡ 15340cba-5e0f-4632-afaa-0df4956c5086
VENDER = "Stanford-Motion"

# ╔═╡ f1972836-28aa-4d7d-b514-0a733279cdba
root_path = joinpath(BASE_PATH, VENDER)

# ╔═╡ 7d783972-56b3-4bb9-b330-6c2424c69d51
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ cb3ab063-91b9-4cb8-88bf-e6b428b8d362
pth = dcm_path_list[68]

# ╔═╡ 8389df5e-7bb0-4fd0-927c-53b25b5c68bb
pth

# ╔═╡ de39da66-713e-4236-9284-32e7e11881f9
scan = basename(pth)

# ╔═╡ 22e5d932-cd55-4d5e-a08f-6f75a2db21ad
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 1d3401c1-e5f2-4a76-b623-bb177fdf8570
md"""
## Helper Functions
"""

# ╔═╡ a797214c-7e65-4408-89fc-fcd780604723
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 8eed3b5e-12cc-4a48-b842-f18f7e009b3b
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 9f928c54-911f-4563-9d46-e8ac73eb4b74
function overlay_mask_plot(array, mask, var, title::AbstractString)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	indices_lbl = findall(x -> x == zs[var], label_array[:,3])
	
	fig = Figure()
	ax = Makie.Axis(fig[1, 1])
	ax.title = title
	heatmap!(array[:, :, zs[var]], colormap=:grays)
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=0.6, color=:red)
	fig
end

# ╔═╡ c28446f5-e9cd-4051-beaf-dc49e4285044
md"""
## Segment Heart
"""

# ╔═╡ 000e728a-180d-41f3-863a-34e8d7f69d5b
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2; radius_val=88);

# ╔═╡ 46fb26a4-a4f3-4f28-83bc-05f3ea3ccd1a
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 978cff4e-d217-44c9-a42f-1c6d59770b36
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ c50b0b96-ae55-449d-bd27-b3e44ea1060f
let
	fig = CairoMakie.Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 10]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ 8c46d0fd-f3c0-446e-9163-98a93e7a89f5
let
	fig = CairoMakie.Figure()
	
	ax2 = Makie.Axis(fig[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ 5ede534f-5bab-476c-ab28-e1f4fff62f30
let
	fig = CairoMakie.Figure()
	
	ax3 = Makie.Axis(fig[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 10]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ 8eec4a39-1b44-42ac-862a-b52356dad5fa
md"""
## Segment Calcium Inserts
"""

# ╔═╡ a49a3561-3047-4bf5-b881-ccf4aadac550
function find_slices(filename)
	if occursin("rod1", filename)
		ROD1 = true
		ROD2 = false
	elseif occursin("rod2", filename)
		ROD1 = false
		ROD2 = true
	end

	if occursin("pos1", filename)
		POS1 = true
		POS2 = false
		POS3 = false
		POS4 = false
		POS5 = false
	elseif occursin("pos2", filename)
		POS1 = false
		POS2 = true
		POS3 = false
		POS4 = false
		POS5 = false
	elseif occursin("pos3", filename)
		POS1 = false
		POS2 = false
		POS3 = true
		POS4 = false
		POS5 = false
	elseif occursin("pos4", filename)
		POS1 = false
		POS2 = false
		POS3 = false
		POS4 = true
		POS5 = false
	elseif occursin("pos5", filename)
		POS1 = false
		POS2 = false
		POS3 = false
		POS4 = false
		POS5 = true
	end
	
	if ROD1
		if POS1 || POS4 || POS5
			slices1 = [4, 7]
			slices2 = [14, 17]
		elseif POS2 || POS3
			slices1 = [6, 10]
			slices2 = [16, 19]
		end

	elseif ROD2
		if POS1 || POS2 || POS5
			slices1 = [2, 5]
			slices2 = [12, 15]
		elseif POS3
			slices1 = [3, 6]
			slices2 = [13, 16]
		elseif POS4
			slices1 = [6, 9]
			slices2 = [16, 20]
		end
	end
	return slices1, slices2
end

# ╔═╡ 25ddf768-c878-4e7f-b16b-a0f9e8a5e416
slices1, slices2 = find_slices(scan)

# ╔═╡ dc525947-a3d8-4127-8e05-e1c3bed7975d
mask1_small, mask2_small = mask_inserts_motion(masked_array, header, slices1, slices2; threshold=115, radius=3);

# ╔═╡ 0dd4e077-118c-4a88-af9a-224a651d1c13
heatmap(mask2_small, colormap=:grays)

# ╔═╡ e4e59956-9ae1-4506-8ed5-2b0f7066f6a4
begin
	msk = mask2_small
	idxs = findall(isone, msk)
	idxs = getindex.(idxs, [1 2])
end;

# ╔═╡ 08d50be1-df5f-4c3f-80be-4fcbfe4bcfe8
@bind m1 PlutoUI.Slider(axes(masked_array, 3); default=15, show_value=true)

# ╔═╡ 37562cf1-3d3e-406b-9acb-e49fe200418c
let
	f = Figure()

	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(masked_array[:, :, m1], colormap=:grays)
	scatter!(idxs, color=:red, markersize=2)

	f
end

# ╔═╡ e2f85089-9292-4fa4-9f10-5e1738a5280a
begin
	masked_slice = masked_array[:, :, m1]
	hu_calcium_1 = mean(masked_slice[msk])
end

# ╔═╡ daf98a91-9533-4741-a450-0fe7a758739d
rod1_ca = Dict(
	800 => [956, 936, 910, 937, 917],
	408 => [506, 488, 483, 501, 501]
)

# ╔═╡ c3b24c77-8f9f-4979-b333-d5066d2c93d0
maximum(rod1_ca[408])

# ╔═╡ 0b34a376-d97e-4ba8-80b3-0619429d528a
rod2_ca = Dict(
	380 => [434, 436, 415, 379, 439],
	196 => [201, 218, 239, 216, 203]
)

# ╔═╡ 0fb6db5b-44e4-44e7-be49-eb1fc7cbb458
begin
	slice_cents_1 = Dict(6 => 956, 15 => 506) # slice center => intensity
	slice_cents_8 = Dict(8 => 936, 18 => 488)
	slice_cents_15 = Dict(8 => 910, 18 => 483)
	slice_cents_22 = Dict(5 => 937, 15 => 501)
	slice_cents_28 = Dict(5 => 917, 15 => 501)
	slice_cents_35 = Dict(5 => 434, 14 => 201)
	slice_cents_42 = Dict(4 => 436, 13 => 218)
	slice_cents_49 = Dict(4 => 415, 16 => 239)
	slice_cents_55 = Dict(8 => 379, 18 => 216)
	slice_cents_62 = Dict(3 => 439, 14 => 203)
end

# ╔═╡ 39c5dc12-ac4b-4c3e-8305-1366923d0168
md"""
## Calibration Prep
"""

# ╔═╡ e1338891-0f12-4851-890f-dc37cf6b5143
cal_root = joinpath(dirname(dirname(pwd())), "cac-simulation", "output_new")

# ╔═╡ 1f7d75b1-4f41-4f75-bb1a-8b4c63365859
cal_path = joinpath(cal_root, "calibrations.csv")

# ╔═╡ 3cf9fd09-b14e-4390-8e75-4d0fb7531674
cal_df = CSV.read(cal_path, DataFrame)

# ╔═╡ 58a6415a-d129-4a89-8fb1-ea7c076f3fc0
kvp = header[(0x0018, 0x0060)]

# ╔═╡ eefde8ed-89c6-464c-9d34-a8190f517eb1
begin
	if kvp == 80.0
		arr = Array(cal_df[1, 2:end])
		intensity_array3 = vcat(0, arr)
	elseif kvp == 100.0
		arr = Array(cal_df[2, 2:end])
		intensity_array3 = vcat(0, arr)
	elseif kvp == 120.0
		arr = Array(cal_df[3, 2:end])
		intensity_array3 = vcat(0, arr)
	end
end

# ╔═╡ 580fac3d-465f-4ce0-ab11-e8d0ebfc6d75
density_array_calc3 = [
		0
		25
		50
		100
		200
		400
		800
	]

# ╔═╡ a9c5bad0-727d-4705-9447-3d11840bd274
begin
	df_cal = DataFrame(:density => density_array_calc3, :intensity => intensity_array3)
	linearRegressor = lm(@formula(intensity ~ density), df_cal)
	linearFit = predict(linearRegressor)
	m = linearRegressor.model.pp.beta0[2]
	b = linearRegressor.model.rr.mu[1]
	density(intensity) = (intensity - b) / m
	intensity(ρ) = m*ρ + b
end

# ╔═╡ 35ac3beb-0867-4e7a-9453-15e9ab4b6502
md"""
# Score Inserts
"""

# ╔═╡ 1347149c-e9a4-4696-ac60-ca7a9797bdd3
md"""
## Insert 1
"""

# ╔═╡ 566e0fcf-2436-42e4-b5be-c611e4b1ddae
slices1

# ╔═╡ 2774a80b-c590-4d71-8e39-bdcaa7178dcc
insert1_array = masked_array[:, :, slices1[1]-1:slices1[2]+1];

# ╔═╡ 16390207-8bce-4d24-9875-2e54edad5c0f
md"""
#### Dilated mask
"""

# ╔═╡ c87ed077-d127-4342-9804-b9df2a1771ea
md"""
#### Ring (background) mask
"""

# ╔═╡ 02702431-09eb-4089-a41e-546a6459d148
S_Obj1 = intensity(800)

# ╔═╡ 971a04f2-c3eb-4f3b-819a-94acdd18e075
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ a01ab968-2704-4b95-a559-8ade6ebfbfcc
md"""
## Insert 2
"""

# ╔═╡ dede7b0d-3461-4ade-b2c3-140c2662b5bd
slices2

# ╔═╡ 9bf7461f-25c2-4789-b49a-97eebdd3d531
insert2_array = masked_array[:, :, slices2[1]-1:slices2[2]+1];

# ╔═╡ 1e0c8be0-bb80-4212-9494-a470639da812
md"""
#### Dilated mask
"""

# ╔═╡ b64dc067-db33-4fbb-9043-e87cffee03b5
md"""
#### Ring (background) mask
"""

# ╔═╡ 60662f6e-456c-4e4f-ac97-32eb35524e52
S_Obj2 = intensity(800)

# ╔═╡ 7a83f6a3-0911-45ed-9001-d9971d362bca
md"""
# Results
"""

# ╔═╡ 0ca4d31e-5a88-4ea0-a436-4f0ec4127697
dens = [
	0.196
	0.380
	0.408
	0.800
] # mg/mm^3

# ╔═╡ d97f7848-25ae-4012-ba21-d053689b8b95
vol = 196.35 # mm^3

# ╔═╡ c2fff2e2-e7e0-430e-b1f7-e90e4b0d4a90
masses = vol .* dens

# ╔═╡ 963b8e26-1cc6-4c79-9680-62dbb4c6e86c
inserts = [
	"insert 1",
	"insert 2"
]

# ╔═╡ 77e95317-4da5-4ee1-8b30-73b2e338ca17
if occursin("rod1", scan)
	rod1 = true
	ground_truth_mass = masses[3:4]
elseif occursin("rod2", scan)
	rod1 = false
	ground_truth_mass = masses[1:2]
end

# ╔═╡ 1dca349f-1cd2-4f18-abe5-0b83c80a4450
ground_truth_mass[2]

# ╔═╡ 458e7816-0066-4962-b432-ffae502b4e0c
ground_truth_mass[1]

# ╔═╡ 25ac6562-39ea-4ff0-b0bf-746b1a1a3fa2
pth

# ╔═╡ ef800e1b-5913-465a-8153-dd0e9c077778
densities = dens .* 1000

# ╔═╡ 22c7ecaf-0ed2-4c05-b80f-3e7b397b7fad
md"""
## Save Results
"""

# ╔═╡ 62e6d211-48ca-48ce-881b-90bc901bc49f
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER))
end

# ╔═╡ 0eaafda4-bce0-414f-b214-f9fdd1208666
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "/", scan, ".csv")

# ╔═╡ d9ed1edf-b64a-494a-8e6d-4359e0d0c48f
# CSV.write(output_path, df)

# ╔═╡ 768d5e24-e519-4adb-8143-22ab0b9c3561
md"""
## Save full df
"""

# ╔═╡ b008d1c1-486b-493a-afaf-671fa2dbb73d
SCAN_NUMBER = 7

# ╔═╡ d9767fd8-852f-4f08-bc77-48324492979d
pth

# ╔═╡ e3182af9-6483-4cdb-b4b6-886d5b624d30
rad = 22

# ╔═╡ cb632ff9-be56-442b-a1bf-73aad48add49
mask1, mask2 = mask_inserts_motion(masked_array, header, slices1, slices2; threshold=115, radius=rad);

# ╔═╡ 053a12a8-4857-4629-a8ae-d9c86575016f
heatmap(mask2, colormap=:grays)

# ╔═╡ 333b91d8-8cb7-4b7e-b07a-9423bdf6243f
begin
	mask1_3D = Array{Bool}(undef, size(insert1_array))
	for z in 1:size(insert1_array, 3)
		mask1_3D[:, :, z] = mask1
	end
end;

# ╔═╡ c7d14e38-c813-436f-a37d-aa2163265901
dilated_mask1 = dilate(dilate(mask1_3D));

# ╔═╡ 019103ef-06ae-4c85-8fe7-6dc7f601fae5
@bind g2 overlay_mask_bind(dilated_mask1)

# ╔═╡ 4e58d18d-242b-4427-8b4e-5901784f78aa
overlay_mask_plot(insert1_array, dilated_mask1, g2, "dilated mask")

# ╔═╡ 9c50ceb4-bdc1-4743-b9ff-414fbc94b478
ring_mask1 = dilate(dilate(dilated_mask1)) - dilated_mask1;

# ╔═╡ 7bbe4340-53af-43de-bccb-2ad255a82b4b
@bind g4 overlay_mask_bind(ring_mask1)

# ╔═╡ 39608ff6-0f80-44e1-8d18-a00b1e7d3108
overlay_mask_plot(insert1_array, ring_mask1, g4, "ring mask")

# ╔═╡ 2a139bac-8201-4d92-add2-5bfcd6dff3d3
s_bkg1 = mean(insert1_array[Bool.(ring_mask1)])

# ╔═╡ 7a6eb472-e18a-4f20-94f9-d41ec701a12d
begin
	alg1 = Integrated(insert1_array[Bool.(mask1_3D)])
	ρ1 = 0.8 # mg/mm^3
	mass1 = score(s_bkg1, S_Obj1, pixel_size, ρ1, alg1)
end

# ╔═╡ 6fe8ea23-7c76-4ee9-8284-e7a68de07d0b
begin
	mask2_3D = Array{Bool}(undef, size(insert2_array))
	for z in 1:size(insert2_array, 3)
		mask2_3D[:, :, z] = mask2
	end
end;

# ╔═╡ 9760db44-abb7-4833-9676-62851eaade66
dilated_mask2 = dilate(dilate(mask2_3D));

# ╔═╡ c3a7a6e6-f342-4579-b895-c5f68a18c363
@bind g3 overlay_mask_bind(dilated_mask2)

# ╔═╡ 0d680ef5-1fd8-44ac-ace4-9aede54a2fb1
overlay_mask_plot(insert2_array, dilated_mask2, g3, "dilated mask")

# ╔═╡ 90313d8b-01e1-410f-a902-aa2028330d3c
ring_mask2 = dilate(dilate(dilated_mask2)) - dilated_mask2;

# ╔═╡ b36c02e1-7496-45e3-9c4e-2c8a181719e7
@bind g5 overlay_mask_bind(ring_mask2)

# ╔═╡ e87d84cc-0b17-43e2-a12d-31c71a6ed1dc
overlay_mask_plot(insert2_array, ring_mask2, g5, "ring mask")

# ╔═╡ 8b4b7587-96a7-4646-a6e8-8688e5c693a3
s_bkg2 = mean(insert2_array[Bool.(ring_mask2)])

# ╔═╡ f6ae439a-7657-4ac2-92b1-d7dab557e391
begin
	alg2 = Integrated(insert2_array[Bool.(mask2_3D)])
	ρ2 = 0.8 # mg/mm^3
	mass2 = score(s_bkg2, S_Obj2, pixel_size, ρ2, alg2)
end

# ╔═╡ 9e99a632-8903-4382-8278-ac718e9f24a1
calculated_mass = [
	mass2,
	mass1,
]

# ╔═╡ 3c0e83cd-c5f4-4f3b-9a46-0c55a45cdb41
df = DataFrame(
	scan = scan,
	rad = rad,
	inserts = inserts,
	ground_truth_mass = ground_truth_mass,
	calculated_mass = calculated_mass,
)

# ╔═╡ 82e3af51-3fbc-4be5-8c09-4b4a67ab1e85
let
	f = Figure()
	ax1 = Axis(f[1, 1])
	ax1.title = "Mass Measurements"
	ax1.ylabel = "Mass (mg)"
	ax1.xlabel = "Density (mg/mm^3)"
	if rod1
		scatter!(densities[3:4], df[!, :ground_truth_mass], label="ground_truth_mass")
		scatter!(densities[3:4], df[!, :calculated_mass], label="calculated_mass")
		
		
	else
		scatter!(densities[1:2], df[!, :ground_truth_mass], label="ground_truth_mass")
		scatter!(densities[1:2], df[!, :calculated_mass], label="calculated_mass")
	end

	f[1, 2] = Legend(f, ax1, framevisible = false)
	f
end

# ╔═╡ 9df0362e-af34-4d70-ae3e-740778d70bc0
dfs = []

# ╔═╡ 953e4214-019a-4d22-b064-46ed33cde266
push!(dfs, df)

# ╔═╡ e11301e4-51f5-4589-99e7-c84716d6899e
# begin
# 	dfa = dfs[3]
# 	dfb = dfs[6]
# 	rms1 = rmsd(dfa[!, :ground_truth_mass], dfa[!, :calculated_mass])
# 	rms2 = rmsd(dfb[!, :ground_truth_mass], dfb[!, :calculated_mass])
# 	rms1, rms2
# end

# ╔═╡ 828062e1-e637-48dd-8083-96eb954979fa
# if length(dfs) == 10
# 	global new_df = vcat(dfs[1:10]...)
# 	output_path_new = string(cd(pwd, "..") , "/data/output/", VENDER, "/full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═66611812-518f-4e1d-9123-8e33c5f8c16f
# ╠═c8357f16-7ea1-4871-a7d5-f823eb5e3c18
# ╟─6bcbc55d-607b-4d06-95f8-20b5451e8e25
# ╟─c0d9d335-b943-4b5d-8941-8fbdfcb9ad5a
# ╠═6a193644-fdb2-40f0-972f-ac0e053c261c
# ╠═15340cba-5e0f-4632-afaa-0df4956c5086
# ╠═f1972836-28aa-4d7d-b514-0a733279cdba
# ╠═7d783972-56b3-4bb9-b330-6c2424c69d51
# ╠═cb3ab063-91b9-4cb8-88bf-e6b428b8d362
# ╠═8389df5e-7bb0-4fd0-927c-53b25b5c68bb
# ╠═de39da66-713e-4236-9284-32e7e11881f9
# ╠═22e5d932-cd55-4d5e-a08f-6f75a2db21ad
# ╟─1d3401c1-e5f2-4a76-b623-bb177fdf8570
# ╟─a797214c-7e65-4408-89fc-fcd780604723
# ╟─8eed3b5e-12cc-4a48-b842-f18f7e009b3b
# ╟─9f928c54-911f-4563-9d46-e8ac73eb4b74
# ╟─c28446f5-e9cd-4051-beaf-dc49e4285044
# ╠═000e728a-180d-41f3-863a-34e8d7f69d5b
# ╟─46fb26a4-a4f3-4f28-83bc-05f3ea3ccd1a
# ╟─978cff4e-d217-44c9-a42f-1c6d59770b36
# ╟─c50b0b96-ae55-449d-bd27-b3e44ea1060f
# ╟─8c46d0fd-f3c0-446e-9163-98a93e7a89f5
# ╟─5ede534f-5bab-476c-ab28-e1f4fff62f30
# ╟─8eec4a39-1b44-42ac-862a-b52356dad5fa
# ╟─a49a3561-3047-4bf5-b881-ccf4aadac550
# ╠═25ddf768-c878-4e7f-b16b-a0f9e8a5e416
# ╠═cb632ff9-be56-442b-a1bf-73aad48add49
# ╠═053a12a8-4857-4629-a8ae-d9c86575016f
# ╠═dc525947-a3d8-4127-8e05-e1c3bed7975d
# ╠═0dd4e077-118c-4a88-af9a-224a651d1c13
# ╠═e4e59956-9ae1-4506-8ed5-2b0f7066f6a4
# ╟─08d50be1-df5f-4c3f-80be-4fcbfe4bcfe8
# ╟─37562cf1-3d3e-406b-9acb-e49fe200418c
# ╠═e2f85089-9292-4fa4-9f10-5e1738a5280a
# ╠═daf98a91-9533-4741-a450-0fe7a758739d
# ╠═c3b24c77-8f9f-4979-b333-d5066d2c93d0
# ╠═0b34a376-d97e-4ba8-80b3-0619429d528a
# ╠═0fb6db5b-44e4-44e7-be49-eb1fc7cbb458
# ╟─39c5dc12-ac4b-4c3e-8305-1366923d0168
# ╠═e1338891-0f12-4851-890f-dc37cf6b5143
# ╠═1f7d75b1-4f41-4f75-bb1a-8b4c63365859
# ╠═3cf9fd09-b14e-4390-8e75-4d0fb7531674
# ╠═58a6415a-d129-4a89-8fb1-ea7c076f3fc0
# ╠═eefde8ed-89c6-464c-9d34-a8190f517eb1
# ╠═580fac3d-465f-4ce0-ab11-e8d0ebfc6d75
# ╠═a9c5bad0-727d-4705-9447-3d11840bd274
# ╟─35ac3beb-0867-4e7a-9453-15e9ab4b6502
# ╟─1347149c-e9a4-4696-ac60-ca7a9797bdd3
# ╠═566e0fcf-2436-42e4-b5be-c611e4b1ddae
# ╠═2774a80b-c590-4d71-8e39-bdcaa7178dcc
# ╠═333b91d8-8cb7-4b7e-b07a-9423bdf6243f
# ╟─16390207-8bce-4d24-9875-2e54edad5c0f
# ╠═c7d14e38-c813-436f-a37d-aa2163265901
# ╟─019103ef-06ae-4c85-8fe7-6dc7f601fae5
# ╠═4e58d18d-242b-4427-8b4e-5901784f78aa
# ╟─c87ed077-d127-4342-9804-b9df2a1771ea
# ╠═9c50ceb4-bdc1-4743-b9ff-414fbc94b478
# ╟─7bbe4340-53af-43de-bccb-2ad255a82b4b
# ╠═39608ff6-0f80-44e1-8d18-a00b1e7d3108
# ╠═2a139bac-8201-4d92-add2-5bfcd6dff3d3
# ╠═02702431-09eb-4089-a41e-546a6459d148
# ╠═971a04f2-c3eb-4f3b-819a-94acdd18e075
# ╠═7a6eb472-e18a-4f20-94f9-d41ec701a12d
# ╠═1dca349f-1cd2-4f18-abe5-0b83c80a4450
# ╟─a01ab968-2704-4b95-a559-8ade6ebfbfcc
# ╠═dede7b0d-3461-4ade-b2c3-140c2662b5bd
# ╠═9bf7461f-25c2-4789-b49a-97eebdd3d531
# ╠═6fe8ea23-7c76-4ee9-8284-e7a68de07d0b
# ╟─1e0c8be0-bb80-4212-9494-a470639da812
# ╠═9760db44-abb7-4833-9676-62851eaade66
# ╟─c3a7a6e6-f342-4579-b895-c5f68a18c363
# ╠═0d680ef5-1fd8-44ac-ace4-9aede54a2fb1
# ╟─b64dc067-db33-4fbb-9043-e87cffee03b5
# ╠═90313d8b-01e1-410f-a902-aa2028330d3c
# ╟─b36c02e1-7496-45e3-9c4e-2c8a181719e7
# ╠═e87d84cc-0b17-43e2-a12d-31c71a6ed1dc
# ╠═8b4b7587-96a7-4646-a6e8-8688e5c693a3
# ╠═60662f6e-456c-4e4f-ac97-32eb35524e52
# ╠═f6ae439a-7657-4ac2-92b1-d7dab557e391
# ╠═458e7816-0066-4962-b432-ffae502b4e0c
# ╟─7a83f6a3-0911-45ed-9001-d9971d362bca
# ╠═0ca4d31e-5a88-4ea0-a436-4f0ec4127697
# ╠═d97f7848-25ae-4012-ba21-d053689b8b95
# ╠═c2fff2e2-e7e0-430e-b1f7-e90e4b0d4a90
# ╠═963b8e26-1cc6-4c79-9680-62dbb4c6e86c
# ╠═77e95317-4da5-4ee1-8b30-73b2e338ca17
# ╠═9e99a632-8903-4382-8278-ac718e9f24a1
# ╠═25ac6562-39ea-4ff0-b0bf-746b1a1a3fa2
# ╠═3c0e83cd-c5f4-4f3b-9a46-0c55a45cdb41
# ╠═ef800e1b-5913-465a-8153-dd0e9c077778
# ╟─82e3af51-3fbc-4be5-8c09-4b4a67ab1e85
# ╟─22c7ecaf-0ed2-4c05-b80f-3e7b397b7fad
# ╠═62e6d211-48ca-48ce-881b-90bc901bc49f
# ╠═0eaafda4-bce0-414f-b214-f9fdd1208666
# ╠═d9ed1edf-b64a-494a-8e6d-4359e0d0c48f
# ╟─768d5e24-e519-4adb-8143-22ab0b9c3561
# ╠═b008d1c1-486b-493a-afaf-671fa2dbb73d
# ╠═d9767fd8-852f-4f08-bc77-48324492979d
# ╠═e3182af9-6483-4cdb-b4b6-886d5b624d30
# ╠═9df0362e-af34-4d70-ae3e-740778d70bc0
# ╠═953e4214-019a-4d22-b064-46ed33cde266
# ╠═e11301e4-51f5-4589-99e7-c84716d6899e
# ╠═828062e1-e637-48dd-8083-96eb954979fa
