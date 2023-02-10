### A Pluto.jl notebook ###
# v0.19.14

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

# ╔═╡ 795a7fe3-22ea-4367-b7db-3941fd9a4888
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")
	
	using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
	using StatsBase: quantile!
	
end

# ╔═╡ a1d626e2-fc5b-484e-9838-5a8463d4617a
TableOfContents()

# ╔═╡ 2bfa5dfc-472d-40fa-86e0-754a3e8abb38
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ f16eecad-68d2-42ae-a4dc-9b574984fd26
begin
	SCAN_NUMBER = 7
	VENDER = "Stanford-Motion"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ bb021d41-db85-4311-b16d-e4303bc38f63
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ d574f97c-bad2-413c-ac73-18363cc559fe
root_path = string(BASE_PATH, VENDER)

# ╔═╡ e0467b68-aecf-4ff8-a221-59b8d7d1a10a
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 75d8379a-1ddb-4127-a4f6-bc88d6ce6d14
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ ca7fd376-b7eb-43cd-8b74-dec6431d1823
scan = basename(pth)

# ╔═╡ 4dee1592-8848-4fb7-801f-c56404018046
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ ea20f143-6cfe-4662-bddb-97c5a15fe81f
md"""
## Helper Functions
"""

# ╔═╡ cb24aac3-1ce5-40b6-a843-00de6c998f2d
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 1eca10c3-2259-4fec-a4e4-f52b7be45b02
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 2c72f868-dda1-4baa-a200-511cbaa2755a
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
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=1, color=:red)
	fig
end

# ╔═╡ c5e319a1-ac8b-4595-b249-b94815116c20
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ d068bbe8-f2e4-4697-bcfd-b0a9333a7957
md"""
## Segment Heart
"""

# ╔═╡ 6d646443-5d78-4711-baac-25bf60fe1515
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2; radius_val=88);

# ╔═╡ 15897036-6264-4451-a1b1-5e42c7c54a8e
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 97f143cf-4725-493d-a6c6-a9693baddf3b
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 36c46879-c037-4236-b106-45012797f1d2
let
	f = Figure()
	ax = Makie.Axis(f[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	f
end

# ╔═╡ 55dcd1d5-1cdb-41be-9c05-66395905a880
let
	f = Figure()
	ax = Makie.Axis(f[1, 1])
	ax.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	f
end

# ╔═╡ c73647f9-489b-4af4-ab78-3e2c85925a11
let
	f = Figure()
	
	ax = Makie.Axis(f[1, 1])
	ax.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 10]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	f
end

# ╔═╡ b49c9980-6648-49b0-8a2f-82b2e180cb1d
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 330b7534-eefa-4946-a33c-efdac78bb213
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

# ╔═╡ 7d9a6306-38ed-467d-b55c-a70ce9f8242a
slices1, slices2 = find_slices(scan)

# ╔═╡ 61268b09-b7db-44d7-8e75-47f7fe4ff5e9
mask1, mask2 = mask_inserts_motion(masked_array, header, slices1, slices2; threshold=115, radius=10);

# ╔═╡ 34710f8d-822c-49a1-a4f9-4cf42a0a00e8
heatmap(mask2, colormap=:grays)

# ╔═╡ da315970-7db8-4070-9db2-e3d5b21b3a84
md"""
## Mass cal factor
"""

# ╔═╡ e4334739-f94b-473f-88bd-d58767fb9fb0
output = calc_output_motion(masked_array, header, slices1);

# ╔═╡ fe9dc076-bf45-44d4-8b6a-a74b839dfd98
center_large_LD = output[4][1]

# ╔═╡ 9d024419-81c2-42c0-8615-918f9a5ba2a9
rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])

# ╔═╡ fb9de4d3-4e19-4c87-80e3-4086908ca48a
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 7a9f5566-0439-4005-a26b-e992613480b1
avg_mass_cals = mean([
0.000931147497044372
0.000931147497044372
0.000931147497044372
0.0009211153369385057
0.0009211153369385057
0.0009211153369385057
0.0009096665764441353
0.0009096665764441353
0.0009096665764441353
0.0009031593242901845
0.0009031593242901845
0.0009031593242901845
0.0009012665464632671
0.0009012665464632671
0.0009012665464632671
0.0008848386753622677
0.0008848386753622677
0.0008848386753622677
0.0008846639896996343
0.0008846639896996343
0.0008846639896996343
0.0008830535691512911
0.0008830535691512911
0.0008830535691512911
0.0008828619418716314
0.0008828619418716314
0.0008828619418716314
0.000871091095403193
0.000871091095403193
0.000871091095403193
0.0008232363272292931
0.0008232363272292931
0.0008232363272292931
0.0008222808145711048
0.0008222808145711048
0.0008222808145711048
0.0008205738803591024
0.0008205738803591024
0.0008205738803591024
0.0008189027838316996
0.0008189027838316996
0.0008189027838316996
0.0008151719733384375
0.0008151719733384375
0.0008151719733384375
0.0008076950603540785
0.0008076950603540785
0.0008076950603540785
0.0008000344342467945
0.0008000344342467945
0.0008000344342467945
0.000795033977195256
0.000795033977195256
0.000795033977195256
0.0007848850854924014
0.0007848850854924014
0.0007848850854924014
0.0007833934362090658
0.0007833934362090658
0.0007833934362090658
0.0007821124363623475
0.0007821124363623475
0.0007821124363623475
0.0007812337722716096
0.0007812337722716096
0.0007812337722716096
0.0007780860432477867
0.0007780860432477867
0.0007780860432477867
0.0007778239848230719
0.0007778239848230719
0.0007778239848230719
0.0007755188620679913
0.0007755188620679913
0.0007755188620679913
0.0007743691127089254
0.0007743691127089254
0.0007743691127089254
0.0007542379553525793
0.0007542379553525793
0.0007542379553525793
0.0007509459790839494
0.0007509459790839494
0.0007509459790839494
0.0007503656580249782
0.0007503656580249782
0.0007503656580249782
0.0007476133079294162
0.0007476133079294162
0.0007476133079294162
0.0007292010274957114
0.0007292010274957114
0.0007292010274957114
0.0007228701895993029
0.0007228701895993029
0.0007228701895993029
0.0007226216721221612
0.0007226216721221612
0.0007226216721221612
0.0007222519504066924
0.0007222519504066924
0.0007222519504066924
0.0007190124559613215
0.0007190124559613215
0.0007190124559613215
0.0007177458439735389
0.0007177458439735389
0.0007177458439735389
0.0007141625582687574
0.0007141625582687574
0.0007141625582687574
0.000713788198460867
0.000713788198460867
0.000713788198460867
0.0007129340558334588
0.0007129340558334588
0.0007129340558334588
0.0007112866926356895
0.0007112866926356895
0.0007112866926356895
])

# ╔═╡ 94a0f103-7157-4116-8d8f-b01e1cf3c099
function mass_calibration(
    dcm_array, center_large_LD, center, cal_rod_slice, rows, cols, spacing
)
    center_LD = center_large_LD
    dist_x = abs(center_LD[1] - center[1])
    dist_y = abs(center_LD[2] - center[2])

    if dist_x == 0
        mass_center_x = center[1]
        if center_LD[2] > center[2]
            mass_center_y = round(center[1] - round(23 / spacing[1], RoundUp), RoundUp)
        else
            mass_center_y = round(center[1] + round(23 / spacing[1], RoundUp), RoundUp)
        end
    elseif dist_y == 0
        mass_center_y = center[2]
        if center_LD[1] > center[1]
            mass_center_x = round(center[1] - round(23 / spacing[1], RoundUp), RoundUp)
        else
            mass_center_x = round(center[0] + round(23 / spacing[1], RoundUp), RoundUp)
        end

    else
        mass_angle = atan(dist_y / dist_x)
        dist_x = (23 / spacing[1]) * cos(mass_angle)
        dist_y = (23 / spacing[1]) * sin(mass_angle)

        if (center_LD[1] < center[1] && center_LD[2] < center[2])
            mass_center_x = round(center[1] + dist_x, RoundUp)
            mass_center_y = round(center[2] + dist_y, RoundUp)
        elseif (center_LD[1] < center[1] && center_LD[2] > center[2])
            mass_center_x = round(center[1] + dist_x, RoundUp)
            mass_center_y = round(center[2] - dist_y, RoundUp)
        elseif (center_LD[1] > center[1] && center_LD[2] < center[2])
            mass_center_x = round(center[1] - dist_x, RoundUp)
            mass_center_y = round(center[2] + dist_y, RoundUp)
        elseif (center_LD[1] > center[1] && center_LD[2] > center[2])
            mass_center_x = round(center[1] - dist_x, RoundUp)
            mass_center_y = round(center[2] - dist_y, RoundUp)
        end
    end

    mass_cal_center = [mass_center_y, mass_center_x]
    x_distance = abs(center[1] - mass_cal_center[2])
    angled_distance = sqrt(
        (center[1] - mass_cal_center[2])^2 + (center[2] - mass_cal_center[2])^2
    )
    angle_0_200HA = acos(x_distance / angled_distance) * 180 / π
    mask_0HU = PhantomSegmentation.create_circular_mask(
        cols, rows, mass_cal_center, Int(round(6.9 / spacing[1]))
    )
	# @info size(mask_0HU)
    masked_0HU = mask_0HU .* dcm_array[:, :, cal_rod_slice]
    nonzero_count = length(findall(x -> x != 0, masked_0HU))
    mean_0HU = sum(masked_0HU) / nonzero_count
    std_0HU = 0

    for voxel in vec(masked_0HU)
        if voxel != 0
            std_0HU += (voxel - mean_0HU)^2
        end
    end

    std_0HU = sqrt(std_0HU / nonzero_count)
    mask_200HU = PhantomSegmentation.create_circular_mask(
        cols, rows, (center[2], center[1]), Int(round(6.9 / spacing[1]))
    )
    masked_200HU = mask_200HU .* dcm_array[:, :, cal_rod_slice]
    nonzero_count_200HU = length(findall(x -> x != 0, masked_200HU))
    mean_200HU = sum(masked_200HU) / nonzero_count_200HU
    mass_cal_factor = 0.2 / (mean_200HU - mean_0HU)
    water_rod_metrics = mean_0HU, std_0HU

    return mass_cal_factor, angle_0_200HA, water_rod_metrics
end

# ╔═╡ f45c6429-333a-4096-ab9e-86df847d5e32
mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, center_large_LD, center_insert, 9, rows, cols, pixel_size)

# ╔═╡ 16edd9f1-dffa-44ec-8e86-eb63281a903c


# ╔═╡ 6cc4f6b2-0517-45a4-827f-8002a5da2fab


# ╔═╡ d6a497af-05cb-47a0-8562-976229ffa9d5


# ╔═╡ bda24523-ad64-4449-870d-71de2318020a


# ╔═╡ 977e0878-0a3c-4671-81e0-12ccfb7a034b


# ╔═╡ 8760d77b-821e-4755-b486-71000ecfaabc


# ╔═╡ 1b861fdf-4c35-4428-9e5d-64053f019fb0


# ╔═╡ cd5cb82d-3a09-4e65-b44a-05233cb6b7ad
md"""
# Score Inserts
"""

# ╔═╡ cb3910a5-f0ec-42f1-862a-86afe8381fbb
md"""
## Insert 1
"""

# ╔═╡ c622b37d-619f-41eb-9697-51ac6355a552
insert1_array = masked_array[:, :, slices1[1]-1:slices1[2]+1];

# ╔═╡ 9cc803fc-da75-4e3e-916f-3fb10a07063a
begin
	mask1_3D = Array{Bool}(undef, size(insert1_array))
	for z in 1:size(insert1_array, 3)
		mask1_3D[:, :, z] = mask1
	end
end;

# ╔═╡ d72c7049-ae03-4559-9948-9bada486e171
md"""
#### Dilated mask
"""

# ╔═╡ 8fadd7bf-1b5f-457d-9895-a8633c7ae12b
dilated_mask1 = dilate(dilate(mask1_3D));

# ╔═╡ d1eb5d2a-7f54-4635-87f0-3a3bcb9ee0c9
@bind g2 overlay_mask_bind(dilated_mask1)

# ╔═╡ e97e9f87-0650-4360-9096-1619257895b9
overlay_mask_plot(insert1_array, dilated_mask1, g2, "dilated mask")

# ╔═╡ 1621eeb9-183f-4f3b-b1ac-17fad532fe39
alg = Agatston()

# ╔═╡ 6d43f921-9b68-49a0-8ad9-2e040d452eb6
overlayed_mask1 = create_mask(insert1_array, mask1_3D);

# ╔═╡ df26ce2c-7e51-4966-9207-9c39501eea69
mass_cal_factor2 = 0.000749257

# ╔═╡ 56972bd7-50f7-4814-8208-ca3c5c55e63a
agat1, mass1 = score(overlayed_mask1, pixel_size, avg_mass_cals, alg)

# ╔═╡ 5c3a2b75-96c6-4f78-bb9e-ac246ecea3ea
md"""
## Insert 2
"""

# ╔═╡ e3f40211-409e-428a-b7bf-da211ff427f4
insert2_array = masked_array[:, :, slices2[1]-1:slices2[2]+1];

# ╔═╡ e08949ed-278c-49b1-ac64-a042da7985d8
begin
	mask2_3D = Array{Bool}(undef, size(insert2_array))
	for z in 1:size(insert2_array, 3)
		mask2_3D[:, :, z] = mask2
	end
end;

# ╔═╡ dd4cbaab-e8d7-4106-a4f3-ef0babff9bee
md"""
#### Dilated mask
"""

# ╔═╡ 6b3107de-b260-45da-97d9-e1f2bbf742a3
dilated_mask2 = dilate(dilate(mask2_3D));

# ╔═╡ cb9972c2-6be8-4a69-9c44-155f1224e90a
@bind g3 overlay_mask_bind(dilated_mask2)

# ╔═╡ 7fa966a2-7f2a-42b4-a7c1-b0c763db5680
overlay_mask_plot(insert2_array, dilated_mask2, g3, "dilated mask")

# ╔═╡ 565e0f0c-6ceb-4794-ab88-527df0093f50
overlayed_mask2 = create_mask(insert2_array, mask2_3D);

# ╔═╡ e06ae422-3b23-4cc6-b66a-46bcdb6eb2ea
agat2, mass2 = score(overlayed_mask2, pixel_size, avg_mass_cals, alg)

# ╔═╡ 4fb7e802-6b65-46ca-84a4-3b14e2f5856a
md"""
# Results
"""

# ╔═╡ e1b9e039-2a3d-4e78-95f5-9cbe5720d3d3
dens = [
	0.196
	0.380
	0.408
	0.800
] # mg/mm^3

# ╔═╡ 75c41d50-8fd3-45a7-8638-bf95a07f8d98
vol = 190.63 # mm^3

# ╔═╡ 624be6a8-d8b0-4a26-a27a-572b450c4313
masses = vol .* dens

# ╔═╡ 71318d50-c416-4308-b8d0-2e22d76c9c71
inserts = [
	"insert 1",
	"insert 2"
]

# ╔═╡ 8f7c6f13-8aff-4c9d-b5ab-79b3302d45fb
if occursin("rod1", scan)
	rod1 = true
	ground_truth_mass = masses[3:4]
elseif occursin("rod2", scan)
	rod1 = false
	ground_truth_mass = masses[1:2]
end

# ╔═╡ 882d2ac2-75dc-448c-904c-72aaa62af50b
calculated_mass = [
	mass2,
	mass1,
]

# ╔═╡ 198f49e3-06fe-4047-a478-5642a66b26d2
calculated_agat = [
	agat1,
	agat2,
]

# ╔═╡ 05af3e51-0e69-40dc-a087-9c9a342af509
pth

# ╔═╡ 2b0a45a1-81cf-49ff-bee1-893e14a51dbc
if rod1
	rod = "rod1"
else
	rod = "rod2"
end

# ╔═╡ 1f2b43d8-b3e3-41c4-8b2c-c22a46a52445
df = DataFrame(
	vender = VENDER,
	scan = scan,
	rod = rod,
	calculated_agat = calculated_agat,
	ground_truth_mass = ground_truth_mass,
	calculated_mass = calculated_mass,
)

# ╔═╡ 10034f35-aeeb-48b5-913b-5eeb747ffb31
densities = dens .* 10e2

# ╔═╡ c86b720e-56c1-4ad9-93d6-68a3607a8589
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

# ╔═╡ a02d4461-0740-485d-8ff7-f26e91c08ff9
md"""
## Save Results
"""

# ╔═╡ faa2e9ef-2433-4283-b373-3632ba524589
# if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER))
# 	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER))
# end

# ╔═╡ 25a170e8-ec5c-49c8-bb2e-afd86ec6c929
# output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "/", scan, ".csv")

# ╔═╡ 1f77205e-54b3-4a61-b823-3f6b9ca303b3
# CSV.write(output_path, df)

# ╔═╡ 762b3014-f0a8-409e-8499-ba6e0cd4ecaa
md"""
## Save full df
"""

# ╔═╡ 2284c1af-d098-498c-9b75-b218985e0e6d
dfs = []

# ╔═╡ 3b3e773d-b426-4edd-a32d-1662332d2332
push!(dfs, df)

# ╔═╡ 08320750-5c21-4554-96e5-bc98d6c3d812
# if length(dfs) == 10
# 	global new_df = vcat(dfs[1:10]...)
# 	output_path_new = string(cd(pwd, "..") , "/data/output/", VENDER, "/full.csv")
# 	CSV.write(output_path_new, new_df)
# end

# ╔═╡ Cell order:
# ╠═795a7fe3-22ea-4367-b7db-3941fd9a4888
# ╠═a1d626e2-fc5b-484e-9838-5a8463d4617a
# ╟─2bfa5dfc-472d-40fa-86e0-754a3e8abb38
# ╠═f16eecad-68d2-42ae-a4dc-9b574984fd26
# ╟─bb021d41-db85-4311-b16d-e4303bc38f63
# ╠═d574f97c-bad2-413c-ac73-18363cc559fe
# ╠═e0467b68-aecf-4ff8-a221-59b8d7d1a10a
# ╠═75d8379a-1ddb-4127-a4f6-bc88d6ce6d14
# ╠═ca7fd376-b7eb-43cd-8b74-dec6431d1823
# ╠═4dee1592-8848-4fb7-801f-c56404018046
# ╠═ea20f143-6cfe-4662-bddb-97c5a15fe81f
# ╟─cb24aac3-1ce5-40b6-a843-00de6c998f2d
# ╟─1eca10c3-2259-4fec-a4e4-f52b7be45b02
# ╟─2c72f868-dda1-4baa-a200-511cbaa2755a
# ╟─c5e319a1-ac8b-4595-b249-b94815116c20
# ╟─d068bbe8-f2e4-4697-bcfd-b0a9333a7957
# ╠═6d646443-5d78-4711-baac-25bf60fe1515
# ╟─15897036-6264-4451-a1b1-5e42c7c54a8e
# ╠═97f143cf-4725-493d-a6c6-a9693baddf3b
# ╠═36c46879-c037-4236-b106-45012797f1d2
# ╠═55dcd1d5-1cdb-41be-9c05-66395905a880
# ╠═c73647f9-489b-4af4-ab78-3e2c85925a11
# ╟─b49c9980-6648-49b0-8a2f-82b2e180cb1d
# ╠═330b7534-eefa-4946-a33c-efdac78bb213
# ╠═7d9a6306-38ed-467d-b55c-a70ce9f8242a
# ╠═61268b09-b7db-44d7-8e75-47f7fe4ff5e9
# ╠═34710f8d-822c-49a1-a4f9-4cf42a0a00e8
# ╟─da315970-7db8-4070-9db2-e3d5b21b3a84
# ╠═e4334739-f94b-473f-88bd-d58767fb9fb0
# ╠═fe9dc076-bf45-44d4-8b6a-a74b839dfd98
# ╠═9d024419-81c2-42c0-8615-918f9a5ba2a9
# ╠═fb9de4d3-4e19-4c87-80e3-4086908ca48a
# ╠═7a9f5566-0439-4005-a26b-e992613480b1
# ╠═94a0f103-7157-4116-8d8f-b01e1cf3c099
# ╠═f45c6429-333a-4096-ab9e-86df847d5e32
# ╠═16edd9f1-dffa-44ec-8e86-eb63281a903c
# ╠═6cc4f6b2-0517-45a4-827f-8002a5da2fab
# ╠═d6a497af-05cb-47a0-8562-976229ffa9d5
# ╠═bda24523-ad64-4449-870d-71de2318020a
# ╠═977e0878-0a3c-4671-81e0-12ccfb7a034b
# ╠═8760d77b-821e-4755-b486-71000ecfaabc
# ╠═1b861fdf-4c35-4428-9e5d-64053f019fb0
# ╟─cd5cb82d-3a09-4e65-b44a-05233cb6b7ad
# ╟─cb3910a5-f0ec-42f1-862a-86afe8381fbb
# ╠═c622b37d-619f-41eb-9697-51ac6355a552
# ╠═9cc803fc-da75-4e3e-916f-3fb10a07063a
# ╟─d72c7049-ae03-4559-9948-9bada486e171
# ╠═8fadd7bf-1b5f-457d-9895-a8633c7ae12b
# ╟─d1eb5d2a-7f54-4635-87f0-3a3bcb9ee0c9
# ╠═e97e9f87-0650-4360-9096-1619257895b9
# ╠═1621eeb9-183f-4f3b-b1ac-17fad532fe39
# ╠═6d43f921-9b68-49a0-8ad9-2e040d452eb6
# ╠═df26ce2c-7e51-4966-9207-9c39501eea69
# ╠═56972bd7-50f7-4814-8208-ca3c5c55e63a
# ╟─5c3a2b75-96c6-4f78-bb9e-ac246ecea3ea
# ╠═e3f40211-409e-428a-b7bf-da211ff427f4
# ╠═e08949ed-278c-49b1-ac64-a042da7985d8
# ╟─dd4cbaab-e8d7-4106-a4f3-ef0babff9bee
# ╠═6b3107de-b260-45da-97d9-e1f2bbf742a3
# ╟─cb9972c2-6be8-4a69-9c44-155f1224e90a
# ╠═7fa966a2-7f2a-42b4-a7c1-b0c763db5680
# ╠═565e0f0c-6ceb-4794-ab88-527df0093f50
# ╠═e06ae422-3b23-4cc6-b66a-46bcdb6eb2ea
# ╟─4fb7e802-6b65-46ca-84a4-3b14e2f5856a
# ╠═e1b9e039-2a3d-4e78-95f5-9cbe5720d3d3
# ╠═75c41d50-8fd3-45a7-8638-bf95a07f8d98
# ╠═624be6a8-d8b0-4a26-a27a-572b450c4313
# ╠═71318d50-c416-4308-b8d0-2e22d76c9c71
# ╠═8f7c6f13-8aff-4c9d-b5ab-79b3302d45fb
# ╠═882d2ac2-75dc-448c-904c-72aaa62af50b
# ╠═198f49e3-06fe-4047-a478-5642a66b26d2
# ╠═05af3e51-0e69-40dc-a087-9c9a342af509
# ╠═2b0a45a1-81cf-49ff-bee1-893e14a51dbc
# ╠═1f2b43d8-b3e3-41c4-8b2c-c22a46a52445
# ╠═10034f35-aeeb-48b5-913b-5eeb747ffb31
# ╟─c86b720e-56c1-4ad9-93d6-68a3607a8589
# ╟─a02d4461-0740-485d-8ff7-f26e91c08ff9
# ╠═faa2e9ef-2433-4283-b373-3632ba524589
# ╠═25a170e8-ec5c-49c8-bb2e-afd86ec6c929
# ╠═1f77205e-54b3-4a61-b823-3f6b9ca303b3
# ╟─762b3014-f0a8-409e-8499-ba6e0cd4ecaa
# ╠═2284c1af-d098-498c-9b75-b218985e0e6d
# ╠═3b3e773d-b426-4edd-a32d-1662332d2332
# ╠═08320750-5c21-4554-96e5-bc98d6c3d812
