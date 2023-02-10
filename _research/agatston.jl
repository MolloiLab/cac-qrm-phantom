### A Pluto.jl notebook ###
# v0.19.11

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

# ╔═╡ 651419e4-4c13-4ddd-8e69-5cdf2b995a51
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("Statistics")
		Pkg.add("StatsBase")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageFiltering")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("GLM")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI
	using CairoMakie
	using Statistics
	using StatsBase: quantile!
	using ImageMorphology
	using ImageFiltering
	using CSV
	using DataFrames
	using GLM
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
	using CalciumScoring
end

# ╔═╡ d9401636-eb09-45b2-9663-b8ab615531f9
TableOfContents()

# ╔═╡ ceb200df-c17b-476f-ae95-a13e572a243e
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 35c2524a-28a4-49eb-9fe2-e4232084af1a
begin
	SCAN_NUMBER = 7
	VENDER = "Siemens_SOMATOM_Force_extra"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ 7c5f2ec6-26bb-4a2a-9dd1-a7f3b0de10b0
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 7dab92d6-13c3-4c1e-8b51-ad80a709a80b
root_path = string(BASE_PATH, VENDER)

# ╔═╡ 9bcfc64a-170b-461c-baa1-ac5c1ab73579
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 041192c8-875e-43f3-881e-dc993955c91f
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 57cd84bb-9307-44eb-87c1-9a663e10e5c2
scan = basename(pth)

# ╔═╡ 807cfdd5-379a-4c65-ab15-534f1e28b65c
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 948b9234-f845-48f5-a246-18cbc961ad67
md"""
## Helper Functions
"""

# ╔═╡ baa55458-c8e2-4852-b486-ef8c854f77bf
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 7434ae0d-06f0-44fc-b9e2-198fb9afb97c
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 7325a8d4-6e0a-4087-ad64-d0874fb26aad
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

# ╔═╡ 63c84140-6510-4186-a0bf-dab85d19e2cc
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ a9516b63-80c0-42ff-9bac-a9ad311fbe9c
md"""
## Segment Heart
"""

# ╔═╡ fb4542a8-3406-48bb-b255-ebde03adb951
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 835490d8-332c-400a-a2f7-2dcff791d3d6
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ a2172056-1f1c-40c9-ac6f-816c80dd09ef
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ 38042368-0be5-4c64-a4e6-c596201cdd2e
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig
end

# ╔═╡ a467ddd4-e2d3-445e-808b-298fea9b927c
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig2
end

# ╔═╡ dfb18ee0-7aa7-43c3-bf7e-914b2ee4b729
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2]+1, center_insert[1]:center_insert[1]+1, markersize=10, color=:red)
	fig3
end

# ╔═╡ 65b0a9dc-efc2-4eb7-9995-1b8bd8113a63
md"""
## Segment Calcium Rod
"""

# ╔═╡ 8ae2afdc-ae49-49e6-ba59-d20373907ba7
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ d635ed56-ef46-428a-ac7e-f2fc5c10e4d6
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ e4052629-fe1c-4b9f-8d74-fe7dce34de06
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ 630ec973-b16d-46b3-bb3b-229d407087b1
cal_rod_slice

# ╔═╡ bdb61513-e02a-43df-be3e-164f415e061b
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 5ec27948-098e-4acf-961a-0bc36307f807
mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
            dcm_array, masked_array, header, slice_CCI, center_insert
);

# ╔═╡ 480f30c9-b0cb-4f49-8767-717d6213902e
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 02cd946e-98a3-4043-b482-7a19dcd1f157
heatmap(masks, colormap=:grays)

# ╔═╡ 84f78e6d-c73f-4e74-bfb2-1de22990d1c8
md"""
## Mass cal factor
"""

# ╔═╡ fe8c8d86-816d-46f8-9068-663c4a45ce97
begin
	output = calc_output(masked_array, header, slice_CCI, 130, trues(3, 3))
    insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
end

# ╔═╡ 0604f5cf-a731-405e-9c3c-a9f918352f96
center_large_LD = insert_centers[:Large_LD]

# ╔═╡ fd78b7f0-a83e-4365-a770-16483010376a
rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])

# ╔═╡ 8e1ca1c1-ea27-4f0e-89cd-b466053e9485
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 3c53f529-9fe6-47e9-86f9-41de9ee92b87
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
	@info size(mask_0HU)
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

# ╔═╡ 64d8f848-cc5b-4b51-90bd-aacca098593b
mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, center_large_LD, center_insert, cal_rod_slice, rows, cols, pixel_size)

# ╔═╡ 11e23023-d88b-4509-a562-51d5255c6b76
md"""
# Score Large Inserts
"""

# ╔═╡ 9c14d354-9f0d-45d5-8dc2-98a719b0465c
arr = masked_array[:, :, slice_CCI-2:slice_CCI+2];

# ╔═╡ 5ee01db8-2a9d-491d-8305-adff24c6dd09
single_arr = masked_array[:, :, slice_CCI];

# ╔═╡ c0f2a223-d737-4823-80c0-76801ae26745
md"""
## High Density
"""

# ╔═╡ c986c5a2-dffb-4200-af4b-fde0e94c5361
begin
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
end;

# ╔═╡ cc598eba-d013-4984-aa5a-6f94f8f03f19
md"""
#### Dilated mask
"""

# ╔═╡ 630c79f5-16c3-4ce7-b197-16348783c2a0
dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D));

# ╔═╡ cf80857d-f58d-46a8-b7fc-f44dfbd8d2a6
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 5e82dcc1-431f-4ba7-9adf-df62beaa0ca9
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 53b8bada-351b-42d1-8807-79ed375b8186
alg = Agatston()

# ╔═╡ dc41660b-6620-4cd6-859b-b5fc3acc7063
overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD);

# ╔═╡ fe4f4945-d3e8-48c8-a83b-f35554de86e3
heatmap(overlayed_mask_l_hd[:, :, 3])

# ╔═╡ 4cd5d9bd-1b56-4cd4-b46c-168cb3da6b4d
agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg)

# ╔═╡ c164f2e6-edea-45f6-8cde-4da1ba1820a0
md"""
## Medium Density
"""

# ╔═╡ c1ba29be-fad6-40f5-a828-122ea002ee74
begin
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
end;

# ╔═╡ 4b2fe0ac-f265-4a29-9ef0-b452e5a66d22
md"""
#### Dilated mask
"""

# ╔═╡ 4493b19e-bfb7-42da-9d5d-457b30707f14
dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D));

# ╔═╡ 52242114-a86b-43ff-a10b-580805e48b27
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 3e627949-fe0f-4f3a-8f94-9b8b8eb04736
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ 9049f063-bcc7-4ddf-90db-adf813b2e645
overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD);

# ╔═╡ a6810274-a229-4e81-ba76-86a116c2cec1
agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_cal_factor, alg)

# ╔═╡ 03d7f139-8283-42fc-b93c-cf52a2f1a171
md"""
## Low Density
"""

# ╔═╡ 04511af5-6cba-4d2e-83c5-423347e9d12a
begin
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
end;

# ╔═╡ b2a7843a-08e1-4fb7-95c7-4c30c66b2403
md"""
#### Dilated mask
"""

# ╔═╡ 149d69f2-2280-4925-8c73-47499a6cbf01
dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D));

# ╔═╡ 797cab84-8ea9-4ec0-8ce9-13e0de6371ba
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ d1cd6b6f-fd66-4255-aa30-eb9a2ba08809
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ 063add05-86dc-44a6-8b7a-0eba4f04da3f
overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD);

# ╔═╡ 9c8e037d-f157-4b19-9e91-88ad2a219e0e
agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg)

# ╔═╡ 4633f872-8ef1-42fd-baa8-3b7a32583dd7
md"""
# Score Medium Inserts
"""

# ╔═╡ d443e864-7eaa-4f4a-b2b9-642c8f7bf165
md"""
## High Density
"""

# ╔═╡ 8e96f01a-5cd3-4dde-a926-d81f851900fa
begin
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
end;

# ╔═╡ 7155e9f1-73a2-4d37-a7f3-92d1676f8058
md"""
#### Dilated mask
"""

# ╔═╡ 91820ed4-6a31-4132-92f3-b138fdebc482
dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 5fe3c5f4-ed73-4a35-b71a-dd8829b59da4
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ 52c7be1a-2daa-4b7c-89c4-d72463535aa9
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 74636d4e-aaef-4fa0-8c55-5afc49418397
overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD);

# ╔═╡ 3fefa81c-89d4-416e-be01-174a8a543dd6
agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg)

# ╔═╡ f7dda15e-b22d-4a33-9efa-da167abe6538
md"""
## Medium Density
"""

# ╔═╡ 1bea5ec6-a17a-4640-89ea-c6666499b2fc
begin
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
end;

# ╔═╡ da9714c0-fead-4343-a50c-a8a2c7f1d9f4
md"""
#### Dilated mask
"""

# ╔═╡ d602f6af-0130-4a51-be9e-465627f6730d
dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))));

# ╔═╡ e36c2ef8-260b-418f-b8a8-52ebf731ebe3
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ efbd8bf7-ce73-463b-a4ca-820f7213b3ee
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 5e154bdf-84b1-4c97-bc71-1c632ce8d5bd
overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD);

# ╔═╡ 7c9ca207-83b0-48e9-982f-f2408a956d9c
agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_cal_factor, alg)

# ╔═╡ f8e48d57-c674-4915-8c79-589bde4cc677
md"""
## Low Density
"""

# ╔═╡ c6bbd448-0e32-4eeb-b7d4-f9718ec18304
begin
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
end;

# ╔═╡ 79fea025-1886-4929-8d96-d9f0505c432d
md"""
#### Dilated mask
"""

# ╔═╡ 28ba6334-8adb-4a55-a48d-19f8e4f5f728
dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))));

# ╔═╡ 81ed8294-fa65-4164-86a8-1c219ab321ef
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 2c0861f2-7d97-4b75-a889-c66e27b09f13
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ 7bfd57f2-2027-49cb-abab-8ef11409645c
overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD);

# ╔═╡ f4f21b51-6b56-4b97-812e-b9b095c22505
agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg)

# ╔═╡ f02905f4-f15f-4876-b26e-47d70b9c488e
md"""
# Score Small Inserts
"""

# ╔═╡ 8df6e736-67af-4fdf-8eeb-18680b703d8a
md"""
## High Density
"""

# ╔═╡ d9695ec4-ad31-47ae-9baa-aa411ac46e1c
begin
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
end;

# ╔═╡ e574c7b9-09c6-4fc5-820c-78ceaf5a3c80
md"""
#### Dilated mask
"""

# ╔═╡ f0abacd4-eb6b-494d-ac93-342c1e29bd91
dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))));

# ╔═╡ 5b43337b-1ad0-4640-b285-1cd0bd032410
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ 5072be8f-34ed-46e7-aad3-9315525f7232
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 9bc5292c-d669-46b4-8827-8ff950e44c6c
overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD);

# ╔═╡ e9e619f8-5b94-45e9-80d3-cdf8faa4523a
agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg)

# ╔═╡ 49d9fb49-5a12-4b14-b759-5cb6263eafed
md"""
## Medium Density
"""

# ╔═╡ 551eb59d-79d5-4098-93ff-c624fee56119
begin
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
end;

# ╔═╡ 31803d25-bc2f-44b3-84b2-88291d2138ae
md"""
#### Dilated mask
"""

# ╔═╡ c91a36be-99bc-495a-9969-954918164375
dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))));

# ╔═╡ b4af4774-916b-4712-bdf8-811c289d316c
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ e66946b7-fcd4-4cfd-9500-ca0b0ececf26
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ fa8ce4b9-898a-4706-97c1-9398b15779ea
overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD);

# ╔═╡ 70581d09-b3eb-4c79-bc83-89945b22881e
agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_cal_factor, alg)

# ╔═╡ 16585130-17c0-46bd-9b35-c55dad77819f
md"""
## Low Density
"""

# ╔═╡ 98b22d21-4a54-4b2b-9ca7-7f32c0442195
begin
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in 1:size(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
end;

# ╔═╡ e45981da-de12-4484-9a07-f6ca3f1c9cf1
md"""
#### Dilated mask
"""

# ╔═╡ 2cff2e76-a209-45d1-8970-bb7b101b91c0
dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))));

# ╔═╡ 3495f22c-8883-48ef-81d2-fcc670b654be
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 8856ad47-91ac-4f59-89a5-26e7106f16c7
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ c3b274e2-7bd5-440e-a977-c29e1c20b5ce
overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD);

# ╔═╡ cb920bb8-7c28-48d5-af49-8f94caa61f94
# agat_s_ld = score(overlayed_mask_s_ld, pixel_size, alg)

# ╔═╡ 996d788e-9d3e-46a6-9532-59002db192d0
agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg)

# ╔═╡ be6ac424-cd70-4dda-99ee-6c08c131a56b
md"""
# Results
"""

# ╔═╡ 3d021782-5514-4266-bfc5-a42b4ab5ec66
density_array = [0, 200, 400, 800]

# ╔═╡ cdf5731b-091f-47d9-a992-8a3c5ef3f878
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ fa458085-3a0c-4d8d-b3d3-9a3cb6f4b36f
md"""
## Agatston
"""

# ╔═╡ 59474596-ba9c-470a-8007-c656623b07fb
calculated_agat_large = [
	agat_l_ld,
	agat_l_md,
	agat_l_hd
]

# ╔═╡ b3909180-04f4-471a-b1e3-457c5431bc19
calculated_agat_medium = [
	agat_m_ld,
	agat_m_md,
	agat_m_hd
]

# ╔═╡ 3dfc2f6e-7ba2-4482-9240-8a6800fa4a2d
calculated_agat_small = [
	agat_s_ld,
	agat_s_md,
	agat_s_hd
]

# ╔═╡ 6652c87e-44dc-4a4b-b291-c5978c97a616
md"""
## Mass
"""

# ╔═╡ 9968b6e4-e852-4766-a30a-43d23be32b7b
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ b8bf59f2-24fa-43fe-8260-e85ed3028b93
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ a48c296d-1e1b-4ecd-89d9-645577bc877d
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 84b7fdd1-76cc-4556-bf99-b93b1d93d555
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ da4af805-291d-4a52-b3c8-f1f8d98bc1d8
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ 749a5fe6-d063-4728-ab2b-530e694a9f17
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 35175bc9-1fec-4bd2-b9a8-381b69996f4a
df = DataFrame(
	inserts = inserts,
	calculated_agat_large = calculated_agat_large,
	calculated_agat_medium = calculated_agat_medium,
	calculated_agat_small = calculated_agat_small,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small,
	mass_cal_factor = mass_cal_factor
)

# ╔═╡ bcacf403-5859-4b30-b8f5-0106ca01cd2c
begin
	fmass22 = Figure()
	axmass22 = Axis(fmass22[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass22.title = "Mass Measurements (Large)"
	axmass22.ylabel = "Mass (mg)"
	axmass22.xlabel = "Density (mg/cm^3)"

	xlims!(axmass22, 0, 850)
	ylims!(axmass22, 0, 100)
	
	fmass22[1, 2] = Legend(fmass22, axmass22, framevisible = false)
	
	fmass22
end

# ╔═╡ b83fc531-4867-4056-bd54-af8eeab50496
begin
	fmass32 = Figure()
	axmass32 = Axis(fmass32[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass32.title = "Mass Measurements (Medium)"
	axmass32.ylabel = "Mass (mg)"
	axmass32.xlabel = "Density (mg/cm^3)"

	xlims!(axmass32, 0, 850)
	ylims!(axmass32, 0, 25)
	
	fmass32[1, 2] = Legend(fmass32, axmass32, framevisible = false)
	
	fmass32
end

# ╔═╡ 3d20bc0e-9156-41d9-93b3-d044ce50b0ef
begin
	fmass42 = Figure()
	axmass42 = Axis(fmass42[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass42.title = "Mass Measurements (Small)"
	axmass42.ylabel = "Mass (mg)"
	axmass42.xlabel = "Density (mg/cm^3)"

	xlims!(axmass42, 0, 850)
	ylims!(axmass42, 0, 1.5)
	
	fmass42[1, 2] = Legend(fmass42, axmass42, framevisible = false)
	
	fmass42
end

# ╔═╡ c9e8865a-9143-4a5d-8b0c-aae193168d0d
md"""
### Save Results
"""

# ╔═╡ 2629fd24-e10e-4c42-a388-9f4b680bf550
if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER, "_agatston"))
	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER, "_agatston",))
end

# ╔═╡ 618cf15a-4ccc-48a3-b122-e346d99d0635
output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "_agatston", "/", scan, ".csv")

# ╔═╡ 62d446bf-6e50-4a6b-aa3c-f8e97e99490f
CSV.write(output_path, df)

# ╔═╡ Cell order:
# ╠═651419e4-4c13-4ddd-8e69-5cdf2b995a51
# ╠═d9401636-eb09-45b2-9663-b8ab615531f9
# ╟─ceb200df-c17b-476f-ae95-a13e572a243e
# ╠═35c2524a-28a4-49eb-9fe2-e4232084af1a
# ╟─7c5f2ec6-26bb-4a2a-9dd1-a7f3b0de10b0
# ╠═7dab92d6-13c3-4c1e-8b51-ad80a709a80b
# ╠═9bcfc64a-170b-461c-baa1-ac5c1ab73579
# ╠═041192c8-875e-43f3-881e-dc993955c91f
# ╠═57cd84bb-9307-44eb-87c1-9a663e10e5c2
# ╠═807cfdd5-379a-4c65-ab15-534f1e28b65c
# ╟─948b9234-f845-48f5-a246-18cbc961ad67
# ╟─baa55458-c8e2-4852-b486-ef8c854f77bf
# ╟─7434ae0d-06f0-44fc-b9e2-198fb9afb97c
# ╟─7325a8d4-6e0a-4087-ad64-d0874fb26aad
# ╟─63c84140-6510-4186-a0bf-dab85d19e2cc
# ╟─a9516b63-80c0-42ff-9bac-a9ad311fbe9c
# ╠═fb4542a8-3406-48bb-b255-ebde03adb951
# ╟─835490d8-332c-400a-a2f7-2dcff791d3d6
# ╟─a2172056-1f1c-40c9-ac6f-816c80dd09ef
# ╟─38042368-0be5-4c64-a4e6-c596201cdd2e
# ╟─a467ddd4-e2d3-445e-808b-298fea9b927c
# ╟─dfb18ee0-7aa7-43c3-bf7e-914b2ee4b729
# ╟─65b0a9dc-efc2-4eb7-9995-1b8bd8113a63
# ╠═8ae2afdc-ae49-49e6-ba59-d20373907ba7
# ╟─d635ed56-ef46-428a-ac7e-f2fc5c10e4d6
# ╠═e4052629-fe1c-4b9f-8d74-fe7dce34de06
# ╠═630ec973-b16d-46b3-bb3b-229d407087b1
# ╟─bdb61513-e02a-43df-be3e-164f415e061b
# ╠═5ec27948-098e-4acf-961a-0bc36307f807
# ╠═480f30c9-b0cb-4f49-8767-717d6213902e
# ╠═02cd946e-98a3-4043-b482-7a19dcd1f157
# ╟─84f78e6d-c73f-4e74-bfb2-1de22990d1c8
# ╠═fe8c8d86-816d-46f8-9068-663c4a45ce97
# ╠═0604f5cf-a731-405e-9c3c-a9f918352f96
# ╠═fd78b7f0-a83e-4365-a770-16483010376a
# ╠═8e1ca1c1-ea27-4f0e-89cd-b466053e9485
# ╠═3c53f529-9fe6-47e9-86f9-41de9ee92b87
# ╠═64d8f848-cc5b-4b51-90bd-aacca098593b
# ╟─11e23023-d88b-4509-a562-51d5255c6b76
# ╠═9c14d354-9f0d-45d5-8dc2-98a719b0465c
# ╠═5ee01db8-2a9d-491d-8305-adff24c6dd09
# ╟─c0f2a223-d737-4823-80c0-76801ae26745
# ╠═c986c5a2-dffb-4200-af4b-fde0e94c5361
# ╟─cc598eba-d013-4984-aa5a-6f94f8f03f19
# ╠═630c79f5-16c3-4ce7-b197-16348783c2a0
# ╟─cf80857d-f58d-46a8-b7fc-f44dfbd8d2a6
# ╠═5e82dcc1-431f-4ba7-9adf-df62beaa0ca9
# ╠═53b8bada-351b-42d1-8807-79ed375b8186
# ╠═dc41660b-6620-4cd6-859b-b5fc3acc7063
# ╠═fe4f4945-d3e8-48c8-a83b-f35554de86e3
# ╠═4cd5d9bd-1b56-4cd4-b46c-168cb3da6b4d
# ╟─c164f2e6-edea-45f6-8cde-4da1ba1820a0
# ╠═c1ba29be-fad6-40f5-a828-122ea002ee74
# ╟─4b2fe0ac-f265-4a29-9ef0-b452e5a66d22
# ╠═4493b19e-bfb7-42da-9d5d-457b30707f14
# ╟─52242114-a86b-43ff-a10b-580805e48b27
# ╠═3e627949-fe0f-4f3a-8f94-9b8b8eb04736
# ╠═9049f063-bcc7-4ddf-90db-adf813b2e645
# ╠═a6810274-a229-4e81-ba76-86a116c2cec1
# ╟─03d7f139-8283-42fc-b93c-cf52a2f1a171
# ╠═04511af5-6cba-4d2e-83c5-423347e9d12a
# ╟─b2a7843a-08e1-4fb7-95c7-4c30c66b2403
# ╠═149d69f2-2280-4925-8c73-47499a6cbf01
# ╟─797cab84-8ea9-4ec0-8ce9-13e0de6371ba
# ╠═d1cd6b6f-fd66-4255-aa30-eb9a2ba08809
# ╠═063add05-86dc-44a6-8b7a-0eba4f04da3f
# ╠═9c8e037d-f157-4b19-9e91-88ad2a219e0e
# ╟─4633f872-8ef1-42fd-baa8-3b7a32583dd7
# ╟─d443e864-7eaa-4f4a-b2b9-642c8f7bf165
# ╠═8e96f01a-5cd3-4dde-a926-d81f851900fa
# ╟─7155e9f1-73a2-4d37-a7f3-92d1676f8058
# ╠═91820ed4-6a31-4132-92f3-b138fdebc482
# ╟─5fe3c5f4-ed73-4a35-b71a-dd8829b59da4
# ╠═52c7be1a-2daa-4b7c-89c4-d72463535aa9
# ╠═74636d4e-aaef-4fa0-8c55-5afc49418397
# ╠═3fefa81c-89d4-416e-be01-174a8a543dd6
# ╟─f7dda15e-b22d-4a33-9efa-da167abe6538
# ╠═1bea5ec6-a17a-4640-89ea-c6666499b2fc
# ╟─da9714c0-fead-4343-a50c-a8a2c7f1d9f4
# ╠═d602f6af-0130-4a51-be9e-465627f6730d
# ╟─e36c2ef8-260b-418f-b8a8-52ebf731ebe3
# ╠═efbd8bf7-ce73-463b-a4ca-820f7213b3ee
# ╠═5e154bdf-84b1-4c97-bc71-1c632ce8d5bd
# ╠═7c9ca207-83b0-48e9-982f-f2408a956d9c
# ╟─f8e48d57-c674-4915-8c79-589bde4cc677
# ╠═c6bbd448-0e32-4eeb-b7d4-f9718ec18304
# ╟─79fea025-1886-4929-8d96-d9f0505c432d
# ╠═28ba6334-8adb-4a55-a48d-19f8e4f5f728
# ╟─81ed8294-fa65-4164-86a8-1c219ab321ef
# ╠═2c0861f2-7d97-4b75-a889-c66e27b09f13
# ╠═7bfd57f2-2027-49cb-abab-8ef11409645c
# ╠═f4f21b51-6b56-4b97-812e-b9b095c22505
# ╟─f02905f4-f15f-4876-b26e-47d70b9c488e
# ╟─8df6e736-67af-4fdf-8eeb-18680b703d8a
# ╠═d9695ec4-ad31-47ae-9baa-aa411ac46e1c
# ╟─e574c7b9-09c6-4fc5-820c-78ceaf5a3c80
# ╠═f0abacd4-eb6b-494d-ac93-342c1e29bd91
# ╟─5b43337b-1ad0-4640-b285-1cd0bd032410
# ╠═5072be8f-34ed-46e7-aad3-9315525f7232
# ╠═9bc5292c-d669-46b4-8827-8ff950e44c6c
# ╠═e9e619f8-5b94-45e9-80d3-cdf8faa4523a
# ╟─49d9fb49-5a12-4b14-b759-5cb6263eafed
# ╠═551eb59d-79d5-4098-93ff-c624fee56119
# ╟─31803d25-bc2f-44b3-84b2-88291d2138ae
# ╠═c91a36be-99bc-495a-9969-954918164375
# ╟─b4af4774-916b-4712-bdf8-811c289d316c
# ╠═e66946b7-fcd4-4cfd-9500-ca0b0ececf26
# ╠═fa8ce4b9-898a-4706-97c1-9398b15779ea
# ╠═70581d09-b3eb-4c79-bc83-89945b22881e
# ╟─16585130-17c0-46bd-9b35-c55dad77819f
# ╠═98b22d21-4a54-4b2b-9ca7-7f32c0442195
# ╟─e45981da-de12-4484-9a07-f6ca3f1c9cf1
# ╠═2cff2e76-a209-45d1-8970-bb7b101b91c0
# ╟─3495f22c-8883-48ef-81d2-fcc670b654be
# ╠═8856ad47-91ac-4f59-89a5-26e7106f16c7
# ╠═c3b274e2-7bd5-440e-a977-c29e1c20b5ce
# ╠═cb920bb8-7c28-48d5-af49-8f94caa61f94
# ╠═996d788e-9d3e-46a6-9532-59002db192d0
# ╟─be6ac424-cd70-4dda-99ee-6c08c131a56b
# ╠═3d021782-5514-4266-bfc5-a42b4ab5ec66
# ╠═cdf5731b-091f-47d9-a992-8a3c5ef3f878
# ╟─fa458085-3a0c-4d8d-b3d3-9a3cb6f4b36f
# ╠═59474596-ba9c-470a-8007-c656623b07fb
# ╠═b3909180-04f4-471a-b1e3-457c5431bc19
# ╠═3dfc2f6e-7ba2-4482-9240-8a6800fa4a2d
# ╠═35175bc9-1fec-4bd2-b9a8-381b69996f4a
# ╟─6652c87e-44dc-4a4b-b291-c5978c97a616
# ╠═9968b6e4-e852-4766-a30a-43d23be32b7b
# ╠═b8bf59f2-24fa-43fe-8260-e85ed3028b93
# ╠═a48c296d-1e1b-4ecd-89d9-645577bc877d
# ╠═84b7fdd1-76cc-4556-bf99-b93b1d93d555
# ╠═da4af805-291d-4a52-b3c8-f1f8d98bc1d8
# ╠═749a5fe6-d063-4728-ab2b-530e694a9f17
# ╟─bcacf403-5859-4b30-b8f5-0106ca01cd2c
# ╟─b83fc531-4867-4056-bd54-af8eeab50496
# ╟─3d20bc0e-9156-41d9-93b3-d044ce50b0ef
# ╟─c9e8865a-9143-4a5d-8b0c-aae193168d0d
# ╠═2629fd24-e10e-4c42-a388-9f4b680bf550
# ╠═618cf15a-4ccc-48a3-b122-e346d99d0635
# ╠═62d446bf-6e50-4a6b-aa3c-f8e97e99490f
