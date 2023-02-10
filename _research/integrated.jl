### A Pluto.jl notebook ###
# v0.19.18

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

# ╔═╡ 3c162d31-407a-4c50-a8a7-30153192f207
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, CairoMakie
    using StatsBase: quantile!

end

# ╔═╡ 0aeb13cf-6299-473a-908e-5b0be1b23c86
TableOfContents()

# ╔═╡ c4bba726-9f22-4f76-a4ae-966efd13afcf
md"""
## Load DICOMS

All you need to do is set `base_path` once and leave it. After that, the only thing that should change is the `VENDER`, once for every set, and the `SCAN_NUMBER`, once for each scan.
"""

# ╔═╡ 80ea14ec-0128-4899-bcae-5d75d6654a89
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ cb6af1f1-9e54-4172-9e65-bcf534840716
begin
	SCAN_NUMBER = 2
	VENDER = "Canon_Aquilion_One_Vision"
	BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"
end

# ╔═╡ 47504ab4-a52f-418d-b9de-3468304c9578
root_path = string(BASE_PATH, VENDER)

# ╔═╡ d41f610c-928b-4a36-a7f8-ac50ff81f6db
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 7efcfa5c-8d99-4e26-9538-efa1d2bfeea6
pth = dcm_path_list[SCAN_NUMBER]

# ╔═╡ 017121d9-a3bd-481a-ba89-7cd261f0ffe7
pth

# ╔═╡ 49fe69bc-7c19-4dc6-86df-ceaccbdd65f2
scan = basename(pth)

# ╔═╡ e95b1cbc-82a1-47db-8026-e2773f0d0cff
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ a4f1ebf0-6f99-42fe-96f3-765caf5f207d
md"""
## Helper Functions
"""

# ╔═╡ 8e874d2a-95bc-414b-91cd-1e6df8220c18
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 1fec01bc-445a-4961-b518-f25a0de737d3
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=3, show_value=true)
end

# ╔═╡ 4494460e-4ffc-4317-9a30-21ac8e005bea
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
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=0.75, color=:red)
	fig
end

# ╔═╡ 3bf98e23-b918-4584-bbb1-d6e125831b16
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 8246c747-7f18-4a67-a861-05f9c0a4049e
function create_circular_mask(h, w, center_circle, radius_circle)
    Y, X = collect(1:h), collect(1:w)'
    dist_from_center = sqrt.((X .- center_circle[1]) .^ 2 .+ (Y .- center_circle[2]) .^ 2)
    mask = dist_from_center .<= radius_circle
    return mask
end

# ╔═╡ 452858f9-6719-46d4-bfbb-c6dbcf402ad6
function dilate_mask_large(mask)
    return dilate(mask)
end

# ╔═╡ 6220010e-ce94-4c14-abe5-8a987501f7b9
function ring_mask_large(dilated_mask)
    return Bool.(((dilate(dilate(dilate(dilate(dilated_mask)))))) - dilate(dilated_mask))
end

# ╔═╡ 0f4ed85d-4141-4f81-b415-6d21535a6584
function dilate_mask_medium(mask)
    return dilate(dilate(dilate(mask)))
end

# ╔═╡ 9e0509b5-a505-4fcd-8095-25f8684f40f3
function ring_mask_medium(dilated_mask)
    return Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilated_mask)))))) - dilated_mask)
end

# ╔═╡ c833d5e2-e046-4472-a072-1bf0417d9bca
function dilate_mask_small(mask)
    return dilate(dilate(dilate(dilate(dilate(mask)))))
end

# ╔═╡ 26f006e4-d84a-4d08-819d-cf35d7596e5f
function ring_mask_small(dilated_mask)
    return Bool.((((dilate(dilate(dilate(dilated_mask)))))) - dilate(dilated_mask))
end

# ╔═╡ 9ce07038-fbf6-4a05-b64d-24f70dfefc9c
function dilate_mask_large_bkg(mask)
    return dilate(dilate(mask))
end

# ╔═╡ aed31300-433e-44c6-a5de-f723bc96ca4f
function dilate_mask_medium_bkg(mask)
    return dilate(mask)
end

# ╔═╡ 29e33799-5513-44dc-a6e2-3c8acc967066
function dilate_mask_small_bkg(mask)
    return (mask)
end

# ╔═╡ 06537b73-9579-4ab9-b829-a4bfa1a8614a
md"""
# Segmentations
"""

# ╔═╡ 9d4357e3-5ee5-4bf9-8879-473a4da59812
md"""
## Heart
"""

# ╔═╡ 4a947165-c565-4c6d-b63c-8cc714d0d039
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ f1339241-169d-40c1-becc-8f230c12e835
@bind a PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ 641043a6-6d00-4f20-96ad-4b588bda2893
heatmap(masked_array[:, :, a], colormap=:grays)

# ╔═╡ b6d3fb6f-fb65-474b-9ee4-29368ed58372
let
    fig = Figure()

    ax = Makie.Axis(fig[1, 1])
	heatmap!(transpose(dcm_array[:, :, 25]), colormap=:grays)

	save("/Users/daleblack/Google Drive/Research/Papers/My Papers/cac-simulation/figures-review/physical_phantom.png", fig)
    fig
end

# ╔═╡ ede1723d-ac8b-43a5-a6d8-6a92fcbaa2fa
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 25]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ 13671de8-9397-422f-81a4-57781e881178
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig2
end

# ╔═╡ c07994fb-87c8-4a61-8dd9-7fe50b087b13
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig3
end

# ╔═╡ 2988d4aa-a2f2-45a5-9e96-fb9c233e9430
md"""
## Calcium Rod
"""

# ╔═╡ 06c46a2f-9f03-46ec-88ed-6219bd2a5f8d
calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header);

# ╔═╡ 5c4cb5db-83f3-4a0f-b2a2-be28dc151b81
slice_CCI

# ╔═╡ c91551f2-38e6-4d2e-b286-7624405339d0
@bind slice PlutoUI.Slider(axes(masked_array, 3); default=slice_CCI, show_value=true)

# ╔═╡ 8233c7a6-138f-4259-9b93-3324640d9820
heatmap(masked_array[:, :, slice], colormap=:grays)

# ╔═╡ 66c50557-49d0-4d9f-beae-404b6784f079
@bind c PlutoUI.Slider(1:size(calcium_image, 3), default=cal_rod_slice, show_value=true)

# ╔═╡ 8806c9d6-c416-47be-981e-a03d35d45602
heatmap(transpose(calcium_image[:, :, c]), colormap=:grays)

# ╔═╡ f30cd88e-c6ca-495d-8ef9-6f029182b5a0
md"""
## Calcium Inserts
"""

# ╔═╡ 3bb994fa-e4e8-4e98-bf1b-5e896f21e96e
begin
	thresh = 115
	mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=0, calcium_threshold=thresh)
end

# ╔═╡ 0aece8c7-c8f7-49ca-b68e-b2fda92060c0
masks = mask_L_HD + mask_M_HD + mask_S_HD + mask_L_MD + mask_M_MD + mask_S_MD + mask_L_LD + mask_M_LD + mask_S_LD;

# ╔═╡ 5eb2cb0e-5353-4bd2-86ac-62f17be96f4c
heatmap(masks, colormap=:grays)

# ╔═╡ a914f5b8-84ab-4b82-886a-7a0fba603252
md"""
## Dilated/Eroded Masks
"""

# ╔═╡ a5189ede-552e-48ea-9660-b5a655ed23f0
begin
	arr = masked_array[:, :, slice_CCI-1:slice_CCI+2]
	single_arr = masked_array[:, :, slice_CCI]
	pixel_size = DICOMUtils.get_pixel_size(header)
end

# ╔═╡ 66e1c059-6ae2-45e2-9deb-1242d8d86999
begin
	# Segmentations
	## Background
	### Large
	mask_bkg_L = create_circular_mask(size(masked_array, 1), size(masked_array, 2), (center_insert[2], center_insert[1]), 10)
	mask_bkg_L_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_bkg_L_3D[:, :, z] = mask_bkg_L
	end
	dilated_mask_L_bkg = dilate_mask_large(mask_bkg_L_3D)
	ring_mask_L_bkg = ring_mask_large(dilated_mask_L_bkg)

	### Medium
	mask_bkg_M = create_circular_mask(size(masked_array, 1), size(masked_array, 2), (center_insert[2], center_insert[1]), 8)
	mask_bkg_M_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_bkg_M_3D[:, :, z] = mask_bkg_M
	end
	dilated_mask_M_bkg = dilate_mask_large(mask_bkg_M_3D)
	ring_mask_M_bkg = ring_mask_large(dilated_mask_M_bkg)

	### Small
	mask_bkg_S = create_circular_mask(size(masked_array, 1), size(masked_array, 2), (center_insert[2], center_insert[1]), 6)
	mask_bkg_S_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_bkg_S_3D[:, :, z] = mask_bkg_S
	end
	dilated_mask_S_bkg = dilate_mask_large(mask_bkg_S_3D)
	ring_mask_S_bkg = ring_mask_large(dilated_mask_S_bkg)
	
	## Large
	### High Density
	mask_L_HD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
	end
	dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D)
	ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)

	### Medium Density
	mask_L_MD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_L_MD_3D[:, :, z] = mask_L_MD
	end
	dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D)
	ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)

	### Low Density
	mask_L_LD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_L_LD_3D[:, :, z] = mask_L_LD
	end
	dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D)
	ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)


	## Medium 
	### High Density
	mask_M_HD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_M_HD_3D[:, :, z] = mask_M_HD
	end
	dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D)
	ring_mask_M_HD = ring_mask_medium(dilated_mask_M_HD)

	### Medium Density
	mask_M_MD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_M_MD_3D[:, :, z] = mask_M_MD
	end
	dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
	ring_mask_M_MD = ring_mask_medium(dilated_mask_M_MD)

	### Low Density
	mask_M_LD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_M_LD_3D[:, :, z] = mask_M_LD
	end
	dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D)
	ring_mask_M_LD = ring_mask_medium(dilated_mask_M_LD)

	## Small
	### High Density
	mask_S_HD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_S_HD_3D[:, :, z] = mask_S_HD
	end
	dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D)
	ring_mask_S_HD = ring_mask_small(dilated_mask_S_HD)

	### Medium Density
	mask_S_MD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_S_MD_3D[:, :, z] = mask_S_MD
	end
	dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D)
	ring_mask_S_MD = ring_mask_small(dilated_mask_S_MD)

	### Low Density
	mask_S_LD_3D = Array{Bool}(undef, size(arr))
	for z in axes(arr, 3)
		mask_S_LD_3D[:, :, z] = mask_S_LD
	end
	dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D)
	ring_mask_S_LD = ring_mask_small(dilated_mask_S_LD)
end;

# ╔═╡ a488da37-5fd7-4d23-a958-d2fe3f1249b1
md"""
# Prepare Calibration
"""

# ╔═╡ 2c5b2279-5c5b-4c97-9feb-3e1346da978c
function erode_array(array, num_erosions, i=1)
	while i <= num_erosions
		@info i
		return erode_array(erode(array), num_erosions, i + 1) 
	end
end

# ╔═╡ e5681b3b-acbf-4e65-9679-e50dcddf83ab
function erode_array2(array, num_erosions, i)
	for i in 1:num_erosions
		@info i
		erode_array2(erode(array), num_erosions, i+1)
	end
end

# ╔═╡ c7ca2d20-33d5-4b45-b936-6198e33a27f4
function erode_array_test(array::AbstractArray, boolean_array::Bool, delta_std)
	eroded_array = erode(boolean_array)
	std1, std2 = std(array[boolean_array], array[eroded_array])
	if std2 - std1 > delta_std
		erode_array2(array, eroded_array, delta_std)
	else
		return eroded_array
	end
end

# ╔═╡ d003e211-3fbb-4c53-8ef4-4c60b40c9f82


# ╔═╡ 390b256d-015e-4581-b72d-c0920291c576
begin
	cal_path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/calibrations.csv"
	
	cal_df = CSV.read(cal_path, DataFrame)
	
	c_img = calcium_image[:, :, cal_rod_slice]
	array_filtered = abs.(mapwindow(median, c_img, (3, 3)))
	bool_arr = array_filtered .> 0
	bool_arr_erode = (((((repeat(bool_arr))))))
	# bool_arr_erode = erode_array(bool_arr, 4)

	hu_calcium = mean(c_img[bool_arr_erode])
	hu_heart_tissue = mean(single_arr[dilated_mask_L_bkg[:, :, 1]])
	# arr_cal = Array(cal_df[3, 2:end])
	# intensity_array3 = vcat(hu_heart_tissue, arr_cal)
	# density_array_calc3 = [
	# 		0
	# 		25
	# 		50
	# 		100
	# 		200
	# 		400
	# 		800
	# 	]

	# intensity_array3 = [40, 296, 1100]
	# density_array_calc3 = [0, 200, 800]
	intensity_array3 = [hu_heart_tissue, hu_calcium]
	density_array_calc3 = [0, 200]
	df_cal = DataFrame(:density => density_array_calc3, :intensity => intensity_array3)
end

# ╔═╡ a50ca0be-2613-498a-bcb8-5bdf10b4eeed
heatmap(bool_arr, colormap=:grays)

# ╔═╡ a2ea95bc-8840-4fbe-a29d-3618c629e9dc
heatmap(bool_arr_erode, colormap=:grays)

# ╔═╡ 07f287df-b706-47b4-9fb9-1df0d6837c4e
cal_df

# ╔═╡ 283eeacd-80e8-43ef-acdc-a0a72e5b8af5
erode_array(bool_arr, 4)

# ╔═╡ 8c4039c8-85a4-40f2-bbc2-eb4da8f84925
erode_array2(bool_arr, 4, 1)

# ╔═╡ 1b3d70c8-2498-402a-a5fd-69d199254280
mean(c_img[bool_arr_erode]), std(c_img[bool_arr_erode])

# ╔═╡ e97d0ae4-58a6-4f27-8a75-585fe3602e79
mean(c_img[bool_arr_erode_new]), std(c_img[bool_arr_erode_new])

# ╔═╡ 06cbd5ff-c18d-4554-8450-f4f88f3d9cd2
hu_calcium

# ╔═╡ 4c5b5a15-a456-4f67-977d-55a727f6ac7f
begin

	linearRegressor = lm(@formula(intensity ~ density), df_cal)
	linearFit = predict(linearRegressor)
	m = linearRegressor.model.pp.beta0[2]
	b = linearRegressor.model.rr.mu[1]
	density(intensity) = (intensity - b) / m
	intensity(ρ) = m*ρ + b
end;

# ╔═╡ 9bdf9ab6-ad78-405a-8a60-71d16a0d6aec
begin
	s_obj = intensity(200)
	dens = 0.2
end

# ╔═╡ af1baf67-bc44-4ed4-87ab-6775a8df1cfc
md"""
# Score Background
"""

# ╔═╡ f5d6851f-0892-4c96-840d-98002bb4baff
md"""
## Large
"""

# ╔═╡ 1a0121f0-5b52-4293-813d-11fcf2400438
@bind bg_l overlay_mask_bind(dilated_mask_L_bkg)

# ╔═╡ 3804dfde-57bc-40e0-994c-fccd37cbe549
overlay_mask_plot(arr, dilated_mask_L_bkg, bg_l, "dilated mask")

# ╔═╡ a019363e-4303-47f2-8e98-3ca11e05013c
overlay_mask_plot(arr, ring_mask_L_bkg, bg_l, "dilated mask")

# ╔═╡ 5c93816d-2d73-46a1-95fe-d2327328ba44
begin
	S_Obj_bkg = s_obj
	ρ_bkg = dens

	S_bkg_large = mean(arr[ring_mask_L_bkg])
	alg_bkg_large = Integrated(arr[dilated_mask_L_bkg])
	mass_large_bkg = score(S_bkg_large, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg_large)
end

# ╔═╡ 5a9bcd46-d89d-4e08-81a1-3de7c17791d8
md"""
## Medium
"""

# ╔═╡ 56a95f11-0f9f-465a-a63c-43dff6d81b79
@bind bg_m overlay_mask_bind(dilated_mask_M_bkg)

# ╔═╡ ebc56920-edaf-4bbc-8a99-8fa9a9a44052
overlay_mask_plot(arr, dilated_mask_M_bkg, bg_m, "dilated mask")

# ╔═╡ 375aae93-204c-49cf-9d3c-8391715c1b83
overlay_mask_plot(arr, ring_mask_M_bkg, bg_m, "dilated mask")

# ╔═╡ 70410004-1507-432f-bc1b-adf7c0499bc7
begin
	S_bkg_medium = mean(arr[ring_mask_M_bkg])
	alg_bkg_medium = Integrated(arr[dilated_mask_M_bkg])
	mass_medium_bkg = score(S_bkg_medium, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg_medium)
end

# ╔═╡ c11cc26b-0e79-4634-8c5f-d1583eb192d6
md"""
## Small
"""

# ╔═╡ 032a4c83-4390-4b23-b4aa-c1d06f1f0fcd
@bind bg_s overlay_mask_bind(dilated_mask_S_bkg)

# ╔═╡ 4b9a1cdd-54ab-4ab7-9e19-ccaaf980cbdf
overlay_mask_plot(arr, dilated_mask_S_bkg, bg_s, "dilated mask")

# ╔═╡ bbe51fc4-5988-4022-a110-8c894511de5f
overlay_mask_plot(arr, ring_mask_S_bkg, bg_s, "dilated mask")

# ╔═╡ 3f73a09a-7b27-484a-940a-da2247e5fda9
begin
	# S_Obj_bkg = intensity(200)
	# ρ_bkg = 0.2 # mg/mm^3

	S_bkg_small = mean(arr[ring_mask_S_bkg])
	alg_bkg_small = Integrated(arr[dilated_mask_S_bkg])
	mass_small_bkg = score(S_bkg_small, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg_small)	
end

# ╔═╡ 2903b7ba-f028-40e4-81ce-3da32c6553b6
mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]

# ╔═╡ db84d9bb-e0da-4f81-88c8-c2a20c3ce006
md"""
# Score Large Inserts
"""

# ╔═╡ 2d83f341-24e4-4f65-885b-e5852e1c5d96
md"""
## High Density
"""

# ╔═╡ 4d08c2b9-b373-4fc4-ad4f-af7621f8e55f
md"""
#### Dilated mask
"""

# ╔═╡ 4da96d00-fc05-4b05-b41e-ff6527c2de36
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 258d3494-a618-4db6-84a2-189b24bb15db
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 337109bb-42b2-4f6f-83f6-e0e98e8838c3
md"""
#### Ring (background) mask
"""

# ╔═╡ 5faa2499-1017-4182-92bc-cb6f18401617
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ 6af0d533-6eb4-4129-8b12-a5eefe82b2c0
overlay_mask_plot(arr, ring_mask_L_HD, g4, "ring mask")

# ╔═╡ 45557e45-f916-43a9-a2c5-944b154ee264
begin
	single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
	s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
end

# ╔═╡ 1a7c8177-6f86-465e-b0dd-2756ff0db034
begin
	# S_Obj_HD = intensity(200)
	# ρ_hd = 0.2 # mg/mm^3
	S_Obj_HD = s_obj
	ρ_hd = dens
end

# ╔═╡ 8305356e-9b8e-4ccf-83bd-2258b4bfa2d6
begin
	alg_L_HD = Integrated(arr[mask_L_HD_3D])
	mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
end

# ╔═╡ 29d535c8-994d-42e1-ad34-48453a3f92f4
md"""
## Medium Density
"""

# ╔═╡ ec4c1348-6eac-4372-88b9-7323d20a80f9
md"""
#### Dilated mask
"""

# ╔═╡ 1e5bb04d-ecd8-4036-85e9-05ba522e65fd
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 5b01eec5-b967-4e52-ada2-d7e6f98e8743
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ bd50ca64-f6eb-4827-b8eb-a4186648d3c7
md"""
#### Ring (background) mask
"""

# ╔═╡ 370c4549-69f5-4835-ae55-fc19e61c6c01
@bind h4 overlay_mask_bind(ring_mask_L_MD)

# ╔═╡ 835b7446-c0dc-4c31-a172-d52af2a81b89
overlay_mask_plot(arr, ring_mask_L_MD, h4, "ring mask")

# ╔═╡ d189b76f-d6ca-455a-9432-ba054377557a
begin
	single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
	s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
end

# ╔═╡ 16e21e05-b582-4986-b2e6-f7e40ea77f39
begin
	# S_Obj_MD = intensity(200)
	# ρ_md = 0.2 # mg/mm^3
	S_Obj_MD = s_obj
	ρ_md = dens
end

# ╔═╡ b7a1a6d7-82a6-40aa-afcb-5323a2642c56
begin
	alg_L_MD = Integrated(arr[mask_L_MD_3D])
	mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
end

# ╔═╡ a42e982c-16e5-45b3-a9a5-859375b72f05
md"""
## Low Density
"""

# ╔═╡ 6d2174b3-16ce-4fe2-9ec1-9c945fa6629e
md"""
#### Dilated mask
"""

# ╔═╡ ba2c6db8-b246-45ba-b667-e7b1883894db
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ a05fb79c-3782-460a-a4ee-54b943198e2e
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ fc0f8bd0-40dc-410b-a980-252179ecd3fc
md"""
#### Ring (background) mask
"""

# ╔═╡ e9aa29a4-4480-44da-aeb8-1de14cf382b2
@bind i4 overlay_mask_bind(ring_mask_L_LD)

# ╔═╡ c3638503-8a36-4416-8842-88d4c5bda21c
overlay_mask_plot(arr, ring_mask_L_LD, i4, "ring mask")

# ╔═╡ d1e942c5-8ef5-443d-ab5e-083e0d81973e
begin	
	single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
	s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
end

# ╔═╡ 54ba65d1-e7e8-40cb-b016-43686789e0e6
begin
	# S_Obj_LD = intensity(200)
	# ρ_ld = 0.2 # mg/mm^3
	S_Obj_LD = s_obj
	ρ_ld = dens
end

# ╔═╡ dc424553-ab65-4013-aa1e-0d8163ac0c6f
begin
	alg_L_LD = Integrated(arr[mask_L_LD_3D])
	mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_ld, alg_L_LD)
end

# ╔═╡ 24b943da-9319-4360-9d41-992f41752bb8
md"""
# Score Medium Inserts
"""

# ╔═╡ ea1e7093-51c3-4058-ad9b-9a43e67b019f
md"""
## High Density
"""

# ╔═╡ fa5e0bef-3991-4b42-921d-60015da8000a
md"""
#### Dilated mask
"""

# ╔═╡ 21f76e13-8748-4eaf-a38e-4fd8a0e9d861
# dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))));

# ╔═╡ 344026d1-d9b5-4dd0-bce5-5b668e9545cc
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ f6e68cba-833a-467e-af86-658dcc5bf6c7
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 33bf73da-96c1-4ea3-a6a2-c022c7b64bb9
md"""
#### Ring (background) mask
"""

# ╔═╡ 894f3319-b750-490e-8e45-bc3710efb0ce
@bind j4 overlay_mask_bind(ring_mask_M_HD)

# ╔═╡ 82acb4b7-85ae-42c7-bf42-fe2ae1abdaf1
overlay_mask_plot(arr, ring_mask_M_HD, j4, "ring mask")

# ╔═╡ 3163aaaa-2463-4664-aa6c-6e44988aa337
begin
	single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
	s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
end

# ╔═╡ 8af7d835-c43b-40b7-a934-4f6ce7661535
begin
	alg_M_HD = Integrated(arr[mask_M_HD_3D])
	mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
end

# ╔═╡ 18eae514-f9f8-4cd3-9626-0d702a8fda80
md"""
## Medium Density
"""

# ╔═╡ 43432b7f-cb36-4fc7-9b2a-97920880004a
md"""
#### Dilated mask
"""

# ╔═╡ 5a320923-1966-4a68-b90d-af1da17e79ea
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 439dd6bb-ae38-4cd0-ad89-002c9d91a7c1
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 63abca9e-2111-4799-be13-ab3f86d5b082
md"""
#### Ring (background) mask
"""

# ╔═╡ 495cc420-a24f-41d2-a0cb-116f778c0bff
@bind k4 overlay_mask_bind(ring_mask_M_MD)

# ╔═╡ 2e52ff34-db03-404c-ad98-c8f62df0d181
overlay_mask_plot(arr, ring_mask_M_MD, k4, "ring mask")

# ╔═╡ 0687ccfd-d58a-4fa4-aa00-815101997a75
begin
	single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
	s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
end

# ╔═╡ 980785a4-0c19-4208-89eb-2b2b4e85523d
begin
	alg_M_MD = Integrated(arr[mask_M_MD_3D])
	mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
end

# ╔═╡ 48400ead-acd0-4c5b-8ff7-87ee2b07b9d1
md"""
## Low Density
"""

# ╔═╡ c8a9335b-c271-4097-be57-455119cd461e
md"""
#### Dilated mask
"""

# ╔═╡ 223efc65-64a0-4d27-b829-3288445967d4
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ e2a2251d-9964-4fa3-becc-09d0b12cd7f0
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ f7c88984-fcaf-4173-8f8b-023979aa472f
md"""
#### Ring (background) mask
"""

# ╔═╡ b03ba157-8b69-4616-8723-f0396bf184c6
@bind l4 overlay_mask_bind(ring_mask_M_LD)

# ╔═╡ f58e9581-8f7d-477e-8820-9c13e210b385
overlay_mask_plot(arr, ring_mask_M_LD, l4, "ring mask")

# ╔═╡ ceecd152-e1f6-497a-872b-c4de20345499
begin
	single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
	s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
end

# ╔═╡ 1a280788-84f9-4a3d-9084-b1ab9c2ba514
begin
	alg_M_LD = Integrated(arr[mask_M_LD_3D])
	mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_ld, alg_M_LD)
end

# ╔═╡ 16d56ff5-238e-4916-bac6-56fe70164b54
md"""
# Score Small Inserts
"""

# ╔═╡ 37acb2c8-30e2-441d-90d2-884ea0fd2d49
md"""
## High Density
"""

# ╔═╡ fdbf5185-6a1f-440b-aa81-81529e045b73
md"""
#### Dilated mask
"""

# ╔═╡ b1a6f986-4d86-4c62-a662-552bdf9421f9
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ de05e5fd-3972-41e9-bdaa-eba6631d8ecb
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 339617d5-49c4-4bad-8c03-579438b55799
md"""
#### Ring (background) mask
"""

# ╔═╡ b712ebf8-bf99-40c8-b861-19e6dfd968c0
@bind m4 overlay_mask_bind(ring_mask_S_HD)

# ╔═╡ d7bbcda2-37bf-43c2-9a4e-e6027d868d26
overlay_mask_plot(arr, ring_mask_S_HD, m4, "ring mask")

# ╔═╡ 10550a99-76b9-487e-922b-3239e7f811f4
begin
	single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
	s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
end

# ╔═╡ 933791e4-934d-40fd-8928-eccb68273059
begin
	alg_S_HD = Integrated(arr[mask_S_HD_3D])
	mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
end

# ╔═╡ cdc09f53-4f83-46ff-8124-bdf54cea5097
md"""
## Medium Density
"""

# ╔═╡ d27098ee-43ae-4072-9843-11934f377df6
md"""
#### Dilated mask
"""

# ╔═╡ 6467393e-8590-45db-924b-0033e9ddb724
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 522d78ed-ee09-4c72-a8d3-4678d4677e1b
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "ring mask")

# ╔═╡ f1057c02-9060-4d44-a0d9-6b04245c9ed0
md"""
#### Ring (background) mask
"""

# ╔═╡ 9b254664-72a6-4030-83b0-27db11c222ba
@bind n4 overlay_mask_bind(ring_mask_S_MD)

# ╔═╡ 60370f6d-51f3-436f-ab32-400de8d98f2b
overlay_mask_plot(arr, ring_mask_S_MD, n4, "ring mask")

# ╔═╡ d62d50a8-dd17-4d25-b353-fdfedb456ff0
begin
	single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
	s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
end

# ╔═╡ a77205bb-eab2-4a2a-90b7-2f4844f039de
begin
	alg_S_MD = Integrated(arr[mask_S_MD_3D])
	mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
end

# ╔═╡ b3455451-651e-4f70-87f4-08f9db13d9b5
md"""
## Low Density
"""

# ╔═╡ 35b79339-bb2e-4987-a0f4-583a732c71fd
md"""
#### Dilated mask
"""

# ╔═╡ e35a7298-79fa-42b3-a2c7-faec6b5c9af4
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 86fac660-6ef4-4041-80d8-61eb6edc2a39
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ 005f8151-6b2b-4476-b92d-e63d0d216325
md"""
#### Ring (background) mask
"""

# ╔═╡ b71e334a-8452-453c-a7a4-6d5f6eafd0bc
@bind o4 overlay_mask_bind(ring_mask_S_LD)

# ╔═╡ 13806613-6eca-4b9b-8216-894b25ff1b01
overlay_mask_plot(arr, ring_mask_S_LD, o4, "ring mask")

# ╔═╡ 9a8f53e1-4900-4b64-a09f-5a9b8fd5a0c5
begin
	single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
	s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
end

# ╔═╡ 05233816-9613-4cb6-90f0-bac0f01816c0
begin
	alg_S_LD = Integrated(arr[mask_S_LD_3D])
	mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_ld, alg_S_LD)
end

# ╔═╡ b7f35164-9e08-49ea-aa40-f968fbbb92c2
md"""
# Results
"""

# ╔═╡ b5e14014-12bf-4004-bccc-2f3c759200b5
density_array = [0, 200, 400, 800]

# ╔═╡ 88edf2c7-3c5a-4b90-9c50-5d448a358ce8
inserts = [
	"Low Density",
	"Medium Density",
	"High Density"
]

# ╔═╡ 786eb037-e39e-41e1-9a63-6106800cb306
ground_truth_mass_large = [
	19.6,
	39.3,
	78.5
] # mg

# ╔═╡ fe51b3b3-31bc-454f-bc09-db5a7e6c6941
calculated_mass_large = [
	mass_l_ld,
	mass_l_md,
	mass_l_hd
]

# ╔═╡ 2867365f-f7fe-460c-b8fd-eb8640f7dc73
ground_truth_mass_medium = [
	4.2,
	8.5,
	17.0
]

# ╔═╡ 61cffd44-a5e8-4223-85ef-eb4f2b22de2b
calculated_mass_medium = [
	mass_m_ld,
	mass_m_md,
	mass_m_hd
]

# ╔═╡ 60324b4b-e880-4ea3-9831-fca775eb22ab
ground_truth_mass_small = [
	0.2,
	0.3,
	0.6
]

# ╔═╡ 7c4fe694-7343-4e87-aed4-8895802e34a9
calculated_mass_small = [
	mass_s_ld,
	mass_s_md,
	mass_s_hd
]

# ╔═╡ 9f12b7bd-ccc2-4108-9171-1902a0c75c30
df = DataFrame(
	scan = scan,
	inserts = inserts,
	mass_bkg = mass_bkg,
	ground_truth_mass_large = ground_truth_mass_large,
	calculated_mass_large = calculated_mass_large,
	ground_truth_mass_medium = ground_truth_mass_medium,
	calculated_mass_medium = calculated_mass_medium,
	ground_truth_mass_small = ground_truth_mass_small,
	calculated_mass_small = calculated_mass_small
)

# ╔═╡ 3c68054c-604b-4060-bf00-37cb36fa30d4
df

# ╔═╡ 2d55c0c1-9292-40bc-8189-f73783b28d47
mn, st = mean(df[!, :mass_bkg]), std(df[!, :mass_bkg])

# ╔═╡ 129baab5-4aa3-4b8b-9f5f-1842e558cd67
df[!, :calculated_mass_small] .< mn + st

# ╔═╡ 958a9d50-c61f-4d8e-89e7-7b37f0fd4246
df[!, :calculated_mass_small] .< mn + st

# ╔═╡ fef2fcda-cc64-469e-9f17-946f06a3e4ae
begin
	fmass2 = Figure()
	axmass2 = Axis(fmass2[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_large], label="ground_truth_mass_large")
	scatter!(density_array[2:end], df[!, :calculated_mass_large], label="calculated_mass_large")
	
	axmass2.title = "Mass Measurements (Large)"
	axmass2.ylabel = "Mass (mg)"
	axmass2.xlabel = "Density (mg/cm^3)"

	xlims!(axmass2, 0, 850)
	ylims!(axmass2, 0, 100)
	
	fmass2[1, 2] = Legend(fmass2, axmass2, framevisible = false)
	
	fmass2
end

# ╔═╡ 4f9c5043-3485-45b8-b493-d62780134251
begin
	fmass3 = Figure()
	axmass3 = Axis(fmass3[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_medium], label="ground_truth_mass_medium")
	scatter!(density_array[2:end], df[!, :calculated_mass_medium], label="calculated_mass_medium")
	
	axmass3.title = "Mass Measurements (Medium)"
	axmass3.ylabel = "Mass (mg)"
	axmass3.xlabel = "Density (mg/cm^3)"

	xlims!(axmass3, 0, 850)
	ylims!(axmass3, 0, 25)
	
	fmass3[1, 2] = Legend(fmass3, axmass3, framevisible = false)
	
	fmass3
end

# ╔═╡ 197e0479-f15a-43f9-a7cf-5b5d58b5a901
begin
	fmass4 = Figure()
	axmass4 = Axis(fmass4[1, 1])
	
	scatter!(density_array[2:end], df[!, :ground_truth_mass_small], label="ground_truth_mass_small")
	scatter!(density_array[2:end], df[!, :calculated_mass_small], label="calculated_mass_small")
	
	axmass4.title = "Mass Measurements (Small)"
	axmass4.ylabel = "Mass (mg)"
	axmass4.xlabel = "Density (mg/cm^3)"

	xlims!(axmass4, 0, 850)
	ylims!(axmass4, 0, 2)
	
	fmass4[1, 2] = Legend(fmass4, axmass4, framevisible = false)
	
	fmass4
end

# ╔═╡ b030dcf1-fc55-4e32-be9b-6cba1d3456f0
md"""
### Save Results
"""

# ╔═╡ 5a9c7e99-0e6e-4162-8ff7-bff1935fd318
# if ~isdir(string(cd(pwd, "..") , "/data/output/", VENDER))
# 	mkdir(string(cd(pwd, "..") , "/data/output/", VENDER))
# end

# ╔═╡ 7976fbad-ae1f-40e8-961f-2ea68f752caa
# output_path = string(cd(pwd, "..") , "/data/output/", VENDER, "/", scan, ".csv")

# ╔═╡ d8aeb68a-697d-4477-8b6c-14712f50a945
# CSV.write(output_path, df)

# ╔═╡ Cell order:
# ╠═3c162d31-407a-4c50-a8a7-30153192f207
# ╠═0aeb13cf-6299-473a-908e-5b0be1b23c86
# ╟─c4bba726-9f22-4f76-a4ae-966efd13afcf
# ╟─80ea14ec-0128-4899-bcae-5d75d6654a89
# ╠═cb6af1f1-9e54-4172-9e65-bcf534840716
# ╠═47504ab4-a52f-418d-b9de-3468304c9578
# ╠═d41f610c-928b-4a36-a7f8-ac50ff81f6db
# ╠═7efcfa5c-8d99-4e26-9538-efa1d2bfeea6
# ╠═017121d9-a3bd-481a-ba89-7cd261f0ffe7
# ╠═49fe69bc-7c19-4dc6-86df-ceaccbdd65f2
# ╠═e95b1cbc-82a1-47db-8026-e2773f0d0cff
# ╟─a4f1ebf0-6f99-42fe-96f3-765caf5f207d
# ╟─8e874d2a-95bc-414b-91cd-1e6df8220c18
# ╟─1fec01bc-445a-4961-b518-f25a0de737d3
# ╟─4494460e-4ffc-4317-9a30-21ac8e005bea
# ╟─3bf98e23-b918-4584-bbb1-d6e125831b16
# ╟─8246c747-7f18-4a67-a861-05f9c0a4049e
# ╟─452858f9-6719-46d4-bfbb-c6dbcf402ad6
# ╟─6220010e-ce94-4c14-abe5-8a987501f7b9
# ╟─0f4ed85d-4141-4f81-b415-6d21535a6584
# ╟─9e0509b5-a505-4fcd-8095-25f8684f40f3
# ╟─c833d5e2-e046-4472-a072-1bf0417d9bca
# ╟─26f006e4-d84a-4d08-819d-cf35d7596e5f
# ╟─9ce07038-fbf6-4a05-b64d-24f70dfefc9c
# ╟─aed31300-433e-44c6-a5de-f723bc96ca4f
# ╟─29e33799-5513-44dc-a6e2-3c8acc967066
# ╟─06537b73-9579-4ab9-b829-a4bfa1a8614a
# ╟─9d4357e3-5ee5-4bf9-8879-473a4da59812
# ╠═4a947165-c565-4c6d-b63c-8cc714d0d039
# ╟─f1339241-169d-40c1-becc-8f230c12e835
# ╟─641043a6-6d00-4f20-96ad-4b588bda2893
# ╟─b6d3fb6f-fb65-474b-9ee4-29368ed58372
# ╟─ede1723d-ac8b-43a5-a6d8-6a92fcbaa2fa
# ╟─13671de8-9397-422f-81a4-57781e881178
# ╟─c07994fb-87c8-4a61-8dd9-7fe50b087b13
# ╟─2988d4aa-a2f2-45a5-9e96-fb9c233e9430
# ╠═06c46a2f-9f03-46ec-88ed-6219bd2a5f8d
# ╠═5c4cb5db-83f3-4a0f-b2a2-be28dc151b81
# ╟─c91551f2-38e6-4d2e-b286-7624405339d0
# ╟─8233c7a6-138f-4259-9b93-3324640d9820
# ╟─66c50557-49d0-4d9f-beae-404b6784f079
# ╠═8806c9d6-c416-47be-981e-a03d35d45602
# ╟─f30cd88e-c6ca-495d-8ef9-6f029182b5a0
# ╠═3bb994fa-e4e8-4e98-bf1b-5e896f21e96e
# ╠═0aece8c7-c8f7-49ca-b68e-b2fda92060c0
# ╟─5eb2cb0e-5353-4bd2-86ac-62f17be96f4c
# ╟─a914f5b8-84ab-4b82-886a-7a0fba603252
# ╠═a5189ede-552e-48ea-9660-b5a655ed23f0
# ╠═66e1c059-6ae2-45e2-9deb-1242d8d86999
# ╟─a488da37-5fd7-4d23-a958-d2fe3f1249b1
# ╠═a50ca0be-2613-498a-bcb8-5bdf10b4eeed
# ╠═a2ea95bc-8840-4fbe-a29d-3618c629e9dc
# ╠═07f287df-b706-47b4-9fb9-1df0d6837c4e
# ╠═2c5b2279-5c5b-4c97-9feb-3e1346da978c
# ╠═e5681b3b-acbf-4e65-9679-e50dcddf83ab
# ╠═c7ca2d20-33d5-4b45-b936-6198e33a27f4
# ╠═283eeacd-80e8-43ef-acdc-a0a72e5b8af5
# ╠═8c4039c8-85a4-40f2-bbc2-eb4da8f84925
# ╠═d003e211-3fbb-4c53-8ef4-4c60b40c9f82
# ╠═390b256d-015e-4581-b72d-c0920291c576
# ╠═1b3d70c8-2498-402a-a5fd-69d199254280
# ╠═e97d0ae4-58a6-4f27-8a75-585fe3602e79
# ╠═06cbd5ff-c18d-4554-8450-f4f88f3d9cd2
# ╠═3c68054c-604b-4060-bf00-37cb36fa30d4
# ╠═129baab5-4aa3-4b8b-9f5f-1842e558cd67
# ╠═9bdf9ab6-ad78-405a-8a60-71d16a0d6aec
# ╠═4c5b5a15-a456-4f67-977d-55a727f6ac7f
# ╟─af1baf67-bc44-4ed4-87ab-6775a8df1cfc
# ╟─f5d6851f-0892-4c96-840d-98002bb4baff
# ╟─3804dfde-57bc-40e0-994c-fccd37cbe549
# ╟─1a0121f0-5b52-4293-813d-11fcf2400438
# ╟─a019363e-4303-47f2-8e98-3ca11e05013c
# ╠═5c93816d-2d73-46a1-95fe-d2327328ba44
# ╟─5a9bcd46-d89d-4e08-81a1-3de7c17791d8
# ╟─ebc56920-edaf-4bbc-8a99-8fa9a9a44052
# ╟─56a95f11-0f9f-465a-a63c-43dff6d81b79
# ╟─375aae93-204c-49cf-9d3c-8391715c1b83
# ╠═70410004-1507-432f-bc1b-adf7c0499bc7
# ╟─c11cc26b-0e79-4634-8c5f-d1583eb192d6
# ╟─4b9a1cdd-54ab-4ab7-9e19-ccaaf980cbdf
# ╟─032a4c83-4390-4b23-b4aa-c1d06f1f0fcd
# ╟─bbe51fc4-5988-4022-a110-8c894511de5f
# ╠═3f73a09a-7b27-484a-940a-da2247e5fda9
# ╠═2903b7ba-f028-40e4-81ce-3da32c6553b6
# ╟─db84d9bb-e0da-4f81-88c8-c2a20c3ce006
# ╟─2d83f341-24e4-4f65-885b-e5852e1c5d96
# ╟─4d08c2b9-b373-4fc4-ad4f-af7621f8e55f
# ╟─4da96d00-fc05-4b05-b41e-ff6527c2de36
# ╠═258d3494-a618-4db6-84a2-189b24bb15db
# ╟─337109bb-42b2-4f6f-83f6-e0e98e8838c3
# ╟─5faa2499-1017-4182-92bc-cb6f18401617
# ╠═6af0d533-6eb4-4129-8b12-a5eefe82b2c0
# ╠═45557e45-f916-43a9-a2c5-944b154ee264
# ╠═1a7c8177-6f86-465e-b0dd-2756ff0db034
# ╠═8305356e-9b8e-4ccf-83bd-2258b4bfa2d6
# ╟─29d535c8-994d-42e1-ad34-48453a3f92f4
# ╟─ec4c1348-6eac-4372-88b9-7323d20a80f9
# ╟─1e5bb04d-ecd8-4036-85e9-05ba522e65fd
# ╠═5b01eec5-b967-4e52-ada2-d7e6f98e8743
# ╟─bd50ca64-f6eb-4827-b8eb-a4186648d3c7
# ╟─370c4549-69f5-4835-ae55-fc19e61c6c01
# ╠═835b7446-c0dc-4c31-a172-d52af2a81b89
# ╠═d189b76f-d6ca-455a-9432-ba054377557a
# ╠═16e21e05-b582-4986-b2e6-f7e40ea77f39
# ╠═b7a1a6d7-82a6-40aa-afcb-5323a2642c56
# ╟─a42e982c-16e5-45b3-a9a5-859375b72f05
# ╟─6d2174b3-16ce-4fe2-9ec1-9c945fa6629e
# ╟─ba2c6db8-b246-45ba-b667-e7b1883894db
# ╠═a05fb79c-3782-460a-a4ee-54b943198e2e
# ╟─fc0f8bd0-40dc-410b-a980-252179ecd3fc
# ╟─e9aa29a4-4480-44da-aeb8-1de14cf382b2
# ╠═c3638503-8a36-4416-8842-88d4c5bda21c
# ╠═d1e942c5-8ef5-443d-ab5e-083e0d81973e
# ╠═54ba65d1-e7e8-40cb-b016-43686789e0e6
# ╠═dc424553-ab65-4013-aa1e-0d8163ac0c6f
# ╟─24b943da-9319-4360-9d41-992f41752bb8
# ╟─ea1e7093-51c3-4058-ad9b-9a43e67b019f
# ╟─fa5e0bef-3991-4b42-921d-60015da8000a
# ╠═21f76e13-8748-4eaf-a38e-4fd8a0e9d861
# ╟─344026d1-d9b5-4dd0-bce5-5b668e9545cc
# ╠═f6e68cba-833a-467e-af86-658dcc5bf6c7
# ╟─33bf73da-96c1-4ea3-a6a2-c022c7b64bb9
# ╟─894f3319-b750-490e-8e45-bc3710efb0ce
# ╠═82acb4b7-85ae-42c7-bf42-fe2ae1abdaf1
# ╠═3163aaaa-2463-4664-aa6c-6e44988aa337
# ╠═8af7d835-c43b-40b7-a934-4f6ce7661535
# ╟─18eae514-f9f8-4cd3-9626-0d702a8fda80
# ╟─43432b7f-cb36-4fc7-9b2a-97920880004a
# ╟─5a320923-1966-4a68-b90d-af1da17e79ea
# ╠═439dd6bb-ae38-4cd0-ad89-002c9d91a7c1
# ╟─63abca9e-2111-4799-be13-ab3f86d5b082
# ╟─495cc420-a24f-41d2-a0cb-116f778c0bff
# ╠═2e52ff34-db03-404c-ad98-c8f62df0d181
# ╠═0687ccfd-d58a-4fa4-aa00-815101997a75
# ╠═980785a4-0c19-4208-89eb-2b2b4e85523d
# ╟─48400ead-acd0-4c5b-8ff7-87ee2b07b9d1
# ╟─c8a9335b-c271-4097-be57-455119cd461e
# ╟─223efc65-64a0-4d27-b829-3288445967d4
# ╠═e2a2251d-9964-4fa3-becc-09d0b12cd7f0
# ╟─f7c88984-fcaf-4173-8f8b-023979aa472f
# ╟─b03ba157-8b69-4616-8723-f0396bf184c6
# ╠═f58e9581-8f7d-477e-8820-9c13e210b385
# ╠═ceecd152-e1f6-497a-872b-c4de20345499
# ╠═1a280788-84f9-4a3d-9084-b1ab9c2ba514
# ╟─16d56ff5-238e-4916-bac6-56fe70164b54
# ╟─37acb2c8-30e2-441d-90d2-884ea0fd2d49
# ╟─fdbf5185-6a1f-440b-aa81-81529e045b73
# ╟─b1a6f986-4d86-4c62-a662-552bdf9421f9
# ╠═de05e5fd-3972-41e9-bdaa-eba6631d8ecb
# ╟─339617d5-49c4-4bad-8c03-579438b55799
# ╟─b712ebf8-bf99-40c8-b861-19e6dfd968c0
# ╠═d7bbcda2-37bf-43c2-9a4e-e6027d868d26
# ╠═10550a99-76b9-487e-922b-3239e7f811f4
# ╠═933791e4-934d-40fd-8928-eccb68273059
# ╟─cdc09f53-4f83-46ff-8124-bdf54cea5097
# ╟─d27098ee-43ae-4072-9843-11934f377df6
# ╟─6467393e-8590-45db-924b-0033e9ddb724
# ╟─522d78ed-ee09-4c72-a8d3-4678d4677e1b
# ╟─f1057c02-9060-4d44-a0d9-6b04245c9ed0
# ╟─9b254664-72a6-4030-83b0-27db11c222ba
# ╠═60370f6d-51f3-436f-ab32-400de8d98f2b
# ╠═d62d50a8-dd17-4d25-b353-fdfedb456ff0
# ╠═a77205bb-eab2-4a2a-90b7-2f4844f039de
# ╟─b3455451-651e-4f70-87f4-08f9db13d9b5
# ╟─35b79339-bb2e-4987-a0f4-583a732c71fd
# ╟─e35a7298-79fa-42b3-a2c7-faec6b5c9af4
# ╠═86fac660-6ef4-4041-80d8-61eb6edc2a39
# ╟─005f8151-6b2b-4476-b92d-e63d0d216325
# ╟─b71e334a-8452-453c-a7a4-6d5f6eafd0bc
# ╠═13806613-6eca-4b9b-8216-894b25ff1b01
# ╠═9a8f53e1-4900-4b64-a09f-5a9b8fd5a0c5
# ╠═05233816-9613-4cb6-90f0-bac0f01816c0
# ╟─b7f35164-9e08-49ea-aa40-f968fbbb92c2
# ╠═b5e14014-12bf-4004-bccc-2f3c759200b5
# ╠═88edf2c7-3c5a-4b90-9c50-5d448a358ce8
# ╠═786eb037-e39e-41e1-9a63-6106800cb306
# ╠═fe51b3b3-31bc-454f-bc09-db5a7e6c6941
# ╠═2867365f-f7fe-460c-b8fd-eb8640f7dc73
# ╠═61cffd44-a5e8-4223-85ef-eb4f2b22de2b
# ╠═60324b4b-e880-4ea3-9831-fca775eb22ab
# ╠═7c4fe694-7343-4e87-aed4-8895802e34a9
# ╠═9f12b7bd-ccc2-4108-9171-1902a0c75c30
# ╠═2d55c0c1-9292-40bc-8189-f73783b28d47
# ╠═958a9d50-c61f-4d8e-89e7-7b37f0fd4246
# ╟─fef2fcda-cc64-469e-9f17-946f06a3e4ae
# ╟─4f9c5043-3485-45b8-b493-d62780134251
# ╟─197e0479-f15a-43f9-a7cf-5b5d58b5a901
# ╟─b030dcf1-fc55-4e32-be9b-6cba1d3456f0
# ╠═5a9c7e99-0e6e-4162-8ff7-bff1935fd318
# ╠═7976fbad-ae1f-40e8-961f-2ea68f752caa
# ╠═d8aeb68a-697d-4477-8b6c-14712f50a945
