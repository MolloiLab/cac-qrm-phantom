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

# ╔═╡ e6bf9991-b066-4d05-8de7-32822573c725
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("Revise")
		Pkg.add("PlutoUI")
		Pkg.add("ImageFiltering")
		Pkg.add("Images")
		Pkg.add("ImageMorphology")
		Pkg.add("ImageSegmentation")
		Pkg.add("ImageComponentAnalysis")
		Pkg.add("DataFrames")
		Pkg.add("Statistics")
		Pkg.add("CairoMakie")
		Pkg.add("DataStructures")
		Pkg.add("LinearAlgebra")
		Pkg.add(url="https://github.com/JuliaHealth/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/PhantomSegmentation.jl")
	end

	using Revise
	using PlutoUI
	using ImageFiltering
	using Images
	using ImageMorphology
	using ImageSegmentation
	using ImageComponentAnalysis
	using DataFrames
	using Statistics
	using CairoMakie
	using DataStructures
	using LinearAlgebra
	using DICOM
	using DICOMUtils
	using PhantomSegmentation
end

# ╔═╡ f203a4e1-03ef-482e-b70f-a989e94c1f08
TableOfContents()

# ╔═╡ e14721a3-25a6-487f-ac02-82fd56c43003
md"""
## `dcm_list_builder`
"""

# ╔═╡ 95da9cb9-03f5-4ef5-8588-69968638e96b
function dcm_list_builder(path)
    dcm_path_list = []
    for (dirpath, dirnames, filenames) in walkdir(path, topdown=true)
        if (dirpath in dcm_path_list) == false
            for filename in filenames
                try
                    tmp_str = string(dirpath, "/", filename)
                    ds = dcm_parse(tmp_str)
                    if (dirpath in dcm_path_list) == false
                        push!(dcm_path_list, dirpath)
					end
				catch
                    nothing
				end
			end
        else
                nothing
		end
	end
    return dcm_path_list
end

# ╔═╡ fc092ed9-e06c-47ec-a994-7b5626733aef
md"""
## `dcm_reader`
"""

# ╔═╡ b1a70a44-cc19-4531-9883-e3159fca9cdc
function dcm_reader(dcm_path)
    dcm_files = []
    for (dirpath, dirnames, filenames) in walkdir(dcm_path, topdown=false)
        for filename in filenames
            try
                if (filename == "DIRFILE") == false   
                    dcm_file = string(dirpath, "/", filename)
                    dcm_parse(dcm_file)
                    push!(dcm_files, dcm_file)
				end
			catch
				nothing
			end
		end
	end

    read_RefDs = true
	local RefDs
    while read_RefDs
        for index in range(1, length(dcm_files))
            try
                RefDs = dcm_parse(dcm_files[index])
                read_RefDs = false
                break
			catch
                nothing
			end
		end
	end

	header = RefDs.meta
	slice_thick_ori = header[(0x0018, 0x0050)]
	rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    
    ConstPixelDims = (rows, cols, length(dcm_files))
    dcm_array = zeros(ConstPixelDims...)

    instances = []    
    for filenameDCM in dcm_files
        try
            ds = dcm_parse(filenameDCM)
			head = ds.meta
			InstanceNumber = head[(0x0020, 0x0013)]
            push!(instances, InstanceNumber)
		catch
            nothing
		end
	end
    sort!(instances)
	instances = unique(instances)
    index = 0
    for filenameDCM in dcm_files
        try
            ds = dcm_parse(filenameDCM)
			head = ds.meta
			InstanceNumber = head[(0x0020, 0x0013)]
			index = findall(x -> x==InstanceNumber, instances)
			pixel_array = head[(0x7fe0, 0x0010)]
            dcm_array[:, :, index] = Int64.(pixel_array)
            index += 1
		catch
            nothing
		end
	end
	
    RescaleSlope = header[(0x0028, 0x1053)]
	RescaleIntercept = header[(0x0028, 0x1052)]
    dcm_array = dcm_array .* RescaleSlope .+ RescaleIntercept
    return RefDs.meta, dcm_array, slice_thick_ori
end

# ╔═╡ dac055b7-e90b-434a-bfb7-a8c24fcdfb28
md"""
## Load DICOMs
"""

# ╔═╡ 9f634fec-f4c3-48ca-8bfd-53f786e9c3e4
root_path = "/Users/daleblack/Google Drive/Datasets/Stanford-Motion/"

# ╔═╡ d8f0ac79-d5d0-4ce2-9d63-478e2750e109
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 70280077-b2f3-4eca-bf64-1514d80ab05d
md"""
# Whole heart mask
"""

# ╔═╡ 62286abc-b4c4-4bd4-9624-b2124fd0a55a
md"""
## `find_circle`
"""

# ╔═╡ 7dfadf75-ab5a-4328-be65-f6047323b158
function find_circle(point_1, point_2, point_3)
    x1, y1 = point_1
    x2, y2 = point_2
    x3, y3 = point_3
    
    x12 = x1 - x2 
    x13 = x1 - x3  
    y12 = y1 - y2  
    y13 = y1 - y3 
    y31 = y3 - y1  
    y21 = y2 - y1
    x31 = x3 - x1  
    x21 = x2 - x1 
 
    sx13 = x1^2 - x3^2  
    sy13 = y1^2 - y3^2
    sx21 = x2^2 - x1^2  
    sy21 = y2^2 - y1^2  
  
    f = (((sx13) * (x12) + (sy13) * (x12) + (sx21) * (x13) + (sy21) * (x13)) ÷ (2 * ((y31) * (x12) - (y21) * (x13)))) 
              
    g = (((sx13) * (y12) + (sy13) * (y12) + (sx21) * (y13) + (sy21) * (y13)) ÷ (2 * ((x31) * (y12) - (x21) * (y13))))  
  
    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0 where center is (h = -g, k = -f)  
    center_insert = [-g, -f]

    return center_insert
end

# ╔═╡ 9a422b5a-62cc-44e6-af52-7065d09643ad
find_circle([309, 309], [312, 200], [155, 155])

# ╔═╡ 902a5fbf-f56b-4644-8d61-e3df1eab6818
md"""
## `mask_heart`
"""

# ╔═╡ 066847a3-ee35-4361-aa96-b6483026adfb
"""
    mask_heart(
        header; 
        array_used=nothing, 
        radius_val=88, 
        slice_used_center=nothing
        )
Given a QRM Phantom with heart insert, this function will create a mask of the whole heart for
image processing purposes.
"""
function mask_heart(header, array_used, slice_used_center; radius_val=88)
    pixel_size = PhantomSegmentation.get_pixel_size(header)

    radius = (radius_val / 2) / pixel_size[1]
    central_image = copy(array_used[:, :, slice_used_center])
    central_image = Int.(central_image .< -200)
    kern = Int.(round(5 / pixel_size[1]))
    if kern % 2 == 0
        kern += 1
    end
    central_image = mapwindow(median, central_image, (kern, kern))
    center = [size(central_image, 1) ÷ 2, size(central_image, 2) ÷ 2]
    a = copy(central_image)
    local point_1
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] + index, center[2] + index] == 1 &&
            central_image[center[1] + index, center[2] + index + 5] == 1
        )
            point_1 = [center[1] + index, center[2] + index]
            break
        else
            a[center[1] + index, center[2] + index] = 2
        end
    end

    local point_2
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] + index, center[2] - index] == 1 &&
            central_image[center[1] + index, center[2] - index - 5] == 1
        )
            point_2 = [center[1] + index, center[2] - index]
            break
        else
            a[center[1] + index, center[2] - index] = 2
        end
    end

    local point_3
    for index in 1:(size(central_image, 2) ÷ 2)
        if (
            central_image[center[1] - index, center[2] - index] == 1 &&
            central_image[center[1] - index, center[2] - index - 5] == 1
        )
            point_3 = [center[1] - index, center[2] - index]
            break
        else
            a[center[1] - index, center[2] - index] = 2
        end
    end

    center_insert = find_circle(point_1, point_2, point_3)
    rows, cols = Int(header[(0x0028, 0x0010)]), Int(header[(0x0028, 0x0011)])
    Y, X = collect(1:rows), collect(1:cols)'
    dist_from_center = @. sqrt((X - center_insert[2])^2 + (Y - center_insert[1])^2)

    mask = dist_from_center .<= radius[1]
    masked_array = zeros(size(array_used))
    for index in 1:size(array_used, 3)
        masked_array[:, :, index] = array_used[:, :, index] .* mask
    end

    return masked_array, center_insert, mask
end

# ╔═╡ ca6a4423-70c0-4948-ba85-ae8aab9eaeb0
pth = dcm_path_list[62]

# ╔═╡ e3ec043f-ac96-47ef-b1ca-6dd80996dd6c
dcm_reader(pth)[2]

# ╔═╡ 7fe7b265-900b-4dd3-b212-e37445153d2c
header, dcm_array, slice_thick_ori = dcm_reader(pth);

# ╔═╡ 6ca97337-1abb-4ae1-b7e5-b75ca82f2757
heatmap(transpose(dcm_array[:, :, 3]), colormap=:grays)

# ╔═╡ 2ddb096f-c2e8-44d5-9e5f-1cf4f8393221
maximum(dcm_array)

# ╔═╡ 982e7f8e-d057-4d67-9580-ed927295623a
minimum(dcm_array)

# ╔═╡ 5a81360d-82e7-424e-9aba-c976b7e347b9
masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2);

# ╔═╡ 45e24bcc-fbc2-4bc1-ad3b-970d93138bd9
center_insert

# ╔═╡ e9a8bf5d-5fcf-4fcf-9079-01f38970515a
begin
	fig2 = Figure()
	
	ax2 = Makie.Axis(fig2[1, 1])
	ax2.title = "Mask Array"
	heatmap!(transpose(mask), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig2
end

# ╔═╡ bef2fd59-68f9-4c98-a27b-cc0d0628eef2
begin
	fig3 = Figure()
	
	ax3 = Makie.Axis(fig3[1, 1])
	ax3.title = "Masked DICOM Array"
	heatmap!(transpose(masked_array[:, :, 23]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig3
end

# ╔═╡ c16f1b87-9d6a-4778-bc1e-46caa060cdee
begin
	fig = Figure()
	
	ax = Makie.Axis(fig[1, 1])
	ax.title = "Raw DICOM Array"
	heatmap!(transpose(dcm_array[:, :, 15]), colormap=:grays)
	scatter!(center_insert[2]:center_insert[2], center_insert[1]:center_insert[1], markersize=10, color=:red)
	fig
end

# ╔═╡ 4a9ebfb1-dbe1-417e-8193-d29de0b02579
@bind a2 PlutoUI.Slider(1:size(masked_array, 3), default=10, show_value=true)

# ╔═╡ cc6e3e96-ff16-4180-90bd-3a9c43c967fe
heatmap(transpose(masked_array[:, :, a2]), colormap=:grays)

# ╔═╡ 068fddda-7351-4d40-9ea7-90cc63a7fdf3
md"""
# Calcium inserts mask
"""

# ╔═╡ d8396c9e-4641-44b9-a3cf-37274201a04d
md"""
## `angle_calc`
"""

# ╔═╡ e470ba06-e8a2-4ab8-b95f-50c41e9e5d6d
function angle_calc(side1, side2)
    #Calculate angle between two sides of rectangular triangle
    if side1 == 0
        angle = 0
	elseif side2 == 0
        angle = π / 2
    else
        angle = atan(side1 / side2)
	end
    
    return angle
end

# ╔═╡ 04cf1650-3430-4d1f-ad9e-bffae4f78f5a
angle_calc(4, 3)

# ╔═╡ d96bf3a8-7b71-4f7d-89b8-7466e91f6b18
md"""
## `create_circular_mask`
"""

# ╔═╡ 32e06755-79a4-4298-aeee-1353814dbb36
function create_circular_mask(h, w, center_circle, radius_circle)
	Y, X = collect(1:h), collect(1:w)'
    dist_from_center = sqrt.((X .- center_circle[1]).^2 .+ (Y .- center_circle[2]).^2)

    mask = dist_from_center .<= radius_circle
    
    return mask
end

# ╔═╡ 67ace303-7a6c-4128-8b82-09d521c4d466
mask1 = create_circular_mask(40, 40, [20, 20], 1);

# ╔═╡ 2ad3c110-6f51-4639-a3af-bbcf03b39993
heatmap(mask1, colormap=:grays)

# ╔═╡ 98b5f041-9f42-47a8-aa3a-5200b3b23308
md"""
## Find key slices
"""

# ╔═╡ d64bcc7b-c8ab-4350-a631-f71fa1657365
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

# ╔═╡ 96978b75-6d6c-4c5a-8adc-55c9f04929e8
filenm = splitdir(pth)[2]

# ╔═╡ 291f98e8-f8bc-4800-b3df-ca0b44edbe95
slices1, slices2 = find_slices(filenm)

# ╔═╡ 4edd69db-2df1-4574-b8bf-88bb432c0e25
md"""
## `calc_output`
"""

# ╔═╡ 8e3a1552-371d-4c6e-846c-10c7f40cb1af
"""
	calc_output(dcm_array, CCI_slice, calcium_threshold=130, comp_connect=trues(3, 3))
Calculate the output of a dcm_array
"""
function calc_output(
    dcm_array, header, slices, calcium_threshold=130, comp_connect=trues(3, 3)
)
    # Actual scoring for CCI insert
    # First step is to remove slices without calcium from arrays
    PixelSpacing = PhantomSegmentation.get_pixel_size(header)
    SliceThickness = header[(0x0018, 0x0050)]
    CCI_min = slices[1] - 1
    CCI_max = slices[2] + 2
    central_CCI = Int(round(CCI_max - CCI_min) / 2)

    if CCI_min ≤ 0
        CCI_min = 1
    end
    if CCI_max > size(dcm_array, 3)
        CCI_max = size(dcm_array, 3)
    end

    CCI_array = copy(dcm_array[:, :, CCI_min:CCI_max])

    image_kernel = Int(round(3 / PixelSpacing[1]))
    if image_kernel % 2 == 0
        image_kernel += 1
    end

    CCI_array_binary = copy(CCI_array)
    CCI_array_binary = Int.(CCI_array_binary .> 1.0 * calcium_threshold)
    inp =
        CCI_array_binary[:, :, central_CCI - 1] +
        CCI_array_binary[:, :, central_CCI] +
        CCI_array_binary[:, :, central_CCI + 1]
    components = ImageComponentAnalysis.label_components(inp, comp_connect)
    a1 = analyze_components(components, BasicMeasurement(; area=true, perimeter=true))
    a2 = analyze_components(components, BoundingBox(; box_area=true))
    df = leftjoin(a1, a2; on=:l)
    centroids = []
    for row in eachrow(df)
        indices = row[:box_indices]
        x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
        y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
        push!(centroids, (x_point, y_point))
    end

    centroids = deleteat!(centroids, 1)

    i1 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI - 1], (image_kernel, image_kernel)
    )
    i2 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI], (image_kernel, image_kernel)
    )
    i3 = mapwindow(
        median, CCI_array_binary[:, :, central_CCI + 1], (image_kernel, image_kernel)
    )

    image_for_center = i1 + i2 + i3

    components2 = ImageComponentAnalysis.label_components(image_for_center, comp_connect)
    components2 = Int.(components2 .> 0)
    components2 = ImageComponentAnalysis.label_components(components2, comp_connect)

    b1 = analyze_components(components2, BasicMeasurement(; area=true, perimeter=true))
    b2 = analyze_components(components2, BoundingBox(; box_area=true))
    df2 = leftjoin(b1, b2; on=:l)
    centroids2 = []
    for row in eachrow(df2)
        indices = row[:box_indices]
        x_point = ((indices[1][end] - indices[1][1]) ÷ 2) + indices[1][1]
        y_point = ((indices[2][end] - indices[2][1]) ÷ 2) + indices[2][1]
        push!(centroids2, (y_point, x_point))
    end

    output = length(unique(components2)), components2, df2, centroids2
    return output
end

# ╔═╡ 4e01a4fa-7c8b-4d96-bf53-02c0e071395c
output = calc_output(masked_array, header, slices1, 115);

# ╔═╡ 2ecf6bd2-edbc-42d4-9cc0-af5d42db9f90
output1 = calc_output(masked_array, header, slices1, 115);

# ╔═╡ bcd2eee1-44f6-4722-b380-f9a11394bf3d
output2= calc_output(masked_array, header, slices2, 115);

# ╔═╡ 31f12f7e-437f-45b7-b9e4-e662db2cdaff
md"""
## Mask Inserts
"""

# ╔═╡ 652ac4a2-a8a6-4900-bb1d-5a6af581d1c4
function mask_inserts(dcm_array, slices1, slices2; threshold=115, radius=5)
	output1 = calc_output(masked_array, header, slices1, threshold)
	output1 = calc_output(masked_array, header, slices1, threshold)
	center1 = output1[4][1]
	center2 = output2[4][1]

	mask1 = create_circular_mask(size(dcm_array)[1:2]..., center1, radius)
	mask2 = create_circular_mask(size(dcm_array)[1:2]..., center2, radius)
	return mask1, mask2
end

# ╔═╡ f5c0b46e-4fc4-492b-9f3e-1246eff4d604
mask11, mask22 = mask_inserts(dcm_array, slices1, slices2);

# ╔═╡ cfff084b-db02-41ac-80a9-3b356c6351fc
md"""
## Visualize
"""

# ╔═╡ 6be0bf36-54e7-4095-ac70-20e719dee0d6
function collect_tuple(tuple_array)
	row_num = size(tuple_array)
	col_num = length(tuple_array[1])
	container = zeros(Int64, row_num..., col_num)
	for i in 1:length(tuple_array)
		container[i,:] = collect(tuple_array[i])
	end
	return container
end

# ╔═╡ 07b6c9aa-d4aa-4bd6-9888-c1bef55146fd
function overlay_mask_bind(mask)
	indices = findall(x -> x == 1, mask)
	indices = Tuple.(indices)
	label_array = collect_tuple(indices)
	zs = unique(label_array[:,3])
	return PlutoUI.Slider(1:length(zs), default=4, show_value=true)
end

# ╔═╡ 774d8b58-b386-402e-906e-ca64d97ae514
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
	scatter!(label_array[:, 1][indices_lbl], label_array[:, 2][indices_lbl], markersize=0.5, color=:red)
	fig
end

# ╔═╡ fba07324-68f8-4b22-bcb8-98bb407f5cee
begin
	masks_3D_1 = Array{Bool}(undef, size(dcm_array))
	for z in 1:size(dcm_array, 3)
		masks_3D_1[:, :, z] = mask11
	end
end;

# ╔═╡ 0bc80c07-1695-47c4-9dba-410bf1fa59ef
@bind v1 overlay_mask_bind(masks_3D_1)

# ╔═╡ d043938d-c39b-4a1f-9afc-4c8f99926c43
overlay_mask_plot(dcm_array, masks_3D_1, v1, "mask 1 overlayed")

# ╔═╡ 5d7c1abe-8ee5-41e5-817e-04f5156b5041
begin
	masks_3D_2 = Array{Bool}(undef, size(dcm_array))
	for z in 1:size(dcm_array, 3)
		masks_3D_2[:, :, z] = mask22
	end
end;

# ╔═╡ f0078851-b87e-4fea-8e3e-8fc35c7c791f
@bind v2 overlay_mask_bind(masks_3D_2)

# ╔═╡ 512e15b8-6336-4883-828e-ba4ec620e979
overlay_mask_plot(dcm_array, masks_3D_2, v2, "mask 2 overlayed")

# ╔═╡ Cell order:
# ╠═e6bf9991-b066-4d05-8de7-32822573c725
# ╠═f203a4e1-03ef-482e-b70f-a989e94c1f08
# ╟─e14721a3-25a6-487f-ac02-82fd56c43003
# ╠═95da9cb9-03f5-4ef5-8588-69968638e96b
# ╟─fc092ed9-e06c-47ec-a994-7b5626733aef
# ╠═b1a70a44-cc19-4531-9883-e3159fca9cdc
# ╠═e3ec043f-ac96-47ef-b1ca-6dd80996dd6c
# ╟─dac055b7-e90b-434a-bfb7-a8c24fcdfb28
# ╠═9f634fec-f4c3-48ca-8bfd-53f786e9c3e4
# ╠═d8f0ac79-d5d0-4ce2-9d63-478e2750e109
# ╠═7fe7b265-900b-4dd3-b212-e37445153d2c
# ╠═6ca97337-1abb-4ae1-b7e5-b75ca82f2757
# ╠═2ddb096f-c2e8-44d5-9e5f-1cf4f8393221
# ╠═982e7f8e-d057-4d67-9580-ed927295623a
# ╟─70280077-b2f3-4eca-bf64-1514d80ab05d
# ╟─62286abc-b4c4-4bd4-9624-b2124fd0a55a
# ╠═7dfadf75-ab5a-4328-be65-f6047323b158
# ╠═9a422b5a-62cc-44e6-af52-7065d09643ad
# ╟─902a5fbf-f56b-4644-8d61-e3df1eab6818
# ╠═066847a3-ee35-4361-aa96-b6483026adfb
# ╠═5a81360d-82e7-424e-9aba-c976b7e347b9
# ╠═45e24bcc-fbc2-4bc1-ad3b-970d93138bd9
# ╟─c16f1b87-9d6a-4778-bc1e-46caa060cdee
# ╟─e9a8bf5d-5fcf-4fcf-9079-01f38970515a
# ╠═bef2fd59-68f9-4c98-a27b-cc0d0628eef2
# ╠═ca6a4423-70c0-4948-ba85-ae8aab9eaeb0
# ╟─4a9ebfb1-dbe1-417e-8193-d29de0b02579
# ╠═cc6e3e96-ff16-4180-90bd-3a9c43c967fe
# ╟─068fddda-7351-4d40-9ea7-90cc63a7fdf3
# ╟─d8396c9e-4641-44b9-a3cf-37274201a04d
# ╠═e470ba06-e8a2-4ab8-b95f-50c41e9e5d6d
# ╠═04cf1650-3430-4d1f-ad9e-bffae4f78f5a
# ╟─d96bf3a8-7b71-4f7d-89b8-7466e91f6b18
# ╠═32e06755-79a4-4298-aeee-1353814dbb36
# ╠═67ace303-7a6c-4128-8b82-09d521c4d466
# ╠═2ad3c110-6f51-4639-a3af-bbcf03b39993
# ╟─98b5f041-9f42-47a8-aa3a-5200b3b23308
# ╠═d64bcc7b-c8ab-4350-a631-f71fa1657365
# ╠═96978b75-6d6c-4c5a-8adc-55c9f04929e8
# ╠═291f98e8-f8bc-4800-b3df-ca0b44edbe95
# ╟─4edd69db-2df1-4574-b8bf-88bb432c0e25
# ╠═8e3a1552-371d-4c6e-846c-10c7f40cb1af
# ╠═4e01a4fa-7c8b-4d96-bf53-02c0e071395c
# ╠═2ecf6bd2-edbc-42d4-9cc0-af5d42db9f90
# ╠═bcd2eee1-44f6-4722-b380-f9a11394bf3d
# ╟─31f12f7e-437f-45b7-b9e4-e662db2cdaff
# ╠═652ac4a2-a8a6-4900-bb1d-5a6af581d1c4
# ╠═f5c0b46e-4fc4-492b-9f3e-1246eff4d604
# ╟─cfff084b-db02-41ac-80a9-3b356c6351fc
# ╟─6be0bf36-54e7-4095-ac70-20e719dee0d6
# ╟─07b6c9aa-d4aa-4bd6-9888-c1bef55146fd
# ╟─774d8b58-b386-402e-906e-ca64d97ae514
# ╠═fba07324-68f8-4b22-bcb8-98bb407f5cee
# ╟─0bc80c07-1695-47c4-9dba-410bf1fa59ef
# ╠═d043938d-c39b-4a1f-9afc-4c8f99926c43
# ╠═5d7c1abe-8ee5-41e5-817e-04f5156b5041
# ╟─f0078851-b87e-4fea-8e3e-8fc35c7c791f
# ╠═512e15b8-6336-4883-828e-ba4ec620e979
