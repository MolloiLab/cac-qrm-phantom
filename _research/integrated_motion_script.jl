### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ e4940e54-49f7-461f-8356-09ea1a67f1bc
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
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

# ╔═╡ 3d805dd9-90aa-4487-b531-8d89a1b8099c
TableOfContents()

# ╔═╡ 142ffc04-a70b-46d8-9f14-2a54d99595b5
function find_attributes(filename)
	if occursin("rod1", filename)
		rod = "rod1"
	elseif occursin("rod2", filename)
		rod = "rod2"
	end

	if occursin("pos1", filename)
		pos = "pos1"
	elseif occursin("pos2", filename)
		pos = "pos2"
	elseif occursin("pos3", filename)
		pos = "pos3"
	elseif occursin("pos4", filename)
		pos = "pos4"
	elseif occursin("pos5", filename)
		pos = "pos5"
	end

	if occursin("10mm_s", filename)
		speed = "10"
	elseif occursin("20mm_s", filename)
		speed = "20"
	elseif occursin("30mm_s", filename)
		speed = "30"
	elseif occursin("40mm_s", filename)
		speed = "40"
	elseif occursin("50mm_s", filename)
		speed = "50"
	elseif occursin("60mm_s", filename)
		speed = "60"
	elseif occursin("0mm_s", filename)
		speed = "0"
	end
	return rod, pos, speed
end

# ╔═╡ 7c96198b-81c3-403a-96da-8e81d8698238
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

# ╔═╡ f5a027fd-4d16-4c1d-b5e5-4538c6249f62
cal_path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac_simulation/output_new/calibrations.csv"

# ╔═╡ 14b07274-62fe-4d10-9fbd-71ff09252a4d
cal_df = CSV.read(cal_path, DataFrame)

# ╔═╡ 6b68b10d-04da-4c61-bd80-dd11c905de4f
BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"

# ╔═╡ 3bbdbc3f-1c0e-445f-b287-10954fbbd4ee
venders = ["Stanford-Motion"]

# ╔═╡ 9d2df848-28eb-4a9f-99e5-cb2f06d4724e
scans = collect(1:68)

# ╔═╡ e80c4efc-d4f9-46ac-b84a-baa898475d85
# rads = [10, 12, 14, 16, 18, 20, 22, 24, 26]

# ╔═╡ d10dce44-f539-4bd9-a4d9-d7033423b6eb
begin
	dfs = []
	for VENDER in venders
		for s in scans
			# for rad in rads
				@info s
				
				SCAN_NUMBER = s
				root_path = string(BASE_PATH, VENDER)
				dcm_path_list = dcm_list_builder(root_path)
				pth = dcm_path_list[SCAN_NUMBER]
				scan = basename(pth)
				header, dcm_array, slice_thick_ori1 = dcm_reader(pth)
				if occursin("rod1", scan)
					rod1 = true
				elseif occursin("rod2", scan)
					rod1 = false
				end

				rod, pos, speed = find_attributes(scan)
				if rod == "rod1"
					if pos == "pos1"
						rad = 22
					elseif pos == "pos2"
						rad = 14
					elseif pos == "pos3"
						rad = 18
					elseif pos == "pos4"
						rad = 14
					elseif pos == "pos5"
						rad = 20
					end
				elseif rod == "rod2"
					if pos == "pos1"
						rad = 26
					elseif pos == "pos2"
						rad = 12
					elseif pos == "pos3"
						rad = 16
					elseif pos == "pos4"
						rad = 10
					elseif pos == "pos5"
						rad = 14
					end
				end
				
				# Segment Heart
				masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2; radius_val=88)
	
				# Segment Calcium Inserts
				slices1, slices2 = find_slices(scan)
				mask1, mask2 = mask_inserts_motion(masked_array, header, slices1, slices2; threshold=115, radius=rad);
	
				# Calibration Prep
				kvp = header[(0x0018, 0x0060)]
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
							
				density_array_calc3 = [
						0
						25
						50
						100
						200
						400
						800
					]
				df_cal = DataFrame(:density => density_array_calc3, :intensity => intensity_array3)
				linearRegressor = lm(@formula(intensity ~ density), df_cal)
				linearFit = predict(linearRegressor)
				m = linearRegressor.model.pp.beta0[2]
				b = linearRegressor.model.rr.mu[1]
				density(intensity) = (intensity - b) / m
				intensity(ρ) = m*ρ + b
				pixel_size = DICOMUtils.get_pixel_size(header)
				
				# Score Inserts
				## Insert 1
				insert1_array = masked_array[:, :, slices1[1]-1:slices1[2]+1]
				mask1_3D = Array{Bool}(undef, size(insert1_array))
				for z in 1:size(insert1_array, 3)
					mask1_3D[:, :, z] = mask1
				end
				dilated_mask1 = dilate(dilate(mask1_3D))
				ring_mask1 = dilate(dilate(dilated_mask1)) - dilated_mask1;
				s_bkg1 = mean(insert1_array[Bool.(ring_mask1)])
				S_Obj1 = intensity(800)
				pixel_size = DICOMUtils.get_pixel_size(header)
				alg1 = Integrated(insert1_array[Bool.(mask1_3D)])
				ρ1 = 0.8 # mg/mm^3
				mass1 = score(s_bkg1, S_Obj1, pixel_size, ρ1, alg1)
				
				## Insert 2
				insert2_array = masked_array[:, :, slices2[1]-1:slices2[2]+1];
				mask2_3D = Array{Bool}(undef, size(insert2_array))
				for z in 1:size(insert2_array, 3)
					mask2_3D[:, :, z] = mask2
				end
				dilated_mask2= dilate(dilate(mask2_3D))
				ring_mask2 = dilate(dilate(dilated_mask2)) - dilated_mask2;
				s_bkg2 = mean(insert2_array[Bool.(ring_mask2)])
				S_Obj2 = intensity(800)
				alg2 = Integrated(insert2_array[Bool.(mask2_3D)])
				ρ2 = 0.8 # mg/mm^3
				mass2 = score(s_bkg2, S_Obj2, pixel_size, ρ2, alg2)
	
				# Results
				dens = [
					0.196
					0.380
					0.408
					0.800
				] # mg/mm^3
				vol = 196.35 # mm^3
				masses = vol .* dens
				
				inserts = [
					"insert 1",
					"insert 2"
				]
				if rod1
					ground_truth_mass = masses[3:4]
				else
					ground_truth_mass = masses[1:2]
				end
				
				calculated_mass = [
					mass2,
					mass1,
				]
				df = DataFrame(
					rod = rod,
					pos = pos,
					speed = speed,
					rad = rad,
					inserts = inserts,
					ground_truth_mass = ground_truth_mass,
					calculated_mass = calculated_mass,
				)
				push!(dfs, df)
			# end
		end
	end
end

# ╔═╡ f864bef2-34c5-47e3-9eb1-aa6f92998052
md"""
# Save Results
"""

# ╔═╡ d54f7159-bfa8-4320-8bdf-0f9b6b1f1fbe
new_df = vcat(dfs[1:length(dfs)]...)

# ╔═╡ 5cbec599-1d86-4d3d-80a0-f2cc85ee2051
if ~isdir(string(cd(pwd, "..") , "/data/output"))
	mkdir(string(cd(pwd, "..") , "/data/output"))
end

# ╔═╡ 321f83d4-6478-4e66-9ff6-fd873dcf743f
output_path = string(cd(pwd, "..") , "/data/output", "/integrated_motion_specific_radius.csv")

# ╔═╡ d4b11d06-01e8-4530-a902-70f23f7c8646
CSV.write(output_path, new_df)

# ╔═╡ Cell order:
# ╠═e4940e54-49f7-461f-8356-09ea1a67f1bc
# ╠═3d805dd9-90aa-4487-b531-8d89a1b8099c
# ╟─142ffc04-a70b-46d8-9f14-2a54d99595b5
# ╟─7c96198b-81c3-403a-96da-8e81d8698238
# ╠═f5a027fd-4d16-4c1d-b5e5-4538c6249f62
# ╠═14b07274-62fe-4d10-9fbd-71ff09252a4d
# ╠═6b68b10d-04da-4c61-bd80-dd11c905de4f
# ╠═3bbdbc3f-1c0e-445f-b287-10954fbbd4ee
# ╠═9d2df848-28eb-4a9f-99e5-cb2f06d4724e
# ╠═e80c4efc-d4f9-46ac-b84a-baa898475d85
# ╠═d10dce44-f539-4bd9-a4d9-d7033423b6eb
# ╟─f864bef2-34c5-47e3-9eb1-aa6f92998052
# ╠═d54f7159-bfa8-4320-8bdf-0f9b6b1f1fbe
# ╠═5cbec599-1d86-4d3d-80a0-f2cc85ee2051
# ╠═321f83d4-6478-4e66-9ff6-fd873dcf743f
# ╠═d4b11d06-01e8-4530-a902-70f23f7c8646
