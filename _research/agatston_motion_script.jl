### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 4e81a5d2-7135-4bb2-b304-dbf24557c0da
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")
	
	using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
	using StatsBase: quantile!
	
end

# ╔═╡ 70cbbd1b-b196-45ed-b1de-33ed13b70a47
TableOfContents()

# ╔═╡ 764a9432-22cd-4ca6-bbb6-a204b428b0a2
md"""
## Helper Functions
"""

# ╔═╡ 39c4d479-4ee5-461c-946e-8f6da5ec536d
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

# ╔═╡ bba19976-7ffd-4c66-9fb7-27f96284790b
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

# ╔═╡ ccabe5c2-31c1-4f83-a271-0a4eeb29c215
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

# ╔═╡ f5472805-c853-4066-bfb5-fad6a365fac9
function create_mask(array, mask)
	@assert size(array) == size(mask)
	idxs = findall(x -> x == true, mask)
	overlayed_mask = zeros(size(array))
	for idx in idxs
		overlayed_mask[idx] = array[idx]
	end
	return overlayed_mask
end

# ╔═╡ a3cdd8cd-7ae6-48ec-bbf2-040f5c6965a9
md"""
## Script
"""

# ╔═╡ fd792846-bb96-4f4d-9adb-7ddacab37de7
BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"

# ╔═╡ d4a3be67-1711-4ba0-8ede-fcdee8ec8898
venders = ["Stanford-Motion"]

# ╔═╡ c928de10-5f15-4289-a9ac-01989f52a592
scans = collect(1:68)

# ╔═╡ eda8e29c-1636-4eb5-8ff9-f99100dc49cd
begin
	dfs = []
	for VENDER in venders
		for s in scans
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
			mask1, mask2 = mask_inserts_motion(masked_array, header, slices1, slices2; threshold=115, radius=8);
			
			# Score Inserts
			## Insert 1
			insert1_array = masked_array[:, :, slices1[1]-1:slices1[2]+1]
			mask1_3D = Array{Bool}(undef, size(insert1_array))
			for z in 1:size(insert1_array, 3)
				mask1_3D[:, :, z] = mask1
			end
			dilated_mask1 = dilate(dilate(mask1_3D))
			pixel_size = DICOMUtils.get_pixel_size(header)
			alg = Agatston()
			overlayed_mask1 = create_mask(insert1_array, mask1_3D)
			agat1, mass1 = score(overlayed_mask1, pixel_size, avg_mass_cals, alg)
			
			## Insert 2
			insert2_array = masked_array[:, :, slices2[1]-1:slices2[2]+1];
			mask2_3D = Array{Bool}(undef, size(insert2_array))
			for z in 1:size(insert2_array, 3)
				mask2_3D[:, :, z] = mask2
			end
			dilated_mask2= dilate(dilate(mask2_3D))
			overlayed_mask2 = create_mask(insert2_array, mask2_3D)
			agat2, mass2 = score(overlayed_mask2, pixel_size, avg_mass_cals, alg)

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
			calculated_agat = [
				agat1,
				agat2,
			]
			
			df = DataFrame(
				vender = VENDER,
				rod = rod,
				pos = pos,
				speed = speed,
				calculated_agat = calculated_agat,
				ground_truth_mass = ground_truth_mass,
				calculated_mass = calculated_mass,
			)
			push!(dfs, df)
		end
	end
end

# ╔═╡ e4a8db82-f1da-42d3-a341-c96d99fbfb77
md"""
## Save Results
"""

# ╔═╡ 6ca49361-faf7-44d2-82f2-9516cc4a13d8
new_df = vcat(dfs[1:length(dfs)]...)

# ╔═╡ 85db3e70-21b2-4293-beb0-aa5ca40232d2
if ~isdir(string(cd(pwd, "..") , "/data/output"))
	mkdir(string(cd(pwd, "..") , "/data/output"))
end

# ╔═╡ 268f5d10-51f8-4265-aa9f-7a0e26ce0cd2
output_path = string(cd(pwd, "..") , "/data/output", "/agatston_motion_specific_radius.csv")

# ╔═╡ 9976fcf0-e998-4b39-9e64-13e15565f8bb
CSV.write(output_path, new_df)

# ╔═╡ Cell order:
# ╠═4e81a5d2-7135-4bb2-b304-dbf24557c0da
# ╠═70cbbd1b-b196-45ed-b1de-33ed13b70a47
# ╟─764a9432-22cd-4ca6-bbb6-a204b428b0a2
# ╟─39c4d479-4ee5-461c-946e-8f6da5ec536d
# ╟─bba19976-7ffd-4c66-9fb7-27f96284790b
# ╟─ccabe5c2-31c1-4f83-a271-0a4eeb29c215
# ╟─f5472805-c853-4066-bfb5-fad6a365fac9
# ╟─a3cdd8cd-7ae6-48ec-bbf2-040f5c6965a9
# ╠═fd792846-bb96-4f4d-9adb-7ddacab37de7
# ╠═d4a3be67-1711-4ba0-8ede-fcdee8ec8898
# ╠═c928de10-5f15-4289-a9ac-01989f52a592
# ╠═eda8e29c-1636-4eb5-8ff9-f99100dc49cd
# ╟─e4a8db82-f1da-42d3-a341-c96d99fbfb77
# ╠═6ca49361-faf7-44d2-82f2-9516cc4a13d8
# ╠═85db3e70-21b2-4293-beb0-aa5ca40232d2
# ╠═268f5d10-51f8-4265-aa9f-7a0e26ce0cd2
# ╠═9976fcf0-e998-4b39-9e64-13e15565f8bb
