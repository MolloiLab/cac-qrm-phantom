### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 0b5306c7-9cc9-47df-a06b-5a579d6b3763
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-qrm-phantom"
	using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, CSVFiles, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
	using StatsBase: quantile!

	include(srcdir("masks.jl"))
end

# ╔═╡ 02f12246-0a64-42e9-aabb-3028247ae230
TableOfContents()

# ╔═╡ f05b56d3-08ba-4dfd-ab49-7564ea380fa2
BASE_PATH = "/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/Datasets/CAC Data";

# ╔═╡ f1fe3b41-0876-487f-8c68-49991b3bc091
VENDORS = ["Stanford-Motion"];

# ╔═╡ 596b838c-c1ba-4b12-becc-be0f531530e4
md"""
## Main Loop
"""

# ╔═╡ 3a766ef0-c422-412a-847e-21f68fb7ac66
md"""
## Save Data
"""

# ╔═╡ 9978f7fc-6136-42ba-b199-be5ecceef846
md"""
## Appendix
"""

# ╔═╡ 2461da34-ccf0-41f1-b21e-d48c9757d474
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

# ╔═╡ 3851936a-5ce5-46a6-9580-08d3c9fa1ee2
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

# ╔═╡ bf32c1c6-6217-4f4e-8a84-9f4d4291c45e
function set_calibration(filename)
	rod1_ca = Dict(
		800 => [956, 936, 910, 937, 917],
		408 => [506, 488, 483, 501, 501]
	)
	rod2_ca = Dict(
		380 => [434, 436, 415, 379, 439],
		196 => [201, 218, 239, 216, 203]
	)
	local densities, intensities
	if occursin("rod1", filename)
		densities = [0, 408, 800]
		intensities = [0, maximum(rod1_ca[408]), maximum(rod1_ca[800])]
	elseif occursin("rod2", filename)
		densities = [0, 196, 380]
		intensities = [0, maximum(rod2_ca[196]), maximum(rod2_ca[380])]
	end
	return densities, intensities
end

# ╔═╡ 037f6fdc-4072-4276-8fa8-a28a86335608
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

# ╔═╡ 2fc1defb-fb2f-468e-a4ac-03280069d661
begin
	dfs_i = []
	dfs_a = []
    for VENDOR in VENDORS
		root_path = joinpath(BASE_PATH, VENDOR)
		dcm_path_list = dcm_list_builder(root_path)
		for path in dcm_path_list
			#---------------- Reusable Pieces ----------------#
			scan = basename(path)
			header, dcm_array, slice_thick_ori1 = dcm_reader(path)
			kV = header[tag"KVP"]
			
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
			@info VENDOR, scan
			
		
			# Segment Heart
			masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2; radius_val=88)
			
			# Segment Calcium Inserts
			slices1, slices2 = find_slices(scan)
			mask1, mask2 = mask_inserts_motion(masked_array, header, slices1, slices2; threshold=115, radius=rad);

			# Calibration Prep
			kvp = header[(0x0018, 0x0060)]
			dens, ints = set_calibration(scan)
			df_cal = DataFrame(:density => dens, :intensity => ints)
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
				vendor = VENDOR,
				rod = rod,
				pos = pos,
				speed = speed,
				rad = rad,
				inserts = inserts,
				ground_truth_mass = ground_truth_mass,
				calculated_mass = calculated_mass,
			)
			push!(dfs_i, df)
	
			#---------------- Agatston ----------------#
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
				vendor = VENDOR,
				rod = rod,
				pos = pos,
				speed = speed,
				calculated_agat = calculated_agat,
				ground_truth_mass = ground_truth_mass,
				calculated_mass = calculated_mass,
			)
			push!(dfs_a, df)
		end
    end
end

# ╔═╡ ded3e355-6737-49a5-9303-6d2aca7e76d7
begin
	dfs_i_tot = vcat(dfs_i[1:length(dfs_i)]...)
	save(datadir("output", "integrated", "motion.csv"), dfs_i_tot)

	dfs_a_tot = vcat(dfs_a[1:length(dfs_a)]...)
	save(datadir("output", "agatston", "motion.csv"), dfs_a_tot)
end

# ╔═╡ Cell order:
# ╠═0b5306c7-9cc9-47df-a06b-5a579d6b3763
# ╠═02f12246-0a64-42e9-aabb-3028247ae230
# ╠═f05b56d3-08ba-4dfd-ab49-7564ea380fa2
# ╠═f1fe3b41-0876-487f-8c68-49991b3bc091
# ╟─596b838c-c1ba-4b12-becc-be0f531530e4
# ╠═2fc1defb-fb2f-468e-a4ac-03280069d661
# ╟─3a766ef0-c422-412a-847e-21f68fb7ac66
# ╠═ded3e355-6737-49a5-9303-6d2aca7e76d7
# ╟─9978f7fc-6136-42ba-b199-be5ecceef846
# ╠═2461da34-ccf0-41f1-b21e-d48c9757d474
# ╠═3851936a-5ce5-46a6-9580-08d3c9fa1ee2
# ╠═bf32c1c6-6217-4f4e-8a84-9f4d4291c45e
# ╠═037f6fdc-4072-4276-8fa8-a28a86335608
