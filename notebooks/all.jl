### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 9a571f4b-c7e0-4281-bd7c-91f1bdad2ada
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-qrm-phantom"
	using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
	using StatsBase: quantile!
end

# ╔═╡ 65c8d3bc-d4ec-44f3-b20a-0cf2f931a1a9
include(srcdir("masks.jl"))

# ╔═╡ da866599-e426-49d8-bf8b-0038aa0db264
BASE_PATH = "/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/Datasets/CAC Data"

# ╔═╡ 58a98b14-521f-48fa-bebd-2275f660d9d8
VENDORS = ["Canon_Aquilion_One_Vision", "GE_Revolution", "Philips_Brilliance_iCT", "Siemens_SOMATOM_Force"]

# ╔═╡ 94218a58-c19c-42d1-b22f-2aa39394e5f3
# cal_root = joinpath(dirname(dirname(pwd())), "cac-simulation", "output_new")

# ╔═╡ 7d901382-ffe6-4d77-8622-a602977528d6
# cal_path = joinpath(cal_root, "calibrations.csv")

# ╔═╡ 67f2154c-2c14-492a-a3dc-64a1201ddc7b
# cal_df = CSV.read(cal_path, DataFrame)

# ╔═╡ f044de90-0427-4164-98bc-cf828f5f6eb4
# OUTPUT = "output_new"

# ╔═╡ 8288b069-8e99-4147-90c5-b804c5867bab
# SAVE_DF = "physical.csv"

# ╔═╡ ab4af475-a031-42d0-b8b0-f20d06b1ac57
begin
	dfs_i = []
	dfs_a = []
    high_dens = []
    med_dens = []
    low_dens = []
    high_dens100 = []
    med_dens50 = []
    low_dens25 = []
    for VENDOR in VENDORS
		root_path = joinpath(BASE_PATH, VENDOR)
		dcm_path_list = dcm_list_builder(root_path)
		for path in dcm_path_list
			#---------------- Reusable Pieces ----------------#
			scan = basename(path)
			header, dcm_array, slice_thick_ori1 = dcm_reader(path)
			kV = header[tag"KVP"]

			@info VENDOR, scan
	
			# Segment Heart
			masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
		
			# Segment Calcium Rod
			calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)
			
			
			thresh = 115
			mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
				dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=0, calcium_threshold=thresh)
		
			arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
			single_arr = masked_array[:, :, slice_CCI]
			pixel_size = DICOMUtils.get_pixel_size(header)
			
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
			mask_bkg_S = create_circular_mask(size(masked_array, 1), size(masked_array, 2), (center_insert[2], center_insert[1]), 8)
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
				
			c_img = calcium_image[:, :, cal_rod_slice]
			array_filtered = abs.(mapwindow(median, c_img, (3, 3)))
			bool_arr = array_filtered .> 0
			bool_arr_erode = erode(erode(erode(erode(erode(bool_arr)))))
			hu_calcium = mean(c_img[bool_arr_erode])

			# Calibration Prep
			kvp = header[(0x0018, 0x0060)]
			# if kvp == 80.0
			# 	arr_cal = Array(cal_df[1, 2:end])
			# 	intensity_array3 = vcat(0, arr_cal)
			# elseif kvp == 100.0
			# 	arr_cal = Array(cal_df[2, 2:end])
			# 	intensity_array3 = vcat(0, arr_cal)
			# elseif kvp == 120.0
			# 	arr_cal = Array(cal_df[3, 2:end])
			# 	intensity_array3 = vcat(0, arr_cal)
			# end
						
			# density_array_calc3 = [
			# 		0
			# 		25
			# 		50
			# 		100
			# 		200
			# 		400
			# 		800
			# 	]
			# df_cal = DataFrame(:density => density_array_calc3, :intensity => intensity_array3)
			df_cal = DataFrame(:density => [0, 200], :intensity => [0, hu_calcium])
			linearRegressor = lm(@formula(intensity ~ density), df_cal)
			linearFit = predict(linearRegressor)
			m = linearRegressor.model.pp.beta0[2]
			b = linearRegressor.model.rr.mu[1]
			density(intensity) = (intensity - b) / m
			intensity(ρ) = m*ρ + b
			pixel_size = DICOMUtils.get_pixel_size(header)

			# Score background
			S_Obj_bkg = intensity(200)
			ρ_bkg = 0.2 # mg/mm^3
			
			S_bkg_large = mean(arr[ring_mask_L_bkg])
			alg_bkg_large = Integrated(arr[dilated_mask_L_bkg])
			mass_large_bkg = score(S_bkg_large, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg_large)

			S_bkg_medium = mean(arr[ring_mask_M_bkg])
			alg_bkg_medium = Integrated(arr[dilated_mask_M_bkg])
			mass_medium_bkg = score(S_bkg_medium, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg_medium)

			S_bkg_small = mean(arr[ring_mask_S_bkg])
			alg_bkg_small = Integrated(arr[dilated_mask_S_bkg])
			mass_small_bkg = score(S_bkg_small, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg_small)

			mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]
			
			## High Density
			single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
			s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
			S_Obj_HD = intensity(800)
			ρ_hd = 0.8 # mg/mm^3
			alg_L_HD = Integrated(arr[mask_L_HD_3D])
			mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
		
			## Medium Density
			single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
			s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
			S_Obj_MD = intensity(400)
			ρ_md = 0.4 # mg/mm^3
			alg_L_MD = Integrated(arr[mask_L_MD_3D])
			mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
		
			## Low Density
			single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
			s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
			S_Obj_LD = intensity(200)
			ρ_ld = 0.2 # mg/mm^3
			alg_L_LD = Integrated(arr[mask_L_LD_3D])
			mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_ld, alg_L_LD)
			
			# Score Medium Inserts
			## High Density
			single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
			s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
			alg_M_HD = Integrated(arr[mask_M_HD_3D])
			mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
			
			## Medium Density
			single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
			s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
			alg_M_MD = Integrated(arr[mask_M_MD_3D])
			mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
		
			## Low Density
			single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
			s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
			alg_M_LD = Integrated(arr[mask_M_LD_3D])
			mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_ld, alg_M_LD)
			
		
			# Score Small Inserts
			## High Density
			single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
			s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
			alg_S_HD = Integrated(arr[mask_S_HD_3D])
			mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
		
			## Medium Density
			single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
			s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
			alg_S_MD = Integrated(arr[mask_S_MD_3D])
			mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
			
			## Low Density
			single_ring_mask_S_LD = Bool.(ring_mask_S_LD[:, :, 3])
			s_bkg_S_LD = mean(single_arr[single_ring_mask_S_LD])
			alg_S_LD = Integrated(arr[mask_S_LD_3D])
			mass_s_ld = score(s_bkg_S_LD, S_Obj_LD, pixel_size, ρ_ld, alg_S_LD)
			
			# Results
			density_array = [0, 200, 400, 800]
			inserts = [
				"Low Density",
				"Medium Density",
				"High Density"
			]
			ground_truth_mass_large = [
				19.6,
				39.3,
				78.5
			] # mg
			
			calculated_mass_large = [
				mass_l_ld,
				mass_l_md,
				mass_l_hd
			]
			ground_truth_mass_medium = [
				4.2,
				8.5,
				17.0
			]
			calculated_mass_medium = [
				mass_m_ld,
				mass_m_md,
				mass_m_hd
			]
			ground_truth_mass_small = [
				0.2,
				0.3,
				0.6
			]
			calculated_mass_small = [
				mass_s_ld,
				mass_s_md,
				mass_s_hd
			]
			
			df = DataFrame(
				vendor = VENDOR,
				scan = scan,
				inserts = inserts,
				ground_truth_mass_large = ground_truth_mass_large,
				calculated_mass_large = calculated_mass_large,
				ground_truth_mass_medium = ground_truth_mass_medium,
				calculated_mass_medium = calculated_mass_medium,
				ground_truth_mass_small = ground_truth_mass_small,
				calculated_mass_small = calculated_mass_small,
				mass_bkg = mass_bkg
			)
			push!(dfs_i, df)
	
			#---------------- Agatston ----------------#
			# Mask Calibration Factor
            output = calc_output(masked_array, header, slice_CCI, thresh, trues(3, 3))
            insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
            rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
            pixel_size = DICOMUtils.get_pixel_size(header)
            mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, cal_rod_slice, rows, cols, pixel_size)

			# Background
			alg = Agatston()
			overlayed_bkg_mask_L = create_mask(arr, dilated_mask_L_bkg)
			overlayed_bkg_mask_M = create_mask(arr, dilated_mask_M_bkg)
			overlayed_bkg_mask_S = create_mask(arr, dilated_mask_S_bkg)

			agat_bkg, mass_bkg_large = score(
				overlayed_bkg_mask_L,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			agat_bkg, mass_bkg_medium = score(
				overlayed_bkg_mask_M,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			agat_bkg, mass_bkg_small = score(
				overlayed_bkg_mask_S,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			mass_bkg = [mass_bkg_large, mass_bkg_medium, mass_bkg_small]

            # Score Large Inserts
            ## High Density
            
            overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
            agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg)

            ## Medium Density
            overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
            agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_cal_factor, alg)

            ## Low Density
            overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
            agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg)


            # Score Medium Inserts
            ## High Density
            overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
            agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg)

            ## Medium Density
            overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
            agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_cal_factor, alg)

            ## Low Density
            overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
            agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg)

            # Score Small Inserts
            ## High Density
            overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
            agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg)

            ## Medium Density
            overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
            agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_cal_factor, alg)

            ## Low Density
            overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
            agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg)

            # Results

            ## Agatston
            calculated_agat_large = [
                agat_l_ld,
                agat_l_md,
                agat_l_hd
            ]
            calculated_agat_medium = [
                agat_m_ld,
                agat_m_md,
                agat_m_hd
            ]
            calculated_agat_small = [
                agat_s_ld,
                agat_s_md,
                agat_s_hd
            ]

            ## Mass
            calculated_mass_large = [
                mass_l_ld,
                mass_l_md,
                mass_l_hd
            ]
            calculated_mass_medium = [
                mass_m_ld,
                mass_m_md,
                mass_m_hd
            ]
            calculated_mass_small = [
                mass_s_ld,
                mass_s_md,
                mass_s_hd
            ]

            df = DataFrame(
                vendor=VENDOR,
                scan=scan,
                inserts=inserts,
                calculated_agat_large=calculated_agat_large,
                calculated_agat_medium=calculated_agat_medium,
                calculated_agat_small=calculated_agat_small,
                ground_truth_mass_large=ground_truth_mass_large,
                calculated_mass_large=calculated_mass_large,
                ground_truth_mass_medium=ground_truth_mass_medium,
                calculated_mass_medium=calculated_mass_medium,
                ground_truth_mass_small=ground_truth_mass_small,
                calculated_mass_small=calculated_mass_small,
				mass_bkg = mass_bkg,
                mass_cal_factor=mass_cal_factor
            )
            push!(dfs_a, df)
		end
    end
end

# ╔═╡ 9eb979d2-45f5-4960-bd73-32c469281e67
md"""
# Save Results
"""

# ╔═╡ b91bbb38-3bd8-4a03-9710-8883ac9480f8
dfs_i

# ╔═╡ 594c8fdf-9bc2-4f19-a560-d30bfd436ef9
dfs_a

# ╔═╡ 4d0acfea-1ee6-4662-8f9e-af09b80f5fc9
begin
	i_path = joinpath(cd(pwd, ".."), OUTPUT, "integrated")
	if ~isdir(i_path)
	    mkpath(i_path)
	end

	a_path = joinpath(cd(pwd, ".."), OUTPUT, "agatston")
	if ~isdir(a_path)
	    mkpath(a_path)
	end
end

# ╔═╡ eacaa8e3-615d-4a98-aa5f-68524f3de424
begin
	new_df = vcat(dfs_i[1:length(dfs_i)]...)
    output_path_new = joinpath(i_path, SAVE_DF)
    CSV.write(output_path_new, new_df)

	new_df = vcat(dfs_a[1:length(dfs_a)]...)
    output_path_new = joinpath(a_path, SAVE_DF)
    CSV.write(output_path_new, new_df)
end

# ╔═╡ 75beba3e-b838-451a-b257-02ed3d26417c


# ╔═╡ Cell order:
# ╠═9a571f4b-c7e0-4281-bd7c-91f1bdad2ada
# ╠═65c8d3bc-d4ec-44f3-b20a-0cf2f931a1a9
# ╠═da866599-e426-49d8-bf8b-0038aa0db264
# ╠═58a98b14-521f-48fa-bebd-2275f660d9d8
# ╠═94218a58-c19c-42d1-b22f-2aa39394e5f3
# ╠═7d901382-ffe6-4d77-8622-a602977528d6
# ╠═67f2154c-2c14-492a-a3dc-64a1201ddc7b
# ╠═f044de90-0427-4164-98bc-cf828f5f6eb4
# ╠═8288b069-8e99-4147-90c5-b804c5867bab
# ╠═ab4af475-a031-42d0-b8b0-f20d06b1ac57
# ╟─9eb979d2-45f5-4960-bd73-32c469281e67
# ╠═b91bbb38-3bd8-4a03-9710-8883ac9480f8
# ╠═594c8fdf-9bc2-4f19-a560-d30bfd436ef9
# ╠═4d0acfea-1ee6-4662-8f9e-af09b80f5fc9
# ╠═eacaa8e3-615d-4a98-aa5f-68524f3de424
# ╠═75beba3e-b838-451a-b257-02ed3d26417c
