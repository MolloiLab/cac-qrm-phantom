### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 015f30f1-7bde-40e7-9294-2bd495a7962c
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
    using StatsBase: quantile!

end

# ╔═╡ 7e617d7e-dd3a-4fd4-a77b-bf2608629d36
TableOfContents()

# ╔═╡ 8b755f00-8176-4da8-9ab5-7d1585a6d4fb
BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"

# ╔═╡ a2871448-7b2d-41b5-9cf5-9a711049e16d
venders = ["Canon_Aquilion_One_Vision", "GE_Revolution", "Philips_Brilliance_iCT", "Siemens_SOMATOM_Force"]

# ╔═╡ ad398e5a-c7fb-4e20-8cfc-3cf50e86fd21
scans = collect(1:10)

# ╔═╡ 731e7e8d-4883-4f1f-b93a-2534b9962771
cal_path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/calibrations.csv"

# ╔═╡ 359a88ad-174e-4e90-a1c8-3d9321279b22
cal_df = CSV.read(cal_path, DataFrame)

# ╔═╡ e4f24c82-60ce-4933-b441-d504fdee078a
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
		
			# Segment Heart
			masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3)÷2)
		
			# Segment Calcium Rod
			calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)
		
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

			
			thresh = 115
			mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
				dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=0, calcium_threshold=thresh)
		
			arr = masked_array[:, :, slice_CCI-3:slice_CCI+3]
			single_arr = masked_array[:, :, slice_CCI]
			pixel_size = DICOMUtils.get_pixel_size(header)
			
			# Score Large Inserts
			
			## Pure Background
			background_mask = zeros(size(arr)...)
			background_mask[center_insert[1]-5:center_insert[1]+5, center_insert[2]-5:center_insert[2]+5, 2] .= 1
			background_mask = Bool.(background_mask)
			background_mask_dil = dilate(dilate(background_mask))
			
			single_bkg_center = Bool.(background_mask[:, :, 2])
			S_Obj_bkg = intensity(200)
			
			ring_background = background_mask_dil - background_mask
			single_ring_bkg = Bool.(ring_background[:, :, 2])
			s_bkg = mean(single_arr[single_ring_bkg])
			
			alg_bkg = Integrated(arr[background_mask])
			ρ_bkg = 0.2 # mg/mm^3
			mass_bkg = score(s_bkg, S_Obj_bkg, pixel_size, ρ_bkg, alg_bkg)
			
			## High Density
			mask_L_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_HD_3D[:, :, z] = mask_L_HD
			end
			dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
			ring_mask_L_HD = dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D)))
			single_ring_mask_L_HD = Bool.(ring_mask_L_HD[:, :, 3])
			s_bkg_L_HD = mean(single_arr[single_ring_mask_L_HD])
			S_Obj_HD = intensity(800)
			ρ_hd = 0.8 # mg/mm^3
			alg_L_HD = Integrated(arr[mask_L_HD_3D])
			mass_l_hd = score(s_bkg_L_HD, S_Obj_HD, pixel_size, ρ_hd, alg_L_HD)
		
			## Medium Density
			mask_L_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_MD_3D[:, :, z] = mask_L_MD
			end
			dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
			ring_mask_L_MD = dilate(dilate(dilate(dilate(mask_L_MD_3D)))) - dilate(dilate(dilate(mask_L_MD_3D)))
			single_ring_mask_L_MD = Bool.(ring_mask_L_MD[:, :, 3])
			s_bkg_L_MD = mean(single_arr[single_ring_mask_L_MD])
			S_Obj_MD = intensity(400)
			ρ_md = 0.4 # mg/mm^3
			alg_L_MD = Integrated(arr[mask_L_MD_3D])
			mass_l_md = score(s_bkg_L_MD, S_Obj_MD, pixel_size, ρ_md, alg_L_MD)
		
			## Low Density
			mask_L_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_LD_3D[:, :, z] = mask_L_LD
			end
			dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
			ring_mask_L_LD = dilate(dilate(dilate(dilate(mask_L_LD_3D)))) - dilate(dilate(dilate(mask_L_LD_3D)))
			single_ring_mask_L_LD = Bool.(ring_mask_L_LD[:, :, 3])
			s_bkg_L_LD = mean(single_arr[single_ring_mask_L_LD])
			S_Obj_LD = intensity(200)
			ρ_ld = 0.2 # mg/mm^3
			alg_L_LD = Integrated(arr[mask_L_LD_3D])
			mass_l_ld = score(s_bkg_L_LD, S_Obj_LD, pixel_size, ρ_ld, alg_L_LD)
			
			# Score Medium Inserts
			## High Density
			mask_M_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_HD_3D[:, :, z] = mask_M_HD
			end
			dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
			ring_mask_M_HD = dilate(dilate(dilate(dilate(dilate(mask_M_HD_3D))))) - dilate(dilate(dilate(dilate(mask_M_HD_3D))))
			single_ring_mask_M_HD = Bool.(ring_mask_M_HD[:, :, 3])
			s_bkg_M_HD = mean(single_arr[single_ring_mask_M_HD])
			alg_M_HD = Integrated(arr[mask_M_HD_3D])
			mass_m_hd = score(s_bkg_M_HD, S_Obj_HD, pixel_size, ρ_hd, alg_M_HD)
			
			## Medium Density
			mask_M_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_MD_3D[:, :, z] = mask_M_MD
			end
			dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
			ring_mask_M_MD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_MD_3D)))))
			single_ring_mask_M_MD = Bool.(ring_mask_M_MD[:, :, 3])
			s_bkg_M_MD = mean(single_arr[single_ring_mask_M_MD])
			alg_M_MD = Integrated(arr[mask_M_MD_3D])
			mass_m_md = score(s_bkg_M_MD, S_Obj_MD, pixel_size, ρ_md, alg_M_MD)
		
			## Low Density
			mask_M_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_LD_3D[:, :, z] = mask_M_LD
			end
			dilated_mask_M_LD = dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
			ring_mask_M_LD = dilate(dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))) - dilate(dilate(dilate(dilate(dilate(mask_M_LD_3D)))))
			single_ring_mask_M_LD = Bool.(ring_mask_M_LD[:, :, 3])
			s_bkg_M_LD = mean(single_arr[single_ring_mask_M_LD])
			alg_M_LD = Integrated(arr[mask_M_LD_3D])
			mass_m_ld = score(s_bkg_M_LD, S_Obj_LD, pixel_size, ρ_ld, alg_M_LD)
			
		
			# Score Small Inserts
			## High Density
			mask_S_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_HD_3D[:, :, z] = mask_S_HD
			end
			dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))))
			ring_mask_S_HD = dilate(dilate(dilate(dilate(dilate(mask_S_HD_3D))))) - dilate(dilate(dilate(dilate(mask_S_HD_3D))))
			single_ring_mask_S_HD = Bool.(ring_mask_S_HD[:, :, 3])
			s_bkg_S_HD = mean(single_arr[single_ring_mask_S_HD])
			alg_S_HD = Integrated(arr[mask_S_HD_3D])
			mass_s_hd = score(s_bkg_S_HD, S_Obj_HD, pixel_size, ρ_hd, alg_S_HD)
		
			## Medium Density
			mask_S_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_MD_3D[:, :, z] = mask_S_MD
			end
			dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))))
			ring_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D))))) - dilate(dilate(dilate(dilate(mask_S_MD_3D))))
			single_ring_mask_S_MD = Bool.(ring_mask_S_MD[:, :, 3])
			s_bkg_S_MD = mean(single_arr[single_ring_mask_S_MD])
			alg_S_MD = Integrated(arr[mask_S_MD_3D])
			mass_s_md = score(s_bkg_S_MD, S_Obj_MD, pixel_size, ρ_md, alg_S_MD)
			
			## Low Density
			mask_S_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_LD_3D[:, :, z] = mask_S_LD
			end
			dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))))
			ring_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D))))) - dilate(dilate(dilate(dilate(mask_S_LD_3D))))
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
				vender = VENDER,
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
			push!(dfs, df)
		end
	end
end

# ╔═╡ d265f1e9-4f2d-4342-99e3-1ee74d60305d
md"""
# Save Results
"""

# ╔═╡ cf2eab60-7418-49fc-a2ac-ba9397a2d8e2
new_df = vcat(dfs[1:length(dfs)]...)

# ╔═╡ a5091569-f68b-48c4-a5d4-f81c6daad5c1
if ~isdir(string(cd(pwd, "..") , "/data/output"))
	mkdir(string(cd(pwd, "..") , "/data/output"))
end

# ╔═╡ 152b5890-eb8e-4f2b-a91e-107416b920b9
output_path = string(cd(pwd, "..") , "/data/output", "/integrated.csv")

# ╔═╡ 2c94457f-fa67-4b2e-9756-68d42e6ac06b
CSV.write(output_path, new_df)

# ╔═╡ Cell order:
# ╠═015f30f1-7bde-40e7-9294-2bd495a7962c
# ╠═7e617d7e-dd3a-4fd4-a77b-bf2608629d36
# ╠═8b755f00-8176-4da8-9ab5-7d1585a6d4fb
# ╠═a2871448-7b2d-41b5-9cf5-9a711049e16d
# ╠═ad398e5a-c7fb-4e20-8cfc-3cf50e86fd21
# ╠═731e7e8d-4883-4f1f-b93a-2534b9962771
# ╠═359a88ad-174e-4e90-a1c8-3d9321279b22
# ╠═e4f24c82-60ce-4933-b441-d504fdee078a
# ╟─d265f1e9-4f2d-4342-99e3-1ee74d60305d
# ╠═cf2eab60-7418-49fc-a2ac-ba9397a2d8e2
# ╠═a5091569-f68b-48c4-a5d4-f81c6daad5c1
# ╠═152b5890-eb8e-4f2b-a91e-107416b920b9
# ╠═2c94457f-fa67-4b2e-9756-68d42e6ac06b
