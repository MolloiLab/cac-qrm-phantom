### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 71942d1f-004f-4edc-af24-9a3d151679f3
# ╠═╡ show_logs = false
using DrWatson

# ╔═╡ 05e5f4f8-9a11-40fb-9e4b-e15312ee6d59
# ╠═╡ show_logs = false
@quickactivate "cac-qrm-phantom"

# ╔═╡ 800ad242-ede0-41ae-96b2-919fd369c302
# ╠═╡ show_logs = false
begin
	using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, CSVFiles, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
	using StatsBase: quantile!
end

# ╔═╡ b5eb5c10-b8ab-4e75-8901-bf29f3eefeeb
include(srcdir("masks.jl"));

# ╔═╡ a0175742-d5d6-40a1-af64-f291ba962bb5
TableOfContents()

# ╔═╡ 8932ed56-e47b-4f01-a417-f68a07ac2f47
BASE_PATH = "/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/Datasets/CAC Data";

# ╔═╡ 31480286-9ca0-468f-893d-3f59600d2e16
VENDORS = ["Canon_Aquilion_One_Vision", "GE_Revolution", "Philips_Brilliance_iCT", "Siemens_SOMATOM_Force"];

# ╔═╡ 624dec30-100b-4ebf-9b69-8507635ab1f1
begin
    dfs_vf = []
    dfs_a = []
	dfs_s = []
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
            masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2)

            # Segment Calcium Rod
            calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)


            thresh = 115
            mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
                dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=0, calcium_threshold=thresh)

            arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
            single_arr = masked_array[:, :, slice_CCI]
            pixel_size = DICOMUtils.get_pixel_size(header)
			voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

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
			ρ_calcium = 0.2

			# Background
			hu_heart_tissue_large_bkg = mean(arr[ring_mask_L_bkg])
			mass_large_bkg = score(arr[dilated_mask_L_bkg], hu_calcium, hu_heart_tissue_large_bkg, voxel_size, ρ_calcium, VolumeFraction())

			hu_heart_tissue_medium_bkg = mean(arr[ring_mask_M_bkg])
			mass_medium_bkg = score(arr[dilated_mask_M_bkg], hu_calcium, hu_heart_tissue_medium_bkg, voxel_size, ρ_calcium, VolumeFraction())

			hu_heart_tissue_small_bkg = mean(arr[ring_mask_S_bkg])
			mass_small_bkg = score(arr[dilated_mask_S_bkg], hu_calcium, hu_heart_tissue_small_bkg, voxel_size, ρ_calcium, VolumeFraction())

			mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]

			# Score Large Inserts
			## High Density
			hu_heart_tissue_large_hd = mean(arr[ring_mask_L_HD])
			mass_l_hd = score(arr[dilated_mask_L_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

			## Medium Density
			hu_heart_tissue_large_md = mean(arr[ring_mask_L_MD])
			mass_l_md = score(arr[dilated_mask_L_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

			## Low Density
			hu_heart_tissue_large_ld = mean(arr[ring_mask_L_LD])
			mass_l_ld = score(arr[dilated_mask_L_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

			# Score Medium Inserts
			## High Density
			hu_heart_tissue_medium_hd = mean(arr[ring_mask_M_HD])
			mass_m_hd = score(arr[dilated_mask_M_HD], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())

			## Medium Density
			hu_heart_tissue_medium_md = mean(arr[ring_mask_M_MD])
			mass_m_md = score(arr[dilated_mask_M_MD], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())

			## Low Density
			hu_heart_tissue_medium_ld = mean(arr[ring_mask_M_LD])
			mass_m_ld = score(arr[dilated_mask_M_LD], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())

			# Score Small Inserts
			## High Density
			hu_heart_tissue_small_hd = mean(arr[ring_mask_S_HD])
			mass_s_hd = score(arr[dilated_mask_S_HD], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

			## Medium Density
			hu_heart_tissue_small_md = mean(arr[ring_mask_S_MD])
			mass_s_md = score(arr[dilated_mask_S_MD], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

			## Low Density
			hu_heart_tissue_small_ld = mean(arr[ring_mask_S_LD])
			mass_s_ld = score(arr[dilated_mask_S_LD], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

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
                vendor=VENDOR,
                scan=scan,
                inserts=inserts,
                ground_truth_mass_large=ground_truth_mass_large,
                calculated_mass_large=calculated_mass_large,
                ground_truth_mass_medium=ground_truth_mass_medium,
                calculated_mass_medium=calculated_mass_medium,
                ground_truth_mass_small=ground_truth_mass_small,
                calculated_mass_small=calculated_mass_small,
                mass_bkg=mass_bkg
            )
            push!(dfs_vf, df)

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
                mass_bkg=mass_bkg,
                mass_cal_factor=mass_cal_factor
            )
            push!(dfs_a, df)

			#---------------- SWCS ----------------#
			μ, σ = mean(c_img[bool_arr_erode]) / 2, std(c_img[bool_arr_erode])
			
			# Background
			alg2 = SpatiallyWeighted()

			overlayed_mask_l_bkg = create_mask(arr, dilated_mask_L_bkg)
			overlayed_mask_m_bkg = create_mask(arr, dilated_mask_M_bkg)
			overlayed_mask_s_bkg = create_mask(arr, dilated_mask_S_bkg)

			swcs_bkg_large = score(overlayed_mask_l_bkg, μ, σ, alg2)
			swcs_bkg_medium = score(overlayed_mask_m_bkg, μ, σ, alg2)
			swcs_bkg_small = score(overlayed_mask_s_bkg, μ, σ, alg2)

			swcs_bkg = [swcs_bkg_large, swcs_bkg_medium, swcs_bkg_small]

			# Score Large Inserts
			## High Density
			overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
			swcs_l_hd = score(overlayed_mask_l_hd, μ, σ, alg2)

			## Medium Density
			overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
			swcs_l_md = score(overlayed_mask_l_md, μ, σ, alg2)

			## Low Density
			overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
			swcs_l_ld = score(overlayed_mask_l_ld, μ, σ, alg2)

			# Score Medium Inserts
			## High Density
			overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
			swcs_m_hd = score(overlayed_mask_m_hd, μ, σ, alg2)

			## Medium Density
			overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
			swcs_m_md = score(overlayed_mask_m_md, μ, σ, alg2)

			## Low Density
			overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
			swcs_m_ld = score(overlayed_mask_m_ld, μ, σ, alg2)

			# Score Small Inserts
			## High Density
			overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
			swcs_s_hd = score(overlayed_mask_s_hd, μ, σ, alg2)

			## Medium Density
			overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
			swcs_s_md = score(overlayed_mask_s_md, μ, σ, alg2)

			## Low Density
			overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
			swcs_s_ld = score(overlayed_mask_s_ld, μ, σ, alg2)

			# Results
			calculated_swcs_large = [swcs_l_ld, swcs_l_md, swcs_l_hd]
			calculated_swcs_medium = [swcs_m_ld, swcs_m_md, swcs_m_hd]
			calculated_swcs_small = [swcs_s_ld, swcs_s_md, swcs_s_hd]

			df = DataFrame(;
                vendor=VENDOR,
                scan=scan,
                inserts=inserts,
				ground_truth_mass_large=ground_truth_mass_large,
				ground_truth_mass_medium=ground_truth_mass_medium,
				ground_truth_mass_small=ground_truth_mass_small,
				calculated_swcs_large=calculated_swcs_large,
				calculated_swcs_medium=calculated_swcs_medium,
				calculated_swcs_small=calculated_swcs_small,
				swcs_bkg=swcs_bkg
			)
			push!(dfs_s, df)
        end
    end
end

# ╔═╡ 92401b3b-9c53-444d-aada-f444d64a291f
md"""
# Save Results
"""

# ╔═╡ 1ba26f74-5228-4ce5-b60a-c15c75c2f26b
dfs_vf

# ╔═╡ 6bc03a29-50a3-4203-ba77-9817fd9ce419
dfs_a

# ╔═╡ 40c32556-c3a7-4fda-8af7-f38360eca3fd
dfs_s

# ╔═╡ db469f7e-6f7d-4137-96cf-20c01a581770
begin
    dfs_vf_tot = vcat(dfs_vf[1:length(dfs_vf)]...)
    save(datadir("rsna", "volume_fraction", "physical.csv"), dfs_vf_tot)

    dfs_a_tot = vcat(dfs_a[1:length(dfs_a)]...)
    save(datadir("rsna", "agatston", "physical.csv"), dfs_a_tot)
	
	dfs_s_tot = vcat(dfs_s[1:length(dfs_s)]...)
    save(datadir("rsna", "swcs", "physical.csv"), dfs_s_tot)
end

# ╔═╡ Cell order:
# ╠═71942d1f-004f-4edc-af24-9a3d151679f3
# ╠═05e5f4f8-9a11-40fb-9e4b-e15312ee6d59
# ╠═800ad242-ede0-41ae-96b2-919fd369c302
# ╠═b5eb5c10-b8ab-4e75-8901-bf29f3eefeeb
# ╠═a0175742-d5d6-40a1-af64-f291ba962bb5
# ╠═8932ed56-e47b-4f01-a417-f68a07ac2f47
# ╠═31480286-9ca0-468f-893d-3f59600d2e16
# ╠═624dec30-100b-4ebf-9b69-8507635ab1f1
# ╟─92401b3b-9c53-444d-aada-f444d64a291f
# ╠═1ba26f74-5228-4ce5-b60a-c15c75c2f26b
# ╠═6bc03a29-50a3-4203-ba77-9817fd9ce419
# ╠═40c32556-c3a7-4fda-8af7-f38360eca3fd
# ╠═db469f7e-6f7d-4137-96cf-20c01a581770
