### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 2df34aa4-a3c7-11ec-3508-61f14f475394
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
    using StatsBase: quantile!

end

# ╔═╡ 61d478bf-7bde-4184-9231-fd4a8bb1e6a1
BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"

# ╔═╡ b4fd6633-56e2-4d0d-ac13-b305b0e2d1c8
venders = ["Canon_Aquilion_One_Vision", "GE_Revolution", "Philips_Brilliance_iCT", "Siemens_SOMATOM_Force"]

# ╔═╡ 64eb1515-c69e-41f0-8fdb-f93389a7003f
md"""
# Save results
"""

# ╔═╡ f5509a75-4a0d-484d-9b9d-7fbe72876cfe
if ~isdir(string(cd(pwd, ".."), "/data/output/"))
    mkdir(string(cd(pwd, ".."), "/data/output/"))
end

# ╔═╡ 09b33976-a24e-4bd7-b8a6-7fc47efb6f94
output_path = string(cd(pwd, ".."), "/data/output/agatston.csv")

# ╔═╡ a1c582c1-85c5-498f-9ecb-e70cd7819715
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ b8afde9f-a02a-4165-908e-41afe98109e8
begin
    dfs = []
    for VENDER in venders
        for s in 1:10
            SCAN_NUMBER = s
            @info s
            root_path = string(BASE_PATH, VENDER)
            dcm_path_list = dcm_list_builder(root_path)
            pth = dcm_path_list[SCAN_NUMBER]
            scan = basename(pth)
            header, dcm_array, slice_thick_ori1 = dcm_reader(pth)

            thresh = 115
            # Segment Heart
            masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2)

            # Segment Calcium Rod
            calcium_image, slice_CCI, quality_slice, cal_rod_slice = mask_rod(masked_array, header)

            # Segment Calcium Inserts
            mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
                dcm_array, masked_array, header, slice_CCI, center_insert;
                calcium_threshold=thresh
            )

            # Mask Calibration Factor
            output = calc_output(masked_array, header, slice_CCI, thresh, trues(3, 3))
            insert_centers = calc_centers(dcm_array, output, header, center_insert, slice_CCI)
            rows, cols = Int(header[tag"Rows"]), Int(header[tag"Columns"])
            pixel_size = DICOMUtils.get_pixel_size(header)
            mass_cal_factor, angle_0_200HA, water_rod_metrics = mass_calibration(masked_array, insert_centers[:Large_LD], center_insert, cal_rod_slice, rows, cols, pixel_size)

            # Score Large Inserts
            arr = masked_array[:, :, slice_CCI-2:slice_CCI+2]
            single_arr = masked_array[:, :, slice_CCI]

            ## High Density
            mask_L_HD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_L_HD_3D[:, :, z] = mask_L_HD
            end
            dilated_mask_L_HD = dilate(dilate(mask_L_HD_3D))
            alg = Agatston()
            overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD)
            agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_cal_factor, alg)

            ## Medium Density
            mask_L_MD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_L_MD_3D[:, :, z] = mask_L_MD
            end
            dilated_mask_L_MD = dilate(dilate(mask_L_MD_3D))
            overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD)
            agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_cal_factor, alg)

            ## Low Density
            mask_L_LD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_L_LD_3D[:, :, z] = mask_L_LD
            end
            dilated_mask_L_LD = dilate(dilate(mask_L_LD_3D))
            overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD)
            agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_cal_factor, alg)


            # Score Medium Inserts
            ## High Density
            mask_M_HD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_M_HD_3D[:, :, z] = mask_M_HD
            end
            dilated_mask_M_HD = dilate(dilate(dilate(dilate(mask_M_HD_3D))))
            overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD)
            agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_cal_factor, alg)

            ## Medium Density
            mask_M_MD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_M_MD_3D[:, :, z] = mask_M_MD
            end
            dilated_mask_M_MD = dilate(dilate(dilate(dilate(mask_M_MD_3D))))
            overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD)
            agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_cal_factor, alg)

            ## Low Density
            mask_M_LD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_M_LD_3D[:, :, z] = mask_M_LD
            end
            dilated_mask_M_LD = dilate(dilate(dilate(dilate(mask_M_LD_3D))))
            overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD)
            agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_cal_factor, alg)

            # Score Small Inserts
            ## High Density
            mask_S_HD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_S_HD_3D[:, :, z] = mask_S_HD
            end
            dilated_mask_S_HD = dilate(dilate(dilate(dilate(dilate((mask_S_HD_3D))))))
            overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD)
            agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_cal_factor, alg)

            ## Medium Density
            mask_S_MD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_S_MD_3D[:, :, z] = mask_S_MD
            end
            dilated_mask_S_MD = dilate(dilate(dilate(dilate(dilate(mask_S_MD_3D)))))
            overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD)
            agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_cal_factor, alg)

            ## Low Density
            mask_S_LD_3D = Array{Bool}(undef, size(arr))
            for z in 1:size(arr, 3)
                mask_S_LD_3D[:, :, z] = mask_S_LD
            end
            dilated_mask_S_LD = dilate(dilate(dilate(dilate(dilate(mask_S_LD_3D)))))
            overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD)
            agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_cal_factor, alg)

            # Results
            density_array = [0, 200, 400, 800]
            inserts = [
                "Low Density",
                "Medium Density",
                "High Density"
            ]

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
                vender=VENDER,
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
                mass_cal_factor=mass_cal_factor
            )
            push!(dfs, df)
        end
    end
end

# ╔═╡ 94bd4f44-cd03-4602-8dd5-b3f0ba99e2ed
new_df = vcat(dfs[1:length(dfs)]...)

# ╔═╡ ef536a67-d206-4c42-ae09-bfc8c9e1abfc
CSV.write(output_path, new_df)

# ╔═╡ Cell order:
# ╠═2df34aa4-a3c7-11ec-3508-61f14f475394
# ╠═61d478bf-7bde-4184-9231-fd4a8bb1e6a1
# ╠═b4fd6633-56e2-4d0d-ac13-b305b0e2d1c8
# ╠═b8afde9f-a02a-4165-908e-41afe98109e8
# ╟─64eb1515-c69e-41f0-8fdb-f93389a7003f
# ╠═94bd4f44-cd03-4602-8dd5-b3f0ba99e2ed
# ╠═f5509a75-4a0d-484d-9b9d-7fbe72876cfe
# ╠═09b33976-a24e-4bd7-b8a6-7fc47efb6f94
# ╠═ef536a67-d206-4c42-ae09-bfc8c9e1abfc
# ╟─a1c582c1-85c5-498f-9ecb-e70cd7819715
