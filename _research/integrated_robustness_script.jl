### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 2bac4509-d9b1-4e45-a70f-3aebc343f3c3
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring
    using StatsBase: quantile!

end

# ╔═╡ 060fc83c-6dbc-48f3-bdb6-eaec53adc0ac
TableOfContents()

# ╔═╡ 6df7f84a-629a-4878-b5e1-2ea2fea61b22
cal_path = "/Users/daleblack/Google Drive/dev/MolloiLab/cac-simulation/output_new/calibrations.csv"

# ╔═╡ 1e86177c-e701-4dda-b77f-60d0dfcd6924
cal_df = CSV.read(cal_path, DataFrame)

# ╔═╡ ce74c835-7f05-4e82-88f2-3fc74ad94917
arr = Array(cal_df[1, 2:end])

# ╔═╡ cc58784a-9ec4-4054-825f-edc993721898
intens = [0]

# ╔═╡ 6198f33e-5249-4812-b5dd-0033ce7d241b
vcat(intens, arr)

# ╔═╡ ee19154c-7af9-44b2-a5f0-cf31e5a2ae4a
BASE_PATH = "/Users/daleblack/Google Drive/Datasets/"

# ╔═╡ 2c3ca39a-f9ac-49f6-afba-4e2e027d370e
venders = ["Siemens_SOMATOM_Force_extra"]

# ╔═╡ 6f2a5bdc-e868-45c7-a34c-87c4aa5e6ecc
scans = collect(1:12)

# ╔═╡ 6becd184-323c-4d45-8f64-7d2c8d8e78d2
dir = readdir(string(BASE_PATH, "/", venders[1]))

# ╔═╡ 1d01b89b-728f-410a-8914-e96931ea90e8
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
            masked_array, center_insert, mask = mask_heart(header, dcm_array, size(dcm_array, 3) ÷ 2)

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
            intensity(ρ) = m * ρ + b
            pixel_size = DICOMUtils.get_pixel_size(header)

            thresh = 115
            mask_L_HD, mask_M_HD, mask_S_HD, mask_L_MD, mask_M_MD, mask_S_MD, mask_L_LD, mask_M_LD, mask_S_LD = mask_inserts(
                dcm_array, masked_array, header, slice_CCI, center_insert; angle_factor=0, calcium_threshold=thresh)

            arr = masked_array[:, :, slice_CCI-3:slice_CCI+3]
            single_arr = masked_array[:, :, slice_CCI]
            pixel_size = DICOMUtils.get_pixel_size(header)

            # Score Large InsertS
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
                vender=VENDER,
                scan=scan,
                inserts=inserts,
                ground_truth_mass_large=ground_truth_mass_large,
                calculated_mass_large=calculated_mass_large,
                ground_truth_mass_medium=ground_truth_mass_medium,
                calculated_mass_medium=calculated_mass_medium,
                ground_truth_mass_small=ground_truth_mass_small,
                calculated_mass_small=calculated_mass_small
            )
            push!(dfs, df)
        end
    end
end

# ╔═╡ 5da17496-f4f4-42b3-916b-421a1d88ba97
md"""
# Save Results
"""

# ╔═╡ de14cc24-6bcc-4dc0-aeba-ee170b3b45ae
new_df = vcat(dfs[1:length(dfs)]...)

# ╔═╡ 6397ec48-a14f-4fd9-8b94-a1498ce2f642
if ~isdir(string(cd(pwd, ".."), "/data/output"))
    mkdir(string(cd(pwd, ".."), "/data/output"))
end

# ╔═╡ a63c5b5f-e6be-420f-9ffd-9154e4459e78
output_path = string(cd(pwd, ".."), "/data/output", "/integrated_robustness.csv")

# ╔═╡ 7cc04dd6-179d-4b51-995d-6c1a84f56493
CSV.write(output_path, new_df)

# ╔═╡ Cell order:
# ╠═2bac4509-d9b1-4e45-a70f-3aebc343f3c3
# ╠═060fc83c-6dbc-48f3-bdb6-eaec53adc0ac
# ╠═6df7f84a-629a-4878-b5e1-2ea2fea61b22
# ╠═1e86177c-e701-4dda-b77f-60d0dfcd6924
# ╠═ce74c835-7f05-4e82-88f2-3fc74ad94917
# ╠═cc58784a-9ec4-4054-825f-edc993721898
# ╠═6198f33e-5249-4812-b5dd-0033ce7d241b
# ╠═ee19154c-7af9-44b2-a5f0-cf31e5a2ae4a
# ╠═2c3ca39a-f9ac-49f6-afba-4e2e027d370e
# ╠═6f2a5bdc-e868-45c7-a34c-87c4aa5e6ecc
# ╠═6becd184-323c-4d45-8f64-7d2c8d8e78d2
# ╠═1d01b89b-728f-410a-8914-e96931ea90e8
# ╟─5da17496-f4f4-42b3-916b-421a1d88ba97
# ╠═de14cc24-6bcc-4dc0-aeba-ee170b3b45ae
# ╠═6397ec48-a14f-4fd9-8b94-a1498ce2f642
# ╠═a63c5b5f-e6be-420f-9ffd-9154e4459e78
# ╠═7cc04dd6-179d-4b51-995d-6c1a84f56493
