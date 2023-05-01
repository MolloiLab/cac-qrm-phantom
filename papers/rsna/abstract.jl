### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 48f3957b-fdcb-4530-a864-8a737323310b
using DrWatson

# ╔═╡ 96fd9ec5-2829-401c-a6ac-bbfe36850424
# ╠═╡ show_logs = false
@quickactivate "cac-qrm-phantom"

# ╔═╡ 41f64723-b7ae-4d6f-a645-54d7c6cdf672
# ╠═╡ show_logs = false
using Images, PlutoUI, CairoMakie

# ╔═╡ 800e4e22-a925-47fa-a7d9-278c09c9f4b4
TableOfContents()

# ╔═╡ 356001e5-bf85-4862-a098-76ecd4f131c8
md"""
# Abstract
"""

# ╔═╡ 5dd3df83-99ee-4134-a9ed-717cdeec4427
md"""
## Purpose
Agatston scoring, the traditional method for measuring coronary artery calcium, is limited in its ability to accurately quantify low-density calcifications, among other things. The inaccuracy of Agatston scoring is likely due partly to the arbitrary thresholding requirement of Agatston scoring. A calcium quantification technique that removes the need for arbitrary thresholding and is more sensitive is needed. Improvements to calcium scoring will likely improve patient risk stratification and outcome.

## Materials and Methods
Volume fraction calcium mass requires no thresholding and includes all calcium information within an image. This study utilized phantom images acquired by G van Praagh et al., with calcium hydroxyapatite (HA) densities in the range of 200-800 mgHAcm-3 to measure calcium according to integrated calcium mass and Agatston scoring. The calcium mass was known, which allowed for accuracy and sensitivity comparisons between volume fraction calcium mass, Agatston scoring, and spatially weighted calcium scoring. Multiple CT vendors (Canon, GE, Philips, Siemens) were used during the image acquisition phase, which provided a more robust comparison between the two calcium scoring techniques. Three calcification inserts of different diameters (1, 3, and 5 mm) and different HA densities (200, 400, and 800 mgHAcm-3) were placed within the phantom.

## Results
Volume fraction mass was less accurate than Agatston scoring for stationary scans (RMSE_{VolumeFraction} = 2.87, RMSE_{Agatson} = 4.07). Spatially weighted calcium scoring produces an arbitrary score and was, therefore, removed from the accuracy analysis. The percentage of false-negative and false-positive calcium scores was lower for volume fraction calcium mass (13.61%, 0.00%) than Agatston scoring (28.33%, 6.67%) and spatially weighted calcium scoring (25.83%, 15.83%). Volume fraction calcium mass produced 49 false-negative zero-CAC scores out of 360 total measurements. Spatially weighted calcium scoring produced 93 false-negative zero-CAC scores out of 360 total measurements. Agatston scoring produced 102 false-negative zero-CAC scores out of 360 total measurements. Volume fraction calcium mass produced zero false-positive (CAC=0) scores out of 120 total measurements. Spatially weighted calcium scoring produced 19 false-negative zero-CAC scores out of 120 total measurements. Agatston scoring produced eight false-negative zero-CAC scores out of 120 total measurements.

## Conclusions
The results of this study indicate that volume fraction calcium mass is more sensitive than Agatston and spatially weighted calcium scoring on a variety of different CT vendors.

## Clinical Significance
The substantial reduction in false-negative scores for volume fraction calcium mass is likely to improve risk-stratification for patients undergoing calcium scoring and their potential outcome.

"""

# ╔═╡ edc73ad9-9462-46fd-a807-7b057407cfd1
md"""
# Figures
"""

# ╔═╡ b6149579-0641-457c-999f-39d4b9097f02
load(plotsdir("rsna", "combined.png"))

# ╔═╡ 6571cbcf-0ff9-4b78-afd8-600ae3efc328
md"""
Fig. 1 Shows (A) sketch of the small QRM-Thorax phantom with the cardiac calcification insert phantom located in the center. (B) Axial and lateral sketch of the cardiac calcification insert phantom. (C) Slice containing the calcification inserts on the large phantom acquired by the Canon scanner. (D) Slice containing the calibration rod on the small phantom acquired by the Philips scanner.
"""

# ╔═╡ 40e13289-486a-41fe-8e90-12227563e8d7
load(plotsdir("rsna", "sensitivity_specificity.png"))

# ╔═╡ 9e7504d4-dacd-4a65-a156-9233cff16133
md"""
Fig. 2 Shows the percentage of false-negative (CAC=0) scores (A) and false-positive (CAC>0) scores (B), computed by volume fraction calcium mass (left), spatially weighted calcium scoring (middle), and Agatston scoring (right).
"""

# ╔═╡ 2ae6e432-58e1-4e32-89ca-879a5552d296
md"""
# Appendix
"""

# ╔═╡ 872a034f-4bca-42e7-928a-f661d1a3a5d0
let
	img1 = load(plotsdir("rsna", "qrm_phantom.png"))
	img2 = load(plotsdir("rsna", "qrm_insert.png"))
	img3 = load(plotsdir("rsna", "slice_Large_rep1.png"))
	img4 = load(plotsdir("rsna", "slice_iCT_small_SER_0001.png"))

	f = Figure(resolution = (1400, 1400))
	
	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(rotr90(img1), colormap=:grays)
	hidedecorations!(ax)
	
	ax = CairoMakie.Axis(f[1, 2])
	heatmap!(rotr90(img2), colormap=:grays)
	hidedecorations!(ax)

	ax = CairoMakie.Axis(f[2, 1])
	heatmap!(rotr90(img3), colormap=:grays)
	hidedecorations!(ax)

	ax = CairoMakie.Axis(f[2, 2])
	heatmap!(rotr90(img4), colormap=:grays)
	hidedecorations!(ax)

	for (label, layout) in zip(["A", "B", "C", "D"], [f[1, 1], f[1, 2], f[2, 1], f[2, 2]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            fontsize=45,
            padding=(0, 0, 0, 0),
            halign=:right,
        )
    end

	save(plotsdir("rsna", "combined.png"), f, resolution=(1400, 1400))

	f
end;

# ╔═╡ Cell order:
# ╠═800e4e22-a925-47fa-a7d9-278c09c9f4b4
# ╟─356001e5-bf85-4862-a098-76ecd4f131c8
# ╟─5dd3df83-99ee-4134-a9ed-717cdeec4427
# ╟─edc73ad9-9462-46fd-a807-7b057407cfd1
# ╟─b6149579-0641-457c-999f-39d4b9097f02
# ╟─6571cbcf-0ff9-4b78-afd8-600ae3efc328
# ╟─40e13289-486a-41fe-8e90-12227563e8d7
# ╟─9e7504d4-dacd-4a65-a156-9233cff16133
# ╟─2ae6e432-58e1-4e32-89ca-879a5552d296
# ╠═48f3957b-fdcb-4530-a864-8a737323310b
# ╠═96fd9ec5-2829-401c-a6ac-bbfe36850424
# ╠═41f64723-b7ae-4d6f-a645-54d7c6cdf672
# ╠═872a034f-4bca-42e7-928a-f661d1a3a5d0
