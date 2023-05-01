### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ d94e58a8-1179-43f5-b8d1-1d5d9e092334
using DrWatson;

# ╔═╡ ffa96149-4d40-4548-ba65-d4bfe2345ba0
# ╠═╡ show_logs = false
@quickactivate "cac-qrm-phantom"

# ╔═╡ 635bafeb-43f8-4a40-8895-bde0d5938e3b
# ╠═╡ show_logs = false
begin
	using PlutoUI, Statistics, ImageMorphology, ImageFiltering, CSV, CSVFiles, DataFrames, GLM, DICOM, DICOMUtils, PhantomSegmentation, CalciumScoring, CairoMakie, Images
	using StatsBase: quantile!
end

# ╔═╡ b4bad713-df42-424d-91cc-dc95eda4d864
include(srcdir("masks.jl"));

# ╔═╡ f1788432-1cc2-4bb8-8f6c-c35a25864ddf
TableOfContents()

# ╔═╡ 60b8b0a7-35d0-4422-aacc-620802a907ef
BASE_PATH = "/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/Datasets/CAC Data";

# ╔═╡ 3594ec52-e98b-4ff8-af03-63ad23780a36
VENDORS = ["Canon_Aquilion_One_Vision", "GE_Revolution", "Philips_Brilliance_iCT", "Siemens_SOMATOM_Force"];

# ╔═╡ 13c0d284-0081-4194-88ef-6f53387b9fa0
VENDOR = VENDORS[3]

# ╔═╡ ba6e00e4-31c3-4a1f-ad3b-52fd0e92a6ef
begin
	root_path = joinpath(BASE_PATH, VENDOR)
	dcm_path_list = dcm_list_builder(root_path)
end

# ╔═╡ e056d922-6b16-4bf7-9b6d-e6b288c5be96
begin
	path = dcm_path_list[6]
	header, dcm_array, slice_thick_ori1 = dcm_reader(path)
end

# ╔═╡ 836c2da3-96b5-47bd-9c9f-2720fd04cc8f
@bind a PlutoUI.Slider(axes(dcm_array, 3); default=8, show_value=true)

# ╔═╡ 175c8d99-f8c8-4afb-8570-f867f05085e2
# let
# 	f = CairoMakie.Figure()
# 	ax = CairoMakie.Axis(f[1, 1])
# 	heatmap!(rotr90(dcm_array[:, :, a]); colormap=:grays)
# 	hidedecorations!(ax)
# 	save(plotsdir("rsna", "slice_$(basename(path)).png"), f)
# 	f
# end

# ╔═╡ Cell order:
# ╠═d94e58a8-1179-43f5-b8d1-1d5d9e092334
# ╠═ffa96149-4d40-4548-ba65-d4bfe2345ba0
# ╠═635bafeb-43f8-4a40-8895-bde0d5938e3b
# ╠═b4bad713-df42-424d-91cc-dc95eda4d864
# ╠═f1788432-1cc2-4bb8-8f6c-c35a25864ddf
# ╠═60b8b0a7-35d0-4422-aacc-620802a907ef
# ╠═3594ec52-e98b-4ff8-af03-63ad23780a36
# ╠═13c0d284-0081-4194-88ef-6f53387b9fa0
# ╠═ba6e00e4-31c3-4a1f-ad3b-52fd0e92a6ef
# ╠═e056d922-6b16-4bf7-9b6d-e6b288c5be96
# ╟─836c2da3-96b5-47bd-9c9f-2720fd04cc8f
# ╠═175c8d99-f8c8-4afb-8570-f867f05085e2
