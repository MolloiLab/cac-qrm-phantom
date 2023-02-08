"""
    create_textbox(loc, coef, r, rms)

Create a textbox, inside a Makie.jl plot that accounts for the sign of the coefficient

#### Inputs
- `loc`: figure location, e.g. `loc=f[1, 1]`
- `coef`: coefficients for the best fit line
- `r`: r squared values
- `rms`: rms values
"""
function create_textbox(loc, coef, r, rms)
    local tb
    if coef[1] > 0
        tb = Textbox(
            loc;
            placeholder="y = $(trunc(coef[2]; digits=2))x+$(trunc(coef[1]; digits=2)) \nr = $(trunc(r; digits=2)) \nRMSE: $(trunc(rms[1]; digits=2)) \nRMSD: $(trunc(rms[2]; digits=2))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    else
        tb = Textbox(
            loc,
            placeholder="y = $(trunc(coef[2]; digits=2))x$(trunc(coef[1]; digits=2)) \nr = $(trunc(r; digits=2)) \nRMSE: $(trunc(rms[1]; digits=2)) \nRMSD: $(trunc(rms[2]; digits=2))",
            tellheight=false,
            tellwidth=false,
            boxcolor=:white,
            halign=:left,
            valign=:top,
            fontsize=12
        )
    end
    return tb
end