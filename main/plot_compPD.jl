
"""
Author: Sergio Pérez Montero\n
Date: 26.07.2022\n

Aim: This script plots important variables (wrt pd) from Yelmo in 3-var columns\n

"""
# Packages
using Pkg
Pkg.activate("ice_env")

# Dependencies
include("/home/sergio/entra/tools/ice_tools/src/ice_calcs.jl")
include("/home/sergio/entra/tools/ice_tools/src/ice_plots.jl")

# Variables
locdata = "/home/sergio/entra/models/yelmo_vers/v1.753_sergio-test/yelmox/output/ismip6/d03_LateralBC/dtt_03/"           # path to locate the data
locplot = "/home/sergio/entra/proyects/d03_LateralBC/plots/dtt_03/"                                         # path to save plots
expnames = ["spinup_dtt_03_m", "spinup_dtt_03_md0.1", "spinup_dtt_03_md1.0"] # each experiment is a column
plot_name = "spinup_dtt_03_marine"

# Load
data = []
varnames = ["H_ice_pd_err", "uxy_s_pd_err", "mask_bed", "rmse_H", "rmse_uxy", "V_sle"]
varlevels = Dict("H_ice_pd_err" => -2000:100:2000, "uxy_s_pd_err" => -5000:100:5000)
varmaps = Dict("H_ice_pd_err" => cgrad([:firebrick, :red, :white, :deepskyblue3, :navyblue]), "uxy_s_pd_err" => cgrad([:navyblue, :deepskyblue3, :white, :red, :firebrick]))

for i in 1:length(varnames)
    data_i, max_i = [], []
    for j in 1:length(expnames)
        if varnames[i] == "V_sle"
            data_ij = ncread(locdata * expnames[j] * "/yelmo1D.nc", varnames[i])[end]
        elseif varnames[i] in ["H_ice_pd_err", "uxy_s_pd_err", "mask_bed"]
            data_ij = ncread(locdata * expnames[j] * "/yelmo2D.nc", varnames[i])[:, :, end]
        else
            data_ij = ncread(locdata * expnames[j] * "/yelmo2D.nc", varnames[i])[end]
        end
        push!(data_i, data_ij)
    end
    push!(data, data_i)
end

# mask
for i in 1:length(expnames)
    mask = copy(data[3][i][:, :])
    data[1][i][mask.==0] .= -9999999999
    data[2][i][mask.==0] .= -9999999999
end

# Plot
fgsz = (1200 * length(expnames), 2400)
fig = Figure(resolution=fgsz)
fntsz1 = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)
fntsz2 = 0.02 * min(fgsz[1], fgsz[2])
fntsz = min(fntsz1, fntsz2)
fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
set_theme!(fontsize_theme)

maps1, maps2 = [], []
for j in 1:length(expnames)
    ax1 = Axis(fig[1, j], title=expnames[j], titlesize=fntsz)
    ax2 = Axis(fig[2, j])
    update_theme!()

    hidespines!(ax1, :t, :b, :r, :l)
    hidedecorations!(ax1, label=false)
    hidespines!(ax2, :t, :b, :r, :l)
    hidedecorations!(ax2, label=false)

    c1 = contourf!(ax1, data[1][j], colormap=varmaps[varnames[1]], levels=varlevels[varnames[1]])
    contour!(ax1, data[3][j], levels=[0, 3], color="black", linewidth=2)
    c2 = contourf!(ax2, data[2][j], colormap=varmaps[varnames[2]], levels=varlevels[varnames[2]])
    contour!(ax2, data[3][j], levels=[0, 3], color="black", linewidth=2)
    text!(ax2, "RMSE (H) = " * string(round(Int, data[4][j])) * "m", position=(10, 220), textsize=0.9 * fntsz)
    text!(ax2, "RMSE (u) = " * string(round(Int, data[5][j])) * "m/yr", position=(10, 210), textsize=0.9 * fntsz)
    text!(ax2, "V = " * string(round(Int, data[end][j])) * "m SLE", position=(10, 200), textsize=0.9 * fntsz)

    push!(maps1, c1)
    push!(maps2, c2)

    H_ticks, u_ticks = [-2000, -1000, 0, 1000, 2000], [-5000, -2500, 0, 2500, 5000]

    (j == length(expnames)) && (Colorbar(fig[1, j+1], maps1[end], height=Relative(1 / 3), width=30, label=varnames[1] * " (m)", ticklabelsize=0.8 * fntsz, ticks=H_ticks))
    (j == length(expnames)) && (Colorbar(fig[2, j+1], maps2[end], height=Relative(1 / 3), width=30, label=varnames[2] * " (m/yr)", ticklabelsize=0.8 * fntsz, ticks=u_ticks))
end

save(locplot * plot_name * "-compPD_plot.png", fig)