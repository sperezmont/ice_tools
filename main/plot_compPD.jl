
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
locdata = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d03_LateralBC/MARINE_PD/v1.0/"           # path to locate the data
locplot = "/home/sergio/entra/proyects/d03_LateralBC/plots/MARINE_PD/"                                         # path to save plots
expnames = ["spinup_MARINE_PD_d0.1", "spinup_MARINE_PD_d0.5", "spinup_MARINE_PD_d1.0"] # each experiment is a column
plot_name = "spinup_MARINE_PD_v1.0"

ylimits = [(0, 4),       # iter_redo
    (0, 1),                # dt_now
    (0, 20)]               # pc_eta          

# Load
data, max_vals = [], []
varnames = ["time", "H_ice_pd_err", "uxy_s_pd_err", "rmse_H", "rmse_uxy", "V_sle"]
varunits = ["yrs", "m", "m/yr", "m", "m/yr", "m SLE"]
varlevels = Dict("H_ice_pd_err" => -2000:100:2000, "uxy_s_pd_err" => -5000:500:5000)
varmaps = Dict("H_ice_pd_err" => :redsblues, "uxy_s_pd_err" => :bluesreds)

for i in 1:length(varnames)
    data_i, max_i = [], []
    for j in 1:length(expnames)
        if varnames[i] != "V_sle"
            data_ij = ncread(locdata * expnames[j] * "/yelmo2D.nc", varnames[i])
        else
            data_ij = ncread(locdata * expnames[j] * "/yelmo1D.nc", varnames[i])
        end
        push!(data_i, data_ij)
        push!(max_i, maximum(data_ij))
    end
    push!(data, data_i)
    push!(max_vals, max_i)
end

# Plot
fgsz = (1200 * length(expnames), 2400)
fig = Figure(resolution=fgsz)
fntsz = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)
fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
set_theme!(fontsize_theme)

for i in 1:2 # i = row, j = columns
    for j in 1:length(expnames)
        if i == 1
            ax = Axis(fig[i, j], title=expnames[j], titlesize=0.8 * fntsz)
        elseif i == 2
            ax = Axis(fig[i, j])
        end
        update_theme!()
        hidespines!(ax, :t, :b, :r, :l)
        hidedecorations!(ax, label=false)

        contourf!(ax, data[i+1][j], colormap=varmaps[varnames[i+1]], levels=varlevels[varnames[i+1]])
        
        (j == length(expnames)) && (Colorbar(fig[i, end+1], maps[end], height=Relative(2 / 3), width=30, label=varnames[i+1], ticklabelsize=0.8 * fntsz, ticks=(varlevels[varnames[i+1]], string.(varlevels[varnames[i+1]]))))

    end
end

save(locplot * plot_name * "-compPD_plot.png", fig)