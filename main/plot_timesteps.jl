
"""
Author: Sergio Pérez Montero\n
Date: 22.07.2022\n

Aim: This script plots timesteps.nc important variables from Yelmo in 3-var columns\n

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
expnames = ["dtt_03_mm", "dtt_03_mmd0.1", "dtt_03_mmd1.0"] # each experiment is a column
plot_name = "dtt_03_mm"

ylimits = [(0, 4),       # iter_redo
    (0, 1),                # dt_now
    (0, 20)]               # pc_eta          

# Load
display("Plotting " * plot_name * "...")
data, max_vals = [], []
varnames = ["time", "iter_redo", "dt_now", "pc_eta"]
varunits = ["yrs", "n", "yrs", "m/yr"]
for i in 1:length(varnames)
    data_i, max_i = [], []
    for j in 1:length(expnames)
        data_ij = ncread(locdata * expnames[j] * "/timesteps.nc", varnames[i])
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
colors = [:red, :blue, :green, :orange]

for i in 1:length(varnames)-1 # i = row, j = columns
    for j in 1:length(expnames)
        if i == 1
            if j == 1
                ax = Axis(fig[i, j], title=expnames[j], titlesize=0.8 * fntsz,
                    ylabel=varnames[i+1] * " (" * varunits[i+1] * ")",
                    xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz)
            else
                ax = Axis(fig[i, j], title=expnames[j], titlesize=0.8 * fntsz,
                    xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz)
            end
        elseif i == 3
            if j == 1
                ax = Axis(fig[i, j],
                    ylabel=varnames[i+1] * " (" * varunits[i+1] * ")", xlabel="Time (yrs)",
                    xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz)
            else
                ax = Axis(fig[i, j],
                    xlabel="Time (yrs)",
                    xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz)
            end
        else
            if j == 1
                ax = Axis(fig[i, j],
                    ylabel=varnames[i+1] * " (" * varunits[i+1] * ")",
                    xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz)
            else
                ax = Axis(fig[i, j],
                    xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz)
            end
        end
        update_theme!()
        lines!(ax, data[1][j], data[i+1][j], color=colors[j], linewidth=3)

        maxy = maximum(max_vals[i+1])
        if maxy != 0
            ylims!(ax, 0, 1.05 * maxy)
        end
    end
end

save(locplot * plot_name * "-timesteps_plot.png", fig)