
"""
Author: Sergio Pérez Montero\n
Date: 27.06.2022\n

Aim: Easy and quick script to plot Yelmo (Robinson et al. 2020) results\n

"""

# Packages
using Pkg
Pkg.activate("ice_env")

# -- ice_tools
include("../src/ice_calcs.jl")
include("../src/ice_plots.jl")

# User parameters
locdata = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d03_LateralBC/dtt_02/"           # path to locate the data
locplot = "/home/sergio/entra/proyects/d03_LateralBC/plots/dtt_02/"                                         # path to save plots
expnames = ["dtt_02_akad0.1", "dtt_02_akad0.5", "dtt_02_akad0.8", "dtt_02_akad1.0",
            "dtt_02_akmd0.1", "dtt_02_akmd0.5", "dtt_02_akmd0.8", "dtt_02_akmd1.0",
            "dtt_02_akfd0.1", "dtt_02_akfd0.5", "dtt_02_akfd0.8", "dtt_02_akfd1.0"]       # vector with the names of the experiments
explabels = copy(expnames)     # vector with the names you want to show (nothing/label)
xpars, ypars = ["dtt = 0.1", "dtt = 0.5", "dtt = 0.8", "dtt = 1.0"],              # vector with x and y labels mor par vs par plot
                ["ALL", "MARINE", "FLOATING"]
plot_name = "dtt_02_abuk"

varnames = ["SLR"]      # vector with variables to plot, if ["all"] it will plot all the implemented variables 
units_time = "yrs"
times2D = "end" # "end"/index

colors, lstyles, lwidths = [repeat(["red"], 4); repeat(["blue"], 4); repeat(["green"], 4)], repeat([:dash, :solid, :dot, :dashdot], 3), repeat([3], 24)
layout_2D = (3, 4) # rows, cols

# Some implemented variables and dictionaries
yelmo1D_vars = ["V_sle", "A_ice", "dHicedt"]
yelmo2D_vars = ["H_ice", "H_grnd", "z_srf", "uxy_s", "pc_tau_max"]
composite_vars = ["SLR", "diffH_ice", "diffuxy_s"]

var1D_ylimits = Dict("V_sle" => [35, 57], "A_ice" => [7, 15], "dHicedt" => [-30, 30])
var1D_xlimits = [0, 500]    # set to [] if you don't want a predefined
var2D_levels = Dict("H_ice" => 0:100:4500, "H_grnd" => -4500:200:4500, "z_srf" => 0:100:4500, "uxy_s" => 0:10:1e4, "pc_tau_max" => 0:0.1:4)
composite_levels = Dict("SLR" => [0, 35], "diffH_ice" => -1000:100:1000, "diffuxy_s" => -1000:100:1000)
clrmps = Dict("H_ice" => :Blues, "H_grnd" => :curl, "z_srf" => :dense, "uxy_s" => :BuPu_7, "pc_tau_max" => :matter)
units = Dict("V_sle" => "m SLE", "A_ice" => "1e6 km^2", "dHicedt" => "m/a",
    "H_ice" => "m", "H_grnd" => "m", "z_srf" => "m", "uxy_s" => "m/a", "pc_tau_max" => "m/a")

#### SCRIPT ####
println("Plotting " * plot_name)

# Check if paths exist 
isdir(locplot) || throw(ErrorException(string(locplot, " does not exist")))
isdir(locdata) || throw(ErrorException(string(locdata, " does not exist")))

# Check if layout and number of experiments match
(layout_2D[1] * layout_2D[2] == length(expnames)) || throw(ErrorException("Layout and number of experiments don't match!"))

## Load
# Check variables to plot, then load them
(varnames[1] == "all") && (varnames = [yelmo1D_vars; yelmo2D_vars; composite_vars])
varlabels = copy(varnames)

yelmo_file = String[]
vars_1D, vars_2D, vars_C = String[], String[], String[]
labels_1D, labels_2D, labels_C = String[], String[], String[]
units_1D, units_2D = String[], String[]
for i in 1:length(varnames)
    if varnames[i] in yelmo1D_vars
        push!(vars_1D, varnames[i])
        push!(labels_1D, varlabels[i])
        push!(units_1D, units[varnames[i]])
    elseif varnames[i] in yelmo2D_vars
        push!(vars_2D, varnames[i])
        push!(labels_2D, varlabels[i])
        push!(units_2D, units[varnames[i]])
    elseif varnames[i] in composite_vars
        push!(vars_C, varnames[i])
        push!(labels_C, varlabels[i])
    else
        throw(ErrorException(string(varnames[i], " is not implemented yet")))
    end
end

# Check, load and plot
if ~isempty(vars_1D)    # Check
    for v in 1:length(vars_1D)
        time_data = []
        data_array_1D = []
        for j in 1:length(expnames)
            push!(data_array_1D, ncread(locdata * expnames[j] * "/yelmo1D.nc", vars_1D[v]))
            push!(time_data, ncread(locdata * expnames[j] * "/yelmo1D.nc", "time"))
        end
        xlab, ylab = "Time (" * units_time * ")", labels_1D[v] .* " (" .* units_1D[v] .* ")"
        plot_name_1D = locplot * vars_1D[v] * "_" * plot_name * ".png"
        plot_lines(time_data, data_array_1D, xlab, ylab, colors, lstyles, lwidths, explabels, ylimits=var1D_ylimits[vars_1D[v]], xlimits=var1D_xlimits,
                fntsz=nothing, ptsv=plot_name_1D)
        display(vars_1D[v] * " plotted")
    end
end
if ~isempty(vars_2D) # Check
    # Load
    yaxis, xaxis = ncread(locdata * expnames[1] * "/yelmo2D.nc", "yc"), ncread(locdata * expnames[1] * "/yelmo2D.nc", "xc")
    for v in 1:length(vars_2D)
        data_array_2D = zeros(length(expnames), length(yaxis), length(xaxis))
        contour_array_2D = zeros(length(expnames), length(yaxis), length(xaxis))
        time_labels = String[]
        for i in 1:length(expnames)
            if times2D == "end"
                data_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", vars_2D[v])[:, :, end]
                contour_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
                time_lbl = ncread(locdata * expnames[i] * "/yelmo2D.nc", "time")[end]
            else
                try
                    data_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", vars_2D[v])[:, :, times2D]
                    contour_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", "mask_bed")[:, :, times2D]
                    time_lbl = ncread(locdata * expnames[i] * "/yelmo2D.nc", "time")[times2D]
                catch
                    data_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", vars_2D[v])[:, :, end]
                    contour_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
                    time_lbl = ncread(locdata * expnames[i] * "/yelmo2D.nc", "time")[end]
                end
            end
            push!(time_labels, ", t = " * string(round(Int, time_lbl)) * "yrs")
        end
        # Now we modify title labels
        explabels2plot = explabels.*time_labels

        # Plot
        if var2D_levels[vars_2D[v]] == []
            n2approx = 500
            min_lvl, max_lvl = round(Int, minimum(data_array_2D) / n2approx) * n2approx, round(Int, maximum(data_array_2D) / n2approx) * n2approx
            distance, boundary = Int(max_lvl - min_lvl), max(abs(max_lvl), abs(min_lvl))
            step = round(Int, 0.1 * boundary / n2approx) * n2approx
            (min_lvl * max_lvl >= 0) ? (levels2D = min_lvl:step:max_lvl) : (levels2D = -boundary:step:boundary) # Check if symmetric
            (min_lvl * max_lvl >= 0) && (levels2D = round.(Int, levels2D))
        else
            levels2D = var2D_levels[vars_2D[v]]
        end

        (vars_2D[v] == "uxy_s") ? (ifuxy = true) : (ifuxy = false)
        xlab, ylab = "xc (km)", "yc (km)"
        figure_size = (800 * layout_2D[2], 750 * layout_2D[1])
        plot_name_2D = locplot * vars_2D[v] * "_" * plot_name * ".png"
        plot_maps(yaxis, xaxis, data_array_2D, xlab, ylab, labels_2D[v] * " (" * units_2D[v] * ")", levels2D, explabels2plot,
            log_scale=ifuxy, clrmp=clrmps[vars_2D[v]], lyout=layout_2D, fgsz=figure_size, fntsz=nothing, ptsv=plot_name_2D, hide_axis=true, cont=contour_array_2D, cont_lvls=[0, 3])
        display("$(vars_2D[v]) plotted")
    end
end
if ~isempty(vars_C)
    for v in vars_C
        if v == "SLR"
            include("../src/ice_SLR-plots.jl")
            calc_and_plot_SLR(locdata, expnames, explabels, colors, lstyles, lwidths, units_time, locplot, plot_name, ylimits=composite_levels["SLR"], xlimits=var1D_xlimits)
            display("SLR plotted")
        end
        if v == "diffH_ice"
            include("../src/ice_dif-plots.jl")
            calc_and_plot_diffH_ice(locdata, expnames, xpars, ypars, times2D, layout_2D, locplot, plot_name, lvls2plot=composite_levels["diffH_ice"])
            display("diffH_ice plotted")
        end
        if v == "diffuxy_s"
            include("../src/ice_dif-plots.jl")
            calc_and_plot_diffuxy_s(locdata, expnames, xpars, ypars, times2D, layout_2D, locplot, plot_name, lvls2plot=composite_levels["diffuxy_s"])
            display("diffuxy_s plotted")
        end
    end
end

