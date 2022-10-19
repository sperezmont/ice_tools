
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
locdata = "/home/sergio/entra/models/yelmo_vers/v1.753_sergio-test/yelmox/output/ismip6/d03_LateralBC/ismip6_dtt_01/"           # path to locate the data
locplot = "/home/sergio/entra/proyects/d03_LateralBC/plots/ismip6_dtt_01/"                                         # path to save plots
expnames = ["ctrl_m", "ctrl_md0.1", "ctrl_md0.5", "ctrl_md1.0",
            #"exp05_m", "exp05_md0.1", "exp05_md0.5", "exp05_md1.0",
            #"exp09_m", "exp09_md0.1", "exp09_md0.5", "exp09_md1.0",
            #"exp10_m", "exp10_md0.1", "exp10_md0.5", "exp10_md1.0",
            ]#"exp13_m", "exp13_md0.1", "exp13_md0.5", "exp13_md1.0"]       # vector with the names of the experiments
explabels = copy(expnames)     # vector with the names you want to show (nothing/label)
xpars, ypars = ["m", "md0.1", "md0.5", "md1.0"],              # vector with x and y labels mor par vs par plot
                ["ctrl", "exp05", "exp09", "exp10", "exp13"]
plot_name = "v1.753_sergio-test_ismip6_dtt_01-ctrl_m"

varnames = ["convergence"]      # vector with variables to plot, if ["all"] it will plot all the implemented variables 
units_time = ""
times2D = [201-190, "end"] # [initial index, "end"/index] the first element is to compare plots

colors, lstyles, lwidths = [:black, :blue, :green, :red],repeat([:solid],5), repeat([2], 25)
layout_2D = (1, 4) # rows, cols

# Some implemented variables and dictionaries
yelmo1D_vars = ["V_sle", "A_ice", "dHicedt"]
yelmo2D_vars = ["H_ice", "H_grnd", "z_srf", "uxy_s", "pc_tau_max"]
composite_vars = ["SLR", "diffH_ice", "diffuxy_s", "convergence", "correspondence"]

conv_idx2cmpare = 1 # default index to compare with in convergence calculations

var1D_ylimits = Dict("V_sle" => [], "A_ice" => [7, 15], "dHicedt" => [-30, 30])
var1D_xlimits = [1900, 2500] #[0, 30000]    # set to [] if you don't want a predefined
var2D_levels = Dict("H_ice" => 0:100:4500, "H_grnd" => -4500:200:4500, "z_srf" => 0:100:4500, "uxy_s" => 0:10:1e4, "pc_tau_max" => 0:0.1:4)
composite_levels = Dict("SLR" => [-0.5, 2.5], "diffH_ice" => -1000:100:1000, "diffuxy_s" => -500:20:500, "convergence" => [-0.1, 0.4], "correspondence" => [-0.5, 1.5])
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
(varnames[1] == "composite") && (varnames = composite_vars)

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
        if units_time == ""
            xlab, ylab = "", labels_1D[v] .* " (" .* units_1D[v] .* ")"
        else
            xlab, ylab = "Time (" * units_time * ")", labels_1D[v] .* " (" .* units_1D[v] .* ")"
        end
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
            if times2D[2] == "end"
                data_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", vars_2D[v])[:, :, end]
                contour_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
                time_lbl = ncread(locdata * expnames[i] * "/yelmo2D.nc", "time")[end]
            else
                try
                    data_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", vars_2D[v])[:, :, times2D[2]]
                    contour_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", "mask_bed")[:, :, times2D[2]]
                    time_lbl = ncread(locdata * expnames[i] * "/yelmo2D.nc", "time")[times2D[2]]
                catch
                    data_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", vars_2D[v])[:, :, end]
                    contour_array_2D[i, :, :] = ncread(locdata * expnames[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
                    time_lbl = ncread(locdata * expnames[i] * "/yelmo2D.nc", "time")[end]
                end
            end
            if units_time == ""
                push!(time_labels, ", t = " * string(round(Int, time_lbl)))
            else
                push!(time_labels, ", t = " * string(round(Int, time_lbl)) * units_time)
            end
        end
        # Now we modify title labels
        explabels2plot = explabels .* time_labels

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
            plot_SLR(locdata, expnames, explabels, colors, lstyles, lwidths, units_time, locplot, plot_name, ylimits=composite_levels["SLR"], xlimits=var1D_xlimits)
            display("SLR plotted")
        end
        if v == "diffH_ice"
            include("../src/ice_dif-plots.jl")
            calc_and_plot_diffH_ice(locdata, expnames, xpars, ypars, times2D, units_time, layout_2D, locplot, plot_name, lvls2plot=composite_levels["diffH_ice"])
            display("diffH_ice plotted")
        end
        if v == "diffuxy_s"
            include("../src/ice_dif-plots.jl")
            calc_and_plot_diffuxy_s(locdata, expnames, xpars, ypars, times2D, units_time, layout_2D, locplot, plot_name, lvls2plot=composite_levels["diffuxy_s"])
            display("diffuxy_s plotted")
        end
        if v == "convergence"   # this part takes into account the number of 2D rows as number of clusters and ncols as mebers of each cluster
            include("../src/ice_dif-plots.jl")
            nclust, nmemb = layout_2D
            calc_plot_convergence(locdata, expnames, explabels, conv_idx2cmpare, nclust, nmemb, colors, lstyles, lwidths, units_time, locplot, plot_name, ylimits=composite_levels["convergence"], xlimits=var1D_xlimits)
            display("convergence plotted")
        end
        if v == "correspondence"   # this part takes into account the number of 2D rows as number of clusters and ncols as mebers of each cluster
            include("../src/ice_dif-plots.jl")
            nclust, nmemb = layout_2D
            calc_plot_correspondence(locdata, expnames, explabels, conv_idx2cmpare, nclust, nmemb, colors, lstyles, lwidths, units_time, locplot, plot_name, ylimits=composite_levels["correspondence"], xlimits=composite_levels["correspondence"])
            display("correspondence plotted")
        end
    end
end

