
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
locplot = "/home/sergio/entra/proyects/d03_LateralBC/plots/taubc_01/"        # path to save plots
locdata = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d03_LateralBC/taubc_01/"        # path to locate the data
expnames = ["abuc_fbc0", "abum-mod_fbc0", "abum_fbc0", "abuk_fbc0",
            "abuc_fbch", "abum-mod_fbch", "abum_fbch", "abuk_fbch",
            "abuc_fbco", "abum-mod_fbco", "abum_fbco", "abuk_fbco",
            "abuc_fbcs", "abum-mod_fbcs", "abum_fbcs", "abuk_fbcs"]       # vector with the names of the experiments
explabels = ["abuc f0", "abum-mod f0", "abum f0", "abuk f0",
             "abuc fh", "abum-mod fh", "abum fh", "abuk fh",
             "abuc fo", "abum-mod fo", "abum fo", "abuk fo",
             "abuc fs", "abum-mod fs", "abum fs", "abuk fs"]      # vector with the names you want to show (nothing/label)
xpars, ypars = ["ABUC", "ABUM-mod", "ABUM", "ABUK"],              # vector with x and y labels for par vs par plot
               ["zero", "ice", "ocean", "ice+ocean"]
                

varnames = ["SLR", "diffH_grnd", "diffuxy_s"]       # vector with variables to plot 
units_time = "yrs"
times2D = "end" # "end"/index

plot_name = "taubc_01f"
colors, lstyles, lwidths = repeat(["red", "blue", "green", "purple"], 4),    
                           [repeat([:solid], 4); repeat([:dash], 4); repeat([:dot], 4); repeat([:dashdot], 4)],
                           repeat([3], 16)
layout_1D, layout_2D = (1, 1), (4, 4)

#### SCRIPT ####
yelmo1D_vars = ["V_sle"]
yelmo2D_vars = ["H_ice", "H_grnd", "z_srf", "uxy_s"]
composite_vars = ["SLR", "diffH_grnd", "diffuxy_s"]

clrmps = Dict("H_ice"=>:Blues, "H_grnd"=>:curl, "z_srf"=>:dense, "uxy_s"=>:rainbow)
units = Dict("H_ice"=>"m", "H_grnd"=>"m", "z_srf"=>"m", "uxy_s"=>"m/a")

# Check if paths exist 
isdir(locplot) || throw(ErrorException(string(locplot," does not exist")))
isdir(locdata) || throw(ErrorException(string(locdata," does not exist")))

# Check if layout and number of experiments match
(layout_2D[1]*layout_2D[2] == length(expnames)) || throw(ErrorException("Layout and number of experiments don't match!"))

# Convert to LaTeX font
for i in 1:length(explabels)
    explabels[i] = replace(explabels[i], "_" => "~")
    explabels[i] = replace(explabels[i], " " => "~")
end
explabels = [L"%$(i)" for i in explabels]

varlabels = copy(varnames)
#for i in 1:length(varlabels)
    #varlabels[i] = replace(varlabels[i], "_" => "")
    #varlabels[i] = replace(varlabels[i], "_" => "_{")
    #varlabels[i] = varlabels[i]*"}"
#end
#varlabels = [L"%$(i)" for i in varlabels]
#units = [L"%$(i)" for i in units]

## Load
# Check variables to plot, then load them
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
        throw(ErrorException(string(varnames[i]," is not implemented yet")))
    end
end

# Check, load and plot
if ~isempty(vars_1D)    # Check
    # Load
    time_data = ncread(locdata*expnames[1]*"/yelmo1D.nc", "time")  
    data_array_1D = zeros(length(vars_1D), length(expnames), length(time_data))
    for i in 1:length(vars_1D), j in 1:length(expnames)
        if ~isfile(locdata*expnames[j]*"/yelmo_killed.nc")
            data_array_1D[i, j, :] = ncread(locdata*expnames[j]*"/yelmo1D.nc", vars_1D[i])
        else
            data_array_1D[i, j, :] .= Inf
        end
    end

    # Plot
    xlab, ylab = "Time ("*units_time*")", labels_1D.*" (".*units_1D.*")"
    figure_size = (800*layout_1D[1], 600*layout_1D[2])  
    plot_name_1D = locplot*"1D_"*plot_name*".png"
    plot_lines(time_data, data_array_1D, xlab, ylab, colors, lstyles, lwidths, explabels,
                         lyout=layout_1D, fgsz=figure_size, fntsz=nothing, ptsv=plot_name_1D)
    display("1D variables plotted")
end
if ~isempty(vars_2D) # Check
    # Load
    yaxis, xaxis = ncread(locdata*expnames[1]*"/yelmo2D.nc", "yc"), ncread(locdata*expnames[1]*"/yelmo2D.nc", "xc")
    for v in 1:length(vars_2D)
        data_array_2D = zeros(length(expnames), length(yaxis), length(xaxis))
        for i in 1:length(expnames)
            if times2D == "end"
                if ~isfile(locdata*expnames[i]*"/yelmo_killed.nc")
                    data_array_2D[i, :, :] = ncread(locdata*expnames[i]*"/yelmo2D.nc", vars_2D[v])[:, :, end]
                else
                    replace!(data_array_2D[i, :, :], 0=>Inf)
                end
            else
                if ~isfile(locdata*expnames[i]*"/yelmo_killed.nc")
                    data_array_2D[i, :, :] = ncread(locdata*expnames[i]*"/yelmo2D.nc", vars_2D[v])[:, :, times2D]
                else
                    replace!(data_array_2D[i, :, :], 0=>Inf)
                end
            end
        end

        # Plot
        n2approx = 500
        min_lvl, max_lvl = round(Int, minimum(data_array_2D)/n2approx)*n2approx, round(Int, maximum(data_array_2D)/n2approx)*n2approx
        distance, boundary = Int(max_lvl - min_lvl), max(abs(max_lvl), abs(min_lvl))
        step = round(Int, 0.1*boundary/n2approx)*n2approx
        (min_lvl*max_lvl >= 0) ? (levels2D = min_lvl:step:max_lvl) : (levels2D = -boundary:step:boundary) # Check if symmetric
        (min_lvl*max_lvl >= 0) && (levels2D = round.(Int, levels2D))

        (vars_2D[v] == "uxy_s") ? (ifuxy = true) : (ifuxy = false)
        xlab, ylab = "xc (km)", "yc (km)"
        figure_size = (800*layout_2D[2], 750*layout_2D[1])  
        plot_name_2D = locplot*vars_2D[v]*"_"*plot_name*".png"
        plot_maps(yaxis, xaxis, data_array_2D, xlab, ylab, labels_2D[v]*" ("*units_2D[v]*")", levels2D, explabels, 
                  log_scale=ifuxy, clrmp=clrmps[vars_2D[v]], lyout=layout_2D, fgsz=figure_size, fntsz=nothing, ptsv=plot_name_2D, hide_axis=true)
        display("$(vars_2D[v]) plotted")
    end
end
if ~isempty(vars_C)
    for v in vars_C
        if v == "SLR"
            include("../src/ice_SLR-plots.jl")
            calc_and_plot_SLR(locdata, expnames, explabels, colors, lstyles, lwidths, units_time, layout_1D, locplot, plot_name)
            display("SLR plotted")
        end
        if v == "diffH_grnd"
            include("../src/ice_dif-plots.jl")
            calc_and_plot_diffH_grnd(locdata, expnames, xpars, ypars, times2D, layout_2D,locplot, plot_name)
            display("diffH_grnd plotted")
        end
        if v == "diffuxy_s"
            include("../src/ice_dif-plots.jl")
            calc_and_plot_diffuxy_s(locdata, expnames, xpars, ypars, times2D, layout_2D,locplot, plot_name)
            display("diffuxy_s plotted")
        end
    end
end

