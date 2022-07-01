# Packages
using Pkg
Pkg.activate("ice_env")

# -- ice_tools
include("../src/ice_calcs.jl")
include("../src/ice_plots.jl")

# User parameters
locplot = "/home/sergio/entra/proyects/d03_LateralBC/plots/taubc_01/"        # path to save plots
locdata = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d03_LateralBC/taubc_01/"        # path to locate the data
expnames = ["abuc_bc0", "abuk_bc0", "abum_bc0", "abum-mod_bc0",
            "abuc_bch", "abuk_bch", "abum_bch", "abum-mod_bch",
            "abuc_bco", "abuk_bco", "abum_bco", "abum-mod_bco"]       # vector with the names of the experiments
explabels = ["abuc_bc0", "abuk_bc0", "abum_bc0", "abum-mod_bc0",
            "abuc_bch", "abuk_bch", "abum_bch", "abum-mod_bch",
            "abuc_bco", "abuk_bco", "abum_bco", "abum-mod_bco"]      # vector with the names you want to show (nothing/label)
xpars, ypars = ["ABUC", "ABUK", "ABUM", "ABUM-mod"],
                ["zero", "ice", "ocean"]

varnames = ["difH_grnd"]       # vector with variables to plot "V_sle","H_ice", "H_grnd", "z_srf", "uxy_s", 
units = [""]     # if composite_var, ""         "m SLE", "m", "m", "m", "m/a", 
units_time = "yrs"
times2D = "end" # "end"/index

plot_name = "taubc_01"
colors, lstyles, lwidths = ["red", "blue", "green","orange","red", "blue", "green","orange","red", "blue", "green","orange"],     # ME QUEDO AQUI, NO FUNCIONA
                           [:solid,:solid,:solid,:solid,:dash,:dash,:dash,:dash,:dot,:dot,:dot,:dot,:dashdot,:dashdot,:dashdot,:dashdot],
                           [2, 3, 4, 6,2, 3, 4, 6,2, 3, 4, 6,2, 3, 4, 6]
layout_1D, layout_2D = (1, 1), (3, 4)

#### SCRIPT ####
yelmo1D_vars = ["V_sle"]
yelmo2D_vars = ["H_ice", "H_grnd", "z_srf", "uxy_s"]
composite_vars = ["difH_grnd"]

clrmps = Dict("H_ice"=>:Blues, "H_grnd"=>:curl, "z_srf"=>:dense, "uxy_s"=>:rainbow)

# Check if paths exist 
isdir(locplot) || throw(ErrorException(string(locplot," does not exist")))
isdir(locdata) || throw(ErrorException(string(locdata," does not exist")))

# Check if layout and number of experiments match
(layout_2D[1]*layout_2D[2] == length(expnames)) || throw(ErrorException("Layout and number of experiments don't match!"))

# Convert to LaTeX font
for i in 1:length(explabels)
    explabels[i] = replace(explabels[i], "_" => "~")
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
        push!(units_1D, units[i])
    elseif varnames[i] in yelmo2D_vars
        push!(vars_2D, varnames[i])
        push!(labels_2D, varlabels[i])
        push!(units_2D, units[i])
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
            continue
        end
    end

    # Plot
    xlab, ylab = "Time ("*units_time*")", labels_1D.*" (".*units_1D.*")"
    figure_size = (800*layout_1D[1], 600*layout_1D[2])  
    plot_name_1D = locplot*plot_name*"_p1D.png"
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
        if v == "difH_grnd"
            include("../src/ice_dif-plots.jl")
            calc_and_plot_difH_grnd(locdata, expnames, xpars, ypars, times2D, layout_2D,locplot, plot_name)
            display("difH_grnd plotted")
        end
    end
end

