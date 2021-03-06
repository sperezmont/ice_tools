
"""
Author: Sergio Pérez Montero\n
Date: 30.06.2022\n

Aim: This script calculates and plots the differences between two fields\n

"""
# Dependencies
include("../src/ice_calcs.jl")
include("../src/ice_plots.jl")

# calc_and_plot_difH_ice
function calc_and_plot_diffH_ice(path2files, file_names, xvals, yvals, time_index, layout2plot, path2save, plot2save)

    # Load
    yaxis, xaxis = ncread(path2files * file_names[1] * "/yelmo2D.nc", "yc"), ncread(path2files * file_names[1] * "/yelmo2D.nc", "xc")
    d_array = zeros(length(file_names), 2, length(yaxis), length(xaxis))
    c_array = zeros(length(file_names), length(yaxis), length(xaxis))
    for i in 1:length(file_names)
        if time_index == "end"
            if ~isfile(locdata * file_names[i] * "/yelmo_killed.nc")
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, 1]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, end]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
            else
                replace!(d_array[i, :, :, :], 0 => Inf)
            end
        else
            if ~isfile(locdata * file_names[i] * "/yelmo_killed.nc")
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, 1]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, time_index]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, time_index]
            else
                d_array[i, :, :, :] = replace!(d_array[i, :, :, :], 0 => Inf)
            end
        end
    end

    # Calculations
    dif_array = d_array[:, 2, :, :] - d_array[:, 1, :, :]
    dif_array[c_array .== 0] .= 0

    # Plot
    min_lvl, max_lvl = round(Int, minimum(dif_array) / 100) * 100, round(Int, maximum(dif_array) / 100) * 100
    boundary = max(abs(max_lvl), abs(min_lvl))
    step = round(Int, 0.1 * boundary / 100) * 100
    levels2D = -boundary:step:boundary # assume symmetric
    levels2D = round.(Int, levels2D)
    #display(levels2D)
    #levels2D = -1500:100:1500
    dif_array[c_array .== 0] .= Inf

    figure_size = (750 * layout2plot[2], 700 * layout2plot[1])
    plot_name_2D = path2save * "diffH_ice_" * plot2save * ".png"
    plot_multivar(dif_array, xvals, yvals, "H_ice change (m)", levels2D, clrmp=:redsblues, fgsz=figure_size, ptsv=plot_name_2D, cont=c_array, cont_lvls=[0, 3])
end

# calc_and_plot_diffuxy_s
function calc_and_plot_diffuxy_s(path2files, file_names, xvals, yvals, time_index, layout2plot, path2save, plot2save)

    # Load
    yaxis, xaxis = ncread(path2files * file_names[1] * "/yelmo2D.nc", "yc"), ncread(path2files * file_names[1] * "/yelmo2D.nc", "xc")
    d_array = zeros(length(file_names), 2, length(yaxis), length(xaxis))
    c_array = zeros(length(file_names), length(yaxis), length(xaxis))
    for i in 1:length(file_names)
        if time_index == "end"
            if ~isfile(locdata * file_names[i] * "/yelmo_killed.nc")
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, 1]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, end]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
            else
                replace!(d_array[i, :, :, :], 0 => Inf)
            end
        else
            if ~isfile(locdata * file_names[i] * "/yelmo_killed.nc")
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, 1]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, time_index]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, time_index]
            else
                replace!(d_array[i, :, :, :], 0 => Inf)
            end
        end
    end

    # Calculations
    dif_array = d_array[:, 2, :, :] - d_array[:, 1, :, :]

    # Plot
    min_lvl, max_lvl = round(Int, minimum(dif_array) / 100) * 100, round(Int, maximum(dif_array) / 100) * 100
    boundary = max(abs(max_lvl), abs(min_lvl))
    step = round(Int, 0.1 * boundary / 100) * 100
    (min_lvl * max_lvl >= 0) ? (levels2D = min_lvl:step:max_lvl) : (levels2D = -boundary:step:boundary) # Check if symmetric
    (min_lvl * max_lvl >= 0) && (levels2D = round.(Int, levels2D))

    levels2D = -1500:100:1500 # [-1500, -500, -100, -50, -5, 5, 50, 100, 500, 1500]
    dif_array[c_array .== 0] .= Inf

    figure_size = (750 * layout2plot[2], 700 * layout2plot[1])
    plot_name_2D = path2save * "diffuxy_s_" * plot2save * ".png"
    plot_multivar(dif_array, xvals, yvals, "uxy_s change (m/a)", levels2D, clrmp=:bluesreds, fgsz=figure_size, ptsv=plot_name_2D, cont=c_array, cont_lvls=[0, 3])
end