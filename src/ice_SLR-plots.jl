
"""
Author: Sergio Pérez Montero\n
Date: 01.07.2022\n

Aim: This script calculates SLR for different experiments\n

"""
# Dependencies
include("../src/ice_calcs.jl")
include("../src/ice_plots.jl")

# calc_and_plot_SLR
function calc_and_plot_SLR(path2files, file_names, file_labels, cs, lss, lws, time_units, path2save, plot2save)

    # Load
    time_data = ncread(path2files * file_names[1] * "/yelmo1D.nc", "time")
    d_array = zeros(length(file_names), length(time_data))
    for i in 1:length(file_names)
        if ~isfile(path2files * file_names[i] * "/yelmo_killed.nc")
            d_array[i, :] = ncread(path2files * file_names[i] * "/yelmo1D.nc", "V_sle")
        else
            d_array[i, :] .= Inf
        end
    end

    #Calculations
    slr_array = SLR(d_array)
    slr_array = replace!(slr_array, NaN => Inf)

    # Plot
    xlab, ylab = "Time (" * time_units * ")", "SLR (m SLE)"
    plot2save = path2save * "SLR_" * plot2save * ".png"
    plot_lines(time_data, slr_array, xlab, ylab, cs, lss, lws, file_labels,
               fntsz=nothing, ptsv=plot2save)
end