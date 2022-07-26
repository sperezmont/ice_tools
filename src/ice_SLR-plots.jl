
"""
Author: Sergio Pérez Montero\n
Date: 01.07.2022\n

Aim: This script calculates SLR for different experiments\n

"""
# Dependencies
include("../src/ice_calcs.jl")
include("../src/ice_plots.jl")

# calc_and_plot_SLR
function calc_and_plot_SLR(path2files, file_names, file_labels, cs, lss, lws, time_units, path2save, plot2save; ylimits=[], xlimits=[])

    # Load
    time_data = []
    slr_array = []

    for i in 1:length(file_names)
        # Load
        d_arrayi = ncread(path2files * file_names[i] * "/yelmo1D.nc", "V_sle")
        push!(time_data, ncread(path2files * file_names[i] * "/yelmo1D.nc", "time"))

        #Calculations
        slr_arrayi = SLR(d_arrayi)
        push!(slr_array, slr_arrayi)
    end

    # Plot
    xlab, ylab = "Time (" * time_units * ")", "SLR (m SLE)"
    plot2save = path2save * "SLR_" * plot2save * ".png"
    plot_lines(time_data, slr_array, xlab, ylab, cs, lss, lws, file_labels, ylimits=ylimits, xlimits=xlimits,
        fntsz=nothing, ptsv=plot2save)
end