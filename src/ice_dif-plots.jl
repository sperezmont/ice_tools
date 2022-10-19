
"""
Author: Sergio Pérez Montero\n
Date: 30.06.2022\n

Aim: This script calculates and plots the differences between two fields\n

"""
# Dependencies
include("../src/ice_calcs.jl")
include("../src/ice_plots.jl")

# calc_and_plot_difH_ice
function calc_and_plot_diffH_ice(path2files, file_names, xvals, yvals, time_index, time_units, layout2plot, path2save, plot2save; lvls2plot=[])

    # Load
    yaxis, xaxis = ncread(path2files * file_names[1] * "/yelmo2D.nc", "yc"), ncread(path2files * file_names[1] * "/yelmo2D.nc", "xc")
    d_array = zeros(length(file_names), 2, length(yaxis), length(xaxis))
    c_array = zeros(length(file_names), length(yaxis), length(xaxis))
    time_labels = String[]
    for i in 1:length(file_names)
        if time_index[2] == "end"
            d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, time_index[1]]
            d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, end]
            c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
            time_lbl = ncread(locdata * file_names[i] * "/yelmo2D.nc", "time")[end]
        else
            try
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, time_index[1]]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, time_index[2]]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, time_index[2]]
                time_lbl = ncread(locdata * file_names[i] * "/yelmo2D.nc", "time")[time_index[2]]
            catch
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, time_index[1]]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "H_ice")[:, :, end]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
                time_lbl = ncread(locdata * file_names[i] * "/yelmo2D.nc", "time")[end]
            end
        end
        if time_units == ""
            push!(time_labels, "t = " * string(round(Int, time_lbl)))
        else
            push!(time_labels, "t = " * string(round(Int, time_lbl)) * time_units)
        end
    end

    # Calculations
    dif_array = d_array[:, 2, :, :] - d_array[:, 1, :, :]
    dif_array[c_array.==0] .= 0

    # Plot
    if lvls2plot == []
        min_lvl, max_lvl = round(Int, minimum(dif_array) / 100) * 100, round(Int, maximum(dif_array) / 100) * 100
        boundary = max(abs(max_lvl), abs(min_lvl))
        step = round(Int, 0.1 * boundary / 100) * 100
        levels2D = -boundary:step:boundary # assume symmetric
        levels2D = round.(Int, levels2D)
    else
        levels2D = lvls2plot
    end
    dif_array[c_array.==0] .= Inf

    figure_size = (750 * layout2plot[2], 700 * layout2plot[1])
    plot_name_2D = path2save * "diffH_ice_" * plot2save * ".png"
    plot_multivar(dif_array, xvals, yvals, "H_ice change "*string(time_index[1])*"-"*String(time_index[2])*" (m)", levels2D, clrmp=cgrad([:red4, :firebrick, :red3, :indianred, :white, :skyblue2, :deepskyblue3, :royalblue3, :navyblue]), fgsz=figure_size, ptsv=plot_name_2D, cont=c_array, cont_lvls=[0, 3], lbls=time_labels)
end

# calc_and_plot_diffuxy_s
function calc_and_plot_diffuxy_s(path2files, file_names, xvals, yvals, time_index, time_units, layout2plot, path2save, plot2save; lvls2plot=[])

    # Load
    yaxis, xaxis = ncread(path2files * file_names[1] * "/yelmo2D.nc", "yc"), ncread(path2files * file_names[1] * "/yelmo2D.nc", "xc")
    d_array = zeros(length(file_names), 2, length(yaxis), length(xaxis))
    c_array = zeros(length(file_names), length(yaxis), length(xaxis))
    time_labels = String[]
    for i in 1:length(file_names)
        if time_index[2] == "end"
            d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, time_index[1]]
            d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, end]
            c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
            time_lbl = ncread(locdata * file_names[i] * "/yelmo2D.nc", "time")[end]
        else
            try
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, time_index[1]]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, time_index[2]]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, time_index[2]]
                time_lbl = ncread(locdata * file_names[i] * "/yelmo2D.nc", "time")[time_index]
            catch
                d_array[i, 1, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, time_index[1]]
                d_array[i, 2, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "uxy_s")[:, :, end]
                c_array[i, :, :] = ncread(locdata * file_names[i] * "/yelmo2D.nc", "mask_bed")[:, :, end]
                time_lbl = ncread(locdata * file_names[i] * "/yelmo2D.nc", "time")[end]
            end
        end
        if time_units == ""
            push!(time_labels, "t = " * string(round(Int, time_lbl)))
        else
            push!(time_labels, "t = " * string(round(Int, time_lbl)) * time_units)
        end
    end

    # Calculations
    dif_array = d_array[:, 2, :, :] - d_array[:, 1, :, :]

    # Plot
    if lvls2plot == []
        min_lvl, max_lvl = round(Int, minimum(dif_array) / 100) * 100, round(Int, maximum(dif_array) / 100) * 100
        boundary = max(abs(max_lvl), abs(min_lvl))
        step = round(Int, 0.1 * boundary / 100) * 100
        (min_lvl * max_lvl >= 0) ? (levels2D = min_lvl:step:max_lvl) : (levels2D = -boundary:step:boundary) # Check if symmetric
        (min_lvl * max_lvl >= 0) && (levels2D = round.(Int, levels2D))
    else
        levels2D = lvls2plot
    end
    dif_array[c_array.==0] .= Inf

    figure_size = (750 * layout2plot[2], 700 * layout2plot[1])
    plot_name_2D = path2save * "diffuxy_s_" * plot2save * ".png"

    plot_multivar(dif_array, xvals, yvals, "uxy_s change "*string(time_index[1])*"-"*String(time_index[2])*" (m/a)", levels2D, clrmp=cgrad([:navyblue, :royalblue3, :deepskyblue3, :skyblue2, :white, :indianred, :red3, :firebrick, :red4]), fgsz=figure_size, ptsv=plot_name_2D, cont=c_array, cont_lvls=[0, 3], lbls=time_labels)
end

# calc_plot_convergence
@doc """
Calculates and plots the convergence of the SLR results
"""
function calc_plot_convergence(path2files, file_names, file_labels, idx2c, clusters, members, cs, lss, lws, time_units, path2save, plot2save; ylimits=[], xlimits=[])
    # Load data
    v_vector, time_data = load_var1D("V_sle", path2files, file_names)

    # Compute clusters (assume equal time length)
    clust_lengths = [length(time_data[i]) for i in 1:length(file_names)]
    lengths_arr = reshape(clust_lengths, members, clusters)
    selected_length = maximum(lengths_arr) # store max length

    conv_array, ref_idx = [], []
    for i in 1:clusters
        reference_values = v_vector[members*(i-1)+idx2c]
        k = 0
        if length(reference_values) < selected_length
            while length(reference_values) < selected_length
                reference_values = v_vector[members*(i-1)+idx2c+k]
                k = k + 1
            end
        end
        
        #max_distancesi = []
        for j in 1:members
            push!(ref_idx, members*(i-1)+idx2c+k)
            if length(v_vector[members*(i-1)+j]) == selected_length
                distance_values = (v_vector[members*(i-1)+j] .- reference_values)
                push!(conv_array, distance_values)
            else
                conv_arrayij = Array{Float64}(undef, selected_length)
                conv_arrayij[:] .= -9999
                time_data[members*(i-1)+j] = Array{Float64}(undef, selected_length)
                push!(conv_array, conv_arrayij)
            end
        end

    end

    # Plot
    if time_units == ""
        xlab, ylab = "", "Distance (V_sle, m SLE)"
    else
        xlab, ylab = "Time (" * time_units * ")", "Distance (V_sle, m SLE)"
    end

    new_file_labels = []
    for i in 1:length(file_labels)
        push!(new_file_labels, file_labels[i]*" - "*file_labels[ref_idx[i]])
    end

    plot2save = path2save * "convergence_" * plot2save * ".png"
    plot_lines(time_data, conv_array, xlab, ylab, cs, lss, lws, new_file_labels, ylimits=ylimits, xlimits=xlimits,
        fntsz=nothing, ptsv=plot2save)

end

# calc_plot_correspondence
function calc_plot_correspondence(path2files, file_names, file_labels, idx2c, clusters, members, cs, lss, lws, time_units, path2save, plot2save; ylimits=[], xlimits=[])
    # Compute SLR
    exp_array, time_data = lc_SLR(path2files, file_names)

    # extract reference data
    ref_data = []
    for i in 1:clusters
        j = 1
        while j <= members
            push!(ref_data, exp_array[members*(i-1)+idx2c])
            j = j + 1
        end
    end

    # eliminate crashing simulations
    clust_lengths = [length(time_data[i]) for i in 1:length(file_names)]
    lengths_arr = reshape(clust_lengths, members, clusters)
    selected_length = maximum(lengths_arr) # store max length
    for i in 1:length(file_names)
        if length(exp_array[i]) != selected_length
            exp_array[i], ref_data[i] = Array{Float64}(undef, selected_length), Array{Float64}(undef, selected_length)
        end
    end

    # Plot
    xlab, ylab = "reference SLR (" * string(idx2c) * ", m SLE)", "SLR (m SLE)"
    plot2save = path2save * "correspondence_" * plot2save * ".png"
    plot_lines(ref_data, exp_array, xlab, ylab, cs, lss, lws, file_labels, ylimits=ylimits, xlimits=xlimits,
        fntsz=nothing, ptsv=plot2save)

end
# # old code for different time points
# function calc_plot_convergence(path2files, file_names, file_labels, idx2c, clusters, members, cs, lss, lws, time_units, path2save, plot2save; ylimits=[], xlimits=[])
#     # Compute SLR
#     slr_array, time_data = lc_SLR(path2files, file_names)

#     # Compute x-length of each cluster taken as a "chief" the first element of each cluster
#     clust_lengths = [length(time_data[i]) for i in 1:length(file_names)]
#     lengths_arr = reshape(clust_lengths, members, clusters)
#     selected_lengths = [findmax(lengths_arr[:, i]) for i in 1:clusters] # store max and position

#     # Interpolate
#     clust_arr, clust_times = [], []

#     k = 1
#     for i in 1:clusters
#         ref_times = copy(time_data[members*(i-1)+selected_lengths[i][2]])
#         ref_dt = ref_times[2] - ref_times[1]    # constant time-store step is presumed

#         arr_i, times_i = [], []
#         for j in 1:members
#             if lengths_arr[j, i] == selected_lengths[i][1]
#                 push!(arr_i, slr_array[k])
#                 push!(times_i, time_data[k])
#             else
#                 old_data, old_times = copy(slr_array[k]), copy(time_data[k])
#                 old_dt = old_times[2] - old_times[1]
#             end

#         end
#         k = k + 1
#     end

# end