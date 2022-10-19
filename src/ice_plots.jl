
"""
Author: Sergio Pérez Montero\n
Date: 22.06.2022\n

Aim: Some functions to make plots\n

"""

# Packages
#using Pkg
#Pkg.activate("ice_env")

using NetCDF
using CairoMakie
using ColorSchemes

# 1D variables

#-- plot_lines
@doc """
        Aim: Plots 1D data                                 
        Parameters
        ----------                                                              
        x              -> list (with each x series)                                            
        data           -> list, [nseries, nx]            
        xname          -> string                                                
        yname          -> vector with strings, one for each panel               
        clrs           -> vector with colornames or colormaps, one for each series
        lnstls         -> style, one for each series
        lnswths        -> width, one for each series
        lbls           -> labeling, one for each series

        Optional
        ----------
        fgsz           -> figure size, tuple
        fntsz          -> font size
        ptsv           -> path/name to save plot
    """
function plot_lines(x, data, xname, yname, clrs, lnstls, lnswdths, lbls; ylimits=[], xlimits=[], fgsz=(800, 600), fntsz=nothing, ptsv="./plot_lines.png")
    # Check fntsz and theming
    isnothing(fntsz) && (fntsz = 0.025 * sqrt(fgsz[1]^2 + fgsz[2]^2)) # isnothing(fntsz) && (fntsz = 0.05 * min(fgsz[1], fgsz[2]))
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Create figure and layout
    fig = Figure(resolution=fgsz)

    # Plotting
    ax = Axis(fig[1, 1], grid=true,
        xlabelsize=fntsz, ylabelsize=fntsz,
        xlabel=xname, ylabel=yname)
    update_theme!()
    min_vs, max_vs = [], []
    for k in 1:length(data)
        all(y -> y == Inf, data[k]) || lines!(ax, x[k], data[k], color=clrs[k], linestyle=lnstls[k], linewidth=lnswdths[k], label=lbls[k])
        all(y -> y == Inf, data[k]) || (push!(min_vs, minimum(data[k])))
        all(y -> y == Inf, data[k]) || (push!(max_vs, maximum(data[k])))
    end

    if length(lbls) != 1
        Legend(fig[1, 2], ax, framevisible=false, labelsize=0.8 * fntsz)
    end
    if ylimits == []
        new_min_vs, new_max_vs = [], []
        for i in 1:length(min_vs)
            (min_vs[i] != -9999) && push!(new_min_vs, min_vs[i])
            (max_vs[i] != -9999) && push!(new_max_vs, max_vs[i])
        end
        miny, maxy = minimum(new_min_vs), maximum(new_max_vs)
    else
        miny, maxy = ylimits
    end
    if xlimits == []
        minx, maxx = x[1], x[end]
    else
        minx, maxx = xlimits
    end
    limits!(ax, minx, maxx, miny, maxy)

    # Resizing
    resize_to_layout!(fig)

    # Saving
    save(ptsv, fig)
end

# 2D variables

#-- plot_maps
@doc """
        Aim: Plots 2D data in different panels                                  
        Parameters
        ----------                                                              
        y,x            -> 1D vector                                             
        data           -> 2D - 3D array, [npanels, ny, nx]            
        xname          -> string                                                
        yname          -> vector with strings, one for each panel  
        varname        -> variable name             
        lvls           -> vector or sequence of levels to plot
        lbls           -> labeling, one for each panel

        Optional
        ----------
        clrmp          -> colormap to be used
        lyout          -> layout, (nrows, ncols)
        fgsz           -> figure size, tuple
        fntsz          -> font size
        ptsv           -> path/name to save plot
        show_axis      -> bool, true=show, false=don't show
    """
function plot_maps(y, x, data, xname, yname, varname, lvls, lbls; log_scale=false, clrmp=:thermal, lyout=(1, 1), fgsz=(800, 600), fntsz=nothing, ptsv="./plot_maps.png", hide_axis=false, cont=[], cont_lvls=[])
    # Check fntsz and theming
    isnothing(fntsz) && (fntsz = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)) #isnothing(fntsz) && (fntsz = 0.02 * min(fgsz[1], fgsz[2]))
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Check log scale 
    if log_scale    # By: Jan Swierczek-Jereczek 
        data2plot = log10.(data .+ 1e-8) # to deal with zeros
        new_lvls = [10.0^i for i in [-2, -1, 0, 1, 2, 3, 4]]
        new_lvls = log10.(new_lvls)
        ticks_val = copy(new_lvls)
        ticks_str = string.([L"$0$", L"$0.1$", L"$1$",
            L"$10$", L"$100$", L"$1000$", L"$10000$"])   # Here we can even write by hand a nicer string that has exponents
    else
        data2plot = copy(data)
        new_lvls = copy(lvls)
        ticks_val = new_lvls[1:round(Int, length(new_lvls) / 5):end]
        ticks_str = string.(ticks_val)
    end

    # Create figure and layout
    fig = Figure(resolution=fgsz)
    nrows, ncols = lyout
    if ndims(data) == 2
        npanels, ny, nx = 1, size(data)[1], size(data)[2]
    elseif ndims(data) == 3
        npanels, ny, nx = size(data)
    end

    # Plotting
    p, maps = 1, []
    for i in 1:nrows, j in 1:ncols
        ax = Axis(fig[i, j], title=lbls[p], titlesize=0.7 * fntsz,
            xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz,
            xlabel=xname, ylabel=yname)
        update_theme!()
        limits!(ax, x[1], x[end], y[1], y[end])
        hide_axis && hidedecorations!(ax)
        if ndims(data) == 2
            if ~all(y -> y in [0, -8, -9999], data2plot)
                c = contourf!(ax, x, y, data2plot, levels=new_lvls, colormap=clrmp)
                (cont != []) && contour!(ax, x, y, cont[:, :], levels=cont_lvls, color="black", linewidth=1.5, linestyle="--")
                c.extendlow = :auto
                c.extendhigh = :auto
                push!(maps, c)
            end
        elseif ndims(data) == 3
            if ~all(y -> y in [0, -8, -9999], data2plot[p, :, :])
                c = contourf!(ax, x, y, data2plot[p, :, :], levels=new_lvls, colormap=clrmp)
                (cont != []) && contour!(ax, x, y, cont[p, :, :], levels=cont_lvls, color="black", linewidth=1.5, linestyle="--")
                c.extendlow = :auto
                c.extendhigh = :auto
                push!(maps, c)
            end
        end
        p = p + 1
        (i * j == nrows * ncols) && (Colorbar(fig[:, end+1], maps[end], height=Relative(2 / 3), width=30, label=varname, ticklabelsize=0.8 * fntsz, ticks=(ticks_val, ticks_str)))
    end

    # Resizing
    resize_to_layout!(fig)

    # Saving
    save(ptsv, fig)
end

#-- plot_multivar
@doc """
        Aim: Plots 2D data as a function of two variables/parameters                                
        Parameters
        ----------                                                              
        data           -> 3D array, [nexps, ny, nx]            
        xnames         -> string vector with values and names of each column (min --> max)                                                
        ynames         -> string vector with values and names of each row (min --> max)
        varname        -> variable name            
        lvls           -> vector or sequence of levels to plot

        Optional
        ----------
        clrmp          -> colormap to be used
        fgsz           -> figure size, tuple
        fntsz          -> font size
        ptsv           -> path/name to save plot
    """
function plot_multivar(data, xnames, ynames, varname, lvls; log_scale=false, clrmp=:thermal, fgsz=(800, 600), fntsz=nothing, ptsv="./plot_multivar.png", cont=[], cont_lvls=[], lbls=[])
    # Check fntsz and theming
    isnothing(fntsz) && (fntsz = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)) #isnothing(fntsz) && (fntsz = 0.05 * min(fgsz[1], fgsz[2]))
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    if lbls == []
        lbls = repeat([""], size(data)[1])
    end

    # Check log scale 
    if log_scale    # By: Jan Swierczek-Jereczek 
        data2plot = log10.(data .+ 1e-8) # to deal with zeros
        new_lvls = [10.0^i for i in [-2, -1, 0, 1, 2, 3, 4]]
        new_lvls = log10.(new_lvls)
        ticks_val = copy(new_lvls)
        ticks_str = string.([L"$0$", L"$0.1$", L"$1$",
            L"$10$", L"$100$", L"$1000$", L"$10000$"])   # Here we can even write by hand a nicer string that has exponents
    else
        data2plot = copy(data)
        new_lvls = copy(lvls)
        ticks_val = new_lvls[1:round(Int, length(new_lvls) / 5):end]
        ticks_str = string.(ticks_val)
    end

    # Create figure and layout
    fig = Figure(resolution=fgsz)
    nrows, ncols = length(ynames), length(xnames)

    # Plotting
    p, maps = 1, []
    for i in 1:nrows, j in 1:ncols
        if j == 1
            if i == nrows
                ax = Axis(fig[i, j], title=lbls[p], titlesize=0.7 * fntsz, xlabelsize=fntsz, ylabelsize=fntsz, ylabel=ynames[i], xlabel=xnames[j])
                update_theme!()
            else
                ax = Axis(fig[i, j], title=lbls[p], titlesize=0.7 * fntsz, xlabelsize=fntsz, ylabelsize=fntsz, ylabel=ynames[i])
                update_theme!()
            end
        else
            if i == nrows
                ax = Axis(fig[i, j], title=lbls[p], titlesize=0.7 * fntsz, xlabelsize=fntsz, ylabelsize=fntsz, xlabel=xnames[j])
                update_theme!()
            else
                ax = Axis(fig[i, j], title=lbls[p], titlesize=0.7 * fntsz, xlabelsize=fntsz, ylabelsize=fntsz)
                update_theme!()
            end
        end
        hidespines!(ax, :t, :b, :r, :l)
        hidedecorations!(ax, label=false)

        if ~all(y -> y in [0, -8, -9999], data2plot)
            try
                c = contourf!(ax, data2plot[p, :, :], levels=new_lvls, colormap=clrmp)
                c.extendlow = :auto
                c.extendhigh = :auto
                push!(maps, c)
                (cont != []) && contour!(ax, cont[p, :, :], levels=cont_lvls, color="black", linewidth=2)
            catch
                (cont != []) && contour!(ax, cont[p, :, :], levels=cont_lvls, color="black", linewidth=2)
                p = p + 1
                (i * j == nrows * ncols) && (Colorbar(fig[:, end+1], maps[end], height=Relative(2 / 3), width=30, label=varname, ticklabelsize=0.8 * fntsz, ticks=(ticks_val, ticks_str)))
                continue
            end
        end

        (i * j == nrows * ncols) && (Colorbar(fig[:, end+1], maps[end], height=Relative(2 / 3), width=30, label=varname, ticklabelsize=0.8 * fntsz, ticks=(ticks_val, ticks_str)))
        p = p + 1
    end

    # Resizing
    resize_to_layout!(fig)

    # Saving
    save(ptsv, fig)
end
# locdata = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d01_TFM/bmb-sliding-fmb_yelmo-v1.75/abum_pmp-m3q0.5f10.0/yelmo1D.nc"
# locdata2 = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d01_TFM/bmb-sliding-fmb_yelmo-v1.75/abum_pmp-m3q0.5f10.0/yelmo2D.nc"

# d = ncread(locdata, "V_sle");
# dd = zeros(2, 2, 501);
# dd[1, 1, :] = d;
# dd[1, 2, :] = d .- 15;
# dd[2, 1, :] = d .+ 30;
# dd[2, 2, :] = d .* 2;

# d2 = ncread(locdata2, "H_ice")
# x, y = ncread(locdata2, "xc"), ncread(locdata2, "yc")
# d2 = d2[:, :, 1]
# dd2 = zeros(6, 191, 191);
# dd2[1, :, :], dd2[2, :, :] = d2, d2
# plot_maps(y, x, dd2[:, :,:], "xc (km)", "yc (km)", "H_ice (m)", 0:500:4500, ["abum_pmp"], lyout=(1, 1), fgsz=(800, 600), hide_axis=true);
# plot_multivar(dd2[:, :,:], ["1", "2", "3"], ["f2", "f3"], "H_ice (m)", 0:500:4500, fgsz=(800, 600));

# plot_lines(1:501, dd, "t", ["V_sle (m SLE)", "v-5"],
#     [:tomato, :cyan, :blue, :green],
#     [:solid, :dash, :dotted, :solid],
#     [5, 3],
#     ["v", "v - 5"], lyout=(1, 2));