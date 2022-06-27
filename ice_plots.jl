
"""
Author: Sergio Pérez Montero\n
Date: 22.06.2022\n

Aim: Some functions to make plots\n

"""

# Packages
using Pkg
Pkg.activate("ice_env")
using NetCDF
using CairoMakie
using ColorSchemes

# 1D variables

#-- plot_lines
@doc """
        Aim: Plots 1D data in different panels                                  
        Parameters
        ----------                                                              
        x              -> 1D vector                                             
        data           -> 1D to 3D array, [npanels, nseries, nx]            
        xname          -> string                                                
        yname          -> vector with strings, one for each panel               
        clrs           -> vector with colornames or colormaps, one for each series
        lnstls         -> style, one for each series
        lnswths        -> width, one for each series
        lbls           -> labeling, one for each series

        Optional
        ----------
        lyout          -> layout, (nrows, ncols)
        fgsz           -> figure size, tuple
        fntsz          -> font size
        ptsv           -> path/name to save plot
    """
function plot_lines(x, data, xname, yname, clrs, lnstls, lnswdths, lbls; lyout=(1, 1), fgsz=(800, 600), fntsz=nothing, ptsv="./plot_lines.png")
    # Check fntsz and theming
    isnothing(fntsz) && (fntsz = 0.05 * min(fgsz[1], fgsz[2]))
    fontsize_theme = Theme(font="Dejavu Sans", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Create figure and layout
    fig = Figure(resolution=fgsz)
    nrows, ncols = lyout
    if ndims(data) == 1
        npanels, nseries, nx = 1, 1, size(data)
    elseif ndims(data) == 2
        npanels, nseries, nx = 1, size(data)[1], size(data)[2]
    elseif ndims(data) == 3
        npanels, nseries, nx = size(data)
    end
    
    miny, maxy = minimum(data), maximum(data)

    # Plotting
    p = 1
    for i in 1:nrows, j in 1:ncols
        ax = Axis(fig[i, j], grid=true,
            xlabelsize=fntsz, ylabelsize=fntsz,
            xlabel=xname, ylabel=yname[p])
        for k in 1:nseries
            if ndims(data) == 1
                lines!(ax, x, data, color=clrs[k], linestyle=lnstls[k], linewidth=lnswdths[k], label=lbls[k])
            elseif ndims(data) == 2
                lines!(ax, x, data[k, :], color=clrs[k], linestyle=lnstls[k], linewidth=lnswdths[k], label=lbls[k])
            elseif ndims(data) == 3
                lines!(ax, x, data[p, k, :], color=clrs[k], linestyle=lnstls[k], linewidth=lnswdths[k], label=lbls[k])
            end
        end
        p = p + 1
        axislegend(ax)
        limits!(ax, x[1], x[end], miny, maxy) 

    end

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
function plot_maps(y, x, data, xname, yname, varname, lvls, lbls; clrmp=:thermal, lyout=(1, 1), fgsz=(800, 600), fntsz=nothing, ptsv="./plot_maps.png", hide_axis=false)
    # Check fntsz and theming
    isnothing(fntsz) && (fntsz = 0.05 * min(fgsz[1], fgsz[2]))
    fontsize_theme = Theme(font="Dejavu Sans", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Create figure and layout
    fig = Figure(resolution=fgsz)
    nrows, ncols = lyout
    if ndims(data) == 2
        npanels, ny, nx = 1, size(data)[1], size(data)[2]
    elseif ndims(data) == 3
        npanels, ny, nx = size(data)
    end

    # Plotting
    p = 1
    for i in 1:nrows, j in 1:ncols
        ax = Axis(fig[i, j], title=lbls[p],
            xlabelsize=fntsz, ylabelsize=fntsz,
            xlabel=xname, ylabel=yname)
        if ndims(data) == 2
            c = contourf!(ax, x, y, data, levels=lvls, colormap=clrmp)
        elseif ndims(data) == 3
            c = contourf!(ax, x, y, data[p, :, :], levels=lvls, colormap=clrmp)
        end
        limits!(ax, x[1], x[end], y[1], y[end])
        hide_axis && hidedecorations!(ax)
        p = p + 1
        (i * j == nrows * ncols) && (Colorbar(fig[:, end+1], c, height=Relative(1.5 / 3), label=varname, ticks=lvls))
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
        data           -> 3D array, [nvariables, ny, nx]            
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
function plot_multivar(data, xnames, ynames, varname, lvls; clrmp=:thermal, fgsz=(800, 600), fntsz=nothing, ptsv="./plot_multivar.png")
    # Check fntsz and theming
    isnothing(fntsz) && (fntsz = 0.05 * min(fgsz[1], fgsz[2]))
    fontsize_theme = Theme(font="Dejavu Sans", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Create figure and layout
    fig = Figure(resolution=fgsz)
    nrows, ncols = length(ynames), length(xnames)

    # Plotting
    p = 1
    for i in 1:nrows, j in 1:ncols
        if j == 1
            if i == nrows
                ax = Axis(fig[i, j], xlabelsize=fntsz, ylabelsize=fntsz, ylabel=ynames[i], xlabel=xnames[j])
            else
                ax = Axis(fig[i, j], xlabelsize=fntsz, ylabelsize=fntsz, ylabel=ynames[i])
            end
        else
            if i == nrows
                ax = Axis(fig[i, j], xlabelsize=fntsz, ylabelsize=fntsz, xlabel=xnames[j])
            else
                ax = Axis(fig[i, j], xlabelsize=fntsz, ylabelsize=fntsz)
            end
        end
        c = contourf!(ax, data[p, :, :], levels=lvls, colormap=clrmp)
        (i * j == nrows * ncols) && (Colorbar(fig[:, end+1], c, height=Relative(1.5 / 3), label=varname, ticks=lvls))
        hidespines!(ax, :t, :b, :r, :l)
        hidedecorations!(ax, label=false)
        p = p + 1
    end

    # Resizing
    resize_to_layout!(fig)

    # Saving
    save(ptsv, fig)
end

locdata = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d01_TFM/bmb-sliding-fmb_yelmo-v1.75/abum_pmp-m3q0.5f10.0/yelmo1D.nc"
locdata2 = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/d01_TFM/bmb-sliding-fmb_yelmo-v1.75/abum_pmp-m3q0.5f10.0/yelmo2D.nc"

d = ncread(locdata, "V_sle");
dd = zeros(2, 2, 501);
dd[1, 1, :] = d;
dd[1, 2, :] = d .- 15;
dd[2, 1, :] = d .+ 30;
dd[2, 2, :] = d .* 2;

d2 = ncread(locdata2, "H_ice")
x, y = ncread(locdata2, "xc"), ncread(locdata2, "yc")
d2 = d2[:, :, 1]
dd2 = zeros(6, 191, 191);
dd2[1, :, :], dd2[2, :, :] = d2, d2
plot_maps(y, x, dd2[:, :,:], "xc (km)", "yc (km)", "H_ice (m)", 0:500:4500, ["abum_pmp"], lyout=(1, 1), fgsz=(800, 600), hide_axis=true);
plot_multivar(dd2[:, :,:], ["1", "2", "3"], ["f2", "f3"], "H_ice (m)", 0:500:4500, fgsz=(800, 600));

plot_lines(1:501, dd, "t", ["V_sle (m SLE)", "v-5"],
    [:tomato, :cyan, :blue, :green],
    [:solid, :dash, :dotted, :solid],
    [5, 3],
    ["v", "v - 5"], lyout=(1, 2));