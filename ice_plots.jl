
"""
Author: Sergio Pérez Montero\n
Date: 22.06.2022\n

Aim: Some functions to make plots\n

"""

# Packages
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
        lyout          -> layout, (nrows, ncols)
        fgsz           -> figure size, tuple
        fntsz          -> font size
        ptsv           -> path/name to save plot
    """
function plot_lines(x, data, xname, yname, clrs, lnstls, lnswdths, lbls; lyout=(1,1), fgsz=(800, 600), fntsz=48, ptsv="./plot_lines.png")
    fontsize_theme = Theme(font = "Dejavu Sans", fontsize = fntsz)
    set_theme!(fontsize_theme)
    fig = Figure(resolution = fgsz)
    nrows, ncols = lyout
    if ndims(data) == 1
        npanels, nseries, nx = 1, 1, size(data)
    elseif ndims(data) == 2
        npanels, nseries, nx = 1, size(data)[1], size(data)[2]
    elseif ndims(data) == 3
        npanels, nseries, nx = size(data)
    end

    p = 1
    for i in 1:nrows, j in 1:ncols
        ax = Axis(fig[i, j],
                xlabelsize=fntsz, ylabelsize=fntsz,
                xlabel = xname, ylabel = yname[p]) 
        for k in 1:nseries
            if ndims(data) == 1
                lines!(ax, x, data, color = clrs[k], linestyle = lnstls[k], linewidth = lnswdths[k], label = lbls[k])
            elseif ndims(data) == 2 
                lines!(ax, x, data[k, :], color = clrs[k], linestyle = lnstls[k], linewidth = lnswdths[k], label = lbls[k])
            elseif ndims(data) == 3
                lines!(ax, x, data[p, k, :], color = clrs[k], linestyle = lnstls[k], linewidth = lnswdths[k], label = lbls[k])
            end
        end
        p = p + 1
        axislegend(ax)
        
    end
    save(ptsv, fig)
end

locdata = "/home/sergio/entra/models/yelmo_vers/v1.75/yelmox/output/ismip6/bmb-sliding-fmb_yelmo-v1.75/abum_pmp-m3q0.5f10.0/yelmo1D.nc"

d = ncread(locdata,"V_sle");
d2 = zeros(2, 2, 501);
d2[1, 1, :] = d;
d2[1, 2, :] = d .- 15;
d2[2, 1, :] = d .+ 30;
d2[2, 2, :] = d .*2;

plot_lines(1:501, d2, "t", ["V_sle (m SLE)", "v-5"], [:tomato, :cyan], [:solid, :dash], [5, 3], ["v", "v - 5"], lyout=(2,1))