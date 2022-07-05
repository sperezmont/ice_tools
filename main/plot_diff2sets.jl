
"""
Author: Sergio Pérez Montero\n
Date: 04.07.2022\n

Aim: This script plots the difference between two sets of simulations (loading method = NetCDF)\n
     As a default case, it plots differences in H_ice for Yelmo\n

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
expnames = [["abuc_fbc0", "abum-mod_fbc0", "abum_fbc0", "abuk_fbc0"],
            ["abuc_fbco", "abum-mod_fbco", "abum_fbco", "abuk_fbco"]]       # vector with the names of the experiments divided [[set 1], [set 2]]

xpars, ypars = ["ABUC", "ABUM-mod", "ABUM", "ABUK"],              # vector with x and y labels for par vs par plot
               ["zero - ocean"]
plot_name = "taubc_01f_0-o"

#### SCRIPT ####
# Load
yaxis, xaxis = ncread(locdata*expnames[1][1]*"/yelmo2D.nc", "yc"), ncread(locdata*expnames[1][1]*"/yelmo2D.nc", "xc")
d_array = zeros(length(expnames), length(expnames[1]), length(yaxis), length(xaxis))
m_array = zeros(length(expnames), length(expnames[1]), length(yaxis), length(xaxis))
for i in 1:length(expnames), j in 1:length(expnames[1])
    if ~isfile(locdata*expnames[i][j]*"/yelmo_killed.nc")
        d_array[i, j, :, :] = ncread(locdata*expnames[i][j]*"/yelmo2D.nc", "H_ice")[:, :, end]
        m_array[i, j, :, :] = ncread(locdata*expnames[i][j]*"/yelmo2D.nc", "mask_bed")[:, :, end]
    else
        d_array[i, j, :, :] .= 0
        m_array[i, j, :, :] .= 0
    end
end
d_array[m_array .== 0] .= NaN

# Calculations
diff_array = d_array[1, :, :, :] .- d_array[2, :, :, :]
for i in 1:length(expnames), j in 1:length(expnames[1])
    all(y->y==0, d_array[i, j, :, :]) && (diff_array[j, :, :] .= 0)   # Just to deal with killed simulations
end
diff_array = replace!(diff_array, NaN=>0)

mask_array = m_array[1, :, :, :] + m_array[2, :, :, :]

# Plots
n2approx = 500
min_lvl, max_lvl = round(Int, minimum(diff_array)/n2approx)*n2approx, round(Int, maximum(diff_array)/n2approx)*n2approx
distance, boundary = Int(max_lvl - min_lvl), max(abs(max_lvl), abs(min_lvl))
step = max(1, round(Int, 0.05*boundary/n2approx))*n2approx
(min_lvl*max_lvl >= 0) ? (levels2D = min_lvl:step:max_lvl) : (levels2D = -boundary:step:boundary) # Check if symmetric
(min_lvl*max_lvl >= 0) && (levels2D = round.(Int, levels2D))

d_array = replace!(d_array, NaN=>0)
for i in 1:length(expnames), j in 1:length(expnames[1])
    all(y->y==0, d_array[i, j, :, :]) && (diff_array[j, :, :] .= -9999999999) # Just to deal with killed simulations  
    all(y->y==0, d_array[i, j, :, :]) && (mask_array[j, :, :] .= -9999999999)
end

diff_array[mask_array .== 0] .= -9999999999 # Just to deal with ocean

figure_size = (800*length(xpars), 750*length(ypars))  
plot_multivar(diff_array, xpars, ypars, "Difference H_grnd (m)", levels2D, clrmp=:balance, fgsz=figure_size, ptsv=locplot*"setDiffs_"*plot_name*".png", cont=mask_array, cont_lvls=[0])


