
"""
Author: Sergio Pérez Montero\n
Date: 22.06.2022\n

Aim: Some functions to make calculations\n

"""

# Packages
#using Pkg
#Pkg.activate("ice_env")

# Functions

#-- SLE
@doc """ Transforms volume data to m SLE \n
        [data] = km^3 \n
        rhow = 1 Gt/m3 -> "disregarding the minor salinity/density effects of mixing fresh meltwater with seawater"
            More about: https://sealevel.info/conversion_factors.html
    """
function SLE(data; rhoi=0.9167, rhow=1.0, Aoc=3.618*10^8)
    SLE = rhoi/rhow .* 1e3 ./ Aoc .* data
    return SLE
end

#-- SLR
@doc """ Transforms volume data to m SLE and calculates the SLR for each time step \n
        [data] = km**3 \n
        rhow = 1 Gt/m3 -> "disregarding the minor salinity/density effects of mixing fresh meltwater with seawater"
            More about: https://sealevel.info/conversion_factors.html 
    """
function SLR(data; rhoi=0.9167, rhow=1, Aoc=3.618*10^8, v2s=false)
    (v2s) ? sle = SLE(data, rhoi=rhoi, rhow=rhow, Aoc=Aoc) : sle = copy(data)
    if ndims(data) == 1             # 1 time series    
        SLR = sle[1] .- sle
    elseif ndims(data) == 2         # >1 time series
        SLR = sle[:, 1] .- sle
    elseif ndims(data) == 3         # 1 2D field
        SLR = sle[1, :, :] .- sle
    elseif ndims(data) == 4         # >1 2D field
        SLR = sle[:, 1, :, :] .- sle
    end
    return SLR
end

#-- isfloating
@doc """ Calculates if ice is floating \n
        thick_data  -> ice thickness data
        slvl        -> sea level height
        bdh         -> bed height
    """
function is_floating(thick_data, slvl, bdh; rhoi=0.9167, rhosw=1.02)
    floating_mask = thick_data - rhosw/rhoi * max(slvl - bdh, 0)
    floating_mask[floating_mask <= 0] = true
    floating_mask[floating_mask > 0] = false
    return floating_mask
end

# Load 1D
@doc """
Loads variable var
"""
function load_var1D(var, path2files, file_names)
    # Load
    time_data = []
    var_array = []

    for i in 1:length(file_names)
        # Load
        d_arrayi = ncread(path2files * file_names[i] * "/yelmo1D.nc", var)
        push!(time_data, ncread(path2files * file_names[i] * "/yelmo1D.nc", "time"))
        push!(var_array, d_arrayi)
    end
    return var_array, time_data
end


# locdata = "/home/sergio/entra/ice_data/Antarctica/ANT-32KM/ANT-32KM_TOPO-BedMachine.nc"

# d = ncread(locdata,"H_ice");
# vol = sum(d/1000*(32)^2);
# sle = SLE(vol)

# size(d)

