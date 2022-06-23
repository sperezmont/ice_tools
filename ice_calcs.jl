
"""
Author: Sergio Pérez Montero\n
Date: 22.06.2022\n

Aim: Some functions to make calculations\n

"""

# Packages
using Pkg
Pkg.activate("ice_env")
using NetCDF

# Functions

function SLE(data, rhoi=0.9167, rhow=1.0, Aoc=3.618*10^8)
    """ Transforms volume data to m SLE \n
        [data] = km^3 \n
        rhow = 1 Gt/m3 -> "disregarding the minor salinity/density effects of mixing fresh meltwater with seawater"
            More about: https://sealevel.info/conversion_factors.html
    """
    SLE = rhoi/rhow * 1e3 / Aoc * data
    return SLE
end

function SLR(data, rhoi=0.9167, rhow=1, Aoc=3.618*10^8)
    """ Transforms volume data to m SLE and calculates the SLR for each time step \n
        [data] = km**3 \n
        rhow = 1 Gt/m3 -> "disregarding the minor salinity/density effects of mixing fresh meltwater with seawater"
            More about: https://sealevel.info/conversion_factors.html 
    """
    sle = SLE(data, rhoi, rhow, Aoc)
    SLR = sle - sle[0]
    return SLR
end

locdata = "/home/sergio/entra/ice_data/Antarctica/ANT-32KM/ANT-32KM_TOPO-BedMachine.nc"

d = ncread(locdata,"H_ice");
vol = sum(d/1000*(32)^2);
sle = SLE(vol)

