
"""
Author: Sergio Pérez Montero\n
Date: 22.06.2022\n

Aim: This script configures the dependencies for using ice_tools\n

"""

# Environment generation
using Pkg
#Pkg.generate("ice_env")
Pkg.activate("ice_env")

# Adding dependencies ... 
## ice_calcs
Pkg.add("NetCDF")

## ice_plots
Pkg.add("CairoMakie")
Pkg.add("ColorSchemes")

# Check status
Pkg.status()

