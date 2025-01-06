using ClimaOcean
using Oceananigans
using SpeedyWeather
using Oceananigans.Units

include("method_extension.jl")

#####
##### Ocean model (idealized quarter degree model)
#####

Nx = 90
Ny = 40
Nz = 10

depth = 5000
φmax  = 75

arch    = Oceananigans.CPU()
r_faces = ClimaOcean.exponential_z_faces(; Nz, depth)

# A quarter degree ocean model (idealized)
grid = LatitudeLongitudeGrid(arch, 
                             size=(Nx, Ny, Nz), 
                             latitude=(-φmax, φmax), 
                             longitude=(0, 360),
                             halo = (6, 6, 2),
                             z = r_faces)

ocean = ClimaOcean.ocean_simulation(grid)

# Initial conditions 
# - parabolic temperature profile with a stratification
# - constant salinity

Tᵢ(λ, φ, z) = 30.0 * cosd(φ)^2 * (1 + z / depth)
Sᵢ(λ, φ, z) = 35.0 

Oceananigans.set!(ocean.model, T=Tᵢ, S=Sᵢ)

#####
##### Atmospheric model
#####

spectral_grid = SpectralGrid(trunc=121, nlayers=8, Grid=FullClenshawGrid)
model         = PrimitiveWetModel(spectral_grid)
atmosphere    = initialize!(model)

#=
We probably need to do SpeedyWeather.first_timesteps!(atmosphere) here
to initialize the leapfrog scheme correctly
=# 

#####
##### Coupled model
##### 

radiation   = Radiation() # Need to fix this to be consistent with SpeedyWeather? 
earth_model = OceanSeaIceModel(ocean; atmosphere, radiation)

# We eventually need to put some checks in `OceanSeaIceModel` to make sure that the atmosphere
# and the ocean have consistent time-step sizes
Δt = atmosphere.model.time_stepping.Δt_sec
earth = ClimaOcean.Simulation(earth_model; Δt)

#####
##### Progress function
#####

function progress(earth)

    # Some outputs from speedyweather?


    # Some outputs from oceananigans?

end

#####
##### Run the coupled model
##### 

time_step!(earth) # If no issue we can think about the details