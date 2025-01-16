# This file contains all the extensions needed from SpeedyWeather to run the coupled model.

using ClimaOcean
using SpeedyWeather

using ClimaOcean.OceanSeaIceModels.Atmospheres: PrescribedAtmosphereThermodynamicsParameters

import Oceananigans: time_step!
import Oceananigans.Models: update_model_field_time_series!

#####
##### Extending the time_step! function and the `update_model_field_time_series!` function
#####

update_model_field_time_series!(::SpeedyWeather.Simulation, time) = nothing

time_step!(atmos::SpeedyWeather.Simulation) = SpeedyWeather.timestep!(atmos)

####
#### Extending inputs to flux computation
####

include("speedy_weather_parameters.jl")

# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.Atmospheres: thermodynamics_parameters, 
                                                 boundary_layer_height, 
                                                 surface_layer_height


# This should be the height of the surface layer in the atmospheric model
# We use a constant surface_layer_height even if a Speedy model has σ coordinates
function surface_layer_height(s::SpeedyWeather.Simulation) 
    T = s.model.atmosphere.temp_ref
    g = s.model.planet.gravity
    Φ = s.model.geopotential.Δp_geopot_full

    return Φ[end] * T / g
end

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::SpeedyWeather.Simulation) = 600

using SpeedyWeather: EarthAtmosphere

Base.eltype(::EarthAtmosphere{FT}) where FT = FT

# This is a _hack_!! The parameters should be consistent with what is specified in SpeedyWeather
thermodynamics_parameters(atmos::SpeedyWeather.Simulation) = 
    PrescribedAtmosphereThermodynamicsParameters(eltype(atmos.model.atmosphere))

#####
##### Extensions for interpolation between the ocean/sea-ice model and the atmospheric model
#####

import ClimaOcean.OceanSeaIceModels.Atmospheres: 
                    regrid_fluxes_to_atmospheric_model!, 
                    interpolate_atmospheric_state!

# Interpolate the atmospheric surface fields to the ocean/sea-ice model grid
function interpolate_atmospheric_state!(surface_atmosphere_state, 
                                        interpolated_prescribed_freshwater_flux, 
                                        atmos::SpeedyWeather.Simulation, 
                                        grid, clock)

    # Get the atmospheric state on the ocean grid
    ua = on_architecture(CPU(), surface_atmosphere_state.u)
    va = on_architecture(CPU(), surface_atmosphere_state.v)
    Ta = on_architecture(CPU(), surface_atmosphere_state.T)
    qa = on_architecture(CPU(), surface_atmosphere_state.q)
    pa = on_architecture(CPU(), surface_atmosphere_state.p)
    Qs = on_architecture(CPU(), surface_atmosphere_state.Qs)
    Qℓ = on_architecture(CPU(), surface_atmosphere_state.Qℓ)

    λ,  φ,  _ = Oceananigans.Grids.nodes(grid, Center(), Center(), Center(), with_halos=true) 

    λ = Array(vec(on_architecture(CPU(), λ)))
    φ = Array(vec(on_architecture(CPU(), φ)))

    spectral_grid = atmos.model.spectral_grid
    interpolator = RingGrids.AnvilInterpolator(Float32, spectral_grid.Grid, spectral_grid.nlat_half, length(λ))
    RingGrids.update_locator!(interpolator, λ, φ)

    RingGrids.interpolate!(vec(view(ua, :, :, 1)), atmos.diagnostic_variables.grid.u_grid[:, end],            interpolator)
    RingGrids.interpolate!(vec(view(va, :, :, 1)), atmos.diagnostic_variables.grid.v_grid[:, end],            interpolator)
    RingGrids.interpolate!(vec(view(Ta, :, :, 1)), atmos.diagnostic_variables.grid.temp_grid[:, end],         interpolator)
    RingGrids.interpolate!(vec(view(qa, :, :, 1)), atmos.diagnostic_variables.grid.humid_grid[:, end],        interpolator)
    RingGrids.interpolate!(vec(view(pa, :, :, 1)), exp.(atmos.diagnostic_variables.grid.pres_grid[:, end]),   interpolator)
    RingGrids.interpolate!(vec(view(Qs, :, :, 1)), atmos.diagnostic_variables.physics.surface_shortwave_down, interpolator)
    RingGrids.interpolate!(vec(view(Qℓ, :, :, 1)), atmos.diagnostic_variables.physics.surface_longwave_down,  interpolator)

    surface_atmosphere_state.u  .= ua 
    surface_atmosphere_state.v  .= va 
    surface_atmosphere_state.T  .= Ta 
    surface_atmosphere_state.q  .= qa 
    surface_atmosphere_state.p  .= pa 
    surface_atmosphere_state.Qs .= Qs
    surface_atmosphere_state.Qℓ .= Qℓ

    return nothing
end

using OrthogonalSphericalShellGrids: InterpolationWeights

arch       = Oceananigans.architecture(grid)
tmp_grid   = LatitudeLongitudeGrid(arch; size=(360, 179, 1), latitude=(-89.5, 89.5), longitude=(0, 360), z = (0, 1))
tmp_field  = Field{Center, Center, Nothing}(tmp_grid)
tmp_ongrid = Field{Center, Center, Nothing}(grid)

weights = InterpolationWeights(tmp_field, tmp_ongrid)

using Oceananigans.Fields: interpolate!

# Regrid the fluxes from the ocean/sea-ice grid to the atmospheric model grid
# This is a _hack_!! (We are performing a double interpolation)
function regrid_fluxes_to_atmospheric_model!(atmos::SpeedyWeather.Simulation, similarity_theory_fields)
    
    spectral_grid = atmos.model.spectral_grid

    Qs = similarity_theory_fields.sensible_heat
    Mv = similarity_theory_fields.water_vapor

    interpolator = RingGrids.AnvilInterpolator(Float32, FullClenshawGrid, 90, spectral_grid.npoints)
    londs, latds = RingGrids.get_londlatds(spectral_grid.Grid, spectral_grid.nlat_half)
    RingGrids.update_locator!(interpolator, londs, latds)

    regrid_flux!(atmos.diagnostic_variables.physics.sensible_heat_flux, Qs, interpolator, weights)
    regrid_flux!(atmos.diagnostic_variables.physics.evaporative_flux,   Mv, interpolator, weights)

    return nothing
end

function regrid_flux!(speedy_flux, climaocean_flux, interpolator, weights)
    interpolate!(tmp_field, climaocean_flux, weights) 
    tmp_field = on_architecture(CPU(), interior(tmp_field, :, :, 1))
    tmp_field[isnan.(tmp_field)] .= 0 # We do not have antarctic
    flux = FullClenshawGrid(vec(reverse(tmp_field, dims=2)))
    RingGrids.interpolate!(speedy_flux, flux, interpolator)
    return nothing
end
