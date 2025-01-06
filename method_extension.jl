# This file contains all the extensions needed from SpeedyWeather to run the coupled model.

using ClimaOcean
using SpeedyWeather

using ClimaOcean.OceanSeaIceModels.Atmospheres: HeatCapacityParameters, 
                                                ConstitutiveParameters

import Oceananigans: time_step!

#####
##### Extending the time_step! function and the `update_model_field_time_series!` function
#####

# Out-source the time_step! to the prognostic atmosphere model
function time_step!(atmos::SpeedyWeather.Simulation) 
    progn = atmos.simulation.prognostic_variables
    diagn = atmos.simulation.diagnostic_variables
    model = atmos.simulation.model
    
    (; clock) = progn
    (; Δt, Δt_millisec) = model.time_stepping

    (; output, feedback) = model
    
    SpeedyWeather.timestep!(progn, diagn, 2Δt, model) # calculate tendencies and leapfrog forward
    SpeedyWeather.timestep!(clock, Δt_millisec)       # time of lf=2 and diagn after timestep!

    SpeedyWeather.progress!(feedback, progn)          # updates the progress meter bar
    SpeedyWeather.output!(output, progn, diagn, model)
    SpeedyWeather.callback!(model.callbacks, progn, diagn, model)

    return nothing
end

update_model_field_time_series!(::SpeedyWeather.Simulation, time) = nothing

####
#### Extending inputs to flux computation
####

# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.Atmospheres: thermodynamics_parameters, 
                                                 boundary_layer_height, 
                                                 surface_layer_height


# This should be the height of the surface layer in the atmospheric model
surface_layer_height(atmos::SpeedyWeather.Simulation) = 0

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::SpeedyWeather.Simulation) = 600

using SpeedyWeather: EarthAtmosphere

Base.eltype(::EarthAtmosphere{FT}) where FT = FT

thermodynamics_parameters(atmos::SpeedyWeather.Simulation) = SpeedyWeatherParameters{FT}(atmos.model.atmosphere)

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
    nothing
end

# Regrid the fluxes from the ocean/sea-ice grid to the atmospheric model grid
function regrid_fluxes_to_atmospheric_model!(atmos::SpeedyWeather.Simulation, net_tracer_fluxes, centered_velocity_fluxes)
    nothing
end
