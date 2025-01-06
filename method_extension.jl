# This file contains all the extensions needed from SpeedyWeather to run the coupled model.

using ClimaOcean
using SpeedyWeather

using ClimaOcean.OceanSeaIceModels.Atmospheres: HeatCapacityParameters, 
                                                ConstitutiveParameters

import Oceananigans: time_step!

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

# Extend the PrognosticAtmosphere functions to work with SpeedyWeather
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

# Make sure the atmospheric parameters from SpeedyWeather can be used in the compute fluxes function
import ClimaOcean.OceanSeaIceModels.Atmospheres: thermodynamics_parameters, 
                                                 boundary_layer_height, 
                                                 surface_layer_height

function thermodynamics_parameters(atmos::SpeedyWeather.Simulation)
    return HeatCapacityParameters()
end