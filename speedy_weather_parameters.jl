using Thermodynamics.Parameters: AbstractThermodynamicsParameters

import Thermodynamics.Parameters:
    gas_constant,   #
    molmass_dryair, # Molar mass of dry air (without moisture)
    molmass_water,  # Molar mass of gaseous water vapor
    molmass_ratio,  # Ratio of the molar masses of dry air to water vapor
    R_v,            # Specific gas constant for water vapor
    R_d,            # Specific gas constant for dry air
    kappa_d,        # Ideal gas adiabatic exponent for dry air
    T_0,            # Enthalpy reference temperature
    LH_v0,          # Vaporization enthalpy at the reference temperature
    LH_s0,          # Sublimation enthalpy at the reference temperature
    LH_f0,          # Fusion enthalpy at the reference temperature
    cp_d,           # Heat capacity of dry air at constant pressure
    cp_v,           # Isobaric specific heat capacity of gaseous water vapor
    cp_l,           # Isobaric specific heat capacity of liquid water
    cp_i,           # Isobaric specific heat capacity of water ice
    cv_v,           # Heat capacity of dry air at constant volume
    cv_l,           # Isobaric specific heat capacity of liquid water
    cv_i,           # Isobaric specific heat capacity of liquid water
    e_int_v0,       # what? someting about reference internal energy of water vapor
    T_freeze,       # Freezing temperature of _pure_ water
    T_triple,       # Triple point temperature of _pure_ water
    press_triple,   # Triple point pressure of pure water
    T_icenuc,       # Lower temperature limit for the presence of liquid condensate
                    # (below which homogeneous ice nucleation occurs)
    pow_icenuc      # "Power parameter" that controls liquid/ice condensate partitioning
                    # during partial ice nucleation

struct SpeedyWeatherParameters{FT, P} <: AbstractThermodynamicsParameters{FT}
    parameters :: P
    SpeedyWeatherParameters{FT}(params::P) = new{FT, P}(params)
end

# Extend Thermodynamics functions to work with SpeedyWeather
const SWP = SpeedyWeatherParameters

@inline R_d(p::SWP)            = R_d(s.parameters.R_dry)
@inline R_v(p::SWP)            = R_v(s.parameters.R_vapour)
@inline gas_constant(p::SWP)   = gas_constant(p.constitutive)
@inline molmass_dryair(p::SWP) = molmass_dryair(p.constitutive)
@inline molmass_water(p::SWP)  = molmass_water(p.constitutive)
@inline molmass_ratio(p::SWP)  = molmass_ratio(p.constitutive)
@inline LH_v0(p::SWP)          = LH_v0(p.parameters.latent_heat_condensation)
@inline LH_s0(p::SWP)          = LH_s0(p.parameters.latent_heat_sublimation)
@inline LH_f0(p::SWP)          = LH_f0(p.phase_transitions)
@inline T_freeze(p::SWP)       = T_freeze(p.phase_transitions)
@inline T_triple(p::SWP)       = T_triple(p.phase_transitions)
@inline T_icenuc(p::SWP)       = T_icenuc(p.phase_transitions)
@inline pow_icenuc(p::SWP)     = pow_icenuc(p.phase_transitions)
@inline press_triple(p::SWP)   = press_triple(p.phase_transitions)
@inline T_0(p::SWP)            = T_0(p.phase_transitions)

@inline e_int_v0(p::SWP)       = LH_v0(p) - R_v(p) * T_0(p)

@inline cp_v(p::SWP)           = cp_v(p.heat_capacity)
@inline cp_l(p::SWP)           = cp_l(p.heat_capacity)
@inline cp_i(p::SWP)           = cp_i(p.heat_capacity)

@inline cv_l(p::SWP)           = cv_l(p.heat_capacity)
@inline cv_i(p::SWP)           = cv_i(p.heat_capacity)

@inline kappa_d(p::SWP)        = kappa_d(p.heat_capacity)
@inline cp_d(p::SWP)           = R_d(p) / kappa_d(p)
@inline cv_d(p::SWP)           = cp_d(p) - R_d(p)
@inline cv_v(p::SWP)           = cp_v(p) - R_v(p)
