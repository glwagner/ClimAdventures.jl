#= Figure this out later

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

# Extend Thermodynamics functions to work with SpeedyWeather
const SWP = SpeedyWeatherParameters

@inline R_d(p::SWP)            = p.parameters.R_dry
@inline R_v(p::SWP)            = p.parameters.R_vapour
@inline LH_v0(p::SWP)          = p.parameters.latent_heat_condensation
@inline LH_s0(p::SWP)          = p.parameters.latent_heat_sublimation

@inline e_int_v0(p::SWP)       = LH_v0(p) - R_v(p) * T_0(p)

@inline cp_d(p::SWP)           = R_d(p) / kappa_d(p)
@inline cv_d(p::SWP)           = cp_d(p) - R_d(p)
@inline cv_v(p::SWP)           = cp_v(p) - R_v(p)

# Fix all these guys over here
@inline gas_constant(p::SWP)   = gas_constant(p.parameters)
@inline molmass_dryair(p::SWP) = 28.7
@inline molmass_water(p::SWP)  = 18.0
@inline molmass_ratio(p::SWP)  = 28.7 / 18.0

@inline LH_f0(p::SWP)          = LH_v0(p)
@inline T_freeze(p::SWP)       = 0
@inline T_triple(p::SWP)       = 0
@inline T_icenuc(p::SWP)       = 0
@inline pow_icenuc(p::SWP)     = 1
@inline press_triple(p::SWP)   = 1
@inline T_0(p::SWP)            = 0

@inline cp_v(p::SWP)           = 1000
@inline cp_l(p::SWP)           = 1000
@inline cp_i(p::SWP)           = 1000

@inline cv_l(p::SWP)           = 1000
@inline cv_i(p::SWP)           = 1000
@inline kappa_d(p::SWP)        = 1.0

=#