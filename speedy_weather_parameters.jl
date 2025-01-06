using Thermodynamics.Parameters: AbstractThermodynamicsParameters

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
