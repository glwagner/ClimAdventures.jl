using ClimaOcean
using Oceananigans
using Oceananigans.Units
using OrthogonalSphericalShellGrids
using SpeedyWeather
using CFTime
using Dates
using Printf
using CairoMakie

####
#### A near-global ocean
####

arch  = Oceananigans.CPU()
depth = 5000meters
Nz    = 40
h     = 30 # e-folding length of the exponential

r_faces = ClimaOcean.exponential_z_faces(; Nz, h, depth)
z_faces = ZStarVerticalCoordinate(r_faces)

Nx = 1440 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 700 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz   = length(r_faces) - 1
grid = TripolarGrid(arch; size=(Nx, Ny, Nz), z=z_faces, halo=(4, 4, 4))

url = "https://www.dropbox.com/scl/fi/zy1cu64ybn93l67rjgiq0/Downsampled_ETOPO_2022.nc?rlkey=5upqvoxrnljj205amqf663vcw&st=ou8b32tt&dl=0"
filename = isfile("Downsampled_ETOPO_2022.nc") ? "Downsampled_ETOPO_2022.nc" : download(url, "Downsampled_ETOPO_2022.nc")
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, filename, dir="./")

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

momentum_advection = WENOVectorInvariant() 
tracer_advection   = WENO()

free_surface = SplitExplicitFreeSurface(grid; substeps=30) 

using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity, 
                                       ExplicitTimeDiscretization, 
                                       DiffusiveFormulation

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity

numerical_closure = HorizontalScalarDiffusivity(ν=5e3)
eddy_closure      = IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3, skew_flux_formulation=DiffusiveFormulation())
vertical_mixing   = CATKEVerticalDiffusivity() 

closure = (eddy_closure, numerical_closure, vertical_mixing) 

ocean = ocean_simulation(grid; 
                         momentum_advection, 
                         tracer_advection, 
                         closure, 
                         free_surface)

# Set up initial conditions for temperature and salinity
ClimaOcean.set!(ocean.model, T=ECCOMetadata(:temperature),
                             S=ECCOMetadata(:salinity))

#####
##### Adding an atmosphere
#####

include("method_extension.jl")

spectral_grid = SpectralGrid(trunc=31, nlayers=20)
surface_heat_flux = PrescribedSurfaceHeatFlux()
surface_evaporation = PrescribedSurfaceEvaporation()
model = PrimitiveWetModel(spectral_grid; surface_heat_flux, surface_evaporation)
atmosphere = initialize!(model)
model.output.active = true
add!(model.output, SpeedyWeather.SurfaceFluxesOutput()...)

# Initializing the atmosphere (necessary to set up correct initial conditions)
SpeedyWeather.initialize!(atmosphere)
SpeedyWeather.set_period!(atmosphere.prognostic_variables.clock, Day(400))
SpeedyWeather.initialize!(atmosphere.prognostic_variables.clock, atmosphere.model.time_stepping)
SpeedyWeather.set!(atmosphere.model.time_stepping, Δt=Minute(10))

#####
##### Coupling the two
#####

Δt = atmosphere.model.time_stepping.Δt_sec

radiation = ClimaOcean.Radiation(ocean_albedo=0)
similarity_theory = SimilarityTheoryTurbulentFluxes(grid; maxiter=10)
sea_ice = ClimaOcean.FreezingLimitedOceanTemperature()
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, similarity_theory)
earth = ClimaOcean.Simulation(earth_model; Δt, stop_time=400days)

# ### Adding some diagnostics
#
# We add a callback to save surface fields as well as surface fluxes, every 6 hours

u, v, _ = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
s = sqrt(u^2 + v^2)

η = ocean.model.free_surface.η 

earth.output_writers[:surface_tracers] = JLD2OutputWriter(ocean.model, (; T, S, s),
                                                          schedule = TimeInterval(6hours),
                                                          indices = (:, :, grid.Nz),
                                                          overwrite_existing = true,
                                                          filename = "surface_fields.jld2")

earth.output_writers[:free_surface] = JLD2OutputWriter(ocean.model, (; η),
                                                       schedule = TimeInterval(6hours),
                                                       overwrite_existing = true,
                                                       filename = "free_surface.jld2")

Q  = earth.model.fluxes.total.ocean.heat
τx = earth.model.fluxes.total.ocean.momentum.u
τy = earth.model.fluxes.total.ocean.momentum.v
PE = earth.model.fluxes.total.ocean.tracers.S

earth.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, (; Q, τx, τy, PE),
                                                 schedule = TimeInterval(6hours),
                                                 overwrite_existing = true,
                                                 filename = "surface_fluxes.jld2")

# Also, we add a callback to print a message about how the simulation is going

wall_time = [time_ns()]

function progress(earth)
    clock = earth.model.clock
    atmos = earth.model.atmosphere
    time  = Second(atmos.prognostic_variables.clock.time - atmos.prognostic_variables.clock.start).value

    maxu = maximum(abs, earth.model.ocean.model.velocities.u)
    maxv = maximum(abs, earth.model.ocean.model.velocities.v)
    maxT = maximum(earth.model.ocean.model.tracers.T)
    minS = minimum(earth.model.ocean.model.tracers.S)
    
    @info @sprintf("Iteration: %d, time O: %s, time A: %s, wall_time: %s, max(|u|, |v|): %.2e %.2e max(T): %.2e, min(S): %.2e\n",
                   clock.iteration, prettytime(clock.time), prettytime(time), prettytime(1e-9 * (time_ns() - wall_time[1])), maxu, maxv, maxT, minS)

    wall_time[1] = time_ns()
end

add_callback!(earth, progress, IterationInterval(10))

# ### Running the simulation
#
# quite simply

ClimaOcean.run!(earth)
