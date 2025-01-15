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
Nz    = 10
h     = 3 # e-folding length of the exponential

r_faces = ClimaOcean.exponential_z_faces(; Nz, h, depth)
z_faces = ZStarVerticalCoordinate(r_faces)

Nx = 256 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 128 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz   = length(r_faces) - 1
grid = TripolarGrid(arch; size=(Nx, Ny, Nz), z=z_faces, halo=(4, 4, 4))

url = "https://www.dropbox.com/scl/fi/zy1cu64ybn93l67rjgiq0/Downsampled_ETOPO_2022.nc?rlkey=5upqvoxrnljj205amqf663vcw&st=ou8b32tt&dl=0"
filename = isfile("Downsampled_ETOPO_2022.nc") ? "Downsampled_ETOPO_2022.nc" : download(url, "Downsampled_ETOPO_2022.nc")
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, filename, dir="./")

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

momentum_advection = WENOVectorInvariant(order=3) 
tracer_advection   = Centered()

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

spectral_grid = SpectralGrid(trunc=31, nlayers=8)
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
SpeedyWeather.set!(atmosphere.model.time_stepping, Δt=Minute(30))

#####
##### Coupling the two
#####

Δt = atmosphere.model.time_stepping.Δt_sec

radiation = ClimaOcean.Radiation(ocean_albedo=0)
similarity_theory = SimilarityTheoryTurbulentFluxes(grid; maxiter=5)
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

# ## Visualizing the results
#
# We can visualize the results using CairoMakie. We record a video of surface variables and fluxes.
# To load the data we can use Oceananigans' `FieldTimeSeries` object.

# using JLD2
# using Oceananigans
# using Oceananigans.Grids: halo_size
# using CairoMakie 
# using Statistics: mean

# file  = jldopen("free_surface.jld2")
# iters = keys(file["timeseries/t"]) 

# Hx, Hy, _ = halo_size(η.grid)
# T  = FieldTimeSeries("surface_fields.jld2", "T")
# S  = FieldTimeSeries("surface_fields.jld2", "S")
# s  = FieldTimeSeries("surface_fields.jld2", "s")

# n  = Observable(1)
# Tn = @lift(interior(T[$n], :, :, 1))
# Sn = @lift(interior(S[$n], :, :, 1))
# sn = @lift(interior(s[$n], :, :, 1))
# ηn = @lift(file["timeseries/η/" * iters[$n]][Hx+1:end-Hx, Hy+1:end-Hy, 1])

# fig = Figure(size = (1800, 800))
# axT = Axis(fig[1, 1], title="Surface temperature ᵒC")
# axS = Axis(fig[1, 2], title="Surface salinity psu")
# axs = Axis(fig[2, 1], title="Surface speed ms⁻¹")
# axη = Axis(fig[2, 2], title="Sea surface height m")

# λ, φ, z = nodes(T[1])

# hmT = heatmap!(axT, Tn, colormap=:magma,  colorrange=(-1, 30))
# hmS = heatmap!(axS, Sn, colormap=:haline, colorrange=(25, 40))
# hms = heatmap!(axs, sn, colormap=:deep,   colorrange=( 0, 0.8))
# hmη = heatmap!(axη, ηn, colormap=:bwr,    colorrange=(-1, 1))

# CairoMakie.record(fig, "surface_fields.mp4", 1:length(T.times); framerate=5) do i 
#     @info "doing $i of $(length(T.times))"
#     n[] = i
# end

# # let's also visualize the surface fluxes that force the model

# Q  = FieldTimeSeries("surface_fluxes.jld2", "Q")
# τx = FieldTimeSeries("surface_fluxes.jld2", "τx")
# τy = FieldTimeSeries("surface_fluxes.jld2", "τy")
# PE = FieldTimeSeries("surface_fluxes.jld2", "PE")

# Qn  = @lift(interior(Q[$n],  :, :, 1))
# τxn = @lift(interior(τx[$n], :, :, 1))
# τyn = @lift(interior(τy[$n], :, :, 1))
# PEn = @lift(interior(PE[$n], :, :, 1))

# fig  = Figure(size = (1800, 800))
# axQ  = Axis(fig[1, 1], title="Net heat flux Wm⁻²")
# axPE = Axis(fig[1, 2], title="Net salt flux psu m s⁻¹")
# axτx = Axis(fig[2, 1], title="Zonal wind stress Nm⁻²")
# axτy = Axis(fig[2, 2], title="Meridional wind stress Nm⁻²")

# hmQ  = heatmap!(axQ,  Qn,  colormap=:magma,   colorrange=(-800,  800))
# hmPE = heatmap!(axPE, PEn, colormap=:haline,  colorrange=(-1e-5, 5e-5))
# hmτx = heatmap!(axτx, τxn, colormap=:balance, colorrange=(-5e-4, 5e-4))
# hmτy = heatmap!(axτy, τyn, colormap=:balance, colorrange=(-5e-4, 5e-4))

# CairoMakie.record(fig, "surface_fluxes.mp4", 1:length(Q.times); framerate=5) do i 
#     @info "doing $i of $(length(Q.times))"
#     n[] = i
# end

# # Let's visualize the internal structure of temperature and salinity
# x, y, z = nodes(ocean.model.tracers.T)

# fig = Figure(size = (1200, 400))
# axT = Axis(fig[1, 1], title="Internal temperature structure ᵒC")
# axS = Axis(fig[1, 2], title="Internal salinity structure psu")

# contourf!(axT, 1:48, z, interior(mean(ocean.model.tracers.T, dims=1), 1, :, :), colormap=:magma)
# contourf!(axS, 1:48, z, interior(mean(ocean.model.tracers.S, dims=1), 1, :, :), colormap=:haline)

# # Running a high-resolution simulation
#
# What are the steps to modify the above script to run an eddying (quarter degree) simulation?
