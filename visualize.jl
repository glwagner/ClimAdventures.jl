using GLMakie
using Oceananigans
using ClimaOcean
using OrthogonalSphericalShellGrids
using NCDatasets

s  = FieldTimeSeries("surface_fields.jld2", "s";  backend=InMemory(10))
Qs = FieldTimeSeries("surface_fluxes.jld2", "Q";  backend=InMemory(10))
τx = FieldTimeSeries("surface_fluxes.jld2", "τx"; backend=InMemory(10))
τy = FieldTimeSeries("surface_fluxes.jld2", "τy"; backend=InMemory(10))
ds = Dataset("run_0016/output.nc", "r")

n  = Observable(1)
ut = XFaceField(s.grid)
vt = YFaceField(s.grid)
τ  = Field(sqrt(ut^2 + vt^2))

τn = @lift begin
    Oceananigans.set!(ut, τx[$n])
    Oceananigans.set!(vt, τy[$n])
    compute!(τ)
    interior(τ, :, :, 1)
end

τxn = @lift interior(τx[$n], :, :, 1)
sn  = @lift interior(s[$n],  :, :, 1)
ζn  = @lift reverse(ds["vor"].var[:, :, 8, $n], dims=2)
Qn  = @lift reverse(ds["sensible_heat_flux"].var[:, :, $n], dims=2)

fig = Figure()
axs = Axis(fig[1, 1], title="ClimaOcean surface speed ms⁻¹")
axζ = Axis(fig[2, 1], title="SpeedyWeather surface vorticity s⁻¹")
axQ = Axis(fig[1, 2], title="Sensible heat flux Wm⁻²")
axτ = Axis(fig[2, 2], title="Surface wind stress Nm⁻²")

hm = heatmap!(axs, sn,  colormap=:deep,    colorrange=(0, 0.8))
hm = heatmap!(axζ, ζn,  colormap=:viridis, colorrange=(-7e-5, 7e-5))
hm = heatmap!(axQ, Qn,  colormap=:magma,   colorrange=(-10, 10))
hm = heatmap!(axτ, τxn, colormap=:bwr,     colorrange=(-2.5, 2.5).*1e-5)

record(fig, "surface_fields.mp4", 1:length(s.times)-1; framerate=5) do i
    @info "doing $i of $(length(s.times))"
    n[] = i
end
