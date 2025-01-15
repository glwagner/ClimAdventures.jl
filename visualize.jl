using GLMakie
using Oceananigans
using ClimaOcean
using OrthogonalSphericalShellGrids
using NCDatasets

s  = FieldTimeSeries("surface_fields.jld2", "s";  backend=InMemory(10))
Qs = FieldTimeSeries("surface_fluxes.jld2", "Q";  backend=InMemory(10))
τx = FieldTimeSeries("surface_fluxes.jld2", "τx"; backend=InMemory(10))
τy = FieldTimeSeries("surface_fluxes.jld2", "τy"; backend=InMemory(10))
ds = Dataset("run_0013/output.nc", "r")

n  = Observable(1)
ut = XFaceField(s.grid)
vt = XFaceField(s.grid)
τ  = Field(sqrt(ut^2 + vt^2))

τxn = @lift begin
    set!(ut, τx[$n])
    set!(vt, τy[$n])
    compute!(τ)
    interior(τ, :, :, 1)
end

sn  = @lift interior(s[$n],  :, :, 1)
Qn  = @lift interior(Qs[$n], :, :, 1)
ζn  = @lift ds["zeta"].var[:, :, $n]

fig = Figure()
axs = Axis(fig[1, 1], title="ClimaOcean surface speed ms⁻¹")
axζ = Axis(fig[1, 1], title="SpeedyWeather surface vorticity s⁻¹")
axQ = Axis(fig[1, 2], title="Surface heat flux Wm⁻²")
axτ = Axis(fig[2, 1], title="Surface wind stress Nm⁻²")

hm = heatmap!(axs, sn, colormap=:deep,    colorrange=(0, 0.8))
hm = heatmap!(axζ, ζn, colormap=:viridis, colorrange=(-1e-4, 1e-4))
hm = heatmap!(axQ, Qn, colormap=:magma,   colorrange=(-800, 800))
hm = heatmap!(axQ, τn, colormap=:solar,   colorrange=(-1, 1))

record(fig, "surface_fields.mp4", 1:length(s.times); framerate=5) do i
    @info "doing $i of $(length(s.times))"
    n[] = i
end
