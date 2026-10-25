# Verification plots for temperature_dependent_kappa (Whittington et al., 2009)
# usage: julia plot_verification.jl <scratch_dir>
using LaMEM, CairoMakie, Printf

S = ARGS[1]

# Whittington et al. (2009), T in K, kappa in mm^2/s
κ(TK) = TK < 846 ? 567.3/TK - 0.062 : 0.732 - 0.000135TK
# Kirchhoff integral F(T) = ∫ κ dT from 273.15 K
lo(t) = 567.3log(t) - 0.062t
hi(t) = 0.732t - 0.0000675t^2
F(TK) = TK < 846 ? lo(TK) - lo(273.15) : lo(846) - lo(273.15) + hi(TK) - hi(846)

# steady state of ∇·(k(T)∇T) = 0 with k = κ(T)ρCp: F(T(z)) is linear in z
function T_exact(z; Ttop=273.15, Tbot=1273.15, L=10.0)
    target = F(Tbot) * (-z) / L
    a, b = Ttop, Tbot
    for _ in 1:100
        m = (a + b) / 2
        F(m) < target ? (a = m) : (b = m)
    end
    return (a + b) / 2 - 273.15
end

ink, muted = "#1f1f1e", "#6b6a63"
c_td, c_ck = "#2a78d6", "#eb6834"      # categorical slots 1 and 2
set_theme!(fontsize = 15, Axis = (xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
           leftspinecolor = muted, bottomspinecolor = muted, topspinevisible = false,
           rightspinevisible = false))

#-------------------------------------------------------------------------------
# Figure 1: 1D steady-state conduction vs exact solution
#-------------------------------------------------------------------------------
data, t = read_LaMEM_timestep("steady", 82, joinpath(S, "val"))
z  = data.z.val[1, 1, :]
T  = data.fields.temperature[1, 1, :]
kc = data.fields.conductivity[1, 1, :]

zf  = range(-10, 0, length = 400)
Tex = T_exact.(zf)
Tl  = 1000 .* (-zf) ./ 10
err = maximum(abs.(T .- T_exact.(z)))
kex = [κ(Ti + 273.15) * 1e-6 * 2700 * 1000 for Ti in T]
kerr = maximum(abs.(kc[2:end-1] .- kex[2:end-1]) ./ kex[2:end-1])
@printf("1D: t = %.2f Myr, max |T - T_exact| = %.3f C, max |T_exact - linear| = %.1f C, max rel k error = %.2e\n",
        t[1], err, maximum(abs.(Tex .- Tl)), kerr)

fig1 = Figure(size = (900, 470))
ax1 = Axis(fig1[1, 1], xlabel = "Temperature [°C]", ylabel = "z [km]",
           title = "Steady-state temperature", titlealign = :left)
lines!(ax1, Tl, collect(zf), color = c_ck, linewidth = 2, linestyle = :dash, label = "constant κ (linear)")
lines!(ax1, Tex, collect(zf), color = ink, linewidth = 2, label = "exact, κ(T)")
scatter!(ax1, T, z, color = c_td, markersize = 8, strokecolor = :white, strokewidth = 1, label = "LaMEM, κ(T)")
axislegend(ax1, position = :rt, framevisible = false)
text!(ax1, 20, -9.4, text = @sprintf("max |T − T_exact| = %.2f °C", err), color = muted, fontsize = 13)

ax2 = Axis(fig1[1, 2], xlabel = "Conductivity [W/m/K]", ylabel = "z [km]",
           title = "Conductivity k = κ(T)·ρ·Cp", titlealign = :left)
Tg = range(0, 1000, length = 400)
lines!(ax2, [κ(Ti + 273.15) * 1e-6 * 2700 * 1000 for Ti in T_exact.(zf)], collect(zf),
       color = ink, linewidth = 2, label = "κ(T_exact)·ρ·Cp")
scatter!(ax2, kc, z, color = c_td, markersize = 8, strokecolor = :white, strokewidth = 1, label = "LaMEM output")
axislegend(ax2, position = :rb, framevisible = false)
Label(fig1[0, :], "1D column, 10 km, T = 0 °C top / 1000 °C bottom, ρ = 2700 kg/m³, Cp = 1000 J/kg/K, 64 cells",
      color = muted, fontsize = 13, halign = :left, tellwidth = false)
save(joinpath(S, "fig_1D_steady.png"), fig1, px_per_unit = 2)

#-------------------------------------------------------------------------------
# Figure 2: t43 test model, magma sphere (700 C) in host crust (500 C)
#-------------------------------------------------------------------------------
function centre_history(dir, name)
    times, Tc = Float64[], Float64[]
    for step in 0:10
        d, tt = read_LaMEM_timestep(name, step, dir)
        x, y, zz = d.x.val[:, 1, 1], d.y.val[1, :, 1], d.z.val[1, 1, :]
        i, j, k = argmin(abs.(x)), argmin(abs.(y)), argmin(abs.(zz .+ 6))
        push!(times, tt[1]); push!(Tc, d.fields.temperature[i, j, k])
    end
    return times, Tc
end
function profile(dir, name)
    d, _ = read_LaMEM_timestep(name, 10, dir)
    x, y = d.x.val[:, 1, 1], d.y.val[1, :, 1]
    i, j = argmin(abs.(x)), argmin(abs.(y))
    return d.z.val[1, 1, :], d.fields.temperature[i, j, :], d.fields.conductivity[i, j, :]
end
ttd, Ttd = centre_history(joinpath(S, "sphere_td"), "t43_TempDepKappa")
tck, Tck = centre_history(joinpath(S, "sphere_ck"), "t43_TempDepKappa")
ztd, Ptd, Ktd = profile(joinpath(S, "sphere_td"), "t43_TempDepKappa")
zck, Pck, Kck = profile(joinpath(S, "sphere_ck"), "t43_TempDepKappa")
@printf("sphere: centre T at %.2f Myr: kappa(T) %.1f C, constant k %.1f C\n", ttd[end], Ttd[end], Tck[end])
@printf("sphere: host k along centre line: kappa(T) %.2f..%.2f W/m/K\n",
        minimum(Ktd[abs.(ztd .+ 6) .> 3]), maximum(Ktd[abs.(ztd .+ 6) .> 3]))

fig2 = Figure(size = (900, 470))
ax3 = Axis(fig2[1, 1], xlabel = "Time [Myr]", ylabel = "Temperature [°C]",
           title = "Temperature at sphere centre", titlealign = :left)
lines!(ax3, tck, Tck, color = c_ck, linewidth = 2, label = "host k = 2.5 W/m/K (constant)")
scatter!(ax3, tck, Tck, color = c_ck, markersize = 8, strokecolor = :white, strokewidth = 1)
lines!(ax3, ttd, Ttd, color = c_td, linewidth = 2, label = "host κ(T), Whittington et al. (2009)")
scatter!(ax3, ttd, Ttd, color = c_td, markersize = 8, strokecolor = :white, strokewidth = 1)
axislegend(ax3, position = :rt, framevisible = false)

ax4 = Axis(fig2[1, 2], xlabel = "Temperature [°C]", ylabel = "z [km]",
           title = @sprintf("Vertical profile through centre, t = %.1f Myr", ttd[end]), titlealign = :left)
hspan!(ax4, -8.5, -3.5, color = (:black, 0.05))
text!(ax4, 502, -3.9, text = "magma sphere", color = muted, fontsize = 12)
lines!(ax4, Pck, zck, color = c_ck, linewidth = 2, label = "host k = 2.5 (constant)")
lines!(ax4, Ptd, ztd, color = c_td, linewidth = 2, label = "host κ(T)")
axislegend(ax4, position = :rb, framevisible = false)
Label(fig2[0, :], "t43_TempDepKappa: 700 °C sphere (r = 2.5 km, k = 1.5 W/m/K) in 500 °C host crust, 16³ cells, 10 steps",
      color = muted, fontsize = 13, halign = :left, tellwidth = false)
save(joinpath(S, "fig_sphere_t43.png"), fig2, px_per_unit = 2)
