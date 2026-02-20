# =======================================================
# === simulation_mpd_base.jl (UNIFICADO Y CORREGIDO MPD) ===
# =======================================================

using DifferentialEquations
using LinearAlgebra
using Plots
using Printf

# =========================================================
# === CONSTANTES GLOBALES =================================
# =========================================================
const M  = 1.0
const mu = 1.0
const q  = 0.0
const Phi = 0.0

# =========================================================
# === COMPONENTES DE LA MÉTRICA SCHWARZSCHILD =============
# =========================================================
g_tt(r) = -(1.0 - 2.0*M/r)
g_rr(r) = 1.0 / (1.0 - 2.0*M/r)
g_thth(r) = r^2
g_phph(r, th) = (r*sin(th))^2
grr(r) = (1.0 - 2.0*M/r)
dg_tt(r) = 2.0*M/r^2

# =========================================================
# === POTENCIAL EFECTIVO Y AUXILIARES =====================
# =========================================================

function shX_m(r, th, J, s)
    f = grr(r)
    sint = sin(th)
    
    epsilon = 1e-12 * (1 + 10*s^2)
    denom = (mu*r)^2 - s^2 * f + epsilon

    term_A = (mu * J * r * sint) / denom
    term_B_numer = (J^2 - s^2) * f + (2.0*M/r) * J^2 * sint^2
    term_B = term_B_numer / denom
    inside_sqrt = term_A^2 - term_B

    if inside_sqrt < 0.0
        return term_A
    end
    
    return term_A - sqrt(inside_sqrt)
end

chX_m(r, th, J, s) = sqrt(shX_m(r, th, J, s)^2 + 1.0)

function Veff_m(r, th, J, s)
    f = grr(r)
    sint = sin(th)
    f_sqrt = sqrt(abs(f))

    shX = shX_m(r, th, J, s)
    chX = chX_m(r, th, J, s)

    t1 = mu * f_sqrt * chX
    bracket = (J * sint / (mu * r)) - shX
    fraction = (M * shX) / (f_sqrt * r * chX)
    t2 = mu * fraction * bracket

    return t1 + t2
end

# =========================================================
# === ECUACIONES DE MOVIMIENTO MPD - COMPLETAS (ODE!) =====
# =========================================================

function ODE!(dq, q, p, τ)
    t, r, th, ph, pt, pr, pth, pph, Str, Stth, Stph, Srth, Srph, Sthph = q
    J, s = p
    
    # Regularización
    r = max(r, 2.001 * M)
    th = max(1e-4, min(π - 1e-4, th))

    sin_th = sin(th)
    cos_th = cos(th)
    sin_th_sq = sin_th^2
    
    f_metric = 1.0 - 2.0*M/r
    
    # Conexiones (Christoffel) y Riemann Tensor (R) como en EoM.py
    f_eom_py = 2.0*M - r 
    
    # Christoffel Symbols (Grtt, Grrr, Grthth, etc.)
    Gttr = -M / (f_eom_py * r); Grtt = -(f_eom_py * M) / r^3
    Grrr = M / (f_eom_py * r); Grthth = f_eom_py; Grphph = f_eom_py * sin_th_sq
    Gthrth = 1.0 / r; Gthphph = -cos_th * sin_th
    Gphrph = 1.0 / r; Gphthph = cos_th / sin_th
    
    # Riemann Tensor dddd (R_mu_nu_rho_sigma)
    R_trtr = -2.0*M / r^3; R_tthtth = -f_eom_py*M / r^2
    R_tphtph = -f_eom_py*M * sin_th_sq / r^2; R_rtrt = -2.0*M / r^3 
    R_rthrth = M / f_metric; R_rphrph = M * sin_th_sq / f_metric
    R_thttht = -f_eom_py*M / r^2; R_thrthr = M / f_metric
    R_thphthph = 2.0*M * r * sin_th_sq; R_phtpht = -f_eom_py*M * sin_th_sq / r^2
    R_phrphr = M * sin_th_sq / f_metric; R_phthphth = 2.0*M * r * sin_th_sq

    # Riemann Tensor Uddd (R^mu_nu_rho_sigma)
    Rt_rtr = -2.0*M / (f_eom_py * r^2); Rt_thtth = -M / r; Rt_phtph = -M * sin_th_sq / r
    Rr_trt = 2.0*f_eom_py*M / r^4; Rr_thrth = -M / r; Rr_phrph = -M * sin_th_sq / r
    Rth_ttht = -f_eom_py*M / r^4; Rth_rthr = M / (f_eom_py * r^2); Rth_phthph = 2.0*M * sin_th_sq / r 
    Rph_tpht = -f_eom_py*M * sin_th_sq / r^4; Rph_rphr = M / (f_eom_py * r^2); Rph_thphth = 2.0*M / r 

    # Coeficientes Auxiliares (a, b, c, d)
    a = - Str*(R_rtrt*pt*Str - R_rthrth*pth*Srth - R_rphrph*pph*Srph) - Stth*(R_thttht*pt*Stth + R_thrthr*pr*Srth - R_thphthph*pph*Sthph) - Stph*(R_tphtph*pt*Stph + R_phrphr*pr*Srph + R_phthphth*pth*Sthph)
    b = - Str*(R_rtrt*pr*Str + R_tthtth*pth*Stth + R_tphtph*pph*Stph) - Srth*(R_thttht*pt*Stth + R_thrthr*pr*Srth - R_thphthph*pph*Sthph) - Srph*(R_tphtph*pt*Stph + R_phrphr*pr*Srph + R_phthphth*pth*Sthph) 
    c = - Stth*(R_rtrt*pr*Str + R_tthtth*pth*Stth + R_tphtph*pph*Stph) + Srth*(R_rtrt*pt*Str - R_thrthr*pth*Srth - R_rphrph*pph*Srph) - Sthph*(R_tphtph*pt*Stph + R_phrphr*pr*Srph + R_phthphth*pth*Sthph)
    d = - Stph*(R_rtrt*pr*Str + R_tthtth*pth*Stth + R_tphtph*pph*Stph) + Srph*(R_rtrt*pt*Str - R_thrthr*pth*Srth - R_rphrph*pph*Srph) + Sthph*(R_tthtth*pt*Stth + R_thrthr*pr*Srth + R_phthphth*pph*Sthph)
    
    # Coeficientes A y Delta
    A = g_tt(r)*a^2 + g_rr(r)*b^2 + g_thth(r)*c^2 + g_phph(r,th)*d^2
    Delta = 1.0 + (R_trtr*Str^2 + R_tthtth*Stth^2 + R_tphtph*Stph^2 + R_rthrth*Srth^2 + R_rphrph*Srph^2 + R_thphthph*Sthph^2) / mu^2
    
    # Factor de Normalización N/mu
    mu_sq = mu^2; Delta_sq = Delta^2
    denom_sq = mu_sq - A/(4.0 * mu_sq * Delta_sq)
    N_mu = (denom_sq <= 1e-12) ? 1.0 / sqrt(1e-12) : 1.0 / sqrt(denom_sq)
    
    # 4-Velocidad (v^mu)
    term_S_R_P_S = 1.0 / (2.0 * mu_sq * Delta)
    vt = N_mu * (pt + term_S_R_P_S * a); vr = N_mu * (pr + term_S_R_P_S * b)
    vth = N_mu * (pth + term_S_R_P_S * c); vph = N_mu * (pph + term_S_R_P_S * d)

    # Derivadas del Momento (dP^mu/dtau)
    dpt = -Gttr*(vt*pr + vr*pt) - Rt_rtr*vr*Str - Rt_thtth*vth*Stth - Rt_phtph*vph*Stph 
    dpr = -Grtt*vt*pt - Grrr*vr*pr - Grthth*vth*pth - Grphph*vph*pph  + Rr_trt*vt*Str - Rr_thrth*vth*Srth - Rr_phrph*vph*Srph
    dpth = -Gthrth*(vr*pth + vth*pr) - Gthphph*vph*pph + Rth_ttht*vt*Stth + Rth_rthr*vr*Srth - Rth_phthph*vph*Sthph
    dpph = -Gphrph*(vr*pph + vph*pr) - Gphthph*(vth*pph + vph*pth) + Rph_tpht*vt*Stph + Rph_rphr*vr*Srph + Rph_thphth*vth*Sthph

    # Derivadas del Spin (dS^mu_nu/dtau)
    dStr = -Gttr*vr*Str - Grrr*vr*Str - Grthth*vth*Stth - Grphph*vph*Stph + pt*vr - pr*vt
    dStth = -Gttr*(vt*Srth + vr*Stth) - Gthrth*(vr*Stth + vth*Str) - Gthphph*vph*Stph + pt*vth - pth*vt
    dStph = -Gttr*(vt*Srph + vr*Stph) - Gphrph*(vr*Stph + vph*Str) - Gphthph*(vth*Stph + vph*Stth) + pt*vph - pph*vt
    dSthph = -Gthrth*(vr*Sthph + vth*Srph) - Gphrph*(vr*Sthph - vph*Srth) - Gphthph*vth*Sthph + pth*vph - pph*vth
    dSrth = -Grtt*vt*Stth - Grrr*vr*Srth + Grphph*vph*Sthph - Gthrth*vr*Srth - Gthphph*vph*Srph  + pr*vth - pth*vr
    dSrph = -Grtt*vt*Stph - Grrr*vr*Srph - Grthth*vth*Sthph - Gphrph*vr*Srph - Gphthph*(vth*Srph + vph*Srth)  + pr*vph - pph*vr


    dq[1] = vt; dq[2] = vr; dq[3] = vth; dq[4] = vph
    dq[5] = dpt; dq[6] = dpr; dq[7] = dpth; dq[8] = dpph
    dq[9] = dStr; dq[10] = dStth; dq[11] = dStph; dq[12] = dSrth
    dq[13] = dSrph; dq[14] = dSthph
    
    return nothing
end

# =========================================================
# === CONDICIONES INICIALES (IC) ==========================
# =========================================================

function initial_conditions(t0, r0, th0, ph0, E, J, s, alpha)
    # Caso geodésico simple para s=0 (blindaje)
    if abs(s) < 1e-12
        pt0 = E / g_tt(r0); pr0 = 0.0; pth0 = 0.0
        pph0 = J / g_phph(r0, th0)
        return [t0, r0, th0, ph0, pt0, pr0, pth0, pph0, zeros(6)...]
    end

    # ESTABILIZACIÓN PARA SPINES GRANDES
    f = max((J * sin(alpha))^2, 1e-12) + (mu * r0)^2
    
    numerator = E + q * Phi + J * s * cos(alpha) * dg_tt(r0) / (2 * sqrt(max(-g_tt(r0) * g_rr(r0) * f, 1e-12)))
    denominator = -1 + r0 * s^2 * dg_tt(r0) / (2 * g_tt(r0) * g_rr(r0) * max(f, 1e-12))
    
    p_t0 = numerator / denominator
    pt0 = p_t0 / g_tt(r0)

    f_sqrt = sqrt(max(f, 1e-12))
    Srth0 = -p_t0 * s * sin(alpha) / f_sqrt
    Srph0 = p_t0 * s * cos(alpha) / f_sqrt
    Sthph0 = J * cos(th0) / max(sin(th0) * r0^2, 1e-12)

    p_th0 = -Srth0 * r0
    pth0 = p_th0 / g_thth(r0)
    p_ph0 = (J - Srph0 * r0) * sin(th0)^2
    pph0 = p_ph0 / g_phph(r0, th0)

    # Cálculo ROBUSTO de pr0 (normalización del 4-momento)
    energy_term = -(mu^2 + g_tt(r0)*pt0^2 + g_thth(r0)*pth0^2 + g_phph(r0, th0)*pph0^2)
    pr_squared = energy_term / g_rr(r0)
    
    pr0 = (pr_squared < 0) ? 0.01 : sqrt(pr_squared)
    p_r0 = pr0 * g_rr(r0)

    # Componentes restantes del spin tensor
    Str0  = -(p_th0^2 + (p_ph0 / max(sin(th0), 1e-12))^2 - J * p_ph0) / (r0 * p_t0)
    Stth0 = (p_r0 * p_th0 + J * p_ph0 * cos(th0) / (r0 * max(sin(th0), 1e-12))) / (r0 * p_t0)
    Stph0 = -(p_r0 * (J - p_ph0 / max(sin(th0)^2, 1e-12)) + J * p_th0 * cos(th0) / (r0 * max(sin(th0), 1e-12))) / (r0 * p_t0)

    # Límite heurístico para evitar NaN si los componentes son demasiado grandes
    clamp_s(val) = max(min(val, s + 0.1), -(s + 0.1))

    return [t0, r0, th0, ph0, pt0, pr0, pth0, pph0, 
            clamp_s(Str0), clamp_s(Stth0), clamp_s(Stph0),
            clamp_s(Srth0), clamp_s(Srph0), clamp_s(Sthph0)]
end

# =========================================================
# === SIMULACIÓN BASE Y SECCIÓN DE POINCARÉ ================
# =========================================================

function run_base_simulation(J, s; r0=6.0, alpha=pi-0.7, tmax=10000000.0)
    
    th0 = pi/2; t0, ph0 = 0.0, 0.0
    r_pts, pr_pts = Float64[], Float64[]

    @printf "INICIANDO: J=%.1f, s=%.2f, r0=%.1f, τ_max=%.0e\n" J s r0 tmax

    # 1. Calcular Energía e IC
    #E = Veff_m(r0, th0, J, s)
    #E = 0.97698396 #This is for s = 0.4
    #E = 0.96730999 #This is for s = 0.6
    E = 0.94738162 #This is for s = 1.0
    #E = 0.93545565 #This is for s = 1.2
    #E = 0.92292941 #This is for s = 1.4
    q0 = initial_conditions(t0, r0, th0, ph0, E, J, s, alpha)
    
    # 2. Seleccionar Integrador (Robusto)
    integrator = (s <= 0.6) ? Rodas5P() : KenCarp4()
    
    # 3. Callbacks para Poincaré y Horizonte
    condition(u, t, integrator) = u[3] - pi/2 # Cruce theta=pi/2
    affect!(integrator) = (push!(r_pts, integrator.u[2]); push!(pr_pts, integrator.u[6]))
    cb_poincare = ContinuousCallback(condition, affect!, nothing, save_positions=(false, false), abstol=1e-12)

    cb_stop = ContinuousCallback((u,t,i) -> u[2] - 2.0001*M, (i) -> terminate!(i), save_positions=(false, false))
    full_callback = CallbackSet(cb_poincare, cb_stop)

    # 4. Solución
    prob = ODEProblem(ODE!, q0, (0.0, tmax), (J, s)) 
    sol = solve(prob, integrator, 
                reltol=1e-9, abstol=1e-13, 
                callback=full_callback, 
                maxiters=1e9) 

    @printf "Puntos de Poincaré obtenidos: %d\n" length(r_pts)
    return r_pts, pr_pts, sol
end

# =========================================================
# === EJECUCIÓN DEL CASO DEMO =============================
# =========================================================

J_DEMO = 4.0
S_DEMO = 1.0  # Spin alto para demostrar la estabilidad de KenCarp4
R0_DEMO = 6.0 # Radio estable

r_pts, pr_pts, sol = run_base_simulation(J_DEMO, S_DEMO; r0=R0_DEMO)

# =========================================================
# === BLOQUE DE GRÁFICAS 2D Y SECUENCIALES ================
# =========================================================

# NOTA: Asegúrate de que J_DEMO y S_DEMO sean los valores usados para la simulación
J_plot = J_DEMO
S_plot = S_DEMO

# 1. Calcular Coordenadas Cartesianas
x = sol[2,:].*sin.(sol[3,:]) .* cos.(sol[4,:])
y = sol[2,:].*sin.(sol[3,:]) .* sin.(sol[4,:])
z = sol[2,:].*cos.(sol[3,:]) # Z ya no se usa para 2D

# --- GRÁFICA 1: TRAYECTORIA 2D (PLANO X-Y) ---
plot2d = plot(x, y, 
              xlabel="x", 
              ylabel="y", 
              title=@sprintf("Trayectoria MPD 2D (J=%.1f, s=%.2f)", J_plot, S_plot),
              aspect_ratio=:equal, # Crucial para ver la órbita correctamente
              legend=false)

# Mostrar la gráfica 2D (Plano Ecuatorial)
display(plot2d)

# --- GRÁFICA 3D: TRAYECTORIA COMPLETA ---
plot3d_traj = plot(
    x, y, z,
    xlabel = "x",
    ylabel = "y",
    zlabel = "z",
    title = @sprintf("Trayectoria MPD 3D (J=%.1f, s=%.2f)", J_plot, S_plot),
    linewidth = 1.5,
    legend = false,
    aspect_ratio = :equal
)

display(plot3d_traj)

# Opcional: guardar la figura 3D
savefig(plot3d_traj, "trayectoria_3d_mpd.png")



# --- GRÁFICA 2: SECCIÓN DE POINCARÉ (r - p_r) ---
scatter_poincare = scatter(r_pts, pr_pts, 
                           xlabel="r", 
                           ylabel="\$p_r\$", # Uso de LaTeX para p_r
                           title=@sprintf("Sección de Poincaré (J=%.1f, s=%.2f)", J_plot, S_plot), 
                           markersize=2.5, 
                           markerstrokewidth=0, 
                           legend=false,
                           color=:orange,
                           alpha=0.7)

# Mostrar la gráfica de Poincaré
display(scatter_poincare)

# Opcional: Si quieres guardar la gráfica 2D por separado
savefig(plot2d, "trayectoria_2d_mpd.png")

