# =======================================================
# === comparacion_s08.jl (ANÁLISIS DE ESTABILIDAD EN R) ===
# =======================================================

using Plots
using Printf

# !!! Importa tu archivo principal (donde está run_base_simulation) !!!
include("Untitled-1.jl") 

# =======================================================
# === PARÁMETROS SOLICITADOS ============================
# =======================================================

const J_VAL = 4.0      # Momento angular fijo (asumiendo tu valor por defecto)
const S_VAL = 1.4     # Spin fijo solicitado
const T_MAX = 1.0e3    # Tiempo de simulación solicitado (1e4)

# Radios iniciales a evaluar
radios = [6.0, 7.0, 8.0]
colores = [:orange, :blue, :green] # Colores para distinguir cada R

# Inicializar la gráfica (vacía al inicio)
p = plot(title=@sprintf("Secciones de Poincaré para s=%.2f (τ=%.0e)", S_VAL, T_MAX),
         xlabel="r", 
         ylabel="\$p_r\$", # Usando LaTeX
         #legend=:topright, 
         legend = false,
         size=(800, 600))

println("INICIANDO ANÁLISIS DE ESTABILIDAD MPD EN FUNCIÓN DEL RADIO...")

# =======================================================
# === BUCLE PRINCIPAL DE SIMULACIÓN =====================
# =======================================================

for (i, r0) in enumerate(radios)
    color_actual = colores[i]

    @printf "\n--- Ejecutando r0 = %.1f ---\n" r0

    # Llama a la función de simulación con los parámetros fijos y el radio actual
    # La función debe devolver (r_pts, pr_pts, sol)
    r_pts, pr_pts, sol = run_base_simulation(J_VAL, S_VAL; r0=r0, tmax=T_MAX)

    n_pts = length(r_pts)

    if n_pts > 0
        @printf "✅ ÉXITO: %d puntos obtenidos.\n" n_pts
        
        # Añadir los puntos de Poincaré a la misma gráfica
        scatter!(p, r_pts, pr_pts,
                 color=color_actual, 
                 markersize=3.0,
                 label=@sprintf("r₀=%.1f (%d pts)", r0, n_pts),
                 markerstrokewidth=0,
                 alpha=0.7)
    else
        @printf "❌ FALLA: No se obtuvieron puntos de Poincaré o la órbita cayó al horizonte.\n"
        # Opcional: Puede agregar un marcador para indicar dónde falló la órbita (e.g., r0=3.7)
    end
end

# Mostrar y guardar la gráfica final con todas las secciones
display(p)
savefig(p, @sprintf("poincare_comparacion_s%.2f_Rset.png", S_VAL))

println("\nGráfica final guardada como 'poincare_comparacion_s0.80_Rset.png'")

# =======================================================
# === FIGURA ZOOM (DETALLE DE POINCARÉ) =================
# =======================================================

p_zoom = deepcopy(p)

# Límites del zoom
xlims!(p_zoom, (7.0, 10.0))     # r
ylims!(p_zoom, (-0.05, 0.05))      # p_r

# Opcional: mejorar legibilidad en el zoom
plot!(p_zoom,
      title = @sprintf("Zoom – Secciones de Poincaré (s=%.2f)", S_VAL),
      markersize = 4.0)

display(p_zoom)

savefig(p_zoom, @sprintf("poincare_zoom_s%.2f_r11_25.png", S_VAL))
