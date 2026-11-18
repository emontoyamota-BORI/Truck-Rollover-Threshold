# =============================================================================
#  paso2_mesa_basculante.py — Prueba estática + mesa basculante virtual (UMTRI).
#   Parte A: el camión se asienta 3 s frenado en piso plano -> cargas por eje
#            contra la Tabla 1 del p22.
#   Parte B: se inclina la gravedad 1°/s (equivale a inclinar la mesa) hasta que
#            despegan las ruedas. a_y/g = tan θ. Se registra el despegue de cada
#            eje (SRT del primer eje), el balanceo de tractor y remolque y el
#            gradiente de balanceo [°/g]: es la referencia directa para el
#            modelo analítico de balanceo eje por eje (Parte B de la guía).
#  Uso:  python paso2_mesa_basculante.py multihoja      (o aire)
#        python paso2_mesa_basculante.py                 (corre las dos)
#  Salidas: resultados/<susp>/mesa/mesa_basculante.csv y resumen_mesa.json
# =============================================================================
import sys, math, json
import numpy as np
import config as P

if len(sys.argv) < 2:
    import lanzar
    for s in P.LISTA_SUSPENSIONES:
        lanzar.correr("paso2_mesa_basculante.py", s)
    sys.exit(0)

SUSP = sys.argv[1]
USAR_IRRLICHT = "--ver" in sys.argv     # ventana 3D opcional:  ... multihoja --ver

import pychrono as chrono
import pychrono.vehicle as veh
import chrono_modelo as C

DT = 1e-3
R = C.resumen(SUSP)
out = P.resultados_susp(SUSP) / "mesa"; out.mkdir(parents=True, exist_ok=True)
z_ini = R["radio_cargado_del_m"] + 0.03
tractor, remolque = C.crear_tractor_remolque(SUSP, chrono.ChVector3d(0, 0, z_ini), vis=USAR_IRRLICHT)
sis = tractor.GetSystem()
terreno = C.crear_terreno(sis, malla=False, mu=1.0)

vis = None
if USAR_IRRLICHT:
    try:
        import pychrono.irrlicht as irr
        vis = veh.ChWheeledVehicleVisualSystemIrrlicht()
        vis.SetWindowTitle(f"p22 — mesa basculante ({SUSP})")
        vis.SetWindowSize(1280, 800)
        vis.SetChaseCamera(chrono.ChVector3d(-6.0, 0.0, 1.5), 20.0, 0.5)
        vis.Initialize(); vis.AddLightDirectional(); vis.AddSkyBox(); vis.AttachVehicle(tractor)
    except Exception as e:
        print("Sin ventana 3D:", e); vis = None

entradas = veh.DriverInputs()
entradas.m_steering = 0.0; entradas.m_throttle = 0.0; entradas.m_braking = 1.0
n_paso = 0

def paso():
    global n_paso
    t = sis.GetChTime()
    if vis is not None:
        if not vis.Run():
            raise SystemExit("ventana cerrada")
        if n_paso % 20 == 0:
            vis.BeginScene(); vis.Render(); vis.EndScene()
    terreno.Synchronize(t); tractor.Synchronize(t, entradas, terreno); remolque.Synchronize(t, entradas, terreno)
    if vis is not None: vis.Synchronize(t, entradas)
    terreno.Advance(DT); tractor.Advance(DT); remolque.Advance(DT)
    if vis is not None: vis.Advance(DT)
    n_paso += 1

# ------------------------------------------------------------- Parte A
while sis.GetChTime() < 3.0:
    paso()
F = C.fuerzas_verticales(tractor, terreno) + C.fuerzas_verticales(remolque, terreno)
tot = sum(l + r for l, r in F)
print(f"\n=== {SUSP} — Parte A: cargas estáticas ===")
print(f"Peso medido {tot/P.G/1e3:.2f} t (esperado {P.masa_total()/1e3:.2f} t)")
estatica = []
for i, (l, r) in enumerate(F):
    pct = 100 * (l + r) / tot
    estatica.append({"eje": i + 1, "Fz_izq_kN": l / 1e3, "Fz_der_kN": r / 1e3, "porcentaje": pct, "p22_pct": 100 * P.F_EJE[i]})
    print(f"Eje {i+1}: {(l+r)/1e3:6.1f} kN = {pct:5.1f} % (p22 {100*P.F_EJE[i]:.0f} %)")

# ------------------------------------------------------------- Parte B
print(f"\n=== {SUSP} — Parte B: mesa basculante ({P.MESA_GRADOS_POR_S}°/s) ===")
est_eje = [l + r for l, r in F]
reg, despegue = [], [None] * 4
theta, t0 = 0.0, sis.GetChTime()
motivo = "ángulo máximo"
while theta < math.radians(P.MESA_ANGULO_MAX_DEG):
    theta = math.radians(P.MESA_GRADOS_POR_S) * (sis.GetChTime() - t0)
    sis.SetGravitationalAcceleration(chrono.ChVector3d(0, -P.G * math.sin(theta), -P.G * math.cos(theta)))
    paso()
    if n_paso % 50:
        continue
    F = C.fuerzas_verticales(tractor, terreno) + C.fuerzas_verticales(remolque, terreno)
    ltr = [(r - l) / max(r + l, 1.0) for l, r in F]               # la gravedad lateral carga el lado derecho
    der = sum(r for l, r in F); izq = sum(l for l, r in F)
    ltr_tot = (der - izq) / max(der + izq, 1.0)
    phi_t = math.degrees(C.angulo_balanceo(tractor.GetChassisBody()))   # + = lado izquierdo arriba = se inclina hacia el lado cargado
    phi_r = math.degrees(C.angulo_balanceo(remolque.GetChassis().GetBody()))
    ay = math.tan(theta)
    try:
        phi_s = [math.degrees(b) for b in C.balanceo_suspension(tractor, remolque)]
    except Exception:
        phi_s = [float("nan")] * 4
    reg.append([sis.GetChTime() - t0, math.degrees(theta), ay] + [f for lr in F for f in lr] + ltr + [ltr_tot, phi_t, phi_r] + phi_s)
    for i, (l, r) in enumerate(F):
        if despegue[i] is None and l < 0.02 * est_eje[i]:
            despegue[i] = ay
            print(f"  Despega el eje {i+1} a a_y/g = {ay:.3f} (balanceo remolque {phi_r:.1f}°)")
    if all(d is not None for d in despegue):
        motivo = "despegaron los 4 ejes"; break
    if max(abs(phi_t), abs(phi_r)) > 25:
        motivo = "balanceo > 25° (vuelco)"; break

A = np.array(reg)
cab = ("t,theta_deg,ay_g,FzL1,FzR1,FzL2,FzR2,FzL3,FzR3,FzL4,FzR4,"
       "LTR1,LTR2,LTR3,LTR4,LTR,balanceo_tractor_deg,balanceo_remolque_deg,"
       "balanceo_susp_eje1_deg,balanceo_susp_eje2_deg,balanceo_susp_eje3_deg,balanceo_susp_eje4_deg")
np.savetxt(out / "mesa_basculante.csv", A, delimiter=",", header=cab, comments="", fmt="%.5f")

# gradiente de balanceo en la zona lineal (a_y < 0.15 g)
lin = A[:, 2] < 0.15
grad_t = float(np.polyfit(A[lin, 2], A[lin, 16], 1)[0]) if lin.sum() > 3 else float("nan")
grad_r = float(np.polyfit(A[lin, 2], A[lin, 17], 1)[0]) if lin.sum() > 3 else float("nan")
grad_s = []
for i in range(4):
    try:
        grad_s.append(float(np.polyfit(A[lin, 2], A[lin, 18 + i], 1)[0]) if lin.sum() > 3 else float("nan"))
    except Exception:
        grad_s.append(float("nan"))
orden = [i + 1 for i in sorted(range(4), key=lambda i: despegue[i] if despegue[i] is not None else 9)]
validos = [d for d in despegue if d is not None]
res = {"suspension": SUSP, "estatica": estatica, "peso_medido_t": tot / P.G / 1e3,
       "despegue_ay_g_por_eje": despegue,
       "SRT_primer_despegue_g": min(validos) if validos else None,
       "ay_ultimo_despegue_g": max(validos) if len(validos) == 4 else None,
       "ay_max_alcanzado_g": float(A[:, 2].max()),
       "orden_despegue": orden, "fin": motivo,
       "gradiente_balanceo_tractor_deg_g": grad_t, "gradiente_balanceo_remolque_deg_g": grad_r,
       "gradiente_balanceo_suspension_por_eje_deg_g": grad_s,
       "SSF_aprox": (sum(P.VIA) / 4) / (2 * R["h_cg_combinacion_m"])}
with open(out / "resumen_mesa.json", "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, ensure_ascii=False)
print(f"\nFin: {motivo}")
print(f"SRT (primer despegue) = {res['SRT_primer_despegue_g']} g | orden de despegue: {orden}")
print(f"Gradiente de balanceo: tractor {grad_t:.1f} °/g, remolque {grad_r:.1f} °/g")
print(f"Resultados -> {out}")
