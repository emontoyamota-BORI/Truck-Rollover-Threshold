# =============================================================================
#  paso3_rampa.py — UNA corrida del tractocamión en la rampa del p22 / ref. [11].
#  Equivale a la co-simulación TruckSim–Simulink del paper (Sec. 3.1):
#    - el conductor sigue el eje de la rampa (ChPathFollowerDriver)
#    - velocidad objetivo: constante o el perfil v_final(s) de las Ecs. 1–19
#    - registra cada 0.05 s: cargas y fuerzas laterales por rueda, ángulos de
#      deslizamiento por eje, LTR total y por eje, balanceo de tractor y remolque,
#      dirección, articulación, aceleración lateral y desplazamiento lateral
#    - si video = 1: abre la ventana 3D, guarda fotogramas y arma GIF (y MP4)
#
#  Uso:  python paso3_rampa.py <susp> <modo> <v_kmh> <mu> <video 0|1>
#        p. ej.  python paso3_rampa.py aire constante 40 0.8 1
#  (correr_todo.py lo llama con todos los casos de config.CASOS_RAMPA)
#
#  CONVENCIÓN DE SIGNOS en la telemetría (la rampa gira a la DERECHA):
#    lado L = EXTERIOR de la curva, lado R = INTERIOR.
#    a_y, F_y, dirección: + = hacia el CENTRO de la curva.
#    balanceo: + = hacia el EXTERIOR, medido RELATIVO A LA VÍA (sin el peralte).
#    balanceo_susp_eje<i>: masa suspendida relativa al eje i (lo que gira la suspensión).
#    Fz = carga normal a la vía; Fy = fuerza lateral en el plano de la vía.
#    LTR = (ext − int)/(ext + int).
#    desv_lat: + = hacia el exterior (izquierda) del eje de la rampa.
# =============================================================================
import sys, os, math, time, json, shutil
import numpy as np
import config as P

SUSP = sys.argv[1] if len(sys.argv) > 1 else "aire"
MODO = sys.argv[2] if len(sys.argv) > 2 else "constante"
V_KMH = float(sys.argv[3]) if len(sys.argv) > 3 else 40.0
MU = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
VIDEO = (sys.argv[5] == "1") if len(sys.argv) > 5 else True
caso = dict(modo=MODO, v_kmh=V_KMH, mu=MU, video=VIDEO)
NOMBRE = P.nombre_corrida(caso)

import pychrono as chrono
import pychrono.vehicle as veh
import chrono_modelo as C
import video as VID

dir_res = P.resultados_susp(SUSP)
dir_tel = dir_res / "telemetria"; dir_tel.mkdir(parents=True, exist_ok=True)

eje = np.loadtxt(P.datos_susp(SUSP) / "p22" / "terreno" / "eje.csv", delimiter=",", skiprows=1)
S_EJE, X_EJE, Y_EJE, Z_EJE, K_EJE = eje[:, 0], eje[:, 1], eje[:, 2], eje[:, 3], eje[:, 4]
A_EJE, B_EJE = eje[:, 5], eje[:, 6]                       # pendiente y peralte [rad]
PSI_EJE = np.unwrap(np.arctan2(np.gradient(Y_EJE), np.gradient(X_EJE)))   # rumbo del eje [rad]
perfil = np.loadtxt(P.PERFIL_P22, delimiter=",", skiprows=1) if MODO.startswith("perfil") else None


def v_objetivo(s):
    if MODO == "constante":
        return V_KMH / 3.6
    col = 1 if MODO == "perfil_mu_constante" else 2
    v0 = perfil[0, col]
    return float(np.interp(s, perfil[:, 0], perfil[:, col], left=v0, right=perfil[-1, col]))


def estacion(x, y, i_prev):
    a, b = max(i_prev - 30, 0), min(i_prev + 60, len(S_EJE) - 1)
    d2 = (X_EJE[a:b] - x) ** 2 + (Y_EJE[a:b] - y) ** 2
    i = a + int(np.argmin(d2))
    j = min(i + 1, len(S_EJE) - 1)
    tx, ty = X_EJE[j] - X_EJE[i], Y_EJE[j] - Y_EJE[i]
    n = math.hypot(tx, ty) or 1.0
    lateral = (-ty * (x - X_EJE[i]) + tx * (y - Y_EJE[i])) / n
    return S_EJE[i], lateral, i


# ----------------------------------------------------------------------------- modelo
print(f"=== Suspensión: {SUSP} | corrida: {NOMBRE} | video: {VIDEO} ===")
R = C.resumen(SUSP)
i0 = int(np.argmin(abs(S_EJE - P.S_INICIO)))
DZ0 = float(os.environ.get("P22_DZ_INICIAL", P.DZ_INICIAL[0]))
print(f"Altura inicial sobre el piso: {DZ0:.3f} m | asentamiento: {P.T_ASENTAR:.1f} s frenado")
pos0 = chrono.ChVector3d(X_EJE[i0], Y_EJE[i0], Z_EJE[i0] + R["radio_cargado_del_m"] + DZ0)
tractor, remolque = C.crear_tractor_remolque(SUSP, pos0, vis=VIDEO)
sis = tractor.GetSystem()
terreno = C.crear_terreno(sis, malla=True, mu=MU)

puntos = chrono.vector_ChVector3d()
for k in range(i0, len(S_EJE), 2):
    puntos.push_back(chrono.ChVector3d(X_EJE[k], Y_EJE[k], Z_EJE[k]))
trayectoria = chrono.ChBezierCurve(puntos, False)
conductor = veh.ChPathFollowerDriver(tractor, trayectoria, "rampa_p22", max(v_objetivo(P.S_INICIO), 1.0))
conductor.GetSteeringController().SetLookAheadDistance(P.LOOKAHEAD)
conductor.GetSteeringController().SetGains(*P.GANANCIAS_DIRECCION)
conductor.GetSpeedController().SetGains(*P.GANANCIAS_VELOCIDAD)
conductor.Initialize()

vis = None
if VIDEO:
    try:
        import pychrono.irrlicht as irr
        vis = veh.ChWheeledVehicleVisualSystemIrrlicht()
        vis.SetWindowTitle(f"p22 — {P.SUSPENSIONES[SUSP]['descripcion']} — {NOMBRE}")
        vis.SetWindowSize(*P.VENTANA_PX)
        vis.SetChaseCamera(chrono.ChVector3d(*P.CAMARA["punto"]), P.CAMARA["distancia"], P.CAMARA["altura"])
        vis.Initialize()
        vis.AddLightDirectional()
        vis.AddSkyBox()
        vis.AttachVehicle(tractor)
        vis.AttachDriver(conductor)
    except Exception as e:
        print("Sin ventana 3D (no habrá fotogramas):", e)
        vis = None

carpeta_cuadros = None
if vis is not None:
    carpeta_cuadros = P.VIDEOS_ASCII / SUSP / NOMBRE      # Chrono no escribe en rutas con acentos
    carpeta_cuadros.mkdir(parents=True, exist_ok=True)
    for viejo in carpeta_cuadros.glob("cuadro_*.jpg"):
        viejo.unlink()
    print(f"Fotogramas -> {carpeta_cuadros}")

# ----------------------------------------------------------------------------- columnas
COLS = (["t", "s", "x", "y", "z", "desv_lat", "kappa", "v", "v_obj", "ay_geom", "ay_vr", "r_guinada",
         "balanceo_tractor_deg", "balanceo_remolque_deg", "articulacion_deg", "direccion_deg",
         "direccion_cmd", "acelerador", "freno"]
        + [f"Fz{lado}{i}" for i in range(1, 5) for lado in ("L", "R")]
        + [f"Fy{lado}{i}" for i in range(1, 5) for lado in ("L", "R")]
        + [f"alfa{i}_deg" for i in range(1, 5)]
        + ["LTR"] + [f"LTR{i}" for i in range(1, 5)]
        + [f"balanceo_susp_eje{i}_deg" for i in range(1, 5)])
IX = {c: k for k, c in enumerate(COLS)}

# ----------------------------------------------------------------------------- simulación
reg, info_cuadros = [], []
i_eje = i0
cada = max(1, int(round(P.DT_REGISTRO / P.DT)))
cada_render = max(1, int(round(1 / 50 / P.DT)))
cada_cuadro = max(1, int(round(1 / P.FPS_FOTOGRAMAS / P.DT)))
paso, n_cuadro = 0, 0
despegue = [None] * 4          # (t, s, v) del primer LTR_i >= LTR_DESPEGUE sostenido DESPEGUE_MIN_S
racha = [0] * 4
n_racha = max(1, int(round(P.DESPEGUE_MIN_S / P.DT_REGISTRO)))
v_ref = min(P.V_REF_INICIAL, v_objetivo(P.S_INICIO))   # objetivo con rampa de aceleración (evita sobrepaso)
i_rem = i0                     # índice del eje de la rampa más cercano al remolque
quieto = veh.DriverInputs()    # entradas durante el asentamiento: frenado, sin acelerador ni dirección
quieto.m_steering = 0.0; quieto.m_throttle = 0.0; quieto.m_braking = 1.0
motivo = "llegó a S_FIN"
t_reloj = time.time()
ultima = None
while True:
    t = sis.GetChTime()
    if vis is not None:
        if not vis.Run():
            motivo = "ventana cerrada"; break
        grabar = carpeta_cuadros is not None and paso % cada_cuadro == 0 and S_EJE[i_eje] >= P.S_GRABAR_DESDE
        if paso % cada_render == 0 or grabar:
            vis.BeginScene(); vis.Render(); vis.EndScene()
        if grabar:
            vis.WriteImageToFile((carpeta_cuadros / f"cuadro_{n_cuadro:05d}.jpg").as_posix())
            u = ultima if ultima is not None else [0.0] * len(COLS)
            info_cuadros.append([n_cuadro, t, S_EJE[i_eje], tractor.GetSpeed(), u[IX["LTR"]],
                                 u[IX["desv_lat"]], u[IX["balanceo_remolque_deg"]], u[IX["ay_geom"]]])
            n_cuadro += 1

    p = tractor.GetPos()
    s, lat, i_eje = estacion(p.x, p.y, i_eje)
    asentando = t < P.T_ASENTAR
    if asentando:
        entradas = quieto
    else:
        v_ref = min(v_ref + P.ACEL_OBJETIVO * P.DT, v_objetivo(s + P.ANTICIPO))
        conductor.SetDesiredSpeed(v_ref)
        entradas = conductor.GetInputs()
        conductor.Synchronize(t)
    terreno.Synchronize(t)
    tractor.Synchronize(t, entradas, terreno); remolque.Synchronize(t, entradas, terreno)
    if vis is not None:
        vis.Synchronize(t, entradas)
    if not asentando:
        conductor.Advance(P.DT)
    terreno.Advance(P.DT); tractor.Advance(P.DT); remolque.Advance(P.DT)
    if vis is not None:
        vis.Advance(P.DT)
    paso += 1

    if paso % cada == 0:
        ch_t = tractor.GetChassisBody(); ch_r = remolque.GetChassis().GetBody()
        pr = ch_r.GetPos()
        _, _, i_rem = estacion(pr.x, pr.y, i_rem)
        n_t = C.normal_via(PSI_EJE[i_eje], A_EJE[i_eje], B_EJE[i_eje])
        n_r = C.normal_via(PSI_EJE[i_rem], A_EJE[i_rem], B_EJE[i_rem])
        Ft = C.fuerzas_por_eje(tractor, terreno, ch_t, n_t)
        Fr = C.fuerzas_por_eje(remolque, terreno, ch_r, n_r)
        F = Ft + Fr
        ext = sum(f[0] for f in F); inte = sum(f[1] for f in F)
        ltr = (ext - inte) / max(ext + inte, 1.0)
        ltr_i = [(f[0] - f[1]) / max(f[0] + f[1], 1.0) for f in F]
        v = tractor.GetSpeed()
        kap = float(np.interp(s, S_EJE, K_EJE))
        try:
            delta = -math.degrees(C.angulo_direccion(tractor))
        except Exception:
            delta = float("nan")
        r_g = -C.velocidad_guinada(ch_t)
        fila = [t, s, p.x, p.y, p.z, lat, kap, v, v_ref, v * v * kap, v * r_g, r_g,
                -math.degrees(C.balanceo_rel_via(ch_t, n_t)), -math.degrees(C.balanceo_rel_via(ch_r, n_r)),
                math.degrees(C.angulo_articulacion(tractor, remolque)), delta,
                entradas.m_steering, entradas.m_throttle, entradas.m_braking]
        fila += [x for f in F for x in (f[0], f[1])]
        fila += [x for f in F for x in (-f[2], -f[3])]
        fila += [math.degrees(f[4]) for f in F]
        fila += [ltr] + ltr_i
        try:
            fila += [-math.degrees(b) for b in C.balanceo_suspension(tractor, remolque)]
        except Exception:
            fila += [float("nan")] * 4
        reg.append(fila); ultima = fila
        for i in range(4):
            racha[i] = racha[i] + 1 if (ltr_i[i] >= P.LTR_DESPEGUE and s >= 0) else 0
            if despegue[i] is None and racha[i] >= n_racha:
                despegue[i] = (t, s, v * 3.6)
                print(f"  ** Despega el eje {i+1} en s = {s:.1f} m a {v*3.6:.1f} km/h (LTR ≥ {P.LTR_DESPEGUE} durante {P.DESPEGUE_MIN_S} s)")
        if len(reg) % 40 == 0:
            print(f"t = {t:6.1f} s | s = {s:7.1f} m | v = {v*3.6:5.1f} km/h | LTR = {ltr:5.2f} | "
                  f"balanceo rem. = {fila[IX['balanceo_remolque_deg']]:5.2f}° | desv. lat. = {lat:+.2f} m")
        if t < P.T_ASENTAR + 3.0 and abs(delta) > 20.0:
            motivo = "FALLA AL ARRANCAR: ruedas delanteras giradas > 20° sin orden de dirección"; break
        if max(abs(fila[IX["balanceo_tractor_deg"]]), abs(fila[IX["balanceo_remolque_deg"]])) > P.BALANCEO_VUELCO_DEG:
            motivo = "VUELCO"; break
    if s >= P.S_FIN:
        break
    if t > P.T_MAX:
        motivo = "tiempo máximo"; break

print(f"\nFin: {motivo}. Tiempo simulado {sis.GetChTime():.1f} s en {time.time()-t_reloj:.0f} s de reloj.")
A = np.array(reg)
np.savetxt(dir_tel / f"{NOMBRE}.csv", A, delimiter=",", header=",".join(COLS), comments="", fmt="%.6g")
print(f"Telemetría -> {dir_tel / (NOMBRE + '.csv')}")

resumen = {"suspension": SUSP, "corrida": NOMBRE, "caso": caso, "fin": motivo, "dz_inicial_m": DZ0,
           "vuelco": motivo == "VUELCO",
           "s_final_m": float(A[-1, IX["s"]]) if len(A) else None,
           "v_final_kmh": float(A[-1, IX["v"]] * 3.6) if len(A) else None,
           "despegue_por_eje": [None if d is None else {"t": d[0], "s": d[1], "v_kmh": d[2]} for d in despegue],
           "fotogramas": n_cuadro, "carpeta_fotogramas": str(carpeta_cuadros) if carpeta_cuadros else None}
with open(dir_tel / f"{NOMBRE}_resumen.json", "w", encoding="utf-8") as f:
    json.dump(resumen, f, indent=2, ensure_ascii=False)

if carpeta_cuadros is not None and n_cuadro > 0:
    np.savetxt(carpeta_cuadros / "cuadros.csv", np.array(info_cuadros), delimiter=",",
               header="cuadro,t,s,v,LTR,desv_lat,balanceo_remolque_deg,ay_geom", comments="", fmt="%.5f")
    if P.COPIAR_FOTOGRAMAS:
        destino = dir_res / "fotogramas" / NOMBRE
        if destino.exists():
            shutil.rmtree(destino)
        shutil.copytree(carpeta_cuadros, destino)
        print(f"{n_cuadro} fotogramas copiados a {destino}")
    titulo = f"{P.SUSPENSIONES[SUSP]['descripcion']} — {NOMBRE}"
    for archivo in VID.hacer_video(carpeta_cuadros, dir_res / "videos", NOMBRE, titulo, P.SUSPENSIONES[SUSP]["color"]):
        print("Video ->", archivo)
