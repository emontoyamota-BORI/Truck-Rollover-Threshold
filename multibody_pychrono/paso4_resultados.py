# =============================================================================
#  paso4_resultados.py — Gráficas y tablas comparativas multihoja vs aire, para
#  comparar después con el modelo analítico (SRT físico eje por eje + MF).
#  No necesita PyChrono. Se puede correr las veces que quieras.
#
#  Lee:   resultados/<susp>/mesa/*            (paso2)
#         resultados/<susp>/telemetria/*.csv  (paso3)
#  Escribe en resultados/comparacion/:
#    figuras/  mesa_*.png, resumen_*.png y <corrida>_*.png
#    tablas/   tabla_mesa_basculante.csv     SRT, despegue por eje, gradiente de balanceo
#              tabla_mesa_curvas.csv         LTR_i(a_y) y balanceo(a_y) cada 0.01 g
#              tabla_rampa_estacionario.csv  promedios en las curvas de R = 125 m y 60 m
#              tabla_eventos.csv             despegue de ejes, vuelco, fin de cada corrida
#              tabla_validacion_ref11.csv    LTR a 30/40/50 km/h y vuelco vs TruckSim
#              resultados_p22_suspensiones.xlsx  (todas las tablas, si hay openpyxl)
#  y en resultados/<susp>/tablas/<corrida>_cada5m.csv (telemetría cada 5 m de estación)
# =============================================================================
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as P
import tablas as T

FIG = P.RESULTADOS / "comparacion" / "figuras"; FIG.mkdir(parents=True, exist_ok=True)
TAB = P.RESULTADOS / "comparacion" / "tablas"; TAB.mkdir(parents=True, exist_ok=True)
SUSPS = P.LISTA_SUSPENSIONES
COL = {s: P.SUSPENSIONES[s]["color"] for s in SUSPS}
DESC = {s: P.SUSPENSIONES[s]["descripcion"] for s in SUSPS}
C_EJE = ['tab:green', 'tab:orange', 'tab:purple', 'tab:brown']
REF_FIG8A_CTE = np.array([[0, 80], [30, 80], [100, 72], [200, 62], [330, 47], [400, 47.5], [480, 47.5], [600, 59]])


def guardar(nombre):
    plt.tight_layout(); plt.savefig(FIG / nombre, dpi=130); plt.close()
    print("  figura ->", nombre)


def nanmean(x):
    x = np.asarray(x, float)
    return float(np.nanmean(x)) if np.isfinite(x).any() else float("nan")


def nanmax(x):
    x = np.asarray(x, float)
    return float(np.nanmax(x)) if np.isfinite(x).any() else float("nan")


def p95(x):
    """Percentil 95: 'máximo' robusto. Los máximos crudos de Fz/LTR por eje traen picos de
       un solo registro (0.05 s) del contacto llanta–malla; el p95 los ignora."""
    x = np.asarray(x, float)
    return float(np.nanpercentile(x, 95)) if np.isfinite(x).any() else float("nan")


# ============================================================================ 1. MESA
mesa, mesa_res = {}, {}
for s in SUSPS:
    f = P.resultados_susp(s) / "mesa" / "mesa_basculante.csv"
    if f.exists():
        mesa[s] = T.leer_telemetria(f)
        with open(f.parent / "resumen_mesa.json", encoding="utf-8") as h:
            mesa_res[s] = json.load(h)

filas_mesa, filas_curvas = [], []
if mesa:
    print("Mesa basculante:")
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True, sharey=True)
    for i in range(4):
        a = ax.flat[i]
        for s, d in mesa.items():
            a.plot(d["ay_g"], d[f"LTR{i+1}"], color=COL[s], lw=2, label=DESC[s])
            dsp = mesa_res[s]["despegue_ay_g_por_eje"][i]
            if dsp is not None:
                a.axvline(dsp, color=COL[s], ls=":", lw=1)
        a.axhline(1, color="k", lw=.8, ls="--"); a.axhline(0.8, color="r", lw=.8, ls="--")
        a.set_title(f"Eje {i+1}"); a.grid(alpha=.3)
        if i >= 2: a.set_xlabel("a_y / g = tan θ")
        if i % 2 == 0: a.set_ylabel("LTR del eje [-]")
    ax.flat[0].legend(fontsize=8)
    fig.suptitle("Mesa basculante virtual: LTR por eje (líneas punteadas = despegue)")
    guardar("mesa_LTR_por_eje.png")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for s, d in mesa.items():
        ax[0].plot(d["ay_g"], d["balanceo_tractor_deg"], color=COL[s], lw=2, label=DESC[s])
        ax[1].plot(d["ay_g"], d["balanceo_remolque_deg"], color=COL[s], lw=2, label=DESC[s])
        for i in range(4):
            c = f"balanceo_susp_eje{i+1}_deg"
            if c in d:
                ax[2].plot(d["ay_g"], d[c], color=COL[s], lw=1.5, ls=["-", "--", "-.", ":"][i], label=f"{s} eje {i+1}")
    for a, t in zip(ax, ["Tractor (total, respecto al piso)", "Semirremolque (total)", "Suspensión: chasis − eje"]):
        a.set_title(t); a.set_xlabel("a_y / g"); a.grid(alpha=.3); a.set_ylabel("Balanceo [°]")
    ax[0].legend(fontsize=8); ax[2].legend(fontsize=7, ncol=2)
    guardar("mesa_balanceo.png")

    fig, ax = plt.subplots(figsize=(8, 4.3))
    w = 0.35
    for k, s in enumerate(mesa):
        dsp = [np.nan if v is None else v for v in mesa_res[s]["despegue_ay_g_por_eje"]]
        ax.bar(np.arange(1, 5) + (k - 0.5) * w, dsp, w, color=COL[s], label=DESC[s])
        for i, v in enumerate(dsp):
            if np.isfinite(v):
                ax.text(i + 1 + (k - 0.5) * w, v + 0.005, f"{v:.3f}", ha="center", fontsize=8)
    ax.set_xticks(range(1, 5)); ax.set_xticklabels([f"Eje {i}" for i in range(1, 5)])
    ax.set_ylabel("a_y / g al despegar"); ax.set_title("SRT por eje (mesa basculante)"); ax.grid(alpha=.3, axis="y")
    ax.legend(fontsize=8)
    guardar("mesa_SRT_por_eje.png")

    grid = np.round(np.arange(0.0, 0.501, 0.01), 3)
    for s, d in mesa.items():
        r = mesa_res[s]
        fila = {"suspension": s, "descripcion": DESC[s],
                "SRT_primer_despegue_g": r["SRT_primer_despegue_g"],
                "ay_ultimo_despegue_g": r["ay_ultimo_despegue_g"], "ay_max_alcanzado_g": r["ay_max_alcanzado_g"],
                "orden_despegue": "-".join(str(x) for x in r["orden_despegue"]), "fin": r["fin"],
                "gradiente_balanceo_tractor_deg_g": r["gradiente_balanceo_tractor_deg_g"],
                "gradiente_balanceo_remolque_deg_g": r["gradiente_balanceo_remolque_deg_g"],
                "SSF_aprox": r["SSF_aprox"], "peso_medido_t": r["peso_medido_t"]}
        for i, gs in enumerate(r.get("gradiente_balanceo_suspension_por_eje_deg_g", [np.nan] * 4)):
            fila[f"gradiente_balanceo_susp_eje{i+1}_deg_g"] = gs
        for i in range(4):
            fila[f"despegue_eje{i+1}_g"] = r["despegue_ay_g_por_eje"][i]
            fila[f"carga_estatica_eje{i+1}_pct"] = r["estatica"][i]["porcentaje"]
        filas_mesa.append(fila)
        o = np.argsort(d["ay_g"])
        for g in grid:
            if g > d["ay_g"].max():
                break
            fc = {"suspension": s, "ay_g": g}
            for i in range(4):
                fc[f"LTR{i+1}"] = float(np.interp(g, d["ay_g"][o], d[f"LTR{i+1}"][o]))
            fc["LTR_total"] = float(np.interp(g, d["ay_g"][o], d["LTR"][o]))
            fc["balanceo_tractor_deg"] = float(np.interp(g, d["ay_g"][o], d["balanceo_tractor_deg"][o]))
            fc["balanceo_remolque_deg"] = float(np.interp(g, d["ay_g"][o], d["balanceo_remolque_deg"][o]))
            for i in range(4):
                c = f"balanceo_susp_eje{i+1}_deg"
                if c in d:
                    fc[c] = float(np.interp(g, d["ay_g"][o], d[c][o]))
            for i in range(4):
                fc[f"FzL{i+1}_kN"] = float(np.interp(g, d["ay_g"][o], d[f"FzL{i+1}"][o])) / 1e3
                fc[f"FzR{i+1}_kN"] = float(np.interp(g, d["ay_g"][o], d[f"FzR{i+1}"][o])) / 1e3
            filas_curvas.append(fc)
    T.escribir_csv(TAB / "tabla_mesa_basculante.csv", filas_mesa)
    T.escribir_csv(TAB / "tabla_mesa_curvas.csv", filas_curvas)
else:
    print("(sin resultados de mesa basculante: corre paso2)")

# ============================================================================ 2. RAMPA
corridas = {}       # nombre -> {susp: (datos, resumen)}
for s in SUSPS:
    for f in sorted((P.resultados_susp(s) / "telemetria").glob("rampa_*.csv")):
        d = T.leer_telemetria(f)
        if "balanceo_susp_eje1_deg" not in d:      # archivo de la versión anterior: volver a correr paso3
            print(f"  AVISO: {s}/{f.name} es de la versión anterior; vuelve a correr esa corrida.")
            continue
        rj = f.with_name(f.stem + "_resumen.json")
        res = json.load(open(rj, encoding="utf-8")) if rj.exists() else {}
        corridas.setdefault(f.stem, {})[s] = (d, res)

filas_est, filas_ev = [], []
perfil = np.loadtxt(P.PERFIL_P22, delimiter=",", skiprows=1) if P.PERFIL_P22.exists() else None

for nombre, por_s in corridas.items():
    print(f"Corrida {nombre}:")
    # --- velocidad
    plt.figure(figsize=(10, 4.3))
    for s, (d, _) in por_s.items():
        m = d["s"] >= -20
        plt.plot(d["s"][m], d["v"][m] * 3.6, color=COL[s], lw=2, label=f"{DESC[s]} (real)")
        plt.plot(d["s"][m], d["v_obj"][m] * 3.6, color=COL[s], lw=.8, ls="--", label=f"{DESC[s]} (objetivo)")
    if "perfil_mu_constante" in nombre:
        plt.plot(*REF_FIG8A_CTE.T, "ko", mfc="w", label="p22 Fig. 8a, μ constante (leído)")
    plt.xlabel("Estación s [m]"); plt.ylabel("Velocidad [km/h]"); plt.grid(alpha=.3); plt.legend(fontsize=8)
    plt.title(f"{nombre}: velocidad"); guardar(f"{nombre}_velocidad.png")

    # --- LTR total y por eje
    fig, ax = plt.subplots(5, 1, figsize=(10, 11), sharex=True)
    for s, (d, _) in por_s.items():
        m = d["s"] >= 0
        ax[0].plot(d["s"][m], d["LTR"][m], color=COL[s], lw=2, label=DESC[s])
        for i in range(4):
            ax[i + 1].plot(d["s"][m], d[f"LTR{i+1}"][m], color=COL[s], lw=1.5, label=DESC[s])
    for k, a in enumerate(ax):
        a.axhline(0.8, color="r", ls="--", lw=.8); a.axhline(1.0, color="k", ls="--", lw=.8)
        a.set_ylabel("LTR total (Ec. 10)" if k == 0 else f"LTR eje {k}"); a.grid(alpha=.3)
        for (s0, s1, _) in P.VENTANAS.values():
            a.axvspan(s0, s1, color="0.9", zorder=0)
    ax[0].legend(fontsize=8); ax[-1].set_xlabel("Estación s [m] (zonas grises = curvas de R constante)")
    ax[0].set_title(f"{nombre}: transferencia de carga")
    guardar(f"{nombre}_LTR.png")

    # --- balanceo
    fig, ax = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)
    for s, (d, _) in por_s.items():
        m = d["s"] >= 0
        ax[0].plot(d["s"][m], d["balanceo_tractor_deg"][m], color=COL[s], lw=1.5, label=DESC[s])
        ax[1].plot(d["s"][m], d["balanceo_remolque_deg"][m], color=COL[s], lw=1.5, label=DESC[s])
        for i in range(4):
            c = f"balanceo_susp_eje{i+1}_deg"
            if c in d:
                ax[2].plot(d["s"][m], d[c][m], color=COL[s], lw=1.2, ls=["-", "--", "-.", ":"][i],
                           label=f"{s} eje {i+1}")
    ax[0].set_ylabel("Tractor [°]\n(relativo a la vía)"); ax[1].set_ylabel("Remolque [°]\n(relativo a la vía)")
    ax[2].set_ylabel("Suspensión [°]\n(chasis − eje)"); ax[2].legend(fontsize=7, ncol=4)
    ax[2].set_xlabel("Estación s [m]"); ax[0].legend(fontsize=8); ax[0].set_title(f"{nombre}: balanceo (+ = hacia el exterior)")
    for a in ax: a.grid(alpha=.3)
    guardar(f"{nombre}_balanceo.png")

    # --- cargas por rueda
    fig, ax = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    for s, (d, _) in por_s.items():
        m = d["s"] >= 0
        for i in range(4):
            ax[i].plot(d["s"][m], d[f"FzL{i+1}"][m] / 1e3, color=COL[s], lw=1.5, label=f"{DESC[s]} ext. (L)")
            ax[i].plot(d["s"][m], d[f"FzR{i+1}"][m] / 1e3, color=COL[s], lw=1.5, ls="--", label=f"{DESC[s]} int. (R)")
    for i, a in enumerate(ax):
        a.set_ylabel(f"Fz eje {i+1} [kN]"); a.grid(alpha=.3)
    ax[0].legend(fontsize=7, ncol=2); ax[-1].set_xlabel("Estación s [m]")
    ax[0].set_title(f"{nombre}: carga vertical por lado de cada eje (cf. Fig. 7a del p22)")
    guardar(f"{nombre}_cargas.png")

    # --- deslizamiento, fuerza lateral
    fig, ax = plt.subplots(4, 2, figsize=(13, 10), sharex=True)
    for s, (d, _) in por_s.items():
        m = d["s"] >= 0
        for i in range(4):
            ax[i, 0].plot(d["s"][m], d[f"alfa{i+1}_deg"][m], color=COL[s], lw=1.3, label=DESC[s])
            ax[i, 1].plot(d["s"][m], (d[f"FyL{i+1}"][m] + d[f"FyR{i+1}"][m]) / 1e3, color=COL[s], lw=1.3, label=DESC[s])
    for i in range(4):
        ax[i, 0].set_ylabel(f"α eje {i+1} [°]"); ax[i, 1].set_ylabel(f"Fy eje {i+1} [kN]")
        ax[i, 0].grid(alpha=.3); ax[i, 1].grid(alpha=.3)
    ax[0, 0].legend(fontsize=8); ax[-1, 0].set_xlabel("Estación s [m]"); ax[-1, 1].set_xlabel("Estación s [m]")
    ax[0, 0].set_title("Ángulo de deslizamiento medio del eje (TMeasy)")
    ax[0, 1].set_title("Fuerza lateral del eje (+ = hacia el centro)")
    guardar(f"{nombre}_deslizamiento_fuerza_lateral.png")

    # --- dirección, articulación, a_y, desviación
    fig, ax = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
    for s, (d, _) in por_s.items():
        m = d["s"] >= 0
        ax[0].plot(d["s"][m], d["direccion_deg"][m], color=COL[s], lw=1.5, label=DESC[s])
        ax[1].plot(d["s"][m], d["articulacion_deg"][m], color=COL[s], lw=1.5)
        ax[2].plot(d["s"][m], d["ay_geom"][m] / P.G, color=COL[s], lw=1.5, label=f"{DESC[s]}: v²/R")
        ax[2].plot(d["s"][m], d["ay_vr"][m] / P.G, color=COL[s], lw=.8, ls=":", label=f"{DESC[s]}: v·r (guiñada medida)")
        ax[3].plot(d["s"][m], d["desv_lat"][m], color=COL[s], lw=1.5)
    ax[0].set_ylabel("Dirección ruedas [°]"); ax[1].set_ylabel("Articulación Γ [°]")
    ax[2].set_ylabel("a_y [g]"); ax[3].set_ylabel("Desv. lateral [m]\n(+ = exterior)")
    ax[0].legend(fontsize=8); ax[2].legend(fontsize=7)
    ax[-1].set_xlabel("Estación s [m]"); ax[0].set_title(f"{nombre}: dirección, articulación, aceleración lateral y trayectoria")
    for a in ax: a.grid(alpha=.3)
    guardar(f"{nombre}_direccion_articulacion_ay.png")

    # --- tablas por corrida
    for s, (d, res) in por_s.items():
        # cada 5 m
        m = d["s"] >= 0
        if m.sum() > 2:
            ss = np.arange(0, d["s"][m].max(), 5.0)
            o = np.argsort(d["s"][m])
            filas5 = []
            for sv in ss:
                fila = {"s_m": sv}
                for c in d:
                    if c in ("s",):
                        continue
                    fila[c] = float(np.interp(sv, d["s"][m][o], d[c][m][o]))
                fila["v_kmh"] = fila["v"] * 3.6
                filas5.append(fila)
            T.escribir_csv(P.resultados_susp(s) / "tablas" / f"{nombre}_cada5m.csv", filas5)
        # ventanas de curva
        for vn, (s0, s1, Rv) in P.VENTANAS.items():
            w = (d["s"] >= s0) & (d["s"] <= s1)
            if w.sum() < 5:
                continue
            fila = {"suspension": s, "corrida": nombre, "curva": vn, "R_m": Rv, "muestras": int(w.sum()),
                    "v_media_kmh": nanmean(d["v"][w]) * 3.6,
                    "ay_geom_media_g": nanmean(d["ay_geom"][w]) / P.G,
                    "ay_vr_media_g": nanmean(d["ay_vr"][w]) / P.G,
                    "LTR_medio": nanmean(d["LTR"][w]), "LTR_p95": p95(d["LTR"][w])}
            for i in range(4):
                fila[f"LTR{i+1}_medio"] = nanmean(d[f"LTR{i+1}"][w])
                fila[f"LTR{i+1}_p95"] = p95(d[f"LTR{i+1}"][w])
            fila["balanceo_tractor_medio_deg"] = nanmean(d["balanceo_tractor_deg"][w])
            fila["balanceo_remolque_medio_deg"] = nanmean(d["balanceo_remolque_deg"][w])
            for i in range(4):
                fila[f"balanceo_susp_eje{i+1}_medio_deg"] = nanmean(d.get(f"balanceo_susp_eje{i+1}_deg", np.full(len(w), np.nan))[w])
            fila["direccion_media_deg"] = nanmean(d["direccion_deg"][w])
            fila["articulacion_media_deg"] = nanmean(d["articulacion_deg"][w])
            fila["guinada_media_deg_s"] = np.degrees(nanmean(d["r_guinada"][w]))
            fila["desv_lat_max_abs_m"] = nanmax(np.abs(d["desv_lat"][w]))
            for i in range(4):
                fila[f"alfa{i+1}_medio_deg"] = nanmean(d[f"alfa{i+1}_deg"][w])
            for i in range(4):
                fila[f"FzL{i+1}_medio_kN"] = nanmean(d[f"FzL{i+1}"][w]) / 1e3
                fila[f"FzR{i+1}_medio_kN"] = nanmean(d[f"FzR{i+1}"][w]) / 1e3
            for i in range(4):
                fila[f"Fy{i+1}_media_kN"] = nanmean(d[f"FyL{i+1}"][w] + d[f"FyR{i+1}"][w]) / 1e3
            filas_est.append(fila)
        # eventos
        ev = {"suspension": s, "corrida": nombre, "fin": res.get("fin", ""), "vuelco": res.get("vuelco", ""),
              "s_final_m": res.get("s_final_m"), "v_final_kmh": res.get("v_final_kmh"),
              "LTR_p95_rampa": p95(d["LTR"][d["s"] >= 0]),
              "balanceo_remolque_max_deg": nanmax(d["balanceo_remolque_deg"][d["s"] >= 0])}
        for i, dsp in enumerate(res.get("despegue_por_eje", [None] * 4)):
            ev[f"despegue_eje{i+1}_s_m"] = None if dsp is None else dsp["s"]
            ev[f"despegue_eje{i+1}_v_kmh"] = None if dsp is None else dsp["v_kmh"]
        filas_ev.append(ev)

T.escribir_csv(TAB / "tabla_rampa_estacionario.csv", filas_est)
T.escribir_csv(TAB / "tabla_eventos.csv", filas_ev)

# ============================================================================ 3. RESÚMENES
filas_val = []
const = [(n, v) for n, v in corridas.items() if n.startswith("rampa_") and "kmh" in n]
if const:
    def vkmh(n):
        return float(n.split("_")[1].replace("kmh", ""))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for s in SUSPS:
        xs, ltr, ltr_i, phi = [], [], [[] for _ in range(4)], []
        for n, por_s in sorted(const, key=lambda x: vkmh(x[0])):
            if s not in por_s:
                continue
            d, res = por_s[s]
            w = (d["s"] >= P.VENTANAS["R60"][0]) & (d["s"] <= P.VENTANAS["R60"][1])
            if w.sum() < 5:
                continue
            xs.append(nanmean(d["v"][w]) * 3.6); ltr.append(nanmean(d["LTR"][w]))
            for i in range(4):
                ltr_i[i].append(p95(d[f"LTR{i+1}"][w]))
            v_ = vkmh(n)
            fila = {"suspension": s, "v_objetivo_kmh": v_, "v_real_media_R60_kmh": xs[-1],
                    "LTR_medio_R60": ltr[-1], "LTR4_p95_R60": ltr_i[3][-1],
                    "LTR_ref11_TruckSim": P.REF11_LTR.get(int(v_), ""),
                    "vuelco": res.get("vuelco", ""), "v_vuelco_ref11_kmh": P.REF11_V_VUELCO_KMH}
            filas_val.append(fila)
        if xs:
            ax[0].plot(xs, ltr, "o-", color=COL[s], lw=2, label=f"{DESC[s]}: LTR total medio")
            ax[0].plot(xs, ltr_i[3], "s--", color=COL[s], lw=1.2, label=f"{DESC[s]}: LTR eje 4 (p95)")
            for i in range(4):
                ax[1].plot(xs, ltr_i[i], marker="o", color=C_EJE[i], ls="-" if s == SUSPS[0] else "--",
                           label=f"Eje {i+1} ({s})")
    ax[0].plot(list(P.REF11_LTR.keys()), list(P.REF11_LTR.values()), "k^", ms=9, label="Ref. [11] TruckSim")
    ax[0].axvline(P.REF11_V_VUELCO_KMH, color="k", ls=":", lw=1, label="Vuelco ref. [11] (54.9 km/h)")
    for a in ax:
        a.axhline(0.8, color="r", ls="--", lw=.8); a.axhline(1.0, color="k", ls="--", lw=.8)
        a.set_xlabel("Velocidad media en la curva de 60 m [km/h]"); a.grid(alpha=.3); a.legend(fontsize=7)
    ax[0].set_ylabel("LTR"); ax[0].set_title("LTR en la curva de 60 m vs velocidad")
    ax[1].set_ylabel("LTR por eje (percentil 95)"); ax[1].set_title("LTR por eje, p95 (— multihoja, - - aire)")
    guardar("resumen_LTR_vs_velocidad_R60.png")
T.escribir_csv(TAB / "tabla_validacion_ref11.csv", filas_val)

if filas_est or mesa:
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for s in SUSPS:
        if s in mesa:
            ax[0].plot(mesa[s]["ay_g"], mesa[s]["balanceo_tractor_deg"], color=COL[s], lw=1, alpha=.6, label=f"{s}: mesa")
            ax[1].plot(mesa[s]["ay_g"], mesa[s]["balanceo_remolque_deg"], color=COL[s], lw=1, alpha=.6, label=f"{s}: mesa")
        f = [x for x in filas_est if x["suspension"] == s]
        if f:
            ax[0].plot([x["ay_geom_media_g"] for x in f], [x["balanceo_tractor_medio_deg"] for x in f], "o", color=COL[s], label=f"{s}: rampa (curvas)")
            ax[1].plot([x["ay_geom_media_g"] for x in f], [x["balanceo_remolque_medio_deg"] for x in f], "o", color=COL[s], label=f"{s}: rampa (curvas)")
    for a, t in zip(ax, ["Tractor", "Semirremolque"]):
        a.set_xlabel("a_y [g]"); a.set_title(f"Balanceo vs aceleración lateral — {t}"); a.grid(alpha=.3); a.legend(fontsize=7)
    ax[0].set_ylabel("Balanceo [°]")
    guardar("resumen_balanceo_vs_ay.png")

ok = T.escribir_excel(TAB / "resultados_p22_suspensiones.xlsx", {
    "parametros": T.leer_csv_dict(TAB / "tabla_parametros_suspension.csv") if (TAB / "tabla_parametros_suspension.csv").exists() else [],
    "masas_cg": T.leer_csv_dict(TAB / "tabla_masas_y_cg.csv") if (TAB / "tabla_masas_y_cg.csv").exists() else [],
    "mesa_basculante": filas_mesa, "mesa_curvas": filas_curvas,
    "rampa_estacionario": filas_est, "eventos": filas_ev, "validacion_ref11": filas_val})
print(f"\nTablas -> {TAB}" + ("  (+ Excel)" if ok else "  (sin Excel: pip install openpyxl)"))
print(f"Figuras -> {FIG}")
