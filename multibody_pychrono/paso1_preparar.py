# =============================================================================
#  paso1_preparar.py — Genera los datos de AMBAS suspensiones (no necesita PyChrono):
#     datos/<susp>/p22/terreno   malla de la rampa (ref. [11]) y eje central
#     datos/<susp>/p22/vehiculo  JSON del tractor
#     datos/<susp>/p22/remolque  JSON del semirremolque
#     resultados/comparacion/tablas/tabla_parametros_suspension.csv
#     resultados/comparacion/figuras/rampa_geometria.png
#  Correr con:  python paso1_preparar.py
# =============================================================================
import shutil
import config as P
import geometria_rampa as GR
import vehiculo_json as VJ
import tablas as T

fig_dir = P.RESULTADOS / "comparacion" / "figuras"
tab_dir = P.RESULTADOS / "comparacion" / "tablas"
fig_dir.mkdir(parents=True, exist_ok=True); tab_dir.mkdir(parents=True, exist_ok=True)

if not P.PERFIL_P22.exists():
    viejo = P.CARPETA.parent.parent / "p22_chrono" / "data" / "p22" / "perfil_v_p22.csv"
    if viejo.exists():
        P.PERFIL_P22.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(viejo, P.PERFIL_P22)
        print(f"Perfil de velocidad del p22 copiado de {viejo}")
    else:
        print("AVISO: falta datos/perfil_v_p22.csv (sólo se usa en los casos 'perfil_*').")

filas = []
for susp in P.LISTA_SUSPENSIONES:
    E, nv = GR.generar(P.datos_susp(susp) / "p22" / "terreno")
    r = VJ.generar(susp)
    print(f"\n=== {susp}: {r['descripcion']} ===")
    print(f"Rampa: {nv} vértices | masa total {r['masa_total_kg']/1e3:.1f} t | perno rey {r['perno_rey_kg']/1e3:.2f} t | "
          f"h_cg combinación {r['h_cg_combinacion_m']:.2f} m")
    for i in range(4):
        print(f"  Eje {i+1}: k = {r['k_resorte_N_m'][i]/1e3:5.0f} kN/m, c = {r['c_amort_Ns_m'][i]/1e3:4.1f} kN·s/m, "
              f"y = ±{r['semi_sep_resortes_m'][i]:.2f} m | Kφ resortes {r['K_balanceo_resortes_Nm_rad'][i]/1e3:4.0f} "
              f"+ auxiliar {r['K_balanceo_auxiliar_Nm_rad'][i]/1e3:4.0f} = {r['K_balanceo_total_Nm_rad'][i]/1e3:4.0f} kN·m/rad")
        filas.append({"suspension": susp, "eje": i + 1,
                      "carga_estatica_kg": round(r["cargas_eje_kg"][i], 1),
                      "masa_no_suspendida_kg": r["masa_no_suspendida_kg"][i],
                      "via_m": r["via_m"][i],
                      "frecuencia_Hz": P.SUSPENSIONES[susp]["FREC_HZ"][i],
                      "k_resorte_kN_m": round(r["k_resorte_N_m"][i] / 1e3, 1),
                      "c_amort_kNs_m": round(r["c_amort_Ns_m"][i] / 1e3, 2),
                      "semi_sep_resortes_m": r["semi_sep_resortes_m"][i],
                      "Kphi_resortes_kNm_rad": round(r["K_balanceo_resortes_Nm_rad"][i] / 1e3, 1),
                      "Kphi_auxiliar_kNm_rad": round(r["K_balanceo_auxiliar_Nm_rad"][i] / 1e3, 1),
                      "Kphi_total_kNm_rad": round(r["K_balanceo_total_Nm_rad"][i] / 1e3, 1),
                      "fraccion_auxiliar": round(r["K_balanceo_auxiliar_Nm_rad"][i] / max(r["K_balanceo_total_Nm_rad"][i], 1), 3),
                      "Fz_estatica_por_llanta_kN": round(r["fz_estatica_por_llanta_N"][i] / 1e3, 2)})
    for a in r["avisos"]:
        print("  AVISO:", a)
    # masas y CG (iguales en ambas)
    if susp == P.LISTA_SUSPENSIONES[0]:
        T.escribir_csv(tab_dir / "tabla_masas_y_cg.csv", [
            {"concepto": "masa total [kg]", "valor": r["masa_total_kg"]},
            {"concepto": "carga perno rey [kg]", "valor": round(r["perno_rey_kg"], 1)},
            {"concepto": "tractor masa suspendida [kg]", "valor": round(r["tractor"]["masa_susp_kg"], 1)},
            {"concepto": "tractor h_cg suspendido [m]", "valor": r["tractor"]["h_cg_susp_m"]},
            {"concepto": "tractor x_cg suspendido desde eje 1 [m]", "valor": round(r["tractor"]["x_cg_susp_m"], 3)},
            {"concepto": "remolque masa suspendida [kg]", "valor": round(r["remolque"]["masa_susp_kg"], 1)},
            {"concepto": "remolque h_cg suspendido [m]", "valor": r["remolque"]["h_cg_susp_m"]},
            {"concepto": "remolque x_cg suspendido desde perno [m]", "valor": round(r["remolque"]["x_cg_susp_m"], 3)},
            {"concepto": "h_cg combinación [m]", "valor": round(r["h_cg_combinacion_m"], 3)},
            {"concepto": "radio cargado llanta delantera [m]", "valor": round(r["radio_cargado_del_m"], 4)},
            {"concepto": "radio cargado llanta remolque [m]", "valor": round(r["radio_cargado_rem_m"], 4)},
            {"concepto": "quinta rueda x [m] / h [m]", "valor": f"{r['x_quinta_rueda_m']} / {r['h_quinta_rueda_m']}"},
        ])
GR.figura(E, fig_dir / "rampa_geometria.png")
T.escribir_csv(tab_dir / "tabla_parametros_suspension.csv", filas)
print(f"\nTablas -> {tab_dir}\nSiguiente: paso2_mesa_basculante.py (o correr_todo.py)")
