# =============================================================================
#  vehiculo_json.py — Genera los JSON de Chrono::Vehicle del tractocamión 4x2 +
#  semirremolque de 2 ejes del p22, para UNA suspensión ("multihoja" o "aire").
#
#  Cálculo (equilibrio estático y supuestos documentados en config.py):
#    - carga en el perno rey para que cada eje tenga el % de peso de la Tabla 1
#    - CG longitudinal de las masas suspendidas de tractor y remolque
#    - rigidez vertical k de cada resorte con la frecuencia de rebote; amortiguamiento ζ
#    - rigidez al balanceo de los resortes 2·k·y² y la auxiliar que falta para K_PHI
#    - parámetros TMeasy de las llantas (receta GuessTruck80Par de Chrono 10)
#  Escribe en datos/<susp>/p22/vehiculo, datos/<susp>/p22/remolque y resumen_vehiculo.json.
#  No necesita PyChrono. Plantillas base: ejemplo "truck" de Chrono (BSD).
# =============================================================================
import json, math
import config as P

G = P.G
REL_V, REL_R = "p22/vehiculo/", "p22/remolque/"
M_RIN, M_MAZA = 40.0, 60.0          # [sup]
L_DIS = 0.495                        # longitud de diseño del resorte (plantilla truck) [m]


def escribir(carpeta, nombre, d):
    carpeta.mkdir(parents=True, exist_ok=True)
    with open(carpeta / nombre, "w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)


def carga_por_indice(LI):
    tabla = {140: 2500, 141: 2575, 142: 2650, 143: 2725, 144: 2800, 145: 2900, 146: 3000, 147: 3075,
             148: 3150, 149: 3250, 150: 3350, 151: 3450, 152: 3550, 153: 3650, 154: 3750, 155: 3875,
             156: 4000, 157: 4125, 158: 4250, 159: 4375, 160: 4500}
    return tabla[LI] * G


def tmeasy_json(L, c_sens=None, mu_ref=P.MU_TMEASY_REF, f_rr=0.008, zeta=0.5):
    secth = L['ancho'] * L['relacion']
    rr = L['rin_pulg'] * 0.0254 / 2
    R = rr + secth
    Fli = carga_por_indice(L['LI'])
    pn = 0.5 * Fli
    CZ = Fli / (0.16 * secth)
    DZ = 2 * zeta * math.sqrt(CZ * L['masa'])
    k = {'fxm': (0.89965, 0.77751), 'fxs': (0.46183, 0.42349),
         'fym': (0.75876, 0.72628), 'fys': (0.68276, 0.65319)}
    if c_sens is not None:
        k = {n: (a, a * 2 ** (c_sens - 1)) for n, (a, b) in k.items()}
    par = {
        "Vertical": {"Nominal Vertical Force [N]": pn, "Maximum Vertical Force [N]": 3.5 * pn,
                     "Vertical Tire Stiffness [N/m]": CZ, "Tire Vertical Damping [Ns/m]": DZ},
        "Longitudinal": {"Initial Slopes dFx/dsx [N]": [17.7764 * pn, 14.5301 * 2 * pn],
                         "Maximum Fx Load [N]": [k['fxm'][0] * pn, k['fxm'][1] * 2 * pn],
                         "Sliding Fx Load [N]": [k['fxs'][0] * pn, k['fxs'][1] * 2 * pn],
                         "Slip sx at Maximum Fx": [0.10811, 0.12389],
                         "Slip sx where sliding begins": [0.66667, 0.66667]},
        "Lateral": {"Initial Slopes dFy/dsy [N]": [7.4013 * pn, 6.8505 * 2 * pn],
                    "Maximum Fy Load [N]": [k['fym'][0] * pn, k['fym'][1] * 2 * pn],
                    "Sliding Fy Load [N]": [k['fys'][0] * pn, k['fys'][1] * 2 * pn],
                    "Slip sy at Maximum Fy": [0.33167, 0.33216],
                    "Slip sy where sliding begins": [1.0296, 1.0296]},
        "Aligning": {"Normalized Trail at Zero Slip sy": [0.178, 0.19],
                     "Slip sy where Trail Changes Sign": [0.36484, 0.36538],
                     "Slip sy where Trail Tends to Zero": [1.0296, 1.0296]},
    }
    m, w = L['masa'], L['ancho']
    Iyy = m * (R**2 + rr**2) / 2
    Ixx = m * (3 * (R**2 + rr**2) + w**2) / 12
    J = {"Name": f"TMeasy {L['nombre']}", "Type": "Tire", "Template": "TMeasyTire",
         "Design": {"Unloaded Radius [m]": R, "Mass [kg]": m, "Inertia [kg.m2]": [Ixx, Iyy, Ixx],
                    "Width [m]": w, "Rim Radius [m]": rr, "Roundness of Cross Section": 0.1},
         "Coefficient of Friction": mu_ref, "Rolling Resistance Coefficient": f_rr,
         "Parameters": par}
    return J, R, CZ, pn


def generar(susp):
    S = P.SUSPENSIONES[susp]
    raiz = P.datos_susp(susp) / "p22"
    OUT_V, OUT_R = raiz / "vehiculo", raiz / "remolque"

    # --- equilibrio estático
    m_tot = P.masa_total()
    N = [f * m_tot for f in P.F_EJE]
    M_rem = P.M_REMOLQUE_VACIO + P.M_CARGA
    K = M_rem - (N[2] + N[3])                       # carga en el perno rey [kg]
    x1, x2 = P.X_EJE_TRACTOR
    xfw = P.X_QUINTA_RUEDA
    x3, x4 = -P.DX_PERNO_EJE3, -P.DX_PERNO_EJE3 - P.DX_TANDEM
    m_t = N[0] + N[1] - K
    x_t = (N[1] * x2 - K * xfw) / m_t
    ms_t = m_t - P.M_NOSUSP[0] - P.M_NOSUSP[1]
    xs_t = (m_t * x_t - P.M_NOSUSP[1] * x2) / ms_t
    x_r = (N[2] * x3 + N[3] * x4) / M_rem
    ms_r = M_rem - P.M_NOSUSP[2] - P.M_NOSUSP[3]
    xs_r = (M_rem * x_r - P.M_NOSUSP[2] * x3 - P.M_NOSUSP[3] * x4) / ms_r

    J_del, R_del, CZ_del, pn_del = tmeasy_json(P.LLANTA_TRACTOR_DEL, P.C_SENSIBILIDAD)
    J_tra, R_tra, CZ_tra, pn_tra = tmeasy_json(P.LLANTA_TRACTOR_TRAS, P.C_SENSIBILIDAD)
    J_rem, R_rem, CZ_rem, pn_rem = tmeasy_json(P.LLANTA_REMOLQUE, P.C_SENSIBILIDAD)
    r_carg_del = R_del - (N[0] * G / 2) / CZ_del
    r_carg_rem = R_rem - (N[2] * G / 4) / CZ_rem
    z0 = r_carg_del
    z_eje_rem = r_carg_rem - z0
    h_t, h_r = P.H_CG_TRACTOR_SUSP, P.H_CG_REMOLQUE_SUSP
    radios = [r_carg_del, r_carg_del, r_carg_rem, r_carg_rem]
    h_comb = (ms_t * h_t + ms_r * h_r + sum(m * r for m, r in zip(P.M_NOSUSP, radios))) / m_tot
    L_, W_, H_ = P.CARGA_CAJA
    I_carga = [P.M_CARGA / 12 * (W_**2 + H_**2), P.M_CARGA / 12 * (L_**2 + H_**2), P.M_CARGA / 12 * (L_**2 + W_**2)]
    I_rem = [a + b for a, b in zip(P.I_REMOLQUE_VACIO, I_carga)]

    # --- suspensión
    F_lado = [(N[i] - P.M_NOSUSP[i]) * G / 2 for i in range(4)]
    k_res = [(2 * math.pi * S["FREC_HZ"][i]) ** 2 * F_lado[i] / G for i in range(4)]
    c_amo = [2 * S["ZETA"] * math.sqrt(k_res[i] * F_lado[i] / G) for i in range(4)]
    y_res = list(S["Y_RES"])
    K_eje = [k * 1e3 for k in S["K_PHI_kNm_rad"]]
    K_res = [2 * k_res[i] * y_res[i] ** 2 for i in range(4)]
    K_aux = [max(K_eje[i] - K_res[i], 0.0) for i in range(4)]
    K_real = [K_res[i] + K_aux[i] for i in range(4)]
    avisos = [f"eje {i+1}: los resortes solos ({K_res[i]/1e3:.0f} kN·m/rad) superan K_PHI" for i in range(4) if K_res[i] > K_eje[i]]
    # cargas por llanta (duales en ejes 2–4)
    fz_llanta = [N[0] * G / 2, N[1] * G / 4, N[2] * G / 4, N[3] * G / 4]

    def resorte(i):
        y = y_res[i]
        return ({"Location Chassis": [0.0, y, 0.585], "Location Axle": [0.0, y, 0.09],
                 "Free Length": L_DIS, "Spring Coefficient": k_res[i], "Preload": F_lado[i],
                 "Minimum Length": L_DIS - 0.16, "Maximum Length": L_DIS + 0.14},
                {"Location Chassis": [0.05, y, 0.585], "Location Axle": [0.05, y, 0.09],
                 "Damping Coefficient": c_amo[i]})

    def eje_trasero_json(i, track, n_llantas, m_llanta):
        m_tubo = P.M_NOSUSP[i] - 2 * M_MAZA - n_llantas * (M_RIN + m_llanta)
        sp, sh = resorte(i)
        return {"Name": P.NOMBRES_EJE[i], "Type": "Suspension", "Template": "SolidThreeLinkAxle",
                "Spindle": {"Mass": M_MAZA, "COM": [0, track / 2, 0], "Inertia": [0.27, 0.48, 0.27], "Radius": 0.10, "Width": 0.06},
                "Axle Tube": {"Mass": m_tubo, "COM": [0, 0, 0], "Inertia": [m_tubo * track**2 / 12, 5.0, m_tubo * track**2 / 12], "Radius": 0.06},
                "Triangular Link": {"Mass": 50.0, "Inertia": [0.2, 0.2, 0.2], "Location Chassis": [-0.762, 0.570, 0.100], "Location Axle": [0.0, 0.0, 0.260]},
                "Longitudinal Link": {"Mass": 25.0, "Inertia": [0.5, 0.1, 0.1], "Location Chassis": [1.139, 0.550, 0.100], "Location Axle": [0.115, 0.700, -0.090]},
                "Spring": sp, "Shock": sh, "Axle": {"Inertia": 0.4}}

    def eje_delantero_json():
        r = (P.VIA[0] / 2) / 1.185
        m_kn = 60.0
        m_tubo = P.M_NOSUSP[0] - 2 * M_MAZA - 2 * m_kn - 2 * (M_RIN + P.LLANTA_TRACTOR_DEL['masa'])
        sp, sh = resorte(0)
        return {"Name": P.NOMBRES_EJE[0], "Type": "Suspension", "Template": "SolidBellcrankThreeLinkAxle",
                "Spindle": {"Mass": M_MAZA, "COM": [0, P.VIA[0] / 2, 0], "Inertia": [0.27, 0.48, 0.27], "Radius": 0.10, "Width": 0.06},
                "Axle Tube": {"Mass": m_tubo, "COM": [0, 0, 0], "Inertia": [21.23, 8.12, 21.31], "Radius": 0.0476},
                "Knuckle": {"Mass": m_kn, "COM": [0, 1.087 * r, 0], "Inertia": [1.0, 1.6, 1.0],
                            "Location Lower": [0, 0.853748866 * r, -0.10], "Location Upper": [0, 0.836251134 * r, 0.10],
                            "Location Tierod": [-0.236, 0.800 * r, 0.153]},
                "Bellcrank": {"Mass": 24.4, "Inertia": [0.05, 0.29, 0.30], "Location Axle": [-0.023, 0.0, 0.250],
                              "Location Tierod": [-0.273, 0.192, 0.153]},
                "Draglink": {"Mass": 10.3, "Inertia": [0.29, 0.67, 0.95], "Location Steering": [0.741, -0.067, 0.089],
                             "Location Bellcrank": [0.045, 0.406, 0.153]},
                "Triangular Link": {"Mass": 50.0, "Inertia": [0.2, 0.2, 0.2], "Location Chassis": [0.871, 0.570, 0.165], "Location Axle": [0.125, 0.0, 0.260]},
                "Longitudinal Link": {"Mass": 25.0, "Inertia": [0.05, 0.1, 0.1], "Location Chassis": [-1.138, 0.420, 0.115], "Location Axle": [-0.115, 0.640 * r, -0.090]},
                "Tierod": {"Mass": 10.0, "Inertia": [0.05, 0.1, 0.105]},
                "Spring": sp, "Shock": sh, "Axle": {"Inertia": 0.4}}

    # --- llantas, rines, frenos
    escribir(OUT_V, "llanta_del.json", J_del)
    escribir(OUT_V, "llanta_tras.json", J_tra)
    escribir(OUT_R, "llanta_rem.json", J_rem)
    rin = {"Name": "rin 22.5", "Type": "Wheel", "Template": "Wheel", "Mass": M_RIN, "Inertia": [0.6, 1.0, 0.6],
           "Visualization": {"Radius": 0.2858, "Width": 0.23}}
    escribir(OUT_V, "rin.json", rin); escribir(OUT_R, "rin.json", rin)
    freno = {"Name": "freno", "Type": "Brake", "Template": "BrakeSimple", "Maximum Torque": 15000}
    escribir(OUT_V, "freno.json", freno); escribir(OUT_R, "freno.json", freno)

    # --- tractor
    escribir(OUT_V, "chasis.json", {
        "Name": "p22 tractor chasis", "Type": "Chassis", "Template": "RigidChassis",
        "Components": [{"Centroidal Frame": {"Location": [xs_t, 0, h_t - z0], "Orientation": [1, 0, 0, 0]},
                        "Mass": ms_t, "Moments of Inertia": P.I_TRACTOR, "Products of Inertia": [0, 0, 0], "Void": False}],
        "Driver Position": {"Location": [0.8, 0.5, 1.6], "Orientation": [1, 0, 0, 0]},
        "Rear Connector Location": [xfw, 0, P.H_QUINTA_RUEDA - z0],
        "Visualization": {"Primitives": [
            {"Type": "BOX", "Location": [0.6, 0, 1.4], "Orientation": [1, 0, 0, 0], "Dimensions": [2.2, 2.45, 2.6]},
            {"Type": "BOX", "Location": [-1.9, 0, 0.45], "Orientation": [1, 0, 0, 0], "Dimensions": [5.2, 0.9, 0.3]}]}})
    escribir(OUT_V, "eje_delantero.json", eje_delantero_json())
    escribir(OUT_V, "eje_motriz.json", eje_trasero_json(1, P.VIA[1], 4, P.LLANTA_TRACTOR_TRAS['masa']))
    escribir(OUT_V, "direccion.json", {
        "Name": "p22 direccion", "Type": "Steering", "Template": "RotaryArm",
        "Pitman Arm": {"Mass": 1.605, "Inertia": [0.00638, 0.00756, 0.0015], "Inertia Products": [0, 0, 0], "Radius": 0.02,
                       "Axis of Rotation": [0, 0, 1], "Point of Rotation": [0.850, -0.028, 0.166],
                       "Point to Draglink": [0.741, -0.217, 0.089], "Maximum Angle (deg)": 40}})
    escribir(OUT_V, "transmision_eje.json", {"Name": "p22 driveline", "Type": "Driveline",
                                             "Template": "SimpleDrivelineXWD", "Differential Max Bias": 2.0})
    escribir(OUT_V, "motor.json", {"Name": "p22 motor", "Type": "Engine", "Template": "EngineSimple",
                                   "Maximum Engine Torque": P.T_MAX_MOTOR, "Maximum Engine Power": P.P_MAX,
                                   "Maximum Engine Speed": P.W_MAX_MOTOR})
    escribir(OUT_V, "caja.json", {"Name": "p22 caja", "Type": "Transmission", "Template": "AutomaticTransmissionSimpleMap",
        "Gear Box": {"Reverse Gear Ratio": -0.151,
                     "Forward Gear Ratios": [0.1708, 0.2791, 0.4218, 0.6223, 1.0173, 1.5361],
                     "Shift Points Map RPM": [[1000, 2226], [1000, 2225], [1000, 2210], [1000, 2226], [1000, 2225], [1000, 2700]]}})
    escribir(OUT_V, "tractor.json", {
        "Name": "p22 tractor 4x2", "Type": "Vehicle", "Template": "WheeledVehicle",
        "Chassis": {"Input File": REL_V + "chasis.json"},
        "Axles": [
            {"Suspension Input File": REL_V + "eje_delantero.json", "Suspension Location": [x1, 0, 0], "Steering Index": 0,
             "Left Wheel Input File": REL_V + "rin.json", "Right Wheel Input File": REL_V + "rin.json",
             "Left Brake Input File": REL_V + "freno.json", "Right Brake Input File": REL_V + "freno.json"},
            {"Suspension Input File": REL_V + "eje_motriz.json", "Suspension Location": [x2, 0, 0],
             "Wheel Separation": P.SEP_DUALES,
             "Left Inside Wheel Input File": REL_V + "rin.json", "Right Inside Wheel Input File": REL_V + "rin.json",
             "Left Outside Wheel Input File": REL_V + "rin.json", "Right Outside Wheel Input File": REL_V + "rin.json",
             "Left Brake Input File": REL_V + "freno.json", "Right Brake Input File": REL_V + "freno.json"}],
        "Steering Subsystems": [{"Input File": REL_V + "direccion.json", "Location": [0, 0, 0], "Orientation": [1, 0, 0, 0]}],
        "Wheelbase": -x2, "Minimum Turning Radius": 7.5, "Maximum Steering Angle (deg)": 40.0,
        "Driveline": {"Input File": REL_V + "transmision_eje.json", "Suspension Indexes": [1]}})

    # --- remolque
    escribir(OUT_R, "chasis.json", {
        "Name": "p22 remolque chasis", "Type": "ChassisRear", "Template": "RigidChassisRear",
        "Components": [{"Centroidal Frame": {"Location": [xs_r, 0, h_r - z0], "Orientation": [1, 0, 0, 0]},
                        "Mass": ms_r, "Moments of Inertia": I_rem, "Products of Inertia": [0, 0, 0], "Void": False}],
        "Front Connector Location": [0, 0, P.H_QUINTA_RUEDA - z0],
        "Visualization": {"Primitives": [
            {"Type": "BOX", "Location": [-L_ / 2 + 1.0, 0, 1.0 + H_ / 2 - z0], "Orientation": [1, 0, 0, 0],
             "Dimensions": [L_, 2.5, H_]}]}})
    escribir(OUT_R, "quinta_rueda.json", {"Name": "p22 quinta rueda", "Type": "ChassisConnector",
                                         "Template": "ChassisConnectorFifthWheel"})
    escribir(OUT_R, "eje_remolque_1.json", eje_trasero_json(2, P.VIA[2], 4, P.LLANTA_REMOLQUE['masa']))
    escribir(OUT_R, "eje_remolque_2.json", eje_trasero_json(3, P.VIA[3], 4, P.LLANTA_REMOLQUE['masa']))
    ejes_r = []
    for xa, f in ((x3, "eje_remolque_1.json"), (x4, "eje_remolque_2.json")):
        ejes_r.append({"Suspension Input File": REL_R + f, "Suspension Location": [xa, 0, z_eje_rem],
                       "Wheel Separation": P.SEP_DUALES,
                       "Left Inside Wheel Input File": REL_R + "rin.json", "Right Inside Wheel Input File": REL_R + "rin.json",
                       "Left Outside Wheel Input File": REL_R + "rin.json", "Right Outside Wheel Input File": REL_R + "rin.json",
                       "Left Brake Input File": REL_R + "freno.json", "Right Brake Input File": REL_R + "freno.json"})
    escribir(OUT_R, "remolque.json", {"Name": "p22 semirremolque 2 ejes", "Type": "Trailer", "Template": "WheeledTrailer",
                                     "Chassis": {"Input File": REL_R + "chasis.json"},
                                     "Connector": {"Input File": REL_R + "quinta_rueda.json"},
                                     "Axles": ejes_r})

    res = {"suspension": susp, "descripcion": S["descripcion"],
           "masa_total_kg": m_tot, "cargas_eje_kg": N, "perno_rey_kg": K,
           "tractor": {"masa_kg": m_t, "x_cg_m": x_t, "masa_susp_kg": ms_t, "x_cg_susp_m": xs_t, "h_cg_susp_m": h_t},
           "remolque": {"masa_kg": M_rem, "x_cg_desde_perno_m": x_r, "masa_susp_kg": ms_r, "x_cg_susp_m": xs_r, "h_cg_susp_m": h_r},
           "h_cg_combinacion_m": h_comb, "radio_cargado_del_m": r_carg_del, "radio_cargado_rem_m": r_carg_rem,
           "fz_estatica_por_llanta_N": fz_llanta, "carga_nominal_tmeasy_N": [pn_del, pn_tra, pn_rem, pn_rem],
           "masa_no_suspendida_kg": P.M_NOSUSP, "via_m": P.VIA,
           "k_resorte_N_m": k_res, "c_amort_Ns_m": c_amo, "semi_sep_resortes_m": y_res,
           "K_balanceo_objetivo_Nm_rad": K_eje, "K_balanceo_resortes_Nm_rad": K_res,
           "K_balanceo_auxiliar_Nm_rad": K_aux, "K_balanceo_total_Nm_rad": K_real,
           "x_ejes_tractor_m": [x1, x2], "x_ejes_remolque_desde_perno_m": [x3, x4],
           "x_quinta_rueda_m": xfw, "h_quinta_rueda_m": P.H_QUINTA_RUEDA,
           "nombres_suspension": P.NOMBRES_EJE, "avisos": avisos}
    escribir(raiz, "resumen_vehiculo.json", res)
    return res
