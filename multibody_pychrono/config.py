# =============================================================================
#  config.py — TODOS los parámetros de la simulación p22 + suspensiones
#
#  Tractocamión 4x2 + semirremolque de 2 ejes del p22 (Yoo, Kim & Lee 2025) en la
#  rampa del accidente de la ref. [11] (Zhang et al. 2024), con DOS suspensiones:
#     "multihoja" : muelles de acero multihoja (eje motriz y tándem del remolque)
#     "aire"      : suspensión neumática de brazo arrastrado (eje motriz y remolque)
#  El eje delantero (muelle de una hoja) es igual en ambos casos.
#
#  Si cambias algo aquí, vuelve a correr  paso1_preparar.py.
#
#  Convención de Chrono::Vehicle: X adelante, Y izquierda, Z arriba. Unidades SI.
#  Fuentes:  [p22] Yoo et al. 2025, Tabla 1        [r11] Zhang et al. 2024, Tablas 1–4
#            [p35] Winkler, Karamihas & Bogard 1992, SAE 922426, tabla B-1 (UMTRI)
#            [p36] Fu & Cebon 2002, IJHVS 9(4)        [sup] supuesto propio
# =============================================================================
from pathlib import Path

G = 9.81
CARPETA = Path(__file__).resolve().parent
DATOS = CARPETA / "datos"            # datos/<suspension>/p22/{vehiculo,remolque,terreno}
RESULTADOS = CARPETA / "resultados"  # resultados/<suspension>/... y resultados/comparacion/...
PERFIL_P22 = DATOS / "perfil_v_p22.csv"   # perfil de velocidad de las Ecs. 1–19 (Colab)


def carpeta_ascii(nombre):
    """Carpeta SIN acentos. Chrono (C++) en Windows no abre ni escribe rutas con
       acentos (Maestría, DINÁMICA, VEHÍCULOS), así que lo que lee o escribe Chrono
       va a C:\\Users\\<usuario>\\<nombre>."""
    for base in (Path.home(), Path("C:/")):
        d = base / nombre
        if d.as_posix().isascii():
            return d
    return Path(nombre).resolve()


DATOS_ASCII = carpeta_ascii("p22_susp_datos")     # copia de datos/ que lee Chrono
VIDEOS_ASCII = carpeta_ascii("p22_susp_videos")   # fotogramas que escribe Chrono

# -----------------------------------------------------------------------------
# 1. SUSPENSIONES A COMPARAR                                       [p35][p36]
# -----------------------------------------------------------------------------
# K_PHI: rigidez al balanceo POR EJE [kN·m/rad], centro de los rangos medidos por
# UMTRI (p35; 1 kN·m/deg = 57.3 kN·m/rad):
#   dirección hoja simple 70–200 | motriz 4-spring 315–625 | motriz aire 150–330
#   remolque 4-spring 460–1030   | remolque aire 690–1090
# FREC: frecuencia de rebote de la masa suspendida, fija la rigidez vertical de
#   cada resorte (acero ≈ 2 Hz, aire ≈ 1.3 Hz cargado) [sup].
# Y_RES: semi-separación lateral de los resortes [m] (2p = 0.98 m en p36) [sup].
# Lo que falta para llegar a K_PHI (K_PHI − 2·k·y²) se agrega como rigidez auxiliar
# (torsión de eje/brazo) con un resorte de torsión entre eje y chasis. En aire es
# > 70 % del total (p36).
# LIMITACIÓN: la altura del centro de balanceo NO cambia entre suspensiones; la
# fija la geometría de la plantilla de eje de Chrono (p35: acero 600–790 mm).
SUSPENSIONES = {
    "multihoja": dict(
        descripcion="Muelles de acero multihoja (4-spring)",
        K_PHI_kNm_rad=[150, 470, 750, 750],
        FREC_HZ=[1.6, 2.0, 2.0, 2.0],
        Y_RES=[0.45, 0.49, 0.49, 0.49],
        ZETA=0.25,
        color="tab:red"),
    "aire": dict(
        descripcion="Neumática de brazo arrastrado",
        K_PHI_kNm_rad=[150, 240, 890, 890],
        FREC_HZ=[1.6, 1.3, 1.3, 1.3],
        Y_RES=[0.45, 0.55, 0.55, 0.55],
        ZETA=0.25,
        color="tab:blue"),
}
LISTA_SUSPENSIONES = ["multihoja", "aire"]

# -----------------------------------------------------------------------------
# 2. CORRIDAS EN LA RAMPA (las ejecuta correr_todo.py)
# -----------------------------------------------------------------------------
# modo: "constante" (velocidad fija) | "perfil_mu_constante" | "perfil_mu_N"
#       (perfil v_final(s) de las Ecs. 1–19 del p22, calculado en Colab)
# video: True = ventana 3D + fotogramas + GIF (más lento); False = sin ventana (rápido)
# mu: fricción del pavimento (0.8 seco [r11]; 0.3 = caso de las Figs. 6–8 del p22)
CASOS_RAMPA = [
    dict(modo="constante", v_kmh=30, video=False, mu=0.8),
    dict(modo="constante", v_kmh=40, video=True,  mu=0.8),   # ref. [11] Tabla 4: LTR 0.467
    dict(modo="constante", v_kmh=44, video=False, mu=0.8),   # para ubicar la velocidad de vuelco
    dict(modo="constante", v_kmh=47, video=False, mu=0.8),
    dict(modo="constante", v_kmh=50, video=True,  mu=0.8),   # ref. [11] Tabla 4: LTR 0.796
    dict(modo="constante", v_kmh=55, video=False, mu=0.8),   # ref. [11]: vuelco a 54.9 km/h
    dict(modo="constante", v_kmh=60, video=False, mu=0.8),
    dict(modo="perfil_mu_constante", video=True, mu=0.3),    # réplica de la Fig. 8a del p22
]

# Valores de la ref. [11] (TruckSim) para validar
REF11_LTR = {30: 0.225, 40: 0.467, 50: 0.796}     # Tabla 4, R = 60 m
REF11_V_VUELCO_KMH = 54.9

# Ventanas de curva de radio constante para promediar (estación del eje delantero)
VENTANAS = {"R125": (170.0, 290.0, 125.0), "R60": (375.0, 430.0, 60.0)}

# -----------------------------------------------------------------------------
# 3. OPCIONES DE SIMULACIÓN
# -----------------------------------------------------------------------------
DT = 1e-3                    # paso de integración [s]
T_MAX = 400.0                # tiempo máximo simulado [s]
S_INICIO = -480.0            # estación inicial del eje delantero [m] (recta de arranque)
S_FIN = 630.0                # estación final [m]
DT_REGISTRO = 0.05           # registrar telemetría cada 0.05 s
ANTICIPO = 10.0              # la velocidad objetivo se lee ANTICIPO m adelante [m]
LOOKAHEAD = 8.0              # anticipación del control de dirección [m]
GANANCIAS_DIRECCION = (0.5, 0.0, 0.0)
GANANCIAS_VELOCIDAD = (0.6, 0.05, 0.0)   # las de p22_chrono (probadas)
T_ASENTAR = 2.0              # [s] al inicio el camión se asienta frenado y sin dirección (como en la
                             # mesa basculante, que nunca falla); después toma el control el conductor
DZ_INICIAL = [0.03, 0.01, 0.06]   # altura inicial sobre el piso [m]; correr_todo reintenta con la
                                  # siguiente si una corrida reporta FALLA AL ARRANCAR
V_REF_INICIAL = 5.0          # [m/s] velocidad objetivo al arrancar (como en p22_chrono: el camión
                             # arranca con acelerador, nunca con freno, mientras se asienta)
ACEL_OBJETIVO = 1.0          # [m/s²] rampa de la velocidad objetivo al acelerar: evita que el
                             # integrador del control se 'cargue' y el camión llegue pasado de velocidad
BALANCEO_VUELCO_DEG = 30.0   # se declara vuelco y se detiene la simulación
LTR_DESPEGUE = 0.98          # LTR de eje a partir del cual se registra "despegue"
DESPEGUE_MIN_S = 0.3         # ...sostenido al menos este tiempo [s] (filtra picos de contacto)

# Video
FPS_FOTOGRAMAS = 10          # fotogramas guardados por segundo SIMULADO
S_GRABAR_DESDE = -40.0       # empieza a grabar en esta estación [m]
VENTANA_PX = (1280, 720)
CAMARA = dict(punto=(-6.0, 0.0, 2.0), distancia=22.0, altura=0.6)   # cámara de persecución
GIF_ANCHO = 640              # ancho del GIF [px] (alto proporcional)
GIF_CADA = 2                 # usa 1 de cada GIF_CADA fotogramas -> GIF a velocidad ×GIF_CADA
COPIAR_FOTOGRAMAS = True     # copia los fotogramas a resultados/<susp>/fotogramas/<corrida>/

# Mesa basculante virtual
MESA_GRADOS_POR_S = 1.0
MESA_ANGULO_MAX_DEG = 35.0

# -----------------------------------------------------------------------------
# 4. MASAS Y REPARTO DE CARGA                                   [p22] Tabla 1
# -----------------------------------------------------------------------------
M_TRACTOR = 6800.0
M_REMOLQUE_VACIO = 10700.0          # [p22]; la [r11] da 7050 kg
M_CARGA = 17200.0
F_EJE = [0.21, 0.30, 0.25, 0.24]    # fracción del peso total por eje
M_NOSUSP = [700.0, 1200.0, 900.0, 900.0]   # masas no suspendidas por eje [kg] [sup]

def masa_total():
    return M_TRACTOR + M_REMOLQUE_VACIO + M_CARGA

# -----------------------------------------------------------------------------
# 5. GEOMETRÍA                                                     [r11][p22]
# -----------------------------------------------------------------------------
VIA = [2.04, 1.82, 1.82, 1.82]
SEP_DUALES = 0.33                   # [sup]
X_EJE_TRACTOR = [0.0, -3.80]        # [r11]
X_QUINTA_RUEDA = -3.20              # [sup]
H_QUINTA_RUEDA = 1.20               # [sup]
DX_PERNO_EJE3 = 9.545               # [r11]
DX_TANDEM = 1.31                    # [r11]
H_CG_TRACTOR_SUSP = 1.10            # [sup]
H_CG_REMOLQUE_SUSP = 2.725          # [r11] CG del remolque cargado
I_TRACTOR = [2283.9, 35402.8, 34802.6]           # [r11] Tabla 1
I_REMOLQUE_VACIO = [9959.7, 171336.0, 179992.0]  # [r11] Tabla 2
CARGA_CAJA = [15.75, 2.40, 3.45]                 # L, ancho, alto de la carga [m]

# -----------------------------------------------------------------------------
# 6. LLANTAS (TMeasy)                                                [r11][sup]
# -----------------------------------------------------------------------------
LLANTA_TRACTOR_DEL = dict(nombre="295_80R22.5_del", ancho=0.295, relacion=0.80, rin_pulg=22.5, LI=152, masa=65.0)
LLANTA_TRACTOR_TRAS = dict(nombre="295_80R22.5_tras", ancho=0.295, relacion=0.80, rin_pulg=22.5, LI=148, masa=65.0)
LLANTA_REMOLQUE = dict(nombre="275_70R22.5", ancho=0.275, relacion=0.70, rin_pulg=22.5, LI=145, masa=55.0)
MU_TMEASY_REF = 0.8
C_SENSIBILIDAD = None       # None = TMeasy por defecto; 0.4 = sensibilidad a la carga tipo Ec. 6 del p22

# -----------------------------------------------------------------------------
# 7. TREN MOTRIZ                                                         [p22]
# -----------------------------------------------------------------------------
P_MAX = 300e3
T_MAX_MOTOR = 38270.0
W_MAX_MOTOR = 300.0

# -----------------------------------------------------------------------------
# 8. RAMPA                                                               [r11]
# -----------------------------------------------------------------------------
# (s_inicio, s_fin, kappa_inicio, kappa_fin, pendiente %, peralte_inicio %, peralte_fin %)
R1, R2 = 125.0, 60.0
TRAMOS = [
    (0.0,     50.0,    0.0,  0.0,  -2.45, 0, 0),
    (50.0,    146.8,   0.0,  1/R1, -2.45, 0, 6),
    (146.8,   308.032, 1/R1, 1/R1, -2.45, 6, 6),
    (308.032, 363.469, 1/R1, 1/R2, -2.45, 6, 7),
    (363.469, 437.268, 1/R2, 1/R2, -2.45, 7, 7),
    (437.268, 543.935, 1/R2, 0.0,   0.893, 7, 0),
    (543.935, 643.935, 0.0,  0.0,   0.0,   0, 0),
]
GIRO = -1                 # -1 = gira a la derecha: lado EXTERIOR = IZQUIERDO (L)
L_ARRANQUE = 500.0        # recta plana para acelerar antes de s = 0 [m]
ANCHO_CARRIL = 14.0
DS_MALLA = 1.0

NOMBRES_EJE = ["p22 eje delantero", "p22 eje motriz", "p22 eje remolque 1", "p22 eje remolque 2"]


def datos_susp(susp):
    return DATOS / susp


def resultados_susp(susp):
    return RESULTADOS / susp


def nombre_corrida(caso):
    if caso["modo"] == "constante":
        base = f"rampa_{caso['v_kmh']:.0f}kmh"
    else:
        base = f"rampa_{caso['modo']}"
    return base + f"_mu{caso.get('mu', 0.8):.1f}".replace(".", "")
