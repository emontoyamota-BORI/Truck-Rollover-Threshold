# =============================================================================
#  chrono_modelo.py — Arma el modelo en PyChrono para UNA suspensión:
#  tractor (WheeledVehicle) + semirremolque (WheeledTrailer) desde los JSON de
#  vehiculo_json.py, llantas TMeasy, tren motriz, resortes de torsión con la
#  rigidez al balanceo auxiliar, terreno (malla de la rampa o plano) y funciones
#  de telemetría (cargas, fuerzas laterales, deslizamiento, balanceo, ángulos).
#  Escrito contra PyChrono 10.0 (las llamadas nuevas van protegidas: si alguna
#  no existe en tu versión, esa columna sale como NaN y la simulación sigue).
# =============================================================================
import json, math, shutil
import pychrono as chrono
import pychrono.vehicle as veh
import config as P

_set_path = getattr(veh, "SetVehicleDataPath", None) or getattr(veh, "SetDataPath")
_get_file = getattr(veh, "GetVehicleDataFile", None) or getattr(veh, "GetDataFile")
NAN = float("nan")


# ----------------------------------------------------------------------------- rutas
def configurar_rutas(susp):
    """Copia datos/<susp> a una carpeta sin acentos y apunta Chrono ahí."""
    origen = P.datos_susp(susp)
    destino = P.DATOS_ASCII / susp
    if origen.as_posix().isascii():
        destino = origen
    else:
        shutil.copytree(origen, destino, dirs_exist_ok=True)
    _set_path(destino.as_posix() + "/")
    return destino


def archivo(rel):
    return _get_file(rel)


def resumen(susp):
    with open(P.datos_susp(susp) / "p22" / "resumen_vehiculo.json", encoding="utf-8") as f:
        return json.load(f)


def ruedas_de(eje):
    w = eje.m_wheels
    return [w[i] for i in range(len(w))]


def ejes_de(sistema_vehicular):
    a = sistema_vehicular.GetAxles()
    return [a[i] for i in range(len(a))]


# ----------------------------------------------------------------------------- modelo
def crear_tractor_remolque(susp, pos_ini, vis=True):
    configurar_rutas(susp)
    vt = chrono.VisualizationType_PRIMITIVES if vis else chrono.VisualizationType_NONE
    tractor = veh.WheeledVehicle(archivo("p22/vehiculo/tractor.json"), chrono.ChContactMethod_NSC)
    tractor.Initialize(chrono.ChCoordsysd(pos_ini, chrono.QUNIT))
    tractor.SetChassisVisualizationType(vt)
    tractor.SetSuspensionVisualizationType(vt)
    tractor.SetSteeringVisualizationType(vt)
    tractor.SetWheelVisualizationType(vt)
    llantas_tractor = ["p22/vehiculo/llanta_del.json", "p22/vehiculo/llanta_tras.json"]
    for i, eje in enumerate(ejes_de(tractor)):
        for rueda in ruedas_de(eje):
            tractor.InitializeTire(veh.TMeasyTire(archivo(llantas_tractor[i])), rueda, vt)
    motor = veh.ReadEngineJSON(archivo("p22/vehiculo/motor.json"))
    caja = veh.ReadTransmissionJSON(archivo("p22/vehiculo/caja.json"))
    tractor.InitializePowertrain(veh.ChPowertrainAssembly(motor, caja))

    remolque = veh.WheeledTrailer(tractor.GetSystem(), archivo("p22/remolque/remolque.json"))
    remolque.Initialize(tractor.GetChassis())
    remolque.SetChassisVisualizationType(vt)
    remolque.SetSuspensionVisualizationType(vt)
    remolque.SetWheelVisualizationType(vt)
    for eje in ejes_de(remolque):
        for rueda in ruedas_de(eje):
            remolque.InitializeTire(veh.TMeasyTire(archivo("p22/remolque/llanta_rem.json")), rueda, vt)

    tractor.GetSystem().SetCollisionSystemType(chrono.ChCollisionSystem.Type_BULLET)
    agregar_torsion_auxiliar(susp, tractor, remolque)
    return tractor, remolque


def _buscar_cuerpo(sistema, nombre):
    cuerpos = sistema.GetBodies()
    for i in range(len(cuerpos)):
        if cuerpos[i].GetName() == nombre:
            return cuerpos[i]
    raise RuntimeError(f"No encontré el cuerpo '{nombre}'")


def agregar_torsion_auxiliar(susp, tractor, remolque):
    """Resorte de torsión (ChLinkRSDA) entre cada tubo de eje y su chasis, alrededor
       del eje X, con la rigidez al balanceo auxiliar K_aux de la suspensión."""
    R = resumen(susp)
    sis = tractor.GetSystem()
    chasis = [tractor.GetChassisBody(), tractor.GetChassisBody(),
              remolque.GetChassis().GetBody(), remolque.GetChassis().GetBody()]
    rot = chrono.QuatFromAngleY(math.pi / 2)        # eje Z del resorte -> eje X del vehículo
    links = []
    for i, nombre in enumerate(R["nombres_suspension"]):
        k = R["K_balanceo_auxiliar_Nm_rad"][i]
        if k <= 0:
            continue
        tubo = _buscar_cuerpo(sis, nombre + "_axleTube")
        rsda = chrono.ChLinkRSDA()
        rsda.Initialize(tubo, chasis[i], chrono.ChFramed(tubo.GetPos(), rot))
        rsda.SetSpringCoefficient(k)
        rsda.SetDampingCoefficient(0.02 * k)
        sis.AddLink(rsda)
        links.append(rsda)
    return links


def crear_terreno(sistema, malla=True, mu=0.8):
    terreno = veh.RigidTerrain(sistema)
    mat = chrono.ChContactMaterialNSC()
    mat.SetFriction(mu)
    mat.SetRestitution(0.01)
    if malla:
        parche = terreno.AddPatch(mat, chrono.CSYSNORM, archivo("p22/terreno/rampa.obj"), True, 0.01)
    else:
        parche = terreno.AddPatch(mat, chrono.CSYSNORM, 400.0, 60.0)
    parche.SetColor(chrono.ChColor(0.55, 0.55, 0.55))
    terreno.Initialize()
    return terreno


# ----------------------------------------------------------------------------- telemetría
def _v(a, b):
    return a.x * b.x + a.y * b.y + a.z * b.z


def ejes_locales(cuerpo):
    q = cuerpo.GetRot()
    return (q.Rotate(chrono.ChVector3d(1, 0, 0)), q.Rotate(chrono.ChVector3d(0, 1, 0)),
            q.Rotate(chrono.ChVector3d(0, 0, 1)))


def angulo_balanceo(cuerpo):
    """Balanceo [rad]: inclinación del eje Y local respecto al horizontal (+ = lado izquierdo arriba)."""
    _, y, _ = ejes_locales(cuerpo)
    return math.asin(max(-1.0, min(1.0, y.z)))


def rumbo(cuerpo):
    x, _, _ = ejes_locales(cuerpo)
    return math.atan2(x.y, x.x)


def angulo_articulacion(tractor, remolque):
    """Γ = rumbo tractor − rumbo remolque [rad], en (−π, π]."""
    d = rumbo(tractor.GetChassisBody()) - rumbo(remolque.GetChassis().GetBody())
    return (d + math.pi) % (2 * math.pi) - math.pi


def velocidad_guinada(cuerpo):
    try:
        return cuerpo.GetAngVelParent().z
    except Exception:
        return NAN


def aceleracion_lateral(cuerpo):
    """Aceleración del CG proyectada en el eje Y del chasis [m/s²] (+ = hacia la izquierda)."""
    try:
        a = cuerpo.GetPosDt2()
    except Exception:
        try:
            a = cuerpo.GetPos_dtdt()
        except Exception:
            return NAN
    _, y, _ = ejes_locales(cuerpo)
    return _v(a, y)


def _cruz(a, b):
    return chrono.ChVector3d(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)


def _unit(a):
    n = math.sqrt(_v(a, a)) or 1.0
    return chrono.ChVector3d(a.x / n, a.y / n, a.z / n)


def normal_via(psi, alfa, beta):
    """Normal de la superficie de la vía a partir del eje de la rampa (rumbo psi,
       pendiente alfa, peralte beta; el peralte levanta el lado exterior = izquierdo
       con GIRO = -1). Se calcula con la geometría y NO con terreno.GetNormal(): esa
       consulta al sistema de colisiones durante la simulación alteró los contactos
       (v2: ruedas delanteras giradas 35° y picos falsos de balanceo)."""
    t = chrono.ChVector3d(math.cos(psi), math.sin(psi), math.tan(alfa))
    l = chrono.ChVector3d(-math.sin(psi), math.cos(psi), (-P.GIRO) * math.tan(beta))
    return _unit(_cruz(t, l))


def fuerzas_por_eje(sistema_vehicular, terreno, chasis, n=None):
    """Para cada eje: (Fn_izq, Fn_der, Fy_izq, Fy_der, alfa_medio, n).
       Fn [N] carga NORMAL a la vía; Fy [N] fuerza lateral EN EL PLANO de la vía
       (+ = izquierda del chasis); alfa [rad] deslizamiento medio de las llantas (TMeasy);
       n = normal de la vía (por defecto Z global: piso plano).
       Orden de ruedas en Chrono: [L, R] o [L_int, R_int, L_ext, R_ext]."""
    xch, _, _ = ejes_locales(chasis)
    n = n if n is not None else chrono.ChVector3d(0, 0, 1)
    ylat = _unit(_cruz(n, xch))                      # lateral en el plano de la vía
    out = []
    for eje in ejes_de(sistema_vehicular):
        fzl = fzr = fyl = fyr = 0.0
        alfas = []
        for k, rueda in enumerate(ruedas_de(eje)):
            llanta = rueda.GetTire()
            f = llanta.ReportTireForce(terreno).force
            fn, fy = _v(f, n), _v(f, ylat)
            if k % 2 == 0:
                fzl += fn; fyl += fy
            else:
                fzr += fn; fyr += fy
            try:
                alfas.append(llanta.GetSlipAngle())
            except Exception:
                pass
        alfa = sum(alfas) / len(alfas) if alfas else NAN
        out.append((fzl, fzr, fyl, fyr, alfa, n))
    return out


def balanceo_rel_via(cuerpo, n):
    """Balanceo [rad] del cuerpo RELATIVO a la superficie de la vía (+ = lado izquierdo arriba)."""
    _, y, _ = ejes_locales(cuerpo)
    return math.asin(max(-1.0, min(1.0, _v(y, n))))


_cache_tubos = {}
def tubos_de_eje(tractor):
    sis = tractor.GetSystem()
    if id(sis) not in _cache_tubos:
        _cache_tubos[id(sis)] = [_buscar_cuerpo(sis, n + "_axleTube") for n in P.NOMBRES_EJE]
    return _cache_tubos[id(sis)]


def balanceo_suspension(tractor, remolque):
    """Balanceo de la masa suspendida RELATIVO a cada eje [rad] (lo que deforma la
       suspensión; es el φ del modelo analítico con Kφ). + = lado izquierdo arriba."""
    ch = [tractor.GetChassisBody()] * 2 + [remolque.GetChassis().GetBody()] * 2
    out = []
    for c, tubo in zip(ch, tubos_de_eje(tractor)):
        _, yc, _ = ejes_locales(c)
        _, yt, _ = ejes_locales(tubo)
        xt, _, _ = ejes_locales(tubo)
        zt = _cruz(xt, yt)
        out.append(math.asin(max(-1.0, min(1.0, _v(yc, zt)))))
    return out


def fuerzas_verticales(sistema_vehicular, terreno):
    out = []
    for eje in ejes_de(sistema_vehicular):
        L = R = 0.0
        for k, rueda in enumerate(ruedas_de(eje)):
            fz = rueda.GetTire().ReportTireForce(terreno).force.z
            if k % 2 == 0: L += fz
            else: R += fz
        out.append((L, R))
    return out


_cache_mazas = {}
def angulo_direccion(tractor):
    """Ángulo de dirección medio de las ruedas delanteras [rad] respecto al chasis,
       a partir de la orientación de las mazas (eje Y de la maza)."""
    sis = tractor.GetSystem()
    clave = id(sis)
    if clave not in _cache_mazas:
        cuerpos = sis.GetBodies()
        mazas = []
        for i in range(len(cuerpos)):
            n = cuerpos[i].GetName()
            if n.startswith(P.NOMBRES_EJE[0]) and "spindle" in n.lower():
                mazas.append(cuerpos[i])
        _cache_mazas[clave] = mazas
    mazas = _cache_mazas[clave]
    if not mazas:
        return NAN
    xc, yc, _ = ejes_locales(tractor.GetChassisBody())
    angs = []
    for m in mazas:
        _, ys, _ = ejes_locales(m)
        angs.append(math.atan2(-_v(ys, xc), _v(ys, yc)) if _v(ys, yc) >= 0 else math.atan2(_v(ys, xc), -_v(ys, yc)))
    return sum(angs) / len(angs)
