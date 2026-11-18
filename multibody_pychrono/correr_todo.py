# =============================================================================
#  correr_todo.py — Corre TODO el estudio de suspensiones de principio a fin:
#    1. paso1_preparar.py        datos de rampa y vehículo de las dos suspensiones
#    2. paso2_mesa_basculante.py SRT y despegue por eje (multihoja y aire)
#    3. paso3_rampa.py           cada caso de config.CASOS_RAMPA con cada suspensión
#                                (los casos con video abren la ventana 3D y generan
#                                 fotogramas + GIF/MP4)
#    4. video.py                 GIF comparativos lado a lado
#    5. paso4_resultados.py      gráficas y tablas comparativas
#  Uso: abre Spyder DESDE el entorno "chrono" y corre este archivo (F5),
#       o en el Anaconda Prompt:  conda activate chrono  ->  python correr_todo.py
#  Para repetir sólo una parte, cambia las banderas de abajo.
#  Tiempo aprox.: 3–10 min por corrida con video, 1–3 min sin video (depende de tu PC).
# =============================================================================
import json, os, time
import config as P
import lanzar


def correr_rampa(s, caso):
    """Corre un caso; si reporta FALLA AL ARRANCAR, reintenta con otra altura inicial."""
    nombre = P.nombre_corrida(caso)
    res = P.resultados_susp(s) / "telemetria" / f"{nombre}_resumen.json"
    for k, dz in enumerate(P.DZ_INICIAL):
        os.environ["P22_DZ_INICIAL"] = str(dz)
        rc = lanzar.correr("paso3_rampa.py", s, caso["modo"], caso.get("v_kmh", 0), caso.get("mu", 0.8),
                           1 if caso.get("video") else 0)
        fin = ""
        if res.exists():
            fin = json.load(open(res, encoding="utf-8")).get("fin", "")
        if rc == 0 and not fin.startswith("FALLA"):
            return True
        print(f"*** {s} {nombre}: '{fin or rc}'. Reintento {k+1} de {len(P.DZ_INICIAL)-1}...")
    return False

PREPARAR = True
MESA = False            # la mesa basculante ya está corrida y es válida
RAMPA = True
COMPARATIVOS = True
RESULTADOS = True
SOLO_SUSPENSIONES = P.LISTA_SUSPENSIONES     # p. ej. ["aire"] para correr sólo una

t0 = time.time()
fallas = []
if PREPARAR and lanzar.correr("paso1_preparar.py") != 0:
    raise SystemExit("Falló paso1_preparar.py; revisa el mensaje de arriba.")
if MESA:
    for s in SOLO_SUSPENSIONES:
        if lanzar.correr("paso2_mesa_basculante.py", s) != 0:
            fallas.append(f"mesa {s}")
if RAMPA:
    for caso in P.CASOS_RAMPA:
        for s in SOLO_SUSPENSIONES:
            if not correr_rampa(s, caso):
                fallas.append(f"rampa {s} {P.nombre_corrida(caso)}")
if COMPARATIVOS:
    lanzar.correr("video.py")
if RESULTADOS:
    lanzar.correr("paso4_resultados.py")

print("\n" + "#" * 78)
print(f"Terminado en {(time.time()-t0)/60:.1f} min.")
print("Fallas: " + (", ".join(fallas) if fallas else "ninguna"))
print(f"Resultados en: {P.RESULTADOS}")
print("  comparacion/figuras, comparacion/tablas, comparacion/videos (GIF lado a lado)")
print("  <suspension>/videos (GIF y MP4), <suspension>/fotogramas, <suspension>/telemetria")
