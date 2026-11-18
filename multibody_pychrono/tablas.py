# =============================================================================
#  tablas.py — Utilidades para leer/escribir tablas CSV (y Excel si hay openpyxl).
#  Sin pandas, para que funcione en el entorno "chrono" tal cual.
# =============================================================================
import csv
import numpy as np


def escribir_csv(ruta, filas, columnas=None):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    if not filas:
        ruta.write_text("", encoding="utf-8")
        return
    columnas = columnas or list(dict.fromkeys(k for f in filas for k in f))
    with open(ruta, "w", newline="", encoding="utf-8-sig") as f:   # utf-8-sig: Excel abre bien los acentos
        w = csv.DictWriter(f, fieldnames=columnas, extrasaction="ignore")
        w.writeheader()
        for fila in filas:
            w.writerow({k: _fmt(fila.get(k, "")) for k in columnas})


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        if np.isnan(v):
            return ""
        return f"{v:.5g}" if abs(v) < 1e5 else f"{v:.6g}"
    return v


def leer_telemetria(ruta):
    """CSV numérico con encabezado -> dict columna: np.array."""
    with open(ruta, encoding="utf-8") as h:
        cab = h.readline().strip().split(",")
    A = np.loadtxt(ruta, delimiter=",", skiprows=1, ndmin=2)
    return {c: A[:, i] for i, c in enumerate(cab)}


def leer_csv_dict(ruta):
    with open(ruta, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def escribir_excel(ruta, hojas):
    """hojas: dict nombre -> lista de filas (dict). Devuelve True si se pudo."""
    try:
        import openpyxl
    except ImportError:
        return False
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for nombre, filas in hojas.items():
        ws = wb.create_sheet(nombre[:31])
        if not filas:
            continue
        cols = list(dict.fromkeys(k for f in filas for k in f))
        ws.append(cols)
        for f in filas:
            fila = []
            for c in cols:
                v = f.get(c, "")
                if isinstance(v, (np.floating, float)):
                    v = None if np.isnan(v) else float(v)
                elif isinstance(v, np.integer):
                    v = int(v)
                fila.append(v)
            ws.append(fila)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    wb.save(ruta)
    return True
