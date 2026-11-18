# =============================================================================
#  paso0_verificar.py — Comprueba el entorno ANTES de simular.
#  Correr en Spyder abierto desde el entorno conda "chrono" (o: python paso0_verificar.py).
# =============================================================================
import sys, importlib

print("Python:", sys.version.split()[0], "|", sys.executable)
ok = True
for mod in ["numpy", "matplotlib", "PIL", "pychrono", "pychrono.vehicle", "pychrono.irrlicht"]:
    try:
        m = importlib.import_module(mod)
        print(f"  OK  {mod} {getattr(m, '__version__', '')}")
    except Exception as e:
        print(f"  --  {mod}: {e}")
        if mod in ("numpy", "matplotlib", "PIL", "pychrono", "pychrono.vehicle"):
            ok = False
for mod in ["imageio", "imageio_ffmpeg", "openpyxl"]:
    try:
        importlib.import_module(mod)
        print(f"  OK  {mod} (opcional)")
    except Exception:
        print(f"  --  {mod} (opcional: " + {"imageio": "MP4", "imageio_ffmpeg": "MP4",
                                            "openpyxl": "tablas en Excel"}[mod] + ")")
if not ok:
    print("\nFalta algo esencial. En el Anaconda Prompt:  conda activate chrono  y abre Spyder desde ahí.")
    print("Si falta Pillow:  conda install -c conda-forge pillow")
else:
    print("\nTODO LISTO. Siguiente: paso1_preparar.py")
print("Opcionales:  pip install imageio imageio-ffmpeg openpyxl")
