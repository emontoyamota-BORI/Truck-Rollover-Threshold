# =============================================================================
#  lanzar.py — Corre un script de este proyecto en un proceso de Python aparte
#  (el mismo intérprete del entorno "chrono") y muestra su salida en vivo.
#  Cada simulación va en su propio proceso: así Chrono e Irrlicht arrancan
#  limpios en cada corrida y un error en una no tumba a las demás.
# =============================================================================
import subprocess, sys, time
import config as P


def correr(script, *args):
    cmd = [sys.executable, "-u", str(P.CARPETA / script)] + [str(a) for a in args]
    print("\n" + "=" * 78 + f"\n>>> {' '.join(cmd[2:])}\n" + "=" * 78, flush=True)
    t0 = time.time()
    proc = subprocess.Popen(cmd, cwd=str(P.CARPETA), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace", bufsize=1)
    for linea in proc.stdout:
        print(linea, end="", flush=True)
    proc.wait()
    print(f"<<< terminó con código {proc.returncode} en {time.time()-t0:.0f} s", flush=True)
    return proc.returncode
