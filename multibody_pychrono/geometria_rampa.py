# =============================================================================
#  geometria_rampa.py — Rampa del accidente (ref. [11] del p22) para Chrono:
#     <carpeta>/rampa.obj   malla de triángulos (RigidTerrain)
#     <carpeta>/eje.csv     eje central: s, x, y, z, kappa, pendiente, peralte
#  Tramos (clotoides, curvas de 125 m y 60 m, peralte 6–7 %, pendiente −2.45 %)
#  en config.TRAMOS. No necesita PyChrono.
# =============================================================================
import numpy as np
import config as P


def eje_central(ds=P.DS_MALLA):
    s_fin = P.TRAMOS[-1][1]
    s = np.arange(-P.L_ARRANQUE, s_fin + 1e-9, ds)
    kap = np.zeros_like(s); pend = np.zeros_like(s); per = np.zeros_like(s)
    for (a, b, k0, k1, Gp, e0, e1) in P.TRAMOS:
        m = (s >= a) & (s <= b)
        x = (s[m] - a) / (b - a)
        kap[m] = k0 + (k1 - k0) * x          # clotoide: curvatura lineal en s
        pend[m] = Gp
        per[m] = e0 + (e1 - e0) * x          # peralte interpolado
    alpha = np.arctan(pend / 100.0)
    beta = np.arctan(per / 100.0)
    rumbo = P.GIRO * np.concatenate([[0], np.cumsum(kap[:-1] * ds)])
    x = np.concatenate([[0], np.cumsum(np.cos(rumbo[:-1]) * ds)]) - P.L_ARRANQUE
    y = np.concatenate([[0], np.cumsum(np.sin(rumbo[:-1]) * ds)])
    z = np.concatenate([[0], np.cumsum(np.tan(alpha[:-1]) * ds)])
    z = z - z.min()
    return dict(s=s, x=x, y=y, z=z, kap=kap, alpha=alpha, beta=beta, rumbo=rumbo)


def escribir_obj(E, archivo, n_lat=8):
    """El peralte levanta el lado exterior (izquierdo, +Y, con giro a la derecha)."""
    d = np.linspace(-P.ANCHO_CARRIL / 2, P.ANCHO_CARRIL / 2, n_lat + 1)
    nx, ny = -np.sin(E['rumbo']), np.cos(E['rumbo'])
    V = []
    for i in range(len(E['s'])):
        for dj in d:
            zj = E['z'][i] + (-P.GIRO) * dj * np.tan(E['beta'][i])
            V.append((E['x'][i] + dj * nx[i], E['y'][i] + dj * ny[i], zj))
    V = np.array(V)
    nc = len(d)
    with open(archivo, "w") as f:
        f.write("# Rampa p22 / ref.[11]\n")
        for v in V:
            f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")
        for i in range(len(E['s']) - 1):
            for j in range(nc - 1):
                a = i * nc + j + 1; b = a + 1; c = a + nc; e = c + 1
                f.write(f"f {a} {c} {b}\nf {b} {c} {e}\n")
    return V


def generar(carpeta):
    carpeta.mkdir(parents=True, exist_ok=True)
    E = eje_central()
    V = escribir_obj(E, carpeta / "rampa.obj")
    np.savetxt(carpeta / "eje.csv",
               np.column_stack([E['s'], E['x'], E['y'], E['z'], E['kap'], E['alpha'], E['beta']]),
               delimiter=",", header="s,x,y,z,kappa,alpha_rad,beta_rad", comments="", fmt="%.6f")
    return E, len(V)


def figura(E, archivo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    m = E['s'] >= 0
    ax[0].plot(E['x'], E['y'], color='0.7', label='recta de arranque')
    sc = ax[0].scatter(E['x'][m], E['y'][m], c=np.degrees(E['beta'][m]), s=3, cmap='viridis')
    fig.colorbar(sc, ax=ax[0], label='peralte [°]'); ax[0].set_aspect('equal')
    ax[0].set_title('Rampa de la ref. [11] — planta'); ax[0].set_xlabel('x [m]'); ax[0].set_ylabel('y [m]')
    ax[0].legend(fontsize=8)
    ax[1].plot(E['s'][m], np.where(E['kap'][m] > 0, 1 / np.maximum(E['kap'][m], 1e-9), np.nan))
    ax[1].set_ylim(0, 300); ax[1].set_title('Radio de curva'); ax[1].set_xlabel('s [m]'); ax[1].set_ylabel('R [m]')
    ax[2].plot(E['s'][m], 100 * np.tan(E['beta'][m]), label='peralte [%]')
    ax[2].plot(E['s'][m], 100 * np.tan(E['alpha'][m]), label='pendiente [%]')
    ax[2].set_title('Peralte y pendiente'); ax[2].set_xlabel('s [m]'); ax[2].legend()
    for a in ax: a.grid(alpha=.3)
    plt.tight_layout(); plt.savefig(archivo, dpi=120); plt.close()
