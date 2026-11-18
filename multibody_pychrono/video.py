# =============================================================================
#  video.py — Arma GIF (y MP4 si hay imageio-ffmpeg) con los fotogramas que guarda
#  paso3_rampa.py, con un rótulo de suspensión, tiempo, estación, velocidad, LTR
#  y balanceo. También arma el video COMPARATIVO lado a lado (multihoja | aire),
#  emparejando los fotogramas por estación s para que los dos camiones estén en
#  el mismo punto de la rampa.
#  Uso directo:  python video.py                         (todos los comparativos posibles)
#                python video.py rampa_40kmh_mu08         (uno)
# =============================================================================
import sys
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import config as P


def _fuente(tam):
    for nombre in ("arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(nombre, tam)
        except Exception:
            pass
    return ImageFont.load_default()


def _info(carpeta):
    f = Path(carpeta) / "cuadros.csv"
    return np.loadtxt(f, delimiter=",", skiprows=1, ndmin=2) if f.exists() else None


def rotular(img, titulo, fila, color="white", ancho=None):
    img = img.convert("RGB")
    if ancho and img.width != ancho:
        img = img.resize((ancho, int(img.height * ancho / img.width)))
    esc = img.width / 1280
    d = ImageDraw.Draw(img)
    h = int(56 * esc) + 8
    d.rectangle([0, 0, img.width, h], fill=(0, 0, 0))
    f1, f2 = _fuente(max(int(22 * esc), 10)), _fuente(max(int(19 * esc), 9))
    d.text((8, 3), titulo, fill=color, font=f1)
    if fila is not None:
        _, t, s, v, ltr, lat, phi, ay = fila[:8]
        d.text((8, 3 + int(27 * esc)),
               f"t = {t:5.1f} s   s = {s:6.1f} m   v = {v*3.6:5.1f} km/h   a_y = {ay/9.81:4.2f} g   "
               f"LTR = {ltr:4.2f}   balanceo = {phi:4.1f}°", fill="white", font=f2)
    return img.crop((0, 0, img.width - img.width % 2, img.height - img.height % 2))


def _color(nombre):
    import matplotlib.colors as mc
    return tuple(int(255 * c) for c in mc.to_rgb(nombre))


def hacer_video(carpeta, destino_dir, nombre, titulo, color="white"):
    """GIF (siempre) + MP4 (si se puede). Devuelve la lista de archivos creados."""
    carpeta, destino_dir = Path(carpeta), Path(destino_dir)
    destino_dir.mkdir(parents=True, exist_ok=True)
    cuadros = sorted(carpeta.glob("cuadro_*.jpg"))
    if not cuadros:
        return []
    info = _info(carpeta)
    col = _color(color) if color != "white" else "white"
    creados = []
    # MP4 a tiempo real con todos los fotogramas
    try:
        import imageio.v2 as imageio
        mp4 = destino_dir / f"{nombre}.mp4"
        with imageio.get_writer(mp4, fps=P.FPS_FOTOGRAMAS, codec="libx264", quality=8, macro_block_size=8) as w:
            for k, c in enumerate(cuadros):
                fila = info[k] if info is not None and k < len(info) else None
                w.append_data(np.asarray(rotular(Image.open(c), titulo, fila, col)))
        creados.append(mp4)
    except Exception as e:
        print(f"(sin MP4: {e}. Para MP4: pip install imageio imageio-ffmpeg)")
    # GIF reducido
    sel = list(range(0, len(cuadros), P.GIF_CADA))
    imgs = []
    for k in sel:
        fila = info[k] if info is not None and k < len(info) else None
        im = rotular(Image.open(cuadros[k]), titulo + (f"  (×{P.GIF_CADA})" if P.GIF_CADA > 1 else ""),
                     fila, col, ancho=P.GIF_ANCHO)
        imgs.append(im.convert("P", palette=Image.ADAPTIVE, colors=128))
    gif = destino_dir / f"{nombre}.gif"
    imgs[0].save(gif, save_all=True, append_images=imgs[1:], duration=int(1000 / P.FPS_FOTOGRAMAS), loop=0,
                 optimize=True)
    creados.append(gif)
    return creados


def comparativo(nombre, susps=None):
    """GIF (y MP4) lado a lado de la misma corrida con las dos suspensiones."""
    susps = susps or P.LISTA_SUSPENSIONES
    carpetas = [P.resultados_susp(s) / "fotogramas" / nombre for s in susps]
    if not all((c / "cuadros.csv").exists() for c in carpetas):
        return []
    infos = [_info(c) for c in carpetas]
    listas = [sorted(c.glob("cuadro_*.jpg")) for c in carpetas]
    # eje común: los fotogramas del primero (cada GIF_CADA); el otro se empareja por estación s
    base = infos[0]
    destino = P.RESULTADOS / "comparacion" / "videos"; destino.mkdir(parents=True, exist_ok=True)
    s_max_otros = [inf[:, 2].max() for inf in infos[1:]]
    pares = []
    for k in range(0, len(base), P.GIF_CADA):
        s = base[k, 2]
        idx = [k]
        for inf in infos[1:]:
            idx.append(int(np.argmin(np.abs(inf[:, 2] - s))))
        pares.append(idx)
        if any(s > sm + 5 for sm in s_max_otros):      # el otro ya terminó (p. ej. volcó)
            break
    marcos = []
    w1 = P.GIF_ANCHO
    for idx in pares:
        partes = []
        for j, (s_, k) in enumerate(zip(susps, idx)):
            k = min(k, len(listas[j]) - 1)
            im = rotular(Image.open(listas[j][k]), P.SUSPENSIONES[s_]["descripcion"], infos[j][k],
                         _color(P.SUSPENSIONES[s_]["color"]), ancho=w1)
            partes.append(im)
        h = max(p.height for p in partes)
        lienzo = Image.new("RGB", (w1 * len(partes), h), (0, 0, 0))
        for j, p in enumerate(partes):
            lienzo.paste(p, (j * w1, 0))
        marcos.append(lienzo)
    if not marcos:
        return []
    creados = []
    try:
        import imageio.v2 as imageio
        mp4 = destino / f"comparativo_{nombre}.mp4"
        with imageio.get_writer(mp4, fps=P.FPS_FOTOGRAMAS, codec="libx264", quality=8, macro_block_size=8) as w:
            for m in marcos:
                w.append_data(np.asarray(m))
        creados.append(mp4)
    except Exception:
        pass
    gif = destino / f"comparativo_{nombre}.gif"
    pal = [m.convert("P", palette=Image.ADAPTIVE, colors=128) for m in marcos]
    pal[0].save(gif, save_all=True, append_images=pal[1:], duration=int(1000 / P.FPS_FOTOGRAMAS), loop=0, optimize=True)
    creados.append(gif)
    return creados


if __name__ == "__main__":
    if len(sys.argv) > 1:
        nombres = sys.argv[1:]
    else:
        d0 = P.resultados_susp(P.LISTA_SUSPENSIONES[0]) / "fotogramas"
        nombres = sorted(p.name for p in d0.iterdir() if p.is_dir()) if d0.exists() else []
    if not nombres:
        print("No hay fotogramas en resultados/<susp>/fotogramas/. Corre paso3_rampa.py con video = 1.")
    for n in nombres:
        for a in comparativo(n):
            print("Comparativo ->", a)
