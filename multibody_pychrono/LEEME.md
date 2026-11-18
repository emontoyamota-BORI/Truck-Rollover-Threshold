# p22 en PyChrono — suspensión multihoja vs neumática

Simulación dinámica (multicuerpo) del tractocamión 4x2 + semirremolque del p22
(Yoo, Kim & Lee 2025) en la rampa del accidente de la ref. [11] (Zhang et al. 2024),
con dos suspensiones: **muelles de acero multihoja** y **neumática de brazo arrastrado**.
Es el modelo de referencia ("el más exacto") contra el que se compara el modelo analítico.

## Cómo correrlo

1. Anaconda Prompt → `conda activate chrono` → `spyder` (Spyder debe abrirse desde ese entorno).
2. Abre `paso0_verificar.py` y córrelo (F5): debe decir TODO LISTO.
   Opcional para MP4 y Excel: `pip install imageio imageio-ffmpeg openpyxl`.
3. Abre `correr_todo.py` y córrelo (F5). Hace todo: datos → mesa basculante → corridas en la rampa
   → GIF comparativos → gráficas y tablas. Con video se abre la ventana 3D del camión en la rampa;
   **no la cierres**: se cierra sola al terminar cada corrida.

Cada paso también se puede correr suelto:

| Archivo | Necesita PyChrono | Qué hace |
|---|---|---|
| `paso1_preparar.py` | No | Malla de la rampa y JSON del vehículo para las dos suspensiones; tabla de parámetros |
| `paso2_mesa_basculante.py [susp] [--ver]` | Sí | Cargas estáticas y mesa basculante: SRT, despegue por eje, gradiente de balanceo |
| `paso3_rampa.py susp modo v_kmh mu video` | Sí | Una corrida en la rampa; ej. `python paso3_rampa.py aire constante 40 0.8 1` |
| `video.py [corrida]` | No | GIF/MP4 lado a lado multihoja / aire |
| `paso4_resultados.py` | No | Todas las gráficas y tablas comparativas |

Los casos de la rampa (velocidades, fricción, cuáles llevan video) están en `config.py`, lista `CASOS_RAMPA`.
Las suspensiones (Kφ por eje, frecuencia, separación de resortes) están en `config.py`, `SUSPENSIONES`.

## Qué sale y dónde

```
resultados/
  comparacion/
    figuras/   mesa_LTR_por_eje, mesa_balanceo, mesa_SRT_por_eje,
               resumen_LTR_vs_velocidad_R60, resumen_balanceo_vs_ay,
               <corrida>_velocidad / _LTR / _balanceo / _cargas /
               _deslizamiento_fuerza_lateral / _direccion_articulacion_ay
    tablas/    tabla_parametros_suspension, tabla_masas_y_cg, tabla_mesa_basculante,
               tabla_mesa_curvas, tabla_rampa_estacionario, tabla_eventos,
               tabla_validacion_ref11  (+ resultados_p22_suspensiones.xlsx)
    videos/    comparativo_<corrida>.gif / .mp4  (lado a lado, emparejado por estación)
  multihoja/ y aire/
    mesa/        mesa_basculante.csv, resumen_mesa.json
    telemetria/  <corrida>.csv (cada 0.05 s, 44 columnas), <corrida>_resumen.json
    tablas/      <corrida>_cada5m.csv (telemetría cada 5 m de estación)
    videos/      <corrida>.gif y .mp4
    fotogramas/  <corrida>/cuadro_*.jpg + cuadros.csv
```

### Tablas para comparar con el modelo analítico

- **tabla_mesa_curvas.csv**: LTR de cada eje, LTR total, balanceo de tractor y remolque, y Fz por lado
  cada 0.01 g de a_y. Es la referencia directa de la Parte B (SRT físico eje por eje).
- **tabla_mesa_basculante.csv**: SRT del primer despegue, a_y de despegue de cada eje, orden de
  despegue y gradiente de balanceo [°/g].
- **tabla_rampa_estacionario.csv**: promedios en las curvas de radio constante (R = 125 m y 60 m) de
  cada corrida: v, a_y, LTR total y por eje, balanceo, ángulo de dirección, articulación, ángulo de
  deslizamiento por eje, Fz por lado y Fy por eje. Es la referencia del modelo estacionario con MF (Parte C).
- **tabla_validacion_ref11.csv**: LTR en la curva de 60 m a 30/40/50 km/h contra TruckSim (0.225/0.467/0.796).

## Convención de signos (la rampa gira a la derecha)

Lado L = exterior de la curva, lado R = interior. a_y, F_y y dirección positivas hacia el centro de la
curva; balanceo positivo hacia el exterior; LTR = (ext − int)/(ext + int); desv_lat positiva hacia el exterior.
En la rampa, `balanceo_tractor_deg` y `balanceo_remolque_deg` son RELATIVOS A LA VÍA (ya sin el peralte) y
`balanceo_susp_eje<i>_deg` es el giro del chasis relativo a cada eje (lo que deforma la suspensión: el φ del
modelo con Kφ). Fz es la carga normal a la vía y Fy la fuerza lateral en el plano de la vía.

## Notas

- Chrono no abre rutas con acentos: copia `datos/` a `C:\Users\<usuario>\p22_susp_datos` y escribe los
  fotogramas en `C:\Users\<usuario>\p22_susp_videos`; luego se copian a `resultados/`.
- Kφ por eje viene de los rangos medidos por UMTRI (p35) y Fu & Cebon (p36). La altura del centro de
  balanceo no cambia entre suspensiones (la fija la plantilla de eje de Chrono): limitación a declarar.
- Las columnas `direccion_deg`, `alfa*_deg`, `ay_vr` (= v·r) y `r_guinada` usan funciones de PyChrono 10; si tu
  versión no las tiene salen vacías (NaN) y todo lo demás funciona igual.

## Historial

- v2 (7 nov 2026): control de velocidad con rampa de aceleración (la v1 llegaba 10–25 km/h pasada en las
  corridas de 55 y 60 km/h); despegue de eje sólo si LTR ≥ 0.98 se sostiene 0.3 s; balanceo relativo a la
  vía y por eje; Fz y Fy medidas en el plano de la vía; a_y medida como v·r.
- v3 (7 nov 2026): la normal de la vía se calcula con la geometría de la rampa; ya no se usa
  terreno.GetNormal(), que alteraba los contactos (en v2 las ruedas delanteras de la multihoja quedaron
  giradas 35° y el remolque mostraba picos falsos de balanceo de −39°). Si al arrancar las ruedas
  delanteras se giran > 20° sin orden de dirección, la corrida se detiene y lo reporta.
- v4 (7 nov 2026): la velocidad objetivo arranca en 5 m/s (no en 0). Con objetivo 0 el control frenaba
  mientras el camión se asentaba y, en la multihoja, el mecanismo de dirección quedaba girado ~35°.
- v5 (7 nov 2026): el camión se asienta 2 s frenado y sin dirección antes de que arranque el conductor
  (igual que la mesa basculante, que nunca falló). Si aun así se reporta FALLA AL ARRANCAR, correr_todo
  reintenta con otra altura inicial (3, 1 y 6 cm). Casos nuevos a 44 y 47 km/h para ubicar el vuelco.
