# =============================================================================
#  Generates Figures 3, 4, 7, 8 and 9 of the article from the results of
#  tractor_semitrailer_srt.py (Models A, B and C).
#
#  Code author : Luis Montoya (UANL-FIME, Mexico)
#  License     : MIT
#
#  Usage (from the folder that contains both files):
#      python make_article_figures.py
#  The PNG files are written to the folder figures/ (200 dpi).
#  Figures 1, 2, 5 and 6 are schematics or PyChrono outputs and are not produced here.
# =============================================================================
import os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Run the main script in this namespace without showing its own diagnostic plots
_show = plt.show
plt.show = lambda *a, **k: plt.close('all')
HERE = os.path.dirname(os.path.abspath(__file__))
exec(compile(open(os.path.join(HERE, 'tractor_semitrailer_srt.py'), encoding='utf-8').read(),
             'tractor_semitrailer_srt.py', 'exec'))
plt.show = _show

OUT = os.path.join(HERE, 'figures'); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 10})
YOO = 'Yoo et al. (2025)'
NAME = {'multileaf': 'Multi-leaf steel', 'air': 'Air (trailing arm)'}
TITLE = {'multileaf': 'Multi-leaf steel suspension', 'air': 'Air (trailing arm) suspension'}

# ---- Figure 3: tilt table, axle LTR versus lateral acceleration ---------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for j, n in enumerate(SUSP):
    c, ch = CURVE[n], CHRONO_TILT[n]; i = c['i_peak'] + 1
    for k in range(4):
        ax[j].plot(c['ayg'][:i], c['LTR'][:i, k], '-', color=AXLE_COLORS[k], lw=2, label=f'Axle {k+1}, Model B')
        ax[j].plot(ch['ay'], ch['LTR'][k], 'o', color=AXLE_COLORS[k], ms=4, mfc='w', label=f'Axle {k+1}, Model C (PyChrono)')
        ax[j].plot(ay_p, L_A[:, k], ':', color=AXLE_COLORS[k], lw=1.5, label=f'Axle {k+1}, Model A ({YOO})')
    ax[j].axvline(c['SRT'], color='k', lw=1, ls='--')
    ax[j].text(c['SRT'] - 0.003, 0.03, f"SRT {c['SRT']:.3f} g", fontsize=8, ha='right')
    ax[j].set_title(TITLE[n]); ax[j].set_xlabel('$a_y/g$ [-]'); ax[j].grid(alpha=.3)
    ax[j].set_xlim(0, 0.3); ax[j].set_ylim(0, 1.05)
ax[0].set_ylabel('axle load transfer ratio LTR$_i$ [-]')
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=4, fontsize=7.5)
plt.tight_layout(rect=[0, 0.2, 1, 1]); plt.savefig(os.path.join(OUT, 'fig3_tilt_table.png'), dpi=200); plt.close()

# ---- Figure 7: LTR versus speed in the R = 60 m curve -------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
ax[0].plot(v_kmh, R60['A'][0], 'k:', lw=2, label=f'Model A ({YOO}), SRT 0.5 g')
ax[1].plot(v_kmh, R60['A'][1].max(1), 'k:', lw=2, label=f'Model A ({YOO})')
for n, s in SUSP.items():
    ch = CHRONO_R60[n]
    ax[0].plot(v_kmh, R60[n][0], '-', color=s['color'], lw=2, label=f'Model B, {NAME[n]}')
    ax[0].plot(ch['v'], ch['LTR'], 'o', color=s['color'], mfc='w', ms=6, label=f'Model C (PyChrono), {NAME[n]}')
    ax[1].plot(v_kmh, R60[n][1].max(1), '-', color=s['color'], lw=2, label=f'Model B, {NAME[n]}')
ax[0].plot(TRUCKSIM_R60['v'], TRUCKSIM_R60['LTR'], 'ks', mfc='w', ms=6, label='TruckSim (Zhang et al., 2024)')
for a_ in ax:
    a_.axhline(LTR_MAX, color='gray', ls='--', lw=1); a_.grid(alpha=.3); a_.set_ylim(0, 1.05)
    a_.set_xlabel('speed in the R = 60 m curve [km/h]')
ax[0].legend(fontsize=7.5, loc='lower right'); ax[1].legend(fontsize=7.5, loc='lower right')
ax[0].set_ylabel('total LTR [-]'); ax[0].set_title('(a) Total LTR (criterion of Yoo et al., 2025)')
ax[1].set_ylabel('max$_i$ LTR$_i$ [-]'); ax[1].set_title('(b) Most loaded axle (proposed criterion)')
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'fig7_R60_LTR.png'), dpi=200); plt.close()

# ---- Figure 8: safe-speed profile along the ramp -------------------------------
FIG8_CASES = [(0.3, 1.0, '(a) $\\mu_0$ = 0.3, constant friction (case of Fig. 8a in Yoo et al., 2025)'),
              (0.8, 1.0, '(b) $\\mu_0$ = 0.8, dry pavement (as in the PyChrono runs)')]
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for j, (mu0, c, t_) in enumerate(FIG8_CASES):
    vf, _ = run_profile(None, mu0, c)
    ax[j].plot(S, vf * KMH, 'k:', lw=2, label=f'Model A ({YOO}), SRT 0.5 g')
    for n, m in MODELS.items():
        vfn, _ = run_profile(m, mu0, c)
        ax[j].plot(S, vfn * KMH, '-', color=SUSP[n]['color'], lw=2, label=f'Model B, {NAME[n]}')
    for n in SUSP:
        lo, hi = CHRONO_V_ROLLOVER[n]
        ax[j].fill_between([363.5, 437.3], lo, hi, color=SUSP[n]['color'], alpha=0.15)
    ax[j].text(366, CHRONO_V_ROLLOVER['air'][0] - 5, 'PyChrono rollover band', fontsize=7)
    ax[j].set_title(t_, fontsize=10); ax[j].set_ylabel('$v_{final}$ [km/h]'); ax[j].grid(alpha=.3); ax[j].set_ylim(30, 85)
    ax[j].legend(fontsize=7, loc='lower left')
ax[-1].set_xlabel('station s [m]')
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'fig8_safe_speed.png'), dpi=200); plt.close()

# ---- Figure 9: peak friction versus normalized load ----------------------------
x = np.linspace(0.05, 2.6, 200)
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.axvspan(0.0, 2.2, color='gray', alpha=0.08)
ax.plot(x, power_law(x, mu0_pot, c_pot), 'k-', lw=2, label=f'Power law fitted to the data: $\\mu_0$ = {mu0_pot:.3f}, c = {c_pot:.3f}')
ax.plot(x, np.clip(p_SE + q_SE * x, 0, None), 'g-', lw=2, label='1987 Magic Formula peak, $D/F_z = a_1F_z + a_2$')
ax.plot(x, mu0_pot * x ** (0.4 - 1), ':', color='gray', lw=1.5, label='Power law, c = 0.4 (Yoo et al., 2025)')
for label, col, lab in zip(['295/80R22.5 (axle 2)', '275/70R22.5 (axles 3-4)'], ['purple', 'mediumorchid'], ['295/80R22.5', '275/70R22.5']):
    Fz = x * FZ0_KN
    mu_tm = mf_D(Fz, MF_TM[label]) / (1000 * Fz); mu_tm1 = mf_D(FZ0_KN, MF_TM[label]) / (1000 * FZ0_KN)
    ax.plot(x, mu_tm / mu_tm1 * mu0_pot, '-', color=col, lw=1.2, label=f'TMeasy {lab}, multibody model (reference only)')
ax.plot(x_SE, mu_SE, 'ko', mfc='w', ms=7, label='Truck-tire data (Fig. 5 of Yoo et al., 2025)')
ax.set_xlabel('$F_z/F_{z0}$ [-]  ($F_{z0}$ = 18 kN)'); ax.set_ylabel('peak friction $\\mu$ [-]'); ax.set_ylim(0, 1.6); ax.grid(alpha=.3)
ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'fig9_friction.png'), dpi=200); plt.close()

# ---- Figure 4: roll angle, CG displacement and stability curve ----------------
fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
for n, sp in SUSP.items():
    c, ch, col = CURVE[n], CHRONO_TILT[n], sp['color']; i = c['i_peak'] + 1
    ax[0].plot(c['ayg'][:i], np.degrees(c['phi_s'][:i]), '-', color=col, lw=2, label=f'{NAME[n]}, Model B')
    ax[0].plot(ch['ay'], ch['phi'], 'o', color=col, ms=4, mfc='w', label=f'{NAME[n]}, Model C (PyChrono)')
    ax[1].plot(c['ayg'][:i], 100 * c['dy'][:i], '-', color=col, lw=2, label=f'{NAME[n]}, Model B')
    ax[2].plot(np.degrees(c['phi_s']), c['ayg'], '-', color=col, lw=2, label=f'{NAME[n]}, Model B')
    for k in range(4):
        if np.isfinite(c['ay_liftoff'][k]):
            j = int(np.argmax(c['LTR'][:, k] >= 1.0))
            ax[2].plot(np.degrees(c['phi_s'][j]), c['ayg'][j], 'v', color=AXLE_COLORS[k], ms=8)
ax[2].plot([], [], 'kv', ms=7, label='axle lift-off (color = axle)')
ax[0].set_xlabel('$a_y/g$ [-]'); ax[0].set_ylabel('sprung-mass roll $\\phi_s$ [deg]'); ax[0].set_xlim(0, 0.3); ax[0].set_ylim(0, 12)
ax[0].set_title('(a) Roll angle')
ax[1].set_xlabel('$a_y/g$ [-]'); ax[1].set_ylabel('lateral CG shift $\\Delta y$ [cm]'); ax[1].set_title('(b) Lateral CG displacement')
ax[2].set_xlabel('$\\phi_s$ [deg]'); ax[2].set_ylabel('sustainable $a_y/g$ [-]'); ax[2].set_title('(c) Stability curve (SRT = peak)')
for a_ in ax: a_.grid(alpha=.3); a_.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'fig4_roll_cg.png'), dpi=200); plt.close()
print('Figures written to', OUT)
