# =============================================================================
#  PHYSICS-BASED STATIC ROLLOVER THRESHOLD (SRT) AND SAFE CURVE SPEED OF A
#  TRACTOR-SEMITRAILER: MULTI-LEAF VERSUS AIR SUSPENSION
#
#  Companion code of the article
#    L. Montoya, J. Ibarra, R. A. Ramirez-Mendoza, "A Physics-Based Static
#    Rollover Threshold with Measured Suspension Roll Stiffness for Safe
#    Curve-Speed Estimation of a Tractor-Semitrailer: Multi-Leaf versus Air
#    Suspension" (manuscript submitted for publication).
#
#  Code author : Luis Montoya
#                Facultad de Ingenieria Mecanica y Electrica (FIME),
#                Universidad Autonoma de Nuevo Leon (UANL), Mexico
#  Version     : 1.0 (November 2026)
#  License     : MIT (see LICENSE in the repository root)
#
#  CONTENTS (section numbers refer to the article)
#    Part A  Model A: reproduction of the quasi-static safe-speed model of
#            Yoo, Kim & Lee (2025), Appl. Sci. 15, 12269. The SRT is imposed
#            and the lateral load transfer is distributed among the axles with
#            a fixed roll-stiffness split rho_i (Sec. 2.5).
#    Part B  Model B: proposed axle-by-axle static roll-plane model. The SRT
#            is computed from the measured roll stiffness of each suspension,
#            the tire vertical stiffness and the lateral shift of the centre of
#            gravity (Sec. 2.6).
#    Part C  Results: Model A vs Model B vs Model C (PyChrono multibody
#            simulation) on the tilt table, in the R = 60 m curve of the ramp
#            and along the whole ramp; sensitivity of the SRT; comparison with
#            the ORNL track tests (Secs. 3.1-3.6 and Discussion).
#    Part D  Sensitivity to the tire-road friction law: load-sensitive power
#            law of Yoo et al. (2025) versus the peak of the 1987 Magic Formula,
#            both fitted to the same truck-tire data (Secs. 2.7 and 3.7).
#            Sections D.3 (second half) and D.4 contain a supplementary
#            slip-angle analysis that is NOT reported in the article; it is kept
#            for completeness and does not affect any reported result.
#
#  HOW TO RUN
#    Single self-contained file. Runs in Google Colab (one cell) or in any
#    Python >= 3.9 environment with numpy, scipy, matplotlib and pandas:
#        python tractor_semitrailer_srt.py
#    The results of Model C (PyChrono) are not recomputed here: they are
#    hard-coded in Part B.3 (and D.4) from the output tables of the multibody
#    model (folder multibody_pychrono/resultados/comparacion/tablas/ of this repository).
#
#  TAGS USED IN THE COMMENTS
#    [Yoo Eq. N]           equation N of Yoo et al. (2025)
#    [Author Year, Eq. N]  equation of another reference (list in Part B header)
#    [Implicit]            relation described only in the text of the source
#    [Assumption]          value not given by the source; the choice is explained
#    [Model C]             same input as the PyChrono model, so that the
#                          comparison isolates the method and not the data
#    SI units throughout: kg, N, m, s, rad. Speeds are in m/s internally and
#    are converted to km/h only for printing and plotting.
# =============================================================================

SHOW_MODEL_A_FIGURES = False   # True = also plot the reproduction of Figs. 4-8 of Yoo et al. (about 1 min more)

# #############################################################################
#  PART A - MODEL A: REPRODUCTION OF YOO, KIM & LEE (2025)
# #############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.optimize import curve_fit

g = 9.81            # gravitational acceleration [m/s^2]  [Assumption: not given by the source]
KMH = 3.6           # conversion factor m/s -> km/h

# =============================================================================
# A.1 PARAMETERS
# =============================================================================

# ---- A.1.1 Vehicle (Yoo et al., Table 1) -----------------------------------
# 4x2 tractor + two-axle semitrailer (Yoo et al., Fig. 2). The text calls it a
# "four-axle rigid truck", but the figure and the original accident study
# (Zhang et al., 2024) show a tractor-semitrailer.
M_TRACTOR = 6.8e3        # empty tractor mass [kg]                 (Table 1)
M_TRAILER = 10.7e3       # empty semitrailer mass [kg]             (Table 1)
                         # [Assumption] Zhang et al. (2024) give 7.05 t, which
                         # would give a total of 31.1 t (close to the "32 t" of
                         # the text). 10.7 t is kept because Fig. 7a of Yoo et
                         # al. implies a total of about 35 t.
M_PAYLOAD = 17.2e3       # payload [kg]                            (Table 1)
M_TOTAL   = M_TRACTOR + M_TRAILER + M_PAYLOAD   # total mass m [kg] = 34 700 kg
H_CG_LADEN = 2.725       # CG height, laden [m]                    (Table 1)
H_CG_EMPTY = 1.60        # CG height, empty [m]                    (Table 1)

F_AXLE   = np.array([0.21, 0.30, 0.25, 0.24])  # static weight fraction per axle f_i [-]   (Table 1)
T_AXLE   = np.array([2.04, 1.82, 1.82, 1.82])  # track width t_i [m]                       (Table 1)
RHO_AXLE = np.array([0.15, 0.20, 0.30, 0.35])  # roll-stiffness split rho_i [-]            (Table 1)
N_TIRES_SIDE = np.array([1, 2, 2, 2])          # tires per side of each axle [-]
                                               # [Assumption] Fig. 2: single tires on
                                               # axle 1, dual tires on axles 2-4

# ---- A.1.2 Tire-road friction (Table 1 and Fig. 5) -------------------------
MU0_TRUCK, C_TRUCK = 0.749, 0.821   # power-law fit, truck tire (Fig. 5)
MU0_CAR,   C_CAR   = 0.880, 0.893   # power-law fit, passenger-car tire (Fig. 5)
DELTA   = 0.2      # safety margin delta [-]: only (1 - delta) = 80 % of mu is used (Table 1)
MU_MIN  = 0.05     # lower bound of mu [-]                                         (Table 1)
MU_MAX  = 0.95     # upper bound of mu [-]                                         (Table 1)

N0_PAPER = 4000.0  # reference load N0 [N] as printed in Table 1
N0_CALIB = 18000.0 # [Assumption] reference load N0 [N] calibrated to reproduce
                   # Fig. 8a (minimum speed of about 40 km/h); used in all results.
                   # Alternative 'static' (each tire uses its own static load as
                   # reference): reproduces mu = mu0 on the tangents of Figs. 6b
                   # and 7b better, but gives a minimum of about 44 km/h in Fig. 8a.
                   # With N0 = 4000 N the minimum is about 31 km/h; 4000 N does not
                   # correspond to the load of any tire of this truck (15-45 kN),
                   # and Fig. 5 normalizes with Fz/Fz0 between 0.2 and 1.15, which
                   # indicates that N0 is of the order of the rated tire load.
                   # Both values are plotted when SHOW_MODEL_A_FIGURES = True.
SKID_MODE = 'min'  # [Implicit] skid criterion:
                   #  'min'   -> the wheel with the lowest mu governs (Yoo et al.,
                   #             Sec. 3.3: "the minimum mu among wheels governs")
                   #  'total' -> sum of the capacities of all wheels (Sec. 2.3)

# ---- A.1.3 Stability (Table 1) ----------------------------------------------
LTR_MAX = 0.80     # load-transfer-ratio threshold [-]                      (Table 1)
SRT_TOL = 0.005    # tolerance of the SRT calibration [g]                   (Table 1)

# ---- A.1.4 Longitudinal dynamics (Table 1) ---------------------------------
RHO_AIR  = 1.225   # air density [kg/m^3]
CDA      = 8.5     # drag coefficient x frontal area [m^2]
F_RR     = 0.008   # rolling-resistance coefficient [-]
P_MAX    = 300e3   # maximum engine power [W]
A_CAP    = 0.5     # comfort limit when accelerating a_acc,cap [m/s^2]   (Table 1)
A_COMF   = 0.5     # comfort limit when braking a_comfort [m/s^2]        [Assumption: same value]
EPS_V    = 1.0     # epsilon of Eq. 15 to avoid division by zero [m/s]  [Assumption]
V_CAP    = 80 / KMH  # entry speed and cap on tangents [m/s]             [Assumption: Figs. 6a and 8a start at 80 km/h]
DS       = 1.0     # station step Delta s [m]                           [Assumption]


# =============================================================================
# A.2 RAMP GEOMETRY  (Zhang et al., 2024, Table 3; article Sec. 2.2)
# =============================================================================
# Each segment: (s_start, s_end, kappa_start, kappa_end, grade %, superelevation_start %, superelevation_end %)
#  - kappa = 1/R curvature [1/m]; it varies linearly along the spirals (clothoids).  [Assumption]
#  - Superelevation interpolated linearly along the spirals.                         [Assumption]
#  - The -2 % cross-slope of the tangents is taken as 0 (it is not curve banking).   [Assumption]
#  - The ramp turns to the RIGHT: the left side (L) is the outer side, as in
#    Fig. 7a of Yoo et al., where the L wheels are the loaded ones.
R1, R2 = 125.0, 60.0
SEGMENTS = [
    (0.0,     50.0,    0.0,  0.0,  -2.45, 0, 0),   # entry tangent
    (50.0,    146.8,   0.0,  1/R1, -2.45, 0, 6),   # entry spiral
    (146.8,   308.032, 1/R1, 1/R1, -2.45, 6, 6),   # circular curve R = 125 m
    (308.032, 363.469, 1/R1, 1/R2, -2.45, 6, 7),   # transition spiral
    (363.469, 437.268, 1/R2, 1/R2, -2.45, 7, 7),   # circular curve R = 60 m (accident site)
    (437.268, 543.935, 1/R2, 0.0,   0.893, 7, 0),  # exit spiral
    (543.935, 643.935, 0.0,  0.0,   0.0,   0, 0),  # exit tangent
]

def build_ramp(ds=DS):
    """Return arrays per station s: curvature kappa, grade alpha, superelevation beta and X-Y layout."""
    s = np.arange(0.0, SEGMENTS[-1][1] + 1e-9, ds)      # station [m]
    kap = np.zeros_like(s); grade = np.zeros_like(s); supel = np.zeros_like(s)
    for (a, b, k0, k1, G, e0, e1) in SEGMENTS:
        m = (s >= a) & (s <= b)
        x = (s[m] - a) / (b - a)                         # relative position in the segment [0..1]
        kap[m]   = k0 + (k1 - k0) * x                    # curvature [1/m]
        grade[m] = G                                     # grade [%]
        supel[m] = e0 + (e1 - e0) * x                    # superelevation [%]
    alpha = np.arctan(grade / 100)                       # grade angle alpha [rad] (negative = downhill)
    beta  = np.arctan(supel / 100)                       # superelevation angle beta [rad]
    # Plan-view layout for X-Y plots (right-hand turn => heading decreases)
    heading = -np.concatenate([[0], np.cumsum(kap[:-1] * ds)])
    X = np.concatenate([[0], np.cumsum(np.cos(heading[:-1]) * ds)])
    Y = np.concatenate([[0], np.cumsum(np.sin(heading[:-1]) * ds)])
    return dict(s=s, kap=kap, alpha=alpha, beta=beta, X=X, Y=Y)

ROAD = build_ramp()
S = ROAD['s']


# =============================================================================
# A.3 MODEL (Yoo Eqs. 1-10)
# =============================================================================

def mu_N(N, mu0, c, N0):
    """[Yoo Eq. 6] Load-sensitive friction: mu(N) = mu0*(N/N0)^(c-1), bounded to [mu_min, mu_max].
       N [N] tire load; mu0 [-] friction at the reference load N0; c [-] exponent
       (c = 1 -> constant mu; c < 1 -> mu decreases as the load increases)."""
    if isinstance(c, str) and c == 'MF':
        # Part D: load sensitivity of the 1987 Magic Formula peak, mu = D/F_z = a1*F_z + a2,
        # normalized so that mu(N0) = mu0 (friction scaling factor lambda_mu of Pacejka, 2006).
        # MF_SE = (p, q): mu_rel(x) = (p + q*x)/(p + q), with x = N/N0; defined in Part D.
        x = np.maximum(N, 1.0) / N0
        return np.clip(mu0 * (MF_SE[0] + MF_SE[1] * x) / (MF_SE[0] + MF_SE[1]), MU_MIN, MU_MAX)
    return np.clip(mu0 * (np.maximum(N, 1.0) / N0) ** (c - 1.0), MU_MIN, MU_MAX)


LOAD_MODEL = None   # None = Model A (Yoo Eq. 4); a PhysicalSRTModel instance (Part B) = Model B

def wheel_loads(v, rd, lam, m=M_TOTAL, h=H_CG_LADEN):
    """Vertical loads at speed v [m/s] at every station of the road rd.
       lam = SRT calibration factor (Yoo Eq. 12). Returns a dictionary."""
    kap, al, be = rd['kap'], rd['alpha'], rd['beta']
    # [Yoo Eq. 1] geometric and effective lateral acceleration [m/s^2]
    a_y_geom = v**2 * kap * np.cos(be)                   # (v^2/R)*cos(beta), with kappa = 1/R
    a_y_eff  = a_y_geom - g * np.sin(be)                 # superelevation reduces it by g*sin(beta)
    # [Yoo Eq. 2] longitudinal acceleration imposed by the grade [m/s^2]
    a_x_req = -g * np.sin(al)
    # [Yoo Eq. 3] total forces that the tires must provide [N]
    Fy_need = m * np.abs(a_y_eff)
    Fx_need = m * a_x_req
    # [Implicit] static load per axle [N]: N_i,stat = f_i*m*g
    N_stat = F_AXLE * m * g
    # [Yoo Eq. 5] base load of each side of the axle
    base = N_stat * (np.cos(al) * np.cos(be))[:, None] / 2
    if LOAD_MODEL is None:
        # [Yoo Eq. 4] lateral load transfer per axle [N]: Delta_i = rho_i*m*h/t_i*|a_y,eff|
        #             (lam scales all Delta_i so that the model rolls over at the target SRT)
        Delta = lam * RHO_AXLE * m * h / T_AXLE * np.abs(a_y_eff)[:, None]
    else:
        # Model B: Delta_i comes from the axle-by-axle roll-plane model with the roll stiffness
        # of each suspension; LTR_i = Delta_i/(2*base)  ->  Delta_i = 2*base*LTR_i(|a_y,eff|/g)
        Delta = 2 * base * LOAD_MODEL.ltr(np.abs(a_y_eff) / g)
    # [Yoo Eq. 5] load on each side of the axle [N]: N = N_stat*cos(alpha)*cos(beta)/2 +/- Delta/2
    N_ext = base + Delta / 2      # outer side (left, L)
    N_int = base - Delta / 2      # inner side (right, R); may be < 0 if the axle lifts off
    # [Implicit] load per tire [N]: shared equally among the tires of each side
    Nw_ext = np.clip(N_ext, 0, None) / N_TIRES_SIDE
    Nw_int = np.clip(N_int, 0, None) / N_TIRES_SIDE
    return dict(a_y_eff=a_y_eff, Fy_need=Fy_need, Fx_need=Fx_need,
                N_ext=N_ext, N_int=N_int, Nw_ext=Nw_ext, Nw_int=Nw_int)


def LTR(v, rd, lam, m=M_TOTAL, h=H_CG_LADEN):
    """[Yoo Eq. 10] Whole-vehicle load transfer ratio LTR = sum(N_out - N_in) / sum(N_out + N_in) [-]."""
    C = wheel_loads(v, rd, lam, m, h)
    return (C['N_ext'] - C['N_int']).sum(1) / (C['N_ext'] + C['N_int']).sum(1)


def LTR_axle(v, rd, lam, m=M_TOTAL, h=H_CG_LADEN):
    """LTR of each axle: (N_out,i - N_in,i)/(N_out,i + N_in,i) [-] (axle criterion of Model B)."""
    C = wheel_loads(v, rd, lam, m, h)
    return np.minimum((C['N_ext'] - C['N_int']) / (C['N_ext'] + C['N_int']), 1.0)


def wheel_mu(v, rd, lam, mu0, c, N0, m=M_TOTAL):
    """Friction coefficient of each tire (outer and inner) and the load per tire.
       N0 can be a number [N] or 'static': each tire then uses its own static load
       as reference (N0_w = f_i*m*g/2/n_i), so that mu = mu0 on the tangents."""
    C = wheel_loads(v, rd, lam, m)
    if isinstance(N0, str):
        N0 = F_AXLE * m * g / 2 / N_TIRES_SIDE       # [N] static load of each tire, per axle
    return mu_N(C['Nw_ext'], mu0, c, N0), mu_N(C['Nw_int'], mu0, c, N0), C


def no_skid(v, rd, lam, mu0, c, N0, delta=DELTA, mode=SKID_MODE, m=M_TOTAL):
    """True where the tires can provide the required force.
       [Yoo Eq. 7] friction circle per tire: (Fx/((1-delta)*mu*N))^2 + (Fy/((1-delta)*mu*N))^2 <= 1
       [Yoo Eq. 8] the lateral force is shared in proportion to the load: F_y,w = F_y,need*N_w/sum(N)
       [Implicit]  the longitudinal force due to the grade is shared in the same way."""
    mu_e, mu_i, C = wheel_mu(v, rd, lam, mu0, c, N0, m)
    F_dem = np.hypot(C['Fy_need'], C['Fx_need'])                 # total required force [N]
    if mode == 'total':
        # total capacity = sum of (1-delta)*mu(N_w)*N_w over all tires [N]
        cap = ((1 - delta) * (mu_e * C['Nw_ext'] + mu_i * C['Nw_int']) * N_TIRES_SIDE).sum(1)
        return F_dem <= cap
    # mode 'min': with load-proportional sharing every tire uses the same
    # fraction F_dem/sum(N) of its load; the tire with the lowest mu saturates first.
    N_tot = (np.clip(C['N_ext'], 0, None) + np.clip(C['N_int'], 0, None)).sum(1)
    inner_loaded = C['Nw_int'] > 1.0                              # inner tires still in contact
    mu_min = np.minimum(mu_e.min(1), np.where(inner_loaded, mu_i, np.inf).min(1))
    return F_dem / N_tot <= (1 - delta) * mu_min


# =============================================================================
# A.4 SRT CALIBRATION  [Yoo Eq. 12]
# =============================================================================

def SRT_model(m=M_TOTAL, h=H_CG_LADEN):
    """SRT_model = a_y,crit/g at which LTR(a_y,crit) = 1, on level ground and with lam = 1 [g]."""
    level = dict(kap=np.array([1.0]), alpha=np.zeros(1), beta=np.zeros(1))  # R = 1 m: v^2 = a_y
    lo, hi = 0.0, 3 * g
    while (hi - lo) / g > SRT_TOL / 10:            # bisection on a_y [m/s^2]
        a = (lo + hi) / 2
        if LTR(np.sqrt([a]), level, 1.0, m, h)[0] < 1.0: lo = a
        else: hi = a
    return lo / g

def lambda_factor(SRT_target, m=M_TOTAL, h=H_CG_LADEN):
    """[Implicit] Global factor that scales all Delta_i so that the model rolls over at SRT_target.
       Since the LTR of the model is linear in a_y, lam = SRT_model / SRT_target.
       (The source scales only when SRT_model > SRT_target; here it is always applied.)"""
    return SRT_model(m, h) / SRT_target


# =============================================================================
# A.5 SAFE SPEED  [Yoo Eqs. 11 and 13 + skid criterion]
# =============================================================================

def v_max(condition, n, v_cap=V_CAP, it=40):
    """Largest v in [0, v_cap] that satisfies a condition, by bisection at all stations simultaneously."""
    lo, hi = np.zeros(n), np.full(n, v_cap)
    ok_cap = condition(hi)
    for _ in range(it):
        mid = (lo + hi) / 2
        ok = condition(mid)
        lo, hi = np.where(ok, mid, lo), np.where(ok, hi, mid)
    return np.where(ok_cap, v_cap, lo)

def safe_speed(rd, SRT_target, mu0, c, N0, delta=DELTA, mode=SKID_MODE, m=M_TOTAL, v_cap=V_CAP):
    """Return v_safe, v_roll, v_skid [m/s] per station and the factor lam.
       v_cap [m/s] = upper bound of the search (speed cap on tangents)."""
    n = len(rd['kap'])
    if LOAD_MODEL is None:
        lam = lambda_factor(SRT_target, m)
        # [Yoo Eqs. 11 and 13] v_roll = max{ v | LTR(v) <= LTR_max }   (whole-vehicle LTR)
        v_roll = v_max(lambda v: LTR(v, rd, lam, m) <= LTR_MAX, n, v_cap)
    else:
        lam = 1.0      # Model B: no calibration; the SRT results from the model (SRT_target is not used)
        # axle criterion: v_roll = max{ v | max_i LTR_i(v) <= LTR_max }
        v_roll = v_max(lambda v: LTR_axle(v, rd, lam, m).max(1) <= LTR_MAX, n, v_cap)
    # [Implicit] v_skid = max{ v | the tires can provide the required force }  (Yoo et al., Sec. 2.3)
    v_skid = v_max(lambda v: no_skid(v, rd, lam, mu0, c, N0, delta, mode, m), n, v_cap)
    # [Implicit] the safe speed is the lower of the two (and of the cap on tangents)
    return np.minimum(v_roll, v_skid), v_roll, v_skid, lam


# =============================================================================
# A.6 LONGITUDINAL DYNAMICS AND FINAL PROFILE  [Yoo Eqs. 14-19]
# =============================================================================

def accelerations(i, v, rd, lam, mu0, c, N0, delta=DELTA, m=M_TOTAL):
    """Available acceleration a_acc and deceleration a_brk [m/s^2] at station i and speed v [m/s]."""
    sub = dict(kap=rd['kap'][i:i+1], alpha=rd['alpha'][i:i+1], beta=rd['beta'][i:i+1])
    mu_e, mu_i, C = wheel_mu(np.array([v]), sub, lam, mu0, c, N0, m)
    N_tot = (np.clip(C['N_ext'], 0, None) + np.clip(C['N_int'], 0, None)).sum()
    F_res = 0.0
    for mu, Nw in ((mu_e, C['Nw_ext']), (mu_i, C['Nw_int'])):
        cap = (1 - delta) * mu * Nw                                   # capacity per tire [N]
        Fy_w = C['Fy_need'][0] * Nw / max(N_tot, 1.0)                 # [Yoo Eq. 8] lateral force per tire [N]
        # [Yoo Eq. 9] remaining longitudinal capacity per tire [N]
        Fx_res = cap * np.sqrt(np.clip(1 - (Fy_w / np.maximum(cap, 1e-9))**2, 0, None))
        F_res += (Fx_res * N_TIRES_SIDE).sum()
    a_long = F_res / m                                                # [Yoo Eq. 14] [m/s^2]
    al = rd['alpha'][i]
    F_aero  = 0.5 * RHO_AIR * CDA * v**2                              # [Implicit] aerodynamic drag [N]
    F_roll  = F_RR * m * g * np.cos(al)                               # [Implicit] rolling resistance [N]
    F_grade = m * g * np.sin(al)                                      # [Implicit] grade resistance [N]
    resist = (F_aero + F_roll + F_grade) / m                          # [m/s^2]
    a_acc_tire = a_long - resist                                      # [Implicit] (not defined by Yoo Eq. 16)
    a_pow = P_MAX / (m * max(v, EPS_V)) - resist                      # [Yoo Eq. 15] engine power limit
    a_acc = min(a_acc_tire, a_pow, A_CAP)                             # [Yoo Eq. 16]
    a_brk = min(a_long + g * np.sin(al), A_COMF)                      # [Yoo Eq. 17]
    # [Assumption] if an acceleration is negative (tires saturated in a curve) it is
    # set to 0: v_safe already guarantees equilibrium at that station.
    return max(a_acc, 0.0), max(a_brk, 0.0)

def speed_profile(rd, v_safe, lam, mu0, c, N0, delta=DELTA, m=M_TOTAL,
                  v_entry=V_CAP, tol=1e-4, it_max=50):
    """[Yoo Eqs. 18 and 19] Forward-backward sweep until convergence. Returns v_final [m/s]."""
    v = v_safe.copy()
    v[0] = min(v[0], v_entry)
    for _ in range(it_max):
        v_prev = v.copy()
        # [Yoo Eq. 18] forward pass: do not exceed the speed reachable by accelerating from the previous station
        for i in range(1, len(v)):
            a_acc, _ = accelerations(i - 1, v[i - 1], rd, lam, mu0, c, N0, delta, m)
            v[i] = min(v[i], np.sqrt(v[i - 1]**2 + 2 * a_acc * DS))
        # [Yoo Eq. 19] backward pass: the vehicle must be able to brake before the next station
        for i in range(len(v) - 2, -1, -1):
            _, a_brk = accelerations(i + 1, v[i + 1], rd, lam, mu0, c, N0, delta, m)
            v[i] = min(v[i], np.sqrt(v[i + 1]**2 + 2 * a_brk * DS))
        if np.max(np.abs(v - v_prev)) < tol:
            break
    return v


# =============================================================================
# A.7 REFERENCE VALUES DIGITIZED FROM THE FIGURES OF YOO ET AL. (approx. +/- 1-2 km/h)
# =============================================================================
REF_FIG8_CONST = np.array([[0, 80], [30, 80], [100, 72], [200, 62], [330, 47], [400, 47.5], [480, 47.5], [600, 59]])
REF_FIG8_MUN   = np.array([[0, 73.5], [50, 69], [100, 61], [150, 56], [200, 56.5], [300, 46], [350, 40], [480, 40], [600, 53.5]])
REF_FIG6B_C04  = np.array([[0, 0.30], [100, 0.24], [170, 0.215], [300, 0.29], [400, 0.204], [500, 0.29], [600, 0.27]])
REF_FIG5_TRUCK = np.array([[0.24, 0.965], [0.66, 0.812], [0.875, 0.765], [0.92, 0.758]])
REF_FIG5_CAR   = np.array([[0.58, 0.933], [0.855, 0.895], [0.925, 0.888], [0.95, 0.884]])


if SHOW_MODEL_A_FIGURES:
    # =============================================================================
    # A.8 REPRODUCTION OF THE FIGURES OF YOO ET AL. (2025)
    # =============================================================================
    print(f"Total mass m = {M_TOTAL/1e3:.1f} t ; h_cg = {H_CG_LADEN} m")
    print(f"Uncalibrated SRT_model (Yoo Eq. 12, lam = 1) = {SRT_model():.3f} g "
          f"(SSF = t/2h = {T_AXLE[1]/(2*H_CG_LADEN):.3f}; the factor of about 2 comes from Yoo Eqs. 4-5)")

    # ---------------------------------------------------------------- Figure 5 ---
    # [Yoo Fig. 5] mu versus Fz/Fz0 for a truck tire and a passenger-car tire.
    def power_law_x(x, mu0, c):           # mu = mu0*x^(c-1), with x = Fz/Fz0
        return mu0 * x ** (c - 1)
    p_trk, _ = curve_fit(power_law_x, REF_FIG5_TRUCK[:, 0], REF_FIG5_TRUCK[:, 1], p0=[0.75, 0.8])
    p_car, _ = curve_fit(power_law_x, REF_FIG5_CAR[:, 0], REF_FIG5_CAR[:, 1], p0=[0.88, 0.9])
    print(f"\n[Fig. 5] Fit to digitized points: truck mu0={p_trk[0]:.3f}, c={p_trk[1]:.3f} "
          f"(source 0.749 / 0.821); car mu0={p_car[0]:.3f}, c={p_car[1]:.3f} (source 0.880 / 0.893)")
    x = np.linspace(0.2, 1.15, 100)
    plt.figure(figsize=(6, 4.5))
    plt.plot(x, power_law_x(x, MU0_TRUCK, C_TRUCK), 'k-', label='Truck-tire model (mu0=0.749, c=0.821)')
    plt.plot(x, power_law_x(x, MU0_CAR, C_CAR), 'k--', label='Car-tire model (mu0=0.880, c=0.893)')
    plt.plot(*REF_FIG5_TRUCK.T, 'o', mfc='w', mec='k', label='Truck tire, Segel & Ervin (digitized from Fig. 5)')
    plt.plot(*REF_FIG5_CAR.T, 's', mfc='w', mec='k', label='Car tire, Segel & Ervin (digitized from Fig. 5)')
    plt.xlabel('$F_z/F_{z0}$ [-]'); plt.ylabel('μ [-]')
    plt.title('Model A, reproduction of Fig. 5 of Yoo et al. (2025): friction vs normalized load')
    plt.grid(alpha=.3); plt.legend(fontsize=8); plt.tight_layout(); plt.show()

    # ---------------------------------------------------------------- Figure 4 ---
    # [Yoo Fig. 4] Map of Delta = v_skid - v_roll [km/h] versus station and mu, for SRT = 0.3/0.5/0.7.
    #   Red (Delta > 0): ROLLOVER governs. Blue (Delta < 0): SKID governs. Green: both equal.
    #   [Assumption] constant mu (c = 1) and NO margin delta: with these assumptions the
    #   boundary mu*(s) is about 0.8*SRT, as in the source (about 0.25 / 0.42 / 0.54-0.57).
    #   [Assumption] a speed cap of 150 km/h is used for the map (not 80 km/h), so that both
    #   limits are not truncated in the curves; on the tangents both reach the cap and
    #   Delta = 0 (green), as in the source.
    V_CAP_FIG4 = 150 / KMH
    MU_GRID = np.round(np.arange(0.20, 1.001, 0.025), 3)
    for SRT_target in (0.3, 0.5, 0.7):
        D = np.zeros((len(MU_GRID), len(S)))
        for k, mu in enumerate(MU_GRID):
            _, vr, vs, _ = safe_speed(ROAD, SRT_target, mu, 1.0, N0_CALIB, delta=0.0, mode='total', v_cap=V_CAP_FIG4)
            D[k] = (vs - vr) * KMH
        # boundary mu*(s): first mu of the grid where Delta >= 0 (rollover starts to govern)
        mu_star = np.array([MU_GRID[np.argmax(D[:, j] >= -1e-6)] if (D[:, j] >= -1e-6).any() else np.nan
                            for j in range(len(S))])
        fig = plt.figure(figsize=(13, 4))
        ax = fig.add_axes([0.05, 0.13, 0.5, 0.75])
        im = ax.pcolormesh(S, MU_GRID, D, cmap='jet', norm=TwoSlopeNorm(0, -20, 20), shading='auto')
        ax.plot(S, mu_star, 'k--', label='μ*(s) (Δ = 0)')
        ax.set_xlabel('Station [m]'); ax.set_ylabel('Friction coefficient μ')
        ax.set_title(f'Model A, reproduction of Fig. 4, SRT = {SRT_target}: dominance map'); ax.legend(loc='upper right')
        fig.colorbar(im, ax=ax, label='Δ = v_skid − v_roll [km/h]')
        for p, mu in enumerate((0.3, 0.5, 0.7, 0.9)):
            k = int(np.argmin(np.abs(MU_GRID - mu)))
            axp = fig.add_axes([0.64 + 0.17 * (p % 2), 0.55 - 0.42 * (p // 2), 0.14, 0.33])
            axp.scatter(ROAD['X'], ROAD['Y'], c=D[k], cmap='jet', norm=TwoSlopeNorm(0, -20, 20), s=2)
            axp.set_title(f'μ = {mu}', fontsize=9); axp.set_aspect('equal'); axp.tick_params(labelsize=7)
        plt.show()
        print(f"[Fig. 4] SRT={SRT_target}: mu* in the 60 m curve = {np.nanmedian(mu_star[(S>364)&(S<437)]):.3f} ; "
              f"in the 125 m curve = {np.nanmedian(mu_star[(S>147)&(S<308)]):.3f}  (about 0.8*SRT = {0.8*SRT_target:.2f})")

    # ------------------------------------------------------ Figures 8a and 6a ---
    # [Yoo Fig. 8a] final safe speed with constant mu and with mu(N), c = 0.4.
    # Case of Yoo Figs. 6-8: mu0 = 0.3, SRT = 0.5.
    MU0_CASE, SRT_CASE, C_CASE = 0.3, 0.5, 0.4

    vs_c, vr_c, vk_c, lam_case = safe_speed(ROAD, SRT_CASE, MU0_CASE, 1.0, N0_CALIB)
    v_const = speed_profile(ROAD, vs_c, lam_case, MU0_CASE, 1.0, N0_CALIB)          # constant mu (c = 1)

    vs_n, _, _, _ = safe_speed(ROAD, SRT_CASE, MU0_CASE, C_CASE, N0_CALIB)
    v_muN = speed_profile(ROAD, vs_n, lam_case, MU0_CASE, C_CASE, N0_CALIB)        # mu(N), calibrated N0

    vs_p, _, _, _ = safe_speed(ROAD, SRT_CASE, MU0_CASE, C_CASE, N0_PAPER)
    v_muN_paper = speed_profile(ROAD, vs_p, lam_case, MU0_CASE, C_CASE, N0_PAPER)  # mu(N), N0 = 4000 N

    vs_e, _, _, _ = safe_speed(ROAD, SRT_CASE, MU0_CASE, C_CASE, 'static')
    v_muN_static = speed_profile(ROAD, vs_e, lam_case, MU0_CASE, C_CASE, 'static') # mu(N), N0 = static load of each tire

    plt.figure(figsize=(9, 4.5))
    plt.plot(S, v_const * KMH, 'k-', lw=2, label='constant μ (Model A)')
    plt.plot(S, v_muN * KMH, 'b--', lw=2, label=f'μ(N), c = 0.4, N0 = {N0_CALIB/1e3:.0f} kN (Model A)')
    plt.plot(S, v_muN_static * KMH, 'c-.', lw=1.5, label='μ(N), c = 0.4, N0 = static load of each tire')
    plt.plot(S, v_muN_paper * KMH, 'm:', lw=1.5, label='μ(N), c = 0.4, N0 = 4000 N (Table 1 as printed)')
    plt.plot(*REF_FIG8_CONST.T, 'ko', mfc='w', label='Yoo et al. Fig. 8a, constant μ (digitized)')
    plt.plot(*REF_FIG8_MUN.T, 'bs', mfc='w', label='Yoo et al. Fig. 8a, μ(N) (digitized)')
    plt.xlabel('Station (distance from the start) [m]'); plt.ylabel('Speed [km/h]')
    plt.title('Model A, reproduction of Fig. 8a: safe speed, constant μ vs μ(N)  [Yoo Eqs. 1–19]')
    plt.grid(alpha=.3); plt.legend(fontsize=8); plt.ylim(25, 90); plt.tight_layout(); plt.show()

    dif = (v_const - v_muN)[S <= 600].mean() * KMH
    print(f"\n[Fig. 8a] minimum, constant mu = {v_const.min()*KMH:.1f} km/h (source about 47)")
    print(f"[Fig. 8a] minimum, mu(N)      = {v_muN.min()*KMH:.1f} km/h (source about 40); with N0 = 4000 N: {v_muN_paper.min()*KMH:.1f} km/h; "
          f"with N0 = static load of each tire: {v_muN_static.min()*KMH:.1f} km/h")
    print(f"[Sec. 3.3] mean difference constant mu - mu(N) = {dif:.1f} km/h (source about 7.5)")

    # ---------------------------------------------------------- Figures 6a and 6b --
    # [Yoo Fig. 6] same reference speed (constant-mu profile) and different c:
    #   6a reference speed; 6b minimum mu among all tires.
    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    styles = {1.0: 'k-', 0.8: 'b--', 0.6: 'r:', 0.4: 'g-.'}
    for c in (1.0, 0.8, 0.6, 0.4):
        ax[0].plot(S, v_const * KMH, styles[c], label=f'c = {c}')
        mu_e, mu_i, C = wheel_mu(v_const, ROAD, lam_case, MU0_CASE, c, N0_CALIB)
        mu_i = np.where(C['Nw_int'] > 1.0, mu_i, np.inf)
        ax[1].plot(S, np.minimum(mu_e.min(1), mu_i.min(1)), styles[c], label=f'c = {c}')
    ax[1].plot(*REF_FIG6B_C04.T, 'gs', mfc='w', label='Yoo et al. Fig. 6b, c = 0.4 (digitized)')
    ax[0].set_ylabel('Speed [km/h]'); ax[0].set_title('Model A, reproduction of Fig. 6a: reference safe speed')
    ax[1].set_ylabel('minimum tire μ [-]'); ax[1].set_xlabel('Station [m]')
    ax[1].set_title('Model A, reproduction of Fig. 6b: minimum μ among tires  [Yoo Eq. 6]')
    for a in ax: a.grid(alpha=.3); a.legend(fontsize=8)
    plt.tight_layout(); plt.show()

    # ---------------------------------------------------------- Figures 7a and 7b --
    # [Yoo Fig. 7] load and mu per axle side (L = outer, R = inner), c = 0.4, mu0 = 0.3, SRT = 0.5.
    #   With N0 = N0_CALIB, mu < 0.3 on the tangents; in the source it starts at 0.3 (use N0 = 'static' to see this).
    #   Note: the source plots DYNAMIC loads from TruckSim; here they are the quasi-static loads of Yoo Eq. 5.
    mu_e, mu_i, C = wheel_mu(v_const, ROAD, lam_case, MU0_CASE, C_CASE, N0_CALIB)
    colors = ['tab:blue', 'tab:orange', 'gold', 'tab:purple']
    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for i in range(4):
        ax[0].plot(S, C['N_ext'][:, i] / 1e3, '-', color=colors[i], label=f'L{i+1} (outer)')
        ax[0].plot(S, np.clip(C['N_int'][:, i], 0, None) / 1e3, '--', color=colors[i], label=f'R{i+1} (inner)')
        ax[1].plot(S, mu_e[:, i], '-', color=colors[i], label=f'L{i+1}')
        ax[1].plot(S, mu_i[:, i], '--', color=colors[i], label=f'R{i+1}')
    ax[0].set_ylabel('Load per axle side [kN]')
    ax[0].set_title('Model A, reproduction of Fig. 7a: vertical wheel loads  [Yoo Eqs. 4–5]')
    ax[1].set_ylabel('μ [-]'); ax[1].set_xlabel('Station [m]')
    ax[1].set_title('Model A, reproduction of Fig. 7b: friction coefficient per wheel  [Yoo Eq. 6]')
    for a in ax: a.grid(alpha=.3); a.legend(fontsize=7, ncol=4)
    plt.tight_layout(); plt.show()
    print(f"[Fig. 7] maximum outer load = {C['N_ext'].max()/1e3:.1f} kN (source about 74-78 kN, TruckSim); "
          f"minimum outer mu = {mu_e.min():.3f} (source: text 0.229, figure about 0.27)")

    # =============================================================================
    # A.9 LIMITATIONS OF THE REPRODUCTION
    # =============================================================================
    print("""
    NOT REPRODUCIBLE, OR REPRODUCIBLE ONLY AS A TREND:
     - Figs. 6c and 8b (lateral offset): obtained with TruckSim (multibody dynamic model).
       The model of Eqs. 1-19 is quasi-static and does not compute a trajectory.
     - Fig. 7 (dynamic TruckSim loads): here the quasi-static loads of Eq. 5 are plotted;
       the peaks and oscillations of the figure cannot be obtained with this model.
     - Figs. 1-3: schematics and TruckSim images.
     - N0 and the skid criterion are not fully defined in the source: N0 was calibrated
       against Fig. 8a (see N0_CALIB); N0 = 4000 N (as printed) and N0 = static load of
       each tire are also shown. No option reproduces Fig. 8a and mu = 0.3 on the
       tangents of Figs. 6b/7b at the same time: the figures of the source are not
       mutually consistent.
    """)


# #############################################################################
#
#  PART B - MODEL B: PHYSICAL SRT WITH THE ROLL STIFFNESS OF EACH SUSPENSION
#
#  Replaces Yoo Eq. 4 (Delta_i = lam*rho_i*m*h/t_i*|a_y|, with lam calibrated
#  to impose an SRT) by an axle-by-axle static roll-plane model:
#    - one rigid sprung mass (tractor + semitrailer, because the fifth wheel of
#      the multibody model transmits roll: both units roll together on the
#      tilt table) that rolls by phi_s;
#    - four axles (unsprung masses) that roll by phi_u,k on their tires
#      (tire vertical stiffness);
#    - each suspension with its own roll stiffness K_phi,k and roll-centre
#      height h_rc,k;
#    - the lateral shift of the weight Delta_y (CG displacement) acting as a
#      destabilizing moment;
#    - axle-by-axle lift-off: an axle whose inner tires have lifted no longer
#      contributes roll stiffness; the SRT is the PEAK of the lateral
#      acceleration that the vehicle can sustain (not the first lift-off).
#
#  REFERENCES CITED IN THE COMMENTS
#    Fu & Cebon (2002), Int. J. Heavy Vehicle Systems 9(4): Eq. 1
#        (K_phi = K_aux + 2*k_s*p^2), Eqs. 2-3 (equilibrium of the sprung and
#        unsprung masses), Eq. 6 (lift-off), Eq. 8 (upper bound with
#        K_phi -> infinity), Table 2 (measured suspensions)
#    Winkler & Ervin (1999), UMTRI-99-19: Eq. 1 (W*h*a_y = sum(dF*T/2) - W*dy),
#        Figs. 15-16 (sequential lift-off; the SRT is the peak of the total moment)
#    Kamnik, Boettiger & Hunt (2003): Eqs. 15-17 (roll-plane load transfer with
#        sprung and unsprung masses)
#    Afshari & Garcia-Pozuelo (2025): Eq. 7 (roll angle with rigid tires, a
#        particular case used as a check)
#    Lundahl et al. (2015): Eq. 27 (lift-off acceleration per axle)
#    Winkler, Karamihas & Bogard (1992), SAE 922426: measured K_phi ranges
#    Oak Ridge National Laboratory (2009): measured lift-off accelerations of a
#        tractor-semitrailer (Part C.8)
#    Zhang et al. (2024): ramp geometry and TruckSim results
# #############################################################################

# ---- B.1 Vehicle parameters of the roll-plane model --------------------------
# EXACTLY the data of the PyChrono model are used, so that the comparison
# measures the difference of METHOD and not of input data.
M_AXLE   = F_AXLE * M_TOTAL                          # static mass per axle [kg]           [Yoo Table 1]
M_U      = np.array([700.0, 1200.0, 900.0, 900.0])  # unsprung mass per axle [kg]         [Model C; Fu & Cebon Table 1: 706/1107/792]
M_S_AXLE = M_AXLE - M_U                              # sprung mass carried by each axle [kg]
# sprung mass and CG height of each unit (tractor and semitrailer roll together)     [Model C]
MS_TRACTOR, H_S_TRACTOR = M_TRACTOR - M_U[0] - M_U[1], 1.10                 # [Assumption, same as Model C]
MS_TRAILER, H_S_TRAILER = M_TRAILER + M_PAYLOAD - M_U[2] - M_U[3], 2.725    # [Zhang et al. 2024: CG of the laden trailer]
MS_TOTAL = MS_TRACTOR + MS_TRAILER
H_S      = (MS_TRACTOR * H_S_TRACTOR + MS_TRAILER * H_S_TRAILER) / MS_TOTAL  # CG height of the sprung mass [m]

T_HALF   = T_AXLE / 2                                # half track (to the centre of the dual pair) [m]  [Yoo; Winkler & Ervin, note 9]
E_DUAL   = 0.33 / 2                                  # half spacing of the dual tires [m]              [Model C, Assumption]

# Tires: vertical stiffness and loaded radius with the SAME formula as the multibody model (TMeasy):
#   k_t = F_LI / (0.16 * section height), with F_LI the load of the tire load index.
TIRES = [dict(width=0.295, aspect=0.80, rim=22.5, LI=152),   # axle 1: 295/80R22.5, single
         dict(width=0.295, aspect=0.80, rim=22.5, LI=148),   # axle 2: 295/80R22.5, dual
         dict(width=0.275, aspect=0.70, rim=22.5, LI=145),   # axles 3-4: 275/70R22.5, dual
         dict(width=0.275, aspect=0.70, rim=22.5, LI=145)]
LI_LOAD = {145: 2900, 148: 3150, 152: 3550}           # [kg] ETRTO load-index table
K_T, H_U = np.zeros(4), np.zeros(4)
for k, L in enumerate(TIRES):
    sec = L['width'] * L['aspect']                        # section height [m]
    K_T[k] = LI_LOAD[L['LI']] * g / (0.16 * sec)          # vertical stiffness of ONE tire [N/m]  [Model C]
    free_radius = L['rim'] * 0.0254 / 2 + sec
    H_U[k] = free_radius - (M_AXLE[k] * g / 2 / N_TIRES_SIDE[k]) / K_T[k]   # loaded radius = height of the unsprung mass [m]

# Roll stiffness of the tires of each axle: K_rt = sum(k_t*y^2) over all its tires
#   single tires: 2*k_t*t^2  (Fu & Cebon, note to Eq. 6)      dual tires: 4*k_t*(t^2 + e^2)
K_RT = np.where(N_TIRES_SIDE == 1, 2 * K_T * T_HALF**2, 4 * K_T * (T_HALF**2 + E_DUAL**2))

# Roll centre: in the multibody model it is set by the axle template (SolidThreeLinkAxle),
# whose triangular link is attached 0.26 m above the axle centre. It is the SAME for both
# suspensions (a limitation stated in the documentation of the multibody model).
H_RC_CHRONO = H_U + 0.26                                                   # [Model C]

# ---- B.2 Suspensions (the same two as in the multibody model) ----------------
#   K_PHI: total roll stiffness per axle [N*m/rad] = K_aux + 2*k_s*p^2     [Fu & Cebon, Eq. 1]
#   values at the centre of the ranges measured by UMTRI (1 kN*m/deg = 57.3 kN*m/rad)
#   [Winkler, Karamihas & Bogard 1992; Fu & Cebon, Table 2]
#   K_RES: part provided by the springs (2*k_s*p^2), computed by the multibody model from the bounce frequency
#   H_RC_LIT: measured roll-centre heights, used in the sensitivity study  [Fu & Cebon, Table 2; Winkler et al. 1992]
SUSP = {
    'multileaf': dict(name='Multi-leaf steel springs', color='tab:red',
                      K_PHI=np.array([150, 470, 750, 750]) * 1e3,
                      K_RES=np.array([134.8, 349.2, 294.8, 281.6]) * 1e3,
                      H_RC_LIT=np.array([0.49, 0.75, 0.71, 0.71])),
    'air':       dict(name='Trailing-arm air suspension', color='tab:blue',
                      K_PHI=np.array([150, 240, 890, 890]) * 1e3,
                      K_RES=np.array([134.8, 185.9, 156.9, 149.9]) * 1e3,
                      H_RC_LIT=np.array([0.49, 0.58, 0.58, 0.58])),
}
for s in SUSP.values():
    s['K_AUX'] = s['K_PHI'] - s['K_RES']          # [Fu & Cebon, Eq. 1] part not provided by the springs
    s['H_RC'] = H_RC_CHRONO.copy()                # baseline: same roll centre as Model C

# ---- B.3 Results of Model C (PyChrono output tables) --------------------------
# Tilt table (multibody_pychrono/paso2_mesa_basculante.py; output table
# resultados/comparacion/tablas/tabla_mesa_curvas.csv): a_y/g = tan(theta), LTR per axle,
# whole-vehicle LTR and tractor roll angle [deg]
CHRONO_TILT = {
 'multileaf': dict(
  ay=np.round(np.arange(0, 0.295, 0.01), 2),
  LTR=np.array([
   [-0.00024,0.008005,0.03153,0.05277,0.06409,0.07338,0.09356,0.1145,0.1269,0.138,0.1564,0.1756,0.19,0.2024,0.2187,0.2374,0.2529,0.2666,0.2798,0.2972,0.3133,0.329,0.3474,0.359,0.3764,0.4007,0.4501,0.5263,0.6558,0.8697],
   [3e-05,0.01046,0.05057,0.09161,0.106,0.1191,0.154,0.1921,0.2112,0.2263,0.2583,0.2924,0.3154,0.3332,0.3611,0.3936,0.4194,0.4399,0.4629,0.4918,0.52,0.5454,0.5767,0.595,0.6256,0.6674,0.76,0.9133,1,1],
   [4e-05,0.01504,0.07875,0.1431,0.1661,0.1856,0.2403,0.3003,0.3303,0.3537,0.4026,0.4575,0.4935,0.5213,0.5633,0.6157,0.6555,0.6855,0.723,0.7678,0.8126,0.855,0.8962,0.9285,0.9716,1,1,1,1,1],
   [4e-05,0.01521,0.08022,0.1459,0.1693,0.1891,0.2447,0.3059,0.3365,0.3602,0.41,0.4661,0.5027,0.531,0.5737,0.6271,0.6676,0.698,0.7363,0.7805,0.8271,0.8714,0.9117,0.9437,0.9892,1,1,1,1,1]]),
  LTR_tot=np.array([-2e-05,0.01222,0.06067,0.1092,0.1273,0.1428,0.1844,0.2299,0.2531,0.2714,0.309,0.3504,0.3781,0.3997,0.4323,0.4718,0.5025,0.5263,0.5544,0.5885,0.6227,0.6548,0.688,0.7116,0.7461,0.7738,0.8117,0.8736,0.9293,0.9708]),
  phi=np.array([0.00018,0.07906,0.4781,0.902,1.052,1.164,1.504,1.894,2.09,2.229,2.525,2.882,3.11,3.27,3.53,3.857,4.105,4.283,4.513,4.783,5.048,5.304,5.57,5.765,6.039,6.546,7.644,9.513,12.66,19.41]),
  liftoff=[0.29333, 0.27259, 0.2382, 0.23451], SRT_first=0.23451, SRT_last=0.29333),
 'air': dict(
  ay=np.round(np.arange(0, 0.265, 0.01), 2),
  LTR=np.array([
   [-0.00024,0.008,0.03157,0.05275,0.06417,0.07346,0.09356,0.1145,0.1274,0.1377,0.1557,0.1758,0.1905,0.2037,0.2185,0.2377,0.2528,0.2644,0.2795,0.299,0.3231,0.3314,0.3514,0.3913,0.4739,0.6361,0.899],
   [3e-05,0.00802,0.03562,0.06338,0.07348,0.08317,0.1071,0.1326,0.1467,0.1568,0.1784,0.2029,0.2185,0.2329,0.2503,0.2737,0.29,0.3028,0.3199,0.342,0.3739,0.378,0.4038,0.4466,0.5436,0.7378,1],
   [5e-05,0.01646,0.0874,0.1591,0.1848,0.2064,0.2669,0.3335,0.368,0.3931,0.4467,0.5091,0.5497,0.5814,0.6271,0.6854,0.7292,0.7609,0.8028,0.8635,0.9221,0.9577,1,1,1,1,1],
   [5e-05,0.01681,0.08994,0.1638,0.1902,0.2124,0.2746,0.3432,0.3787,0.4044,0.4596,0.5238,0.5656,0.5981,0.6451,0.7051,0.7499,0.7799,0.8238,0.8884,0.943,0.9815,1,1,1,1,1]]),
  LTR_tot=np.array([-2e-05,0.01223,0.0607,0.1091,0.1273,0.1429,0.1843,0.2294,0.2535,0.2711,0.308,0.3505,0.3785,0.4013,0.4323,0.4723,0.502,0.5233,0.5526,0.5941,0.6367,0.6574,0.6852,0.7056,0.7529,0.8447,0.979]),
  phi=np.array([0.00019,0.07891,0.4775,0.9014,1.053,1.165,1.505,1.894,2.092,2.232,2.527,2.883,3.113,3.275,3.536,3.861,4.109,4.278,4.501,4.836,5.183,5.382,5.719,6.666,8.763,13.05,20.23]),
  liftoff=[0.26231, 0.25672, 0.2107, 0.20523], SRT_first=0.20523, SRT_last=0.26231),
}
# Ramp (multibody_pychrono/paso3_rampa.py; output table tabla_rampa_estacionario.csv),
# R = 60 m curve at constant speed, mu = 0.8:
# mean whole-vehicle LTR and rollover flag
CHRONO_R60 = {
 'multileaf': dict(v=[29.99, 40.05, 44.03, 47.06, 49.87], LTR=[0.1379, 0.4172, 0.5617, 0.6744, 0.8926],
                   rollover=[False, False, False, False, True]),
 'air':       dict(v=[29.99, 40.04, 44.03, 46.94, 49.77], LTR=[0.1359, 0.4300, 0.5589, 0.7805, 0.9075],
                   rollover=[False, False, False, True, True]),
}
CHRONO_V_ROLLOVER = {'multileaf': (47.1, 49.9), 'air': (44.0, 46.9)}  # [km/h] last run without / first run with rollover
TRUCKSIM_R60 = dict(v=[30, 40, 50], LTR=[0.225, 0.467, 0.796], v_rollover=54.9)  # TruckSim, Zhang et al. (2024), Table 4


# ---- B.4 Axle-by-axle roll-plane model ----------------------------------------
def equilibrium(phi_s, susp, lifted, gn=g, h_s=H_S, e_load=0.0):
    """Solve the static equilibrium for a given sprung-mass roll angle phi_s [rad].
       Unknowns x = [a_y, phi_u1, phi_u2, phi_u3, phi_u4]  (a_y [m/s^2], phi_u = roll angle of each axle [rad]).
       lifted[k] = True if the inner tires of axle k are no longer in contact with the ground.
       gn = component of gravity normal to the ground [m/s^2] (g on level ground; g*cos(theta) on the tilt table).
       e_load = lateral offset of the payload towards the outside [m] (off-centre load, Winkler & Ervin, Fig. 17)."""
    K, h_rc = susp['K_PHI'], susp['H_RC']
    ms, mu, hu = M_S_AXLE, M_U, H_U
    d = h_s - h_rc                                           # height of the sprung CG above each roll centre [m]
    A = np.zeros((5, 5)); b = np.zeros(5)
    # (A) sprung mass, moments about the roll centres   [Fu & Cebon Eq. 3 summed over axles; Kamnik et al. Eq. 17]
    #     sum(m_s,k*d_k*(a_y + gn*phi_s)) + M_load = sum(K_phi,k*(phi_s - phi_u,k))
    A[0, 0] = np.sum(ms * d)
    A[0, 1:] = K
    b[0] = (K.sum() - gn * np.sum(ms * d)) * phi_s - M_PAYLOAD * gn * e_load
    for k in range(4):
        # (B) axle k, moments about the centre of the tire contact patch   [Fu & Cebon Eq. 2, written per axle]
        #     M_tire,k - (m_s,k*gn*h_rc,k + m_u,k*gn*h_u,k)*phi_u,k - K_phi,k*(phi_s - phi_u,k) = (m_s,k*h_rc,k + m_u,k*h_u,k)*a_y
        #     in contact: M_tire,k = K_rt,k*phi_u,k
        #     lifted:     M_tire,k = W_k*t_k  (the whole load on the outer tires)   [Winkler & Ervin, Fig. 15]
        A[k + 1, 0] = -(ms[k] * h_rc[k] + mu[k] * hu[k])
        A[k + 1, k + 1] = K[k] - gn * (ms[k] * h_rc[k] + mu[k] * hu[k]) + (0.0 if lifted[k] else K_RT[k])
        b[k + 1] = K[k] * phi_s - (M_AXLE[k] * gn * T_HALF[k] if lifted[k] else 0.0)
    x = np.linalg.solve(A, b)
    a_y, phi_u = x[0], x[1:]
    W = M_AXLE * gn                                          # normal load per axle [N]
    M_tire = np.where(lifted, W * T_HALF, K_RT * phi_u)      # moment of the tires [N*m]
    Delta = M_tire / T_HALF                                  # N_ext - N_int [N]
    LTR_k = Delta / W                                        # [Yoo Eq. 10 applied to each axle]
    # lateral displacement of the CG of the whole vehicle  [Winkler & Ervin Eq. 1: term W*dy]
    dy = (np.sum(ms * (d * phi_s + h_rc * phi_u)) + np.sum(mu * hu * phi_u) + M_PAYLOAD * e_load) / M_TOTAL
    return a_y, phi_u, LTR_k, M_tire, dy


def roll_curve(susp, tilt=False, h_s=H_S, e_load=0.0, dphi=np.radians(0.02), phi_max=np.radians(30)):
    """Sweep phi_s from 0 and track the axle-by-axle lift-off.
       tilt=True reproduces the tilt table of the multibody model: a_y = g*sin(theta),
       gn = g*cos(theta), and the reported a_y/g = tan(theta).
       Returns a dictionary with the full curve, the stable branch (up to the peak) and the SRT."""
    lifted = np.zeros(4, bool)
    out = {k: [] for k in ('ayg', 'phi_s', 'phi_u', 'LTR', 'M', 'dy')}
    ay_liftoff = np.full(4, np.nan)                          # a_y/g at which each axle lifts off
    order = []                                               # actual lift-off sequence (as phi_s increases)
    gn = g
    for phi_s in np.arange(0.0, phi_max, dphi):
        for _ in range(4 if tilt else 1):                    # fixed-point iteration on the tilt table: gn depends on theta
            while True:
                a, pu, ltr, M, dy = equilibrium(phi_s, susp, lifted, gn, h_s, e_load)
                new = (~lifted) & (ltr >= 1.0)               # [Fu & Cebon Eq. 6] the axle lifts off when N_int = 0
                if not new.any(): break
                for k in np.where(new)[0]:
                    ay_liftoff[k] = a / gn
                    order.append(k + 1)
                lifted = lifted | new
            if tilt:
                gn = g / np.sqrt(1 + (a / gn) ** 2)          # gn = g*cos(theta) with tan(theta) = a/gn
        for key, val in zip(('ayg', 'phi_s', 'phi_u', 'LTR', 'M', 'dy'), (a / gn, phi_s, pu, np.minimum(ltr, 1), M, dy)):
            out[key].append(val)
        peak = max(out['ayg'])
        if lifted.all() or (peak > 0 and len(out['ayg']) > 50 and a / gn < 0.5 * peak):
            break
    out = {k: np.array(v) for k, v in out.items()}
    i_peak = int(np.argmax(out['ayg']))                      # [Winkler & Ervin, Figs. 15-16] SRT = peak of the curve
    out.update(SRT=out['ayg'][i_peak], i_peak=i_peak, ay_liftoff=ay_liftoff,
               SRT_first=ay_liftoff[order[0] - 1] if order else np.nan,
               order='-'.join(map(str, order)))
    # lift-offs before the peak (stable) and after it (vehicle already overturning)  [Winkler & Ervin, Fig. 16]
    out['before_peak'] = [k for k in order if np.argmax(out['LTR'][:, k - 1] >= 1.0) <= i_peak]
    # roll gradient in the linear range [deg/g]
    j = out['ayg'] < 0.1
    out['roll_gradient'] = np.degrees(np.polyfit(out['ayg'][j], out['phi_s'][j], 1)[0])
    return out


class PhysicalSRTModel:
    """Couples the roll-plane model (Model B) with the safe-speed framework of Part A:
       returns LTR_k(|a_y,eff|/g) for Yoo Eq. 5 (level ground). Beyond the SRT there is
       no equilibrium and LTR_k = 1."""
    def __init__(self, susp, **kw):
        self.susp = susp
        self.c = roll_curve(susp, tilt=False, **kw)
        i = self.c['i_peak'] + 1
        self.ay, self.LTR = self.c['ayg'][:i], self.c['LTR'][:i]
        self.SRT = self.c['SRT']
    def ltr(self, ayg):
        ayg = np.atleast_1d(np.abs(ayg))
        L = np.column_stack([np.interp(ayg, self.ay, self.LTR[:, k]) for k in range(4)])
        L[ayg >= self.SRT] = 1.0
        return L


def ltr_model_a_level(ayg, SRT_target):
    """Axle and whole-vehicle LTR of Model A on level ground (Yoo Eqs. 4, 5, 10, 12)."""
    lam = lambda_factor(SRT_target)
    Delta = lam * RHO_AXLE * M_TOTAL * H_CG_LADEN / T_AXLE * (np.atleast_1d(ayg)[:, None] * g)
    base = F_AXLE * M_TOTAL * g / 2
    L_axle = np.minimum(Delta / (2 * base), 1.0)
    L_tot = np.minimum(Delta.sum(1) / (2 * base.sum()), 1.0)
    return L_axle, L_tot



# #############################################################################
#
#  PART C - RESULTS: MODEL A (Yoo et al.) | MODEL B (physical SRT) | MODEL C (PyChrono)
#
# #############################################################################
import pandas as pd
pd.set_option('display.width', 200); pd.set_option('display.max_columns', 30)
AXLE_COLORS = ['tab:blue', 'tab:orange', 'gold', 'tab:purple']   # one color per axle in all plots
SRT_MODEL_A = 0.5      # SRT imposed in Model A, case of Yoo Figs. 6-8 (Table 1 gives 0.3 / 0.5 / 0.7)
H_CG_VEHICLE = (MS_TOTAL * H_S + np.sum(M_U * H_U)) / M_TOTAL   # CG height of the whole vehicle [m]

def explain(title, text):
    """Print, BEFORE each plot, what it shows and how to read it."""
    print("\n" + "-" * 78 + f"\nPLOT: {title}\n" + "-" * 78)
    for line in text.strip().splitlines():
        print("  " + line.strip())

# ---- C.1 Parameter tables (article Secs. 2.3-2.4) ---------------------------
print("=" * 78); print("C.1  PARAMETERS OF THE ROLL-PLANE MODEL (same as the multibody model)"); print("=" * 78)
tab_veh = pd.DataFrame({
    'axle': [1, 2, 3, 4],
    'static load [kg]': M_AXLE.round(0),
    'unsprung mass m_u [kg]': M_U,
    'sprung mass m_s [kg]': M_S_AXLE.round(0),
    'track T [m]': T_AXLE,
    'tires per side': N_TIRES_SIDE,
    'k_t per tire [kN/m]': (K_T / 1e3).round(0),
    'loaded radius h_u [m]': H_U.round(3),
    'tire roll stiffness K_rt [MN*m/rad]': (K_RT / 1e6).round(2),
    'roll centre h_rc, Model C [m]': H_RC_CHRONO.round(3)})
print(tab_veh.to_string(index=False))
print(f"Total sprung mass = {MS_TOTAL:.0f} kg (tractor {MS_TRACTOR:.0f} kg at {H_S_TRACTOR} m; "
      f"semitrailer {MS_TRAILER:.0f} kg at {H_S_TRAILER} m) -> h_s = {H_S:.3f} m")
print(f"CG height of the whole vehicle h = {H_CG_VEHICLE:.3f} m (Model C: 2.255 m). "
      f"Model A uses h = {H_CG_LADEN} m for the WHOLE mass.")
rows = []
for n, s in SUSP.items():
    for k in range(4):
        rows.append([n, k + 1, s['K_PHI'][k] / 1e3, s['K_RES'][k] / 1e3, s['K_AUX'][k] / 1e3,
                     s['K_AUX'][k] / s['K_PHI'][k], s['H_RC_LIT'][k]])
tab_susp = pd.DataFrame(rows, columns=['suspension', 'axle', 'total K_phi [kN*m/rad]', 'springs 2*k_s*p^2 [kN*m/rad]',
                                       'auxiliary K_aux [kN*m/rad]', 'auxiliary fraction [-]', 'measured h_rc [m]'])
print("\nSuspensions [Winkler et al. 1992; Fu & Cebon 2002, Eq. 1 and Table 2]:")
print(tab_susp.round(3).to_string(index=False))

# ---- C.2 Consistency checks of the implementation ----------------------------
print("\n" + "=" * 78); print("C.2  CONSISTENCY CHECKS"); print("=" * 78)
s0 = SUSP['multileaf']
a, pu, ltr, M, dy = equilibrium(np.radians(3.0), s0, np.zeros(4, bool))
res = M.sum() - (MS_TOTAL * H_S * a + np.sum(M_U * H_U) * a + M_TOTAL * g * dy)
print(f"1) [Winkler & Ervin Eq. 1] sum of tire moments - (W*h*a_y + W*dy) = {res:.2e} N*m  (must be ~0)")
K_RT_backup = K_RT.copy(); K_RT = K_RT * 1e6                 # (temporarily) rigid tires
a, pu, ltr, M, dy = equilibrium(np.radians(3.0), s0, np.zeros(4, bool))
K_RT = K_RT_backup
d_mean = np.sum(M_S_AXLE * (H_S - s0['H_RC'])) / MS_TOTAL
phi_rigid = MS_TOTAL * a * d_mean / (s0['K_PHI'].sum() - MS_TOTAL * g * d_mean)
print(f"2) [Afshari & Garcia-Pozuelo Eq. 7] rigid tires: phi model = 3.000 deg, phi closed form = {np.degrees(phi_rigid):.3f} deg")
print(f"3) Rigid-body bound t/h [Winkler & Ervin]: {T_HALF[1] / H_CG_VEHICLE:.3f} g; "
      f"the physical SRT must be lower.")
print(f"4) Model A without calibration (lam = 1): SRT = {SRT_model():.3f} g  (about 2*t/h because Yoo Eq. 5 uses +/-Delta/2)")

# ---- C.3 Tilt table: axle LTR versus a_y (article Sec. 3.2) -----------------
print("\n" + "=" * 78); print("C.3  TILT TABLE: SRT AND AXLE LIFT-OFF"); print("=" * 78)
CURVE = {n: roll_curve(s, tilt=True) for n, s in SUSP.items()}   # same test as Model C (gn = g*cos(theta))
ay_p = np.linspace(0, 0.3, 61)
L_A, L_A_tot = ltr_model_a_level(ay_p, SRT_MODEL_A)

explain("C.3  LTR of each axle versus lateral acceleration (tilt table)", """
    Virtual tilt-table test: the ground is tilted until the wheels lift off;
    a_y/g = tan(theta) is the equivalent lateral acceleration. One panel per suspension.
    x axis: a_y/g [-].  y axis: LTR_i = (N_ext - N_int)/(N_ext + N_int) of each axle [-]; 1 = axle lifted off.
    Solid lines: Model B.  Circles: Model C (PyChrono).  Dotted: Model A with SRT = 0.5 g.
    Dashed black line: SRT of Model B (peak).  Grey line: first lift-off in Model C.
    If the lines pass through the circles, the static model reproduces the multibody simulation.
    Model A uses a fixed split rho_i: it does not distinguish between suspensions and lies far below
    (it overestimates stability).""")
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for j, n in enumerate(SUSP):
    c, ch = CURVE[n], CHRONO_TILT[n]
    i = c['i_peak'] + 1
    for k in range(4):
        ax[j].plot(c['ayg'][:i], c['LTR'][:i, k], '-', color=AXLE_COLORS[k], lw=2, label=f'Axle {k+1}, Model B')
        ax[j].plot(ch['ay'], ch['LTR'][k], 'o', color=AXLE_COLORS[k], ms=4, mfc='w', label=f'Axle {k+1}, Model C (PyChrono)')
        ax[j].plot(ay_p, L_A[:, k], ':', color=AXLE_COLORS[k], lw=1.5, label=f'Axle {k+1}, Model A (SRT {SRT_MODEL_A} g)')
    ax[j].axvline(c['SRT'], color='k', lw=1, ls='--')
    ax[j].text(c['SRT'] - 0.003, 0.03, f"Model B SRT\n{c['SRT']:.3f} g", fontsize=8, ha='right')
    ax[j].axvline(ch['SRT_first'], color='gray', lw=1, ls='-.')
    ax[j].text(ch['SRT_first'] + 0.004, 0.03, f"Model C first\nlift-off {ch['SRT_first']:.3f} g", fontsize=8, color='dimgray')
    ax[j].set_title(f"{SUSP[n]['name']}"); ax[j].set_xlabel('$a_y/g$ [-]'); ax[j].grid(alpha=.3)
    ax[j].set_xlim(0, 0.3); ax[j].set_ylim(0, 1.05)
ax[0].set_ylabel('axle LTR [-]')
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc='lower center', ncol=4, fontsize=7.5)
fig.suptitle('Tilt table: axle LTR, Model B (line) vs Model C (circles) vs Model A (dotted)')
plt.tight_layout(rect=[0, 0.17, 1, 0.95]); plt.show()

rows = []
for n in SUSP:
    c, ch = CURVE[n], CHRONO_TILT[n]
    rows.append([n, SRT_MODEL_A, '- (no axle lift-off)',
                 c['SRT_first'], c['SRT'], c['order'],
                 ch['SRT_first'], ch['SRT_last'], '4-3-2-1',
                 100 * (c['SRT'] - ch['SRT_first']) / ch['SRT_first'], c['roll_gradient']])
tab_srt = pd.DataFrame(rows, columns=['suspension', 'Model A imposed SRT [g]', 'Model A lift-off',
                                      'Model B first lift-off [g]', 'Model B SRT (peak) [g]', 'Model B order',
                                      'Model C first lift-off [g]', 'Model C last lift-off [g]', 'Model C order',
                                      'Model B SRT vs Model C first lift-off [%]', 'Model B roll gradient [deg/g]'])
print(tab_srt.round(3).to_string(index=False))
print("Model C: tractor roll gradient = 25.85 deg/g (multi-leaf) and 25.89 deg/g (air).")
print("Note: in Model C the lift-offs after the first one occur with the roll angle already diverging (> 7 deg):")
print("      the table tilts at 1 deg/s and the vehicle is already overturning. Model B places the SRT at the")
print("      peak of the stability curve, as in Winkler & Ervin (1999).")

# ---- C.4 Roll angle, CG displacement and stability curve (Sec. 3.2) ---------
explain("C.4  Roll angle, CG displacement and stability curve", """
    Left: sprung-mass roll angle phi_s [deg] versus a_y/g; line = Model B, circles = Model C.
    Centre: lateral CG displacement dy [cm] versus a_y/g (Winkler & Ervin Eq. 1): how far the weight shifts outwards.
    Right: a_y/g sustained by the vehicle versus phi_s [deg]. It rises while the suspensions resist; each triangle
    marks the lift-off of an axle (which loses its roll stiffness). The PEAK is the SRT (Winkler & Ervin, Fig. 15):
    beyond the peak there is no equilibrium.""")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
for n, s in SUSP.items():
    c, ch, col = CURVE[n], CHRONO_TILT[n], s['color']
    i = c['i_peak'] + 1
    ax[0].plot(c['ayg'][:i], np.degrees(c['phi_s'][:i]), '-', color=col, lw=2, label=f'{n}, Model B')
    ax[0].plot(ch['ay'], ch['phi'], 'o', color=col, ms=4, mfc='w', label=f'{n}, Model C (PyChrono)')
    ax[1].plot(c['ayg'][:i], 100 * c['dy'][:i], '-', color=col, lw=2, label=n)
    ax[2].plot(np.degrees(c['phi_s']), c['ayg'], '-', color=col, lw=2, label=n)
    for k in range(4):
        if np.isfinite(c['ay_liftoff'][k]):
            j = int(np.argmax(c['LTR'][:, k] >= 1.0))
            ax[2].plot(np.degrees(c['phi_s'][j]), c['ayg'][j], 'v', color=AXLE_COLORS[k], ms=8)
            ax[2].annotate(f'axle {k+1}', (np.degrees(c['phi_s'][j]), c['ayg'][j]), fontsize=7,
                           xytext=(4, 4 if n == 'multileaf' else -10), textcoords='offset points')
ax[0].set_xlabel('$a_y/g$'); ax[0].set_ylabel('sprung-mass roll angle φ_s [deg]')
ax[0].set_title('Roll angle: Model B vs Model C'); ax[0].set_xlim(0, 0.3); ax[0].set_ylim(0, 12)
ax[1].set_xlabel('$a_y/g$'); ax[1].set_ylabel('lateral CG shift Δy [cm]')
ax[1].set_title('Lateral CG displacement  [Winkler & Ervin, Eq. 1]')
ax[2].set_xlabel('φ_s [deg]'); ax[2].set_ylabel('sustainable $a_y/g$ [-]')
ax[2].set_title('Stability curve: the SRT is the peak  [Winkler & Ervin, Fig. 15]')
for a_ in ax: a_.grid(alpha=.3); a_.legend(fontsize=8)
plt.tight_layout(); plt.show()
for n in SUSP:
    c = CURVE[n]
    print(f"{n:10s}: CG shift at the SRT = {100*c['dy'][c['i_peak']]:.1f} cm -> loss relative to the rigid body "
          f"= {100*c['dy'][c['i_peak']]/H_CG_VEHICLE/(T_HALF[1]/H_CG_VEHICLE):.0f} % of t/h  "
          f"(Winkler & Ervin: 1 cm of dy = about 1 % of the rigid-body stability)")

# ---- C.5 R = 60 m curve of the ramp (accident site; article Sec. 3.5) --------
print("\n" + "=" * 78); print("C.5  RAMP, R = 60 m CURVE (superelevation 7 %, grade -2.45 %)"); print("=" * 78)
i60 = int(np.argmin(np.abs(S - 400.0)))
st60 = dict(kap=ROAD['kap'][i60:i60+1], alpha=ROAD['alpha'][i60:i60+1], beta=ROAD['beta'][i60:i60+1])
v_kmh = np.arange(25, 70.1, 0.5)
MODELS = {n: PhysicalSRTModel(s) for n, s in SUSP.items()}     # Model B on level ground, used on the ramp

def ltr_R60(model, srt=SRT_MODEL_A):
    """Whole-vehicle LTR and axle LTR versus speed in the R = 60 m curve (model=None -> Model A)."""
    global LOAD_MODEL
    LOAD_MODEL = model
    lam = 1.0 if model is not None else lambda_factor(srt)
    Lt = np.array([LTR(np.array([v / KMH]), st60, lam)[0] for v in v_kmh])
    Lk = np.array([LTR_axle(np.array([v / KMH]), st60, lam)[0] for v in v_kmh])
    LOAD_MODEL = None
    return np.minimum(Lt, 1), Lk

def v_rollover_R60(model, srt=SRT_MODEL_A):
    """Speed above which there is no equilibrium in the curve (all inner tires lifted off)."""
    global LOAD_MODEL
    LOAD_MODEL = model
    lam = 1.0 if model is not None else lambda_factor(srt)
    if model is None:
        ok = lambda v: LTR(v, st60, lam) < 1.0
    else:
        ok = lambda v: LTR_axle(v, st60, lam).min(1) < 1.0
    v = v_max(ok, 1, 100 / KMH)[0] * KMH
    LOAD_MODEL = None
    return v

R60 = {'A': ltr_R60(None)}
R60.update({n: ltr_R60(m) for n, m in MODELS.items()})
explain("C.5  LTR versus speed in the R = 60 m curve of the ramp", """
    The vehicle travels through the accident curve (superelevation 7 %, grade -2.45 %) at constant speed.
    Left: whole-vehicle LTR [-] (criterion of Yoo et al.): Model A, Model B, Model C (circles) and
    TruckSim of Zhang et al. 2024 (squares).
    Right: LTR of the most loaded axle [-] (criterion of Model B, max LTR_i <= 0.8).
    When the curve reaches 1 the vehicle rolls over. Model A rolls over at about 66 km/h; Model C between 44 and 50 km/h.""")
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
ax[0].plot(v_kmh, R60['A'][0], 'k:', lw=2, label=f'Model A (imposed SRT {SRT_MODEL_A} g)')
for n, s in SUSP.items():
    ax[0].plot(v_kmh, R60[n][0], '-', color=s['color'], lw=2, label=f'Model B, {n}')
    ch = CHRONO_R60[n]
    ax[0].plot(ch['v'], ch['LTR'], 'o', color=s['color'], mfc='w', ms=7, label=f'Model C (PyChrono), {n}')
    ax[1].plot(v_kmh, R60[n][1].max(1), '-', color=s['color'], lw=2, label=f'Model B, {n}: max LTR_i')
ax[1].plot(v_kmh, R60['A'][1].max(1), 'k:', lw=2, label='Model A: max LTR_i')
ax[0].plot(TRUCKSIM_R60['v'], TRUCKSIM_R60['LTR'], 'ks', mfc='w', ms=7, label='TruckSim (Zhang et al., 2024)')
for a_ in ax:
    a_.axhline(LTR_MAX, color='gray', ls='--', lw=1); a_.text(26, LTR_MAX + 0.02, 'LTR_max = 0.8', fontsize=8, color='gray')
    a_.set_xlabel('speed in the R = 60 m curve [km/h]'); a_.grid(alpha=.3); a_.legend(fontsize=8); a_.set_ylim(0, 1.05)
ax[0].set_ylabel('whole-vehicle LTR [-]'); ax[0].set_title('Whole-vehicle LTR (criterion of Yoo et al., 2025)')
ax[1].set_ylabel('LTR of the most loaded axle [-]'); ax[1].set_title('Axle LTR (proposed criterion)')
plt.tight_layout(); plt.show()

rows = [['Model A', f'SRT {SRT_MODEL_A} imposed', np.interp(30, v_kmh, R60['A'][0]), np.interp(40, v_kmh, R60['A'][0]),
         np.interp(50, v_kmh, R60['A'][0]), v_rollover_R60(None)]]
for n, m in MODELS.items():
    rows.append([f'Model B, {n}', f'SRT {m.SRT:.3f} computed', np.interp(30, v_kmh, R60[n][0]),
                 np.interp(40, v_kmh, R60[n][0]), np.interp(50, v_kmh, R60[n][0]), v_rollover_R60(m)])
for n in SUSP:
    ch = CHRONO_R60[n]; lo, hi = CHRONO_V_ROLLOVER[n]
    rows.append([f'Model C (PyChrono), {n}', 'multibody', ch['LTR'][0], ch['LTR'][1],
                 ch['LTR'][4] if not ch['rollover'][4] else np.nan, f'{lo}-{hi}'])
rows.append(['TruckSim (Zhang et al., 2024)', '-', *TRUCKSIM_R60['LTR'], TRUCKSIM_R60['v_rollover']])
tab_r60 = pd.DataFrame(rows, columns=['model', 'SRT', 'LTR 30 km/h', 'LTR 40 km/h', 'LTR 50 km/h', 'rollover speed [km/h]'])
print(tab_r60.round(3).to_string(index=False))
print("(Model C already rolls over at 50 km/h in the 60 m curve; its mean LTR is omitted.)")

# ---- C.6 Safe-speed profile along the whole ramp (article Sec. 3.6) ---------
print("\n" + "=" * 78); print("C.6  SAFE-SPEED PROFILE v_final(s)  [Yoo Eqs. 13-19 with each load model]"); print("=" * 78)
def run_profile(model, mu0, c, srt=SRT_MODEL_A):
    """Final speed profile [m/s] and governing limit per station for a load model (None -> Model A)."""
    global LOAD_MODEL
    LOAD_MODEL = model
    vs, vr, vk, lam = safe_speed(ROAD, srt, mu0, c, N0_CALIB)
    vf = speed_profile(ROAD, vs, lam, mu0, c, N0_CALIB)
    LOAD_MODEL = None
    governs = np.where(vr < vk, 'rollover', 'skid')
    return vf, governs

CASES = [(0.3, 1.0, 'mu0 = 0.3, constant mu (case of Fig. 8a of Yoo et al., 2025)'),
         (0.3, 0.4, 'mu0 = 0.3, mu(N) with c = 0.4 (case of Fig. 8a of Yoo et al., 2025)'),
         (0.8, 1.0, 'mu0 = 0.8, dry pavement (as in the PyChrono runs)')]
explain("C.6  Safe-speed profile v_final(s) along the whole ramp", """
    Recommended speed at each station s [m] after the forward-backward sweep (Yoo Eqs. 13-19).
    One panel per friction case. Dotted: Model A (imposed SRT 0.5 g). Red/blue: Model B.
    Shaded bands in the 60 m curve: speed range at which Model C rolls over (red multi-leaf, blue air).
    A 'safe' speed inside or above the band is NOT safe.""")
fig, ax = plt.subplots(len(CASES), 1, figsize=(10, 11), sharex=True)
rows = []
for j, (mu0, c, title) in enumerate(CASES):
    vf, governs = run_profile(None, mu0, c)
    ax[j].plot(S, vf * KMH, 'k:', lw=2, label=f'Model A (SRT {SRT_MODEL_A} g)')
    row = [title, vf[(S > 364) & (S < 437)].min() * KMH, governs[i60]]
    for n, m in MODELS.items():
        vfn, governs_n = run_profile(m, mu0, c)
        ax[j].plot(S, vfn * KMH, '-', color=SUSP[n]['color'], lw=2, label=f'Model B, {n} (SRT {m.SRT:.3f} g)')
        row += [vfn[(S > 364) & (S < 437)].min() * KMH, governs_n[i60]]
    for n in SUSP:
        lo, hi = CHRONO_V_ROLLOVER[n]
        ax[j].fill_between([363.5, 437.3], lo, hi, color=SUSP[n]['color'], alpha=0.15)
    ax[j].text(366, CHRONO_V_ROLLOVER['air'][0] - 4, 'Model C rollover (bands)', fontsize=7)
    ax[j].axvspan(363.5, 437.3, color='gray', alpha=0.06)
    ax[j].set_title(title, fontsize=10); ax[j].set_ylabel('v [km/h]'); ax[j].grid(alpha=.3); ax[j].set_ylim(20, 85)
    ax[j].legend(fontsize=8, loc='lower left')
    rows.append(row)
ax[-1].set_xlabel('station [m]  (shaded: R = 60 m curve)')
plt.tight_layout(); plt.show()
tab_v = pd.DataFrame(rows, columns=['case', 'Model A v_min R60 [km/h]', 'Model A governs',
                                    'Model B multi-leaf v_min [km/h]', 'multi-leaf governs',
                                    'Model B air v_min [km/h]', 'air governs'])
print(tab_v.round(1).to_string(index=False))
print(f"Model C rolls over in the 60 m curve between {CHRONO_V_ROLLOVER['multileaf'][0]}-{CHRONO_V_ROLLOVER['multileaf'][1]} km/h "
      f"(multi-leaf) and {CHRONO_V_ROLLOVER['air'][0]}-{CHRONO_V_ROLLOVER['air'][1]} km/h (air).")

# ---- C.7 Sensitivity: roll centre, payload height and off-centre load (Sec. 3.3) ----
print("\n" + "=" * 78); print("C.7  SENSITIVITY OF THE PHYSICAL SRT (tilt table, peak a_y/g)"); print("=" * 78)
rows = []
for n, s in SUSP.items():
    base = CURVE[n]['SRT']
    s_lit = dict(s, H_RC=s['H_RC_LIT'])
    rows.append([n, 'baseline (roll centre of Model C)', base, 0.0])
    v_ = roll_curve(s_lit, tilt=True)['SRT']; rows.append([n, 'measured h_rc [Fu & Cebon, Table 2]', v_, 100 * (v_ / base - 1)])
    for dh in (-0.3, +0.3):                                   # payload-height variants of Kamnik et al., Table 1
        hs = (MS_TRACTOR * H_S_TRACTOR + MS_TRAILER * (H_S_TRAILER + dh)) / MS_TOTAL
        v_ = roll_curve(s, tilt=True, h_s=hs)['SRT']; rows.append([n, f'semitrailer CG {dh:+.1f} m [Kamnik et al.]', v_, 100 * (v_ / base - 1)])
    v_ = roll_curve(s, tilt=True, e_load=0.3)['SRT']; rows.append([n, 'payload 0.3 m off-centre [Kamnik et al.; Winkler & Ervin]', v_, 100 * (v_ / base - 1)])
tab_sens = pd.DataFrame(rows, columns=['suspension', 'variant', 'SRT [g]', 'change [%]'])
print(tab_sens.round(3).to_string(index=False))
print("References: Fu & Cebon (2002) - 10 % higher payload CG lowers the SRT by 12-14 %; "
      "ORNL (2009, track tests) - 0.06 g per foot (0.30 m) of payload height.")

# ---- C.7b Sensitivity to the limits of the measured roll-stiffness ranges (Sec. 3.3) ----
# Ranges per axle measured by UMTRI (Winkler, Karamihas & Bogard, 1992) / Fu & Cebon (2002) [N*m/rad]:
# drive axle and semitrailer axles at the lower (soft) and upper (stiff) limit. Axle 1 is unchanged.
print("\n" + "=" * 78); print("C.7b  SRT WITH THE LIMITS OF THE MEASURED K_phi RANGES"); print("=" * 78)
KPHI_RANGES = {'multileaf': {'drive': (315e3, 625e3), 'trailer': (460e3, 1030e3)},
               'air':       {'drive': (150e3, 330e3), 'trailer': (690e3, 1090e3)}}
rows = []
for n, s_ in SUSP.items():
    for i2, lbl2 in enumerate(('soft', 'stiff')):
        for i3, lbl3 in enumerate(('soft', 'stiff')):
            K = np.array(s_['K_PHI'], float)
            K[1] = KPHI_RANGES[n]['drive'][i2]; K[2] = K[3] = KPHI_RANGES[n]['trailer'][i3]
            v_ = roll_curve(dict(s_, K_PHI=K), tilt=True)['SRT']
            rows.append([n, f'drive axle {lbl2}, semitrailer axles {lbl3}', K[1] / 1e3, K[2] / 1e3, v_])
tab_ranges = pd.DataFrame(rows, columns=['suspension', 'case', 'drive-axle K_phi [kN*m/rad]',
                                         'semitrailer-axle K_phi [kN*m/rad]', 'SRT [g]'])
print(tab_ranges.round(3).to_string(index=False))

# ---- C.8 Comparison with track tests (Oak Ridge National Laboratory, 2009; Discussion) ----
# ORNL measured the lift-off of the semitrailer axles of a flatbed tractor-semitrailer with
# air suspension (about 36 t) with the PAYLOAD at two heights:
#   9.5 ft = 2.90 m  -> semitrailer lift-off at 0.30 g (step steer) / 0.33 g (constant radius)
#   8.5 ft = 2.59 m  -> semitrailer lift-off at 0.37 g (both tests)
# The 2.725 m of this study is the CG of the WHOLE sprung mass of the semitrailer (structure +
# payload). For the comparison the semitrailer is split into payload (M_PAYLOAD) and empty
# structure (MS_TRAILER - M_PAYLOAD); the payload is placed at the ORNL heights and the CG of
# the structure is assumed between 1.0 and 1.4 m (not reported). Only the CG height changes;
# everything else (masses, stiffnesses, tires) is the vehicle of this study.
#   h_sS = (m_e*h_e + m_p*h_p) / (m_e + m_p)      (article, Sec. 2.9)
print("\n" + "=" * 78); print("C.8  COMPARISON WITH MEASURED TRACK TESTS (ORNL, 2009)"); print("=" * 78)
M_STRUCT = MS_TRAILER - M_PAYLOAD                    # sprung structure of the empty semitrailer [kg]
H_STRUCT = (1.0, 1.2, 1.4)                           # assumed CG height of the structure [m]
H_PAYLOAD_STUDY = (MS_TRAILER * H_S_TRAILER - M_STRUCT * 1.2) / M_PAYLOAD
print(f"Payload CG height implied in this study (structure at 1.2 m): {H_PAYLOAD_STUDY:.2f} m")
ORNL = [('ORNL high CG (9.5 ft)', 2.90, '0.30-0.33'), ('ORNL low CG (8.5 ft)', 2.59, '0.37')]
rows = []
for n, s_ in SUSP.items():
    rows.append([n, 'this study', H_PAYLOAD_STUDY, H_S_TRAILER, H_S_TRAILER, CURVE[n]['SRT'], CURVE[n]['SRT'], '-'])
    for label, h_pay, measured in ORNL:
        h_trl = [(M_STRUCT * he + M_PAYLOAD * h_pay) / MS_TRAILER for he in H_STRUCT]
        srt = [roll_curve(s_, tilt=True, h_s=(MS_TRACTOR * H_S_TRACTOR + MS_TRAILER * hr) / MS_TOTAL)['SRT']
               for hr in h_trl]
        rows.append([n, label, h_pay, min(h_trl), max(h_trl), min(srt), max(srt), measured])
tab_ornl = pd.DataFrame(rows, columns=['suspension', 'case', 'payload CG [m]', 'semitrailer CG min [m]',
                                       'semitrailer CG max [m]', 'SRT min [g]', 'SRT max [g]', 'ORNL measured [g]'])
print(tab_ornl.round(3).to_string(index=False))
print("The ORNL semitrailer has an air suspension: compare with the 'air' rows.")



# #############################################################################
#
#  PART D - SENSITIVITY TO THE TIRE-ROAD FRICTION LAW (1987 MAGIC FORMULA)
#
#  Magic Formula (MF) of Bakker, Nyborg & Pacejka (1987), SAE 870421; see also
#  Pacejka (2006), Tire and Vehicle Dynamics, Ch. 4:
#
#     F_y(alpha) = D*sin( C*atan( B*alpha - E*(B*alpha - atan(B*alpha)) ) )
#     C   = 1.30                                   shape factor (lateral force)      [-]
#     D   = a1*F_z^2 + a2*F_z                      peak value                         [N]
#     BCD = a3*sin( a4*atan(a5*F_z) )              cornering stiffness (slope at 0)   [N/deg]
#     B   = BCD / (C*D)                            stiffness factor                   [1/deg]
#     E   = a6*F_z^2 + a7*F_z + a8                 curvature factor                   [-]
#     F_z in kN and alpha in degrees (convention of SAE 870421).
#
#  No coefficients a1-a8 of a TRUCK tire were available, therefore:
#   - D.2 (reported in the article, Secs. 2.7 and 3.7): the load-sensitive power law of
#     Yoo et al. (Yoo Eq. 6) is compared with the load sensitivity of the MF peak,
#     mu = D/F_z = a1*F_z + a2, BOTH fitted to the same truck-tire data (Segel & Ervin, as
#     reproduced in Fig. 5 of Yoo et al.). Only D is needed: the skid criterion of Yoo et al.
#     uses the peak friction (friction circle).
#   - D.3-D.4 (SUPPLEMENTARY, not reported in the article): the complete MF (D, BCD, C, E)
#     is fitted to the TMeasy tire curves of the multibody model (the same 295/80R22.5 and
#     275/70R22.5 tires) to estimate slip angles and friction utilization per axle with the
#     loads of Model B. Only the TMeasy load-sensitivity curve of D.3 is used in the article
#     (Fig. 9, "reference only").
# #############################################################################
from scipy.optimize import curve_fit, brentq

print("\n" + "#" * 78); print("PART D - SENSITIVITY TO THE FRICTION LAW (1987 MAGIC FORMULA)"); print("#" * 78)

C_MF = 1.30                                     # shape factor for F_y [-]  (SAE 870421)

def mf_fy(alpha_deg, Fz_kN, a, lam_mu=1.0):
    """Lateral force of ONE tire with the 1987 MF [N].
       alpha_deg: slip angle [deg];  Fz_kN: vertical load [kN];  a = (a1..a8);
       lam_mu: friction scaling factor lambda_mu (Pacejka, 2006, Ch. 4): scales the peak D, keeps BCD."""
    a1, a2, a3, a4, a5, a6, a7, a8 = a
    Fz = np.maximum(Fz_kN, 1e-3)
    D = lam_mu * (a1 * Fz**2 + a2 * Fz)                       # [N]
    BCD = a3 * np.sin(a4 * np.arctan(a5 * Fz))                # [N/deg]
    B = BCD / (C_MF * D)                                      # [1/deg]
    E = np.minimum(a6 * Fz**2 + a7 * Fz + a8, 1.0)            # [-]  (E <= 1 by definition of the MF)
    Ba = B * alpha_deg
    return D * np.sin(C_MF * np.arctan(Ba - E * (Ba - np.arctan(Ba))))

def mf_D(Fz_kN, a, lam_mu=1.0):
    """Peak D [N] for F_z [kN]."""
    return lam_mu * (a[0] * Fz_kN**2 + a[1] * Fz_kN)

def mf_BCD(Fz_kN, a):
    """Cornering stiffness BCD [N/deg] for F_z [kN]."""
    return a[2] * np.sin(a[3] * np.arctan(a[4] * Fz_kN))


# =============================================================================
# D.2  POWER LAW OF YOO ET AL. VERSUS MF PEAK, FITTED TO THE SAME DATA
# =============================================================================
print("\n" + "=" * 78); print("D.2  POWER LAW mu(N) OF YOO ET AL. vs LOAD SENSITIVITY OF THE MF PEAK"); print("=" * 78)
# Data: peak mu versus x = F_z/F_z0 of a truck tire (Segel & Ervin), digitized from Fig. 5 of Yoo et al.
#   F_z0 = N0 = 18 kN (the same reference load as Model A, N0_CALIB).
FZ0_KN = N0_CALIB / 1e3                                        # [kN]
x_SE, mu_SE = REF_FIG5_TRUCK[:, 0], REF_FIG5_TRUCK[:, 1]

# (a) Power law of Yoo et al. (Yoo Eq. 6):  mu(x) = mu0*x^(c-1)
def power_law(x, mu0, c):
    return mu0 * x ** (c - 1)
(mu0_pot, c_pot), _ = curve_fit(power_law, x_SE, mu_SE, p0=[0.75, 0.8])

# (b) 1987 MF: D = a1*F_z^2 + a2*F_z  ->  mu = D/(1000*F_z) = a1/1000*F_z + a2/1000   (linear in F_z)
#     in terms of x:  mu(x) = p + q*x,   with  a2 = 1000*p [N/kN],  a1 = 1000*q/F_z0 [N/kN^2]
(q_SE, p_SE) = np.polyfit(x_SE, mu_SE, 1)
MF_SE = (p_SE, q_SE)                                           # used by mu_N(..., c='MF') of Part A
a1_SE, a2_SE = 1000 * q_SE / FZ0_KN, 1000 * p_SE

def rms(y, yhat):
    """Root-mean-square error."""
    return np.sqrt(np.mean((y - yhat) ** 2))
print(f"Power law of Yoo et al. (Eq. 6):  mu0 = {mu0_pot:.3f}, c = {c_pot:.3f}  (source: 0.749 / 0.821);  "
      f"RMS error = {rms(mu_SE, power_law(x_SE, mu0_pot, c_pot)):.4f}")
print(f"1987 MF (peak D):  a1 = {a1_SE:.2f} N/kN^2, a2 = {a2_SE:.1f} N/kN  (with F_z0 = {FZ0_KN:.0f} kN);  "
      f"RMS error = {rms(mu_SE, p_SE + q_SE * x_SE):.4f}")
print("Both laws fit the 4 data points (x = 0.24-0.92) equally well; they differ when EXTRAPOLATED.")

# --- Plot D.2a -------------------------------------------------------------------
explain("D.2a  Peak friction versus normalized load: three laws fitted to the same data", """
    Peak tire mu versus x = F_z/F_z0 (F_z0 = 18 kN). Circles: Segel & Ervin data (Fig. 5 of Yoo et al.).
    Black: power law of Yoo et al., fitted (Eq. 6).  Green: 1987 MF, mu = D/F_z = a1*F_z + a2 (straight line).
    Grey dotted: power law with c = 0.4, the hypothetical case of Fig. 8a of Yoo et al.
    Purple: TMeasy tires of the multibody model (D/F_z of the MF fitted in D.3), normalized at x = 1.
    Grey band: range of tire loads that occurs on the ramp (from the unloaded inner tire to the most
    loaded outer tire). As the load tends to 0 the power law diverges and the MF does not.""")
x = np.linspace(0.05, 2.6, 200)
fig, ax = plt.subplots(figsize=(8, 5))
ax.axvspan(0.0, 2.2, color='gray', alpha=0.08)
ax.plot(x, power_law(x, mu0_pot, c_pot), 'k-', lw=2, label=f'Yoo et al. power law (μ0 = {mu0_pot:.3f}, c = {c_pot:.3f})')
ax.plot(x, np.clip(p_SE + q_SE * x, 0, None), 'g-', lw=2, label=f'1987 MF: a1 = {a1_SE:.2f} N/kN², a2 = {a2_SE:.0f} N/kN')
ax.plot(x, mu0_pot * x ** (0.4 - 1), ':', color='gray', lw=1.5, label='Yoo et al. power law, c = 0.4 (case of Fig. 8a)')
ax.plot(x_SE, mu_SE, 'ko', mfc='w', ms=8, label='Segel & Ervin (Fig. 5 of Yoo et al.)')
ax.set_xlabel('$x = F_z/F_{z0}$ [-]'); ax.set_ylabel('peak μ [-]'); ax.set_ylim(0, 1.6); ax.grid(alpha=.3)
AX_D2A = ax                     # the TMeasy curve is added in D.3, once it has been fitted


# =============================================================================
# D.3  MF FITTED TO THE TMeasy TIRES OF THE MULTIBODY MODEL
# =============================================================================
# TMeasy (tire model of PyChrono, template ChTMeasyTire; parameters in
# multibody_pychrono/vehiculo_json.py):
#   force versus lateral slip s (in free rolling s = tan(alpha))                    [Assumption]
#     0 <= s <= s_M :  F = dF0*s / (1 + sigma*(sigma + dF0*s_M/F_M - 2)),   sigma = s/s_M
#     s_M < s <= s_S:  F = F_M - (F_M - F_S)*sigma^2*(3 - 2*sigma),        sigma = (s - s_M)/(s_S - s_M)
#     s > s_S       :  F = F_S
#   dF0 [N] initial slope, F_M [N] peak, s_M [-] slip at the peak, F_S [N] sliding force.
#   Each parameter is given at the nominal load p_n and at 2*p_n; for another load:
#     forces:   Y(F_z) = (F_z/p_n)*[ 2*Y1 - Y2/2 - (Y1 - Y2/2)*F_z/p_n ]        (passes through Y1 and Y2)
#     s_M, s_S: y(F_z) = y1 + (y2 - y1)*(F_z/p_n - 1)
#   p_n = half of the load of the load index. The coefficients are the generic values of
#   vehiculo_json.py (NOT measurements of these tires): this comparison is model versus model.
TM_TIRES = {}
for label, LI in (('295/80R22.5 (axle 1)', 152), ('295/80R22.5 (axle 2)', 148), ('275/70R22.5 (axles 3-4)', 145)):
    pn = 0.5 * LI_LOAD[LI] * g                                               # nominal load [N]
    TM_TIRES[label] = dict(pn=pn,
        dF0=(7.4013 * pn, 6.8505 * 2 * pn), FM=(0.75876 * pn, 0.72628 * 2 * pn),
        FS=(0.68276 * pn, 0.65319 * 2 * pn), sM=(0.33167, 0.33216), sS=(1.0296, 1.0296))
AXLE_TIRE = ['295/80R22.5 (axle 1)', '295/80R22.5 (axle 2)', '275/70R22.5 (axles 3-4)', '275/70R22.5 (axles 3-4)']

def tm_param(L, key, Fz):
    """TMeasy parameter at load Fz [N], interpolated from its values at p_n and 2*p_n."""
    Y1, Y2 = L[key]; r = Fz / L['pn']
    if key in ('sM', 'sS'):
        return Y1 + (Y2 - Y1) * (r - 1)
    return r * (2 * Y1 - Y2 / 2 - (Y1 - Y2 / 2) * r)

def tm_fy(alpha_deg, Fz, L):
    """TMeasy lateral force of one tire [N] (road mu = 0.8, as in the PyChrono runs)."""
    s = np.abs(np.tan(np.radians(alpha_deg)))
    dF0, FM, FS = tm_param(L, 'dF0', Fz), tm_param(L, 'FM', Fz), tm_param(L, 'FS', Fz)
    sM, sS = tm_param(L, 'sM', Fz), tm_param(L, 'sS', Fz)
    sig1 = s / sM
    F1 = dF0 * s / (1 + sig1 * (sig1 + dF0 * sM / FM - 2))
    sig2 = np.clip((s - sM) / (sS - sM), 0, 1)
    F2 = FM - (FM - FS) * sig2**2 * (3 - 2 * sig2)
    return np.where(s <= sM, F1, np.where(s <= sS, F2, FS))

# Fit of the MF to each tire:
#   D:   TMeasy gives F_M(F_z) quadratic without constant term -> a1, a2 EXACT (same form as the MF).
#   BCD: a4 = 2 is fixed (original form of Bakker et al., 1987) and a3, a5 are fitted to
#        dF0(F_z)*pi/180 [N/deg] between 0.2*p_n and 2.5*p_n (least squares). With a4 free the fit
#        is not unique: the TMeasy curve is almost linear in F_z and a3 tends to infinity.
#   E:   a6, a7, a8 are fitted to the curves F_y(alpha), alpha = 0-20 deg, F_z = 0.5-2.5*p_n (least squares).
MF_TM = {}
rows = []
for label, L in TM_TIRES.items():
    pn_kN = L['pn'] / 1e3
    Y1, Y2 = L['FM']
    a2 = (2 * Y1 - Y2 / 2) / pn_kN                           # [N/kN]
    a1 = -(Y1 - Y2 / 2) / pn_kN**2                           # [N/kN^2]
    Fz_g = np.linspace(0.2, 2.5, 30) * L['pn']
    bcd_tm = tm_param(L, 'dF0', Fz_g) * np.pi / 180          # [N/deg]
    a4 = 2.0
    f_bcd = lambda Fz, a3, a5: a3 * np.sin(a4 * np.arctan(a5 * Fz))
    (a3, a5), _ = curve_fit(f_bcd, Fz_g / 1e3, bcd_tm, p0=[bcd_tm.max() * 1.2, 0.02], maxfev=20000)
    al_g, fz_g = np.meshgrid(np.linspace(0.5, 20, 25), np.array([0.5, 1, 1.5, 2, 2.5]) * L['pn'])
    al_g, fz_g = al_g.ravel(), fz_g.ravel()
    y_tm = tm_fy(al_g, fz_g, L)
    f_E = lambda X, a6, a7, a8: mf_fy(X[0], X[1] / 1e3, (a1, a2, a3, a4, a5, a6, a7, a8))
    (a6, a7, a8), _ = curve_fit(f_E, np.vstack([al_g, fz_g]), y_tm, p0=[0.0, 0.0, -0.5], maxfev=20000)
    MF_TM[label] = (a1, a2, a3, a4, a5, a6, a7, a8)
    err = rms(y_tm, f_E(np.vstack([al_g, fz_g]), a6, a7, a8)) / L['FM'][0] * 100
    rows.append([label, pn_kN, a1, a2, a3, a4, a5, a6, a7, a8, err])
tab_mf = pd.DataFrame(rows, columns=['tire', 'p_n [kN]', 'a1 [N/kN^2]', 'a2 [N/kN]', 'a3 [N/deg]', 'a4 [-]',
                                     'a5 [1/kN]', 'a6 [1/kN^2]', 'a7 [1/kN]', 'a8 [-]', 'RMS error [% of F_M]'])
print("\n" + "=" * 78); print("D.3  1987 MF COEFFICIENTS FITTED TO THE TMeasy TIRES OF THE MULTIBODY MODEL"); print("=" * 78)
print(tab_mf.round(4).to_string(index=False))
print("Shape reference for a truck tire (Lundahl et al. 2015, Table 2, simplified MF): C = 1.19, E = -1.21. Here C = 1.30.")

# Completes plot D.2a with the load sensitivity of the TMeasy tires
for label, col in zip(['295/80R22.5 (axle 2)', '275/70R22.5 (axles 3-4)'], ['purple', 'mediumorchid']):
    Fz = x * FZ0_KN
    mu_tm = mf_D(Fz, MF_TM[label]) / (1000 * Fz)
    mu_tm1 = mf_D(FZ0_KN, MF_TM[label]) / (1000 * FZ0_KN)
    AX_D2A.plot(x, mu_tm / mu_tm1 * mu0_pot, '-', color=col, lw=1.5, label=f'TMeasy, multibody model, {label} (scaled to μ0)')
AX_D2A.legend(fontsize=8); AX_D2A.set_title('Load sensitivity of friction: Yoo et al. vs 1987 MF vs TMeasy')
plt.tight_layout(); plt.show()
print("Drop of mu between x = 0.5 and x = 2: "
      f"power law {100*(power_law(2, mu0_pot, c_pot)/power_law(0.5, mu0_pot, c_pot)-1):.0f} %, "
      f"MF {100*((p_SE+2*q_SE)/(p_SE+0.5*q_SE)-1):.0f} %, "
      f"TMeasy (axle 2) {100*((mf_D(2*FZ0_KN, MF_TM[AXLE_TIRE[1]])/2)/(mf_D(0.5*FZ0_KN, MF_TM[AXLE_TIRE[1]])/0.5)-1):.0f} %")

# --- Plot D.3: quality of the MF -> TMeasy fit (supplementary) ----------------------
explain("D.3  Fit of the 1987 MF to the TMeasy tires of the multibody model (supplementary)", """
    Lateral force of ONE tire [kN] versus slip angle alpha [deg], at 0.5, 1 and 2 times its nominal load p_n.
    Solid line: TMeasy (tire model of PyChrono).  Dashed line: fitted 1987 MF (coefficients of table D.3).
    If they overlap, the MF represents the same tires as the multibody model and the slip angles are comparable.
    Note: the peak of these generic tires lies at alpha = 18 deg, higher than for a measured truck tire.""")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
al = np.linspace(0, 30, 200)
for j, (label, L) in enumerate(TM_TIRES.items()):
    for f, col in zip((0.5, 1.0, 2.0), ('tab:green', 'tab:blue', 'tab:red')):
        ax[j].plot(al, tm_fy(al, f * L['pn'], L) / 1e3, '-', color=col, lw=2, label=f'TMeasy, F_z = {f}·p_n')
        ax[j].plot(al, mf_fy(al, f * L['pn'] / 1e3, MF_TM[label]) / 1e3, '--', color='k', lw=1)
    ax[j].plot([], [], 'k--', label='fitted 1987 MF')
    ax[j].set_title(label); ax[j].set_xlabel('α [deg]'); ax[j].grid(alpha=.3); ax[j].legend(fontsize=8)
ax[0].set_ylabel('F_y per tire [kN]')
plt.tight_layout(); plt.show()


# --- D.2 (continued): effect of the friction law on the safe speed (article Sec. 3.7) ----
# Case of Fig. 8a of Yoo et al.: mu0 = 0.3 (slippery pavement), where skidding governs.
# Laws: 'constant mu' (c = 1) | 'power law, fitted c' (c = c_pot) | 'power law, c = 0.4' (case of the source) | '1987 MF' (c = 'MF')
LAWS = [('constant μ', 1.0, 'k-'), (f'Yoo et al. power law, c = {c_pot:.2f} (fitted)', c_pot, 'b--'),
        ('Yoo et al. power law, c = 0.4 (Fig. 8a)', 0.4, 'c:'), ('1987 MF (fitted to the same data)', 'MF', 'g-.')]
LOAD_MODELS = [('Model A (SRT 0.5 g)', None), ('Model B, multi-leaf', MODELS['multileaf']),
               ('Model B, air', MODELS['air'])]
MU0_D2 = 0.3
rows, PROFILES_D2 = [], {}
for label_m, model in LOAD_MODELS:
    for label_l, c, _ in LAWS:
        vf, governs = run_profile(model, MU0_D2, c)
        PROFILES_D2[(label_m, label_l)] = vf
        rows.append([label_m, label_l, vf[(S > 364) & (S < 437)].min() * KMH, vf[S <= 600].mean() * KMH, governs[i60]])
tab_d2 = pd.DataFrame(rows, columns=['load model', 'friction law', 'v_min in R60 [km/h]',
                                     'mean v 0-600 m [km/h]', 'governs in R60'])
explain("D.2b  Safe speed along the ramp with each friction law (mu0 = 0.3)", """
    Profile v_final(s) [km/h] versus station s [m], as in Fig. 8a of Yoo et al., with four friction laws.
    One panel per load model: Model A, Model B multi-leaf and Model B air.
    Black: constant mu.  Blue: power law fitted to Segel & Ervin.  Cyan: c = 0.4 (case of the source).
    Green: 1987 MF fitted to the same data. Shaded: 60 m curve.
    How much the safe speed decreases when mu falls with load, depending on the law that is used.""")
fig, ax = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
for j, (label_m, _) in enumerate(LOAD_MODELS):
    for label_l, c, style in LAWS:
        ax[j].plot(S, PROFILES_D2[(label_m, label_l)] * KMH, style, lw=2, label=label_l)
    ax[j].axvspan(363.5, 437.3, color='gray', alpha=0.08)
    ax[j].set_title(label_m, fontsize=10); ax[j].set_ylabel('v [km/h]'); ax[j].grid(alpha=.3)
    ax[j].set_ylim(25, 85); ax[j].legend(fontsize=8, loc='lower left')
ax[-1].set_xlabel('station [m]')
plt.tight_layout(); plt.show()
print(tab_d2.round(1).to_string(index=False))
for label_m, _ in LOAD_MODELS:
    base = PROFILES_D2[(label_m, LAWS[0][0])][S <= 600].mean()
    difs = [(base - PROFILES_D2[(label_m, l)][S <= 600].mean()) * KMH for l, _, _ in LAWS[1:]]
    print(f"{label_m:26s}: mean decrease relative to constant mu -> fitted c {difs[0]:.1f}, "
          f"c = 0.4 {difs[1]:.1f}, MF {difs[2]:.1f} km/h  (Yoo et al. report about 7.5 km/h with c = 0.4)")


# =============================================================================
# D.4  SUPPLEMENTARY: SLIP ANGLE AND FRICTION UTILIZATION PER AXLE (not reported in the article)
# =============================================================================
print("\n" + "=" * 78); print("D.4  SUPPLEMENTARY: SLIP ANGLE AND FRICTION UTILIZATION PER AXLE (MF fitted to TMeasy)"); print("=" * 78)
# Assumptions:
#  - Steady-state cornering: each axle provides a lateral force proportional to the load it carries,
#      F_y,i^dem = m_i*a_y   [N]     (Pacejka 2006, Ch. 1, steady-state single-track model; the same
#    idea as Yoo Eq. 8 but at axle level). In a real tandem the two axles do not share equally
#    (tandem scrub): the multibody model shows it; here the tandem average is compared.  [Assumption]
#  - Outer and inner tires with their actual load (Model B or Model A):
#      F_y,i(alpha) = n_i*[ F_y(alpha, N_ext,i/n_i) + F_y(alpha, N_int,i/n_i) ]     (1987 MF, D.3)
#  - alpha_i = slip angle at which F_y,i(alpha_i) = F_y,i^dem   (solved with brentq).
#  - Friction utilization of the axle: eta_i = F_y,i^dem / ( n_i*[D(N_ext/n) + D(N_int/n)] ) [-]; eta_i = 1 -> axle saturates.
#  - Effective cornering stiffness of the axle (Pacejka 2006, Ch. 1): C_eff,i = n_i*[BCD(N_ext/n) + BCD(N_int/n)];
#    since BCD is concave in F_z, load transfer reduces it relative to 2*n_i*BCD(N_stat/2n_i).
#  - Other surface: the peak is scaled with lambda_mu = mu_road / 0.8 (the TMeasy tires of the
#    multibody model are defined for mu = 0.8).

def level_loads(ayg, model):
    """N_ext, N_int [N] of each axle on level ground for a_y/g; model = PhysicalSRTModel, 'A' or 'static'."""
    W = M_AXLE * g
    if model == 'A':
        L = ltr_model_a_level(np.atleast_1d(ayg), SRT_MODEL_A)[0][0]
    elif model == 'static':
        L = np.zeros(4)
    else:
        L = model.ltr(np.atleast_1d(ayg))[0]
    return W / 2 * (1 + L), W / 2 * (1 - L)

def axle_fy(alpha, k, Ne, Ni, lam_mu=1.0):
    """Lateral force of axle k [N] at slip angle alpha [deg] with outer/inner side loads Ne, Ni [N]."""
    a, n = MF_TM[AXLE_TIRE[k]], N_TIRES_SIDE[k]
    return n * (mf_fy(alpha, Ne / n / 1e3, a, lam_mu) + mf_fy(alpha, Ni / n / 1e3, a, lam_mu))

def axle_state(ayg, model, lam_mu=1.0):
    """For a given a_y/g: alpha_i [deg] (NaN if the axle saturates), eta_i [-] and C_eff,i/C_eff,i(no transfer) [-]."""
    srt = SRT_MODEL_A if model == 'A' else (np.inf if model == 'static' else model.SRT)
    if ayg >= srt:                                   # beyond the SRT there is no equilibrium: the vehicle has rolled over
        return np.full(4, np.nan), np.full(4, np.nan), np.full(4, np.nan)
    Ne, Ni = level_loads(ayg, model)
    alpha, eta, ceff = np.full(4, np.nan), np.zeros(4), np.zeros(4)
    for k in range(4):
        a, n = MF_TM[AXLE_TIRE[k]], N_TIRES_SIDE[k]
        dem = M_AXLE[k] * ayg * g
        cap = n * (mf_D(Ne[k] / n / 1e3, a, lam_mu) + mf_D(max(Ni[k], 0) / n / 1e3, a, lam_mu))
        eta[k] = dem / cap
        ceff[k] = n * (mf_BCD(Ne[k] / n / 1e3, a) + mf_BCD(max(Ni[k], 1) / n / 1e3, a)) / \
                  (2 * n * mf_BCD(M_AXLE[k] * g / 2 / n / 1e3, a))
        al_g = np.linspace(0, 40, 801)
        F_g = axle_fy(al_g, k, Ne[k], max(Ni[k], 0), lam_mu)
        i_max = int(np.argmax(F_g))
        if dem <= F_g[i_max] and dem > 0:
            alpha[k] = brentq(lambda x_: axle_fy(x_, k, Ne[k], max(Ni[k], 0), lam_mu) - dem, 0, al_g[i_max])
        elif dem == 0:
            alpha[k] = 0.0
    return alpha, eta, ceff

# --- Model C data: mean slip angle per axle in the constant-radius curves -----------
#   (table tabla_rampa_estacionario.csv of the multibody model; runs without rollover only).
#   a_y,eff/g = a_y,geom/g*cos(beta) - sin(beta).
#   columns: radius [m], a_y,geom/g, alpha1, alpha2, alpha3, alpha4 [deg]
CHRONO_SLIP = {
 'multileaf': np.array([
   [125, 0.102, 0.276, 0.178, -0.076, 0.544], [125, 0.123, 0.393, 0.356, 0.079, 0.691], [125, 0.141, 0.496, 0.532, 0.202, 0.791],
   [125, 0.160, 0.650, 0.699, 0.372, 0.992], [125, 0.197, 1.213, 1.066, 0.764, 1.380], [125, 0.241, 1.721, 1.539, 1.303, 1.924],
   [60, 0.118, 0.393, 0.260, -0.379, 0.927], [60, 0.210, 0.927, 1.053, 0.392, 1.695], [60, 0.254, 1.606, 1.508, 0.878, 2.211]]),
 'air': np.array([
   [125, 0.102, 0.287, 0.193, -0.066, 0.549], [125, 0.123, 0.380, 0.349, 0.093, 0.692], [125, 0.141, 0.490, 0.539, 0.199, 0.784],
   [125, 0.160, 0.672, 0.717, 0.362, 0.977], [125, 0.197, 1.208, 1.065, 0.778, 1.392], [125, 0.241, 1.709, 1.521, 1.336, 1.944],
   [60, 0.118, 0.407, 0.269, -0.371, 0.973], [60, 0.210, 0.981, 1.049, 0.420, 1.708], [60, 0.254, 1.615, 1.498, 0.952, 2.260]]),
}
for n_, A in CHRONO_SLIP.items():
    beta_ = np.arctan(np.where(A[:, 0] == 125, 0.06, 0.07))
    CHRONO_SLIP[n_] = np.column_stack([np.cos(beta_) * A[:, 1] - np.sin(beta_), A[:, 2:]])   # [a_y,eff/g, alpha1..alpha4]

ay_d = np.linspace(0.005, 0.45, 90)
STATE = {}
for n_, m_ in list(MODELS.items()) + [('A', 'A'), ('no transfer', 'static')]:
    STATE[n_] = {lam: np.array([np.concatenate(axle_state(a_, m_, lam)) for a_ in ay_d]) for lam in (1.0, 0.3 / 0.8)}

explain("D.4a  Slip angle per axle: MF + Model B versus Model C (dry road, mu = 0.8)", """
    Slip angle alpha_i [deg] that each axle needs to provide its lateral force m_i*a_y, versus a_y,eff/g.
    Left panel: tractor (axle 1 steer, axle 2 drive). Centre panel: semitrailer tandem (mean of axles 3-4).
    Lines: model (MF fitted to TMeasy + loads of Model B); red multi-leaf, blue air.  Grey: no load transfer.
    Markers: Model C in the 125 m and 60 m curves (circle = multi-leaf, square = air). In the tandem the bar
    joins alpha3 and alpha4 of Model C: the two axles do not slip equally (tandem scrub), an effect
    that this model does not include.
    Right panel: alpha1 - alpha2 [deg] of the tractor. Positive = understeer (delta = L/R + alpha1 - alpha2).""")
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for n_, col in (('multileaf', 'tab:red'), ('air', 'tab:blue')):
    E = STATE[n_][1.0]
    ax[0].plot(ay_d, E[:, 0], '-', color=col, lw=2, label=f'axle 1, {n_}')
    ax[0].plot(ay_d, E[:, 1], '--', color=col, lw=2, label=f'axle 2, {n_}')
    ax[1].plot(ay_d, (E[:, 2] + E[:, 3]) / 2, '-', color=col, lw=2, label=f'tandem, {n_}')
    ax[2].plot(ay_d, E[:, 0] - E[:, 1], '-', color=col, lw=2, label=f'model, {n_}')
    A = CHRONO_SLIP[n_]; mk = 'o' if n_ == 'multileaf' else 's'
    ax[0].plot(A[:, 0], A[:, 1], mk, color=col, mfc='w', ms=6, label=f'Model C axle 1, {n_}')
    ax[0].plot(A[:, 0], A[:, 2], mk, color=col, ms=5, alpha=0.6, label=f'Model C axle 2, {n_}')
    ax[1].plot(A[:, 0], (A[:, 3] + A[:, 4]) / 2, mk, color=col, mfc='w', ms=6, label=f'Model C mean of axles 3-4, {n_}')
    ax[1].vlines(A[:, 0], A[:, 3], A[:, 4], color=col, alpha=0.3)
    ax[2].plot(A[:, 0], A[:, 1] - A[:, 2], mk, color=col, mfc='w', ms=6, label=f'Model C, {n_}')
E0 = STATE['no transfer'][1.0]
ax[0].plot(ay_d, E0[:, 0], '-', color='gray', lw=1, label='axle 1, no load transfer')
ax[0].plot(ay_d, E0[:, 1], '--', color='gray', lw=1, label='axle 2, no load transfer')
ax[1].plot(ay_d, (E0[:, 2] + E0[:, 3]) / 2, '-', color='gray', lw=1, label='no load transfer')
ax[2].plot(ay_d, E0[:, 0] - E0[:, 1], '-', color='gray', lw=1, label='no load transfer')
for a_, t_ in zip(ax, ('Tractor: axles 1 and 2', 'Semitrailer tandem (axles 3-4)', 'Tractor: α1 − α2 (understeer)')):
    a_.set_xlabel('$a_{y,eff}/g$ [-]'); a_.grid(alpha=.3); a_.set_title(t_); a_.legend(fontsize=7); a_.set_xlim(0, 0.25)
ax[0].set_ylabel('α [deg]'); ax[1].set_ylabel('α [deg]'); ax[2].set_ylabel('α1 − α2 [deg]')
ax[0].set_ylim(0, 3); ax[1].set_ylim(-0.5, 3)
plt.tight_layout(); plt.show()

# Error of the model against Model C (mean |d alpha| over the points of each axle / tandem)
for n_ in SUSP:
    A = CHRONO_SLIP[n_]
    mod = np.array([axle_state(a_, MODELS[n_], 1.0)[0] for a_ in A[:, 0]])
    e1, e2 = np.nanmean(np.abs(mod[:, 0] - A[:, 1])), np.nanmean(np.abs(mod[:, 1] - A[:, 2]))
    e34 = np.nanmean(np.abs((mod[:, 2] + mod[:, 3]) / 2 - (A[:, 3] + A[:, 4]) / 2))
    print(f"{n_:10s}: mean |alpha model - alpha Model C|: axle 1 {e1:.2f} deg, axle 2 {e2:.2f} deg, tandem {e34:.2f} deg")

# Criteria with margin, as in Yoo et al.: rollover when max LTR_i = LTR_max = 0.8 (Model A: whole-vehicle LTR = 0.8);
# skid of axle i when it uses (1 - delta) = 80 % of its friction: eta_i/(1 - delta) = 1 (Yoo Eq. 7 applied per axle).
def ay_ltr_max(n_):
    """a_y/g at which LTR_max is reached (level ground)."""
    if n_ == 'A':
        return LTR_MAX * SRT_MODEL_A                           # whole-vehicle LTR is linear in a_y and = 1 at the SRT
    m_ = MODELS[n_]
    return np.interp(LTR_MAX, m_.LTR.max(1), m_.ay)

explain("D.4b  Friction utilization per axle: which axle saturates first (mu = 0.8 and mu = 0.3)", """
    eta_i = lateral force required by the axle / maximum force its tires can provide with their actual load [-].
    Top: dry road (mu = 0.8). Bottom: slippery road (mu = 0.3, lambda_mu = 0.375). One color per axle.
    Columns: Model A (loads with SRT 0.5 g), Model B multi-leaf and Model B air.
    The curves end at the SRT of each model: beyond it there is no equilibrium (the vehicle rolls over).
    Grey dashed line: eta = 1 - delta = 0.8, the design limit of Yoo et al. (Eq. 7).  Black line: a_y at which
    some axle reaches LTR = 0.8 (rollover limit with margin). Whichever is reached first defines the failure mode:
    axle 1 -> the tractor does not turn (understeer); axle 2 -> jackknife tendency; axles 3-4 -> trailer swing;
    LTR -> rollover.""")
fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex=True, sharey=True)
SRT_MOD = {'A': SRT_MODEL_A, 'multileaf': MODELS['multileaf'].SRT, 'air': MODELS['air'].SRT}
rows, rows_eta = [], []
i20 = int(np.argmin(np.abs(ay_d - 0.20)))
for i_, lam in enumerate((1.0, 0.3 / 0.8)):
    for j_, n_ in enumerate(('A', 'multileaf', 'air')):
        E = STATE[n_][lam]
        ay_v = ay_ltr_max(n_)
        ay_sat = []
        for k in range(4):
            ax[i_, j_].plot(ay_d, E[:, 4 + k], '-', color=AXLE_COLORS[k], lw=2, label=f'axle {k+1}')
            eta_d = E[:, 4 + k] / (1 - DELTA)
            ok_ = np.isfinite(eta_d)
            ay_sat.append(np.interp(1.0, eta_d[ok_], ay_d[ok_]) if (eta_d[ok_] >= 1.0).any() else np.nan)
        ax[i_, j_].axhline(1 - DELTA, color='gray', ls='--', lw=1); ax[i_, j_].axvline(ay_v, color='k', lw=1)
        ax[i_, j_].set_title(f"{'Model A' if n_ == 'A' else 'Model B, ' + n_} — μ = {0.8*lam:.1f}", fontsize=10)
        ax[i_, j_].grid(alpha=.3); ax[i_, j_].set_ylim(0, 1.3)
        ay_sat = np.array(ay_sat)
        if np.isfinite(ay_sat).any() and np.nanmin(ay_sat) < ay_v:
            k1 = int(np.nanargmin(ay_sat)); mode = f"skid of axle {k1 + 1} at {ay_sat[k1]:.3f} g"
        else:
            mode = f"rollover (LTR 0.8) at {ay_v:.3f} g"
        rows.append([f"mu = {0.8 * lam:.1f}", n_, ay_v] + list(ay_sat) + [mode])
        rows_eta.append([f"mu = {0.8 * lam:.1f}", n_] + list(E[i20, 4:8]))
ax[0, 0].legend(fontsize=8)
for a_ in ax[1]: a_.set_xlabel('$a_y/g$ [-]')
for a_ in ax[:, 0]: a_.set_ylabel('friction utilization η_i [-]')
plt.tight_layout(); plt.show()
tab_d4 = pd.DataFrame(rows, columns=['road', 'load model', 'a_y at LTR = 0.8 [g]', 'a_y axle 1 at 80 % [g]',
                                     'a_y axle 2 at 80 % [g]', 'a_y axle 3 at 80 % [g]', 'a_y axle 4 at 80 % [g]', 'first limit'])
print(tab_d4.round(3).to_string(index=False))
tab_eta = pd.DataFrame(rows_eta, columns=['road', 'load model', 'eta axle 1', 'eta axle 2', 'eta axle 3', 'eta axle 4'])
print("\nFriction utilization eta_i at a_y = 0.20 g (which axle is closest to saturation):")
print(tab_eta.round(3).to_string(index=False))

# Effective cornering stiffness per axle at 0.15 g (Pacejka 2006, Ch. 1: load transfer reduces it)
i15 = int(np.argmin(np.abs(ay_d - 0.15)))
rows = [[n_] + list(STATE[n_][1.0][i15, 8:12]) for n_ in ('A', 'multileaf', 'air')]
tab_ceff = pd.DataFrame(rows, columns=['load model', 'C_eff axle 1 [-]', 'C_eff axle 2 [-]', 'C_eff axle 3 [-]', 'C_eff axle 4 [-]'])
print("\nEffective cornering stiffness of the axle at a_y = 0.15 g, relative to the axle without load transfer:")
print(tab_ceff.round(3).to_string(index=False))


# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "#" * 78); print("SUMMARY"); print("#" * 78)
print(f"SRT: Model A {SRT_MODEL_A} g (imposed) | Model B multi-leaf {CURVE['multileaf']['SRT']:.3f} g, air "
      f"{CURVE['air']['SRT']:.3f} g | Model C first lift-off {CHRONO_TILT['multileaf']['SRT_first']:.3f} / "
      f"{CHRONO_TILT['air']['SRT_first']:.3f} g")
print(f"Rollover in R60: Model A {v_rollover_R60(None):.1f} km/h | Model B {v_rollover_R60(MODELS['multileaf']):.1f} / "
      f"{v_rollover_R60(MODELS['air']):.1f} km/h | Model C 47.1-49.9 / 44.0-46.9 km/h | TruckSim 54.9 km/h")
print("Friction law: see table D.2.  Supplementary failure mode per axle: see table D.4.")
