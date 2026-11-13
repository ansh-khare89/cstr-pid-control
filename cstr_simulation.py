#!/usr/bin/env python3
"""
cstr_simulation.py
==================
Non-isothermal CSTR with PID temperature control (Ziegler-Nichols tuned)

Process
-------
Exothermic, first-order, irreversible liquid-phase reaction  A -> B  in a
jacketed, perfectly mixed CSTR.  Cooling water flows through a well-mixed
jacket.  The reactor temperature is the controlled variable (CV), the coolant
flow rate is the manipulated variable (MV), and the coolant inlet temperature
is the disturbance variable (DV).

Study
-----
A +10 K step in coolant inlet temperature at t = 10 min is applied to
  (a) the open-loop reactor (coolant flow held constant), and
  (b) the closed-loop reactor under a discrete PID controller whose gains
      are obtained from the Ziegler-Nichols *ultimate-gain* (continuous
      cycling) rules.  Ku and Pu are computed analytically from the
      linearised plant (no trial-and-error), so the tuning is reproducible.

Outputs
-------
  * Performance metrics printed to the terminal
  * cstr_temperature_response.png   (300 dpi)
  * cstr_concentration_control.png  (300 dpi)

Run:  python cstr_simulation.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")                     # headless backend (works on servers/CI)
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, brentq

# =============================================================================
# 1. PROCESS PARAMETERS  (SI-friendly "plant" units: L, min, K, mol, J)
# =============================================================================
R_GAS = 8.314            # J/(mol.K)   universal gas constant

# --- Reactor / feed ---------------------------------------------------------
V      = 500.0           # L           reactor liquid volume
q      = 50.0            # L/min       feed (and product) volumetric flow
CAf    = 2.0             # mol/L       feed concentration of A
Tf     = 320.0           # K           feed temperature
rho    = 1000.0          # g/L         reaction-mixture density
Cp     = 4.18            # J/(g.K)     reaction-mixture heat capacity

# --- Reaction kinetics:  k = k0 * exp(-E / (R T)) ---------------------------
E_ACT  = 7.0e4           # J/mol       activation energy
dHr    = -1.5e5          # J/mol       heat of reaction (exothermic < 0)

# --- Cooling jacket ---------------------------------------------------------
Vj     = 50.0            # L           jacket liquid holdup
rho_c  = 1000.0          # g/L         coolant density
Cp_c   = 4.18            # J/(g.K)     coolant heat capacity
Tc_in0 = 290.0           # K           nominal coolant inlet temperature

# --- Instrumentation --------------------------------------------------------
TAU_SENSOR = 0.2         # min         RTD / thermowell first-order lag
TS         = 0.05        # min         controller sample time (3 s)
QC_MIN, QC_MAX_FACTOR = 0.0, 2.0   # valve sized ~50 % open at nominal duty (0 .. 2x nominal)

# --- Design operating point (used to back-calculate k0 and nominal qc) ------
T_SP   = 350.0           # K           reactor temperature setpoint
X_DES  = 0.75            # -           design conversion at the setpoint
T_JAC_DES = 320.0        # K           design jacket temperature

# --- Scenario ---------------------------------------------------------------
T_DIST      = 10.0       # min         time of disturbance
DT_DIST     = 10.0       # K           magnitude of coolant-inlet step
T_END       = 80.0       # min         simulation horizon
T_SAFE      = 356.0      # K           high-temperature safety alarm / trip limit

# Derived quantities: pre-exponential factor, jacket UA and nominal coolant flow
# are back-calculated so the plant sits exactly at (T_SP, X_DES, T_JAC_DES) at
# steady state (a standard way of 'designing' the nominal operating point).
_k_req = (q / V) * X_DES / (1.0 - X_DES)                    # 1/min at T_SP
K0     = _k_req * np.exp(E_ACT / (R_GAS * T_SP))            # 1/min
_Qrem  = (-dHr) * _k_req * CAf * (1 - X_DES) * V \
         + q * rho * Cp * (Tf - T_SP)                       # J/min heat to remove
UA     = _Qrem / (T_SP - T_JAC_DES)                         # J/(min.K)
QC_NOM = _Qrem / (rho_c * Cp_c * (T_JAC_DES - Tc_in0))      # L/min
QC_MAX = QC_MAX_FACTOR * QC_NOM


# =============================================================================
# 2. PLANT MODEL  (mass + energy balances, reactor + jacket + sensor)
# =============================================================================
def arrhenius(T):
    """First-order rate constant k(T) [1/min]."""
    return K0 * np.exp(-E_ACT / (R_GAS * T))


def cstr_odes(t, x, qc, Tc_in):
    """
    State vector x = [CA, T, Tj, Tm]
        CA : reactor concentration of A      [mol/L]
        T  : reactor temperature             [K]
        Tj : jacket (coolant) temperature    [K]
        Tm : measured temperature (sensor)   [K]

    Mass balance on A        :  dCA/dt = q/V (CAf - CA) - k CA
    Reactor energy balance   :  dT/dt  = q/V (Tf - T) + (-dHr) k CA /(rho Cp)
                                          - UA /(V rho Cp) (T - Tj)
    Jacket energy balance    :  dTj/dt = qc/Vj (Tc_in - Tj)
                                          + UA /(Vj rho_c Cp_c) (T - Tj)
    Sensor dynamics          :  dTm/dt = (T - Tm) / tau_s
    """
    CA, T, Tj, Tm = x
    k = arrhenius(T)
    dCA = q / V * (CAf - CA) - k * CA
    dT = (q / V * (Tf - T)
          + (-dHr) * k * CA / (rho * Cp)
          - UA / (V * rho * Cp) * (T - Tj))
    dTj = qc / Vj * (Tc_in - Tj) + UA / (Vj * rho_c * Cp_c) * (T - Tj)
    dTm = (T - Tm) / TAU_SENSOR
    return [dCA, dT, dTj, dTm]


def steady_state():
    """Solve the nominal steady state (should reproduce T_SP, CA = CAf(1-X))."""
    f = lambda z: cstr_odes(0, [z[0], z[1], z[2], z[1]], QC_NOM, Tc_in0)[:3]
    CA, T, Tj = fsolve(f, [CAf * (1 - X_DES), T_SP, T_JAC_DES], xtol=1e-12)
    return np.array([CA, T, Tj, T])


# =============================================================================
# 3. ZIEGLER-NICHOLS TUNING (ultimate-gain method from the linearised plant)
# =============================================================================
def linearise(x_ss):
    """Numerical Jacobians A = df/dx, B = df/dqc about the steady state."""
    eps = 1e-6
    A = np.zeros((4, 4))
    for i in range(4):
        d = np.zeros(4); d[i] = eps
        A[:, i] = (np.array(cstr_odes(0, x_ss + d, QC_NOM, Tc_in0))
                   - np.array(cstr_odes(0, x_ss - d, QC_NOM, Tc_in0))) / (2 * eps)
    B = ((np.array(cstr_odes(0, x_ss, QC_NOM + eps, Tc_in0))
          - np.array(cstr_odes(0, x_ss, QC_NOM - eps, Tc_in0))) / (2 * eps))
    C = np.array([0.0, 0.0, 0.0, 1.0])          # measured output = Tm
    return A, B.reshape(-1, 1), C


def ziegler_nichols_pid(x_ss):
    """
    Classical Z-N ultimate-gain rules:
        Kp = 0.6 Ku ,  tau_i = Pu / 2 ,  tau_d = Pu / 8
    Ku (ultimate gain) and Pu (ultimate period) are found from the frequency
    response L(jw) = -C (jwI - A)^-1 B of the plant (sign flipped because the
    controller is reverse-acting: more coolant -> lower temperature).
    The ultimate frequency is where phase(L) = -180 deg.
    """
    A, B, C = linearise(x_ss)
    I = np.eye(4)

    def L(w):
        return -(C @ np.linalg.solve(1j * w * I - A, B)).item()

    def phase_unwrapped(w):
        # Sweep from low frequency to unwrap the phase continuously.
        ws = np.logspace(-4, np.log10(w), 4000)
        ph = np.unwrap([np.angle(L(wi)) for wi in ws])
        return ph[-1]

    # Bracket the -180 deg crossing on a log grid, then refine with brentq.
    grid = np.logspace(-2, 2, 800)
    ph = np.unwrap([np.angle(L(w)) for w in grid]) + np.pi
    idx = np.where(np.diff(np.sign(ph)) != 0)[0][0]
    w_u = brentq(lambda w: phase_unwrapped(w) + np.pi, grid[idx], grid[idx + 1])

    Ku = 1.0 / abs(L(w_u))                    # [L/min per K]
    Pu = 2.0 * np.pi / w_u                    # [min]
    return dict(Ku=Ku, Pu=Pu, Kp=0.6 * Ku, tau_i=Pu / 2.0, tau_d=Pu / 8.0)


# =============================================================================
# 4. SIMULATION ENGINE
# =============================================================================
def simulate(x0, controller=None):
    """
    Integrate the plant with solve_ivp over successive sample intervals.
    The coolant flow is held constant (zero-order hold) within each interval,
    which is exactly how a digital controller drives a control valve.

    controller = None   -> open loop (qc fixed at QC_NOM)
    controller = dict   -> discrete PID with derivative filter and
                           anti-windup (conditional integration)
    """
    n = int(round(T_END / TS))
    t_log = np.zeros(n + 1)
    X = np.zeros((n + 1, 4))
    U = np.zeros(n + 1)
    X[0] = x0
    qc = QC_NOM
    integral, d_filt, e_prev = 0.0, 0.0, 0.0
    N_FILT = 10.0                              # derivative filter coefficient

    for i in range(n):
        t = i * TS
        Tc_in = Tc_in0 + (DT_DIST if t >= T_DIST else 0.0)

        if controller is not None:
            # Reverse-acting: error = measurement - setpoint
            e = X[i, 3] - T_SP
            Kp, ti, td = controller["Kp"], controller["tau_i"], controller["tau_d"]
            # Filtered derivative (first-order, Tf = td/N)
            alpha = (td / N_FILT) / (td / N_FILT + TS)
            d_raw = (e - e_prev) / TS
            d_filt = alpha * d_filt + (1 - alpha) * d_raw
            e_prev = e
            u_unsat = QC_NOM + Kp * (e + integral / ti + td * d_filt)
            qc = float(np.clip(u_unsat, QC_MIN, QC_MAX))
            # Anti-windup: only integrate when the valve is not saturated
            if u_unsat == qc:
                integral += e * TS
        U[i] = qc

        sol = solve_ivp(cstr_odes, (t, t + TS), X[i], args=(qc, Tc_in),
                        method="LSODA", rtol=1e-8, atol=1e-10)
        X[i + 1] = sol.y[:, -1]
        t_log[i + 1] = t + TS
    U[-1] = U[-2]
    return t_log, X, U


# =============================================================================
# 5. PERFORMANCE METRICS
# =============================================================================
def metrics(t, T, label, band_K=0.1):
    """
    Peak overshoot : max(T) - T_SP after the disturbance  [K and % of T_SP]
    Steady-state error : T(end) - T_SP                      [K]
    Settling time  : time after the disturbance after which T stays inside
                     +/- band_K of its final value          [min]
    """
    post = t >= T_DIST
    peak = T[post].max()
    over_K = max(peak - T_SP, 0.0)
    over_pct = 100.0 * over_K / T_SP
    sse = T[-1] - T_SP
    outside = np.where(np.abs(T[post] - T[-1]) > band_K)[0]
    ts = (t[post][outside[-1]] + TS - T_DIST) if outside.size else 0.0
    return dict(label=label, peak_K=peak, over_K=over_K, over_pct=over_pct,
                settle=ts, sse=sse)


# =============================================================================
# 6. MAIN
# =============================================================================
def main():
    x_ss = steady_state()
    print("=" * 66)
    print(" CSTR PID CONTROL  -  Ziegler-Nichols tuned temperature loop")
    print("=" * 66)
    print(f" Nominal steady state : CA = {x_ss[0]:.4f} mol/L, T = {x_ss[1]:.2f} K, "
          f"Tj = {x_ss[2]:.2f} K")
    print(f" Nominal coolant flow : {QC_NOM:.2f} L/min   (valve range 0-{QC_MAX:.1f})")

    A, _, _ = linearise(x_ss)
    ev = np.linalg.eigvals(A)
    print(f" Open-loop stability  : max Re(eig) = {ev.real.max():+.4f} 1/min "
          f"({'stable' if ev.real.max() < 0 else 'UNSTABLE'})")

    pid = ziegler_nichols_pid(x_ss)
    print("\n Ziegler-Nichols (ultimate gain) tuning")
    print(f"   Ku = {pid['Ku']:.3f} L/min/K,  Pu = {pid['Pu']:.3f} min")
    print(f"   Kp = {pid['Kp']:.3f} L/min/K,  tau_i = {pid['tau_i']:.3f} min,"
          f"  tau_d = {pid['tau_d']:.3f} min")

    t, X_ol, U_ol = simulate(x_ss, controller=None)
    _, X_cl, U_cl = simulate(x_ss, controller=pid)

    m_ol = metrics(t, X_ol[:, 1], "Open-loop")
    m_cl = metrics(t, X_cl[:, 1], "PID closed-loop")
    red = 100.0 * (m_ol["over_K"] - m_cl["over_K"]) / m_ol["over_K"]

    print("\n Performance metrics (+%.0f K coolant-inlet step at t = %.0f min)"
          % (DT_DIST, T_DIST))
    print("-" * 66)
    print(f" {'Metric':<34}{'Open-loop':>14}{'PID (Z-N)':>16}")
    print("-" * 66)
    print(f" {'Peak temperature (K)':<34}{m_ol['peak_K']:>14.3f}{m_cl['peak_K']:>16.3f}")
    print(f" {'Peak overshoot above SP (K)':<34}{m_ol['over_K']:>14.3f}{m_cl['over_K']:>16.3f}")
    print(f" {'Peak overshoot (% of SP, K basis)':<34}{m_ol['over_pct']:>14.3f}{m_cl['over_pct']:>16.3f}")
    print(f" {'Settling time (min, +/-0.1 K)':<34}{m_ol['settle']:>14.2f}{m_cl['settle']:>16.2f}")
    print(f" {'Steady-state error (K)':<34}{m_ol['sse']:>14.3f}{m_cl['sse']:>16.3f}")
    print("-" * 66)
    print(f" Overshoot reduction vs open-loop : {red:.1f} %")
    print(f" Max reactor T below safety limit ({T_SAFE:.0f} K): "
          f"open-loop {'YES' if X_ol[:,1].max() < T_SAFE else 'NO'}, "
          f"PID {'YES' if X_cl[:,1].max() < T_SAFE else 'NO'}")
    print("=" * 66)

    make_plots(t, X_ol, X_cl, U_ol, U_cl)
    print(" Saved: cstr_temperature_response.png, cstr_concentration_control.png")


# =============================================================================
# 7. PLOTTING
# =============================================================================
def make_plots(t, X_ol, X_cl, U_ol, U_cl):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.grid": True, "grid.alpha": 0.3})
    c_ol, c_cl = "#C0392B", "#1F5FA8"

    # --- Figure 1: temperature response -------------------------------------
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.plot(t, X_ol[:, 1], color=c_ol, lw=2.2, label="Open-loop (uncontrolled)")
    ax.plot(t, X_cl[:, 1], color=c_cl, lw=2.2, label="PID closed-loop (Z-N tuned)")
    ax.axhline(T_SP, color="k", ls=":", lw=1.2, label=f"Setpoint ({T_SP:.0f} K)")
    ax.axhline(T_SAFE, color="darkorange", ls="--", lw=1.8,
               label=f"Safety threshold ({T_SAFE:.0f} K)")
    ax.axvline(T_DIST, color="grey", ls="-.", lw=1.0)
    ax.annotate(f"+{DT_DIST:.0f} K coolant-inlet\nstep @ t = {T_DIST:.0f} min",
                xy=(T_DIST, T_SP), xytext=(T_DIST + 3, T_SP - 1.6),
                arrowprops=dict(arrowstyle="->", color="grey"), fontsize=9.5)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Reactor temperature (K)")
    ax.set_title("CSTR Temperature Response to Coolant Inlet Disturbance")
    ax.set_xlim(0, T_END)
    ax.legend(loc="upper right", framealpha=0.95)
    fig.tight_layout()
    fig.savefig("cstr_temperature_response.png", dpi=300)
    plt.close(fig)

    # --- Figure 2: concentration + manipulated variable ---------------------
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    a1.plot(t, X_ol[:, 0], color=c_ol, lw=2.2, label="Open-loop")
    a1.plot(t, X_cl[:, 0], color=c_cl, lw=2.2, label="PID closed-loop")
    a1.axvline(T_DIST, color="grey", ls="-.", lw=1.0)
    a1.set_ylabel("Reactant concentration $C_A$ (mol/L)")
    a1.set_title("Reactant Concentration and Manipulated Coolant Flow")
    a1.legend(loc="best")
    a2.plot(t, U_ol, color=c_ol, lw=2.2, label="Open-loop (fixed flow)")
    a2.plot(t, U_cl, color=c_cl, lw=2.2, label="PID controller output")
    a2.axvline(T_DIST, color="grey", ls="-.", lw=1.0)
    a2.set_xlabel("Time (min)")
    a2.set_ylabel("Coolant flow rate $q_c$ (L/min)")
    a2.set_xlim(0, T_END)
    a2.legend(loc="best")
    fig.tight_layout()
    fig.savefig("cstr_concentration_control.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
