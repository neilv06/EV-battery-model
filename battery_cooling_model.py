#Lumped-paramater thermal model of an EV race car battery cooling system

"""
System: battery pack --(cold plate)---> water-glycol loop --
(radiator)--> ambient air

State variables
 T_b : battery pack temperature [degC]
 T_c : coolant loop temperature [degC]

Energy balances
 m_b.cp_b.dT_b/dt = Q_gen - Q_cp
 m_c.cp_c.dT_c/dt = Q_cp - Q_rad

Q_gen : Joule heating, I^2 R(T), with I = P/V
Q_cp : cold plate heat transfer, eps-NTU with an isothermal hot side
Q_rad : radiator heat rejection, eps-NTU crossflow (both fluids unmixed)

Run: python battery_cooling_model.py
"""
from dataclasses import dataclass, replace

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# Paramaters (illustrative values for a ~40 kWh, 800V race pack)

@dataclass
class Params:
  #Battery
  m_b: float = 280.0    # pack mass [kg]
  cp_b: float = 1000.0  # pack specific heat [J/kg/K]
  V: float = 800.0      # pack voltage [V]
  R_ref: float = 0.045  # pack resistance at T_ref [ohm]

  T_ref: float = 25.0   # reference temperature [degC]
  B: float = 1800.0     # Arrhenius-type resistance factor [K]


  #Coolant loop (50/50 water-glycol)
  m_c: float = 6.0      # coolant inventory [kg]
  cp_c: float = 3500.0  # [J/kg/K]
  mdot_c: float = 0.30  # coolant mass flow [Kg/s]
  UA_cp: float = 600.0  # cold plate conductance [W/K]

  #radiator
  UA_rad: float = 900.0 # radiator conductance [W/K]
  A_face: float = 0.12  # radiator face area [m^2]
  duct_eff: float = 0.35 # fraction of car speed seen at radiator face [-]

  v_fan: float = 3.0 # minimum face velcoity from fan [m/s]

  rho_air: float = 1.18 # [kg/m^3]
  cp_air: float = 1005.0 # [J/kg/K]

  #Conditions
  T_amb: float = 30.0  # ambient [degC]
  T0: float = 30.0     # initial temperature of everything [degC]

  t_race: float = 2700.0 # race duration [s] (45 min)

LAP_T = np.array([0, 12, 14, 18, 30, 34, 38, 52, 54, 58, 70, 74, 80, 90]) # s
LAP_P = np.array([300, 300, -250, -250, 120, 120, -200, 300, -250, 80, 80, 300, 300, 300]) *1e3 # W
LAP_V = np.array([140, 230, 220, 120, 110, 150, 100, 210, 200, 90, 100, 130, 160, 140]) / 3.6 # m/s
LAP_LEN = LAP_T[-1]

def lap_conditions(t):
  """Pack power [W] and car speed [m/s] at time t (lap repeats). """
  tl = t % LAP_LEN
  return np.interp(tl, LAP_T, LAP_P), np.interp(tl, LAP_T, LAP_V)

def pack_resistance(T_b, p):
  """Internal resistance rises as the pack gets colder."""
  T_k, Tref_k = T_b + 273.15, p.T_ref + 273.15
  return p.R_ref * np.exp(p.B * (1.0 / T_k - 1.0 / Tref_k))

def heat_generation(t, T_b, p):
  P, _ = lap_conditions(t)
  I = P / p.V
  return I**2 * pack_resistance(T_b, p)

def cold_plate_heat(T_b, T_c, p):
  """Battery side treated as isothermal, so eps = 1 - exp(-NTU)."""
  C_c = p.mdot_c * p.cp_c
  NTU = p.UA_cp / C_c
  eps = 1.0 - np.exp(-NTU)
  return eps * C_c * (T_b - T_c)

def radiator_heat(t, T_c, p):
  """Crossflow, both fluids unmixed (standard eps-NTU correlation)."""
  _, v_car = lap_conditions(t)
  v_face = max(p.v_fan, p.duct_eff * v_car)
  C_air = p.rho_air * p.A_face * v_face * p.cp_air
  C_c = p.mdot_c * p.cp_c
  C_min, C_max = min(C_air, C_c), max(C_air, C_c)
  Cr = C_min / C_max
  NTU = p.UA_rad / C_min
  eps = 1.0 - np.exp((NTU**0.22 / Cr) * (np.exp(-Cr * NTU**0.78) -1.0))
  return eps * C_min * (T_c - p.T_amb)

def rhs(t, y, p):
  T_b, T_c = y
  Q_gen = heat_generation(t, T_b, p)
  Q_cp = cold_plate_heat(T_b, T_c, p)
  Q_rad = radiator_heat(t, T_c, p)
  dTb = (Q_gen - Q_cp) / (p.m_b * p.cp_b)
  dTc = (Q_cp - Q_rad) / (p.m_c * p.cp_c)
  return [dTb, dTc]

def simulate(p):
  t_eval = np.arange(0, p.t_race + 0.5, 0.5) # uniform 0.5 s output grid
  sol = solve_ivp(rhs, (0, p.t_race), [p.T0, p.T0], args=(p,), method="LSODA", max_step=0.5, t_eval=t_eval)
  t = sol.t
  T_b, T_c = sol.y
  Q_gen = np.array([heat_generation(ti, Tb, p) for ti, Tb in zip(t, T_b)])
  Q_cp = cold_plate_heat(T_b, T_c, p)
  Q_rad = np.array([radiator_heat(ti, Tc, p) for ti, Tc in zip(t, T_c)])
  return dict(t=t, T_b=T_b, T_c=T_c, Q_gen=Q_gen, Q_cp=Q_cp, Q_rad=Q_rad)

def energy_balance_check(r, p):
  """Heat generated should equal heat rejected + heat stored"""
  E_gen = np.trapezoid(r["Q_gen"], r["t"])
  E_rad = np.trapezoid(r["Q_rad"], r["t"])
  E_stored = (p.m_b * p.cp_b * (r["T_b"][-1] - p.T0) + p.m_c * p.cp_c * (r["T_c"][-1] - p.T0))
  err = (E_gen - E_rad - E_stored) / E_gen *100
  print(f"Heat generated : {E_gen/1e6:7.2f} MJ")
  print(f"Heat rejected : {E_rad/1e6:7.2f} MJ")
  print(f"Heat stored : {E_stored/1e6:7.2f} MJ")
  print(f"Balance error : {err:7.3f} %")

#plots

def plot_baseline(r, p, fname):
  t_min = r["t"] / 60
  fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
  ax1.plot(t_min, r["T_b"], label="Battery pack", lw=2)
  ax1.plot(t_min, r["T_c"], label="Coolant", lw=1.5)
  ax1.axhline(p.T_amb, color="grey", ls=":", label="Ambient")
  ax1.axhline(55, color="red", ls="--", label="Typical cell limit (55 °C)")
  ax1.set_ylabel("Temperature [°C]")
  ax1.set_title("Battery cooling system: 45 min race")
  ax1.legend(loc="upper left")
  ax1.grid(alpha=0.3)

  # Trailing 1-lap average so the trend is readable through the lap spikes
  n = int(LAP_LEN / (r["t"][1] - r["t"][0]))
  def lap_avg(x):
    return np.convolve(x, np.ones(n) / n, mode="valid")
  t_avg = t_min[n - 1:]
  ax2.plot(t_avg, lap_avg(r["Q_gen"]) / 1e3, label="Generated (lap avg)")
  ax2.plot(t_avg, lap_avg(r["Q_cp"]) / 1e3, label="Into coolant (lap avg)")
  ax2.plot(t_avg, lap_avg(r["Q_rad"]) / 1e3, label="Rejected by radiator (lap avg)")
  ax2.set_xlabel("Time [min]")
  ax2.set_ylabel("Heat flow [kW]")
  ax2.legend(loc="lower right")
  ax2.grid(alpha=0.3)

  fig.tight_layout()
  fig.savefig(fname, dpi=150)
  plt.close(fig)

def plot_sensitivity(p, fname):
  flows = np.linspace(0.05, 0.8, 12)
  rad_sizes = [500, 900, 1500]
  fig, ax = plt.subplots(figsize=(8, 5))
  for UA in rad_sizes:
    peak = [simulate(replace(p, mdot_c=m, UA_rad=UA))["T_b"].max() for m in flows]
    ax.plot(flows, peak, "o-", label=f"Radiator UA = {UA} W/K")
  ax.axhline(55, color='red', ls='--', label="Cell limit")
  ax.set_xlabel("Coolant mass flow [kg/s]")
  ax.set_ylabel("Peak battery temperature [°C]")
  ax.set_title("Design sensitivity: pump flow vs radiator size")
  ax.legend()
  ax.grid(alpha=0.3)
  fig.tight_layout()
  fig.savefig(fname, dpi=150)
  plt.close(fig)

if __name__ == "__main__":
  p = Params()
  r = simulate(p)
  print(f"Peak battery temperature: {r['T_b'].max(): .1f} °C")
  print(f"Final coolant temperature: {r['T_c'][-1]: .1f} °C\n")
  energy_balance_check(r, p)
  plot_baseline(r, p, "baseline_race.png")
  plot_sensitivity(p, "sensitivity.png")
  print("\nSaved baseline_race.png and sensitivity.png")
