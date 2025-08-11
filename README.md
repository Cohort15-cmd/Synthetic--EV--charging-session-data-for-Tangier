# Synthetic--EV--charging-session-data-for-Tangier
Synthetic electric vehicle (EV) charging session data for Tangier,

Overview
This project generates synthetic electric vehicle (EV) charging session data for Tangier, Morocco, using a physics-aware, behavior-driven simulation tied to real stations and a fixed fleet of EVs in Tangier.

It produces:
- tangier_ev_charging_synthetic.csv — per-session records (station, EV model, timestamps, SOC in/out, kWh, ambient temperature, queue wait, etc.)
- station_daily_report.csv — daily per-station KPIs (sessions, energy, utilization, peak concurrency, queue stats)

---

## What the Code Does
This README reflects the current MATLAB script in this repo. Where helpful, we mention inspiration; however, the formulas below are exactly what the code implements.

1) Daily driving -> need to charge
Each vehicle draws a daily distance from a mixture of log-normals to create diversity of trip lengths:
- Probabilities 0.45 / 0.40 / 0.15 for short (logN(mu=log 12, sigma=0.5)), medium (logN(log 35, 0.5)), long (logN(log 90, 0.6), capped at 250 km).

Temperature-dependent energy per km:
- e_km(T) = e0 * (1 + a*max(0,20 - T) + b*max(0, T - 30)), with e0=0.16, a=0.008, b=0.004.
Trip energy: E_trip = D * e_km(T). SOC update: dSOC = 100 * E_trip / B (battery B in kWh).

The probability of attempting to charge today is a range-based logistic (different thresholds for AC vs DC urgency):
- p_AC = 1 / (1 + exp((range_km - Rth_AC)/k_AC)), with Rth_AC=40 km, k_AC=20.
- p_DC = 1 / (1 + exp((range_km - Rth_DC)/k_DC)), with Rth_DC=120 km, k_DC=25.
- p_need = max(p_AC, 0.25 * p_DC).

Why: This yields emergent, city-level arrivals without forcing a Poisson process. It couples the charging need to the actual remaining range after daily driving.

2) Arrival time within the day
Given the vehicle decides to charge, we draw one timestamp within the day using a two-peak intraday profile (morning & evening), sharpened by a shape knob lambda0_time_shape:
- w(h) proportional to [0.3*N(h; 9, 1.2) + 0.5*N(h; 19.5, 1.8) + 0.2/24] ^ lambda0.
We normalize w to a PMF, draw an hour by inverse CDF, and add a uniform minute in [0, 60).

3) Station choice & compatibility
We score each real Tangier station by: connector compatibility, site power, operator (Tesla drivers prefer Supercharger), and a DC-capability bonus. We then draw a station with probability proportional to the score.

4) Charging physics (power, taper, efficiency)
For each 1% SOC step from arrival to target, we compute instantaneous power and integrate time/energy.

Power available at SOC:
- P(SOC, T, isDC) = min(P_site, P_vehicle) * derate_T(T) * taper(SOC, isDC).
  - Thermal derate: derate_T(T) = max(0.85, 1 - 0.01*max(0, T - 35)).
  - Taper (no taper <= 80%; linear after 80%): taper = 1 for SOC<=80; else taper = 1 - beta*(SOC - 80)/20, with beta = 0.7 (DC) or 0.5 (AC), floored at 0.2.
  - Efficiency vs temperature: eta(T) = min(0.98, max(0.85, 0.94 - 0.0008*(25 - T)^2)).

Time/energy integration per 1% SOC step:
- dE = B / (100 * eta(T)); dt_min = 60 * dE / P(SOC, T, ...).
- Duration_min = sum dt_min from ceil(SOC_arr) to floor(SOC_target); Energy_kWh = sum dE.

5) Target SOC selection
We draw a plausible departure SOC by charger type and enforce a minimum delta:
- DC: targets {70, 80, 85, 90}% with probs {0.1, 0.5, 0.3, 0.1}; min delta 10%.
- AC: targets {80, 90, 100}% with probs {0.25, 0.5, 0.25}; min delta 5%.
We clip to [arrival, 100].

6) Queueing & stall assignment
Per station, we track next-free time for each stall. If all are busy at arrival, we compute the wait to the soonest-free stall and accept only if:
- wait <= base * (0.8 + 0.002 * min(Duration_min, 180)), with base = 20 min (DC) or 15 min (AC).
If accepted, plug-in = arrival + wait, and stall’s next-free = plug-out.

7) Temperature model (Tangier)
For each day: draw Tmin/Tmax from monthly normals + AR(1) anomaly (summer heatwave chance), then hourly temperature via a coastal diurnal cosine (min ~05:00, max ~15:00):
- T(h) = Tmin + 0.5*(Tmax - Tmin) * [1 - cos(2*pi*(h - 5)/24)].
Ambient temperature used in physics is the value at the arrival time.

---

Why we chose these formulations (and limits)
- Emergent arrivals from driving need tie charging to range/SOC, rather than assuming a Poisson counter.
- 1% SOC integration guarantees physical plausibility (no kW or kWh exceeding caps; taper handled).
- Simple queue tolerance is heuristic; can be calibrated with surveys if available.
- Climate model follows monthly normals; AR(1) adds day-to-day persistence. (Diurnal phase fixed; can be refined.)

---

## Data & Assumptions
- Stations: Real Tangier sites (name, lat/lon, connectors, power, stalls). Sources: PlugShare / Chargemap / Tesla maps (cross-checked manually).
- Fleet: Fixed 300 EVs with brand/model mix you provided; each has battery capacity and AC/DC caps from public specs.
- Parameters (Rth, taper beta, efficiency curve, derate) are engineering assumptions reasonable for Moroccan coastal climate; tune if you have calibration data.

---

## How to Run
```
clear; clc
GenerateTangierEVData_Logic
```
Outputs:
- tangier_ev_charging_synthetic.csv
- station_daily_report.csv

---

