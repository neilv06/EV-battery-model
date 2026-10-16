EV Race Car Battery Cooling Model

A transient thermal model of a liquid-cooled EV race car battery pack, written in Python. It predicts how battery and coolant temperatures change over a 45-minute race and shows which part of the cooling system limits performance.

Contents

	•	Overview
	•	How to run
	•	Assumptions
	•	Theory and equations
	•	Hand calculations
	•	Results
	•	Conclusions
	•	Limitations and future work

Overview

Heat generated in the battery passes through a cold plate into a water–glycol coolant loop, which carries it to a radiator where it is rejected to the air:

Battery pack → cold plate → coolant loop → radiator → ambient air

The model answers three questions:

	1.	How does pack temperature evolve over a 45-minute race?
	2.	Which part of the cooling system limits heat rejection?
	3.	How do coolant flow rate and radiator size affect the peak battery temperature?

How to run

Requires Python 3 with NumPy, SciPy and Matplotlib.

pip install numpy scipy matplotlib
python battery_cooling_model.py

This prints the peak temperatures and energy balance, and saves baseline_race.png and sensitivity.png. Parameters can be changed in the Params class at the top of the script.

Assumptions

	•	Lumped nodes: the battery and coolant each have one uniform temperature.
	•	Constant properties: masses, specific heats and conductances (UA) do not vary.
	•	Joule heating only: heat generation is $I^2R$; reversible (entropic) heat is neglected.
	•	Constant voltage: 800 V, so current is $I = P/V$.
	•	Isothermal battery side in the cold plate.
	•	No direct losses from the pack to the air.
	•	Repeating lap: a 90 s power and speed profile repeated for the whole race.
	•	Illustrative parameters: representative of an 800 V race pack, not a specific car.

Theory and equations

Energy balances

Applying accumulation = in − out + generation to each node gives two coupled ODEs:

$$m_b\, c_{p,b}\,\frac{dT_b}{dt} = Q_{\mathrm{gen}} - Q_{\mathrm{cp}}$$

$$m_c\, c_{p,c}\,\frac{dT_c}{dt} = Q_{\mathrm{cp}} - Q_{\mathrm{rad}}$$

Heat generation

$$I = \frac{P}{V}, \qquad Q_{\mathrm{gen}} = I^2 R(T_b)$$

$$R(T) = R_{\mathrm{ref}} \exp\left[B\left(\frac{1}{T} - \frac{1}{T_{\mathrm{ref}}}\right)\right]$$

Because heating scales with $I^2$, the lap average uses the mean of $P^2$, not the mean power:

$$\bar{Q}_{\mathrm{gen}} = \frac{\overline{P^2}}{V^2}\, R$$

Cold plate (ε-NTU, isothermal hot side)

$$\mathrm{NTU} = \frac{UA_{\mathrm{cp}}}{C_c}, \qquad \varepsilon_{\mathrm{cp}} = 1 - e^{-\mathrm{NTU}}, \qquad Q_{\mathrm{cp}} = \varepsilon_{\mathrm{cp}}\, C_c\, (T_b - T_c)$$

Radiator (ε-NTU, crossflow, both fluids unmixed)

$$v_{\mathrm{face}} = \max\left(v_{\mathrm{fan}},\ \eta_{\mathrm{duct}}\, v_{\mathrm{car}}\right), \qquad C_{\mathrm{air}} = \rho_{\mathrm{air}}\, A_{\mathrm{face}}\, v_{\mathrm{face}}\, c_{p,\mathrm{air}}$$

$$\varepsilon_{\mathrm{rad}} = 1 - \exp\left[\frac{\mathrm{NTU}^{0.22}}{C_r}\left(e^{-C_r\,\mathrm{NTU}^{0.78}} - 1\right)\right], \qquad Q_{\mathrm{rad}} = \varepsilon_{\mathrm{rad}}\, C_{\mathrm{min}}\, (T_c - T_{\mathrm{amb}})$$

Steady state and time constant

$$UA_{\mathrm{overall}} = \left(\frac{1}{\varepsilon_{\mathrm{cp}} C_c} + \frac{1}{\varepsilon_{\mathrm{rad}} C_{\mathrm{min}}}\right)^{-1}, \qquad \tau = \frac{m_b\, c_{p,b}}{UA_{\mathrm{overall}}}$$

The ODEs are solved with SciPy's solve_ivp (LSODA, maximum step 0.5 s).

Hand calculations

Quantity	Value
RMS lap power / mean lap power	205.2 kW / 101.3 kW
Lap-average heat generation at 30 °C → 39 °C	2.68 kW → 2.26 kW
Cold plate conductance, \varepsilon_{\mathrm{cp}} C_c	457 W/K
Radiator conductance, lap average	515 W/K
Overall conductance	242 W/K
Time constant	1157 s (≈ 19 min)
Steady-state battery / coolant temperature	39.3 °C / 34.4 °C

Radiator performance at different speeds:

Condition	C_{\mathrm{min}} (W/K)	NTU	ε	ε·Cmin (W/K)
230 km/h straight	1050	0.857	0.524	551
150 km/h	1050	0.857	0.499	524
90 km/h corner	1050	0.857	0.453	476
Stationary (fan only)	427 (air)	2.108	0.776	331

Full worked calculations are in the project report.

Results

Baseline race

Quantity	Hand calculation	Simulation
Battery temperature at 45 min	38.4 °C	38.8 °C
Coolant temperature at 45 min	≈ 34 °C	34.1 °C
Energy balance error	–	0.004 %

The pack stays well below the 55 °C cell limit. Heat generated (6.47 MJ) equals heat rejected (3.95 MJ) plus heat stored (2.52 MJ).

Design sensitivity

Conclusions

	•	The pack stays within limits: it peaks at 38.8 °C, about 16 K below the cell limit.
	•	A race is too short to reach steady state: with a 19-minute time constant, the pack reaches only about 90% of its steady-state rise.
	•	RMS power drives heating: the 205 kW RMS power produces four times the heat implied by the 101 kW mean, and regen braking adds heat too.
	•	Coolant flow has diminishing returns: above about 0.3 kg/s, extra flow lowers the peak by less than 1 K.
	•	The radiator is the bottleneck, and it loses about 40% of its capacity when the car is stationary.
	•	Cold packs generate more heat: about 19% more at 30 °C than at 39 °C.

Limitations and future work

The model gives an average pack temperature using illustrative parameters, so it is a system-level design tool rather than a prediction of cell temperatures.

Planned improvements:

	1.	Multi-node pack model to capture the hottest cell
	2.	Flow-dependent heat transfer coefficients (Dittus–Boelter)
	3.	State-of-charge-dependent resistance and voltage sag
	4.	Pump and radiator flap control strategy
	5.	Harder scenarios: hot day, qualifying laps, safety car, pit stops
	6.	Real cell parameters and validation against published data
	7.	Parameter uncertainty analysis

Author

Neil Vyawahare, MEng Chemical Engineering, University of Manchester
