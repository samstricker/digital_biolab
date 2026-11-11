# Mammalian Cell Culture Simulation
## User Guide

To use the mammalian cell culture simulator, you can follow these steps:

1. **Define Reactor Type**: Create a JSON file specifying the reactor parameters (e.g., scale factor, temperature and pH parameters, noise parameters).
2. **Define Cell Type**: Create a JSON file specifying the cell parameters (e.g., growth rates, yield coefficients).
3. **Define Recipe**: Create a JSON file specifying the initial conditions, feeding strategy, and runtime.
   The recipe also sets the impeller `stirring_rate` [rpm]; like temperature and pH, a feeding event can switch it.
   Each feeding event adds its `feed_concentrations` (amounts per culture volume before the feed) and then dilutes the whole culture by `1 + volume_fraction`, the feed volume relative to the culture volume. It can optionally set a new `temperature` and/or `ph_value`; the culture runs at that setpoint from the feed on (e.g. a temperature or pH shift). Leave them out (or `null`) to keep the current setpoint.
4. **Run Simulation**: Use the `MammalianCellCultureSimulator` class to run the simulation with the defined reactor, cell type, and recipe.
5. **Analyze Results**: The simulator will return time-series data for all state variables, which can be analyzed to evaluate process performance, optimize conditions, or test control strategies.

You can refer to the `setup/experiment_setup.ipynb` notebook for examples on how to create these JSON files and run the simulation and analyze results.

## Reference
This simulator is based on the publications:

*Stephen Craven, Nishikant Shirsat, Jessica Whelan, and Brian Glennon. Process model comparison and transfer-
ability across bioreactor scales and modes of operation for a mammalian cell bioprocess. Biotechnology Progress,
29(1):186–196, 2013.*


*A. Martens, M. Neufang, A. Butté, M. von Stosch, A. del R. Chanona, and L. M. Helleckes, “Holistic Bioprocess Development Across Scales Using Multi-Fidelity Batch Bayesian Optimization,” Aug. 14, 2025, arXiv.*

and was adapted and extended in 2026 by Sam Stricker at Imperial College London.

## Theoretical Model Summary
Below you can see a summary of the used differential equations to simulate a CHO bioprocess.

### State Variables
- $c_{product}$: Product concentration
- $x_{total}$: Total cell density
- $x_{viable}$: Viable cell density
- $x_{dead}$: Dead cell density
- $c_{glc}$: Glucose concentration
- $c_{gln}$: Glutamine concentration
- $c_{lac}$: Lactate concentration
- $c_{ammonia}$: Ammonia concentration
- $c_{aggregates}$: Product aggregate concentration
- $c_{glycated}$: Glycated product concentration
- $c_{HCP}$: Host cell protein concentration
- $c_{O_2}$: Dissolved oxygen concentration
- $c_{CO_2}$: Dissolved CO2 concentration

### Growth Rate
The specific growth rate is calculated as:

$$\mu = \mu_{max} \cdot \frac{c_{glc}}{k_{sat,glc} + c_{glc}} \cdot \frac{c_{gln}}{k_{sat,gln} + c_{gln}} \cdot \frac{k_{sat,lac}}{k_{sat,lac} + c_{lac}} \cdot \frac{k_{sat,ammonia}}{k_{sat,ammonia} + c_{ammonia}} \cdot f_{temp} \cdot f_{pH}$$

where $f_{temp}$ and $f_{pH}$ are temperature and pH penalty factors.

### Environmental Effects

#### Temperature Effect
The temperature penalty factor uses asymmetric Gaussian functions:

$$f_{temp} = r_{temp} \cdot \begin{cases}
\exp\left(-steepness_{temp,l} \cdot \left(\frac{T - T_{opt}}{smoothness_{temp,l}}\right)^2\right) & \text{if } T < T_{opt} \\
\exp\left(-steepness_{temp,r} \cdot \left(\frac{T - T_{opt}}{smoothness_{temp,r}}\right)^2\right) & \text{if } T \geq T_{opt}
\end{cases}$$

where $T_{opt}$ is the optimal temperature (e_activation) and $r_{temp}$ is the temperature robustness factor. The factor is clipped to a minimum of 0.2.

#### pH Effect
The pH penalty factor similarly uses asymmetric Gaussian functions:

$$f_{pH} = r_{pH} \cdot \begin{cases}
\exp\left(-steepness_{pH,l} \cdot \left(\frac{pH - pH_{opt}}{smoothness_{pH,l}}\right)^2\right) & \text{if } pH < pH_{opt} \\
\exp\left(-steepness_{pH,r} \cdot \left(\frac{pH - pH_{opt}}{smoothness_{pH,r}}\right)^2\right) & \text{if } pH \geq pH_{opt}
\end{cases}$$

where $pH_{opt}$ is the optimal pH and $r_{pH}$ is the pH robustness factor. The factor is clipped to a minimum of 0.2.

Note: product formation uses its own temperature and pH responses (see below).

### Oxygen Transfer and Shear

The reactor is a stirred, air-sparged tank. The impeller power input (turbulent regime) is

$$P = N_p \cdot \rho \cdot N^3 \cdot D^5$$

with power number $N_p$, medium density $\rho$, stirring rate $N$ and impeller diameter $D$. The oxygen transfer coefficient follows the van't Riet (1979) correlation for non-coalescing media (SI units, $k_La$ in s$^{-1}$):

$$k_La = C \cdot \left(\frac{P}{V}\right)^{\alpha} \cdot v_s^{\beta}, \qquad C = 0.002,\ \alpha = 0.7,\ \beta = 0.2$$

and the dissolved oxygen balance is

$$\frac{dc_{O_2}}{dt} = k_La \cdot \left(c^*_{O_2} - c_{O_2}\right) - q_{O_2} \cdot x_{viable} \cdot f_{met} \cdot f_{O_2}, \qquad f_{O_2} = \frac{c_{O_2}}{k_{sat,O_2} + c_{O_2}}$$

$f_{O_2}$ also multiplies the growth rate and the production rate, and a lack of oxygen raises the lactate yield by the factor $1 + k_{lac,hyp} \cdot (1 - f_{O_2})$. Shear damage depends on the local maximum energy dissipation near the impeller, $\varepsilon_{max} = \phi \cdot P / (\rho V)$ (Nienow, 2006), and adds $(\varepsilon_{max} / \varepsilon_{crit})^2$ to the death-rate stress. Too little stirring starves a dense culture of oxygen; too much damages cells.

**Dissolved oxygen control.** With a `do_setpoint` (fraction of air saturation) in the recipe, an ideal controller enriches the sparged gas with pure oxygen. It picks the gas oxygen fraction $y_{O_2} \in [0.21, 1]$ at which transfer matches the uptake with dissolved oxygen at the setpoint, $y_{O_2} = 0.21 \cdot (c_{sp} + OUR / k_La) / c^*_{O_2}$, and the saturation concentration becomes $c^*_{O_2} \cdot y_{O_2} / 0.21$. The fraction is reported as `o2_fraction`. Without a setpoint (`None`) the tank is sparged with air.

**Dissolved CO2.** Cells produce CO2 with the respiratory quotient, $CER = RQ \cdot OUR$, and the gas strips it:

$$\frac{dc_{CO_2}}{dt} = CER - k_{strip} \cdot \left(c_{CO_2} - c_{CO_2,gas}\right), \qquad k_{strip} = \left(\frac{1}{0.89 \, k_La} + \frac{1}{(Q_g/V) \cdot H}\right)^{-1}$$

The liquid-side coefficient is $0.89 \, k_La$ (square root of the ratio of the CO2 and O2 diffusivities), in series with the capacity of the gas: once the bubbles saturate, they carry off at most the gas flow per volume $Q_g/V$ times the gas/liquid partition $H = 1.68$ at 37 C. Large tanks run at low gas flow per volume, so CO2 builds up there. The sparged gas may carry CO2 for pH control (`gas_co2_fraction`); $c_{CO_2,gas}$ is the dissolved CO2 in equilibrium with it (713 mmHg dry gas at 37 C). With $pCO_2 = c_{CO_2} / 0.0307$ mmHg, growth and production are divided by $1 + (pCO_2 / pCO_{2,crit})^3$, so normal levels (30-60 mmHg) are nearly harmless and effects start around 100 mmHg [`respiratory_quotient`, `pco2_crit`].

**Mixing.** The 95 % mixing time follows $\theta = 5.9 \, T^{2/3} \bar\varepsilon^{-1/3} (D/T)^{-1/3}$ (Grenville; Nienow, 1997). With slow mixing, cells pass through zones of concentrated feed and base and of low oxygen: with $g = (\theta / \theta_{crit})^2$, the growth rate is divided by $1 + g$ and $g$ is added to the death-rate stress. Stirring therefore trades off mixing against shear.

*K. van't Riet, "Review of measuring methods and results in nonviscous gas-liquid mass transfer in stirred vessels," Ind. Eng. Chem. Process Des. Dev. 18(3):357–364, 1979.*
*A. W. Nienow, "On impeller circulation and mixing effectiveness in the turbulent flow regime," Chem. Eng. Sci. 52(15):2557–2565, 1997.*
*A. W. Nienow, "Reactor engineering in large scale animal cell culture," Cytotechnology 50:9–33, 2006.*

### Growth Inhibition

Growth inhibition reduces the growth rate:

$$\mu_{inhibited} = \mu \cdot g_{inhibition}$$

This represents the reduction in growth rate due to reactor-specific conditions (e.g., bioreactor scale effects, oxygen limitations, shear stress).

### Death Rate
The specific death rate is:

$$k_{death} = k_{death,max} \cdot \frac{k_{death,\mu}}{\mu + k_{death,\mu}} \cdot s_{scale} \cdot \left(1 + \frac{c_{lac}}{k_{death,lac}} + \frac{c_{ammonia}}{k_{death,ammonia}} + k_{death,stress} \cdot \left(1 - f_{temp} \cdot f_{pH}\right) + \left(\frac{\varepsilon_{max}}{\varepsilon_{crit}}\right)^2\right)$$

Lactate (acidification, osmolality) and ammonia toxicity raise the death rate, and so do temperature and pH far from the growth optima (stress), so a culture cannot be parked at harsh conditions without losing cells.


### Cell Population Dynamics

$$\frac{dx_{total}}{dt} = \mu_{inhibited} \cdot x_{viable} - k_{lysis} \cdot x_{dead}$$

$$\frac{dx_{viable}}{dt} = (\mu_{inhibited} - k_{death}) \cdot x_{viable}$$

$$\frac{dx_{dead}}{dt} = k_{death} \cdot x_{viable} - k_{lysis} \cdot x_{dead}$$

### Metabolite Dynamics

#### Product Formation and Aggregation
$$q_{P} = x_{viable} \cdot \left(y_{product,ng} + y_{product,viable} \cdot \mu\right) \cdot f_{temp,P} \cdot f_{pH,P} \cdot \frac{c_{glc}}{k_{sat,glc} + c_{glc}} \cdot \frac{c_{gln}}{k_{sat,gln} + c_{gln}} \cdot \frac{k_{inh,lac,P}}{k_{inh,lac,P} + c_{lac}} \cdot \frac{k_{inh,amm,P}}{k_{inh,amm,P} + c_{ammonia}}$$

$$\frac{dc_{product}}{dt} = q_{P} - 2 \cdot k_{agg} \cdot c_{product}^2 - k_{glyc} \cdot c_{glc} \cdot c_{product}$$

Only viable cells produce. $f_{temp,P}$ uses the same functional form as $f_{temp}$ with its own production parameters (`t_opt_product`, `temperature_robustness_factor_product`, `steepness_factor_t_{l,r}_product`, `smoothness_factor_t_{l,r}_product`); a production optimum below the growth optimum models the higher specific productivity at mild hypothermia. $f_{pH,P}$ uses the same functional form as $f_{pH}$ but with its own production parameters (`ph_opt_product`, `ph_robustness_factor_product`, `steepness_factor_ph_{l,r}_product`, `smoothness_factor_ph_{l,r}_product`), so the production pH optimum can differ from the growth optimum, and its own lower bound (`ph_min_factor_product`, instead of 0.2) so that little product is made far from the production optimum. Production stops when glucose or glutamine is depleted.

#### Product Aggregation
$$\frac{dc_{aggregates}}{dt} = k_{agg} \cdot c_{product}^2$$

#### Product Glycation
$$\frac{dc_{glycated}}{dt} = k_{glyc} \cdot c_{glc} \cdot c_{product}$$

Glucose reacts non-enzymatically with the product, so glucose left in the culture degrades product quality and overfeeding glucose has a cost.

#### Glucose Consumption
$$\frac{dc_{glc}}{dt} = x_{viable} \cdot \left(-\frac{\mu}{y_{cells,glc}} - m_{glc}\right) - y_{glc,product} \cdot q_{P}$$

#### Glutamine Dynamics
$$q_{gln} = x_{viable} \cdot \left(\frac{\mu}{y_{cells,gln}} + m_{gln}\right) + y_{gln,product} \cdot q_{P}$$

$$\frac{dc_{gln}}{dt} = -q_{gln} - k_{deg,gln} \cdot c_{gln}$$

#### Lactate Production
$$\frac{dc_{lac}}{dt} = -y_{lac}(c_{glc}) \cdot \frac{dc_{glc}}{dt} - q_{lac,uptake} \cdot x_{viable} \cdot \frac{k_{switch}}{k_{switch} + c_{glc}} \cdot \frac{c_{lac}}{k_{sat,lac,up} + c_{lac}}, \qquad y_{lac}(c_{glc}) = y_{lactate} + \left(y_{lactate,max} - y_{lactate}\right) \frac{c_{glc}}{k_{overflow} + c_{glc}}$$

Overflow metabolism: the lactate yield on all glucose taken up (growth, maintenance and product) rises as glucose becomes abundant. Lactate inhibits growth through $k_{sat,lac}$ (realistic values 20-60 mM), production through $k_{inh,lac,P}$ and raises the death rate through $k_{death,lac}$. When glucose runs low (below about $k_{switch}$), cells switch to consuming lactate. A higher pH raises lactate formation: $y_{lac}$ is multiplied by $1 + k_{lac,pH} \cdot (pH - pH_{ref,lac})$ [`k_lactate_ph`, `ph_ref_lactate`]. A lower temperature slows maintenance metabolism and lactate formation: $m_{glc}$, $m_{gln}$ and $y_{lac}$ are multiplied by $Q_{10}^{(T - T_{ref})/10}$ [`q10_metabolism`, `t_ref_metabolism`]. The lactate yield is capped at 2 lactate per glucose (glycolysis stoichiometry).

#### Ammonia Production
$$\frac{dc_{ammonia}}{dt} = y_{ammonia} \cdot q_{gln} + k_{deg,gln} \cdot c_{gln}$$

#### Product Purity

Product purity is calculated as the share of active product among all protein in the culture:

$$\text{purity} = \frac{c_{product}}{c_{product} + c_{aggregates} + c_{glycated} + c_{HCP}}$$

where host cell protein accumulates as dead cells lyse:
$$\frac{dc_{HCP}}{dt} = hcp_{per\_cell} \cdot k_{lysis} \cdot x_{dead}$$

### Parameters
- $c_{product}$: Active product concentration
- $c_{aggregates}$: Product aggregate concentration
- $c_{HCP}$: Host cell protein released from cell lysis
- $c_{glc}$: Glucose concentration
- $c_{gln}$: Glutamine concentration
- $c_{lac}$: Lactate concentration
- $c_{ammonia}$: Ammonia concentration
- $\mu_{max}$: Maximum specific growth rate [`mu_max`]
- $k_{sat,glc}$, $k_{sat,gln}$, $k_{sat,lac}$, $k_{sat,ammonia}$: Monod saturation constants [`k_sat_glc`, `k_sat_gln`, `k_sat_lac`, `k_sat_ammonia`]
- $y_{cells,glc}$, $y_{cells,gln}$: Cell yield coefficients [`y_cells_glc`, `y_cells_gln`]
- $y_{product,ng}$, $y_{product,viable}$: Non-growth- and growth-associated production per viable cell [`y_product_non_growth`, `y_product_from_viable`]
- $y_{lactate}$, $y_{ammonia}$: Byproduct yield coefficients [`y_lactate`, `y_ammonia`]
- $m_{glc}$, $m_{gln}$: Maintenance coefficients [`m_glc`, `m_gln`]
- $k_{death,max}$, $k_{death,\mu}$: Death rate parameters [`k_death_max`, `k_death_mu`]
- $k_{lysis}$: Cell lysis rate [`k_lysis`]
- $k_{deg,gln}$: Glutamine degradation rate [`k_deg_gln`]
- $k_{agg}$: Product aggregation rate constant [`k_agg_product`]
- $k_{glyc}$: Product glycation rate constant by glucose [`k_glycation`]
- $hcp_{per\_cell}$: HCP released per lysed cell [`hcp_per_cell`]
- $q_{O_2}$, $k_{sat,O_2}$, $k_{lac,hyp}$, $\varepsilon_{crit}$: Specific oxygen uptake, oxygen saturation constant, rise of the lactate yield without oxygen, and critical local energy dissipation [`q_o2`, `k_sat_o2`, `k_lactate_hypoxia`, `eps_crit`]
- $\theta_{crit}$: Mixing time at which gradients double the death rate [`t_mix_crit`]
- $T$, $D$, $N_p$, $v_s$, $C$, $\alpha$, $\beta$, $\phi$, $c^*_{O_2}$: Reactor geometry and oxygen transfer [`tank_diameter`, `impeller_diameter`, `power_number`, `superficial_gas_velocity`, `kla_constant`, `kla_power_exponent`, `kla_gas_exponent`, `max_to_mean_dissipation`, `o2_saturation`]
- $s_{scale}$: Reactor scale factor [`scale_factor`]
- $g_{inhibition}$: Growth inhibition exponent [`growth_inhibition`]
- $T_{opt}$: Optimal temperature for growth [`e_activation`]
- $r_{temp}$: Temperature robustness factor [`temperature_robustness_factor`]
- $steepness_{temp,l}$, $steepness_{temp,r}$: Steepness factors for temperature effect on the left and right side of the optimal temperature [`steepness_factor_t_l`, `steepness_factor_t_r`]
- $pH_{opt}$: Optimal pH for growth [`ph_opt`]
- $r_{pH}$: pH robustness factor [`ph_robustness_factor`]
- $steepness_{pH,l}$, $steepness_{pH,r}$: Steepness factors for pH effect on the left and right side of the optimal pH [`steepness_factor_ph_l`, `steepness_factor_ph_r`]
- $y_{glc,product}$: Glucose consumed per unit product formed [`y_glc_product`]
- $y_{lactate}$, $y_{lactate,max}$, $k_{overflow}$: Lactate yield at low and high glucose, and the glucose at which it is halfway [`y_lactate`, `y_lactate_max`, `k_overflow_glc`]
- $k_{inh,lac,P}$: Lactate inhibition constant for production [`k_inh_lac_product`]
- $k_{death,lac}$: Lactate concentration that doubles the death rate [`k_death_lac`]
- $y_{gln,product}$: Glutamine consumed per unit product formed [`y_gln_product`]
- $k_{inh,amm,P}$: Ammonia inhibition constant for production [`k_inh_ammonia_product`]
- $k_{death,ammonia}$: Ammonia concentration that doubles the death rate [`k_death_ammonia`]
- $k_{death,stress}$: Extra death far from the growth temperature and pH optima [`k_death_stress`]
- $T_{opt,P}$, $r_{temp,P}$, steepness/smoothness: temperature parameters for product formation [`t_opt_product`, `temperature_robustness_factor_product`, `steepness_factor_t_l_product`, `steepness_factor_t_r_product`, `smoothness_factor_t_l_product`, `smoothness_factor_t_r_product`]
- $q_{lac,uptake}$, $k_{switch}$, $k_{sat,lac,up}$: Maximum lactate consumption, the glucose below which it switches on, and its saturation constant [`q_lac_uptake`, `k_glc_lac_switch`, `k_sat_lac_uptake`]
- $pH_{opt,P}$, $r_{pH,P}$, steepness/smoothness: pH parameters for product formation [`ph_opt_product`, `ph_robustness_factor_product`, `steepness_factor_ph_l_product`, `steepness_factor_ph_r_product`, `smoothness_factor_ph_l_product`, `smoothness_factor_ph_r_product`]

### Noise in Simulations
Noise enters where it arises physically. Everything defaults to zero: `NoiseParameters()` is the noiseless, deterministic model, and `NoiseParameters.typical()` gives levels typical of CHO fed-batch processes. Every level can be set individually (fractions, 0.05 = 5 %).

| source | parameters | model |
|---|---|---|
| Batch variability | `batch.mu_max`, `productivity`, `q_o2`, `k_death_max`, `y_cells_glc`, `seeding_density`, `initial_glucose`, `initial_glutamine` | once per run, the parameter or starting value is multiplied by a lognormal factor with mean 1 and the given coefficient of variation |
| Process noise | `process.growth`, `production`, `correlation_time`, `update_interval` | growth and production rates are multiplied by $e^{\xi(t)}$, $\xi$ an Ornstein-Uhlenbeck process, piecewise constant over the update interval; acts inside the dynamics, so all states stay consistent |
| Actuator noise | `actuator.temperature` [C], `ph` [pH], `stirring`, `feed` | per run, setpoints are missed by an offset (temperature, pH) or a relative error (stirring); every feed has a relative dosing error |
| Measurement noise | `measurement.relative[var]`, `absolute[var]`, `calibration[var]`, `censor_below_detection` | the true value times a per-run calibration factor and a per-sample lognormal, mean-preserving error, plus a Gaussian absolute error (three times it is the detection limit); never negative. The cell counter measures `x_total` and `viability` (viable = total x viability, dead = total - viable); for `purity` the errors act on the impurity fraction. Optionally, values below the detection limit are reported as NaN |

`MammalianCellcultureReactor.simulate(recipe, seed)` returns a `SimulationResult` with the `true` trajectory (batch, process and actuator noise) and the `measured` one (measurement noise on top); `run_experiment` returns the measured one. Each source draws from its own random stream, so switching one off does not change the others, and every recipe in `run_dataset` gets independent noise. Reasonable ranges for every setting and an analysis of each source are in `research/noise_analysis.ipynb`.

