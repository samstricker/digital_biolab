# Mammalian Cell Culture Simulation
## User Guide

To use the mammalian cell culture simulator, you can follow these steps:

1. **Define Reactor Type**: Create a JSON file specifying the reactor parameters (e.g., scale factor, temperature and pH parameters, noise parameters).
2. **Define Cell Type**: Create a JSON file specifying the cell parameters (e.g., growth rates, yield coefficients).
3. **Define Recipe**: Create a JSON file specifying the initial conditions, feeding strategy, and runtime.
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

### Growth Inhibition

Growth inhibition reduces the growth rate:

$$\mu_{inhibited} = \mu \cdot g_{inhibition}$$

This represents the reduction in growth rate due to reactor-specific conditions (e.g., bioreactor scale effects, oxygen limitations, shear stress).

### Death Rate
The specific death rate is:

$$k_{death} = k_{death,max} \cdot \frac{k_{death,\mu}}{\mu + k_{death,\mu}} \cdot s_{scale} \cdot \left(1 + \frac{c_{lac}}{k_{death,lac}} + \frac{c_{ammonia}}{k_{death,ammonia}} + k_{death,stress} \cdot \left(1 - f_{temp} \cdot f_{pH}\right)\right)$$

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
There are two types of noise that can be added to the simulations to create more realistic and variable outcomes, mimicking the inherent stochasticity in biological systems and measurement processes.

**Process Noise (Brownian Motion):**
Process noise represents intrinsic stochasticity in biological systems. It is modeled as cumulative Brownian motion increments, scaled relative to the current state value:

$$dX_t = X_t \cdot \sigma_{process} \sqrt{dt} \cdot dW_t$$

where $X_t$ is the current state value, $\sigma_{process}$ is the relative noise intensity (dimensionless), $dW_t$ is a Wiener process increment, and $dt$ is the time step. The noise scales with the variable magnitude, creating realistic proportional variability.

**Measurement Noise (Gaussian White Noise):**
Measurement noise represents sensor/analytical inaccuracies. It is modeled as independent Gaussian random values at each timepoint, scaled relative to the measured value:

$$X_{measured} = X_{true} + X_{true} \cdot \sigma_{measurement} \cdot \mathcal{N}(0, 1)$$

where $X_{true}$ is the true value and $\sigma_{measurement}$ is the relative measurement noise standard deviation (dimensionless). This noise is proportional to the measurement magnitude and uncorrelated across timepoints and species.


