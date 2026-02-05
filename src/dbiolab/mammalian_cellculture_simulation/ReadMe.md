# Reference
This simulator is based on the publications:

*Stephen Craven, Nishikant Shirsat, Jessica Whelan, and Brian Glennon. Process model comparison and transfer-
ability across bioreactor scales and modes of operation for a mammalian cell bioprocess. Biotechnology Progress,
29(1):186–196, 2013.*


*A. Martens, M. Neufang, A. Butté, M. von Stosch, A. del R. Chanona, and L. M. Helleckes, “Holistic Bioprocess Development Across Scales Using Multi-Fidelity Batch Bayesian Optimization,” Aug. 14, 2025, arXiv.*


# Differential Equations
Below you can see a summary of the used differential equations to simulate a CHO bioprocess.

## State Variables
- $c_{product}$: Product concentration
- $x_{total}$: Total cell density
- $x_{viable}$: Viable cell density
- $x_{dead}$: Dead cell density
- $c_{glc}$: Glucose concentration
- $c_{gln}$: Glutamine concentration
- $c_{lac}$: Lactate concentration
- $c_{ammonia}$: Ammonia concentration

## Growth Rate
The specific growth rate is calculated as:

$$\mu = \mu_{max} \cdot \frac{c_{glc}}{k_{sat,glc} + c_{glc}} \cdot \frac{c_{gln}}{k_{sat,gln} + c_{gln}} \cdot \frac{k_{sat,lac}}{k_{sat,lac} + c_{lac}} \cdot \frac{k_{sat,ammonia}}{k_{sat,ammonia} + c_{ammonia}} \cdot f_{temp} \cdot f_{pH} \cdot f_{inhibition}$$

where $f_{temp}$ and $f_{pH}$ are temperature and pH penalty factors.

## Environmental Effects

### Temperature Effect
The temperature penalty factor uses asymmetric Gaussian functions:

$$f_{temp} = r_{temp} \cdot \begin{cases}
\exp\left(-0.8 \cdot \left(\frac{T - T_{opt}}{3.0}\right)^2\right) & \text{if } T < T_{opt} \\
\exp\left(-1.0 \cdot \left(\frac{T - T_{opt}}{8.0}\right)^2\right) & \text{if } T \geq T_{opt}
\end{cases}$$

where $T_{opt}$ is the optimal temperature and $r_{temp}$ is the temperature robustness factor. The factor is clipped to a minimum of 0.2.

### pH Effect
The pH penalty factor similarly uses asymmetric Gaussian functions:

$$f_{pH} = r_{pH} \cdot \begin{cases}
\exp\left(-0.8 \cdot \left(\frac{pH - pH_{opt}}{1.5}\right)^2\right) & \text{if } pH < pH_{opt} \\
\exp\left(-1.5 \cdot \left(\frac{pH - pH_{opt}}{0.5}\right)^2\right) & \text{if } pH \geq pH_{opt}
\end{cases}$$

where $pH_{opt}$ is the optimal pH and $r_{pH}$ is the pH robustness factor. The factor is clipped to a minimum of 0.2.

Note: For growth rate calculations, different pH and temperature parameters are used than for product formation.

## Death Rate
The specific death rate is:

$$k_{death} = k_{death,max} \cdot \frac{k_{death,\mu}}{\mu + k_{death,\mu}} \cdot s_{scale}$$

## Cell Population Dynamics

$$\frac{dx_{total}}{dt} = \mu \cdot x_{viable} - k_{lysis} \cdot x_{dead}$$

$$\frac{dx_{viable}}{dt} = (\mu - k_{death}) \cdot \max(x_{viable}, 10^{-6})^{g_{inhibition}}$$

$$\frac{dx_{dead}}{dt} = k_{death} \cdot x_{viable} - k_{lysis} \cdot x_{dead}$$

## Metabolite Dynamics

### Product Formation
$$\frac{dc_{product}}{dt} = \left(x_{total} \cdot y_{product,total} + x_{viable} \cdot y_{product,viable} \cdot \mu\right) \cdot f_{pH, product} \cdot f_{temp, product}$$

### Glucose Consumption
$$\frac{dc_{glc}}{dt} = x_{viable} \cdot \left(-\frac{\mu}{y_{cells,glc}} - m_{glc}\right)$$

### Glutamine Dynamics
$$\frac{dc_{gln}}{dt} = x_{viable} \cdot \left(-\frac{\mu}{y_{cells,gln}} - m_{gln}\right) - k_{deg,gln} \cdot c_{gln}$$

### Lactate Production
$$\frac{dc_{lac}}{dt} = -y_{lactate} \cdot \frac{dc_{glc}}{dt}$$

### Ammonia Production
$$\frac{dc_{ammonia}}{dt} = -y_{ammonia} \cdot x_{viable} \cdot \left(-\frac{\mu}{y_{cells,gln}} - m_{gln}\right) + k_{deg,gln} \cdot c_{gln}$$

## Parameters
- $\mu_{max}$: Maximum specific growth rate [`mu_max`]
- $k_{sat,glc}$, $k_{sat,gln}$, $k_{sat,lac}$, $k_{sat,ammonia}$: Monod saturation constants [`k_sat_glc`, `k_sat_gln`, `k_sat_lac`, `k_sat_ammonia`]
- $y_{cells,glc}$, $y_{cells,gln}$: Cell yield coefficients [`y_cells_glc`, `y_cells_gln`]
- $y_{product,total}$, $y_{product,viable}$: Product yield coefficients [`y_product_from_total`, `y_product_from_viable`]
- $y_{lactate}$, $y_{ammonia}$: Byproduct yield coefficients [`y_lactate`, `y_ammonia`]
- $m_{glc}$, $m_{gln}$: Maintenance coefficients [`m_glc`, `m_gln`]
- $k_{death,max}$, $k_{death,\mu}$: Death rate parameters [`k_death_max`, `k_death_mu`]
- $k_{lysis}$: Cell lysis rate [`k_lysis`]
- $k_{deg,gln}$: Glutamine degradation rate [`k_deg_gln`]
- $s_{scale}$: Reactor scale factor [`scale_factor`]
- $g_{inhibition}$: Growth inhibition exponent [`growth_inhibition`]
- $T_{opt}$: Optimal temperature for growth [`e_activation`]
- $r_{temp}$: Temperature robustness factor [`temperature_robustness_factor`]
- $pH_{opt}$: Optimal pH for growth [`ph_opt`]
- $r_{pH}$: pH robustness factor [`ph_robustness_factor`]
- $T_{opt,product}$, $r_{temp,product}$: Product formation temperature parameters [`temperature_opt_product`, `temperature_opt_product_robustness`]
- $pH_{opt,product}$, $r_{pH,product}$: Product formation pH parameters [`ph_opt_product`, `ph_opt_product_robustness`]
