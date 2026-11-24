# Bioreactor Operator Station

An operator-style dashboard (HMI) for simulated fed-batch runs: one reactor
over time, the way it looks in a control room. It shows

- a vessel schematic with live tags (temperature, pH, DO, pCO2, agitation,
  gas O2 and volume), feed pump activity and the liquid level,
- process values against setpoints, and the latest offline samples,
- trends with setpoints, offline samples and feeds on the timeline,
- the alarm list, the event journal (feeds and setpoint changes), batch KPIs
  and the product titer,
- a run selector and a time cursor with play, to step through the batch,
- a unit detail view (**⤢ Detail** on the vessel panel, Esc to close): the
  vessel large and animated, controller faceplates (PV, SP, deviation,
  output, mode and a 48 h trend) around it, the latest offline analytics and
  the unit's active alarms.

The runs come from the tutorial recipes (best and bad fed-batch), with and
without typical noise, the best recipe scaled to an Ambr 250 and a production
tank from the scale library, and two 20-day perfusion runs from the perfusion
tutorial (`runs.py`):

- **PB-301**: perfusion at 2.5 vvd with a TFF; a PI controller on the
  filtered VCD opens the bleed to hold 50e6 cells/mL (steady state).
- **PB-302**: the perfusion rate follows the cells (CSPR 40 pL/cell/day, no
  bleed) until oxygen and CO2 transfer limit the culture.

Perfusion runs add a TFF retentate loop (culture circulating through the
filter and back into the vessel), the permeate harvest and the bleed (taken
from the loop) to the schematic, with flow animated by rate, perfusion/bleed/harvest tiles, a flow trend, the VCD
setpoint, harvest titer and product collected, and perfusion KPIs (medium
used, volumetric productivity). Alarms are debounced (merged within 3 h,
blips under 1 h ignored), as a control system would.

## Build

From the repository root:

```sh
uv run python hmi_dashboard/build.py
```

This simulates every run and writes
`hmi_dashboard/output/bioreactor_operator_station.html`, a single file with
the data embedded. Open it in a browser. URL parameters select a run and a
time, for example `...station.html?run=bad&t=200`; add `&view=detail` to
open the detail view.

## Files

| file | purpose |
|---|---|
| `runs.py` | run catalogue; simulates a run and exports tags, events and alarms |
| `template.html` | the dashboard page (HTML, CSS, JavaScript, Plotly from a CDN) |
| `build.py` | simulates all runs and fills the template |

To add a run, append a `RunSpec` to `RUNS` in `runs.py`. Alarm limits are in
`ALARMS` in the same file.
