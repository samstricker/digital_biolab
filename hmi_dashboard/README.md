# Bioreactor Operator Station

An operator-style dashboard (HMI) for simulated fed-batch runs: one reactor
over time, the way it looks in a control room. It shows

- a vessel schematic with live tags (temperature, pH, DO, pCO2, agitation,
  gas O2 and volume), feed pump activity and the liquid level,
- process values against setpoints, and the latest offline samples,
- trends with setpoints, offline samples and feeds on the timeline,
- the alarm list, the event journal (feeds and setpoint changes), batch KPIs
  and the product titer,
- a run selector and a time cursor with play, to step through the batch.

The runs come from the tutorial recipes (best and bad fed-batch), with and
without typical noise, and the best recipe scaled to an Ambr 250 and a
production tank from the scale library (`runs.py`).

## Build

From the repository root:

```sh
uv run python hmi_dashboard/build.py
```

This simulates every run and writes
`hmi_dashboard/output/bioreactor_operator_station.html`, a single file with
the data embedded. Open it in a browser. URL parameters select a run and a
time, for example `...station.html?run=bad&t=200`.

## Files

| file | purpose |
|---|---|
| `runs.py` | run catalogue; simulates a run and exports tags, events and alarms |
| `template.html` | the dashboard page (HTML, CSS, JavaScript, Plotly from a CDN) |
| `build.py` | simulates all runs and fills the template |

To add a run, append a `RunSpec` to `RUNS` in `runs.py`. Alarm limits are in
`ALARMS` in the same file.
