# Digital BioLab

A digital twin of CHO cell culture bioprocesses: simulate batch, fed-batch and
perfusion runs in-silico to test optimisation and control strategies before
applying them in the lab.

## Contents

| Path | What it is |
|---|---|
| [`src/dbiolab/mammalian_cellculture_simulation`](src/dbiolab/mammalian_cellculture_simulation/ReadMe.md) | The simulator: cell growth and death, metabolism, product formation and quality, stirred-tank engineering (oxygen, CO2, mixing, shear), batch / fed-batch / perfusion, controllers, realistic noise and a library of reactor scales and cell types. |
| [`hmi_dashboard`](hmi_dashboard/README.md) | An operator-style dashboard (control-room HMI) to replay simulated runs. |
| `tests` | Unit and regression tests (`pytest`). |

## Setup

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```sh
uv sync                 # install the package and its dependencies
uv run pytest           # run the tests
```

## Quick start

```python
from dbiolab.mammalian_cellculture_simulation import MammalianCellcultureReactor

simulator = MammalianCellcultureReactor(reactor_type)    # a ReactorType
result = simulator.simulate(recipe, seed=1)              # a MammalianRecipe
result.true.c_product        # what happened in the reactor
result.measured.c_product    # what the analytics report (with noise)
```

The notebooks in
[`setup_and_tutorial`](src/dbiolab/mammalian_cellculture_simulation/setup_and_tutorial)
show complete examples:

| Notebook | Shows |
|---|---|
| `tutorial.ipynb` | A fed-batch run step by step: reactor, cell type, recipe, results, setpoints; a good and a bad recipe. |
| `perfusion_tutorial.ipynb` | Perfusion with cell retention, the built-in controllers, writing your own control loop, verification against Richelle et al. (2022). |
| `scale_library.ipynb` | Sampling realistic reactors from Ambr to production scale and CHO cell types; the same process at every scale. |
| `setup.ipynb` | Writing the input JSON files. |
| `../research/noise_analysis.ipynb` | What each noise source does, with reasonable ranges for every setting. |

## Operator dashboard

```sh
uv run python hmi_dashboard/build.py
open hmi_dashboard/output/bioreactor_operator_station.html
```
