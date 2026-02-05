"""
This script adds a new cell type to the folder:
src/dbiolab/mammalian_cellculture_simulation/celltypes
"""

import json
from pathlib import Path

from dbiolab.mammalian_cellculture_simulation.data_models import CellType

path = Path("src/dbiolab/mammalian_cellculture_simulation/celltypes")


def save_celltype():
    name = "craven_3L_batch"
    celltype = CellType(
        name=name,
        mu_max=0.035,
        k_lysis=0.04,
        k_deg_gln=0.001,
        k_death_max=0.01,
        k_death_mu=0.01,
        k_sat_lac=150.0,
        k_sat_ammonia=40.0,
        k_sat_glc=1.0,
        k_sat_gln=0.22,
        y_ammonia=0.68,
        y_lactate=1.6,
        y_cells_glc=9.23e7,
        y_cells_gln=8.8e8,
        m_glc=8.0e-13,
        m_gln=3.0e-12,
        ph_opt=7.5,
        e_activation=38.5,
        temperature_robustness_factor=1.0,
        ph_robustness_factor=1.0,
        y_product_from_total=1.0e-15,
        y_product_from_viable=1.0e-9,
    )
    if not path.exists():
        path.mkdir(parents=True)
    with open(path / f"{name}.json", "w", encoding="utf-8") as f:
        json.dump(celltype.to_dict(), f, indent=4)


if __name__ == "__main__":
    save_celltype()
