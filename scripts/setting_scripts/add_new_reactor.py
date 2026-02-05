"""
This script adds a new reactor type to the folder:
src/dbiolab/mammalian_cellculture_simulation/reactor_types
"""

import json
from pathlib import Path

from dbiolab.data_models import ReactorType

path = Path("src/dbiolab/mammalian_cellculture_simulation/reactor_types")


def save_reactor():
    name = "MTP_1"
    reactor = ReactorType(
        name=name,
        volume=3.0,
        scale_factor=5.0,
        growth_inhibition=0.95,
        noise_level=15e-2,
    )
    if not path.exists():
        path.mkdir(parents=True)
    with open(path / f"{name}.json", "w", encoding="utf-8") as f:
        json.dump(reactor.to_dict(), f, indent=4)


if __name__ == "__main__":
    save_reactor()
