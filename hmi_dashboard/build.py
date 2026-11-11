"""Build the operator dashboard: simulate the runs and write one HTML file.

Usage (from the repository root)::

    uv run python hmi_dashboard/build.py [--out PATH]

The run data is embedded in the page, so the file opens directly in a browser.
"""

import argparse
import json
from pathlib import Path

from runs import export_all

HERE = Path(__file__).parent
TEMPLATE = HERE / "template.html"
DEFAULT_OUTPUT = HERE / "output" / "bioreactor_operator_station.html"


def build(output: Path = DEFAULT_OUTPUT) -> Path:
    """Simulate all runs and write the dashboard to ``output``."""
    data = json.dumps(export_all(), separators=(",", ":"), allow_nan=False)
    page = TEMPLATE.read_text(encoding="utf-8").replace("__RUN_DATA__", data)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(f"Wrote {build(args.out)}")


if __name__ == "__main__":
    main()
