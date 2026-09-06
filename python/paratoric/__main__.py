"""Launch the repository CLI with python -m paratoric."""
from __future__ import annotations

import runpy
import sys
from pathlib import Path


def main() -> None:
    """Execute the sibling CLI script with its local imports available."""
    cli_path = Path(__file__).resolve().parents[1] / "cli" / "paratoric.py"
    sys.path.insert(0, str(cli_path.parent))
    runpy.run_path(str(cli_path), run_name="__main__")


if __name__ == "__main__":
    main()
