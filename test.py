"""Run the Adamson preprocessing with default settings (see scripts/preprocess_adamson.py)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))

from preprocess_adamson import main  # noqa: E402

if __name__ == "__main__":
    main()
