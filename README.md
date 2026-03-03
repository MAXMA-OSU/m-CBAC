# mcbac

`mcbac` packages the m-CBAC workflow as an installable Python project.

## What it does

- Bundles the original CBAC lookup databases.
- Provides a Python API and CLI.
- Processes CIF files from an input directory and writes `FINAL_*.cif` outputs.
- Runs fully in Python (no shell/awk/sed/dos2unix/gfortran required at runtime).

## Install

```bash
pip install -e .
```

## CLI

```bash
mcbac run --data-dir DATA --output-dir FINAL_DATA --log-file m-cbac.log
```

Useful options:

- `--work-dir`: use a persistent working directory instead of a temp dir.
- `--keep-work-dir`: keep temporary files for debugging.

## Python API

```python
from pathlib import Path
from mcbac.pipeline import run_pipeline

result = run_pipeline(
    data_dir=Path("DATA"),
    output_dir=Path("FINAL_DATA"),
    log_path=Path("m-cbac.log"),
)
print(result.final_files)
```

## Runtime requirements

- Python 3.9+