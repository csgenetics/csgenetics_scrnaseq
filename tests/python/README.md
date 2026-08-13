# Python unit tests for `bin/`

Behavioural pytest coverage of the deterministic Python logic in the pipeline's
`bin/` scripts (the existing nf-tests only cover Nextflow wiring).

## Running

Always run via pixi, never bare `pytest`. From the repository root:

```bash
pixi run --manifest-path tests/requirements/python/pixi.toml pytest
```

This auto-discovers `pytest.ini` (`testpaths = tests`) and runs the whole
suite. To select by marker:

```bash
# pure-logic tests only
pixi run --manifest-path tests/requirements/python/pixi.toml pytest -m unit
# subprocess / anndata round-trip tests only
pixi run --manifest-path tests/requirements/python/pixi.toml pytest -m integration
```

## Environments

* `tests/requirements/python/pixi.toml` - the main test env (pytest + the modern
  numpy/pandas/scipy/anndata/scanpy/plotly/jinja2/pysam stack the scripts need).
* `tests/requirements/gtfparse/pixi.toml` - a small pinned env mirroring
  `conda_envs/gtfparse.yml` (gtfparse 1.2.1 / pandas <2.0). Only
  `test_features_names.py` uses it, because `features_names.py` depends on the
  pandas DataFrame API that newer polars-based gtfparse releases break. Install
  it once with:

  ```bash
  pixi install --manifest-path tests/requirements/gtfparse/pixi.toml
  ```

  If it is not installed, `test_features_names.py` skips with a clear message.

## Markers

* `unit` - pure-logic tests, no heavy fixtures.
* `integration` - build/round-trip anndata or h5ad, or run a `bin/` script as a
  subprocess on crafted tiny inputs.
