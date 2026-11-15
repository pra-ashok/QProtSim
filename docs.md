# QProtSim documentation

The package documentation lives in [docs/](docs/index.md) and is configured in `mkdocs.yml`.

```bash
python -m pip install -e ".[cpu,docs]"
python -m mkdocs serve
```

Build the static site with `python -m mkdocs build --strict` or `make docs`.
