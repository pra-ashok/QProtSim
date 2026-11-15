# Build and maintain the documentation

## Install the documentation tools

From the repository root, install the package and documentation extra.

```bash
python -m pip install -e ".[cpu,docs]"
python -m mkdocs serve
```

The preview normally opens at `http://127.0.0.1:8000`. The Makefile also provides `make docs-serve`.

For a documentation-only environment, install the tools directly. The site uses hand-maintained API pages and builds without importing the simulation dependencies.

```bash
python -m pip install "mkdocs>=1.6.1,<2" "mkdocs-shadcn>=0.12.3,<0.13" "pymdown-extensions>=10,<11" "pygments>=2.18,<3"
python -m mkdocs build --strict
```

The generated site is written to `site/`, which is ignored by Git. `make docs` performs the same strict build. Serve the result with `python -m http.server --directory site 8000` to inspect the generated HTML.

## Organisation

| Location | Purpose |
| --- | --- |
| `mkdocs.yml` | Grouped sidebar navigation, Shadcn theme, Markdown extensions, and validation |
| `docs/getting-started.md` | A complete first-run example |
| `docs/configuration.md` | Supported YAML keys and CLI behaviour |
| `docs/api.md` | Public functions, examples, and return types |
| `docs/theory.md` | Equations, assumptions, and method references |
| `docs/outputs.md` | Artifact schemas and interpretation |
| `docs/assets/diagrams/` | Reusable academic-style vector figures |
| `docs/assets/stylesheets/extra.css` | Typography, colour, and figure styling |
| `overrides/templates/` | Site branding and author-credit footer |
| `hooks/docs.py` | Support strict builds from source archives without Git history |

Shadcn keeps the documentation sections visible in the left sidebar on desktop and provides a menu on smaller screens. Search, light and dark themes, code copying, and a page contents panel are supplied by the theme. KaTeX, Geist fonts, and syntax styles are included in the generated site. Math options belong to `theme.katex_options` in `mkdocs.yml`.

