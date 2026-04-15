# NESSi Documentation

This directory contains the source for the NESSi documentation, built with
[Doxygen](https://www.doxygen.nl/) and [Sphinx](https://www.sphinx-doc.org/)
using the [Read the Docs](https://sphinx-rtd-theme.readthedocs.io/) theme.
Doxygen parses the C++ source and generates XML, which Sphinx ingests via the
[Breathe](https://breathe.readthedocs.io/) extension to produce the final HTML.

## Prerequisites

### Doxygen

Install Doxygen via your system package manager:

```bash
# Ubuntu / Debian
sudo apt install doxygen

# macOS (Homebrew)
brew install doxygen
```

### Python packages

Python 3.8 or newer is required. Install the dependencies into a virtual
environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The `requirements.txt` installs:

| Package | Purpose |
|---|---|
| `sphinx` | Documentation generator |
| `sphinx-rtd-theme` | Read the Docs HTML theme |
| `breathe` | Bridge between Doxygen XML and Sphinx |
| `myst-parser` | Allows Markdown files to be included in Sphinx |

## Building the documentation

With the virtual environment active, run from this directory:

```bash
make html
```

This will:
1. Run Doxygen on the C++ source in `cntr/`, writing XML and HTML to `doc/`
2. Run Sphinx, reading the RST sources in `sphinx_docs/` and producing the
   final HTML site in `_build/html/`

Open `_build/html/index.html` in a browser to view the result.

### Other targets

```bash
make doxygen        # run Doxygen only (generates XML and HTML)
make doxygen-html   # build and open the pure source-code documentation (classes, functions, etc.) in the browser
make clean          # remove all build artifacts (_build/ and doc/)
```

The `doxygen-html` target produces a standalone API reference generated directly
from the C++ source, including all classes, functions, and their documentation.
The output is written to `doc/html/index.html`.

## Directory structure

```
nessi_doc/
├── Makefile            # top-level build entry point
├── Doxyfile            # Doxygen configuration
├── conf.py             # Sphinx configuration
├── index.rst           # Sphinx root document
├── requirements.txt    # Python dependencies
├── cntr/               # C++ library headers (Doxygen input)
├── sphinx_docs/        # RST source pages and examples
├── _static/            # Images, CSS, and other static assets
└── _templates/         # Sphinx HTML templates
```
