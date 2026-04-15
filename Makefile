# Makefile for NESSi documentation
# Runs Doxygen (XML for Breathe) then Sphinx.

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
DOXYFILE      = Doxyfile
DOXYGEN       ?= doxygen

# Default target: build HTML docs (runs Doxygen first)
html: doxygen
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Run Doxygen to regenerate XML used by Breathe
doxygen:
	@$(DOXYGEN) $(DOXYFILE)

# Open Doxygen HTML output in the browser
doxygen-html: doxygen
	@echo "Doxygen HTML output is in doc/html/index.html"
	@open doc/html/index.html 2>/dev/null || xdg-open doc/html/index.html 2>/dev/null || true

# Remove Sphinx build output and Doxygen XML
clean:
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	rm -rf doc/xml doc/html

.PHONY: help html doxygen doxygen-html clean Makefile

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all: pass any other targets (latexpdf, epub, ...) straight to Sphinx
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
