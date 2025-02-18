SPHINX_OPTS := -j auto

BUILD_DIR := $(CURDIR)/_build

.PHONY: all
all: html pdf

.PHONY: cp
cp:
	mkdir -p src/ src/docs
	cp ../*.md src/docs/
	cp ../../*.md src/

.PHONY: html
html: cp # apidoc  refs.bibtex.bib
	sphinx-build -M html $(CURDIR) $(BUILD_DIR) $(SPHINX_OPTS)

.PHONY: pdf
pdf: cp # apidoc  refs.bibtex.bib
	sphinx-build -M latexpdf $(CURDIR) $(BUILD_DIR) $(SPHINX_OPTS)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) src

.PHONY: distclean
distclean: clean
