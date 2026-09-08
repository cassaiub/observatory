# Convenience wrappers. `make install` is equivalent to ./install.sh.
#
# Everything here assumes the environment is already active, except `install`,
# which builds it.

PYTHON ?= python

.PHONY: help install test test-offline lint format doctor docs clean

help:
	@echo "make install       one-command setup (delegates to ./install.sh)"
	@echo "make test          run the full test suite"
	@echo "make test-offline  run the tests that need no network"
	@echo "make lint          ruff check"
	@echo "make format        ruff check --fix"
	@echo "make doctor        check this machine can run the pipeline"
	@echo "make docs          rebuild the LaTeX manual"
	@echo "make clean         remove build and cache artifacts"

install:
	./install.sh

test:
	$(PYTHON) -m pytest

test-offline:
	$(PYTHON) -m pytest -m "not network"

lint:
	$(PYTHON) -m ruff check src tests

format:
	$(PYTHON) -m ruff check --fix src tests

doctor:
	$(PYTHON) -m cassa_photometry.doctor

docs:
	cd ../docs && latexmk -pdf main.tex

clean:
	rm -rf build dist .pytest_cache .ruff_cache htmlcov .coverage
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.py[co]' -delete
