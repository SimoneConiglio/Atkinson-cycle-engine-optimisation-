# Convenience targets. Everything here is also runnable via tox.
.DEFAULT_GOAL := help
PY ?= python3
VENV := .venv
BIN := $(VENV)/bin
export MPLBACKEND := Agg

.PHONY: help venv install dev test test-all lint format type figures animation clean

help:  ## show this help
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

venv:  ## create a virtualenv with uv (falls back to python -m venv)
	@command -v uv >/dev/null 2>&1 \
		&& uv venv $(VENV) \
		|| $(PY) -m venv $(VENV)

install: venv  ## install the package and its runtime dependencies
	@command -v uv >/dev/null 2>&1 \
		&& uv pip install --python $(BIN)/python -e "." \
		|| $(BIN)/pip install -e "."

dev: venv  ## install with the dev and multi-objective extras
	@command -v uv >/dev/null 2>&1 \
		&& uv pip install --python $(BIN)/python -e ".[dev,moea]" \
		|| $(BIN)/pip install -e ".[dev,moea]"

test:  ## run the fast tests
	$(BIN)/pytest -m "not slow"

test-all:  ## run every test, optimizers included
	$(BIN)/pytest --cov --cov-report=term-missing

lint:  ## check style
	$(BIN)/ruff check src tests examples
	$(BIN)/ruff format --check src tests examples

format:  ## reformat in place
	$(BIN)/ruff format src tests examples
	$(BIN)/ruff check --fix src tests examples

type:  ## static type check
	$(BIN)/mypy src/exlink

figures:  ## regenerate the static figures
	$(BIN)/exlink plot -o figures

animation:  ## regenerate the animated GIF
	$(BIN)/exlink animate -o figures/exlink.gif --frames 120

clean:  ## remove build and test artefacts
	rm -rf build dist .pytest_cache .mypy_cache .ruff_cache .coverage htmlcov
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.egg-info' -type d -prune -exec rm -rf {} +
