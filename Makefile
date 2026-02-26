.PHONY: help install install-dev test test-cov lint format type-check clean build all

# Default target
.DEFAULT_GOAL := help

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install project dependencies
	poetry install --without dev

install-dev:  ## Install project dependencies including dev dependencies
	poetry install

test:  ## Run tests without coverage
	poetry run pytest -v

test-cov:  ## Run tests with coverage report
	poetry run pytest

lint:  ## Run ruff linter
	poetry run ruff check PyPWR test

lint-fix:  ## Run ruff linter and fix issues automatically
	poetry run ruff check --fix PyPWR test

format:  ## Format code with ruff
	poetry run ruff format PyPWR test

format-check:  ## Check code formatting without making changes
	poetry run ruff format --check PyPWR test

type-check:  ## Run mypy type checker
	poetry run mypy PyPWR

check: lint type-check  ## Run all checks (lint + type-check)

check-all: format-check lint type-check test  ## Run all checks and tests

clean:  ## Clean up build artifacts and cache files
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf coverage.xml
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete

build:  ## Build the package
	poetry build

publish:  ## Publish package to PyPI (requires authentication)
	poetry publish

publish-test:  ## Publish package to TestPyPI
	poetry publish -r testpypi

all: clean install-dev check-all build  ## Run full CI pipeline (clean, install, check, test, build)
