.PHONY: help install install-dev test test-cov lint format type-check clean build all

# Default target
.DEFAULT_GOAL := help

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install project dependencies
	uv sync --no-dev

install-dev:  ## Install project dependencies including dev dependencies
	uv sync

test:  ## Run tests without coverage
	uv run pytest -v

test-cov:  ## Run tests with coverage report
	uv run pytest

lint:  ## Run ruff linter
	uv run ruff check PyPWR test

lint-fix:  ## Run ruff linter and fix issues automatically
	uv run ruff check --fix PyPWR test

format:  ## Format code with ruff
	uv run ruff format PyPWR test

format-check:  ## Check code formatting without making changes
	uv run ruff format --check PyPWR test

type-check:  ## Run zuban type checker
	uv run zuban check PyPWR

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
	rm -rf .zuban_cache/
	rm -rf .ruff_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete

build:  ## Build the package
	uv build

publish:  ## Publish package to PyPI (requires authentication)
	uv publish

publish-test:  ## Publish package to TestPyPI
	uv publish --publish-url https://test.pypi.org/legacy/

all: clean install-dev check-all build  ## Run full CI pipeline (clean, install, check, test, build)
