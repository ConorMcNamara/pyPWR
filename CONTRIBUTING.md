# Contributing to pyPWR

Thank you for your interest in contributing to pyPWR! This document provides guidelines and instructions for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Development Workflow](#development-workflow)
- [Code Standards](#code-standards)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Reporting Bugs](#reporting-bugs)
- [Requesting Features](#requesting-features)

## Code of Conduct

By participating in this project, you agree to maintain a respectful and inclusive environment for all contributors.

## Getting Started

Before contributing, please:

1. Check the [existing issues](https://github.com/ConorMcNamara/pyPWR/issues) to see if your concern has already been raised
2. For major changes, open an issue first to discuss what you would like to change
3. Ensure your contribution aligns with the project's goals of providing a Pythonic statistical power analysis package

## Development Setup

### Prerequisites

- Python 3.13 or higher
- Git
- pip

### Installation

1. Fork the repository on GitHub
2. Clone your fork locally:

```bash
git clone https://github.com/YOUR_USERNAME/pyPWR.git
cd pyPWR
```

3. Install the package in development mode with all dependencies:

```bash
pip install -e ".[dev]"
```

This will install:
- Core dependencies (NumPy, SciPy)
- Development dependencies (pytest, ruff, mypy, coverage tools)

4. Verify the installation:

```bash
pytest test/ -v
```

## Development Workflow

### Creating a Branch

Create a new branch for your feature or bugfix:

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bugfix-name
```

### Making Changes

1. Make your changes in the appropriate files
2. Add or update tests as necessary
3. Update documentation if you're changing functionality
4. Ensure your code follows the project's code standards (see below)

### Code Standards

pyPWR maintains high code quality standards. Before committing, ensure your code passes all checks:

#### Formatting

Format your code with ruff:

```bash
ruff format .
```

#### Linting

Check for linting issues:

```bash
ruff check .
```

To automatically fix some linting issues:

```bash
ruff check . --fix
```

#### Type Checking

Run mypy to ensure type correctness:

```bash
mypy PyPWR
```

All public functions should include type hints following PEP 484/585.

#### Running All Checks

Run all quality checks at once:

```bash
ruff format . && ruff check . && mypy PyPWR && pytest --cov=PyPWR
```

## Testing

### Running Tests

Run the full test suite:

```bash
pytest test/ -v
```

### Test Coverage

Check code coverage:

```bash
pytest test/ --cov=PyPWR --cov-report=term-missing
```

Generate an HTML coverage report:

```bash
pytest test/ --cov=PyPWR --cov-report=html
# Open htmlcov/index.html in your browser
```

### Writing Tests

- All new functionality must include tests
- Aim for >85% code coverage
- Tests should validate compatibility with R's pwr library where applicable (differences should be < 1e-03)
- Place tests in the `test/` directory
- Use descriptive test names that explain what is being tested

Example test structure:

```python
def test_pwr_t_test_two_sample():
    """Test two-sample t-test power calculation."""
    result = pwr_t_test(
        d=0.5,
        n=64,
        sig_level=0.05,
        type='two-sample',
        alternative='two-sided'
    )
    assert result['power'] > 0.80
    assert abs(result['power'] - 0.8014) < 1e-03
```

## Submitting Changes

### Commit Messages

Write clear, descriptive commit messages:

- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit the first line to 72 characters or less
- Reference issues and pull requests when relevant

Example:

```
Add power calculation for paired t-test

- Implement pwr_t_test with type='paired'
- Add validation for paired sample sizes
- Include tests comparing results with R's pwr package

Fixes #123
```

### Pull Request Process

1. Update your branch with the latest main:

```bash
git checkout main
git pull upstream main
git checkout your-branch
git rebase main
```

2. Push your branch to your fork:

```bash
git push origin your-branch
```

3. Open a Pull Request on GitHub with:
   - A clear title describing the change
   - A detailed description of what changed and why
   - References to related issues
   - Screenshots or examples if applicable

4. Ensure all CI checks pass:
   - Tests pass on all supported Python versions
   - Code coverage remains >85%
   - Linting and type checking pass

5. Address any review feedback promptly

6. Once approved, a maintainer will merge your PR

## Reporting Bugs

When reporting bugs, please include:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected behavior
- Actual behavior
- Python version and operating system
- Minimal code example demonstrating the issue
- Any relevant error messages or stack traces

Use the [GitHub issue tracker](https://github.com/ConorMcNamara/pyPWR/issues) to report bugs.

## Requesting Features

For feature requests, please:

- Check if the feature has already been requested
- Provide a clear description of the feature
- Explain why this feature would be useful
- Include examples of how the feature would be used
- Consider whether it aligns with the project's goal of matching R's pwr functionality

## Questions?

If you have questions about contributing, feel free to:

- Open an issue on GitHub
- Check existing documentation in the [README](README.md)

Thank you for contributing to pyPWR!
