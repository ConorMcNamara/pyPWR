# pyPWR

[![CI](https://github.com/ConorMcNamara/pyPWR/actions/workflows/ci.yml/badge.svg)](https://github.com/ConorMcNamara/pyPWR/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ConorMcNamara/pyPWR/branch/main/graph/badge.svg)](https://codecov.io/gh/ConorMcNamara/pyPWR)
[![Python 3.13+](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Typed with mypy](https://img.shields.io/badge/typed-mypy-blue.svg)](https://github.com/python/mypy)

A comprehensive Python package for statistical power analysis, providing a Pythonic equivalent to R's popular [pwr](https://github.com/heliosdrm/pwr) library.

## Overview

Despite the proliferation of A/B testing and experimental design in data science, Python lacks a unified package for calculating statistical power across different test types. **pyPWR** fills this gap by providing a complete suite of power analysis tools based on Cohen's seminal work on statistical power analysis.

### Key Features

- 🔬 **Comprehensive Test Coverage**: Power analysis for t-tests, ANOVA, chi-squared, correlation, proportions, and general linear models
- 📊 **Effect Size Calculations**: Built-in functions for computing standardized effect sizes (Cohen's d, h, w, f, f²)
- 🔄 **Flexible Analysis**: Calculate any of the four parameters (effect size, sample size, significance level, or power) given the other three
- ✅ **Validated Results**: Unit tests ensure compatibility with R's pwr library (differences < 1e-03)
- 🐍 **Pythonic API**: Clean, modern Python 3.13+ codebase with type hints and comprehensive documentation
- 🚀 **High Performance**: Leverages NumPy and SciPy for efficient numerical computations

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Available Functions](#available-functions)
- [Usage Examples](#usage-examples)
- [Effect Size Calculations](#effect-size-calculations)
- [Important Notes](#important-notes)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)

## Installation

### From PyPI (when available)

```bash
pip install pypwr
```

### From Source

```bash
git clone https://github.com/ConorMcNamara/pyPWR.git
cd pyPWR
pip install -e .
```

### Requirements

- Python 3.13 or higher
- NumPy >= 2.2.4
- SciPy >= 1.15.2

## Quick Start

```python
from PyPWR import pwr_2p_test

# Calculate power for a two-proportion test
result = pwr_2p_test(
    h=0.3,              # Effect size
    n=200,              # Sample size per group
    sig_level=0.05,     # Significance level
    alternative='greater'
)

print(f"Power: {result['power']:.4f}")
# Output: Power: 0.9123
```

**Key Principle**: Exactly one of `effect_size`, `n` (or `n1`/`n2`), `sig_level`, or `power` must be `None`. The function will solve for the missing parameter.

## Available Functions

### Power Test Functions

| Function | Description | Parameters |
|----------|-------------|------------|
| `pwr_p_test()` | One proportion test | `h`, `n`, `sig_level`, `power`, `alternative` |
| `pwr_2p_test()` | Two proportions (equal n) | `h`, `n`, `sig_level`, `power`, `alternative` |
| `pwr_2p2n_test()` | Two proportions (unequal n) | `h`, `n1`, `n2`, `sig_level`, `power`, `alternative` |
| `pwr_t_test()` | One or two sample t-test | `d`, `n`, `sig_level`, `power`, `type`, `alternative` |
| `pwr_t2n_test()` | Two sample t-test (unequal n) | `d`, `n1`, `n2`, `sig_level`, `power`, `alternative` |
| `pwr_anova_test()` | Balanced one-way ANOVA | `k`, `n`, `f`, `sig_level`, `power` |
| `pwr_r_test()` | Correlation test | `r`, `n`, `sig_level`, `power`, `alternative` |
| `pwr_chisq_test()` | Chi-squared test | `w`, `n`, `df`, `sig_level`, `power` |
| `pwr_f2_test()` | General linear model | `u`, `v`, `f2`, `sig_level`, `power` |
| `pwr_norm_test()` | Mean (normal, known variance) | `d`, `n`, `sig_level`, `power`, `alternative` |

### Effect Size Functions

| Function | Description | Returns |
|----------|-------------|---------|
| `es_h()` | Effect size for proportions | Cohen's h |
| `es_w1()` | Effect size for goodness-of-fit χ² | Cohen's w |
| `es_w2()` | Effect size for association χ² | Cohen's w |
| `cohen_es()` | Conventional effect sizes | Small/medium/large values |

## Usage Examples

### Example 1: Calculate Required Sample Size

Determine how many participants you need for a two-sample t-test:

```python
from PyPWR import pwr_t_test

result = pwr_t_test(
    d=0.5,              # Medium effect size
    power=0.80,         # 80% power
    sig_level=0.05,     # 5% significance
    type='two-sample',
    alternative='two-sided'
)

print(f"Required sample size per group: {result['n']}")
# Output: Required sample size per group: 64
```

### Example 2: Calculate Statistical Power

Check if your current sample size provides adequate power:

```python
from PyPWR import pwr_anova_test

result = pwr_anova_test(
    k=4,                # 4 groups
    n=25,               # 25 per group
    f=0.25,             # Small-to-medium effect
    sig_level=0.05
)

print(f"Statistical power: {result['power']:.2%}")
# Output: Statistical power: 68.71%
```

### Example 3: Determine Minimum Detectable Effect

Find what effect size you can detect with your constraints:

```python
from PyPWR import pwr_r_test

result = pwr_r_test(
    n=100,              # Sample size
    power=0.80,         # Desired power
    sig_level=0.05,
    alternative='two-sided'
)

print(f"Minimum detectable correlation: {result['effect_size']:.3f}")
# Output: Minimum detectable correlation: 0.277
```

### Example 4: Chi-Squared Test with Unequal Sample Sizes

```python
from PyPWR import pwr_2p2n_test, es_h

# Calculate effect size from proportions
p1, p2 = 0.60, 0.45
effect_size = es_h(p1, p2)

result = pwr_2p2n_test(
    h=effect_size,
    n1=100,
    n2=150,
    sig_level=0.05,
    alternative='two-sided'
)

print(f"Effect size h: {effect_size:.3f}")
print(f"Power: {result['power']:.2%}")
```

### Example 5: Using Conventional Effect Sizes

```python
from PyPWR import cohen_es, pwr_t_test

# Get conventional "medium" effect size for t-test
effect = cohen_es(test='t', size='medium')
print(f"Medium effect size for t-test: {effect['effect_size']}")  # 0.5

# Use in power calculation
result = pwr_t_test(
    d=effect['effect_size'],
    n=30,
    sig_level=0.05,
    type='paired',
    alternative='two-sided'
)

print(f"Power with n=30 and d=0.5: {result['power']:.2%}")
```

## Effect Size Calculations

### Proportions (h)

Cohen's h for the difference between two proportions:

```python
from PyPWR import es_h

h = es_h(p1=0.60, p2=0.45)
print(f"Effect size h: {h:.3f}")  # 0.304
```

### Chi-Squared Tests (w)

For goodness-of-fit tests:

```python
from PyPWR import es_w1

# Null hypothesis: equal proportions [0.25, 0.25, 0.25, 0.25]
# Alternative: [0.375, 0.208, 0.208, 0.208]
p0 = [0.25, 0.25, 0.25, 0.25]
p1 = [0.375, 0.208, 0.208, 0.208]

w = es_w1(p0, p1)
print(f"Effect size w: {w:.3f}")
```

For contingency table tests:

```python
from PyPWR import es_w2

# 2x4 contingency table
prob = [[0.225, 0.125, 0.125, 0.125],
        [0.160, 0.160, 0.040, 0.040]]

w = es_w2(prob)
print(f"Effect size w: {w:.3f}")
```

### Conventional Effect Sizes

```python
from PyPWR import cohen_es

# Get conventional effect sizes for different tests
tests = ['t', 'r', 'anova', 'chisq', 'f2']
sizes = ['small', 'medium', 'large']

for test in tests:
    print(f"\n{test.upper()} test:")
    for size in sizes:
        es = cohen_es(test=test, size=size)
        print(f"  {size}: {es['effect_size']}")
```

## Important Notes

### Parameter Requirements

- **Exactly one parameter must be `None`**: The function will solve for the missing value
- All probability values (significance level, power) must be between 0 and 1
- Sample sizes must be positive integers
- Effect sizes are typically positive (direction specified by `alternative`)

### Alternative Hypotheses

Most functions support three alternatives:
- `'two-sided'`: Tests for any difference (default)
- `'greater'`: Tests for positive effect
- `'less'`: Tests for negative effect

### Computational Methods

pyPWR uses SciPy's `brentq` and `bisect` for root-finding, while R's pwr uses `uniroot`. This may result in minor numerical differences (typically < 1e-03), which are within acceptable tolerance for practical applications.

### Type Hints and IDE Support

The package includes full type hints (PEP 484/585) and a `py.typed` marker for excellent IDE support and static type checking.

## Code Coverage

pyPWR maintains high test coverage (>85%) to ensure reliability and correctness. Coverage reports are automatically generated for each commit and uploaded to [Codecov](https://codecov.io/gh/ConorMcNamara/pyPWR).

### Coverage Breakdown

- **PyPWR/__init__.py**: Public API exports
- **PyPWR/effect_size.py**: Effect size calculation functions
- **PyPWR/power_classes.py**: Core power calculation classes
- **PyPWR/pwr_tests.py**: User-facing test functions (>98% coverage)

### Viewing Coverage Reports

```bash
# Run tests with coverage
pytest test/ --cov=PyPWR --cov-report=term-missing

# Generate interactive HTML report
pytest test/ --cov=PyPWR --cov-report=html
open htmlcov/index.html  # macOS
# or: xdg-open htmlcov/index.html  # Linux
# or: start htmlcov/index.html     # Windows
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

### Development Setup

```bash
git clone https://github.com/ConorMcNamara/pyPWR.git
cd pyPWR
pip install -e ".[dev]"  # Install with development dependencies
```

### Running Tests

```bash
# Run tests
pytest test/ -v

# Run tests with coverage report
pytest test/ -v --cov=PyPWR --cov-report=term-missing

# Generate HTML coverage report
pytest test/ --cov=PyPWR --cov-report=html
# Open htmlcov/index.html in your browser
```

### Code Quality

```bash
# Format code
ruff format .

# Lint code
ruff check .

# Type check
mypy PyPWR

# Run all checks (formatting, linting, type checking, tests with coverage)
ruff format . && ruff check . && mypy PyPWR && pytest --cov=PyPWR
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates, Publishers.

## Acknowledgments

This package is inspired by and designed to be compatible with R's [pwr](https://github.com/heliosdrm/pwr) package. While maintaining API similarity, pyPWR has been redesigned to follow Python conventions and best practices.

---

**Questions or Issues?** Please [open an issue](https://github.com/ConorMcNamara/pyPWR/issues) on GitHub.
