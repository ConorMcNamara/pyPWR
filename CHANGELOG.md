# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- CONTRIBUTING.md with comprehensive contribution guidelines
- CHANGELOG.md to track project changes

## [0.1.0] - 2026-02-25

### Added
- Code coverage reporting with codecov integration
- Comprehensive test suite with >85% coverage
- GitHub Actions CI/CD pipeline
- Support for multiple Python versions (3.13+)
- Type hints throughout the codebase (PEP 484/585)
- `py.typed` marker for type checking support

### Changed
- Modernized GitHub Actions workflow
- Updated build configuration for Linux compatibility
- Improved error handling across the package
- Code formatting with ruff and autopep8
- Migrated to pyproject.toml for package configuration

### Fixed
- Linting issues throughout codebase
- Windows-specific error handling in test suite
- Build errors on Linux systems
- T2N test failures on recent Python builds for Windows

### Infrastructure
- Set up automated code formatting with ruff
- Configured docstring formatting with pydocstringformatter (numpydoc style)
- Integrated mypy for static type checking
- Added coverage reporting to CI pipeline
- Configured GitHub Actions for pull request automation

## Initial Release

### Added
- Core power analysis functions:
  - `pwr_p_test()` - One proportion test
  - `pwr_2p_test()` - Two proportions (equal n)
  - `pwr_2p2n_test()` - Two proportions (unequal n)
  - `pwr_t_test()` - One or two sample t-test
  - `pwr_t2n_test()` - Two sample t-test (unequal n)
  - `pwr_anova_test()` - Balanced one-way ANOVA
  - `pwr_r_test()` - Correlation test
  - `pwr_chisq_test()` - Chi-squared test
  - `pwr_f2_test()` - General linear model
  - `pwr_norm_test()` - Mean (normal, known variance)

- Effect size calculation functions:
  - `es_h()` - Cohen's h for proportions
  - `es_w1()` - Cohen's w for goodness-of-fit χ²
  - `es_w2()` - Cohen's w for association χ²
  - `cohen_es()` - Conventional effect sizes

- Comprehensive documentation with usage examples
- Unit tests validating compatibility with R's pwr library
- MIT License
- README with detailed API documentation and examples

[Unreleased]: https://github.com/ConorMcNamara/pyPWR/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/ConorMcNamara/pyPWR/releases/tag/v0.1.0
