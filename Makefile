.PHONY: test test-unit test-integration test-error test-fast test-coverage test-parallel clean install install-dev help

# Default target
.DEFAULT_GOAL := help

help:  ## Show this help message
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install package dependencies
	pip install -e .

install-dev:  ## Install package and development dependencies
	pip install -e ".[dev]"

test:  ## Run all tests
	@pip show qprotsim > /dev/null 2>&1 || (echo "Error: qprotsim not installed. Run 'make install' or 'make install-dev' first." && exit 1)
	python -m pytest tests/ -v

test-unit:  ## Run only unit tests
	python -m pytest tests/ -m unit -v

test-integration:  ## Run only integration tests
	python -m pytest tests/ -m integration -v

test-error:  ## Run only error case tests
	python -m pytest tests/ -m error_case -v

test-fast:  ## Run tests, skip slow ones
	python -m pytest tests/ -m "not slow" -v

test-coverage:  ## Run tests with coverage report
	python -m pytest tests/ --cov=src/qprotsim --cov-report=html --cov-report=term-missing

test-coverage-xml:  ## Run tests with XML coverage (for CI)
	python -m pytest tests/ --cov=src/qprotsim --cov-report=xml

test-parallel:  ## Run tests in parallel
	python -m pytest tests/ -n auto -v

test-verbose:  ## Run tests with maximum verbosity
	python -m pytest tests/ -vv -s

test-failfast:  ## Stop at first test failure
	python -m pytest tests/ -x -v

test-specific:  ## Run specific test file (usage: make test-specific FILE=test_helper.py)
	python -m pytest tests/$(FILE) -v

test-watch:  ## Run tests in watch mode (requires pytest-watch)
	python -m ptw tests/

coverage-report:  ## Open HTML coverage report
	@echo "Opening coverage report..."
	@python -c "import webbrowser; webbrowser.open('htmlcov/index.html')"

clean:  ## Clean up test artifacts
	rm -rf .pytest_cache
	rm -rf htmlcov
	rm -rf .coverage
	rm -rf coverage.xml
	rm -rf *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

lint:  ## Run code linters
	flake8 src/qprotsim tests/
	pylint src/qprotsim

format:  ## Format code with black
	black src/qprotsim tests/

format-check:  ## Check code formatting
	black --check src/qprotsim tests/

type-check:  ## Run type checking with mypy
	mypy src/qprotsim

quality:  ## Run all quality checks
	make format-check
	make lint
	make type-check

ci:  ## Run full CI pipeline locally
	make clean
	make format-check
	make lint
	make test-coverage

benchmark:  ## Run performance benchmarks
	pytest tests/ --benchmark-only

profile:  ## Profile test execution
	pytest tests/ --profile

# Development helpers
shell:  ## Start IPython shell with package loaded
	ipython -i -c "import sys; sys.path.insert(0, 'src'); from qprotsim import helper"

jupyter:  ## Start Jupyter notebook
	jupyter notebook

docs:  ## Build documentation
	python -m mkdocs build --strict

docs-serve:  ## Serve documentation locally
	python -m mkdocs serve
