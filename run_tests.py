#!/usr/bin/env python3
"""
Test runner script for protqc test suite.

Usage:
    python run_tests.py                    # Run all tests
    python run_tests.py --unit             # Run only unit tests
    python run_tests.py --integration      # Run only integration tests
    python run_tests.py --coverage         # Run with coverage report
    python run_tests.py --fast             # Skip slow tests
    python run_tests.py --verbose          # Verbose output
    python run_tests.py --parallel         # Run tests in parallel
"""
import argparse
import subprocess
import sys
from pathlib import Path


def run_tests(args):
    """Run pytest with specified options."""
    cmd = ["pytest"]
    
    # Add test directory
    cmd.append("tests/")
    
    # Verbosity
    if args.verbose:
        cmd.extend(["-v", "-s"])
    else:
        cmd.append("-v")
    
    # Test selection
    if args.unit:
        cmd.extend(["-m", "unit"])
    elif args.integration:
        cmd.extend(["-m", "integration"])
    elif args.fast:
        cmd.extend(["-m", "not slow"])
    elif args.error_cases:
        cmd.extend(["-m", "error_case"])
    
    # Coverage
    if args.coverage:
        cmd.extend([
            "--cov=src/protqc",
            "--cov-report=html",
            "--cov-report=term-missing",
            "--cov-report=xml"
        ])
    
    # Parallel execution
    if args.parallel:
        cmd.extend(["-n", "auto"])
    
    # Specific test file
    if args.file:
        cmd = ["pytest", f"tests/{args.file}"]
    
    # Specific test function
    if args.test:
        cmd.extend(["-k", args.test])
    
    # Fail fast
    if args.failfast:
        cmd.append("-x")
    
    # Show warnings
    if args.warnings:
        cmd.append("--disable-warnings")
    
    # Benchmark
    if args.benchmark:
        cmd.append("--benchmark-only")
    
    # Print command
    print(f"Running: {' '.join(cmd)}")
    print("=" * 70)
    
    # Run tests
    result = subprocess.run(cmd)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(
        description="Run protqc test suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                          Run all tests
  %(prog)s --unit                   Run only unit tests
  %(prog)s --coverage               Run with coverage report
  %(prog)s --fast                   Skip slow tests
  %(prog)s --file test_helper.py    Run specific test file
  %(prog)s --test test_fetch        Run tests matching pattern
  %(prog)s --parallel               Run tests in parallel
        """
    )
    
    # Test selection
    test_group = parser.add_mutually_exclusive_group()
    test_group.add_argument(
        "--unit",
        action="store_true",
        help="Run only unit tests"
    )
    test_group.add_argument(
        "--integration",
        action="store_true",
        help="Run only integration tests"
    )
    test_group.add_argument(
        "--error-cases",
        dest="error_cases",
        action="store_true",
        help="Run only error case tests"
    )
    test_group.add_argument(
        "--fast",
        action="store_true",
        help="Skip slow tests"
    )
    
    # Output options
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--coverage",
        action="store_true",
        help="Generate coverage report"
    )
    parser.add_argument(
        "--warnings",
        action="store_true",
        help="Show warnings"
    )
    
    # Execution options
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Run tests in parallel"
    )
    parser.add_argument(
        "-x", "--failfast",
        action="store_true",
        help="Stop on first failure"
    )
    
    # Test filtering
    parser.add_argument(
        "--file",
        type=str,
        help="Run specific test file (e.g., test_helper.py)"
    )
    parser.add_argument(
        "--test", "-k",
        type=str,
        help="Run tests matching pattern"
    )
    
    # Performance
    parser.add_argument(
        "--benchmark",
        action="store_true",
        help="Run only benchmark tests"
    )
    
    args = parser.parse_args()
    
    # Check if pytest is installed
    try:
        import pytest
    except ImportError:
        print("Error: pytest is not installed")
        print("Install with: pip install -e '.[dev]'")
        return 1
    
    # Run tests
    return run_tests(args)


if __name__ == "__main__":
    sys.exit(main())
