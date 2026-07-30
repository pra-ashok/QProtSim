# setup.py
from setuptools import setup, find_packages

setup(
    name='protqc',
    version='0.5.0',
    author='Pratik Patil, Bhushan Bonde',
    author_email='p.patil@uos.ac.uk, b.bonde@uos.ac.uk',
    description='A package for thermodynamic modeling of proteins with Quantum computers',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/pratikp204/protqc',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={},
    license="MIT",
    license_files=['LICENSE'],
    python_requires='>=3.11',
    install_requires=[
        'numpy',
        'MDAnalysis',
        "biopython>=1.86",
        'pandas>=2.2.2',
        'requests>=2.32.3',
        # 'torch>=2.5.0',
        'OpenMM>=8.3.0',
        'scipy>=1.11.4',
        ##'mdakit-sasa==0.2.7',  removed for now
        'freesasa>=2.2.0',
        'netCDF4>=1.7.1',
        'rich>=13.9.2',
        'plotly>=5.20.0',
        'ipywidgets>=8.1.0',
        'zstandard>=0.25.0',
        'pennylane==0.44.0',
        'pennylane_lightning==0.44.0',
        'pennylane_lightning_gpu==0.44.0',
        # Note: cupy must be installed via conda: conda install cupy cuda-version=13.1 -c conda-forge
        'reportlab',
        'pyyaml',
        'seaborn >=0.13.2',
        'pillow',   
        'pymol-open-source',  # Molecular visualization for structure rendering
        'psutil',
        'pennylane-qiskit'
    ],
    extras_require={
        "dev": [
            # Testing framework
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",          # Coverage reporting
            "pytest-mock>=3.11.1",        # Enhanced mocking
            "pytest-timeout>=2.1.0",      # Test timeout handling
            "pytest-xdist>=3.3.1",        # Parallel test execution
            "pytest-benchmark>=4.0.0",    # Performance benchmarking
            
            # Mocking and fixtures
            "responses>=0.23.0",          # Mock HTTP requests
            "freezegun>=1.2.0",           # Mock datetime
            "faker>=19.0.0",              # Generate fake data
            
            # Code quality
            "black>=23.7.0",              # Code formatting
            "flake8>=6.1.0",              # Linting
            "mypy>=1.5.0",                # Type checking
            "pylint>=2.17.0",             # Code analysis
            "isort>=5.12.0",              # Import sorting
            
            # Documentation
            "sphinx>=7.1.0",              # Documentation generation
            "sphinx-rtd-theme>=1.3.0",    # ReadTheDocs theme
            
            # Development tools
            "ipython>=8.14.0",            # Interactive shell
            "ipdb>=0.13.13",              # Debugging
            "pre-commit>=3.3.0",          # Git hooks
            
            # Performance profiling
            "memory-profiler>=0.61.0",    # Memory usage profiling
            "line-profiler>=4.0.0",       # Line-by-line profiling
        ]
    },
    entry_points={
        'console_scripts': [
            'protqc=protqc.__main__:main',
        ],
    },
)
