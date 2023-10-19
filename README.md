![Python versions](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10%20%7C%203.11-blue)
[![code style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# GetRandomPCMol
Python program for obtaining random PubChem molecules and generating conformer ensembles for them (optionally).

## Dependencies

`getrandompcmol` in its current state depends on existing installations of 
- the `xtb` program package for optimization and MDs/MTDs (for details and installation instructions, see [github.com/grimme-lab/xtb](https://github.com/grimme-lab/xtb))
- the `mctc-convert` executable that can be obtained from the `mctc-lib` library (for details and installation instructions, see [github.com/grimme-lab/mctc-lib](https://github.com/grimme-lab/mctc-lib))
- the `PubGrep` script for downloading PubChem molecules based on their CID (for details, see [github.com/grimme-lab/PubGrep](https://github.com/grimme-lab/PubGrep))
- _if conformer search is desired_: `crest` for creating conformer ensembles from a given structure file (for details, see [github.com/crest-lab/crest](https://github.com/grimme-lab/PubGrep)](https://github.com/crest-lab/crest))

## Installation

After cloning the code via `git clone git@github.com:grimme-lab/NumgradPy.git`, a new virtual `conda` environment with required Python pre-requisites can be set up with:
```
conda env create -f environment.yml
conda activate numgradpy
```
`numgradpy` can be installed into this environment with 
```
pip install -e .
```
The flag `-e` allows modification of the code without the need to reinstall.

## Use

After installation, the package can generate a data set of random molecules. The randomization is based on a random number generator in `NumPy`.
```
getrandompcmol -n 500 --opt --crest --maxnumat 35 --maxcid 1000000 --seed 27051997 --evalconf 5 10
```
The distinct keywords are described as follows:
- `-n/--n <int>`: controls the number of molecules generated.
- `--maxnumat <int>`: Maximally allowed number of atoms per molecule.
- `--maxcid <int>`: Range between 1 and `<maxcid>`, in which random CIDs are generated.
- - `--opt`: Optimize the molecules using GFN2-xTB in `xtb`.
- `--crest`: Generate conformer ensembles from the given structures via `crest`.
- `--seed <int>`: Starting seed for random number generation.
- `--evalconf <int> <int>`: Range of conformer ensemble size allowed for post-processing.

If only the evaluation and post-processing of a previously generated conformer ensemble is desired, use ``--evalconfonly``. Further information on possible input flags is available via `--help`.

CIDs and names of generated molecules are saved into a file `compounds.txt`. Additionally, (if `--crest`) is active, each directory contains a file `conformer.json`, in which relevant properties of the generated ensemble are saved.

The `crest` conformer ensemble generation is executed in parallel. The parallelization adapts to the number of available cores on your machine and the requested number of molecules.

## Source code

The source code is in the [src/getrandompcmol](src/getrandompcmol) directory. Here, also some _dunder_ files can be found:

- [\_\_version\_\_.py](src/getrandompcmol/__version__.py): just the version number as a string, used by config files
- [\_\_init\_\_.py](src/getrandompcmol/__init__.py): entry point for program/library

<be>

## Setup files and Packaging

Packaging is done with [`setuptools`](https://setuptools.pypa.io/en/latest/index.html), which is configured through the `pyproject.toml` and/or `setup.cfg`/`setup.py` files.

<details>
<summary>
  <code>pyproject.toml</code> vs.
  <code>setup.cfg</code> vs
  <code>setup.py</code>
</summary>

The `setup.py` file is a Python script, and configuration is passed through keyword arguments of `setuptools.setup()`. This is not recommended due to possible security and parsing issues. The same setup can be accomplished in a declarative style within `setup.cfg`, and `setup.py` remains mostly empty only calling `setuptools.setup()`.
The `pyproject.toml` file aims to unify configuration files including various tools like black or pytest. For packaging, it is very similar to `setup.cfg`. However, `pyproject.toml` has not been adopted as the default yet, and many projects still use `setup.cfg` to declare the packaging setup. Note that `setup.py` is not necessary if a `pyproject.toml` is present.

</details>

#### `pyproject.toml`

- minimal build specification to use with setuptools
- configuration of other tools (black, pytest, mypy, ...)

[](https://setuptools.pypa.io/en/latest/userguide/declarative_config.html#using-a-src-layout)

#### `setup.cfg`

- declarative configuration for setuptools
- [_metadata_](https://setuptools.pypa.io/en/latest/userguide/declarative_config.html#metadata): must at least contain _name_ and _version_
- [_options_](https://setuptools.pypa.io/en/latest/userguide/declarative_config.html#options): package discovery, dependencies
  - [additional setup](https://setuptools.pypa.io/en/latest/userguide/declarative_config.html#using-a-src-layout) required for `src/` layout
- [_options.extras_require_](https://setuptools.pypa.io/en/latest/userguide/dependency_management.html#optional-dependencies): optional dependencies (dev tools, docs, ...)
- [_options.package_data_](https://setuptools.pypa.io/en/latest/userguide/datafiles.html#package-data): inclusion of other, non-Python files (marker files, data, ...)
  - alternative: `MANIFEST.in`
- [_options.entry_points_](https://setuptools.pypa.io/en/latest/userguide/entry_point.html): entry point for command line interface
- can also hold configuration of other tools

<br>

The package can be installed with `pip install .` or something like `pip install . [dev]` to also install additional dependencies specified in `setup.cfg`'s _options.extras_require_. Pass the `-e` flag for editable mode, which loads the package from the source directory, i.e., changing the source code does not require a new installation.
