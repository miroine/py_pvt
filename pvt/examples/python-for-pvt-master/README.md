![Python package](https://github.com/equinor/python-for-pvt/workflows/Python%20package/badge.svg)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Generic badge](https://img.shields.io/badge/internal_feed-yes-blue.svg)](https://shields.io/)
[![Generic badge](https://img.shields.io/badge/global_analytics-yes-blue.svg)](https://shields.io/)

# python-for-pvt
All the useful PVT formula, equations and empirical relationships you need when working with python.


# Getting started

```
import pvt.equations

bo = pvt.equations.bo()
```

## Installation

There are two approaches to use this package.

Collect it from our internal global analytics python package feed or setup the package yourself.

It is always encoraged to work with environments in python.

### Create environment

In a location that is good to store local environments

```
python -m venv pvt
```

### Activate environment & upgrade pip

Mac / Linux

```
source geochemdb/bin/activate

geochemdb/Scripts/python.exe -m pip install --upgrade pip
```

Windows Git Bash (Use this if working on windows)

```
source ./geochemdb/Scripts/activate

geochemdb/Scripts/python.exe -m pip install --upgrade pip
```

Windows

```
geochemdb\Scripts\activate

\geochemdb\Scripts\python.exe -m pip install --upgrade pip
```

### Connect to local python feed
Follow these steps to connect to the global analytics python feed.

* Step 1
* Step 2

```
pip install pvt
```

### Setup localy using setup.py

Ensure you have cloned and moved into the folder for the repository.

Make the package available locally for this environment using.

```
pip install --upgrade setuptools

python setup.py develop
```

Instead dependancies can be installed

```
pip install -r requirements.txt
```

### Testing

To run the tests

```
pytest tests
```

### Contributions

Contributions are welcome, submit a pull request
