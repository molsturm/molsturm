# molsturm
[![Build Status](https://travis-ci.org/molsturm/molsturm.svg?branch=master)](https://travis-ci.org/molsturm/molsturm)
[![Licence](https://img.shields.io/github/license/molsturm/molsturm.svg)](LICENCE)

``molsturm`` is a modular electronic structure theory program,
which tries to employ concepts like lazy matrices and orthogonal programming
to achieve a versatile structure.
Our motto is *flexibility first, but speed second*, even though we are probably very far
from either of the two still.

Note that this is a very young project.
Things may change in incompatible ways at any time.

## Build dependencies
``molsturm`` consists of a set of compiled ``C++`` libraries
as well as a ``python`` module, which the user may use to setup and drive the calculations.

The build setup configures and builds most of these components automatically
provided that *all* of the following requirements are met:
- ``cmake`` >= 3.0.0
- A compiler supporting ``C++11``: ``clang`` starting from `clang-3.5` and `gcc` starting
  from `gcc-5` should work.
  `gcc-4.8` support is intended, but not quite there yet.
- To build the compiled part of the ``molsturm`` python module:
	- [``swig``](http://swig.org/) >= 2.0.11
	- [``python``](https://www.python.org/) >= 3.4, including the development headers
- To build the [``linalgwrap``](https://linalgwrap.org) linear algebra library
  you will need
	- A BLAS implementation, e.g. [OpenBLAS](https://github.com/xianyi/OpenBLAS/)
	- A LAPACK compatible library, e.g. [LAPACK](http://netlib.org/lapack)
	- [armadillo](http://arma.sourceforge.net/)
  See [github.com/linalgwrap/linalgwrap](https://github.com/linalgwrap/linalgwrap/blob/master/README.md)
  for more details about ``linalgwrap``'s dependencies.

On a recent **Debian/Ubuntu** the following should install the required packages
```
apt-get install cmake swig python3-dev libopenblas-dev liblapack-dev libarmadillo-dev
```
provided that your standard compiler is recent enough.

## Building ``molsturm``
The **recommended build** enables the Gaussian integral backend
via Edward Valeev's [``libint``](https://github.com/evaleev/libint) library.

In order to achieve this, just execute
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON -DGINT_ENABLE_LIBINT=ON ..
cmake --build .
ctest
```
Note that this will automatically download and build ``libint`` as well as
[``rapidcheck``](https://github.com/emil-e/rapidcheck) and
[``Catch``](https://github.com/philsquared/Catch/)
(both required for testing) along with ``molsturm``.

To build without tests (not recommended), run
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON -DGINT_ENABLE_LIBINT=ON \
  -DMOLSTURM_ENABLE_TESTS=OFF -DGINT_ENABLE_TESTS=OFF -DGSCF_ENABLE_TESTS=OFF \
  -DLINALGWRAP_ENABLE_TESTS=OFF -DKRIMS_ENABLE_TESTS=OFF ..
cmake --build .
```

Note that installing `molsturm` is not yet fully functional.

## Using `molsturm` from `python`
In order to use the `molsturm` module the following other `python`
packages are required:
- [h5py](https://pypi.python.org/pypi/h5py)
- [numpy](https://pypi.python.org/pypi/numpy)
- [PyYAML](https://pypi.python.org/pypi/PyYAML)

In a recent **Debian/Ubuntu** the following should do:
```
apt-get install python3-h5py python3-yaml python3-numpy
```

To install them via **pip** run
```
pip3 install h5py numpy pyyaml
```

`molsturm` has been successfully used on systems with
`python3` starting from version 3.4.

### Running calculations
To get started with using molsturm take a look at the [examples](examples/)
subfolder. A good starting point are especially the [single_point](examples/single_point)
scripts. For example run
```
# Setup environment to run molsturm
cd examples
. setup_environment.sh

# Run a water sto-3g calculation
./single_point/water_sto3g.py
```
