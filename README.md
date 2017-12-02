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

## Dependencies
``molsturm`` consists of a set of compiled ``C++`` libraries
as well as a ``python`` module, which the user may use to setup and drive the calculations.

For *building* `molsturm` the following things are required:
- ``cmake`` >= 3.0.0
- A compiler supporting ``C++11``: ``clang`` starting from `clang-3.5` and `gcc` starting
  from `gcc-5` should work.
- [``swig``](http://swig.org/) >= 2.0.11
- [``python``](https://www.python.org/) >= 3.4, including the development headers
- The [regular build process](#building-molsturm) mentioned below
  will automatically build the [``lazyten``](https://lazyten.org) linear algebra library
  as well. This requires further
    - A BLAS implementation, e.g. [OpenBLAS](https://github.com/xianyi/OpenBLAS/)
    - A LAPACK compatible library, e.g. [LAPACK](http://netlib.org/lapack)
    - [armadillo](http://arma.sourceforge.net/)

  See [github.com/lazyten/lazyten](https://github.com/lazyten/lazyten/blob/master/README.md)
  for more details about ``lazyten``'s dependencies.
- If you want to use Gaussian integrals from
  [``libint``](https://github.com/evaleev/libint) you further need
    - Eigen3
    - Autoconf
    - GNU Multiprecision library

In order to actually [use the `molsturm`-module](#using-molsturm-from-python) once
it has been built the following `python` packages are required:
- [h5py](https://pypi.python.org/pypi/h5py)
- [numpy](https://pypi.python.org/pypi/numpy)
- [scipy](https://pypi.python.org/pypi/scipy)
- [PyYAML](https://pypi.python.org/pypi/PyYAML)

On a recent **Debian/Ubuntu** you can install all the aforementioned dependencies by running
```
apt-get install cmake swig python3-dev libopenblas-dev liblapack-dev libarmadillo-dev \
                python3-h5py python3-yaml python3-numpy python3-scipy libeigen3-dev \
                autoconf libgmp-dev
```
as root.

### Optional dependencies
The dependencies in this section are only needed for some extra functionality.
`molsturm` will work without them, but the mentioned features will be unavailable.

- [pyscf](https://github.com/sunqm/pyscf) for full-configuration interaction (FCI)


## Building ``molsturm``
The **recommended build** enables
Edward Valeev's [``libint``](https://github.com/evaleev/libint) library
as an integral backend for Gaussian-type basis functions
and installs ``molsturm`` to `$HOME/opt` of the local user:
```sh
# Fully download molsturm repository and dependent submodules
git clone https://github.com/molsturm/molsturm molsturm
cd molsturm
git submodule update --init --recursive

# Configure inside build directory
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=~/opt -DAUTOCHECKOUT_MISSING_REPOS=ON -DGINT_ENABLE_LIBINT=ON ..
cmake --build .

# Test and install the build
ctest
cmake --build . --target install
```
Note, that this will automatically download and build ``libint`` as well as
[``rapidcheck``](https://github.com/emil-e/rapidcheck) and
[``Catch``](https://github.com/philsquared/Catch/)
as well.
The latter two libraries are only needed for testing `molsturm` and are
not installed to `$HOME/opt`

Installing `molsturm` is **still experimental** and might not work properly
in all cases. Please let us know about any problems.
Every other feedback to improve the process is welcomed as well.

### Configure options
In order to configure the build further there are a number of options,
which can be passed to the first invocation of `cmake` in the instructions above.
For example
- `-DCMAKE_INSTALL_PREFIX=<dir>`: Choose a different installation directory
  for  ``molsturm`` (if you remove the option entirely,
  the default is ``/usr/local``)
- `-DENABLE_TESTS=OFF`: Disable building the unit test executables
- `-DENABLE_EXAMPLES=OFF`: Disable building the example executables
- `-DENABLE_DOCUMENTATION=ON`: Build and install the *sparse*
  [doxygen](http://www.stack.nl/~dimitri/doxygen/index.html)-generated
  in-source documentation.
- `-DGINT_ENABLE_LIBCINT=ON`: Enable Gaussian integrals via the
  [``libcint``](https://github.com/sunqm/libcint) library.
- `-DGINT_ENABLE_STURMINT=ON`: Enable the Coulomb-Sturmian basis functions
  via [`sturmint`](https://molsturm.org/sturmint)
  (This library is not yet publicly available, but will be released soon.)


## Using `molsturm` from `python`
Given that you followed the build procedure sketched above,
you just need to make sure that your `python` interpreter finds `molsturm`
in the installation directory.
This is done by modifying the `PYTHONPATH`.
For example if you chose the `CMAKE_INSTALL_PREFIX` of `~/opt`
as we shown above, you need to type
```
# Find out the major and minor version of your python interpreter:
VERSION=$(python3 -c 'import sys; print(str(sys.version_info.major) + "." + str(sys.version_info.minor))')
export PYTHONPATH="$PYTHONPATH:$HOME/opt/lib/python${VERSION}/site-packages"
```
to setup your environment for `molsturm`.

You might need to change `python3` to `python` in the first line
if your operating system uses `python` to refer to the binary of your
`python` version 3 interpreter.

### Running calculations
To get started with `molsturm` take a look at the [examples](examples/) subfolder.
Good entry points are especially the [single_point](examples/single_point) scripts.

In principle a calculation is a sequence of calls
to `python` functions of the `molsturm` module.
For example the snippet
```python
import molsturm
hfres = molsturm.hartree_fock("Be", basis_type="gaussian",
                              basis_set_name="6-31g")
print("Be HF energy", hfres["energy_ground_state"])
```
just performs a Hartree-Fock calculation on a beryllium atom and
prints the resulting SCF energy.

