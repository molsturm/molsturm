//
// Copyright (C) 2017 by the molsturm authors
//
// This file is part of molsturm.
//
// molsturm is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// molsturm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with molsturm. If not, see <http://www.gnu.org/licenses/>.
//

%module molsturm_iface

%{
#include "interface/python/ExcInvalidParameters.hh"
#include "interface/python/hartree_fock.hh"
#include "interface/python/HfResults.hh"
#include "interface/python/HfParameters.hh"
#include <molsturm/Version.hh>

// Run the %init block (to setup numpy)
#define SWIG_FILE_WITH_INIT
%}

// TODO Extremely rudimentary exception handling
%exception {
  try {
    $action
  } catch (const krims::ExcNotImplemented& e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    return NULL;
  } catch (const molsturm::iface::ExcInvalidParameters& e) {
    PyErr_SetString(PyExc_RuntimeError, e.extra().c_str());
    return NULL;
  } catch (const krims::ExceptionBase& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }
}

%include "numpy.i"
%include "std_string.i"
%include "std_vector.i"

// Instantiate what is needed
namespace std {
  %template(DoubleVector)  vector<double>;
}

// Setup import of numpy array:
%init %{
import_array();
%}

%include "HfParameters.hh"
%include "HfResults.hh"
%include "hartree_fock.hh"
%include "../../molsturm/Version.hh"

// vi: syntax=c
