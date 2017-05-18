// vi: syntax=c

%module molsturm_iface

%{
#include "Parameters.hh"
#include "HfResults.hh"
#include "hartree_fock.hh"
#include <krims/ExceptionSystem/Exceptions.hh>
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_array.i"

// Instantiate what is needed
namespace std {
  %template(DoubleVector)        vector<double>;
  %template(Coord)               array<double,3>;
  %template(CoordVector)         vector<array<double,3>>;
  %template(AtomicNumberVector)  vector<unsigned int>;
  %template(Nlm)                 array<int,3>;
  %template(NlmBasis)            vector<array<int, 3>>;
}

// TODO Extremely rudimentary exception handling
%exception {
  try {
    $action
  } catch (krims::ExcNotImplemented& e) {
    PyErr_SetString(PyExc_NotImplementedError, e.what());
    return NULL;
  } catch (krims::ExceptionBase& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }
}

%include "Parameters.hh"
%include "HfResults.hh"
%include "hartree_fock.hh"
