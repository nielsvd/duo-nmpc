// SQP interface

%module sqp_interface
%{
#include "ocp.hpp"
#include "transverse_ocp.hpp"
#include "tangential_ocp.hpp"
#include "duo_nmpc.hpp"
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const char e []) {
    SWIG_exception(SWIG_RuntimeError, e);
  }
}

%include "stl.i"
%include "std_vector.i"
// Instantiate templates used by example
namespace std {
   %template(VectorInt) vector<int>;
   %template(VectorDouble) vector<double>;
   %template(DoubleVectorVector) vector< vector<double> >;
   %template(StringVector) vector<string>;
}

%include ocp.hpp
%include transverse_ocp.hpp
%include tangential_ocp.hpp
%include duo_nmpc.hpp