#include <pybind11/pybind11.h>
#include "../tests_cfl.h"
#include "../cflobdd_int.h"

namespace py = pybind11;

PYBIND11_MODULE(python_wrapper_cflobdd, m) {
    m.doc() = "python wrapper for CFLOBDD"; // Optional module docstring

    py::class_<CFL_OBDD::CFLTests>(m, "CFLTests")
    	.def("runtest", &CFL_OBDD::CFLTests::runTests, "ada");

    py::class_<CFL_OBDD::CFLOBDD_T<int>>(m, "CFLOBDD")
    	.def("projection", &CFL_OBDD::MkProjection, "projection");
}
