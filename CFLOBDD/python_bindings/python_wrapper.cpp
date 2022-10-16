#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "../tests_cfl.h"
#include "../cflobdd_int.h"

namespace py = pybind11;


namespace CFL_OBDD{

	PYBIND11_MODULE(python_wrapper_cflobdd, m) {
	    m.doc() = "python wrapper for CFLOBDD"; // Optional module docstring

	    py::class_<CFLTests>(m, "CFLTests")
		.def("runtest", &CFLTests::runTests, "runTest")
		.def("init", &CFLTests::InitModules, "init")
		.def("clear", &CFLTests::ClearModules, "clear");

	    py::class_<CFLOBDD_T<int>>(m, "CFLOBDD")
		.def("MkProjection", &MkProjection, "projection")
		.def("MkAnd", py::overload_cast<CFLOBDD, CFLOBDD>(&MkAnd))
		.def("MkOr", py::overload_cast<CFLOBDD, CFLOBDD>(&MkOr))
		.def("MkNot", &MkNot)
		.def("MkFalse", &MkFalse)
		.def("MkTrue", &MkTrue)
		.def("Print", &PrintCFLOBDD)
		.def("compute_prob", &ComputeProbability);
	
	}
}
