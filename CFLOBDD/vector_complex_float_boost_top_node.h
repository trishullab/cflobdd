#ifndef VECTOR_COMPLEX_FLOAT_BOOST_TOP_NODE_GUARD
#define VECTOR_COMPLEX_FLOAT_BOOST_TOP_NODE_GUARD

#include <iostream>
#include <string>
#include <fstream>
#include <boost/multiprecision/cpp_complex.hpp>
#include "vector_complex_float_boost.h"
#include "return_map_T.h"
namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef mp::cpp_complex_100 BIG_COMPLEX_FLOAT;
	typedef ReturnMapBody<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<CFLOBDDReturnMapHandle> Connection;
}

namespace CFL_OBDD {

	typedef CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT> CFLOBDDTopNodeComplexFloatBoost;
	typedef CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeComplexFloatBoostRefPtr;

	namespace VectorComplexFloatBoost {

		// Initialization routine
		extern void VectorInitializerTop();

		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkBasisVectorTop(unsigned int level, unsigned int index); // Representation of basis vector at index
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkBasisVectorTop(unsigned int level, std::string s); // Representation of basis vector at index
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n); // Convert vector to matrix with variables interleaved
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkColumn1MatrixTop(unsigned int level);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixToVectorTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr NoDistinctionNodeTop(unsigned int level, BIG_COMPLEX_FLOAT val);
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr ConvertToDoubleTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c);
		extern long double getNonZeroProbabilityTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n);
		extern unsigned long long int GetPathCountTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, long double p);
		//ifdef PATH_COUNTING_ENABLED
		extern std::string SamplingTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, bool voctwo = false);
		extern std::string SamplingV2Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n);
		//#endif
		extern CFLOBDDTopNodeComplexFloatBoostRefPtr VectorWithAmplitudeTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n);
		extern void VectorPrintColumnMajorTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);
		extern void VectorPrintColumnMajorInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);
	}
}

#endif

