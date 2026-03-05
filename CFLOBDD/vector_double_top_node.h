#ifndef VECTOR_DOUBLE_TOP_NODE_GUARD
#define VECTOR_DOUBLE_TOP_NODE_GUARD

#include <iostream>
#include <string>
#include <fstream>
#include "vector_double.h"
#include "return_map_T.h"

namespace CFL_OBDD {
	typedef ReturnMapBody<double> DoubleReturnMapBody;
	typedef ReturnMapHandle<double> DoubleReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<CFLOBDDReturnMapHandle> Connection;
}

namespace CFL_OBDD {

	typedef CFLOBDDTopNodeT<double> CFLOBDDTopNodeDouble;
	typedef CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeDoubleRefPtr;

	namespace VectorDouble {

		// Initialization routine
		extern void VectorInitializerTop();

		extern CFLOBDDTopNodeDoubleRefPtr MkBasisVectorTop(unsigned int level, unsigned int index); // Representation of basis vector at index
		extern CFLOBDDTopNodeDoubleRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeDoubleRefPtr n); // Convert vector to matrix with variables interleaved
		extern CFLOBDDTopNodeDoubleRefPtr MkColumn1MatrixTop(unsigned int level);
		extern CFLOBDDTopNodeDoubleRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeDoubleRefPtr n);
		extern CFLOBDDTopNodeDoubleRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeDoubleRefPtr n);
		extern CFLOBDDTopNodeDoubleRefPtr MatrixToVectorTop(CFLOBDDTopNodeDoubleRefPtr c);
		extern CFLOBDDTopNodeDoubleRefPtr NoDistinctionNodeTop(unsigned int level);
#ifdef PATH_COUNTING_ENABLED
		extern std::string SamplingTop(CFLOBDDTopNodeDoubleRefPtr n);
#endif
		extern CFLOBDDTopNodeDoubleRefPtr VectorWithAmplitudeTop(CFLOBDDTopNodeDoubleRefPtr n);
		extern void VectorPrintColumnMajorTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out);
		extern void VectorPrintColumnMajorInterleavedTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out);
	}
}

#endif

