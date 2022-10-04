#ifndef VECTOR_TOP_NODE_GUARD
#define VECTOR_TOP_NODE_GUARD


#include <iostream>
#include <fstream>
#include "cflobdd_int.h"

namespace CFL_OBDD {

	namespace VectorInt {

		// Initialization routine
		extern void VectorInitializerTop();

		extern CFLOBDDTopNodeIntRefPtr MkBasisVectorTop(unsigned int level, unsigned int index); // Representation of basis vector at index
		extern CFLOBDDTopNodeIntRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeIntRefPtr n); // Convert vector to matrix with variables interleaved
		extern CFLOBDDTopNodeIntRefPtr MkColumn1MatrixTop(unsigned int level);
		extern CFLOBDDTopNodeIntRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeIntRefPtr n);
		extern CFLOBDDTopNodeIntRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeIntRefPtr n);
		extern CFLOBDDTopNodeIntRefPtr MatrixToVectorTop(CFLOBDDTopNodeIntRefPtr c);
		extern CFLOBDDTopNodeIntRefPtr NoDistinctionNodeTop(unsigned int level);
	}
}

#endif

