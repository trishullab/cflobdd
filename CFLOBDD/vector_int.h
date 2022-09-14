#ifndef VECTOR_INT_GUARD
#define VECTOR_INT_GUARD


#include <iostream>
#include <fstream>
#include "cflobdd_int.h"

namespace CFL_OBDD {

	namespace VectorInt {

		// Initialization routine
		extern void VectorInitializer();

		extern CFLOBDD MkBasisVector(unsigned int level, unsigned int index);         // Representation of basis vector with index i
		extern CFLOBDD VectorToMatrixInterleaved(CFLOBDD c); // Convert vector to matrix with variables interleaved
		extern CFLOBDD MatrixToVector(CFLOBDD c); // Convert vector to matrix with variables interleaved
		extern CFLOBDD MkColumn1Matrix(unsigned int level);
		extern CFLOBDD MkVectorWithVoc12(CFLOBDD c); // convert the vector into another vector with extra volcabularies
		extern CFLOBDD KroneckerProduct(CFLOBDD m1, CFLOBDD m2);
		extern CFLOBDD VectorShiftVocs1To2(CFLOBDD m1);
		extern CFLOBDD VectorPadWithZeros(CFLOBDD c, unsigned int level);
		extern CFLOBDD NoDistinctionNode(unsigned int level);
	}
}

#endif