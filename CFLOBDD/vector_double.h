#ifndef VECTOR_DOUBLE_GUARD
#define VECTOR_DOUBLE_GUARD


#include <iostream>
#include <string>
#include <fstream>
#include "cflobdd_t.h"

namespace CFL_OBDD {

	typedef CFLOBDD_T<double> CFLOBDD_DOUBLE;

	namespace VectorDouble {
		// Initialization routine
		extern void VectorInitializer();

		extern CFLOBDD_DOUBLE MkBasisVector(unsigned int level, unsigned int index);         // Representation of basis vector with index i
		extern CFLOBDD_DOUBLE VectorToMatrixInterleaved(CFLOBDD_DOUBLE c); // Convert vector to matrix with variables interleaved
		extern CFLOBDD_DOUBLE MatrixToVector(CFLOBDD_DOUBLE c); // Convert vector to matrix with variables interleaved
		extern CFLOBDD_DOUBLE MkColumn1Matrix(unsigned int level);
		extern CFLOBDD_DOUBLE MkVectorWithVoc12(CFLOBDD_DOUBLE c); // convert the vector into another vector with extra volcabularies
		extern CFLOBDD_DOUBLE KroneckerProduct(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2);
		extern CFLOBDD_DOUBLE VectorShiftVocs1To2(CFLOBDD_DOUBLE m1);
		extern CFLOBDD_DOUBLE VectorPadWithZeros(CFLOBDD_DOUBLE c, unsigned int level);
		extern CFLOBDD_DOUBLE NoDistinctionNode(unsigned int level);
#ifdef PATH_COUNTING_ENABLED
		extern std::string Sampling(CFLOBDD_DOUBLE c);
#endif
		extern CFLOBDD_DOUBLE VectorWithAmplitude(CFLOBDD_DOUBLE c);
		extern void VectorPrintColumnMajor(CFLOBDD_DOUBLE c, std::ostream & out);
		extern void VectorPrintColumnMajorInterleaved(CFLOBDD_DOUBLE c, std::ostream & out);
	}
}

#endif

