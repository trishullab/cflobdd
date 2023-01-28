#ifndef W_VECTOR_FOURIER_MUL_GUARD
#define W_VECTOR_FOURIER_MUL_GUARD


#include <iostream>
#include <fstream>
#include <random>
#include "fourier_semiring.h"
#include "weighted_cflobdd_t.h"


namespace CFL_OBDD {

	typedef WEIGHTED_CFLOBDD_T<fourierSemiring, std::multiplies<fourierSemiring>> WEIGHTED_CFLOBDD_FOURIER_MUL;

	namespace WeightedVectorFourierMul {
		// Initialization routine
		extern void VectorInitializer();

		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkBasisVector(unsigned int level, unsigned int index);         // Representation of basis vector with index i
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkBasisVector(unsigned int level, std::string s); // index i represented as string
		extern WEIGHTED_CFLOBDD_FOURIER_MUL VectorToMatrixInterleaved(WEIGHTED_CFLOBDD_FOURIER_MUL c); // Convert vector to matrix with variables interleaved
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkColumn1Matrix(unsigned int level);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL KroneckerProduct(WEIGHTED_CFLOBDD_FOURIER_MUL m1, WEIGHTED_CFLOBDD_FOURIER_MUL m2);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL NoDistinctionNode(unsigned int level, fourierSemiring val);
		extern std::string Sampling(WEIGHTED_CFLOBDD_FOURIER_MUL c, bool voctwo, std::string = "");
		extern WEIGHTED_CFLOBDD_FOURIER_MUL VectorWithAmplitude(WEIGHTED_CFLOBDD_FOURIER_MUL c);
		extern void VectorPrintColumnMajor(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out);
		extern void VectorPrintColumnMajorInterleaved(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out);
	}
}

#endif

