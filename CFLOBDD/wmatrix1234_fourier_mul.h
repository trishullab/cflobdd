#ifndef W_MATRIX1234_FOURIER_MUL_GUARD
#define W_MATRIX1234_FOURIER_MUL_GUARD


#include <iostream>
#include <fstream>
#include <random>
#include "fourier_semiring.h"
#include "weighted_cflobdd_t.h"


namespace CFL_OBDD {

	typedef WEIGHTED_CFLOBDD_T<fourierSemiring, std::multiplies<fourierSemiring>> WEIGHTED_CFLOBDD_FOURIER_MUL;

	namespace WeightedMatrix1234FourierMul {
		// Initialization routine
		extern void Matrix1234Initializer();

		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkIdRelationInterleaved(unsigned int i); // Representation of identity relation
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkWalshInterleaved(unsigned int i);              // Representation of Walsh matrix
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkNegationMatrixInterleaved(unsigned int i);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkCNOTInterleaved(unsigned int i);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkExchangeInterleaved(unsigned int i); // Representation of exchange matrix
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled); // Representation of CNOT matrix with index1 as controller and index2 as controlled bits

		extern WEIGHTED_CFLOBDD_FOURIER_MUL KroneckerProduct2Vocs(WEIGHTED_CFLOBDD_FOURIER_MUL m1, WEIGHTED_CFLOBDD_FOURIER_MUL m2);

		extern void MatrixPrintRowMajor(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out);
		extern void MatrixPrintRowMajorInterleaved(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out);

		extern WEIGHTED_CFLOBDD_FOURIER_MUL MatrixMultiplyV4(WEIGHTED_CFLOBDD_FOURIER_MUL m1, WEIGHTED_CFLOBDD_FOURIER_MUL m2);
        extern WEIGHTED_CFLOBDD_FOURIER_MUL CreateBalancedFn(int n, std::mt19937 mt);

		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkSwapGate(unsigned int i, long int c1, long int c2);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkCPGate(unsigned int i, long int c1, long int c2, fourierSemiring theta);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkCSwapGate(unsigned int i, long int c1, long int x1, long int x2);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkRZGate(unsigned int level, fourierSemiring theta);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkCADDGate(unsigned int level, int c, int x, WEIGHTED_CFLOBDD_FOURIER_MUL f);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkCADDGate2(unsigned int level, int x, WEIGHTED_CFLOBDD_FOURIER_MUL f);
		// specifically for shors
		extern bool CheckIfIndexIsNonZero(unsigned int level, int index, WEIGHTED_CFLOBDD_FOURIER_MUL f);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL MkSetBToZero(unsigned int level, WEIGHTED_CFLOBDD_FOURIER_MUL f);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL ComputeIQFT(unsigned int level, WEIGHTED_CFLOBDD_FOURIER_MUL f, BIG_INT N, int n);
		extern std::pair<WEIGHTED_CFLOBDD_FOURIER_MUL, int> MeasureAndReset(unsigned int level, long int n, WEIGHTED_CFLOBDD_FOURIER_MUL f, fourierSemiring R);
		extern WEIGHTED_CFLOBDD_FOURIER_MUL ResetState(unsigned int level, WEIGHTED_CFLOBDD_FOURIER_MUL f);
	}
}

#endif

