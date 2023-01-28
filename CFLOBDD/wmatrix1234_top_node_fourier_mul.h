#ifndef W_MATRIX1234_FOURIER_MUL_TOP_NODE_GUARD
#define W_MATRIX1234_FOURIER_MUL_TOP_NODE_GUARD

#include <iostream>
#include <fstream>
#include "fourier_semiring.h"
#include "wmatrix1234_fourier_mul.h"
#include "return_map_T.h"

namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef ReturnMapBody<fourierSemiring> FourierReturnMapBody;
	typedef ReturnMapHandle<fourierSemiring> FourierReturnMapHandle;
	typedef ReturnMapBody<int> CFLOBDDReturnMapBody;
	typedef ReturnMapHandle<int> CFLOBDDReturnMapHandle;
}

namespace CFL_OBDD {

	typedef WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDTopNodeFourier;
	typedef WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr WeightedCFLOBDDTopNodeFourierRefPtr;
	

	namespace WeightedMatrix1234FourierMul {

		// Initialization routine
		extern void Matrix1234InitializerTop();

		extern WeightedCFLOBDDTopNodeFourierRefPtr MkIdRelationInterleavedTop(unsigned int i); // Representation of identity relation
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkWalshInterleavedTop(unsigned int i); // Representation of Walsh matrix
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkNegationMatrixInterleavedTop(unsigned int i);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkCNOTInterleavedTop(unsigned int i);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkExchangeInterleavedTop(unsigned int i); // Representation of exchange matrix
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled);
		extern WeightedCFLOBDDTopNodeFourierRefPtr KroneckerProduct2VocsTop(WeightedCFLOBDDTopNodeFourierRefPtr m1, WeightedCFLOBDDTopNodeFourierRefPtr m2);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled);

		extern void MatrixPrintRowMajorTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out);
		extern void MatrixPrintRowMajorInterleavedTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out);

		extern WeightedCFLOBDDTopNodeFourierRefPtr MatrixMultiplyV4TopNode(WeightedCFLOBDDTopNodeFourierRefPtr c1, WeightedCFLOBDDTopNodeFourierRefPtr c2);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkSwapGateTop(unsigned int level, long int i, long int j);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkCPGateTop(unsigned int level, long int i, long int j, fourierSemiring theta);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkCSwapGateTop(unsigned int level, long int c, long int i, long int j);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkRZGateTop(unsigned int level, fourierSemiring theta);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkCADDGateTop(unsigned int level, int c, int x, WeightedCFLOBDDTopNodeFourierRefPtr f);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkCADDGate2Top(unsigned int level, int x, WeightedCFLOBDDTopNodeFourierRefPtr f);
		extern bool CheckIfIndexIsNonZeroTop(unsigned int level, int index, WeightedCFLOBDDTopNodeFourierRefPtr f);
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkSetBToZeroTop(unsigned int level, WeightedCFLOBDDTopNodeFourierRefPtr f);
		extern WeightedCFLOBDDTopNodeFourierRefPtr ComputeIQFTTop(unsigned int level, WeightedCFLOBDDTopNodeFourierRefPtr f, BIG_INT N, int n);
		extern std::pair<WeightedCFLOBDDTopNodeFourierRefPtr, int> MeasureAndResetTop(unsigned int level, long int n, WeightedCFLOBDDTopNodeFourierRefPtr f, fourierSemiring R);
		extern WeightedCFLOBDDTopNodeFourierRefPtr ResetStateTop(unsigned int level, WeightedCFLOBDDTopNodeFourierRefPtr f);
	}
}

#endif

