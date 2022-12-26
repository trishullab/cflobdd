#ifndef W_VECTOR_FOURIER_TOP_NODE_GUARD
#define W_VECTOR_FOURIER_TOP_NODE_GUARD

#include <iostream>
#include <string>
#include <fstream>
#include <boost/multiprecision/cpp_complex.hpp>
#include "wvector_fourier_mul.h"
#include "return_map_T.h"
#include "weighted_cflobdd_top_node_t.h"
#include "fourier_semiring.h"

namespace CFL_OBDD {
	typedef ReturnMapBody<fourierSemiring> FourierReturnMapBody;
	typedef ReturnMapHandle<fourierSemiring> FourierReturnMapHandle;
}
namespace CFL_OBDD {

    typedef WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDTopNodeFourier;
	typedef WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr WeightedCFLOBDDTopNodeFourierRefPtr;

	namespace WeightedVectorFourierMul {

		// Initialization routine
		extern void VectorInitializerTop();

		extern WeightedCFLOBDDTopNodeFourierRefPtr MkBasisVectorTop(unsigned int level, unsigned int index); // Representation of basis vector at index
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkBasisVectorTop(unsigned int level, std::string s); // Representation of basis vector at index
		extern WeightedCFLOBDDTopNodeFourierRefPtr VectorToMatrixInterleavedTop(WeightedCFLOBDDTopNodeFourierRefPtr n); // Convert vector to matrix with variables interleaved
		extern WeightedCFLOBDDTopNodeFourierRefPtr MkColumn1MatrixTop(unsigned int level);
		extern WeightedCFLOBDDTopNodeFourierRefPtr NoDistinctionNodeTop(unsigned int level, fourierSemiring val);
		extern std::string SamplingTop(WeightedCFLOBDDTopNodeFourierRefPtr n, bool voctwo = false, std::string = "");
		extern WeightedCFLOBDDTopNodeFourierRefPtr VectorWithAmplitudeTop(WeightedCFLOBDDTopNodeFourierRefPtr n);
		extern void VectorPrintColumnMajorTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out);
		extern void VectorPrintColumnMajorInterleavedTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out);
	}
}

#endif

