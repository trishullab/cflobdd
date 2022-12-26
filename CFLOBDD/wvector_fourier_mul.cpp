#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "wvector_fourier_mul.h"
#include "wvector_top_node_fourier_mul.h"

namespace CFL_OBDD {

	namespace WeightedVectorFourierMul {

		void VectorInitializer()
		{
			VectorInitializerTop();
			return;
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkBasisVector(unsigned int level, unsigned int index)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkBasisVectorTop(level, index));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkBasisVector(unsigned int level, std::string s)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkBasisVectorTop(level, s));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkColumn1Matrix(unsigned int level)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkColumn1MatrixTop(level));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL VectorToMatrixInterleaved(WEIGHTED_CFLOBDD_FOURIER_MUL c)
		{
			WEIGHTED_CFLOBDD_FOURIER_MUL tempNode = WEIGHTED_CFLOBDD_FOURIER_MUL(VectorToMatrixInterleavedTop(c.root));
			WEIGHTED_CFLOBDD_FOURIER_MUL v = MkColumn1Matrix(tempNode.root->level);
			return tempNode * v;
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL KroneckerProduct(WEIGHTED_CFLOBDD_FOURIER_MUL m1, WEIGHTED_CFLOBDD_FOURIER_MUL m2)
		{
			assert(m1.root->level == m2.root->level);
			WEIGHTED_CFLOBDD_FOURIER_MUL m2_1To2 = m2;
			return m1 * m2_1To2;
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL NoDistinctionNode(unsigned int level, fourierSemiring val)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(NoDistinctionNodeTop(level, val));
		}

		std::string Sampling(WEIGHTED_CFLOBDD_FOURIER_MUL c, bool isTwoVoc, std::string func)
		{
			return SamplingTop(c.root, isTwoVoc, func);
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL VectorWithAmplitude(WEIGHTED_CFLOBDD_FOURIER_MUL c)
		{
			return VectorWithAmplitudeTop(c.root);
		}

		void VectorPrintColumnMajor(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out)
		{
			VectorPrintColumnMajorTop(c.root, out);
			return;
		}

		void VectorPrintColumnMajorInterleaved(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out)
		{
			VectorPrintColumnMajorInterleavedTop(c.root, out);
			return;
		}
	}
}

