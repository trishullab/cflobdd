#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "vector_complex_float_boost_top_node.h"
#include "vector_complex_float_boost.h"
#include "matrix1234_complex_float_boost.h"

namespace CFL_OBDD {

	namespace VectorComplexFloatBoost {

		void VectorInitializer()
		{
			VectorInitializerTop();
			return;
		}

		CFLOBDD_COMPLEX_BIG MkBasisVector(unsigned int level, unsigned int index)
		{
			return CFLOBDD_COMPLEX_BIG(MkBasisVectorTop(level, index));
		}

		CFLOBDD_COMPLEX_BIG MkBasisVector(unsigned int level, std::string s)
		{
			return CFLOBDD_COMPLEX_BIG(MkBasisVectorTop(level, s));
		}

		CFLOBDD_COMPLEX_BIG MkColumn1Matrix(unsigned int level)
		{
			return CFLOBDD_COMPLEX_BIG(MkColumn1MatrixTop(level));
		}

		CFLOBDD_COMPLEX_BIG VectorToMatrixInterleaved(CFLOBDD_COMPLEX_BIG c)
		{
			CFLOBDD_COMPLEX_BIG tempNode = CFLOBDD_COMPLEX_BIG(VectorToMatrixInterleavedTop(c.root));
			CFLOBDD_COMPLEX_BIG v = MkColumn1Matrix(tempNode.root->level);
			return tempNode * v;
		}

		CFLOBDD_COMPLEX_BIG MkVectorWithVoc12(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(MkVectorWithVoc12Top(c.root));
		}

		CFLOBDD_COMPLEX_BIG KroneckerProduct(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_COMPLEX_BIG m2_1To2 = VectorShiftVocs1To2(m2);
			return m1 * m2_1To2;
		}

		CFLOBDD_COMPLEX_BIG VectorShiftVocs1To2(CFLOBDD_COMPLEX_BIG m1)
		{
			return CFLOBDD_COMPLEX_BIG(VectorShiftVocs1To2Top(m1.root));
		}

		CFLOBDD_COMPLEX_BIG NoDistinctionNode(unsigned int level, BIG_COMPLEX_FLOAT val)
		{
			return CFLOBDD_COMPLEX_BIG(NoDistinctionNodeTop(level, val));
		}


		//#ifdef PATH_COUNTING_ENABLED
		std::string Sampling(CFLOBDD_COMPLEX_BIG c, bool isTwoVoc)
		{
			return SamplingTop(c.root, isTwoVoc);
		}

		std::string SamplingV2(CFLOBDD_COMPLEX_BIG c)
		{
			return SamplingV2Top(c.root);
		}
		//#endif

		CFLOBDD_COMPLEX_BIG VectorWithAmplitude(CFLOBDD_COMPLEX_BIG c)
		{
			return CFLOBDD_COMPLEX_BIG(VectorWithAmplitudeTop(c.root));
		}

		void VectorPrintColumnMajor(CFLOBDD_COMPLEX_BIG c, std::ostream & out)
		{
			VectorPrintColumnMajorTop(c.root, out);
			return;
		}

		void VectorPrintColumnMajorInterleaved(CFLOBDD_COMPLEX_BIG c, std::ostream & out)
		{
			VectorPrintColumnMajorInterleavedTop(c.root, out);
			return;
		}

		long double getNonZeroProbability(CFLOBDD_COMPLEX_BIG n)
		{
			return getNonZeroProbabilityTop(n.root);
		}

		unsigned long long int GetPathCount(CFLOBDD_COMPLEX_BIG n, long double p)
		{
			return GetPathCountTop(n.root, p);
		}
	}
}

