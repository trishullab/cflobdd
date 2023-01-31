#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "wvector_complex_fb_mul.h"
#include "wvector_top_node_complex_fb_mul.h"

namespace CFL_OBDD {

	namespace WeightedVectorComplexFloatBoostMul {

		void VectorInitializer()
		{
			VectorInitializerTop();
			return;
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkBasisVector(unsigned int level, unsigned int index, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkBasisVectorTop(level, index, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkBasisVector(unsigned int level, std::string s, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkBasisVectorTop(level, s, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkColumn1Matrix(unsigned int level)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkColumn1MatrixTop(level));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL VectorToMatrixInterleaved(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c)
		{
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL tempNode = WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(VectorToMatrixInterleavedTop(c.root));
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL v = MkColumn1Matrix(tempNode.root->level);
			return tempNode * v;
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL KroneckerProduct(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m1, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m2)
		{
			assert(m1.root->level == m2.root->level);
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m2_1To2 = m2;
			return m1 * m2_1To2;
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL NoDistinctionNode(unsigned int level, BIG_COMPLEX_FLOAT val)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(NoDistinctionNodeTop(level, val));
		}

		std::string Sampling(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c, bool isTwoVoc, std::mt19937 mt, std::uniform_real_distribution<double> dis, std::string func)
		{
			return SamplingTop(c.root, mt, dis, isTwoVoc, func);
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL VectorWithAmplitude(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c)
		{
			return VectorWithAmplitudeTop(c.root);
		}

		void VectorPrintColumnMajor(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c, std::ostream & out)
		{
			VectorPrintColumnMajorTop(c.root, out);
			return;
		}

		void VectorPrintColumnMajorInterleaved(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c, std::ostream & out)
		{
			VectorPrintColumnMajorInterleavedTop(c.root, out);
			return;
		}

		long double getNonZeroProbability(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL n)
		{
			return getNonZeroProbabilityTop(n.root);
		}
	}
}

