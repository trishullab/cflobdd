#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "vector_float_boost_top_node.h"
#include "vector_float_boost.h"
#include "matrix1234_float_boost.h"

namespace CFL_OBDD {

	namespace VectorFloatBoost {

		void VectorInitializer()
		{
			VectorInitializerTop();
			return;
		}

		CFLOBDD_FLOAT_BOOST MkBasisVector(unsigned int level, unsigned int index)
		{
			return CFLOBDD_FLOAT_BOOST(MkBasisVectorTop(level, index));
		}

		CFLOBDD_FLOAT_BOOST MkBasisVector(unsigned int level, std::string s)
		{
			return CFLOBDD_FLOAT_BOOST(MkBasisVectorTop(level, s));
		}

		CFLOBDD_FLOAT_BOOST MkColumn1Matrix(unsigned int level)
		{
			return CFLOBDD_FLOAT_BOOST(MkColumn1MatrixTop(level));
		}

		CFLOBDD_FLOAT_BOOST VectorToMatrixInterleaved(CFLOBDD_FLOAT_BOOST c)
		{
			CFLOBDD_FLOAT_BOOST tempNode = CFLOBDD_FLOAT_BOOST(VectorToMatrixInterleavedTop(c.root));
			CFLOBDD_FLOAT_BOOST v = MkColumn1Matrix(tempNode.root->level);
			return tempNode * v;
		}

		CFLOBDD_FLOAT_BOOST MkVectorWithVoc12(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(MkVectorWithVoc12Top(c.root));
		}

		CFLOBDD_FLOAT_BOOST KroneckerProduct(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_FLOAT_BOOST m2_1To2 = VectorShiftVocs1To2(m2);
			return m1 * m2_1To2;
		}

		CFLOBDD_FLOAT_BOOST VectorShiftVocs1To2(CFLOBDD_FLOAT_BOOST m1)
		{
			return CFLOBDD_FLOAT_BOOST(VectorShiftVocs1To2Top(m1.root));
		}

		CFLOBDD_FLOAT_BOOST VectorPadWithZeros(CFLOBDD_FLOAT_BOOST c, unsigned int level)
		{
			assert(c.root->level < level);
			CFLOBDD_FLOAT_BOOST m(c);
			for (unsigned int i = c.root->level; i < level; i++)
			{
				CFLOBDD_FLOAT_BOOST tempV = VectorFloatBoost::MkBasisVector(m.root->level, 0);
				tempV = VectorFloatBoost::MkVectorWithVoc12(tempV);
				CFLOBDD_FLOAT_BOOST m_12 = VectorFloatBoost::MkVectorWithVoc12(m);
				m = VectorFloatBoost::KroneckerProduct(tempV, m_12);
			}

			assert(m.root->level == level);
			return m;
		}

		CFLOBDD_FLOAT_BOOST NoDistinctionNode(unsigned int level, int val)
		{
			return CFLOBDD_FLOAT_BOOST(NoDistinctionNodeTop(level, val));
		}

		CFLOBDD_FLOAT_BOOST ConvertToDouble(CFLOBDD_FLOAT_BOOST n)
		{
			return CFLOBDD_FLOAT_BOOST(ConvertToDoubleTop(n.root));
		}

		CFLOBDD_FLOAT_BOOST MatrixToVector(CFLOBDD_FLOAT_BOOST c)
		{
			CFLOBDD_FLOAT_BOOST temp = NoDistinctionNode(c.root->level, 1);
			CFLOBDD_FLOAT_BOOST v = Matrix1234FloatBoost::MatrixMultiplyV4(c, temp);
			return CFLOBDD_FLOAT_BOOST(MatrixToVectorTop(v.root));
		}

//#ifdef PATH_COUNTING_ENABLED
		std::string Sampling(CFLOBDD_FLOAT_BOOST c, bool isTwoVoc, std::string func)
		{
			return SamplingTop(c.root, isTwoVoc, func);
		}

		std::string SamplingV2(CFLOBDD_FLOAT_BOOST c)
		{
			return SamplingV2Top(c.root);
		}
//#endif

		CFLOBDD_FLOAT_BOOST VectorWithAmplitude(CFLOBDD_FLOAT_BOOST c)
		{
			return VectorWithAmplitudeTop(c.root);
		}

		void VectorPrintColumnMajor(CFLOBDD_FLOAT_BOOST c, std::ostream & out)
		{
			VectorPrintColumnMajorTop(c.root, out);
			return;
		}

		void VectorPrintColumnMajorInterleaved(CFLOBDD_FLOAT_BOOST c, std::ostream & out)
		{
			VectorPrintColumnMajorInterleavedTop(c.root, out);
			return;
		}
	}
}

