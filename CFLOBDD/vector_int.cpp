#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "vector_int_top_node.h"
#include "vector_int.h"
#include "matrix1234_int.h"

namespace CFL_OBDD {

	namespace VectorInt {

		void VectorInitializer()
		{
			VectorInitializerTop();
			return;
		}

		CFLOBDD MkBasisVector(unsigned int level, unsigned int index)
		{
			return CFLOBDD(MkBasisVectorTop(level, index));
		}

		CFLOBDD MkColumn1Matrix(unsigned int level)
		{
			return CFLOBDD(MkColumn1MatrixTop(level));
		}

		CFLOBDD VectorToMatrixInterleaved(CFLOBDD c)
		{
			CFLOBDD tempNode = CFLOBDD(VectorToMatrixInterleavedTop(c.root));
			CFLOBDD v = MkColumn1Matrix(tempNode.root->level);
			return tempNode * v;
		}

		CFLOBDD MkVectorWithVoc12(CFLOBDD c)
		{
			return CFLOBDD(MkVectorWithVoc12Top(c.root));
		}

		CFLOBDD KroneckerProduct(CFLOBDD m1, CFLOBDD m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD m2_1To2 = VectorShiftVocs1To2(m2);
			m2_1To2.PrintYield(&std::cout);
			std::cout << std::endl;
			std::cout << std::endl;
			return m1 * m2_1To2;
		}

		CFLOBDD VectorShiftVocs1To2(CFLOBDD m1)
		{
			return CFLOBDD(VectorShiftVocs1To2Top(m1.root));
		}

		CFLOBDD VectorPadWithZeros(CFLOBDD c, unsigned int level)
		{
			assert(c.root->level < level);
			CFLOBDD m(c);
			for (unsigned int i = c.root->level; i < level; i++)
			{
				CFLOBDD tempV = VectorInt::MkBasisVector(m.root->level, 0);
				tempV = VectorInt::MkVectorWithVoc12(tempV);
				CFLOBDD m_12 = VectorInt::MkVectorWithVoc12(m);
				m = VectorInt::KroneckerProduct(tempV, m_12);
			}

			assert(m.root->level == level);
			return m;
		}

		CFLOBDD NoDistinctionNode(unsigned int level)
		{
			return CFLOBDD(NoDistinctionNodeTop(level));
		}

		// CFLOBDD MatrixToVector(CFLOBDD c)
		// {
		// 	CFLOBDD temp = NoDistinctionNode(c.root->level);
		// 	temp = Matrix1234Int::PromoteInterleavedTo12(temp);
		// 	c = Matrix1234Int::PromoteInterleavedTo12(c);
		// 	CFLOBDD v = Matrix1234Int::MatrixMultiply(c,temp);
		// 	v = Matrix1234Int::Demote12ToInterleaved(v);
		// 	return CFLOBDD(MatrixToVectorTop(v.root));
		// }
	}
}

