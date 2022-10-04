#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "vector_double_top_node.h"
#include "vector_double.h"
#include "matrix1234_double.h"

namespace CFL_OBDD {

	namespace VectorDouble {

		void VectorInitializer()
		{
			VectorInitializerTop();
			return;
		}

		CFLOBDD_DOUBLE MkBasisVector(unsigned int level, unsigned int index)
		{
			return CFLOBDD_DOUBLE(MkBasisVectorTop(level, index));
		}

		CFLOBDD_DOUBLE MkColumn1Matrix(unsigned int level)
		{
			return CFLOBDD_DOUBLE(MkColumn1MatrixTop(level));
		}

		CFLOBDD_DOUBLE VectorToMatrixInterleaved(CFLOBDD_DOUBLE c)
		{
			CFLOBDD_DOUBLE tempNode = CFLOBDD_DOUBLE(VectorToMatrixInterleavedTop(c.root));
			CFLOBDD_DOUBLE v = MkColumn1Matrix(tempNode.root->level);
			return tempNode * v;
		}

		CFLOBDD_DOUBLE MkVectorWithVoc12(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(MkVectorWithVoc12Top(c.root));
		}

		CFLOBDD_DOUBLE KroneckerProduct(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_DOUBLE m2_1To2 = VectorShiftVocs1To2(m2);
			return m1 * m2_1To2;
		}

		CFLOBDD_DOUBLE VectorShiftVocs1To2(CFLOBDD_DOUBLE m1)
		{
			return CFLOBDD_DOUBLE(VectorShiftVocs1To2Top(m1.root));
		}

		CFLOBDD_DOUBLE VectorPadWithZeros(CFLOBDD_DOUBLE c, unsigned int level)
		{
			assert(c.root->level < level);
			CFLOBDD_DOUBLE m(c);
			for (unsigned int i = c.root->level; i < level; i++)
			{
				CFLOBDD_DOUBLE tempV = VectorDouble::MkBasisVector(m.root->level, 0);
				tempV = VectorDouble::MkVectorWithVoc12(tempV);
				CFLOBDD_DOUBLE m_12 = VectorDouble::MkVectorWithVoc12(m);
				m = VectorDouble::KroneckerProduct(tempV, m_12);
			}

			assert(m.root->level == level);
			return m;
		}

		CFLOBDD_DOUBLE NoDistinctionNode(unsigned int level)
		{
			return CFLOBDD_DOUBLE(NoDistinctionNodeTop(level));
		}

		CFLOBDD_DOUBLE MatrixToVector(CFLOBDD_DOUBLE c)
		{
			CFLOBDD_DOUBLE temp = NoDistinctionNode(c.root->level);
			CFLOBDD_DOUBLE v = Matrix1234Double::MatrixMultiplyV4(c, temp);
			return CFLOBDD_DOUBLE(MatrixToVectorTop(v.root));
		}

#ifdef PATH_COUNTING_ENABLED
		std::string Sampling(CFLOBDD_DOUBLE c)
		{
			return SamplingTop(c.root);
		}
#endif

		CFLOBDD_DOUBLE VectorWithAmplitude(CFLOBDD_DOUBLE c)
		{
			return VectorWithAmplitudeTop(c.root);
		}

		void VectorPrintColumnMajor(CFLOBDD_DOUBLE c, std::ostream & out)
		{
			VectorPrintColumnMajorTop(c.root, out);
			return;
		}

		void VectorPrintColumnMajorInterleaved(CFLOBDD_DOUBLE c, std::ostream & out)
		{
			VectorPrintColumnMajorInterleavedTop(c.root, out);
			return;
		}
	}
}

