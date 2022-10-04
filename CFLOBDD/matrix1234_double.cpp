#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "matrix1234_node.h"
#include "matrix1234_double_top_node.h"
#include "matrix1234_double.h"
#include "vector_double.h"
#include "matrix1234_int.h"

namespace CFL_OBDD {

	namespace Matrix1234Double {

		void Matrix1234Initializer()
		{
			Matrix1234InitializerTop();
			return;
		}

		// Create representation of identity relation
		CFLOBDD_DOUBLE MkIdRelationInterleaved(unsigned int i)
		{
			return CFLOBDD_DOUBLE(MkIdRelationInterleavedTop(i));
		}

		CFLOBDD_DOUBLE MkNegationMatrixInterleaved(unsigned int i)
		{
			return CFLOBDD_DOUBLE(MkNegationMatrixInterleavedTop(i));
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDD_DOUBLE MkWalshInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD::maxLevel);
			return CFLOBDD_DOUBLE(MkWalshInterleavedTop(i));
		}

		CFLOBDD_DOUBLE MkCNOTInterleaved(unsigned int i)
		{
			return CFLOBDD_DOUBLE(MkCNOTInterleavedTop(i));
		}


		// ****************************************************************************
		// Matrix-related operations (on matrices with room for two extra vocabularies)
		// ****************************************************************************

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for
		// two extra vocabularies
		CFLOBDD_DOUBLE MkWalshVoc12(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD_DOUBLE(MkWalshVoc12Top(i));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_DOUBLE MatrixShiftVocs12To34(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(MatrixShiftVocs12To34Top(c.root));
		}

		CFLOBDD_DOUBLE PromoteInterleavedTo12(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(PromoteInterleavedTo12Top(c.root));
		}

		extern CFLOBDD_DOUBLE Demote12ToInterleaved(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(Demote12ToInterleavedTop(c.root));
		}

		CFLOBDD_DOUBLE ReverseColumns(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(ReverseColumnsTop(c.root));
		}

		CFLOBDD_DOUBLE MatrixTranspose(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(MatrixTransposeTop(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_DOUBLE MatrixShiftVoc42(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(MatrixShiftVoc42Top(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_DOUBLE KroneckerProduct(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_DOUBLE m2_12To34 = MatrixShiftVocs12To34(m2);
			return m1 * m2_12To34;
		}

		CFLOBDD_DOUBLE MatrixShiftToAConnection(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(MatrixShiftToAConnectionTop(c.root));
		}

		CFLOBDD_DOUBLE MatrixShiftToBConnection(CFLOBDD_DOUBLE c)
		{
			return CFLOBDD_DOUBLE(MatrixShiftToBConnectionTop(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_DOUBLE KroneckerProduct2Vocs(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_DOUBLE m1_A = MatrixShiftToAConnection(m1);
			CFLOBDD_DOUBLE m2_B = MatrixShiftToBConnection(m2);
			return m1_A * m2_B;
		}

		CFLOBDD_DOUBLE ConvertIntToDouble(CFLOBDD c)
		{
			return CFLOBDD_DOUBLE(ConvertIntToDoubleTop(c.root));
		}

		//
		// MatrixDetensor
		//
		// Return the result of detensoring a 4-vocabulary matrix
		//
		// CFLOBDD_DOUBLE MatrixDetensor(CFLOBDD_DOUBLE k)
		// {
		// 	assert(k.root->level >= 2);
		// 	//std::cout << "Kronecker product done.." << std::endl;
		// 	CFLOBDD_DOUBLE e = MkDetensorConstraintInterleaved(k.root->level);
		// 	//std::cout << "Detensor Constraint done.." << std::endl;
		// 	CFLOBDD_DOUBLE m = k * e;
		// 	//std::cout << "Project start.. " << m.root->level << std::endl;
		// 	CFLOBDD_DOUBLE p = MatrixProjectVoc23(m);
		// 	//std::cout << "Project done.." << std::endl;
		// 	CFLOBDD_DOUBLE ans = MatrixShiftVoc42(p);
		// 	//std::cout << "Multiply ans.." << std::endl;
		// 	return ans;
		// }

		//
		// MatrixMultiply
		//
		// Return the matrix product of m1 and m2
		//
		// CFLOBDD_DOUBLE MatrixMultiply(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2)
		// {
		// 	assert(m1.root->level == m2.root->level);
		// 	assert(m1.root->level >= 2);

		// 	CFLOBDD_DOUBLE k = KroneckerProduct(m1, m2);
		// 	CFLOBDD_DOUBLE ans = MatrixDetensor(k);
		// 	return ans;
		// }

		//
		// MatrixMultiplyV4
		// Naive Matrix Multiplication
		//
		// Return the matrix product of m1 and m2
		//
		CFLOBDD_DOUBLE MatrixMultiplyV4(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2)
		{
			return CFLOBDD_DOUBLE(MatrixMultiplyV4TopNode(m1.root, m2.root));
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal.
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDD_DOUBLE MkDetensorConstraintInterleaved(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD_DOUBLE(MkDetensorConstraintInterleavedTop(i));
		}

		// // Vocabulary projection
		// CFLOBDD_DOUBLE MatrixProjectVoc23(CFLOBDD_DOUBLE c)
		// {
		// 	return CFLOBDD_DOUBLE(MatrixProjectVoc23Top(c.root));
		// }

		// Convert Matrix to Matrix with additional vocabularies
		/*CFLOBDD_DOUBLE MatrixConvertVocs1234To12(CFLOBDD_DOUBLE c)
		{
		return CFLOBDD_DOUBLE(MatrixConvertVocs1234To12Top(c.root));
		}*/

		CFLOBDD_DOUBLE MatrixPadWithZeros(CFLOBDD_DOUBLE c, unsigned int level)
		{
			assert(c.root->level < level);
			CFLOBDD_DOUBLE m(c);
			for (unsigned int i = c.root->level; i < level; i++)
			{
				CFLOBDD_DOUBLE tempV = VectorDouble::MkBasisVector(m.root->level - 1, 0);
				CFLOBDD_DOUBLE tempM = VectorDouble::VectorToMatrixInterleaved(tempV);
				m = Matrix1234Double::KroneckerProduct2Vocs(tempM, m);
			}

			assert(m.root->level == level);
			return m;
		}

		void MatrixPrintRowMajor(CFLOBDD_DOUBLE c, std::ostream & out)
		{
			MatrixPrintRowMajorTop(c.root, out);
			return;
		}

		void MatrixPrintRowMajorInterleaved(CFLOBDD_DOUBLE c, std::ostream & out)
		{
			MatrixPrintRowMajorInterleavedTop(c.root, out);
			return;
		}

		CFLOBDD_DOUBLE SMatrix(std::string s)
		{
			return CFLOBDD_DOUBLE(SMatrixTop(s));
		}

		// CFLOBDD_DOUBLE PermuteColumns(CFLOBDD_DOUBLE c, std::string m)
		// {
		// 	if (m.at(m.length() - 1) == '0')
		// 	{
		// 		CFLOBDD p = Matrix1234Int::MkCFLOBDDMatrixEqVoc14(c.root->level) * Matrix1234Int::MkDetensorConstraintInterleaved(c.root->level);
		// 		CFLOBDD_DOUBLE pd = Matrix1234Double::ConvertIntToDouble(p);
		// 		pd = Matrix1234Double::MatrixTranspose(pd);
		// 		pd = PromoteInterleavedTo12(pd);
		// 		c = PromoteInterleavedTo12(c);
		// 		CFLOBDD_DOUBLE ans = Matrix1234Double::MatrixMultiply(c, pd);
		// 		ans = Demote12ToInterleaved(ans);
		// 		return ans;
		// 	}
		// 	else
		// 	{
		// 		CFLOBDD_DOUBLE p = Matrix1234Double::MkCNOTInterleaved(c.root->level);
		// 		p = Matrix1234Double::MatrixTranspose(p);
		// 		p = PromoteInterleavedTo12(p);
		// 		c = PromoteInterleavedTo12(c);
		// 		CFLOBDD_DOUBLE ans = Matrix1234Double::MatrixMultiply(c, p);
		// 		ans = Demote12ToInterleaved(ans);
		// 		return ans;
		// 	}
		// }

		// CFLOBDD_DOUBLE Create2To1Func(CFLOBDD_DOUBLE F1, std::string s1, std::string m1, CFLOBDD_DOUBLE F2, std::string s2, std::string m2)
		// {
		// 	CFLOBDD_DOUBLE S1 = SMatrix(s1);
		// 	CFLOBDD_DOUBLE S2 = SMatrix(s2);
		// 	CFLOBDD_DOUBLE F1prime = F1 * S1;
		// 	CFLOBDD_DOUBLE F2prime = F2 * S2;
		// 	CFLOBDD_DOUBLE F1Permute = F1prime;
		// 	CFLOBDD_DOUBLE F2Permute = F2prime;
		// 	if (!(s1.find('1') == std::string::npos || s2.find('1') == std::string::npos))
		// 	{
		// 		F1Permute = PermuteColumns(F1prime, m1);
		// 		F2Permute = PermuteColumns(F2prime, m2);
		// 	}

		// 	CFLOBDD_DOUBLE t1 = Matrix1234Double::KroneckerProduct2Vocs(F1prime, F2prime);
		// 	//Matrix1234Double::MatrixPrintRowMajor(t1, std::cout);
			
		// 	DoubleReturnMapHandle rt1 = t1.root->rootConnection.returnMapHandle;
		// 	DoubleReturnMapHandle rt1temp;
		// 	for (unsigned int i = 0; i < rt1.Size(); i++)
		// 	{
		// 		if (rt1[i] == -1)
		// 			rt1temp.AddToEnd(0);
		// 		else
		// 			rt1temp.AddToEnd(rt1[i]);
		// 	}
		// 	rt1temp.Canonicalize();
		// 	CFLOBDDNodeHandle nt1 = *(t1.root->rootConnection.entryPointHandle);
		// 	ReductionMapHandle inducedReductionMapHandle;
		// 	DoubleReturnMapHandle inducedReturnMap;
		// 	rt1temp.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
		// 	CFLOBDDNodeHandle reduced_n = nt1.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
		// 	t1 = CFLOBDD_DOUBLE(new CFLOBDDTopNodeDouble(reduced_n, inducedReturnMap));

		// 	CFLOBDD_DOUBLE t2 = Matrix1234Double::KroneckerProduct2Vocs(F1Permute, F2Permute);
		// 	ReductionMapHandle inducedReductionMapHandle_2;
		// 	DoubleReturnMapHandle inducedReturnMap_2;
		// 	DoubleReturnMapHandle rt2 = t2.root->rootConnection.returnMapHandle;
		// 	DoubleReturnMapHandle rt2temp;
		// 	for (unsigned int i = 0; i < rt2.Size(); i++)
		// 	{
		// 		if (rt2[i] == 1)
		// 			rt2temp.AddToEnd(0);
		// 		else if (rt2[i] == -1)
		// 			rt2temp.AddToEnd(1);
		// 		else
		// 			rt2temp.AddToEnd(rt2[i]);
		// 	}
		// 	rt2temp.Canonicalize();
		// 	CFLOBDDNodeHandle nt2 = *(t2.root->rootConnection.entryPointHandle);
		// 	rt2temp.InducedReductionAndReturnMap(inducedReductionMapHandle_2, inducedReturnMap_2);
		// 	CFLOBDDNodeHandle reduced_n_2 = nt2.Reduce(inducedReductionMapHandle_2, inducedReturnMap_2.Size());
		// 	t2 = CFLOBDD_DOUBLE(new CFLOBDDTopNodeDouble(reduced_n_2, inducedReturnMap_2));
		// 	CFLOBDD_DOUBLE ans = t1 + t2;
		// 	return ans;
		// }

		// CFLOBDD_DOUBLE Func2To1CFLOBDDMatrix(std::string s, std::string t)
		// {
		// 	if (s.length() == 2)
		// 	{
		// 		if (s == "00")
		// 			return Matrix1234Double::MkIdRelationInterleaved(2);
		// 		else if (s == "01")
		// 		{
		// 			CFLOBDD I = Matrix1234Int::MkIdRelationInterleaved(2);
		// 			I = MkRestrict(I, 2, false);
		// 			CFLOBDD_DOUBLE Id = Matrix1234Double::ConvertIntToDouble(I);
		// 			return Id;
		// 		}
		// 		else if (s == "10")
		// 		{
		// 			CFLOBDD I = Matrix1234Int::MkIdRelationInterleaved(2);
		// 			I = MkRestrict(I, 0, false);
		// 			CFLOBDD_DOUBLE Id = Matrix1234Double::ConvertIntToDouble(I);
		// 			return Id;
		// 		}
		// 		else if (s == "11")
		// 		{
		// 			CFLOBDD_DOUBLE A = Matrix1234Double::MkIdRelationInterleaved(1);
		// 			CFLOBDD_DOUBLE B = Matrix1234Double::MkNegationMatrixInterleaved(1);
					
		// 			CFLOBDD_DOUBLE v1 = VectorDouble::MkBasisVector(0, 0);
		// 			v1 = VectorDouble::VectorToMatrixInterleaved(v1);
		// 			CFLOBDD_DOUBLE answer = Matrix1234Double::KroneckerProduct2Vocs(v1, A);

		// 			CFLOBDD_DOUBLE v2 = VectorDouble::MkBasisVector(0, 1);
		// 			v2 = VectorDouble::VectorToMatrixInterleaved(v2);
		// 			answer = answer + Matrix1234Double::KroneckerProduct2Vocs(v2, B);
		// 			return answer;
		// 		}
		// 	}
		// 	CFLOBDD_DOUBLE F1 = Func2To1CFLOBDDMatrix(s.substr(0, s.length() / 2), t.substr(0, t.length()/2));
		// 	CFLOBDD_DOUBLE F2 = Func2To1CFLOBDDMatrix(s.substr(s.length() / 2), t.substr(t.length()/2));
		// 	CFLOBDD_DOUBLE ans = Create2To1Func(F1, s.substr(0, s.length() / 2), t.substr(0,t.length()/2), F2, s.substr(s.length() / 2), t.substr(t.length()/2));
		// 	return ans;
		// }

		CFLOBDD_DOUBLE NormalizeOutputTo1(CFLOBDD_DOUBLE c)
		{
			ReductionMapHandle inducedReductionMapHandle_2;
			DoubleReturnMapHandle inducedReturnMap_2;
			DoubleReturnMapHandle rt = c.root->rootConnection.returnMapHandle;
			DoubleReturnMapHandle rttemp;
			for (unsigned int i = 0; i < rt.Size(); i++)
			{
				if (rt[i] != 0)
					rttemp.AddToEnd(1);
				else
					rttemp.AddToEnd(rt[i]);
			}
			rttemp.Canonicalize();
			CFLOBDDNodeHandle nt = *(c.root->rootConnection.entryPointHandle);
			rttemp.InducedReductionAndReturnMap(inducedReductionMapHandle_2, inducedReturnMap_2);
			CFLOBDDNodeHandle reduced_n_2 = nt.Reduce(inducedReductionMapHandle_2, inducedReturnMap_2.Size());
			return CFLOBDD_DOUBLE(new CFLOBDDTopNodeDouble(reduced_n_2, inducedReturnMap_2));
		}
	}
}

