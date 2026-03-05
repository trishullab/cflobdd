#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <cstdarg>
#include <boost/random.hpp>
#include <boost/preprocessor/logical/xor.hpp>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "matrix1234_node.h"
#include "matrix1234_float_boost_top_node.h"
#include "matrix1234_float_boost.h"
#include "vector_float_boost.h"
#include "matrix1234_int.h"

namespace CFL_OBDD {

	typedef boost::random::independent_bits_engine<boost::random::mt19937, 256, boost::multiprecision::cpp_int> generator_type;
	generator_type gen;

	namespace Matrix1234FloatBoost {

		void Matrix1234Initializer()
		{
			Matrix1234InitializerTop();
			return;
		}

		// Create representation of identity relation
		CFLOBDD_FLOAT_BOOST MkIdRelationInterleaved(unsigned int i)
		{
			// TODO: Check for error - "CodeConvert"
			CFLOBDDTopNodeFloatBoostRefPtr tmp = MkIdRelationInterleavedTop(i);
			return CFLOBDD_FLOAT_BOOST(tmp);
		}

		CFLOBDD_FLOAT_BOOST MkNegationMatrixInterleaved(unsigned int i)
		{
			return CFLOBDD_FLOAT_BOOST(MkNegationMatrixInterleavedTop(i));
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDD_FLOAT_BOOST MkWalshInterleaved(unsigned int i)
		{
			assert(i <= CFLOBDD::maxLevel);
			return CFLOBDD_FLOAT_BOOST(MkWalshInterleavedTop(i));
		}

		CFLOBDD_FLOAT_BOOST MkCNOTInterleaved(unsigned int i)
		{
			return CFLOBDD_FLOAT_BOOST(MkCNOTInterleavedTop(i));
		}

		CFLOBDD_FLOAT_BOOST MkExchangeInterleaved(unsigned int i)
		{
			return CFLOBDD_FLOAT_BOOST(MkExchangeInterleavedTop(i));
		}


		// ****************************************************************************
		// Matrix-related operations (on matrices with room for two extra vocabularies)
		// ****************************************************************************

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for
		// two extra vocabularies
		CFLOBDD_FLOAT_BOOST MkWalshVoc12(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD_FLOAT_BOOST(MkWalshVoc12Top(i));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_FLOAT_BOOST MatrixShiftVocs12To34(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(MatrixShiftVocs12To34Top(c.root));
		}

		CFLOBDD_FLOAT_BOOST PromoteInterleavedTo12(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(PromoteInterleavedTo12Top(c.root));
		}

		CFLOBDD_FLOAT_BOOST PromoteInterleavedTo13(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(PromoteInterleavedTo13Top(c.root));
		}

		extern CFLOBDD_FLOAT_BOOST Demote12ToInterleaved(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(Demote12ToInterleavedTop(c.root));
		}

		CFLOBDD_FLOAT_BOOST ReverseColumns(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(ReverseColumnsTop(c.root));
		}

		CFLOBDD_FLOAT_BOOST MatrixTranspose(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(MatrixTransposeTop(c.root));
		}

		// Vocabulary shift in a matrix
		CFLOBDD_FLOAT_BOOST MatrixShiftVoc42(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(MatrixShiftVoc42Top(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_FLOAT_BOOST KroneckerProduct(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_FLOAT_BOOST m2_12To34 = MatrixShiftVocs12To34(m2);
			//std::cout << "voc shift done" << std::endl;
			return m1 * m2_12To34;
		}

		CFLOBDD_FLOAT_BOOST MatrixShiftToAConnection(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(MatrixShiftToAConnectionTop(c.root));
		}

		CFLOBDD_FLOAT_BOOST MatrixShiftToBConnection(CFLOBDD_FLOAT_BOOST c)
		{
			return CFLOBDD_FLOAT_BOOST(MatrixShiftToBConnectionTop(c.root));
		}

		// Return the Kronecker product of two matrices
		CFLOBDD_FLOAT_BOOST KroneckerProduct2Vocs(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2)
		{
			assert(m1.root->level == m2.root->level);
			CFLOBDD_FLOAT_BOOST m1_A = MatrixShiftToAConnection(m1);
			CFLOBDD_FLOAT_BOOST m2_B = MatrixShiftToBConnection(m2);
			CFLOBDD_FLOAT_BOOST c = m1_A * m2_B;
			return c;
		}

		CFLOBDD_FLOAT_BOOST ConvertIntToFloatBoost(CFLOBDD c)
		{
			return CFLOBDD_FLOAT_BOOST(ConvertIntToFloatBoostTop(c.root));
		}

		//
		// MatrixDetensor
		//
		// Return the result of detensoring a 4-vocabulary matrix
		//
		// CFLOBDD_FLOAT_BOOST MatrixDetensor(CFLOBDD_FLOAT_BOOST k)
		// {
		// 	assert(k.root->level >= 2);
		// 	//std::cout << "Kronecker product done.." << std::endl;
		// 	CFLOBDD_FLOAT_BOOST e = MkDetensorConstraintInterleaved(k.root->level);
		// 	//std::cout << "Detensor Constraint done.." << std::endl;
		// 	CFLOBDD_FLOAT_BOOST m = k * e;
		// 	DisposeOfPairProductCache();
		// 	InitPairProductCache();
		// 	CFLOBDDNodeHandle::DisposeOfReduceCache();
		// 	CFLOBDDNodeHandle::InitReduceCache();
		// 	unsigned int nodeCount = 0;
		// 	m.CountNodes(nodeCount);
		// 	//std::cout << "Project start.. " << m.root->level << " " << nodeCount << " " << m.root->rootConnection.returnMapHandle.Size() << std::endl;
		// 	//std::cout << "IsValid: " << m.IsValid() << std::endl;
		// 	CFLOBDD_FLOAT_BOOST p = Matrix1234FloatBoost::MkIdRelationInterleaved(1);
		// 	//std::cout << "Project begin.. " << std::endl;
		// 	try{
		// 		p = MatrixProjectVoc23(m);
		// 		//printMemory();
		// 	}
		// 	catch (std::exception e){
		// 		std::cout << e.what() << std::endl;
		// 		std::cout << "Matrix Detensor" << std::endl;
		// 		// printMemory();
		// 		throw e;
		// 	}
		// 	//std::cout << "Project done.." << std::endl;
		// 	CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MkIdRelationInterleaved(1);
		// 	try{
		// 		ans = MatrixShiftVoc42(p);
		// 	}
		// 	catch (std::exception e){
		// 		std::cout << e.what() << std::endl;
		// 		std::cout << "Matrix ShiftVoc42" << std::endl;
		// 		// printMemory();
		// 		throw e;
		// 	}
		// 	//std::cout << "Multiply ans.." << std::endl;
		// 	return ans;
		// }


		// //
		// // MatrixMultiply
		// //
		// // Return the matrix product of m1 and m2
		// //
		// CFLOBDD_FLOAT_BOOST MatrixMultiply(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2)
		// {
		// 	assert(m1.root->level == m2.root->level);
		// 	assert(m1.root->level >= 2);

		// 	CFLOBDD_FLOAT_BOOST k = KroneckerProduct(m1, m2);
		// 	//std::cout << "Kronecker Product done" << std::endl;
		// 	CFLOBDD_FLOAT_BOOST ans = MatrixDetensor(k);
		// 	return ans;
		// }

		CFLOBDD_FLOAT_BOOST MatrixMultiplyV3(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2)
		{
			assert(m1.root->level == m2.root->level);
			assert(m1.root->level >= 1);
			CFLOBDD_FLOAT_BOOST ans = KroneckerProduct2Vocs(m1, m2);
			Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(ans, std::cout);
			CFLOBDD_FLOAT_BOOST k = MkMatrixMultiplyConstraint(ans.root->level);
			Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(k, std::cout);
			ans = ans * k;
			Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(ans, std::cout);
			std::cout << ans << std::endl;
			ans = MultiplyOperation(ans);
			return ans;
		}

		// Naive matrix multiplication
		CFLOBDD_FLOAT_BOOST MatrixMultiplyV4(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2)
		{
			return CFLOBDD_FLOAT_BOOST(MatrixMultiplyV4TopNode(m1.root, m2.root));
		}

		// Naive matrix multiplication with 0 nodes information
		CFLOBDD_FLOAT_BOOST MatrixMultiplyV4WithInfo(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2)
		{
			return CFLOBDD_FLOAT_BOOST(MatrixMultiplyV4WithInfoTopNode(m1.root, m2.root));
		}

		CFLOBDD_FLOAT_BOOST MkMatrixMultiplyConstraint(unsigned int level)
		{
			return CFLOBDD_FLOAT_BOOST(MkMatrixMultiplyConstraintTopNode(level));
		}

		CFLOBDD_FLOAT_BOOST MultiplyOperation(CFLOBDD_FLOAT_BOOST m)
		{
			return CFLOBDD_FLOAT_BOOST(MultiplyOperationTopNode(m.root));
		}

		CFLOBDD_FLOAT_BOOST MkDistinctionTwoVars(int x, int y, unsigned int var_level, unsigned int total_level)
		{
			return CFLOBDD_FLOAT_BOOST(MkDistinctionTwoVarsTopNode(x, y, var_level, total_level));
		}

		CFLOBDD_FLOAT_BOOST SubMatrix(CFLOBDD_FLOAT_BOOST m, int x, int y, unsigned int matrix_level)
		{
			return m * MkDistinctionTwoVars(x, y, ((unsigned int)1 << matrix_level), m.root->level);
		}

		//
		// MatrixMultiply V2
		//
		// Return the matrix product of m1 and m2
		//
		CFLOBDD_FLOAT_BOOST MatrixMultiplyV2(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2, int matrix_level = 0)
		{
			assert(m1.root->level == m2.root->level);
			assert(m1.root->level >= 2);

			CFLOBDD_FLOAT_BOOST m10 = SubMatrix(m1, 0, 0, matrix_level);
			CFLOBDD_FLOAT_BOOST m11 = SubMatrix(m1, 0, 1, matrix_level);
			CFLOBDD_FLOAT_BOOST m12 = SubMatrix(m1, 1, 0, matrix_level);
			CFLOBDD_FLOAT_BOOST m13 = SubMatrix(m1, 1, 1, matrix_level);

			CFLOBDD_FLOAT_BOOST m20 = SubMatrix(m2, 0, 0, matrix_level);
			CFLOBDD_FLOAT_BOOST m21 = SubMatrix(m2, 0, 1, matrix_level);
			CFLOBDD_FLOAT_BOOST m22 = SubMatrix(m2, 1, 0, matrix_level);
			CFLOBDD_FLOAT_BOOST m23 = SubMatrix(m2, 1, 1, matrix_level);

			// Strassen's multiplication
			CFLOBDD_FLOAT_BOOST M1 = MatrixMultiplyV2(m10 + m13, m20 + m23, matrix_level + 1);
			CFLOBDD_FLOAT_BOOST M2 = MatrixMultiplyV2(m12 + m13, m20, matrix_level + 1);
			CFLOBDD_FLOAT_BOOST M3 = MatrixMultiplyV2(m10, m21 + (-1 * m23), matrix_level + 1);
			CFLOBDD_FLOAT_BOOST M4 = MatrixMultiplyV2(m13, m22 + (-1 * m20), matrix_level + 1);
			CFLOBDD_FLOAT_BOOST M5 = MatrixMultiplyV2(m10 + m11, m23, matrix_level + 1);
			CFLOBDD_FLOAT_BOOST M6 = MatrixMultiplyV2(m12 + (-1 * m10), m20 + m21, matrix_level + 1);
			CFLOBDD_FLOAT_BOOST M7 = MatrixMultiplyV2(m11 + (-1 * m13), m22 + m23, matrix_level + 1);

			return m1;
		}
		
		// CFLOBDD_FLOAT_BOOST AddMatrixRows(CFLOBDD_FLOAT_BOOST m){
		// 	assert(1 <= m.root->level && m.root->level <= CFLOBDD::maxLevel);
		// 	return CFLOBDD_FLOAT_BOOST(AddMatrixRowsTopNode(m.root));
		// }
		
		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal.
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDD_FLOAT_BOOST MkDetensorConstraintInterleaved(unsigned int i)
		{
			assert(2 <= i && i <= CFLOBDD::maxLevel);
			return CFLOBDD_FLOAT_BOOST(MkDetensorConstraintInterleavedTop(i));
		}

		// Vocabulary projection
		// CFLOBDD_FLOAT_BOOST MatrixProjectVoc23(CFLOBDD_FLOAT_BOOST c)
		// {
		// 	return CFLOBDD_FLOAT_BOOST(MatrixProjectVoc23Top(c.root));
		// }

		// Convert Matrix to Matrix with additional vocabularies
		/*CFLOBDD_FLOAT_BOOST MatrixConvertVocs1234To12(CFLOBDD_FLOAT_BOOST c)
		{
		return CFLOBDD_FLOAT_BOOST(MatrixConvertVocs1234To12Top(c.root));
		}*/

		CFLOBDD_FLOAT_BOOST MatrixPadWithZeros(CFLOBDD_FLOAT_BOOST c, unsigned int level)
		{
			assert(c.root->level < level);
			CFLOBDD_FLOAT_BOOST m(c);
			for (unsigned int i = c.root->level; i < level; i++)
			{
				CFLOBDD_FLOAT_BOOST tempV = VectorFloatBoost::MkBasisVector(m.root->level - 1, 0);
				CFLOBDD_FLOAT_BOOST tempM = VectorFloatBoost::VectorToMatrixInterleaved(tempV);
				m = Matrix1234FloatBoost::KroneckerProduct2Vocs(tempM, m);
			}

			assert(m.root->level == level);
			return m;
		}

		void MatrixPrintRowMajor(CFLOBDD_FLOAT_BOOST c, std::ostream & out)
		{
			MatrixPrintRowMajorTop(c.root, out);
			return;
		}

		void MatrixPrintRowMajorInterleaved(CFLOBDD_FLOAT_BOOST c, std::ostream & out)
		{
			MatrixPrintRowMajorInterleavedTop(c.root, out);
			return;
		}

		CFLOBDD_FLOAT_BOOST SMatrix(std::string s)
		{
			return CFLOBDD_FLOAT_BOOST(SMatrixTop(s));
		}

		CFLOBDD_FLOAT_BOOST PermuteColumns(CFLOBDD_FLOAT_BOOST c, std::mt19937 mt)
		{
			
			if (c.root->level <= 4){
				int level = c.root->level;
				std::vector<int> v(pow(2, pow(2, level - 1)), 0);
				for (unsigned int i = 0; i < v.size(); i++)
					v[i] = i;
				std::shuffle(v.begin(), v.end(), mt);
				std::string s(pow(2, level-1), '0');
				CFLOBDD_FLOAT_BOOST p = Matrix1234FloatBoost::Func2To1CFLOBDDMatrix(s, v);
				//p = Matrix1234FloatBoost::MatrixTranspose(p);
				CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(c, p);
				//CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::Demote12ToInterleaved(Matrix1234FloatBoost::MatrixMultiply(
					//Matrix1234FloatBoost::PromoteInterleavedTo12(c), Matrix1234FloatBoost::PromoteInterleavedTo12(p)));
				return ans;
			}
			
			int rand_val = mt() % 3;
			
			if (rand_val == 0)
			{
				CFLOBDD p = Matrix1234Int::MkCFLOBDDMatrixEqVoc14(c.root->level) * Matrix1234Int::MkDetensorConstraintInterleaved(c.root->level);
				CFLOBDD_FLOAT_BOOST pd = Matrix1234FloatBoost::ConvertIntToFloatBoost(p);
				//pd = Matrix1234FloatBoost::MatrixTranspose(pd);
				CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(c, pd);
				//CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::Demote12ToInterleaved(Matrix1234FloatBoost::MatrixMultiply(
					//Matrix1234FloatBoost::PromoteInterleavedTo12(c), Matrix1234FloatBoost::PromoteInterleavedTo12(pd)));
				return ans;
			}
			else if (rand_val == 1)
			{
				CFLOBDD_FLOAT_BOOST p = Matrix1234FloatBoost::MkCNOTInterleaved(c.root->level);
				//p = Matrix1234FloatBoost::MatrixTranspose(p);
				CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(c, p);
				//CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::Demote12ToInterleaved(Matrix1234FloatBoost::MatrixMultiply(
					//Matrix1234FloatBoost::PromoteInterleavedTo12(c), Matrix1234FloatBoost::PromoteInterleavedTo12(p)));
				return ans;
			}
			else
			{
				CFLOBDD_FLOAT_BOOST p = Matrix1234FloatBoost::MkExchangeInterleaved(c.root->level);
				//p = Matrix1234FloatBoost::MatrixTranspose(p);
				CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(c, p);
				//CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::Demote12ToInterleaved(Matrix1234FloatBoost::MatrixMultiply(
					//Matrix1234FloatBoost::PromoteInterleavedTo12(c), Matrix1234FloatBoost::PromoteInterleavedTo12(p)));
				return ans;
			}
		}

		CFLOBDD_FLOAT_BOOST Create2To1Func(CFLOBDD_FLOAT_BOOST F1, std::string s1, CFLOBDD_FLOAT_BOOST F2, std::string s2, std::mt19937 mt)
		{
			CFLOBDD_FLOAT_BOOST S1 = SMatrix(s1);
			CFLOBDD_FLOAT_BOOST S2 = SMatrix(s2);
			CFLOBDD_FLOAT_BOOST F1prime = F1 * S1;
			CFLOBDD_FLOAT_BOOST F2prime = F2 * S2;
			CFLOBDD_FLOAT_BOOST F1Permute = F1prime;
			CFLOBDD_FLOAT_BOOST F2Permute = F2prime;
			if (!(s1.find('1') == std::string::npos || s2.find('1') == std::string::npos))
			{
				F1Permute = PermuteColumns(F1prime, mt);
				F2Permute = PermuteColumns(F2prime, mt);
			}

			CFLOBDD_FLOAT_BOOST t1 = Matrix1234FloatBoost::KroneckerProduct2Vocs(F1prime, F2prime);
			//Matrix1234FloatBoost::MatrixPrintRowMajor(t1, std::cout);

			FloatBoostReturnMapHandle rt1 = t1.root->rootConnection.returnMapHandle;
			FloatBoostReturnMapHandle rt1temp;
			for (unsigned int i = 0; i < rt1.Size(); i++)
			{
				if (rt1[i] == -1)
					rt1temp.AddToEnd(0);
				else
					rt1temp.AddToEnd(rt1[i]);
			}
			rt1temp.Canonicalize();
			CFLOBDDNodeHandle nt1 = *(t1.root->rootConnection.entryPointHandle);
			ReductionMapHandle inducedReductionMapHandle;
			FloatBoostReturnMapHandle inducedReturnMap;
			rt1temp.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_n = nt1.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			//CFLOBDDNodeHandle::DisposeOfReduceCache();
			t1 = CFLOBDD_FLOAT_BOOST(new CFLOBDDTopNodeFloatBoost(reduced_n, inducedReturnMap));

			CFLOBDD_FLOAT_BOOST t2 = Matrix1234FloatBoost::KroneckerProduct2Vocs(F1Permute, F2Permute);
			ReductionMapHandle inducedReductionMapHandle_2;
			FloatBoostReturnMapHandle inducedReturnMap_2;
			FloatBoostReturnMapHandle rt2 = t2.root->rootConnection.returnMapHandle;
			FloatBoostReturnMapHandle rt2temp;
			for (unsigned int i = 0; i < rt2.Size(); i++)
			{
				if (rt2[i] == 1)
					rt2temp.AddToEnd(0);
				else if (rt2[i] == -1)
					rt2temp.AddToEnd(1);
				else
					rt2temp.AddToEnd(rt2[i]);
			}
			rt2temp.Canonicalize();
			CFLOBDDNodeHandle nt2 = *(t2.root->rootConnection.entryPointHandle);
			rt2temp.InducedReductionAndReturnMap(inducedReductionMapHandle_2, inducedReturnMap_2);
			//CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_n_2 = nt2.Reduce(inducedReductionMapHandle_2, inducedReturnMap_2.Size());
			//CFLOBDDNodeHandle::DisposeOfReduceCache();
			t2 = CFLOBDD_FLOAT_BOOST(new CFLOBDDTopNodeFloatBoost(reduced_n_2, inducedReturnMap_2));
			CFLOBDD_FLOAT_BOOST ans = t1 + t2;
			return ans;
		}

		CFLOBDD_FLOAT_BOOST Func2To1CFLOBDDMatrix_DivideAndConquer(std::string s, std::mt19937 mt)
		{
			if (s.length() == 2)
			{
				if (s == "00"){
					std::vector<int> v;
					for (int i = 0; i < 4; i++)
						v.push_back(mt() % 4);
					/*
					v[0] = 0;
					v[1] = 0;
					v[2] = 0;
					v[3] = 0;*/
					
					return Matrix1234FloatBoost::Func2To1CFLOBDDMatrix(s, v);
				}
				else
				{
					std::vector<int> v;
					for (int i = 0; i < 2; i++)
						v.push_back(mt() % 4);
					
					/*
					v[0] = 0;
					v[1] = 1;*/
					
					return Matrix1234FloatBoost::Func2To1CFLOBDDMatrix(s, v);
				}
			} 
			/*
			else if (s.length() <= 16)
			{
				return Matrix1234FloatBoost::Func2To1CFLOBDDMatrix(s, std::vector<int>());
			}*/
			CFLOBDD_FLOAT_BOOST F1 = Func2To1CFLOBDDMatrix_DivideAndConquer(s.substr(0, s.length() / 2), mt);
			CFLOBDD_FLOAT_BOOST F2 = Func2To1CFLOBDDMatrix_DivideAndConquer(s.substr(s.length() / 2), mt);
			CFLOBDD_FLOAT_BOOST ans = Create2To1Func(F1, s.substr(0, s.length() / 2), F2, s.substr(s.length() / 2), mt);
			return ans;
		}

		CFLOBDD_FLOAT_BOOST NormalizeOutputTo1(CFLOBDD_FLOAT_BOOST c)
		{
			ReductionMapHandle inducedReductionMapHandle_2;
			FloatBoostReturnMapHandle inducedReturnMap_2;
			FloatBoostReturnMapHandle rt = c.root->rootConnection.returnMapHandle;
			FloatBoostReturnMapHandle rttemp;
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
			//CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_n_2 = nt.Reduce(inducedReductionMapHandle_2, inducedReturnMap_2.Size());
			//CFLOBDDNodeHandle::DisposeOfReduceCache();
			return CFLOBDD_FLOAT_BOOST(new CFLOBDDTopNodeFloatBoost(reduced_n_2, inducedReturnMap_2));
		}

		boost::multiprecision::cpp_int getNumberFromBits(std::string s){
			boost::multiprecision::cpp_int num = 0;
			for (int i = s.length() - 1; i >= 0; i--){
				num = 2 * num + (s[i] == '1');
			}
			return num;
		}

		boost::multiprecision::cpp_int compute_xor(boost::multiprecision::cpp_int num, std::string s){
			boost::multiprecision::cpp_int ans = 0;
			int i = s.length() - 1;
			while (i >= 0){
				bool v = (num % 2).convert_to<bool>();
				ans =  (v ^ (s[i] - '0'))*boost::multiprecision::pow(boost::multiprecision::cpp_int(2), s.length() - 1 - i) + ans;
				num = (num - num % 2) / 2;
				i--;
			}
			return ans;
		}

		std::vector<CFLOBDD_FLOAT_BOOST> getSizeTwoMatrices(std::vector<CFLOBDD_FLOAT_BOOST> baseMatrices,
			boost::multiprecision::cpp_int row, boost::multiprecision::cpp_int column, boost::multiprecision::cpp_int size){
			std::vector<CFLOBDD_FLOAT_BOOST> sizeTwoMatrices;

			if (size == 2){
				int r = row.convert_to<int>();
				int c = column.convert_to<int>();
				sizeTwoMatrices.push_back(baseMatrices[2 * r + c]);
				return sizeTwoMatrices;
			}

			int r = (row >= (size / 2));
			int c = (column >= (size / 2));
			CFLOBDD_FLOAT_BOOST t1 = baseMatrices[2 * r + c];
			sizeTwoMatrices.push_back(t1);
			boost::multiprecision::cpp_int new_row = (r == 1) ? row - size / 2 : row;
			boost::multiprecision::cpp_int new_column = (c == 1) ? column - size / 2 : column;
			std::vector<CFLOBDD_FLOAT_BOOST> t2 = getSizeTwoMatrices(baseMatrices, new_row, new_column, size / 2);
			sizeTwoMatrices.insert(sizeTwoMatrices.end(), t2.begin(), t2.end());

			return sizeTwoMatrices;
		}

		CFLOBDD_FLOAT_BOOST formMatrixFromSizeTwoMatrices(std::vector<CFLOBDD_FLOAT_BOOST> sizeTwoMatrices, unsigned int start, unsigned int end){
			if (start + 1 == end){
				return Matrix1234FloatBoost::KroneckerProduct2Vocs(sizeTwoMatrices[start], sizeTwoMatrices[end]);
			}
			unsigned int mid = (end - start) / 2 + start;
			CFLOBDD_FLOAT_BOOST t1 = formMatrixFromSizeTwoMatrices(sizeTwoMatrices, start, mid);
			CFLOBDD_FLOAT_BOOST t2 = formMatrixFromSizeTwoMatrices(sizeTwoMatrices, mid + 1, end);
			return Matrix1234FloatBoost::KroneckerProduct2Vocs(t1, t2);
		}

		CFLOBDD_FLOAT_BOOST formMatrixWithOneNonZeroValue(std::vector<CFLOBDD_FLOAT_BOOST> baseMatrices,
			boost::multiprecision::cpp_int row, boost::multiprecision::cpp_int column, boost::multiprecision::cpp_int size){
			std::vector<CFLOBDD_FLOAT_BOOST> sizeTwoMatrices = getSizeTwoMatrices(baseMatrices, row, column, size);
			return formMatrixFromSizeTwoMatrices(sizeTwoMatrices, 0, sizeTwoMatrices.size()-1);
		}

		std::vector<CFLOBDD_FLOAT_BOOST> formMatrixWithOneNonZeroValueUtil(){
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(0, 0);
			CFLOBDD_FLOAT_BOOST e1 = VectorFloatBoost::MkBasisVector(0, 1);
			e0 = VectorFloatBoost::MkVectorWithVoc12(e0);
			e1 = VectorFloatBoost::MkVectorWithVoc12(e1);
			CFLOBDD_FLOAT_BOOST UL = VectorFloatBoost::KroneckerProduct(e0, e0);
			CFLOBDD_FLOAT_BOOST UR = VectorFloatBoost::KroneckerProduct(e0, e1);
			CFLOBDD_FLOAT_BOOST LL = VectorFloatBoost::KroneckerProduct(e1, e0);
			CFLOBDD_FLOAT_BOOST LR = VectorFloatBoost::KroneckerProduct(e1, e1);
			std::vector<CFLOBDD_FLOAT_BOOST> baseMatrices;
			baseMatrices.push_back(UL);
			baseMatrices.push_back(UR);
			baseMatrices.push_back(LL);
			baseMatrices.push_back(LR);
			return baseMatrices;
		}

		std::string to_string(boost::multiprecision::cpp_int num){
			std::string ans = "";
			while (num > 0){
				int v = (num % 10).convert_to<int>();
				ans = std::to_string(v) + ans;
				num = (num - num % 10) / 10;
			}
			return ans;
		}



		CFLOBDD_FLOAT_BOOST Func2To1CFLOBDDMatrix(std::string s, std::vector<int> v = std::vector<int>()){
			assert(v.size() <= 256);
			unsigned int numOfBits = s.length();

			boost::multiprecision::cpp_int sizeofMatrix = boost::multiprecision::pow(boost::multiprecision::cpp_int(2), numOfBits);

			// Used for mapping rows in matrix based on s vector
			std::unordered_map<std::string, boost::multiprecision::cpp_int> halfMatrixValuesMap;

			std::vector<CFLOBDD_FLOAT_BOOST> baseMatrices = formMatrixWithOneNonZeroValueUtil();
			int v_index = 0;
			boost::multiprecision::cpp_int random_value;
			if (v.size() != 0)
				random_value = v[v_index++];
			else
				random_value = gen() % sizeofMatrix;
			halfMatrixValuesMap.insert(std::make_pair(to_string(compute_xor(0, s)), random_value));
			CFLOBDD_FLOAT_BOOST Func2To1Matrix = formMatrixWithOneNonZeroValue(baseMatrices, 0, random_value, sizeofMatrix);
			// generate 2^(n-1) random numbers
			for (boost::multiprecision::cpp_int i = 1; i < sizeofMatrix; i++)
			{
				std::string strIndex = to_string(i);
				if (halfMatrixValuesMap.find(strIndex) == halfMatrixValuesMap.end()){
					boost::multiprecision::cpp_int random_value;
					if (v.size() != 0)
						random_value = v[v_index++];
					else
						random_value = gen() % sizeofMatrix;
					halfMatrixValuesMap.insert(std::make_pair(to_string(compute_xor(i, s)), random_value));
					Func2To1Matrix = Func2To1Matrix + formMatrixWithOneNonZeroValue(baseMatrices, i, random_value, sizeofMatrix);
				}
				else{
					Func2To1Matrix = Func2To1Matrix + formMatrixWithOneNonZeroValue(baseMatrices, i, halfMatrixValuesMap[strIndex], sizeofMatrix);
				}
			}

			return Func2To1Matrix;

		}

		CFLOBDD_FLOAT_BOOST MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled){
			return CFLOBDD_FLOAT_BOOST(MkCNOTTopNode(level, n, controller, controlled));
		}

		CFLOBDD_FLOAT_BOOST ApplyExchangeAndIdentity(std::string s){
			if (s.find('1') == std::string::npos){
				return MkIdRelationInterleaved(ceil(log2(s.length())) + 1);
			}
			else if (s.find('0') == std::string::npos){
				return MkExchangeInterleaved(ceil(log2(s.length())) + 1);
			}
			CFLOBDD_FLOAT_BOOST F1 = ApplyExchangeAndIdentity(s.substr(0, s.length() / 2));
			CFLOBDD_FLOAT_BOOST F2 = ApplyExchangeAndIdentity(s.substr(s.length() / 2));
			return KroneckerProduct2Vocs(F1, F2);
		}

		CFLOBDD_FLOAT_BOOST CreateBalancedFn(int n, std::mt19937 mt){
			std::string s(2*n, '0');
			for (unsigned int i = 0; i < n; i++)
				s[i] = mt() % 2 ? '1' : '0';
			
			unsigned int level = ceil(log2(n)) + 2;
			CFLOBDD_FLOAT_BOOST F = ApplyExchangeAndIdentity(s);
			CFLOBDD_FLOAT_BOOST C = Matrix1234FloatBoost::MkCNOT(level, n, 0, n);
			for (int i = 1; i < n; i++){
				CFLOBDD_FLOAT_BOOST tmp = Matrix1234FloatBoost::MkCNOT(level, n, i, n);
				C = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(C, tmp);
			}
			std::cout << "C created" << std::endl;
			CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(C, F);
			ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(F, ans);
			return ans;
		}

		CFLOBDD_FLOAT_BOOST ComputeShortestPath(CFLOBDD_FLOAT_BOOST c1, CFLOBDD_FLOAT_BOOST c2){
			return CFLOBDD_FLOAT_BOOST(ComputeShortestPathTop(c1.root, c2.root));
		}
	}
}

