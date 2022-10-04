#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <bitset>
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "matrix1234_node.h"
#include "matrix1234_float_boost_top_node.h"
// #include "general_map.h"
#include "zero_indices_map.h"

namespace CFL_OBDD {

	namespace Matrix1234FloatBoost {

		void Matrix1234InitializerTop()
		{
			Matrix1234InitializerNode();
			return;
		}

		// Create representation of identity relation (with interleaved variable order).
		// That is, input (x0,y0,x1,y1,...,xN,yN) yield Id[(x0,x1,...,xN)][(y0,y1,...,yN)]
		// which equals 1 iff xi == yi, for 0 <= i <= N.
		CFLOBDDTopNodeFloatBoostRefPtr MkIdRelationInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m10;

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m10);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkNegationMatrixInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m01;

			tempHandle = MkNegationMatrixInterleavedNode(i);
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkCNOTInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m01;

			tempHandle = MkCNOTInterleavedNode(i);
			m01.AddToEnd(1);
			m01.AddToEnd(0);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkExchangeInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m01;

			tempHandle = MkExchangeInterleavedNode(i);
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr ReverseColumnsTop(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			CFLOBDDNodeHandle tmp = ReverseColumnsNode(*(n->rootConnection.entryPointHandle));
			return new CFLOBDDTopNodeFloatBoost(tmp, n->rootConnection.returnMapHandle);
		}

		CFLOBDDTopNodeFloatBoostRefPtr MatrixTransposeTop(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			std::unordered_map<CFLOBDDNodeHandle, std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>, 
				CFLOBDDNodeHandle::CFLOBDDNodeHandle_Hash> hashMap;
			auto pc = MatrixTransposeNode(hashMap, *(n->rootConnection.entryPointHandle));
			CFLOBDDNodeHandle temp = pc.first;
			ReductionMapHandle reductionMapHandle;
			FloatBoostReturnMapHandle v;
			for (unsigned int i = 0; i < pc.second.Size(); i++){
				reductionMapHandle.AddToEnd(i);
				v.AddToEnd(n->rootConnection.returnMapHandle[pc.second[i]]);
			}
			reductionMapHandle.Canonicalize();
			v.Canonicalize();
			temp = temp.Reduce(reductionMapHandle, v.Size(), true);
			return new CFLOBDDTopNodeFloatBoost(temp, v);
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDDTopNodeFloatBoostRefPtr MkWalshInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeFloatBoost::maxLevel);

			tempHandle = MkWalshInterleavedNode(i);
			/*m.AddToEnd(1.0 / ((1 << (1 << (i - 2)))));
			m.AddToEnd(-1.0 / ((1 << (1 << (i - 2)))));*/
			/*m.AddToEnd(1.0 / (pow(2, pow(2, i - 2))));
			m.AddToEnd(-1.0 / (pow(2, pow(2, i - 2))));*/
			auto val = boost::multiprecision::pow(boost::multiprecision::cpp_int(2), pow(2, i - 1)).convert_to<BIG_FLOAT>();
			m.AddToEnd(1.0/val);
			m.AddToEnd(-1.0/val);
			m.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m);
			return v;
		}


		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for two
		// additional vocabularies
		CFLOBDDTopNodeFloatBoostRefPtr MkWalshVoc12Top(unsigned int i)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeFloatBoost::maxLevel);

			tempHandle = MkWalshVoc12Node(i);
			m.AddToEnd(1 / sqrt((1 << (1 << (i - 2)))));
			m.AddToEnd(-1 / sqrt((1 << (1 << (i - 2)))));
			m.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m);
			return v;
		}


		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs12To34Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}


		CFLOBDDTopNodeFloatBoostRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			assert(1 <= n->rootConnection.entryPointHandle->handleContents->level);
			assert(n->rootConnection.entryPointHandle->handleContents->level < CFLOBDDTopNodeFloatBoost::maxLevel);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = PromoteInterleavedTo12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr PromoteInterleavedTo13Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			assert(1 <= n->rootConnection.entryPointHandle->handleContents->level);
			assert(n->rootConnection.entryPointHandle->handleContents->level < CFLOBDDTopNodeFloatBoost::maxLevel);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = PromoteInterleavedTo13Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = Demote12ToInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc42Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal:
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDDTopNodeFloatBoostRefPtr MkDetensorConstraintInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeFloatBoost::maxLevel);

			tempHandle = MkDetensorConstraintInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m);
			return v;
		}

		// // Vocabulary projection
		// CFLOBDDTopNodeFloatBoostRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		// {
		// 	assert(n->rootConnection.entryPointHandle.handleContents->level >= 2);

		// 	CFLOBDDLinearMapMemoTableRefPtr memoTable = new CFLOBDDLinearMapMemoTable;
		// 	nodeNum = 0;
		// 	cacheNodeNum = 0;
		// 	notL1NodeNum = 0;
		// 	/*CFLOBDDNodeHandle::DisposeOfReduceCache();
		// 	CFLOBDDNodeHandle::InitReduceCache();*/
		// 	CFLOBDDTopNodeLinearMapRefPtr temp = MatrixProjectVoc23Node(memoTable, n->rootConnection.entryPointHandle, TopLevelVisit);

		// 	// Begin construction of the appropriate CFLOBDDReturnMapHandle by applying each LinearMap in temp's
		// 	// CFLOBDDLinearMapTupleHandle to n's CFLOBDDReturnMapHandle
		// 	FloatBoostReturnMapHandle rmh;
		// 	for (unsigned int i = 0; i < temp->rootConnection.returnMapHandle.Size(); i++) {
		// 		rmh.AddToEnd(ApplyLinearMap(temp->rootConnection.returnMapHandle[i], n->rootConnection.returnMapHandle));
		// 	}

		// 	CFLOBDDNodeHandle tempHandle = temp->rootConnection.entryPointHandle;

		// 	// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
		// 	ReductionMapHandle inducedReductionMapHandle;
		// 	FloatBoostReturnMapHandle inducedReturnMap;
		// 	rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
		// 	//CFLOBDDNodeHandle::InitReduceCache();
		// 	CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
		// 	//CFLOBDDNodeHandle::DisposeOfReduceCache();
		// 	// Create and return CFLOBDDTopNode
		// 	return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
		// }
		
		
		
		// CFLOBDDTopNodeFloatBoostRefPtr AddMatrixRowsTopNode(CFLOBDDTopNodeFloatBoostRefPtr n)
		// {
		// 	CFLOBDDGeneralMapNodeMemoTableRefPtr memoTable = new CFLOBDDGeneralMapNodeMemoTable;
		// 	 CFLOBDD_GENERAL temp = AddMatrixRowsNode(memoTable, n->rootConnection.entryPointHandle);
		// 	 FloatBoostReturnMapHandle ret;
		// 	 for (unsigned int i = 0; i < temp.root->rootConnection.returnMapHandle.Size(); i++)
		// 	 {
		// 		 BIG_FLOAT val = 0;
		// 		 for (unsigned int j = 0; j < temp.root->rootConnection.returnMapHandle[i].Size(); j++)
		// 		 {
		// 			 std::pair<unsigned int, unsigned int> retVal = temp.root->rootConnection.returnMapHandle[i][j];
		// 			 val += n->rootConnection.returnMapHandle[retVal.first] * retVal.second;
		// 		 }
		// 		ret.AddToEnd(val);
		// 	 }
		// 	 ret.Canonicalize();

		// 	 CFLOBDDNodeHandle tempHandle = temp.root->rootConnection.entryPointHandle;

		// 	 // Perform reduction on tempHandle, with respect to the common elements that rmh maps together
		// 	 ReductionMapHandle inducedReductionMapHandle;
		// 	 FloatBoostReturnMapHandle inducedReturnMap;
		// 	 ret.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
		// 	 //CFLOBDDNodeHandle::InitReduceCache();
		// 	 CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
		// 	 //CFLOBDDNodeHandle::DisposeOfReduceCache();
		// 	 // Create and return CFLOBDDTopNode
		// 	 return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));

		// }
		
		CFLOBDDTopNodeFloatBoostRefPtr ConvertIntToFloatBoostTop(CFLOBDDTopNodeIntRefPtr c)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m;

			tempHandle = *(c->rootConnection.entryPointHandle);
			CFLOBDDReturnMapHandle r = c->rootConnection.returnMapHandle;
			for (unsigned i = 0; i < r.Size(); i++)
			{
				BIG_FLOAT r_bigfloat(r[i]);
				m.AddToEnd(r_bigfloat);
			}
			m.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m);
			return v;
		}

		//
		// Print the matrix in row-major order
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 2 && level <= 4 || true) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				for (unsigned long int i = 0UL; i < rows; i++) { // Vocs 1 and 3
					// Fill in the even positions of a 
					unsigned long int maskVoc1 = 1UL << (indexBits / 2);
					unsigned long int maskVoc3 = 1UL;
					for (int k = (indexBits / 2 - 1); k >= 0; k--) {
						a[4 * k] = (i & maskVoc1);
						a[4 * k + 2] = (i & maskVoc3);
						maskVoc1 = maskVoc1 << 1;
						maskVoc3 = maskVoc3 << 1;
					}
					for (unsigned long int j = 0UL; j < cols; j++) {  // Vocs 2 and 4
						// Fill in the odd positions of a
						unsigned long int maskVoc2 = 1UL << (indexBits / 2);
						unsigned long int maskVoc4 = 1UL;
						for (int k = (indexBits / 2 - 1); k >= 0; k--) {
							a[4 * k + 1] = (j & maskVoc2);
							a[4 * k + 3] = (j & maskVoc4);
							maskVoc2 = maskVoc2 << 1;
							maskVoc4 = maskVoc4 << 1;
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;
						double b1 = n->EvaluateIteratively(a).convert_to<double>();
						if (b1 != 0)
							out << std::bitset<4>(i) << " " << std::bitset<4>(j) << " " << b1 << std::endl;
						// std::string b;
						// Change precision
						// b = n->EvaluateIteratively(a).str(5);
						//out << b << " ";
					}
					//out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [2 .. 4]" << std::endl;
			}
		}

		// MatrixPrintRowMajorInterleaved
		//
		// Print the matrix in row-major order
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 1 && level <= 4 || true) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				for (unsigned long int i = 0UL; i < rows; i++) {
					// Fill in the even positions of a 
					unsigned long int mask1 = 1UL;
					for (int k = indexBits - 1; k >= 0; k--) {
						a[2 * k] = (i & mask1);
						mask1 = mask1 << 1;
					}
					for (unsigned long int j = 0UL; j < cols; j++) {
						// Fill in the odd positions of a
						unsigned long int mask2 = 1UL;
						for (int k = indexBits - 1; k >= 0; k--) {
							a[2 * k + 1] = (j & mask2);
							mask2 = mask2 << 1;
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;

						double b1 = n->EvaluateIteratively(a).convert_to<double>();
						//if (b1 != 0)
						//	out << std::bitset<2>(i) << " " << std::bitset<2>(j) << " " << b1 << std::endl;

						//std::string b;
						//b = n->EvaluateIteratively(a).str(5);
						out << b1 << " ";
					}
					out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [1 .. 4]" << std::endl;
			}
		}

		CFLOBDDTopNodeFloatBoostRefPtr SMatrixTop(std::string s)
		{
			std::string stemp;
			for (int i = 0; i < s.length(); i++)
			{
				stemp += "0";
			}
			for (int i = 0; i < s.length(); i++)
			{
				if (s[i] == '1')
				{
					stemp[i] = '1';
					break;
				}
			}
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m;

			tempHandle = SMatrixNode(stemp);
			m.AddToEnd(1);
			if (tempHandle.handleContents->numExits == 2)
				m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeFloatBoostRefPtr c)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToAConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeFloatBoostRefPtr c)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToBConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}
		
		CFLOBDDTopNodeFloatBoostRefPtr MkDistinctionTwoVarsTopNode(int x, int y, unsigned int var_level, unsigned int total_level)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle m;
			if (x == 0 && y == 0){
				m.AddToEnd(1);
				m.AddToEnd(0);
			}
			else{
				m.AddToEnd(0);
				m.AddToEnd(1);
			}
			m.Canonicalize();
			tempHandle = MkDistinctionTwoVarsNode(x,y,var_level, total_level);
			v = new CFLOBDDTopNodeFloatBoost(tempHandle, m);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MultiplyOperationTopNode(CFLOBDDTopNodeFloatBoostRefPtr c)
		{
			CFLOBDDTopNodeFloatBoostRefPtr v;
			auto m = MultiplyOperationNode(*(c->rootConnection.entryPointHandle));
			
			FloatBoostReturnMapHandle rmh;
			for (unsigned int i = 0; i < m.second.size(); i++) {
				BIG_FLOAT f = 0;
				for (unsigned int j = 0; j < m.second[i].size(); j++){
					f += m.second[i][j].first * c->rootConnection.returnMapHandle[m.second[i][j].second];
				}
				rmh.AddToEnd(f);
			}
			rmh.Canonicalize();
			CFLOBDDNodeHandle tempHandle = m.first;
			std::cout << tempHandle << std::endl;
			// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
			ReductionMapHandle inducedReductionMapHandle;
			FloatBoostReturnMapHandle inducedReturnMap;
			rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			//CFLOBDDNodeHandle::DisposeOfReduceCache();
			// Create and return CFLOBDDTopNode
			return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkMatrixMultiplyConstraintTopNode(unsigned int level)
		{
			CFLOBDDNodeHandle c = MkMatrixMultiplyConstraintNode(level);
			FloatBoostReturnMapHandle v;
			v.AddToEnd(1);
			v.AddToEnd(0);
			v.Canonicalize();
			return new CFLOBDDTopNodeFloatBoost(c, v);
		}
		
		
		CFLOBDDTopNodeFloatBoostRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeFloatBoostRefPtr c1, CFLOBDDTopNodeFloatBoostRefPtr c2)
		{
			std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash> hashMap;
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4Node(hashMap, *(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle));
			FloatBoostReturnMapHandle v;
			boost::unordered_map<BIG_FLOAT, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				BIG_FLOAT val = 0;
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					BIG_FLOAT factor(j.second);
					val = val + (factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]));
				}
				if (reductionMap.find(val) == reductionMap.end()){
					v.AddToEnd(val);
					reductionMap.insert(std::make_pair(val, v.Size() - 1));
					reductionMapHandle.AddToEnd(v.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[val]);
				}
			}
			
			v.Canonicalize();
			reductionMapHandle.Canonicalize();
			//std::cout << v << std::endl;
			CFLOBDDNodeHandle tempHandle = *(c->rootConnection.entryPointHandle);
			// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
			/*ReductionMapHandle inducedReductionMapHandle;
			FloatBoostReturnMapHandle inducedReturnMap;
			v.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());*/
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			// Create and return CFLOBDDTopNode
			//return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, v));
		}

		CFLOBDDTopNodeFloatBoostRefPtr MatrixMultiplyV4WithInfoTopNode(CFLOBDDTopNodeFloatBoostRefPtr c1, CFLOBDDTopNodeFloatBoostRefPtr c2)
		{
			std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash> hashMap;
			int c1_zero_index = -1, c2_zero_index = -1;
			c1_zero_index = c1->rootConnection.returnMapHandle.LookupInv(0);
			c2_zero_index = c2->rootConnection.returnMapHandle.LookupInv(0);
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4WithInfoNode
				(hashMap, *(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle), c1_zero_index, c2_zero_index);
			FloatBoostReturnMapHandle v;
			boost::unordered_map<BIG_FLOAT, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				BIG_FLOAT val = 0;
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					if (index1 != -1 && index2 != -1){
						BIG_FLOAT factor(j.second);
						val = val + (factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]));
					}
				}
				if (reductionMap.find(val) == reductionMap.end()){
					v.AddToEnd(val);
					reductionMap.insert(std::make_pair(val, v.Size() - 1));
					reductionMapHandle.AddToEnd(v.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[val]);
				}
			}

			v.Canonicalize();
			reductionMapHandle.Canonicalize();
			//std::cout << v << std::endl;
			CFLOBDDNodeHandle tempHandle = *(c->rootConnection.entryPointHandle);
			// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
			/*ReductionMapHandle inducedReductionMapHandle;
			FloatBoostReturnMapHandle inducedReturnMap;
			v.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());*/
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			// Create and return CFLOBDDTopNode
			//return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, v));
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled)
		{
			CFLOBDDNodeHandle c = MkCNOTNode(level, n, controller, controlled);
			FloatBoostReturnMapHandle m;
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			return new CFLOBDDTopNodeFloatBoost(c, m);
		}

		CFLOBDDTopNodeFloatBoostRefPtr ComputeShortestPathTop(CFLOBDDTopNodeFloatBoostRefPtr c1, CFLOBDDTopNodeFloatBoostRefPtr c2)
		{
			std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash> hashMap;
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4Node(hashMap, *(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle));
			FloatBoostReturnMapHandle v;
			boost::unordered_map<BIG_FLOAT, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				BIG_FLOAT val = std::numeric_limits<BIG_FLOAT>::infinity();
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					BIG_FLOAT factor(j.second);
					val = min(val, (factor * (c1->rootConnection.returnMapHandle[index1] + c2->rootConnection.returnMapHandle[index2])));
				}
				if (reductionMap.find(val) == reductionMap.end()){
					v.AddToEnd(val);
					reductionMap.insert(std::make_pair(val, v.Size() - 1));
					reductionMapHandle.AddToEnd(v.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[val]);
				}
			}

			v.Canonicalize();
			reductionMapHandle.Canonicalize();
			CFLOBDDNodeHandle tempHandle = *(c->rootConnection.entryPointHandle);
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, v));
		}

	}
}

