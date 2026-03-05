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

namespace CFL_OBDD {

	namespace Matrix1234Double {

		void Matrix1234InitializerTop()
		{
			Matrix1234InitializerNode();
			return;
		}

		// Create representation of identity relation (with interleaved variable order).
		// That is, input (x0,y0,x1,y1,...,xN,yN) yield Id[(x0,x1,...,xN)][(y0,y1,...,yN)]
		// which equals 1 iff xi == yi, for 0 <= i <= N.
		CFLOBDDTopNodeDoubleRefPtr MkIdRelationInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m10;

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m10);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr MkNegationMatrixInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m01;

			tempHandle = MkNegationMatrixInterleavedNode(i);
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr MkCNOTInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m01;

			tempHandle = MkCNOTInterleavedNode(i);
			m01.AddToEnd(1);
			m01.AddToEnd(0);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr ReverseColumnsTop(CFLOBDDTopNodeDoubleRefPtr n)
		{
			CFLOBDDNodeHandle tmp = ReverseColumnsNode(*(n->rootConnection.entryPointHandle));
			return new CFLOBDDTopNodeDouble(tmp, n->rootConnection.returnMapHandle);
		}

		CFLOBDDTopNodeDoubleRefPtr MatrixTransposeTop(CFLOBDDTopNodeDoubleRefPtr n)
		{
			std::unordered_map<CFLOBDDNodeHandle, std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>, 
				CFLOBDDNodeHandle::CFLOBDDNodeHandle_Hash> hashMap;
			auto pc = MatrixTransposeNode(hashMap, *(n->rootConnection.entryPointHandle));
			CFLOBDDNodeHandle temp = pc.first;
			ReductionMapHandle reductionMapHandle;
			DoubleReturnMapHandle v;
			for (unsigned int i = 0; i < pc.second.Size(); i++){
				reductionMapHandle.AddToEnd(i);
				v.AddToEnd(n->rootConnection.returnMapHandle[pc.second[i]]);
			}
			reductionMapHandle.Canonicalize();
			v.Canonicalize();
			temp = temp.Reduce(reductionMapHandle, v.Size(), true);
			return new CFLOBDDTopNodeDouble(temp, v);
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDDTopNodeDoubleRefPtr MkWalshInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeDouble::maxLevel);

			tempHandle = MkWalshInterleavedNode(i);
			/*m.AddToEnd(1.0 / ((1 << (1 << (i - 2)))));
			m.AddToEnd(-1.0 / ((1 << (1 << (i - 2)))));*/
			/*m.AddToEnd(1.0 / (pow(2, pow(2, i - 2))));
			m.AddToEnd(-1.0 / (pow(2, pow(2, i - 2))));*/
			m.AddToEnd(1.0);
			m.AddToEnd(-1.0);
			m.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m);
			return v;
		}


		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for two
		// additional vocabularies
		CFLOBDDTopNodeDoubleRefPtr MkWalshVoc12Top(unsigned int i)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeDouble::maxLevel);

			tempHandle = MkWalshVoc12Node(i);
			m.AddToEnd(1/sqrt((1 << (1 << (i-2)))));
			m.AddToEnd(-1/sqrt((1 << (1 << (i - 2)))));
			m.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m);
			return v;
		}


		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeDoubleRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeDoubleRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs12To34Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}


		CFLOBDDTopNodeDoubleRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeDoubleRefPtr n)
		{
			assert(1 <= n->rootConnection.entryPointHandle->handleContents->level);
			assert(n->rootConnection.entryPointHandle->handleContents->level < CFLOBDDTopNodeDouble::maxLevel);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = PromoteInterleavedTo12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeDoubleRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = Demote12ToInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeDoubleRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeDoubleRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc42Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal:
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDDTopNodeDoubleRefPtr MkDetensorConstraintInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeDouble::maxLevel);

			tempHandle = MkDetensorConstraintInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m);
			return v;
		}

		// Vocabulary projection
		// CFLOBDDTopNodeDoubleRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeDoubleRefPtr n)
		// {
		// 	assert(n->rootConnection.entryPointHandle.handleContents->level >= 2);

		// 	CFLOBDDLinearMapMemoTableRefPtr memoTable = new CFLOBDDLinearMapMemoTable;

		// 	CFLOBDDTopNodeLinearMapRefPtr temp = MatrixProjectVoc23Node(memoTable, n->rootConnection.entryPointHandle, TopLevelVisit);

		// 	// Begin construction of the appropriate CFLOBDDReturnMapHandle by applying each LinearMap in temp's
		// 	// CFLOBDDLinearMapTupleHandle to n's CFLOBDDReturnMapHandle
		// 	DoubleReturnMapHandle rmh;
		// 	for (unsigned int i = 0; i < temp->rootConnection.returnMapHandle.Size(); i++) {
		// 		rmh.AddToEnd(ApplyLinearMap(temp->rootConnection.returnMapHandle[i], n->rootConnection.returnMapHandle));
		// 	}

		// 	CFLOBDDNodeHandle tempHandle = temp->rootConnection.entryPointHandle;

		// 	// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
		// 	ReductionMapHandle inducedReductionMapHandle;
		// 	DoubleReturnMapHandle inducedReturnMap;
		// 	rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
		// 	//     CFLOBDDNodeHandle::InitReduceCache();
		// 	CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
		// 	//     CFLOBDDNodeHandle::DisposeOfReduceCache();

		// 	// Create and return CFLOBDDTopNode
		// 	return(new CFLOBDDTopNodeDouble(reduced_tempHandle, inducedReturnMap));
		// }

		CFLOBDDTopNodeDoubleRefPtr ConvertIntToDoubleTop(CFLOBDDTopNodeIntRefPtr c)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m;

			tempHandle = *(c->rootConnection.entryPointHandle);
			CFLOBDDReturnMapHandle r = c->rootConnection.returnMapHandle;
			for (unsigned i = 0; i < r.Size(); i++)
			{
				m.AddToEnd((double)r[i]);
			}
			m.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m);
			return v;
		}

		//
		// Print the matrix in row-major order
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 2 && level <= 4) {
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
						double b;
						b = n->EvaluateIteratively(a);
						if (b != 0)
							out << i << " " << j << std::endl;
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
		void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 1 && level <= 4) {
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
						double b;
						b = n->EvaluateIteratively(a);
						if (b != 0)
							out << i << " " << j << " " << b << std::endl;
						//out << b << " ";
					}
					//out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [1 .. 4]" << std::endl;
			}
		}

		CFLOBDDTopNodeDoubleRefPtr SMatrixTop(std::string s)
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
			
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle m;
			tempHandle = SMatrixNode(stemp);
			m.AddToEnd(1);
			if (tempHandle.handleContents->numExits == 2)
				m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNodeDouble(tempHandle, m);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeDoubleRefPtr c)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToAConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeDouble(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeDoubleRefPtr c)
		{
			CFLOBDDTopNodeDoubleRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToBConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeDouble(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeDoubleRefPtr c1, CFLOBDDTopNodeDoubleRefPtr c2)
		{
			std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash> hashMap;
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4Node(hashMap, *(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle));
			DoubleReturnMapHandle v;
			std::unordered_map<double, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				double val = 0;
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					double factor = j.second.convert_to<double>();
					val = val + factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]);
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
			// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
			//ReductionMapHandle inducedReductionMapHandle;
			//DoubleReturnMapHandle inducedReturnMap;
			//v.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			// Create and return CFLOBDDTopNode
			//return(new CFLOBDDTopNodeDouble(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeDouble(reduced_tempHandle, v));
		}

	}
}

