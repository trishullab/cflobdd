#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <bitset>
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include "wmatrix1234_top_node_fb_mul.h"
#include "wmatrix1234_fb_mul_node.h"

namespace CFL_OBDD {

	namespace WeightedMatrix1234FloatBoostMul {

		void Matrix1234InitializerTop()
		{
			Matrix1234InitializerNode();
			return;
		}

		// Create representation of identity relation (with interleaved variable order).
		// That is, input (x0,y0,x1,y1,...,xN,yN) yield Id[(x0,x1,...,xN)][(y0,y1,...,yN)]
		// which equals 1 iff xi == yi, for 0 <= i <= N.
		WeightedCFLOBDDTopNodeFloatBoostRefPtr MkIdRelationInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFloatBoostRefPtr v;
			WeightedCFLOBDDFloatBoostMulNodeHandle tempHandle;
			FloatBoostReturnMapHandle m10;

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFloatBoost(tempHandle, m10);
			return v;
		}

		WeightedCFLOBDDTopNodeFloatBoostRefPtr MkNegationMatrixInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFloatBoostRefPtr v;
			// WeightedCFLOBDDFloatBoostMulNodeHandle tempHandle;
			// FloatBoostReturnMapHandle m01;

			// tempHandle = MkNegationMatrixInterleavedNode(i);
			// m01.AddToEnd(0);
			// m01.AddToEnd(1);
			// m01.Canonicalize();
			// v = new WeightedCFLOBDDTopNodeFloatBoost(tempHandle, m01);
			return v;
		}

		WeightedCFLOBDDTopNodeFloatBoostRefPtr MkCNOTInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFloatBoostRefPtr v;
			// WeightedCFLOBDDFloatBoostMulNodeHandle tempHandle;
			// FloatBoostReturnMapHandle m01;

			// tempHandle = MkCNOTInterleavedNode(i);
			// m01.AddToEnd(1);
			// m01.AddToEnd(0);
			// m01.Canonicalize();
			// v = new WeightedCFLOBDDTopNodeFloatBoost(tempHandle, m01);
			return v;
		}

		WeightedCFLOBDDTopNodeFloatBoostRefPtr MkExchangeInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFloatBoostRefPtr v;
			WeightedCFLOBDDFloatBoostMulNodeHandle tempHandle;
			FloatBoostReturnMapHandle m01;

			tempHandle = MkExchangeInterleavedNode(i);
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFloatBoost(tempHandle, m01);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		WeightedCFLOBDDTopNodeFloatBoostRefPtr MkWalshInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFloatBoostRefPtr v;
			WeightedCFLOBDDFloatBoostMulNodeHandle tempHandle;
			FloatBoostReturnMapHandle m;

			tempHandle = MkWalshInterleavedNode(i);
			// auto val = boost::multiprecision::pow(BIG_FLOAT(sqrt(2)), i).convert_to<BIG_FLOAT>();
			m.AddToEnd(1);
			m.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFloatBoost(tempHandle, m, 1.0);
			return v;
		}

        WeightedCFLOBDDTopNodeFloatBoostRefPtr KroneckerProduct2VocsTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr m1, WeightedCFLOBDDTopNodeFloatBoostRefPtr m2)
        {
            unsigned int level = m1->rootConnection.entryPointHandle->handleContents->level;
            if (m1->rootConnection.factor == 0 || m2->rootConnection.factor == 0)
            {
                FloatBoostReturnMapHandle m;
                m.AddToEnd(0);
                m.Canonicalize();
                return new WeightedCFLOBDDTopNodeFloatBoost(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level + 1], m);
            }
            int zero_index_m1 = m1->rootConnection.returnMapHandle.LookupInv(0);
            int zero_index_m2 = m2->rootConnection.returnMapHandle.LookupInv(0);

            WeightedCFLOBDDTopNodeFloatBoostRefPtr v;
            WeightedCFLOBDDFloatBoostMulNodeHandle tempHandle; 
            FloatBoostReturnMapHandle m;
            tempHandle = KroneckerProduct2VocsNode(*(m1->rootConnection.entryPointHandle), *(m2->rootConnection.entryPointHandle), zero_index_m1, zero_index_m2);
            FloatBoostReturnMapHandle m01, m10, m_1;
            m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
            m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
            m_1.AddToEnd(1); m_1.Canonicalize();

            if (m1->rootConnection.returnMapHandle == m01)
                m = m01;
            else if (m2->rootConnection.returnMapHandle == m01)
                m = m01;
            else
                m = m10;

            if (m1->rootConnection.returnMapHandle == m_1 && m2->rootConnection.returnMapHandle == m_1)
                m = m_1;
            BIG_FLOAT f1 = m1->rootConnection.factor;
            BIG_FLOAT f2 = m2->rootConnection.factor;
            BIG_FLOAT f = f1 * f2;
            v = new WeightedCFLOBDDTopNodeFloatBoost(tempHandle, m, f);
            return v;
        }


		//
		// Print the matrix in row-major order
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out)
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
						double b1 = 0;//n->EvaluateIteratively(a).convert_to<double>();
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
		void MatrixPrintRowMajorInterleavedTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out)
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

						double b1 = 0;//n->EvaluateIteratively(a).convert_to<double>();
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

		
		
		WeightedCFLOBDDTopNodeFloatBoostRefPtr MatrixMultiplyV4TopNode(WeightedCFLOBDDTopNodeFloatBoostRefPtr c1, WeightedCFLOBDDTopNodeFloatBoostRefPtr c2)
		{
            int zero_exit_1 = -1, zero_exit_2 = -1;
            zero_exit_1 = c1->rootConnection.returnMapHandle.LookupInv(0);
            zero_exit_2 = c2->rootConnection.returnMapHandle.LookupInv(0);
			auto c = MatrixMultiplyV4Node(*(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle), zero_exit_1, zero_exit_2);
			FloatBoostReturnMapHandle v;
			boost::unordered_map<BIG_FLOAT, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
            WeightedValuesListHandle<BIG_FLOAT> valList;
            WeightedCFLOBDDFloatBoostMulNodeHandle n = std::get<0>(c);
            CFLOBDDMatMultMapHandle n_return = std::get<1>(c);
            BIG_FLOAT n_factor = std::get<2>(c);
            bool first = true;
            BIG_FLOAT common_f = 1.0;
			for (unsigned int i = 0; i < n_return.Size(); i++){
				WeightedMatMultMapHandle r = n_return[i];
				BIG_FLOAT val = 0;
				for (auto &j : r.mapContents->map){
					long int index1 = j.first.first;
					long int index2 = j.first.second;
					BIG_FLOAT factor(j.second);
                    if (!(index1 == -1 && index2 == -1))
					    val = val + (factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]));
				}
                if (first && val != 0)
                {
                    common_f = val;
                    val = 1.0;
                    first = false;
                }
                else
                {
                    val = val/common_f;
                }
                BIG_FLOAT val_to_check = (val == 0)? 0 : 1;
				if (reductionMap.find(val_to_check) == reductionMap.end()){
					v.AddToEnd(val_to_check);
					reductionMap.insert(std::make_pair(val_to_check, v.Size() - 1));
					reductionMapHandle.AddToEnd(v.Size() - 1);
                    valList.AddToEnd(val);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[val_to_check]);
                    valList.AddToEnd(val);
				}
			}
			v.Canonicalize();
			reductionMapHandle.Canonicalize();
            valList.Canonicalize();
            // n.print(std::cout);
            // std::cout << valList << std::endl;
			auto g = n.Reduce(reductionMapHandle, v.Size(), valList, true);
			// Create and return CFLOBDDTopNode
			//return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
			return(new WeightedCFLOBDDTopNodeFloatBoost(g.first, v, g.second * common_f * c1->rootConnection.factor * c2->rootConnection.factor));
		}

		WeightedCFLOBDDTopNodeFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled)
		{
			WeightedCFLOBDDFloatBoostMulNodeHandle c = MkCNOTNode(level, n, controller, controlled);
			FloatBoostReturnMapHandle m;
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			return new WeightedCFLOBDDTopNodeFloatBoost(c, m, 1.0);
		}

	}
}

