#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <bitset>
#include <unordered_map>
#include <boost/unordered_map.hpp>
#include "wmatrix1234_top_node_fourier_mul.h"
#include "wmatrix1234_fourier_mul_node.h"

namespace CFL_OBDD {

	namespace WeightedMatrix1234FourierMul {

		void Matrix1234InitializerTop()
		{
			Matrix1234InitializerNode();
			return;
		}

		// Create representation of identity relation (with interleaved variable order).
		// That is, input (x0,y0,x1,y1,...,xN,yN) yield Id[(x0,x1,...,xN)][(y0,y1,...,yN)]
		// which equals 1 iff xi == yi, for 0 <= i <= N.
		WeightedCFLOBDDTopNodeFourierRefPtr MkIdRelationInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m10;

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(fourierSemiring(1, 1));
			m10.AddToEnd(fourierSemiring(0, 1));
			m10.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m10, fourierSemiring(1, 1));
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkNegationMatrixInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			// WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			// FourierReturnMapHandle m01;

			// tempHandle = MkNegationMatrixInterleavedNode(i);
			// m01.AddToEnd(0);
			// m01.AddToEnd(1);
			// m01.Canonicalize();
			// v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01);
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkCNOTInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			// WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			// FourierReturnMapHandle m01;

			// tempHandle = MkCNOTInterleavedNode(i);
			// m01.AddToEnd(1);
			// m01.AddToEnd(0);
			// m01.Canonicalize();
			// v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01);
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkExchangeInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m01;

			tempHandle = MkExchangeInterleavedNode(i);
			m01.AddToEnd(fourierSemiring(0, 1));
			m01.AddToEnd(fourierSemiring(1, 1));
			m01.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01, fourierSemiring(1, 1));
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		WeightedCFLOBDDTopNodeFourierRefPtr MkWalshInterleavedTop(unsigned int i)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m;

			tempHandle = MkWalshInterleavedNode(i);
			// auto val = boost::multiprecision::pow(fourierSemiring(sqrt(2)), i).convert_to<fourierSemiring>();
			m.AddToEnd(fourierSemiring(1, 1));
			m.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m, fourierSemiring(1, 1));
			return v;
		}

        WeightedCFLOBDDTopNodeFourierRefPtr KroneckerProduct2VocsTop(WeightedCFLOBDDTopNodeFourierRefPtr m1, WeightedCFLOBDDTopNodeFourierRefPtr m2)
        {
            unsigned int level = m1->rootConnection.entryPointHandle->handleContents->level;
            if (m1->rootConnection.factor == fourierSemiring(0, 1) || m2->rootConnection.factor == fourierSemiring(0, 1))
            {
                FourierReturnMapHandle m;
                m.AddToEnd(fourierSemiring(0, 1));
                m.Canonicalize();
                return new WeightedCFLOBDDTopNodeFourier(WeightedCFLOBDDFourierMulNodeHandle::NoDistinctionNode_Ann[level + 1], m, fourierSemiring(0, 1));
            }
            int zero_index_m1 = m1->rootConnection.returnMapHandle.LookupInv(fourierSemiring(0, 1));
            int zero_index_m2 = m2->rootConnection.returnMapHandle.LookupInv(fourierSemiring(0, 1));
			fourierSemiring zero = fourierSemiring(0, 1);
			fourierSemiring one = fourierSemiring(1, 1);

            WeightedCFLOBDDTopNodeFourierRefPtr v;
            WeightedCFLOBDDFourierMulNodeHandle tempHandle; 
            FourierReturnMapHandle m;
            tempHandle = KroneckerProduct2VocsNode(*(m1->rootConnection.entryPointHandle), *(m2->rootConnection.entryPointHandle), zero_index_m1, zero_index_m2);
            FourierReturnMapHandle m01, m10, m_1;
            m01.AddToEnd(zero); m01.AddToEnd(one); m01.Canonicalize();
            m10.AddToEnd(one); m10.AddToEnd(zero); m10.Canonicalize();
            m_1.AddToEnd(one); m_1.Canonicalize();

            if (m1->rootConnection.returnMapHandle == m01)
                m = m01;
            else if (m2->rootConnection.returnMapHandle == m01)
                m = m01;
            else
                m = m10;

            if (m1->rootConnection.returnMapHandle == m_1 && m2->rootConnection.returnMapHandle == m_1)
                m = m_1;
            fourierSemiring f1 = m1->rootConnection.factor;
            fourierSemiring f2 = m2->rootConnection.factor;
            fourierSemiring f = f1 * f2;
            v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m, f);
            return v;
        }


		//
		// Print the matrix in row-major order
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out)
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
		void MatrixPrintRowMajorInterleavedTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out)
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

		
		
		WeightedCFLOBDDTopNodeFourierRefPtr MatrixMultiplyV4TopNode(WeightedCFLOBDDTopNodeFourierRefPtr c1, WeightedCFLOBDDTopNodeFourierRefPtr c2)
		{
            int zero_exit_1 = -1, zero_exit_2 = -1;
			fourierSemiring zero = fourierSemiring(0, 1);
			fourierSemiring one = fourierSemiring(1, 1);
            zero_exit_1 = c1->rootConnection.returnMapHandle.LookupInv(zero);
            zero_exit_2 = c2->rootConnection.returnMapHandle.LookupInv(zero);
			auto c = MatrixMultiplyV4Node(*(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle), zero_exit_1, zero_exit_2);
			FourierReturnMapHandle v;
			boost::unordered_map<fourierSemiring, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
            WeightedValuesListHandle<fourierSemiring> valList;
            WeightedCFLOBDDFourierMulNodeHandle n = std::get<0>(c);
            CFLOBDDMatMultMapHandle n_return = std::get<1>(c);
            fourierSemiring n_factor = std::get<2>(c);
            bool first = true;
			// n.print(std::cout);
            fourierSemiring common_f = one;
			for (unsigned int i = 0; i < n_return.Size(); i++){
				WeightedMatMultMapHandle r = n_return[i];
				fourierSemiring val = zero;
				for (auto &j : r.mapContents->map){
					long int index1 = j.first.first;
					long int index2 = j.first.second;
					fourierSemiring factor(j.second);
                    if (!(index1 == -1 && index2 == -1))
					    val = val + (factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]));
				}
                if (first && val != zero)
                {
                    common_f = val;
                    val = one;
                    first = false;
                }
                else
                {
                    val = val/common_f;
                }
                fourierSemiring val_to_check = (val == zero)? zero : one;
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
			auto g = n.Reduce(reductionMapHandle, v.Size(), valList, true);
			return(new WeightedCFLOBDDTopNodeFourier(g.first, v, g.second * common_f * c1->rootConnection.factor * c2->rootConnection.factor));
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled)
		{
			WeightedCFLOBDDFourierMulNodeHandle c = MkCNOTNode(level, n, controller, controlled);
			FourierReturnMapHandle m;
			m.AddToEnd(fourierSemiring(1, 1));
			m.AddToEnd(fourierSemiring(0, 1));
			m.Canonicalize();
			return new WeightedCFLOBDDTopNodeFourier(c, m, fourierSemiring(1, 1));
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkSwapGateTop(unsigned int i, long int c1, long int c2)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m01;

			tempHandle = MkSwapGateNode(i, c1, c2, -1);
			m01.AddToEnd(fourierSemiring(1, 1));
			m01.AddToEnd(fourierSemiring(0, 1));
			m01.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01, fourierSemiring(1, 1));
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkCPGateTop(unsigned int i, long int c1, long int c2, fourierSemiring theta)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m10;

			// std::unordered_map<std::string, CFLOBDDNodeHandle> cp_hashMap;
			std::unordered_map<std::string, WeightedCFLOBDDFourierMulNodeHandle> cp_hashMap;
			tempHandle = MkCPGateNode(cp_hashMap, i, c1, c2, theta);
			m10.AddToEnd(fourierSemiring(1, 1));
			m10.AddToEnd(fourierSemiring(0, 1));
			m10.Canonicalize();
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m10, fourierSemiring(1, 1));
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkCSwapGateTop(unsigned int level, long int c, long int i, long int j)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m01;
			tempHandle = MkCSwapGateNode(level, c, i, j, 1);
			m01.AddToEnd(fourierSemiring(1, 1));
			m01.AddToEnd(fourierSemiring(0, 1));
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01, fourierSemiring(1, 1));
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkRZGateTop(unsigned int level, fourierSemiring theta)
		{
			assert(level == 1);
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m01;
			tempHandle = MkRZGateNode(level, theta);
			m01.AddToEnd(fourierSemiring(1, 1));
			m01.AddToEnd(fourierSemiring(0, 1));
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01, fourierSemiring(1, 1));
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkCADDGateTop(unsigned int level, int c, int x, WeightedCFLOBDDTopNodeFourierRefPtr f)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m01;
			tempHandle = MkCADDGateNode(level, c, x, *(f->rootConnection.entryPointHandle), 1);
			m01.AddToEnd(fourierSemiring(1, 1));
			m01.AddToEnd(fourierSemiring(0, 1));
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01, fourierSemiring(1, 1));
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkCADDGate2Top(unsigned int level, int x, WeightedCFLOBDDTopNodeFourierRefPtr f)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle m01;
			tempHandle = MkCADDGate2Node(level, x, *(f->rootConnection.entryPointHandle), 1);
			m01.AddToEnd(fourierSemiring(1, 1));
			m01.AddToEnd(fourierSemiring(0, 1));
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, m01, fourierSemiring(1, 1));
			return v;
		}

		bool CheckIfIndexIsNonZeroTop(unsigned int level, int index, WeightedCFLOBDDTopNodeFourierRefPtr f)
		{
			return CheckIfIndexIsNonZeroNode(level, index, *(f->rootConnection.entryPointHandle), 1);
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkSetBToZeroTop(unsigned int level, WeightedCFLOBDDTopNodeFourierRefPtr f)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			tempHandle = MkSetBToZeroNode(level, *(f->rootConnection.entryPointHandle), 1);
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, f->rootConnection.returnMapHandle, fourierSemiring(1, 1));
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr ComputeIQFTTop(unsigned int level, WeightedCFLOBDDTopNodeFourierRefPtr f, BIG_INT N, int n)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			tempHandle = ComputeIQFTNode(level, *(f->rootConnection.entryPointHandle), N, n, 1, n+1, level).first;
			v = new WeightedCFLOBDDTopNodeFourier(tempHandle, f->rootConnection.returnMapHandle, fourierSemiring(1, 1));
			return v;
		}

		std::pair<WeightedCFLOBDDTopNodeFourierRefPtr, int> MeasureAndResetTop(unsigned int level, long int n, WeightedCFLOBDDTopNodeFourierRefPtr f, fourierSemiring R)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			auto t = MeasureAndResetNode(level, n, *(f->rootConnection.entryPointHandle), R);
			v = new WeightedCFLOBDDTopNodeFourier(t.first, f->rootConnection.returnMapHandle, fourierSemiring(1, 1));
			return std::make_pair(v, t.second);	
		}

		WeightedCFLOBDDTopNodeFourierRefPtr ResetStateTop(unsigned int level, WeightedCFLOBDDTopNodeFourierRefPtr f)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr v;
			auto t = ResetStateNode(level, *(f->rootConnection.entryPointHandle));
			v = new WeightedCFLOBDDTopNodeFourier(t, f->rootConnection.returnMapHandle, fourierSemiring(1, 1));
			return v;	
		}

	}
}

