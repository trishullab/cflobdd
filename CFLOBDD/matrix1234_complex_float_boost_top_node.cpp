//
//    Copyright (c) 2017, 2018 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <complex>
#include <math.h>
#include <iomanip>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sin_pi.hpp>

#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "matrix1234_node.h"
#include "matrix1234_complex_float_boost_top_node.h"
#include "assignment.h"

namespace CFL_OBDD {

	namespace Matrix1234ComplexFloatBoost {

		void Matrix1234InitializerTop()
		{
			Matrix1234InitializerNode();
			return;
		}

		// MatrixPrintRowMajorInterleaved
		//
		// Given a CFLOBDD n that represents a 2-vocabulary matrix M_n(W,X) with the
		// vocabulary interleaving (W |x| X), print M_n in row-major order as if M_n were the matrix M(W, X).
		//
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 1 && level <= 5) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				// For each row, column position (i,j) in matrix M_n, set a to the assignment
				// such that evaluating CFLOBDD n w.r.t. a yields the value of M_n[i,j]
				for (unsigned long int i = 0UL; i < rows; i++) {  // Voc1
					// Fill in the even positions of a
					//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits of i
					//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits of i
					unsigned long int mask1 = 1UL;
					for (int k = indexBits - 1; k >= 0; k--) {  // least-significant to most-significant bits
						a[2 * k] = (i & mask1);
						mask1 = mask1 << 1;  // Prepare mask1 for the next-most-significant bit
					}
					for (unsigned long int j = 0UL; j < cols; j++) {  // Voc2
						// Fill in the odd positions of a
						//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits of j
						//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits of j
						unsigned long int mask2 = 1UL;
						for (int k = indexBits - 1; k >= 0; k--) {  // least-significant to most-significant bits
							a[2 * k + 1] = (j & mask2);
							mask2 = mask2 << 1;  // Prepare mask2 for the next-most-significant bit
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;
						std::complex<double> b;
						b = n->EvaluateIteratively(a).convert_to<std::complex<double>>();
						if (b != std::complex<double>(0))
							out << i << " " << j << " " << b << " " << std::endl;
					}
					// out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [1 .. 4]" << std::endl;
			}
		}

		// MatrixPrintRowMajor
		//
		// Given a CFLOBDD n that represents a 4-vocabulary matrix M_n(W,X,Y,Z) with the
		// vocabulary interleaving (W |x| X) |x| (Y |x| Z) = (W |x| Y |x| X |x| Z),
		// print M_n in row-major order as if M_n were the matrix M(W || X, Y || Z).
		//
		// For each entry, the assignment created has Booleans for the bits in root-to-leaf order.
		// The loop that fills the assignment runs least-significant bit to most-significant bit,
		// interleaving row and column bits.  The most-significant bit (of the row value) corresponds
		// to the root; the least-significant bit (of the col value) corresponds to the leaf.
		//
		void MatrixPrintRowMajorTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 2 && level <= 4) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				// For each row, column position (i,j) in matrix M_n, set a to the assignment
				// such that evaluating CFLOBDD n w.r.t. a yields the value of M_n[i,j]
				for (unsigned long int i = 0UL; i < rows; i++) { // Vocs 1 and 3
					// Fill in the even positions of a:
					//  o Voc1 bits are obtained from the top (most-significant) half of i
					//  o Voc3 bits are obtained from the bottom (least-significant) half of i
					// For each vocabulary,
					//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits
					//    of the appropriate half of i
					//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits
					//    of the appropriate half of i
					unsigned long int maskVoc1 = 1UL << (indexBits / 2);
					unsigned long int maskVoc3 = 1UL;
					for (int k = (indexBits / 2 - 1); k >= 0; k--) {  // least-significant to most-significant bits
						a[4 * k] = (i & maskVoc1);
						a[4 * k + 2] = (i & maskVoc3);
						maskVoc1 = maskVoc1 << 1;  // Prepare mask1 for the next-most-significant bit
						maskVoc3 = maskVoc3 << 1;  // Prepare mask3 for the next-most-significant bit
					}
					for (unsigned long int j = 0UL; j < cols; j++) {  // Vocs 2 and 4
						// Fill in the odd positions of a:
						//  o Voc2 bits are obtained from the top (most-significant) half of j
						//  o Voc4 bits are obtained from the bottom (least-significant) half of j
						// For each vocabulary,
						//  o the least-significant bits of a (= a[large]) are obtained from the least-signicant bits
						//    of the appropriate half of j
						//  o the most-significant bits of a (= a[small]) are obtained from the most-significant bits
						//    of the appropriate half of j
						unsigned long int maskVoc2 = 1UL << (indexBits / 2);
						unsigned long int maskVoc4 = 1UL;
						for (int k = (indexBits / 2 - 1); k >= 0; k--) {  // least-significant to most-significant bits
							a[4 * k + 1] = (j & maskVoc2);
							a[4 * k + 3] = (j & maskVoc4);
							maskVoc2 = maskVoc2 << 1;  // Prepare mask2 for the next-most-significant bit
							maskVoc4 = maskVoc4 << 1;  // Prepare mask4 for the next-most-significant bit
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;
						std::complex<double> b;
						b = n->EvaluateIteratively(a).convert_to<std::complex<double>>();
						out << b << " ";
					}
					out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [2 .. 4]" << std::endl;
			}
		}

		// Create representation of identity relation (with interleaved variable order).
		// That is, input (x0,y0,x1,y1,...,xN,yN) yield Id[(x0,x1,...,xN)][(y0,y1,...,yN)]
		// which equals 1 iff xi == yi, for 0 <= i <= N.
		CFLOBDDTopNodeComplexFloatBoostRefPtr MkIdRelationInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m10;

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m10);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDDTopNodeComplexFloatBoostRefPtr MkWalshInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkWalshInterleavedNode(i);
			auto val = boost::multiprecision::pow(sqrt(2), boost::multiprecision::pow(BIG_COMPLEX_FLOAT(2), i-1));
			m.AddToEnd(1/val);
			m.AddToEnd(-1/val);
			// m.AddToEnd(1);
			// m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkPauliYMatrixInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkPauliYInterleavedNode(i);
			BIG_COMPLEX_FLOAT img_i(0, 1);
			m.AddToEnd(0);
			m.AddToEnd(-1 * boost::multiprecision::pow(img_i, i-1));
			m.AddToEnd(boost::multiprecision::pow(img_i, i-1));
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkPauliZMatrixInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkPauliZInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkSGateInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkSGateInterleavedNode(i);
			if (i >= 1){
				m.AddToEnd(1);
				m.AddToEnd(0);
				m.AddToEnd(BIG_COMPLEX_FLOAT(0, 1));
			}
			if (i >= 2)
			{
				m.AddToEnd(-1);
			}
			if (i >= 3)
			{
				m.AddToEnd(BIG_COMPLEX_FLOAT(0, -1));
			}
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkSdgGateInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i == 1);

			tempHandle = MkSGateInterleavedNode(i);
			if (i >= 1){
				m.AddToEnd(1);
				m.AddToEnd(0);
				m.AddToEnd(BIG_COMPLEX_FLOAT(0, -1));
			}
			// if (i >= 2)
			// {
			// 	m.AddToEnd(-1);
			// }
			// if (i >= 3)
			// {
			// 	m.AddToEnd(BIG_COMPLEX_FLOAT(0, -1));
			// }
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkTGateInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i == 1);

			tempHandle = MkSGateInterleavedNode(i);
			if (i >= 1){
				m.AddToEnd(1);
				m.AddToEnd(0);
				double cos_v = SQRT2_2; // boost::math::cos_pi(0.25);
				double sin_v = SQRT2_2; // boost::math::sin_pi(0.25);
				BIG_COMPLEX_FLOAT val(cos_v, sin_v);
				m.AddToEnd(val);
			}
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkTdgGateInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i == 1);

			tempHandle = MkSGateInterleavedNode(i);
			if (i >= 1){
				m.AddToEnd(1);
				m.AddToEnd(0);
				double cos_v = SQRT2_2; // boost::math::cos_pi(-0.25);
				double sin_v = -SQRT2_2; // boost::math::sin_pi(-0.25);
				BIG_COMPLEX_FLOAT val(cos_v, sin_v);
				m.AddToEnd(val);
			}
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkRZGateTop(unsigned int i, double theta)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i == 1);

			double cos_v = boost::math::cos_pi(theta);
			double sin_v = boost::math::sin_pi(theta);
			BIG_COMPLEX_FLOAT val1(cos_v, sin_v);
			BIG_COMPLEX_FLOAT val2(-cos_v, -sin_v);
			tempHandle = MkSGateInterleavedNode(i);
			m.AddToEnd(val1);
			m.AddToEnd(0);
			m.AddToEnd(val2);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkSXGateTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m10;

			BIG_COMPLEX_FLOAT one_plus_i (1, 1);
			BIG_COMPLEX_FLOAT one_minus_i (1, -1);

			tempHandle = MkIdRelationInterleavedNode(i);
			m10.AddToEnd(one_plus_i);
			m10.AddToEnd(one_minus_i);
			m10.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m10);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkSYGateTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m10;

			BIG_COMPLEX_FLOAT one_plus_i (1/2, 1/2);
			BIG_COMPLEX_FLOAT one_minus_i (-1/2, -1/2);

			tempHandle = MkSYGateNode(i);
			m10.AddToEnd(one_plus_i);
			m10.AddToEnd(one_minus_i);
			m10.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m10);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkPhaseShiftGateInterleavedTop(unsigned int i, double theta)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i == 1);

			double cos_v = boost::math::cos_pi(theta);
			double sin_v = boost::math::sin_pi(theta);
			BIG_COMPLEX_FLOAT val(cos_v, sin_v);
			tempHandle = MkSGateInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.AddToEnd(val);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;	
		}

		// Create representation of the Inverse Reed-Muller matrix IRM(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields IRM[(x0,x1,...,xN)][(y0,y1,...,yN)]
		CFLOBDDTopNodeComplexFloatBoostRefPtr MkInverseReedMullerInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkInverseReedMullerInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for two
		// additional vocabularies
		CFLOBDDTopNodeComplexFloatBoostRefPtr MkWalshVoc13Top(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkWalshVoc13Node(i);
			m.AddToEnd(1);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;
		}

		// Create representation of the Walsh matrix W(2**(i-2))
		// [i.e., a matrix of size 2**(2**(i-2))) x 2**(2**(i-2)))]
		// with interleaved indexing of components and room for two
		// additional vocabularies
		CFLOBDDTopNodeComplexFloatBoostRefPtr MkWalshVoc12Top(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkWalshVoc12Node(i);
			m.AddToEnd(1);
			m.AddToEnd(-1);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;
		}


		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVocs13To24Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs13To24Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVocs12To34Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			assert(1 <= n->rootConnection.entryPointHandle->handleContents->level);
			assert(n->rootConnection.entryPointHandle->handleContents->level < CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = PromoteInterleavedTo12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = Demote12ToInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}


		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVoc43Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc43Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Vocabulary shift in the representation of a matrix
		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 2);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixShiftVoc42Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal:
		// (W,X,Y,Z) s.t. X==Y with interleaved variables
		CFLOBDDTopNodeComplexFloatBoostRefPtr MkDetensorConstraintInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkDetensorConstraintInterleavedNode(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;
		}

		// Create representation of a matrix in which vocabularies 1 and 4 are constrained to be equal:
		// (W,X,Y,Z) s.t. W==Z with interleaved variables
		CFLOBDDTopNodeComplexFloatBoostRefPtr MkCFLOBDDMatrixEqVoc14Top(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;

			assert(2 <= i && i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkCFLOBDDMatrixEqVoc14Node(i);
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m);
			return v;
		}

		// Vocabulary projection
		// CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		// {
		// 	assert(n->rootConnection.entryPointHandle.handleContents->level >= 2);

		// 	CFLOBDDLinearMapMemoTableRefPtr memoTable = new CFLOBDDLinearMapMemoTable;

		// 	CFLOBDDTopNodeLinearMapRefPtr temp = MatrixProjectVoc23Node(memoTable, n->rootConnection.entryPointHandle, TopLevelVisit);

		// 	// Begin construction of the appropriate ComplexFloatBoostReturnMapHandle by applying each LinearMap in temp's
		// 	// CFLOBDDLinearMapTupleHandle to n's ComplexFloatBoostReturnMapHandle
		// 	ComplexFloatBoostReturnMapHandle rmh;
		// 	for (unsigned int i = 0; i < temp->rootConnection.returnMapHandle.Size(); i++) {
		// 		rmh.AddToEnd(ApplyLinearMap(temp->rootConnection.returnMapHandle[i], n->rootConnection.returnMapHandle));
		// 	}

		// 	CFLOBDDNodeHandle tempHandle = temp->rootConnection.entryPointHandle;

		// 	// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
		// 	ReductionMapHandle inducedReductionMapHandle;
		// 	ComplexFloatBoostReturnMapHandle inducedReturnMap;
		// 	// TODO: Uncomment this following line
		// 	// rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
		// 	//     CFLOBDDNodeHandle::InitReduceCache();
		// 	CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
		// 	//     CFLOBDDNodeHandle::DisposeOfReduceCache();

		// 	// Create and return CFLOBDDTopNodeComplexFloatBoost
		// 	return(new CFLOBDDTopNodeComplexFloatBoost(reduced_tempHandle, inducedReturnMap));
		// }

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkCNOTInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m01;

			tempHandle = MkCNOTInterleavedNode(i);
			m01.AddToEnd(1);
			m01.AddToEnd(0);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkCPGateTop(unsigned int i, long int c1, long int c2, double theta)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m012;

			double cos_v = boost::math::cos_pi(theta);
			double sin_v = boost::math::sin_pi(theta);
			BIG_COMPLEX_FLOAT val(cos_v, sin_v);

			// std::unordered_map<std::string, CFLOBDDNodeHandle> cp_hashMap;

			tempHandle = MkCPGateNode(i, c1, c2);
			m012.AddToEnd(1);
			m012.AddToEnd(0);
			m012.AddToEnd(val);
			m012.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m012);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkSwapGateTop(unsigned int i, long int c1, long int c2)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m01;

			tempHandle = MkSwapGateNode(i, c1, c2, -1);
			m01.AddToEnd(1);
			m01.AddToEnd(0);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkiSwapGateTop(unsigned int i, long int c1, long int c2)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m012;

			tempHandle = MkiSwapGateNode(i, c1, c2, -1);
			m012.AddToEnd(1);
			m012.AddToEnd(0);
			m012.AddToEnd(BIG_COMPLEX_FLOAT(0, 1));
			m012.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m012);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkCSwapGateTop(unsigned int level, long int c, long int i, long int j)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m01;
			tempHandle = MkCSwapGate2Node(level, c, i, j, -1);
			m01.AddToEnd(1);
			m01.AddToEnd(0);
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkRestrictTop(unsigned int level, std::string s)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m01;
			auto tmp = MkRestrictNode(level, s);
			if (tmp.second == 0)
			{
				m01.AddToEnd(0);
				m01.AddToEnd(1);
			}
			else if (tmp.second == 1)
			{
				m01.AddToEnd(1);
				m01.AddToEnd(0);
			}
			else if (tmp.second == -1)
			{
				m01.AddToEnd(1);
			}
			m01.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tmp.first, m01);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToAConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;

			tempHandle = MatrixShiftToBConnectionNode(*(c->rootConnection.entryPointHandle));
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, c->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkFourierDiagonalComponentTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m;
			BIG_COMPLEX_FLOAT one(1);

			assert(2 <= i && i <= CFLOBDDTopNodeComplexFloatBoost::maxLevel);

			tempHandle = MkFourierDiagonalComponentNode(i);
			unsigned long int n = 1UL << (1 << (i - 2));
			unsigned long int log_n = 1 << (i - 2);
			unsigned long int nSq = 1UL << (1 << (i - 1));

			// Create a listing of the diagonal entries for a matrix with the interleaved order
			// (Voc1 || Voc4) |x| (Voc2 || Voc3)
			BIG_COMPLEX_FLOAT *diagonalEntries = new (std::nothrow) BIG_COMPLEX_FLOAT[nSq];
			if (diagonalEntries == nullptr) {
				std::cout << "memory allocation failed in MkFourierDiagonalComponentNode" << std::endl;
				abort();
			}
			double temp = static_cast<double>(n)* static_cast<double>(n);
			double temp2 = 1.0 / temp;
			for (unsigned long int j = 0; j < n; j++) {
				for (unsigned long int k = 0; k < n; k++) {
					//std::complex<double> w_jk = std::pow(w, ((j*k) % nSq));   // w ** ((j*k) mod nSq)
					double v = (-2 * ((j*k) % nSq) / static_cast<double>(nSq));
					double cos_v = boost::math::cos_pi(v);
					double sin_v = boost::math::sin_pi(v);
					BIG_COMPLEX_FLOAT w_jk(cos_v, sin_v);
					diagonalEntries[j*n + k] = w_jk;
				}
			}
			// Create ReturnMapHandle m in a way that respects interleaved order
			// Voc1 |x| Voc2 |x| Voc3 |x| Voc4.
			boost::unordered_map<BIG_COMPLEX_FLOAT, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned long int j = 0; j < nSq; j++) {
				// Set adjustedIndex to the unshuffle of the bits of j
				unsigned long int adjustedIndex = 0;
				unsigned long int maskVoc34 = 1UL;
				unsigned long int maskVoc12 = 2UL;
				for (unsigned long int k = 0; k < log_n; k++) {
					if (j & maskVoc34)
						adjustedIndex += (1UL << k);
					if (j & maskVoc12)
						adjustedIndex += (1UL << (log_n + k));
					maskVoc34 = maskVoc34 << 2;
					maskVoc12 = maskVoc12 << 2;
				}
				auto it = reductionMap.find(diagonalEntries[adjustedIndex]);
				if (it == reductionMap.end()){
					m.AddToEnd(diagonalEntries[adjustedIndex]);
					reductionMap.emplace(diagonalEntries[adjustedIndex], m.Size() - 1);
					reductionMapHandle.AddToEnd(m.Size() - 1);
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[diagonalEntries[adjustedIndex]]);
				}
				if (j == 0){
					BIG_COMPLEX_FLOAT zero(0.0, 0.0);
					m.AddToEnd(zero);
					reductionMap.emplace(zero, m.Size() - 1);
					reductionMapHandle.AddToEnd(m.Size() - 1);
				}
			}
			delete[] diagonalEntries;
			m.Canonicalize();
			reductionMapHandle.Canonicalize();

			// Perform reduction on tempHandle, with respect to the common elements that m maps together
			//ReductionMapHandle inducedReductionMapHandle;
			//ComplexFloatBoostReturnMapHandle inducedReturnMap;
			//m.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//     CFLOBDDNodeHandle::InitReduceCache();
			//CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, m.Size(), true);
			//     CFLOBDDNodeHandle::DisposeOfReduceCache();

			// Create and return CFLOBDDTopNodeComplexFloatBoost
			//return(new CFLOBDDTopNodeComplexFloatBoost(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeComplexFloatBoost(reduced_tempHandle, m));
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr ConvertToComplexTop(CFLOBDDTopNodeFourierRefPtr c)
		{
			ComplexFloatBoostReturnMapHandle m;
			for (int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				double v = (((2 * c->rootConnection.returnMapHandle[i].GetVal()) % (c->rootConnection.returnMapHandle[i].GetRingSize())).convert_to<double>() 
					/ static_cast<double>(c->rootConnection.returnMapHandle[i].GetRingSize().convert_to<double>()));
				double cos_v = boost::math::cos_pi(v);
				double sin_v = boost::math::sin_pi(v);
				BIG_COMPLEX_FLOAT w_jk(cos_v, sin_v);
				m.AddToEnd(w_jk);
			}
			m.Canonicalize();
			return (new CFLOBDDTopNodeComplexFloatBoost(*(c->rootConnection.entryPointHandle), m));
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkExchangeInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m01;

			tempHandle = MkExchangeInterleavedNode(i);
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkNegationMatrixInterleavedTop(unsigned int i)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr v;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle m01;

			tempHandle = MkNegationMatrixInterleavedNode(i);
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, m01);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled)
		{
			CFLOBDDNodeHandle c = MkCNOT2Node(level, n, controller, controlled);
			ComplexFloatBoostReturnMapHandle m;
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			return new CFLOBDDTopNodeComplexFloatBoost(c, m);
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkCCNOTTopNode(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled)
		{
			CFLOBDDNodeHandle c = MkCCNOTNode(level, n, controller1, controller2, controlled);
			ComplexFloatBoostReturnMapHandle m;
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			return new CFLOBDDTopNodeComplexFloatBoost(c, m);
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkMCXTopNode(unsigned int level, unsigned int n, std::vector<long int>& controllers, long int controlled)
		{
			CFLOBDDNodeHandle c = MkMCXNode(level, n, controllers, controlled);
			ComplexFloatBoostReturnMapHandle m;
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.Canonicalize();
			return new CFLOBDDTopNodeComplexFloatBoost(c, m);
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkCCPTopNode(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled, double theta)
		{
			
			double cos_v = boost::math::cos_pi(theta);
			double sin_v = boost::math::sin_pi(theta);
			BIG_COMPLEX_FLOAT val(cos_v, sin_v);

			CFLOBDDNodeHandle c = MkCCPNode(level, n, controller1, controller2, controlled);
			ComplexFloatBoostReturnMapHandle m;
			m.AddToEnd(1);
			m.AddToEnd(0);
			m.AddToEnd(val);
			m.Canonicalize();
			return new CFLOBDDTopNodeComplexFloatBoost(c, m);
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixTransposeTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			std::unordered_map<CFLOBDDNodeHandle, std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>, 
				CFLOBDDNodeHandle::CFLOBDDNodeHandle_Hash> hashMap;
			auto pc = MatrixTransposeNode(hashMap, *(n->rootConnection.entryPointHandle));
			CFLOBDDNodeHandle temp = pc.first;
			ReductionMapHandle reductionMapHandle;
			ComplexFloatBoostReturnMapHandle v;
			for (unsigned int i = 0; i < pc.second.Size(); i++){
				reductionMapHandle.AddToEnd(i);
				v.AddToEnd(n->rootConnection.returnMapHandle[pc.second[i]]);
			}
			reductionMapHandle.Canonicalize();
			v.Canonicalize();
			temp = temp.Reduce(reductionMapHandle, v.Size(), true);
			return new CFLOBDDTopNodeComplexFloatBoost(temp, v);
			return n;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr ConjugateTransposeTop(CFLOBDDTopNodeComplexFloatBoostRefPtr c)
		{
			// std::unordered_map<CFLOBDDNodeHandle, std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>, 
			// 	CFLOBDDNodeHandle::CFLOBDDNodeHandle_Hash> hashMap;
			// auto pc = MatrixTransposeNode(hashMap, *(c->rootConnection.entryPointHandle));
			// CFLOBDDNodeHandle temp = pc.first;
			// ReductionMapHandle reductionMapHandle;
			// ComplexFloatBoostReturnMapHandle v;
			// for (unsigned int i = 0; i < pc.second.Size(); i++){
			// 	reductionMapHandle.AddToEnd(i);
			// 	v.AddToEnd(c->rootConnection.returnMapHandle[pc.second[i]]);
			// }
			// reductionMapHandle.Canonicalize();
			// v.Canonicalize();
			// temp = temp.Reduce(reductionMapHandle, v.Size(), true);
			// return new CFLOBDDTopNodeComplexFloatBoost(temp, v);	
			// return c;
			ComplexFloatBoostReturnMapHandle m;
			for (int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				auto val = c->rootConnection.returnMapHandle[i];
				m.AddToEnd(BIG_COMPLEX_FLOAT(val.real(), -val.imag()));
			}
			m.Canonicalize();
			return new CFLOBDDTopNodeComplexFloatBoost(*(c->rootConnection.entryPointHandle), m);
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeComplexFloatBoostRefPtr c1, CFLOBDDTopNodeComplexFloatBoostRefPtr c2)
		{
			std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash> hashMap;
			// if (c1->level >= 5)
			// 	clearMultMap();
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4Node(hashMap, *(c1->rootConnection.entryPointHandle),*(c2->rootConnection.entryPointHandle));
			ComplexFloatBoostReturnMapHandle v;
			// boost::unordered_map<std::complex<double>, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				BIG_COMPLEX_FLOAT val = 0;
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					auto factor = j.second.convert_to<unsigned long long int>();
					val = val + (factor * (c1->rootConnection.returnMapHandle[index1] * c2->rootConnection.returnMapHandle[index2]));
				}
				unsigned int k = v.Size();
				for (k = 0; k < v.Size(); k++)
				{
					if (v[k] == val)
					{
						break;
					}
				}
				if (k < v.Size())
				{
					reductionMapHandle.AddToEnd(k);
				}
				else
				{
					v.AddToEnd(val);
					reductionMapHandle.AddToEnd(v.Size()-1);
				}
				// if (reductionMap.find(val) == reductionMap.end()){
				// 	v.AddToEnd(BIG_COMPLEX_FLOAT(val));
				// 	reductionMap.insert(std::make_pair(val, v.Size() - 1));
				// 	reductionMapHandle.AddToEnd(v.Size() - 1);
				// }
				// else{
				// 	reductionMapHandle.AddToEnd(reductionMap[val]);
				// }
			}

			v.Canonicalize();
			reductionMapHandle.Canonicalize();
			CFLOBDDNodeHandle tempHandle = *(c->rootConnection.entryPointHandle);
			// Perform reduction on tempHandle, with respect to the common elements that rmh maps together
			/*ReductionMapHandle inducedReductionMapHandle;
			FloatBoostReturnMapHandle inducedReturnMap;
			v.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());*/
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			// Create and return CFLOBDDTopNode
			//return(new CFLOBDDTopNodeFloatBoost(reduced_tempHandle, inducedReturnMap));
			return(new CFLOBDDTopNodeComplexFloatBoost(reduced_tempHandle, v));
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixMultiplyV4WithInfoTopNode(CFLOBDDTopNodeComplexFloatBoostRefPtr c1, CFLOBDDTopNodeComplexFloatBoostRefPtr c2)
		{
			// if (c1->level >= 5)
			// 	clearMultMap();
			const BIG_COMPLEX_FLOAT TOLERANCE_LEVEL = BIG_COMPLEX_FLOAT(1e10, 1e10);
			std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash> hashMap;
			int c1_zero_index = -1, c2_zero_index = -1;
			c1_zero_index = c1->rootConnection.returnMapHandle.LookupInv(0);
			c2_zero_index = c2->rootConnection.returnMapHandle.LookupInv(0);
			CFLOBDDTopNodeMatMultMapRefPtr c = MatrixMultiplyV4WithInfoNode
				(hashMap, *(c1->rootConnection.entryPointHandle), *(c2->rootConnection.entryPointHandle),
				c1_zero_index, c2_zero_index);
			ComplexFloatBoostReturnMapHandle v;
			boost::unordered_map<BIG_COMPLEX_FLOAT, unsigned int> reductionMap;
			ReductionMapHandle reductionMapHandle;
			std::cout << std::setprecision(15);
			for (unsigned int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++){
				MatMultMapHandle r = c->rootConnection.returnMapHandle[i];
				BIG_COMPLEX_FLOAT val = 0;
				for (auto &j : r.mapContents->map){
					unsigned int index1 = j.first.first;
					unsigned int index2 = j.first.second;
					if (index1 != -1 && index2 != -1){
						auto factor = BIG_COMPLEX_FLOAT(j.second, 0);
						BIG_COMPLEX_FLOAT previous_val = val;
						BIG_COMPLEX_FLOAT tmp_val = mul(c1->rootConnection.returnMapHandle[index1], c2->rootConnection.returnMapHandle[index2]);
						BIG_COMPLEX_FLOAT added_val = mul(factor, tmp_val);
						val = add(previous_val, added_val);
					}
				}
				val = roundNearBy(val);
				
				if (reductionMap.find(val) == reductionMap.end()){
					v.AddToEnd(val);
					reductionMap.insert(std::make_pair(val, v.Size() - 1));
					reductionMapHandle.AddToEnd(v.Size() - 1);
					// std::cout << "New Value: " << val << " at " << v.Size() - 1 << std::endl;
				}
				else{
					reductionMapHandle.AddToEnd(reductionMap[val]);
					// std::cout << "Existing Value: " << val << " at " << reductionMap[val] << std::endl;
				}
			}
			v.Canonicalize();
			reductionMapHandle.Canonicalize();
			CFLOBDDNodeHandle tempHandle = *(c->rootConnection.entryPointHandle);
			CFLOBDDNodeHandle reduced_tempHandle = tempHandle.Reduce(reductionMapHandle, v.Size(), true);
			return(new CFLOBDDTopNodeComplexFloatBoost(reduced_tempHandle, v));
		}
	}
}

