#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <cstdarg>
#include <boost/random.hpp>
#include <boost/preprocessor/logical/xor.hpp>
#include "wmatrix1234_complex_fb_mul.h"
#include "wmatrix1234_top_node_complex_fb_mul.h"

namespace CFL_OBDD {


	namespace WeightedMatrix1234ComplexFloatBoostMul {

		void Matrix1234Initializer()
		{
			Matrix1234InitializerTop();
			return;
		}

		// Create representation of identity relation
		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkIdRelationInterleaved(unsigned int i, int cflobdd_kind)
		{
			// TODO: Check for error - "CodeConvert"
			WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr tmp = MkIdRelationInterleavedTop(i, cflobdd_kind);
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(tmp);
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkNegationMatrixInterleaved(unsigned int i, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkNegationMatrixInterleavedTop(i, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkPauliYGate(unsigned int i, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkPauliYGateTop(i, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkPauliZGate(unsigned int i, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkPauliZGateTop(i, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkSGate(unsigned int i, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkSGateTop(i, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkRestrictMatrix(unsigned int level, std::string s, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkRestrictTop(level, s, cflobdd_kind));	
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkWalshInterleaved(unsigned int i, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkWalshInterleavedTop(i, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCNOTInterleaved(unsigned int i)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkCNOTInterleavedTop(i));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkExchangeInterleaved(unsigned int i)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkExchangeInterleavedTop(i));
		}

		// Return the Kronecker product of two matrices
		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL KroneckerProduct2Vocs(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m1, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m2)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(KroneckerProduct2VocsTop(m1.root, m2.root)); 
		}

		// Naive matrix multiplication
		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MatrixMultiplyV4(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m1, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m2)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MatrixMultiplyV4TopNode(m1.root, m2.root));
		}

        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled, int cflobdd_kind){
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkCNOTTopNode(level, n, controller, controlled, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCCNOT(unsigned int level, long int controller1, long int controller2, long int controlled, int cflobdd_kind){
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkCCNOTTop(level, controller1, controller2, controlled, cflobdd_kind));
		}

		void MatrixPrintRowMajor(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c, std::ostream & out)
		{
			MatrixPrintRowMajorTop(c.root, out);
			return;
		}

		void MatrixPrintRowMajorInterleaved(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c, std::ostream & out)
		{
			MatrixPrintRowMajorInterleavedTop(c.root, out);
			return;
		}

        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ApplyExchangeAndIdentity(std::string s){
			if (s.find('1') == std::string::npos){
				return MkIdRelationInterleaved(ceil(log2(s.length())) + 1);
			}
			else if (s.find('0') == std::string::npos){
				return MkExchangeInterleaved(ceil(log2(s.length())) + 1);
			}
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL F1 = ApplyExchangeAndIdentity(s.substr(0, s.length() / 2));
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL F2 = ApplyExchangeAndIdentity(s.substr(s.length() / 2));
			return KroneckerProduct2Vocs(F1, F2);
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CreateBalancedFn(int n, std::mt19937 mt){
			std::string s(2*n, '0');
			for (unsigned int i = 0; i < n; i++)
				s[i] = mt() % 2 ? '1' : '0';
			
			unsigned int level = ceil(log2(n)) + 2;
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL F = ApplyExchangeAndIdentity(s);
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL C = WeightedMatrix1234ComplexFloatBoostMul::MkCNOT(level, n, 0, n);
			for (int i = 1; i < n; i++){
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL tmp = WeightedMatrix1234ComplexFloatBoostMul::MkCNOT(level, n, i, n);
				C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, tmp);
			}
			std::cout << "C created" << std::endl;
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ans = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, F);
			ans = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(F, ans);
			return ans;
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkSwapGate(unsigned int i, long c1, long c2, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkSwapGateTop(i, c1, c2, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkiSwapGate(unsigned int i, long c1, long c2, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkiSwapGateTop(i, c1, c2, cflobdd_kind));
		}
		
		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCPGate(unsigned int i, long c1, long c2, double theta, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkCPGateTop(i, c1, c2, theta, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkPhaseShiftGate(unsigned int i, double theta, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkPhaseShiftGateTop(i, theta, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCZGate(unsigned int i, long c1, long c2, double theta, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkCZGateTop(i, c1, c2, theta, cflobdd_kind));
		}

		WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCSwapGate(unsigned int i, long int c1, long int x1, long int x2, int cflobdd_kind)
		{
			return WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(MkCSwapGateTop(i, c1, x1, x2, cflobdd_kind));
		}
	}
}

