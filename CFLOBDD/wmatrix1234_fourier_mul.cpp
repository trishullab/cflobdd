#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <cstdarg>
#include <boost/random.hpp>
#include <boost/preprocessor/logical/xor.hpp>
#include "fourier_semiring.h"
#include "wmatrix1234_fourier_mul.h"
#include "wmatrix1234_top_node_fourier_mul.h"

namespace CFL_OBDD {


	namespace WeightedMatrix1234FourierMul {

		void Matrix1234Initializer()
		{
			Matrix1234InitializerTop();
			return;
		}

		// Create representation of identity relation
		WEIGHTED_CFLOBDD_FOURIER_MUL MkIdRelationInterleaved(unsigned int i)
		{
			// TODO: Check for error - "CodeConvert"
			WeightedCFLOBDDTopNodeFourierRefPtr tmp = MkIdRelationInterleavedTop(i);
			return WEIGHTED_CFLOBDD_FOURIER_MUL(tmp);
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkNegationMatrixInterleaved(unsigned int i)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkNegationMatrixInterleavedTop(i));
		}

		// Create representation of the Walsh matrix W(2**(i-1))
		// [i.e., a matrix of size 2**(2**(i-1))) x 2**(2**(i-1)))]
		// with interleaved indexing of components: that is, input
		// (x0,y0,x1,y1,...,xN,yN) yields W[(x0,x1,...,xN)][(y0,y1,...,yN)]
		WEIGHTED_CFLOBDD_FOURIER_MUL MkWalshInterleaved(unsigned int i)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkWalshInterleavedTop(i));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkCNOTInterleaved(unsigned int i)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkCNOTInterleavedTop(i));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkExchangeInterleaved(unsigned int i)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkExchangeInterleavedTop(i));
		}

		// Return the Kronecker product of two matrices
		WEIGHTED_CFLOBDD_FOURIER_MUL KroneckerProduct2Vocs(WEIGHTED_CFLOBDD_FOURIER_MUL m1, WEIGHTED_CFLOBDD_FOURIER_MUL m2)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(KroneckerProduct2VocsTop(m1.root, m2.root)); 
		}

		// Naive matrix multiplication
		WEIGHTED_CFLOBDD_FOURIER_MUL MatrixMultiplyV4(WEIGHTED_CFLOBDD_FOURIER_MUL m1, WEIGHTED_CFLOBDD_FOURIER_MUL m2)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MatrixMultiplyV4TopNode(m1.root, m2.root));
		}

        WEIGHTED_CFLOBDD_FOURIER_MUL MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled){
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkCNOTTopNode(level, n, controller, controlled));
		}

		void MatrixPrintRowMajor(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out)
		{
			MatrixPrintRowMajorTop(c.root, out);
			return;
		}

		void MatrixPrintRowMajorInterleaved(WEIGHTED_CFLOBDD_FOURIER_MUL c, std::ostream & out)
		{
			MatrixPrintRowMajorInterleavedTop(c.root, out);
			return;
		}

        WEIGHTED_CFLOBDD_FOURIER_MUL ApplyExchangeAndIdentity(std::string s){
			if (s.find('1') == std::string::npos){
				return MkIdRelationInterleaved(ceil(log2(s.length())) + 1);
			}
			else if (s.find('0') == std::string::npos){
				return MkExchangeInterleaved(ceil(log2(s.length())) + 1);
			}
			WEIGHTED_CFLOBDD_FOURIER_MUL F1 = ApplyExchangeAndIdentity(s.substr(0, s.length() / 2));
			WEIGHTED_CFLOBDD_FOURIER_MUL F2 = ApplyExchangeAndIdentity(s.substr(s.length() / 2));
			return KroneckerProduct2Vocs(F1, F2);
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL CreateBalancedFn(int n, std::mt19937 mt){
			std::string s(2*n, '0');
			for (unsigned int i = 0; i < n; i++)
				s[i] = mt() % 2 ? '1' : '0';
			
			unsigned int level = ceil(log2(n)) + 2;
			WEIGHTED_CFLOBDD_FOURIER_MUL F = ApplyExchangeAndIdentity(s);
			WEIGHTED_CFLOBDD_FOURIER_MUL C = WeightedMatrix1234FourierMul::MkCNOT(level, n, 0, n);
			for (int i = 1; i < n; i++){
				WEIGHTED_CFLOBDD_FOURIER_MUL tmp = WeightedMatrix1234FourierMul::MkCNOT(level, n, i, n);
				C = WeightedMatrix1234FourierMul::MatrixMultiplyV4(C, tmp);
			}
			std::cout << "C created" << std::endl;
			WEIGHTED_CFLOBDD_FOURIER_MUL ans = WeightedMatrix1234FourierMul::MatrixMultiplyV4(C, F);
			ans = WeightedMatrix1234FourierMul::MatrixMultiplyV4(F, ans);
			return ans;
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkSwapGate(unsigned int i, long c1, long c2)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkSwapGateTop(i, c1, c2));
		}
		
		WEIGHTED_CFLOBDD_FOURIER_MUL MkCPGate(unsigned int i, long c1, long c2, fourierSemiring theta)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkCPGateTop(i, c1, c2, theta));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkCSwapGate(unsigned int i, long int c1, long int x1, long int x2)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkCSwapGateTop(i, c1, x1, x2));
		}
		
		WEIGHTED_CFLOBDD_FOURIER_MUL MkRZGate(unsigned int level, fourierSemiring theta)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkRZGateTop(level, theta));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkCADDGate(unsigned int level, int c, int x, WEIGHTED_CFLOBDD_FOURIER_MUL f)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkCADDGateTop(level, c, x, f.root));
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkCADDGate2(unsigned int level, int x, WEIGHTED_CFLOBDD_FOURIER_MUL f)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkCADDGate2Top(level, x, f.root));
		}

		bool CheckIfIndexIsNonZero(unsigned int level, int index, WEIGHTED_CFLOBDD_FOURIER_MUL f)
		{
			return CheckIfIndexIsNonZeroTop(level, index, f.root);
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL MkSetBToZero(unsigned int level, WEIGHTED_CFLOBDD_FOURIER_MUL f)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(MkSetBToZeroTop(level, f.root)); 
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL ComputeIQFT(unsigned int level, WEIGHTED_CFLOBDD_FOURIER_MUL f, BIG_INT N, int n)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(ComputeIQFTTop(level, f.root, N, n)); 
		}

		std::pair<WEIGHTED_CFLOBDD_FOURIER_MUL, int> MeasureAndReset(unsigned int level, long int n, WEIGHTED_CFLOBDD_FOURIER_MUL f, fourierSemiring R)
		{
			auto t = MeasureAndResetTop(level, n, f.root, R);
			return std::make_pair(WEIGHTED_CFLOBDD_FOURIER_MUL(t.first), t.second);
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL ResetState(unsigned int level, WEIGHTED_CFLOBDD_FOURIER_MUL f)
		{
			return WEIGHTED_CFLOBDD_FOURIER_MUL(ResetStateTop(level, f.root));
		}
	}
}

