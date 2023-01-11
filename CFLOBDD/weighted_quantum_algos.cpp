#include "weighted_quantum_algos.h"
#include "wvector_fb_mul.h"
#include <ctime>
#include <time.h>
#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>
#include "wmatrix1234_complex_fb_mul.h"
#include "wvector_complex_fb_mul.h"
#include "wvector_fourier_mul.h"
#include "wmatrix1234_fourier_mul.h"
using namespace std;
using namespace chrono;

namespace CFL_OBDD
{
    namespace WeightedQuantumAlgos
    {
        namespace mp = boost::multiprecision;

        std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> GHZ(unsigned long long int n)
        {
            int level = ceil(log2(n)) + 2;
            auto s1 = high_resolution_clock::now();
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F = WeightedMatrix1234FloatBoostMul::MkCNOT(level, 2*n, 0, n);
			std::cout << "Starting loop" << std::endl;
			for (unsigned int i = 1; i < n; i++){
                // auto x = high_resolution_clock::now(); 
				WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL tmp = WeightedMatrix1234FloatBoostMul::MkCNOT(level, 2*n, i, n);
                auto y = high_resolution_clock::now(); 
				F = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(F, tmp);
                auto z = high_resolution_clock::now(); 
                // auto dx = duration_cast<milliseconds>(y - x); 
                auto dy = duration_cast<milliseconds>(z - y); 
                // std::cout << i << " " << dy.count() << std::endl;
			}
            // std::cout << F << std::endl;
            std::cout << "start" << std::endl;
            auto s2 = high_resolution_clock::now();
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL e0 = WeightedVectorFloatBoostMul::NoDistinctionNode(level - 2, 1);
            e0 = WeightedVectorFloatBoostMul::VectorToMatrixInterleaved(e0);
			std::string last_one(2*n, '0');
			last_one[0] = '1';
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL e1 = WeightedVectorFloatBoostMul::MkBasisVector(level - 1, last_one);
            // e1 = WeightedVectorFloatBoostMul::VectorToMatrixInterleaved(e1);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL stateV = WeightedMatrix1234FloatBoostMul::KroneckerProduct2Vocs(e0, e1);
            auto s3 = high_resolution_clock::now();
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL H2 = WeightedMatrix1234FloatBoostMul::MkWalshInterleaved(level);
			stateV = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(F, stateV);
            auto s4 = high_resolution_clock::now();
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 2 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(H2, stateV);
            // auto val = boost::multiprecision::pow(BIG_FLOAT(sqrt(2)), 3*n).convert_to<BIG_FLOAT>();
            // stateV = val * stateV;
            // stateV.print(std::cout);
            auto s5 = high_resolution_clock::now();
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			// stateV = VectorFloatBoost::VectorWithAmplitude(stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 4 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV.ComputeWeightOfPathsAsAmpsToExits();
			std::string ans_s = "";
            ans_s = WeightedVectorFloatBoostMul::Sampling(stateV, true).substr(0, n + 1);
            auto s6 = high_resolution_clock::now();
            auto d1 = duration_cast<milliseconds>(s2 - s1);
            auto d2 = duration_cast<milliseconds>(s3 - s2);
            auto d3 = duration_cast<milliseconds>(s4 - s3);
            auto d4 = duration_cast<milliseconds>(s5 - s4);
            auto d5 = duration_cast<milliseconds>(s6 - s5);
            std::cout << d1.count() << " " << d2.count() << " " << d3.count() << " " << d4.count() << " " << d5.count() << std::endl;
			return std::make_pair(ans_s, stateV);
        }

        std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> BV(long long int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F){
			int level = ceil(log2(n));
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL e = WeightedVectorFloatBoostMul::NoDistinctionNode(level, 1);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL e0 = WeightedVectorFloatBoostMul::MkBasisVector(level + 1, 0);
			e = WeightedVectorFloatBoostMul::VectorToMatrixInterleaved(e);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL stateV = WeightedMatrix1234FloatBoostMul::KroneckerProduct2Vocs(e, e0);
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 0 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(F, stateV);
			// BIG_FLOAT coeff = mp::pow(BIG_FLOAT(2), -n / 2);
			// stateV = coeff * stateV;
            // stateV.print(std::cout);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL H = WeightedMatrix1234FloatBoostMul::MkWalshInterleaved(level + 1);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL I = WeightedMatrix1234FloatBoostMul::MkIdRelationInterleaved(level + 1);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL HI = WeightedMatrix1234FloatBoostMul::KroneckerProduct2Vocs(H, I);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(HI, stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
		    //std::cout << "Step 2 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			// stateV = VectorFloatBoost::VectorWithAmplitude(stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV.ComputeWeightOfPathsAsAmpsToExits();
            // stateV.print(std::cout);
			std::string ans_s = "";
			while (ans_s.find('1') == std::string::npos){
				ans_s = WeightedVectorFloatBoostMul::Sampling(stateV, true).substr(0, n);
			}
			return std::make_pair(ans_s, stateV);
		}


        std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> DeutschJozsaAlgo(unsigned int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F)
		{
			int level = ceil(log2(n));
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL X = WeightedVectorFloatBoostMul::NoDistinctionNode(level, 1);
			std::string last_one(2*n, '0');
			last_one[0] = '1';
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL Y = WeightedVectorFloatBoostMul::MkBasisVector(level+1, last_one);
			// Y = WeightedVectorFloatBoostMul::VectorToMatrixInterleaved(Y);
			X = WeightedVectorFloatBoostMul::VectorToMatrixInterleaved(X);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL H = WeightedMatrix1234FloatBoostMul::MkWalshInterleaved(level + 1);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL I = WeightedMatrix1234FloatBoostMul::MkIdRelationInterleaved(level + 1);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL HI = WeightedMatrix1234FloatBoostMul::KroneckerProduct2Vocs(H, I);
			Y = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(H, Y);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL ans = WeightedMatrix1234FloatBoostMul::KroneckerProduct2Vocs(X, Y);
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(F, ans);
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 2 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(HI, ans);
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 4 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans.ComputeWeightOfPathsAsAmpsToExits();
			std::string ans_s = WeightedVectorFloatBoostMul::Sampling(ans, true).substr(0, n);
			return std::make_pair(ans_s, ans);
		}

        unsigned int getIndexFromBitString(std::string s){
			unsigned int index = 0;
			for (int i = 0; i < s.length(); i++){
				index = index * 2 + ((s[i] == '0') ? 0 : 1);
				index = index * 2 + ((s[i] == '0') ? 0 : 1);
			}
			return index;
		}

		std::string getIndexStringFromBitString(std::string s){
			std::string ans(s.length() * 2, '0');
			unsigned int j = 0;
			for (int i = 0; i < s.length(); i++){
				ans[j] = s[i];
				ans[j + 1] = s[i];
				j+=2;
			}
			return ans;
		}

        WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL MkU_w(std::string s, int n){
			std::string index = getIndexStringFromBitString(s);
			unsigned int level = ceil(log2(n));
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL I = WeightedMatrix1234FloatBoostMul::MkIdRelationInterleaved(level + 1);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F = WeightedVectorFloatBoostMul::MkBasisVector(level + 1, index);
            BIG_FLOAT c = -2;
			F = c * F;
			return I + F;
		}

		WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL MkU_s(int n){
			unsigned int level = ceil(log2(n));
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL I = WeightedMatrix1234FloatBoostMul::MkIdRelationInterleaved(level + 1);
			BIG_FLOAT val = BIG_FLOAT(1.0) / mp::pow(BIG_FLOAT(2), n - 1);
            BIG_FLOAT c = -1;
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL ans = val * WeightedVectorFloatBoostMul::NoDistinctionNode(level+1, 1) + (c *I);
			return ans;
		}

        std::pair<WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL>
			MultiplyRec(WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL M, mp::cpp_int iters,
			boost::unordered_map<mp::cpp_int, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL>& memo, unsigned int n, unsigned int* level){
			//std::cout << iters << std::endl;
			auto it = memo.find(iters);
			if (it != memo.end())
				return std::make_pair(it->second, it->second);
			if (iters == 1){
				BIG_FLOAT val = BIG_FLOAT(1.0);// / mp::pow(BIG_FLOAT(2), n / 4);
				WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL ans = M;
				if ((*level) >= n) {
					memo[iters] = ans;
					//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
					//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
					//std::cout << "Step " << iters << " : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
				}
				if ((*level) < n){
					*level = *level + 1;
					return std::make_pair((BIG_FLOAT(1.0) / sqrt(2)) * ans, ans);
				}
				return std::make_pair(ans, ans);
			}
			mp::cpp_int half_iters = iters / 2;
			auto F1 = MultiplyRec(M, half_iters, memo, n, level);
			auto F2 = MultiplyRec(M, iters - half_iters, memo, n, level);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL ans = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(F1.first, F2.first);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL true_ans = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(F1.second, F2.second);
			if ((*level) >= n) {
				memo[iters] = true_ans;
				//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
				//true_ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
				//std::cout << "Step " << iters << " : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
				//std::cout << "iter: " << iters << " return: " << true_ans.root->rootConnection.returnMapHandle << std::endl;
			}
			/*if ((*level) < n) {
				*level = *level + 1;
				return std::make_pair((BIG_FLOAT(1.0) / sqrt(2)) * ans, true_ans);
			}*/
			return std::make_pair(ans, true_ans);
		}


        std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> GroversAlgoWithV4(int n, std::string s){
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL U_w = MkU_w(s, n);
			std::cout << "U_w created" << std::endl;
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL U_s = MkU_s(n);
			std::cout << "U_s created" << std::endl;
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL M = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(U_s, U_w);
			std::cout << "M created" << std::endl;
			std::cout << "M leaves count: " << M.root->rootConnection.returnMapHandle.Size() << std::endl;
			std::cout << "M leaves: " << M.root->rootConnection.returnMapHandle << std::endl;
			
			double const pi = 4 * std::atan(1);
			BIG_FLOAT tmp_val(mp::pow(BIG_FLOAT(2.0), n / 2) * boost::math::constants::pi<BIG_FLOAT>() * 0.25);
			mp::cpp_int iters = mp::floor(tmp_val).convert_to<mp::cpp_int>();
			unsigned int level = ceil(log2(n));
			BIG_FLOAT val = BIG_FLOAT(1.0);
			if (n - iters > 0)
				val = val / mp::pow(BIG_FLOAT(2), (BIG_FLOAT(n) - BIG_FLOAT(iters)) / 2);
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL ans = val * WeightedVectorFloatBoostMul::NoDistinctionNode(level + 1, 1);
			std::cout << "Iter start, num iters: " << iters << std::endl;
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			boost::unordered_map<mp::cpp_int, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> memo;
			unsigned int ulevel = 0;
			// auto M_rec = MultiplyRec(M, iters, memo, n, &ulevel);
			// WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL M_actual = M_rec.first;
			// std::cout << ulevel << std::endl;
			// std::cout << M_actual.root->rootConnection.returnMapHandle << std::endl;
			// ans = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(M_actual, ans);

            for (mp::cpp_int i = 0; i < iters; i++)
            {
                ans = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(M, ans);
                std::cout << "i: " << i << std::endl;
            }

			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step i : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;

			std::cout << "Iter end" << std::endl;
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step j : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans.ComputeWeightOfPathsAsAmpsToExits();
            // std::cout << ans << std::endl;
			std::string ans_s = WeightedVectorFloatBoostMul::Sampling(ans, true, "Grovers").substr(0, n);
			return std::make_pair(ans_s, ans);
		}

        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL Hadamard(unsigned int n, unsigned int i)
		{
			if (n == 1)
			{
				return WeightedMatrix1234ComplexFloatBoostMul::MkWalshInterleaved(1);
			}
			else {
                int level = ceil(log2(n/2));
				if (i < n/2)
				{
					WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL T = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(level + 1);
					WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = Hadamard(n/2, i);
					return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, T);
				}
				else
				{
					WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL T = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(level + 1);
					return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(T, Hadamard(n/2, i - n/2)); 
				}
			}
		}

        WEIGHTED_CFLOBDD_FOURIER_MUL Hadamard_i(unsigned int n, unsigned int i)
		{
			if (n == 1)
			{
				return WeightedMatrix1234FourierMul::MkWalshInterleaved(1);
			}
			else {
                int level = ceil(log2(n/2));
				if (i < n/2)
				{
					WEIGHTED_CFLOBDD_FOURIER_MUL T = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level + 1);
					WEIGHTED_CFLOBDD_FOURIER_MUL H = Hadamard_i(n/2, i);
					return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(H, T);
				}
				else
				{
					WEIGHTED_CFLOBDD_FOURIER_MUL T = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level + 1);
					return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(T, Hadamard_i(n/2, i - n/2)); 
				}
			}
		}

        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL QFT(long long int n, std::string s)
		{
			unsigned int level = ceil(log2(n));
			
			// std::cout << "s: " << s << std::endl;
			// std::reverse(s.begin(), s.end());
            std::string S(2*s.length(), '0');
            for (int i = 0; i < s.length(); i++)
            {
                S[2*i] = s[i];
            }
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL stateV = WeightedVectorComplexFloatBoostMul::MkBasisVector(level+1, S);
			// stateV = WeightedVectorComplexFloatBoostMul::VectorToMatrixInterleaved(stateV);
			std::cout << "start" << std::endl;


			for (long long int i = 0; i < n/2; i++)
			{
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL SwapM = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(level+1, i, n-i-1);
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(SwapM, stateV);
			}

			std::cout << "loop start" << std::endl;

			for (long long int i = n-1; i >= 0; i--)
			{
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = Hadamard(n, i);
                // BIG_COMPLEX_FLOAT c = 1.0/sqrt(2);
                // H = c * H;
                // auto start = high_resolution_clock::now();
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateV);
                // auto end = high_resolution_clock::now();
                // auto duration = duration_cast<milliseconds>(end - start);
                // std::cout << "(i): " << i << " " << duration.count() << std::endl;
				for (long int j = 0; j < i; j++)
				{
					double theta = std::pow(2, j - i);
					WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CP = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(level+1, j, i, theta);
                    // auto start = high_resolution_clock::now();
					stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CP, stateV);
                    // std::cout << "(i,j): " << i << "," << j << std::endl;
                    // auto end = high_resolution_clock::now();
                    // auto duration = duration_cast<milliseconds>(end - start);
                    // std::cout << "(i,j): " << i << "," << j << " " << duration.count() << std::endl;
				}
			}

            /*
                Inverse QFT

            for (long long int i = 0; i < n; i++)
			{
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = Hadamard(n, i);
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateV);
				for (long int j = i+1; j < n; j++)
				{
					double theta = -1 * std::pow(2, i - j);
					WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CP = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(level+1, i, j, theta);
					stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CP, stateV);
				}
			}

            for (long long int i = 0; i < n/2; i++)
			{
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL SwapM = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(level+1, i, n-i-1);
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(SwapM, stateV);
			}
            */

			std::cout << "done" << std::endl;


			return stateV;
		}

        WEIGHTED_CFLOBDD_FOURIER_MUL QFT_fourier(long long int n, std::string s)
		{
			unsigned int level = ceil(log2(n));
			
			// std::cout << "s: " << s << std::endl;
			// std::reverse(s.begin(), s.end());
            std::string S(2*s.length(), '0');
            for (int i = 0; i < s.length(); i++)
            {
                S[2*i] = s[i];
            }
			WEIGHTED_CFLOBDD_FOURIER_MUL stateV = WeightedVectorFourierMul::MkBasisVector(level+1, S);
			std::cout << "start" << std::endl;


			for (long long int i = 0; i < n/2; i++)
			{
				WEIGHTED_CFLOBDD_FOURIER_MUL SwapM = WeightedMatrix1234FourierMul::MkSwapGate(level+1, i, n-i-1);
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(SwapM, stateV);
			}

			std::cout << "loop start" << std::endl;
			for (long long int i = n-1; i >= 0; i--)
			{
				WEIGHTED_CFLOBDD_FOURIER_MUL H = Hadamard_i(n, i);
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(H, stateV);
				for (long int j = 0; j < i; j++)
				{
					fourierSemiring theta(1, boost::multiprecision::pow(BIG_INT(2), i - j + 1));
					WEIGHTED_CFLOBDD_FOURIER_MUL CP = WeightedMatrix1234FourierMul::MkCPGate(level+1, j, i, theta);
					stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(CP, stateV);
				}
			}

            /* IQFT
            for (long long int i = 0; i < n; i++)
			{
				WEIGHTED_CFLOBDD_FOURIER_MUL H = Hadamard_i(n, i);
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(H, stateV);
				for (long int j = i+1; j < n; j++)
				{
					fourierSemiring theta(-1, boost::multiprecision::pow(BIG_INT(2), -i + j + 1));
                    theta.SetComplexValue();
					WEIGHTED_CFLOBDD_FOURIER_MUL CP = WeightedMatrix1234FourierMul::MkCPGate(level+1, i, j, theta);
					stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(CP, stateV);
				}
			}
            */


		// 	std::cout << "done" << std::endl;


			return stateV;
		}

        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ShorsAlgoNew(int a, int N)
		{
			unsigned int level = ceil(log2(N));
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL allOnes = WeightedVectorComplexFloatBoostMul::NoDistinctionNode(level,1);
            allOnes = WeightedVectorComplexFloatBoostMul::VectorToMatrixInterleaved(allOnes);
			std::string s(2*N, '0');
			s[N-2] = '1';
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL e = WeightedVectorComplexFloatBoostMul::MkBasisVector(level + 1, s);
            // e.print(std::cout);
			// e = WeightedVectorComplexFloatBoostMul::VectorToMatrixInterleaved(e);
            // e.print(std::cout);
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL stateV = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(allOnes, e);

			for (int q = N-1; q >= 0; q--)
			{
				// std::cout << "q: " << q << std::endl;
				unsigned int power = std::pow(2, N-1-q);
				for (unsigned int i = 0; i < power; i++)
				{
					if (a == 2 || a == 13)
					{
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N, N+1);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
						CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N+1, N+2);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
						CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N+2, N+3);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
					}
					if (a == 7 || a == 8)
					{
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N+2, N+3);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
						CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N+1, N+2);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
						CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N, N+1);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
					}
					if (a == 4 || a == 11)
					{
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N+1, N+3);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
						CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N, N+2);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
					}
					if (a == 7 || a == 11 || a == 13)
					{
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL X = WeightedMatrix1234ComplexFloatBoostMul::MkExchangeInterleaved(level);
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL Id = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(level);
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL XI = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(X, Id);
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(level+1);
						WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL IXI = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I, XI);
						stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(IXI, stateV);
					}
				}
			}

			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(level+1);
			for (long long int i = 0; i < N; i++)
			{
				for (long int j = 0; j < i; j++)
				{
					double theta = -1*std::pow(2, j - i - 1);
					// std::cout << j << " " << i << std::endl;
					WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CP = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(level+1, j, i, theta);
					WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CPI = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(CP, I);
					stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CPI, stateV);
				}
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = Hadamard(2*N, i);
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateV);
			}
			for (long long int i = 0; i < N/2; i++)
			{
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL SwapM = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(level+1, i, N-i-1);
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL SwapMI = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(SwapM, I);
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(SwapMI, stateV);
			}

			return stateV;
		}


        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ShorsAlgo(int a, int N, int bits)
		{
			unsigned int level = ceil(log2(N));
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL allOnes = WeightedVectorComplexFloatBoostMul::NoDistinctionNode(level,1);
            allOnes = WeightedVectorComplexFloatBoostMul::VectorToMatrixInterleaved(allOnes);
			std::string s(2*N, '0');
			s[2*bits-2] = '1';
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL e = WeightedVectorComplexFloatBoostMul::MkBasisVector(level + 1, s);
            // e.print(std::cout);
			// e = WeightedVectorComplexFloatBoostMul::VectorToMatrixInterleaved(e);
            // e.print(std::cout);
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL stateV = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(allOnes, e);
            int alpha = std::floor(std::log2(a));
			for (int q = 2*bits-1; q >= 0; q--)
			{
				std::cout << "q: " << q << std::endl;
				unsigned int power = std::pow(2, 2*bits-1-q);
				for (unsigned int i = 0; i < power; i++)
				{
                    for (int k = alpha-1; k >= 0; k--)
                    {
                        for (int j = k; j < bits-1; j++)
                        {
                           WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CSWAP = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(level+2, q, N + j, N+j+1); 
                           stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CSWAP, stateV);
                        }
                    }
				}
			}

			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(level+1);
			for (long long int i = 0; i < 1/*2*bits*/; i++)
			{
				// for (long int j = 0; j < i; j++)
				// {
				// 	double theta = -1*std::pow(2, j - i - 1);
				// 	std::cout << j << " " << i << std::endl;
				// 	WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CP = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(level+1, j, i, theta);
				// 	WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CPI = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(CP, I);
				// 	stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(CPI, stateV);
				// }
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = Hadamard(2*N, i);
                H.print(std::cout);
                stateV.print(std::cout);
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateV);
			}
			for (long long int i = 0; i < bits; i++)
			{
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL SwapM = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(level+1, i, 2*bits-i-1);
				WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL SwapMI = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(SwapM, I);
				stateV = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(SwapMI, stateV);
			}

			return stateV;
		}

        void QFT_helper(int level, int b_vars, int total_vars, WEIGHTED_CFLOBDD_FOURIER_MUL& stateV)
        {
            // level is of the second half
            auto I = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level);
            int n = total_vars / 4;
            for (long long int i = 0; i < b_vars/2; i++)
			{
				WEIGHTED_CFLOBDD_FOURIER_MUL SwapM = WeightedMatrix1234FourierMul::MkSwapGate(level+1, i + n + total_vars/2, b_vars-i-1 + n + total_vars/2);
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(SwapM, stateV);
			}

			for (long long int i = b_vars-1; i >= 0; i--)
			{
				WEIGHTED_CFLOBDD_FOURIER_MUL H = Hadamard_i(std::pow(2, level), i + n + total_vars/2);
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(H, stateV);
				for (long int j = 0; j < i; j++)
				{
					fourierSemiring theta(1, boost::multiprecision::pow(BIG_INT(2), i - j + 1));
					WEIGHTED_CFLOBDD_FOURIER_MUL CP = WeightedMatrix1234FourierMul::MkCPGate(level+1, j + n + total_vars/2, i + n + total_vars/2, theta);
					stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(CP, stateV);
					// stateV.print(std::cout);
				}
			} 
        }

        void IQFT_helper(int level, int b_vars, int total_vars, WEIGHTED_CFLOBDD_FOURIER_MUL& stateV)
        {
			stateV = WeightedMatrix1234FourierMul::ComputeIQFT(level + 1, stateV, boost::multiprecision::pow(BIG_INT(2), b_vars), b_vars-1);
			return ;
            // level is of the second half
            auto I = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level);
            int n = total_vars / 4;

			for (long long int i = 0; i < b_vars; i++)
			{
				WEIGHTED_CFLOBDD_FOURIER_MUL H = Hadamard_i(std::pow(2, level), i + n + total_vars/2);
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(H, stateV);
				for (long int j = i+1; j < b_vars; j++)
				{
					fourierSemiring theta(-1, boost::multiprecision::pow(BIG_INT(2), -i + j + 1));
					WEIGHTED_CFLOBDD_FOURIER_MUL CP = WeightedMatrix1234FourierMul::MkCPGate(level+1, i + n + total_vars/2, j + n + total_vars/2, theta); 
					stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(CP, stateV);
				}
			} 

            for (long long int i = 0; i < b_vars/2; i++)
			{
				WEIGHTED_CFLOBDD_FOURIER_MUL SwapM = WeightedMatrix1234FourierMul::MkSwapGate(level+1, i + n + total_vars/2, b_vars-i-1 + n + total_vars/2);
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(SwapM, stateV);
			}
            stateV.root->rootConnection.factor = fourierSemiring(1, 1);
        }

		WEIGHTED_CFLOBDD_FOURIER_MUL NOT_i(unsigned int n, unsigned int i)
		{
			if (n == 1)
			{
				return WeightedMatrix1234FourierMul::MkExchangeInterleaved(1);
			}
			else {
                int level = ceil(log2(n/2));
				if (i < n/2)
				{
					WEIGHTED_CFLOBDD_FOURIER_MUL T = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level + 1);
					WEIGHTED_CFLOBDD_FOURIER_MUL H = NOT_i(n/2, i);
					return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(H, T);
				}
				else
				{
					WEIGHTED_CFLOBDD_FOURIER_MUL T = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level + 1);
					return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(T, NOT_i(n/2, i - n/2)); 
				}
			}
		}

        WEIGHTED_CFLOBDD_FOURIER_MUL MkADDGate(int level, BIG_INT a, int& count, int n)
        {
            if (count > n)
            {
                return WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level);
            }
            if (level == 1)
            {
               fourierSemiring a_modn(a, (int)std::pow(2, count));
               count += 1;
               return WeightedMatrix1234FourierMul::MkRZGate(level, a_modn); 
            }
            int vars = std::pow(2, level-2);
            auto lhs = MkADDGate(level-1, a, count, n);
            auto rhs = MkADDGate(level-1, a, count, n);
            return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(lhs, rhs);
        }

        void CADDModN(int controller, int x_c, int a, int N, WEIGHTED_CFLOBDD_FOURIER_MUL& stateV, int b_vars, int x_vars)
        {
            int level = stateV.root->level - 2; /* 4 */
            int n = (int)std::pow(2, level-1); // 8
            int count = 1;
            auto ADD_amodn = MkADDGate(level, a, count, b_vars);
            auto CADD_amodn = WeightedMatrix1234FourierMul::MkCADDGate(stateV.root->level, controller, x_c, ADD_amodn);
            stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(CADD_amodn, stateV);
            BIG_INT MinN = BIG_INT(std::pow(2, b_vars) - N);
            count = 1;
            auto ADD_MinNModn = WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level+1), 
                                    WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level), MkADDGate(level, MinN, count, b_vars)));
            stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(ADD_MinNModn, stateV);
            IQFT_helper(level+1, b_vars, 4*n, stateV);
			auto Cij = WeightedMatrix1234FourierMul::MkCNOT(stateV.root->level, 4 * n, x_vars + 2 * n, 3 * n);
			// Cij.print(std::cout);
			auto SW = WeightedMatrix1234FourierMul::MkSwapGate(level + 2, x_vars + 2*n, 3*n);
			// SW.print(std::cout);
			Cij = WeightedMatrix1234FourierMul::MatrixMultiplyV4(Cij, SW);
			Cij = WeightedMatrix1234FourierMul::MatrixMultiplyV4(SW, Cij);
			// Cij.print(std::cout);
			stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(Cij, stateV);
            QFT_helper(level + 1, b_vars, 4 * n, stateV);
			count = 1;
			auto ADD_NModn = WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level + 1),
							WeightedMatrix1234FourierMul::MkCADDGate2(stateV.root->level-1, x_vars, MkADDGate(level, BIG_INT(N), count, b_vars)));
			stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(ADD_NModn, stateV);
			BIG_INT MinA = BIG_INT(std::pow(2, b_vars) - a);
			count = 1;
			auto CADD_Minamodn = WeightedMatrix1234FourierMul::MkCADDGate(stateV.root->level, controller, x_c, MkADDGate(level, MinA, count, b_vars));
			stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(CADD_Minamodn, stateV);
            IQFT_helper(level+1, b_vars, 4*n, stateV);
			WEIGHTED_CFLOBDD_FOURIER_MUL X = NOT_i(std::pow(2, stateV.root->level-1), 3*n);
			stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(X, stateV);
			stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(Cij, stateV);
			stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(X, stateV);
			QFT_helper(level + 1, b_vars, 4 * n, stateV);
			stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(CADD_amodn, stateV);
        }

        void CMULT(int controller, int powerOfa, int N, WEIGHTED_CFLOBDD_FOURIER_MUL& stateV, int N_vars /* 5 */, int b_vars /* 6 */, int total_vars /* 32 */, int index)
        {
            int n = total_vars / 4;
            QFT_helper(stateV.root->level-1, b_vars, total_vars, stateV);
			int a = 1;
            for (int i = 0; i < N_vars; i++)
            {
				if (i == 0)
					a = powerOfa;
				else
					a = (2 * a) % N;
                CADDModN(controller /* c */, (N_vars - 1 - i + total_vars / 2) /* x */, a , N, stateV, b_vars, N_vars);
				// if (index == 12)
				// {
				// 	std::cout << "CMULT i " << i << std::endl;
				// 	stateV.print(std::cout);
				// }
            }
            IQFT_helper(stateV.root->level-1, b_vars, total_vars, stateV);
        }

        // Source: geeks for geeks
        int gcdExtended(int a, int b, int* x, int* y)
        {
            // Base Case
            if (a == 0) {
                *x = 0, *y = 1;
                return b;
            }
        
            // To store results of recursive call
            int x1, y1;
            int gcd = gcdExtended(b % a, a, &x1, &y1);
        
            // Update x and y using results of recursive
            // call
            *x = y1 - (b / a) * x1;
            *y = x1;
        
            return gcd;
        }

        int modInverse(int a, int N)
        {
            int x, y;
            int g = gcdExtended(a, N, &x, &y);
            assert(g == 1);
            int res = (x % N + N) % N;
            return res;
        }
        
        WEIGHTED_CFLOBDD_FOURIER_MUL HadamardMatrix(unsigned int level, int end /*15*/, int n/*9*/)
        {
            if (end == n){
                return WeightedMatrix1234FourierMul::MkWalshInterleaved(level);
            }
            if (n < 0)
                return WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level);
            int mid = (end)/2;
            auto lhs = HadamardMatrix(level-1, mid, std::min(mid, n));
            int new_end = end - mid - 1;
            int new_n = n - mid - 1;
            auto rhs = HadamardMatrix(level-1, new_end, new_n);
            return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(lhs, rhs);
        } 

		WEIGHTED_CFLOBDD_FOURIER_MUL R_i(unsigned int n, unsigned int i, fourierSemiring R)
		{
			if (n == 1)
			{
				return WeightedMatrix1234FourierMul::MkRZGate(1, R);
			}
			else {
                int level = ceil(log2(n/2));
				if (i < n/2)
				{
					WEIGHTED_CFLOBDD_FOURIER_MUL T = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level + 1);
					WEIGHTED_CFLOBDD_FOURIER_MUL H = R_i(n/2, i, R);
					return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(H, T);
				}
				else
				{
					WEIGHTED_CFLOBDD_FOURIER_MUL T = WeightedMatrix1234FourierMul::MkIdRelationInterleaved(level + 1);
					return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(T, R_i(n/2, i - n/2, R)); 
				}
			}
		}

		WEIGHTED_CFLOBDD_FOURIER_MUL ComputeKroneckerProduct(std::vector<WEIGHTED_CFLOBDD_FOURIER_MUL>& m, int start, int end)
		{
			if (end - start == 0)
				return m[start];
			int mid = (end - start)/2 + start;
			auto lhs = ComputeKroneckerProduct(m, mid + 1, end);
			auto rhs = ComputeKroneckerProduct(m, start, mid);
			return WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(rhs, lhs);
		}
        
        
        std::tuple<WEIGHTED_CFLOBDD_FOURIER_MUL, std::string, BIG_INT> ShorsFourier(int a, int N)
        {
            // NEED 2logN+3 BITS. But we will be using 4 * nearestPow2(n)
            // Consider N = 21
            // <c, x, b> == <2*n, n, n + 2>
            int x_vars = ceil(log2(N)); // 5
            int b_vars = x_vars + 1; // 6
            int total_vars = pow(2, ceil(log2(b_vars))); // 8 bits
            total_vars = 2 * total_vars; // 16
            total_vars = 2 * total_vars; // 32

            int level = log2(total_vars); // of the vector (level-1 x level-1) x (level-1 x level-1)

            WEIGHTED_CFLOBDD_FOURIER_MUL e0 = WeightedVectorFourierMul::MkBasisVector(level, 0); // matrix
            // auto H = HadamardMatrix(level, total_vars/2 - 1, 2 * x_vars-1);
			auto H = Hadamard_i(std::pow(2, level-1), 0);
            e0 = WeightedMatrix1234FourierMul::MatrixMultiplyV4(H, e0);
            std::string s(total_vars, '0');
            s[2*x_vars-2] = '1'; // (00001) x (00000)

            WEIGHTED_CFLOBDD_FOURIER_MUL e1 = WeightedVectorFourierMul::MkBasisVector(level, s);
            auto stateV = WeightedMatrix1234FourierMul::KroneckerProduct2Vocs(e0, e1);

			std::vector<long long int> powersOfa;
			for (unsigned int i = 0; i < 2 * x_vars; i++)
			{
				if (i == 0)
					powersOfa.push_back(a % N);
				else
				{
					long long int prev = powersOfa[powersOfa.size()-1];
					powersOfa.push_back((prev * prev) % N);
				}
			}
            // Let's create U_a^{2^i}
            long long int prevPowOfa = 1;
			std::string sampled_string = "";
			BIG_INT sampled_number = 0;
			std::vector<WEIGHTED_CFLOBDD_FOURIER_MUL> m0s, m1s;
			auto no_dist = WeightedVectorFourierMul::MkBasisVector(1, 0) + WeightedVectorFourierMul::MkBasisVector(1, 2);
			for (long int i = 0; i < total_vars; i++)
			{
				if (i == 0)
				{
					m0s.push_back(WeightedVectorFourierMul::MkBasisVector(1, 0));
					m1s.push_back(WeightedVectorFourierMul::MkBasisVector(1, 2));
				}
				else
				{
					m0s.push_back(no_dist);
					m1s.push_back(no_dist);	
				}
			}
			WEIGHTED_CFLOBDD_FOURIER_MUL Measure0 = ComputeKroneckerProduct(m0s, 0, m0s.size() - 1);
			WEIGHTED_CFLOBDD_FOURIER_MUL Measure1 = ComputeKroneckerProduct(m1s, 0, m1s.size() - 1);

            for (int i = 0; i < 2 * x_vars; i++)
            {
                // U_a => CMULT(a)modN SWAP CMULT(a^-1)modN
                // CMULT(a)modN = QFT Add_{2^k a} mod N, k \in 0..N_vars-1 QFT^-1
                std::cout << i << std::endl;
                // int controller = (2 * x_vars - 1 - i);
				// int controller = i;
				int controller = 0;
                long long int powerOfa = powersOfa[2 * x_vars - 1 - i];
				// std::cout << powerOfa << std::endl;
                // if (i == 0)
                //     powerOfa = a % N;
                // else
                //     powerOfa = (prevPowOfa * prevPowOfa) % N;
                // prevPowOfa = powerOfa;
				// if (i == 12)
				// 	stateV.print(std::cout);
                CMULT(controller, powerOfa, N, stateV, x_vars, b_vars, total_vars, i);
				// if (i == 12)
				// 	stateV.print(std::cout);
                // SWAP
                for (int x = 0; x < x_vars; x++)
                {
                    auto SW = WeightedMatrix1234FourierMul::MkCSwapGate(level + 1, controller, total_vars/2 + x, total_vars/2 + total_vars/4 + x + 1);
                    stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(SW, stateV);
                }
                int powOfaInv = modInverse(powerOfa, N);
				int minusInv = N - powOfaInv;
				// stateV.print(std::cout);
				// stateV = WeightedMatrix1234FourierMul::MkSetBToZero(stateV.root->level, stateV);
                CMULT(controller, minusInv, N, stateV, x_vars, b_vars, total_vars, i);
				
				auto stateV_tmp = stateV;
				fourierSemiring R;
				if (i == 0)
					R = fourierSemiring(1, 1);
				else
				{
					R = fourierSemiring( boost::multiprecision::pow(BIG_INT(2), (i + 1)) - sampled_number, boost::multiprecision::pow(BIG_INT(2), (i + 1)));
				}
				// std::cout << R << std::endl;
				auto R_M = R_i(std::pow(2, level), 0, R);
				stateV_tmp = WeightedMatrix1234FourierMul::MatrixMultiplyV4(R_M, stateV_tmp);
				// stateV_tmp.print(std::cout);
				auto H_M = Hadamard_i(std::pow(2, level), 0);
				// H_M.print(std::cout);
				stateV_tmp = WeightedMatrix1234FourierMul::MatrixMultiplyV4(H_M, stateV_tmp);
				// stateV_tmp.print(std::cout);
				// abort();
				auto res = WeightedMatrix1234FourierMul::MeasureAndReset(stateV.root->level, total_vars, stateV_tmp, R);
				// stateV_tmp.print(std::cout);
				if (res.second == 0){
					stateV = Measure0 * stateV;
				}
				else
					stateV = Measure1 * stateV;
				stateV = WeightedMatrix1234FourierMul::MatrixMultiplyV4(H_M, stateV);
				stateV = WeightedMatrix1234FourierMul::ResetState(stateV.root->level, stateV);
				sampled_string = (res.second == 0 ? "0" : "1") + sampled_string;
				sampled_number = boost::multiprecision::pow(BIG_INT(2), i + 1) * res.second + sampled_number;
				// std::cout << "sampled_string: " << sampled_string << " sampled_number: " << sampled_number << std::endl;
				// stateV.print(std::cout);
            }

            std::cout << "partly done" << std::endl;
			std::cout << "sampled_string: " << sampled_string << std::endl;
			return std::make_tuple(stateV, sampled_string, sampled_number);
        }



		
    }
}