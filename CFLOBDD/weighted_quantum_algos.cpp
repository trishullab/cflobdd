#include "weighted_quantum_algos.h"
#include "wvector_fb_mul.h"
#include <ctime>
#include <time.h>
#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>
#include "wmatrix1234_complex_fb_mul.h"
#include "wvector_complex_fb_mul.h"
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
			std::string last_one(n, '0');
			last_one[0] = '1';
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL e1 = WeightedVectorFloatBoostMul::MkBasisVector(level - 2, last_one);
            e1 = WeightedVectorFloatBoostMul::VectorToMatrixInterleaved(e1);
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
			std::string last_one(n, '0');
			last_one[0] = '1';
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL Y = WeightedVectorFloatBoostMul::MkBasisVector(level, last_one);
			Y = WeightedVectorFloatBoostMul::VectorToMatrixInterleaved(Y);
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

        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL QFT(long long int n, std::string s)
		{
			unsigned int level = ceil(log2(n));
			
			// std::cout << "s: " << s << std::endl;
			// std::reverse(s.begin(), s.end());
			WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL stateV = WeightedVectorComplexFloatBoostMul::MkBasisVector(level, s);
			stateV = WeightedVectorComplexFloatBoostMul::VectorToMatrixInterleaved(stateV);
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
                BIG_COMPLEX_FLOAT c = 1.0/sqrt(2);
                H = c * H;
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
			std::cout << "done" << std::endl;


			return stateV;
		}



    }
}