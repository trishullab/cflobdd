#include "weighted_quantum_algos.h"
#include "wvector_fb_mul.h"
#include <ctime>
#include <time.h>
#include <chrono>
using namespace std;
using namespace chrono;

namespace CFL_OBDD
{
    namespace WeightedQuantumAlgos
    {

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
    }
}