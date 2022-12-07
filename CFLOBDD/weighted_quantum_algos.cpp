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
                auto x = high_resolution_clock::now(); 
				WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL tmp = WeightedMatrix1234FloatBoostMul::MkCNOT(level, 2*n, i, n);
                auto y = high_resolution_clock::now(); 
				F = WeightedMatrix1234FloatBoostMul::MatrixMultiplyV4(F, tmp);
                auto z = high_resolution_clock::now(); 
                auto dx = duration_cast<milliseconds>(y - x); 
                auto dy = duration_cast<milliseconds>(z - y); 
                std::cout << dx.count() << " " << dy.count() << std::endl;
			}
            auto s2 = high_resolution_clock::now();
			WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL e0 = WeightedVectorFloatBoostMul::NoDistinctionNode(level - 1, 1);
			std::string last_one(n, '0');
			last_one[n - 1] = '1';
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
            auto s5 = high_resolution_clock::now();
            auto d1 = duration_cast<milliseconds>(s2 - s1);
            auto d2 = duration_cast<milliseconds>(s3 - s2);
            auto d3 = duration_cast<milliseconds>(s4 - s3);
            auto d4 = duration_cast<milliseconds>(s5 - s4);
            std::cout << d1.count() << " " << d2.count() << " " << d3.count() << " " << d4.count() << std::endl;
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			// stateV = VectorFloatBoost::VectorWithAmplitude(stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 4 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			// stateV.CountPaths();
			std::string ans_s = "";
            // ans_s = VectorFloatBoost::Sampling(stateV, true).substr(0, n + 1);
			return std::make_pair(ans_s, stateV);
        }
    }
}