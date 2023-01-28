#ifndef W_FOURIER_VECTOR_NODE_GUARD
#define W_FOURIER_VECTOR_NODE_GUARD

#include <string.h>
#include <unordered_map>
#include "fourier_semiring.h"
#include "weighted_cflobdd_node_t.h"

namespace CFL_OBDD {
	

    namespace WeightedVectorFourierMul {


        typedef WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierMulNodeHandle;
        typedef WeightedCFLOBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierInternalNode;
        typedef WeightedCFLOBDDDontCareNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierDontCareNode;
        typedef WeightedCFLOBDDForkNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierForkNode;
        typedef WeightedCFLOBDDLeafNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierLeafNode;
        typedef WConnection<fourierSemiring, std::multiplies<fourierSemiring>> Connection;
        typedef WeightedCFLOBDDNodeMemoTable<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierMulNodeMemoTable;
        typedef WeightedCFLOBDDNodeMemoTableRefPtr<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierMulNodeMemoTableRefPtr;


        extern WeightedCFLOBDDFourierMulNodeHandle MkBasisVectorNode(unsigned int level, unsigned int index);
        extern WeightedCFLOBDDFourierMulNodeHandle MkBasisVectorNode(unsigned int level, std::string s);
        
        extern void VectorInitializerNode();  // Empty for now

        extern WeightedCFLOBDDFourierMulNodeHandle VectorToMatrixInterleavedNode(std::unordered_map<WeightedCFLOBDDFourierMulNodeHandle, WeightedCFLOBDDFourierMulNodeHandle, WeightedCFLOBDDFourierMulNodeHandle::WeightedCFLOBDDNodeHandleT_Hash> hashMap,
                                                        WeightedCFLOBDDFourierMulNodeHandle& nh);
        extern WeightedCFLOBDDFourierMulNodeHandle MkColumn1MatrixNode(unsigned int level);
        extern std::pair<std::string,std::string> SamplingNode(WeightedCFLOBDDFourierMulNodeHandle nh, unsigned int index, bool voctwo = false, std::string = "");
        // needs to be removed and linked to the one in cflobdd_node.cpp
        extern long double addNumPathsToExit(long double path1, long double path2);
        extern long double addNumPathsToExit(std::vector<long double>& paths);
        extern fourierSemiring addNumPathsToExit(std::vector<fourierSemiring>& paths);
        long double getLogSumNumPaths(std::vector<std::pair<long double, unsigned int>>& numBPaths, unsigned int size);
        fourierSemiring getLogSumNumPaths(std::vector<std::pair<fourierSemiring, unsigned int>>& numBPaths, unsigned int size);
        
        template <typename T>
        bool sortNumPathPairs(const std::pair<T, unsigned int>& p1, const std::pair<T, unsigned int> &p2){
            if (p1.first < p2.first)
                return true;
            else if (p1.first > p2.first)
                return false;
            return p1.second < p2.second;
        }
    }
}

#endif

