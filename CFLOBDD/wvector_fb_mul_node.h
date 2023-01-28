#ifndef W_VECTOR_NODE_GUARD
#define W_VECTOR_NODE_GUARD

#include <string.h>
#include <unordered_map>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "weighted_cflobdd_node_t.h"

namespace CFL_OBDD {
	

    namespace WeightedVectorFloatBoostMul {

	    typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;

        typedef WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostMulNodeHandle;
        typedef WeightedCFLOBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostInternalNode;
        typedef WeightedCFLOBDDDontCareNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostDontCareNode;
        typedef WeightedCFLOBDDForkNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostForkNode;
        typedef WeightedCFLOBDDLeafNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostLeafNode;
        typedef WConnection<BIG_FLOAT, std::multiplies<BIG_FLOAT>> Connection;
        typedef WeightedCFLOBDDNodeMemoTable<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostMulNodeMemoTable;
        typedef WeightedCFLOBDDNodeMemoTableRefPtr<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostMulNodeMemoTableRefPtr;


        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, unsigned int index);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, std::string s);
        
        extern void VectorInitializerNode();  // Empty for now

        extern WeightedCFLOBDDFloatBoostMulNodeHandle VectorToMatrixInterleavedNode(std::unordered_map<WeightedCFLOBDDFloatBoostMulNodeHandle, WeightedCFLOBDDFloatBoostMulNodeHandle, WeightedCFLOBDDFloatBoostMulNodeHandle::WeightedCFLOBDDNodeHandleT_Hash> hashMap,
                                                        WeightedCFLOBDDFloatBoostMulNodeHandle& nh);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkColumn1MatrixNode(unsigned int level);
        extern std::pair<std::string,std::string> SamplingNode(WeightedCFLOBDDFloatBoostMulNodeHandle nh, unsigned int index, bool voctwo = false, std::string = "");
        // needs to be removed and linked to the one in cflobdd_node.cpp
        extern long double addNumPathsToExit(long double path1, long double path2);
        extern long double addNumPathsToExit(std::vector<long double>& paths);
        extern BIG_FLOAT addNumPathsToExit(std::vector<BIG_FLOAT>& paths);
        long double getLogSumNumPaths(std::vector<std::pair<long double, unsigned int>>& numBPaths, unsigned int size);
        BIG_FLOAT getLogSumNumPaths(std::vector<std::pair<BIG_FLOAT, unsigned int>>& numBPaths, unsigned int size);
        
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

