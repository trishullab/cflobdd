#ifndef W_COMPLEX_VECTOR_NODE_GUARD
#define W_COMPLEX_VECTOR_NODE_GUARD

#include <string.h>
#include <unordered_map>
#include <boost/multiprecision/cpp_complex.hpp>
#include "weighted_cflobdd_node_t.h"
#include "wvector_complex_fb_mul_bdd_node.h"

namespace CFL_OBDD {
	

    namespace WeightedVectorComplexFloatBoostMul {

	    typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;

        typedef WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostMulNodeHandle;
        typedef WeightedCFLOBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostInternalNode;
        typedef WeightedBDDTopNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostTopNode;
        typedef WeightedCFLOBDDDontCareNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostDontCareNode;
        typedef WeightedCFLOBDDForkNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostForkNode;
        typedef WeightedCFLOBDDLeafNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostLeafNode;
        typedef WConnection<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> Connection;
        typedef WeightedCFLOBDDNodeMemoTable<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostMulNodeMemoTable;
        typedef WeightedCFLOBDDNodeMemoTableRefPtr<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostMulNodeMemoTableRefPtr;


        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, unsigned int index, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, std::string s, int cflobdd_kind = 1);
        extern long double getNonZeroProbabilityNode(WeightedCFLOBDDComplexFloatBoostMulNodeHandle nh);

        extern void VectorInitializerNode();  // Empty for now

        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle VectorToMatrixInterleavedNode(std::unordered_map<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, WeightedCFLOBDDComplexFloatBoostMulNodeHandle, WeightedCFLOBDDComplexFloatBoostMulNodeHandle::WeightedCFLOBDDNodeHandleT_Hash> hashMap,
                                                        WeightedCFLOBDDComplexFloatBoostMulNodeHandle& nh);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkColumn1MatrixNode(unsigned int level);
        extern std::pair<std::string,std::string> SamplingNode(WeightedCFLOBDDComplexFloatBoostMulNodeHandle nh, unsigned int index, std::mt19937 mt, std::uniform_real_distribution<double> dis, bool voctwo = false, std::string = "");
        // needs to be removed and linked to the one in cflobdd_node.cpp
        extern long double addNumPathsToExit(long double path1, long double path2);
        extern long double addNumPathsToExit(std::vector<long double>& paths);
        extern BIG_COMPLEX_FLOAT addNumPathsToExit(std::vector<BIG_COMPLEX_FLOAT>& paths);
        long double getLogSumNumPaths(std::vector<std::pair<long double, unsigned int>>& numBPaths, unsigned int size);
        BIG_COMPLEX_FLOAT getLogSumNumPaths(std::vector<std::pair<BIG_COMPLEX_FLOAT, unsigned int>>& numBPaths, unsigned int size);
        
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

