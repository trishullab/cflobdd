#ifndef W_COMPLEX_VECTOR_BDD_NODE_GUARD
#define W_COMPLEX_VECTOR_BDD_NODE_GUARD

#include <string.h>
#include <unordered_map>
#include <random>
#include <boost/multiprecision/cpp_complex.hpp>
#include "weighted_bdd_node_t.h"

namespace CFL_OBDD {
	

    namespace WeightedVectorBDDComplexFloatBoostMul {

	    typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;

        typedef WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostMulNodeHandle;
        typedef WeightedBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostInternalNode;
        typedef WeightedBDDLeafNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostLeafNode;


        extern WeightedBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, unsigned int index, long int id = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, std::string s, long int id = 0);
        
        extern std::pair<std::string,std::string> SamplingNode(WeightedBDDComplexFloatBoostMulNodeHandle nh, unsigned int numVars, std::mt19937 mt, std::uniform_real_distribution<double> dis, long int count = 0);
    }
}

#endif

