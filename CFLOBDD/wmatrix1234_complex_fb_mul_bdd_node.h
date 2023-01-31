#ifndef W_MATRIX1234_COMPLEX_FB_MUL_BDD_NODE_GUARD
#define W_MATRIX1234_COMPLEX_FB_MUL_BDD_NODE_GUARD

#include <map>
#include <unordered_map>
#include <boost/multiprecision/cpp_complex.hpp>
#include "weighted_bdd_node_t.h"
#include "return_map_T.h"
#include "weighted_matmult_map.h"

namespace CFL_OBDD {

    namespace WeightedMatrix1234BDDComplexFloatBoostMul {

	    typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
        typedef WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostMulNodeHandle;
        typedef WeightedBDDNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostMulNode;
        typedef WeightedBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostMulInternalNode;
        typedef WeightedBDDLeafNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostMulLeafNode;
        typedef std::tuple<WeightedBDDComplexFloatBoostMulNodeHandle, std::vector<BIG_COMPLEX_FLOAT>, BIG_COMPLEX_FLOAT> BDDMatMultReturnT;

        extern WeightedBDDComplexFloatBoostMulNodeHandle MkIdRelationInterleavedNode(unsigned int numVars, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkWalshInterleavedNode(unsigned int i, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkInverseReedMullerInterleavedNode(unsigned int i);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkNegationMatrixInterleavedNode(unsigned int i, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkPauliYGateNode(unsigned int i, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkPauliZGateNode(unsigned int i, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkSGateNode(unsigned int i, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkCNOTInterleavedNode(unsigned int i);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkExchangeInterleavedNode(unsigned int i);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkCNOTNode(unsigned int numVars, unsigned int n, long int controller, long int controlled, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkCCNOTNode(unsigned int numVars, long int controller1, long int controller2, long int controlled, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkCPGateNode(unsigned int level, long int controller, long int controlled, BIG_COMPLEX_FLOAT theta_val, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkPhaseShiftGateNode(unsigned int level, BIG_COMPLEX_FLOAT theta_val, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkiSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkCSwapGateNode(unsigned int level, long int controller, long int i, long int j, int case_num, unsigned int index = 0);
        extern WeightedBDDComplexFloatBoostMulNodeHandle MkRestrictNode(unsigned int level, std::string s, unsigned int index = 0);

        class RenamePair{
        public:
            WeightedBDDComplexFloatBoostMulNodeHandle m1;
            long int offset;
            RenamePair(WeightedBDDComplexFloatBoostMulNodeHandle p1, long int p2)
            {
                m1 = p1;
                offset = p2;
            }

            struct RenamePairHash {
                size_t operator()(const RenamePair& p) const
                {
                    WeightedBDDComplexFloatBoostMulNodeHandle t1 = p.m1;
                    long int t2 = p.offset;
                    auto hash1 = t1.Hash(997);
                    auto hash2 = t2;
                    return 117 * (hash1 + 1) + hash2;
                }
            };

            bool operator==(const RenamePair& p) const
            {
                WeightedBDDComplexFloatBoostMulNodeHandle m11 = m1;
                long int m12 = offset;
                WeightedBDDComplexFloatBoostMulNodeHandle m21 = p.m1;
                long int m22 = p.offset;
                return (m11 == m21) && (m12 == m22);
            }
        };

        class WeightedBDDMatMultPair{
        public:
            WeightedBDDComplexFloatBoostMulNodeHandle m1;
            WeightedBDDComplexFloatBoostMulNodeHandle m2;
            WeightedBDDMatMultPair(WeightedBDDComplexFloatBoostMulNodeHandle p1, WeightedBDDComplexFloatBoostMulNodeHandle p2)
            {
                m1 = p1;
                m2 = p2;
            }

            struct MatMultPairHash {
                size_t operator()(const WeightedBDDMatMultPair& p) const
                {
                    WeightedBDDComplexFloatBoostMulNodeHandle t1 = p.m1;
                    WeightedBDDComplexFloatBoostMulNodeHandle t2 = p.m2;
                    auto hash1 = t1.Hash(997);
                    auto hash2 = t2.Hash(997);
                    return 117 * (hash1 + 1) + hash2;
                }
            };

            bool operator==(const WeightedBDDMatMultPair& p) const
            {
                WeightedBDDComplexFloatBoostMulNodeHandle m11 = m1;
                WeightedBDDComplexFloatBoostMulNodeHandle m12 = m2;
                WeightedBDDComplexFloatBoostMulNodeHandle m21 = p.m1;
                WeightedBDDComplexFloatBoostMulNodeHandle m22 = p.m2;
                return (m11 == m21) && (m12 == m22);
            }
        };

        extern WeightedBDDComplexFloatBoostMulNodeHandle KroneckerProduct2VocsNode(std::unordered_map<WeightedBDDComplexFloatBoostMulNodeHandle, WeightedBDDComplexFloatBoostMulNodeHandle, WeightedBDDComplexFloatBoostMulNodeHandle::WeightedBDDNodeHandle_Hash>& hashMap, WeightedBDDComplexFloatBoostMulNodeHandle m1, WeightedBDDComplexFloatBoostMulNodeHandle m2, long int numVars); 
        
        extern std::tuple<WeightedBDDComplexFloatBoostMulNodeHandle, std::vector<BIG_COMPLEX_FLOAT>, BIG_COMPLEX_FLOAT>
         MatrixMultiplyV4Node(WeightedBDDComplexFloatBoostMulNodeHandle m1, WeightedBDDComplexFloatBoostMulNodeHandle m2, long int numVars, long int count = 0);
    }
 }

#endif

