#ifndef W_MATRIX1234_FB_MUL_NODE_GUARD
#define W_MATRIX1234_FB_MUL_NODE_GUARD

#include <map>
#include <unordered_map>
#include "weighted_cflobdd_node_t.h"
#include "return_map_T.h"
#include "weighted_matmult_map.h"

namespace CFL_OBDD {

    namespace WeightedMatrix1234FloatBoostMul {

        typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
        // typedef double BIG_FLOAT;
        //typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000> > BIG_FLOAT;
        typedef WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostMulNodeHandle;
        typedef WeightedCFLOBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostInternalNode;
        typedef WeightedCFLOBDDDontCareNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostDontCareNode;
        typedef WeightedCFLOBDDForkNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostForkNode;
        typedef WeightedCFLOBDDLeafNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDFloatBoostLeafNode;
        typedef WConnection<BIG_FLOAT, std::multiplies<BIG_FLOAT>> Connection;
	    typedef ReturnMapHandle<WeightedMatMultMapHandle<BIG_FLOAT>> CFLOBDDMatMultMapHandle;
        typedef std::tuple<WeightedCFLOBDDFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_FLOAT> MatMultReturnT;

        class WeightedMatMultPair{
        public:
            WeightedCFLOBDDFloatBoostMulNodeHandle m1;
            WeightedCFLOBDDFloatBoostMulNodeHandle m2;
            WeightedMatMultPair(WeightedCFLOBDDFloatBoostMulNodeHandle p1, WeightedCFLOBDDFloatBoostMulNodeHandle p2)
            {
                m1 = p1;
                m2 = p2;
            }

            struct MatMultPairHash {
                size_t operator()(const WeightedMatMultPair& p) const
                {
                    WeightedCFLOBDDFloatBoostMulNodeHandle t1 = p.m1;
                    WeightedCFLOBDDFloatBoostMulNodeHandle t2 = p.m2;
                    auto hash1 = t1.Hash(997);
                    auto hash2 = t2.Hash(997);
                    return 117 * (hash1 + 1) + hash2;
                }
            };

            bool operator==(const WeightedMatMultPair& p) const
            {
                WeightedCFLOBDDFloatBoostMulNodeHandle m11 = m1;
                WeightedCFLOBDDFloatBoostMulNodeHandle m12 = m2;
                WeightedCFLOBDDFloatBoostMulNodeHandle m21 = p.m1;
                WeightedCFLOBDDFloatBoostMulNodeHandle m22 = p.m2;
                return (m11 == m21) && (m12 == m22);
            }
        };

        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkIdRelationInterleavedNode(unsigned int level);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkWalshInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkInverseReedMullerInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkNegationMatrixInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkCNOTInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkExchangeInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkCNOTNode(unsigned int level, unsigned int n, long int controller, long int controlled);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkCPGateNode(unsigned int level, long int controller, long int controlled);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num);
        extern WeightedCFLOBDDFloatBoostMulNodeHandle MkCSwapGateNode(unsigned int level, long int controller, long int i, long int j, int case_num);

        extern WeightedCFLOBDDFloatBoostMulNodeHandle KroneckerProduct2VocsNode(WeightedCFLOBDDFloatBoostMulNodeHandle m1, WeightedCFLOBDDFloatBoostMulNodeHandle m2, int zero_index_m1, int zero_index_m2); 
        extern std::tuple<WeightedCFLOBDDFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_FLOAT>
		MatrixMultiplyV4Node(WeightedCFLOBDDFloatBoostMulNodeHandle c1, WeightedCFLOBDDFloatBoostMulNodeHandle c2, int zero_exit_1, int zero_exit_2);
        
        // Initialization routine that needs to be called before any call to MatrixProjectVoc23Node
        extern void Matrix1234InitializerNode();  // Empty for now

        // extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4Node(
        // 	std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash>& hashMap,
        // 	WeightedCFLOBDDFloatBoostMulNodeHandle c1, WeightedCFLOBDDFloatBoostMulNodeHandle c2);
        // extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4WithInfoNode(
        // 	std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash>& hashMap,
        // 	WeightedCFLOBDDFloatBoostMulNodeHandle c1, WeightedCFLOBDDFloatBoostMulNodeHandle c2, int c1_zero_index, int c2_zero_index);
    }
 }

#endif

