#ifndef W_MATRIX1234_COMPLEX_FB_MUL_NODE_GUARD
#define W_MATRIX1234_COMPLEX_FB_MUL_NODE_GUARD

#include <map>
#include <unordered_map>
#include <boost/multiprecision/cpp_complex.hpp>
#include "weighted_cflobdd_node_t.h"
#include "return_map_T.h"
#include "weighted_matmult_map.h"

namespace CFL_OBDD {

    namespace WeightedMatrix1234ComplexFloatBoostMul {

	    typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
        // typedef double BIG_COMPLEX_FLOAT;
        //typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000> > BIG_COMPLEX_FLOAT;
        typedef WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostMulNodeHandle;
        typedef WeightedCFLOBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostInternalNode;
        typedef WeightedCFLOBDDDontCareNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostDontCareNode;
        typedef WeightedCFLOBDDForkNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostForkNode;
        typedef WeightedCFLOBDDLeafNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDComplexFloatBoostLeafNode;
        typedef WeightedBDDTopNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedBDDComplexFloatBoostTopNode;
        typedef WConnection<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> Connection;
	    typedef ReturnMapHandle<WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>> CFLOBDDMatMultMapHandle;
        typedef std::tuple<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_COMPLEX_FLOAT> MatMultReturnT;

        class WeightedMatMultPair{
        public:
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle m1;
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle m2;
            WeightedMatMultPair(WeightedCFLOBDDComplexFloatBoostMulNodeHandle p1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle p2)
            {
                m1 = p1;
                m2 = p2;
            }

            struct MatMultPairHash {
                size_t operator()(const WeightedMatMultPair& p) const
                {
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle t1 = p.m1;
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle t2 = p.m2;
                    auto hash1 = t1.Hash(997);
                    auto hash2 = t2.Hash(997);
                    return 117 * (hash1 + 1) + hash2;
                }
            };

            bool operator==(const WeightedMatMultPair& p) const
            {
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle m11 = m1;
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle m12 = m2;
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle m21 = p.m1;
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle m22 = p.m2;
                return (m11 == m21) && (m12 == m22);
            }
        };

        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkIdRelationInterleavedNode(unsigned int level, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkWalshInterleavedNode(unsigned int i, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkInverseReedMullerInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkNegationMatrixInterleavedNode(unsigned int i, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkPauliYGateNode(unsigned int i, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkPauliZGateNode(unsigned int i, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkSGateNode(unsigned int i, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCNOTInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkExchangeInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCNOTNode(unsigned int level, unsigned int n, long int controller, long int controlled, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCCNOTNode(unsigned int level, long int controller1, long int contoller2, long int controlled, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCPGateNode(std::unordered_map<std::string, WeightedCFLOBDDComplexFloatBoostMulNodeHandle>& cp_hashMap, unsigned int level, long int controller, long int controlled, BIG_COMPLEX_FLOAT theta_val, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkPhaseShiftGateNode(unsigned int level, BIG_COMPLEX_FLOAT theta_val, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkiSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num, int cflobdd_kind = 1);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCSwapGateNode(unsigned int level, long int controller, long int i, long int j, int case_num);
        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCSwapGate2Node(unsigned int level, long int controller, long int i, long int j, int case_num, int cflobdd_kind = 1);

        extern std::pair<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, int> MkRestrictNode(unsigned int level, std::string s, int cflobdd_kind = 1);

        extern WeightedCFLOBDDComplexFloatBoostMulNodeHandle KroneckerProduct2VocsNode(WeightedCFLOBDDComplexFloatBoostMulNodeHandle m1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle m2, int zero_index_m1, int zero_index_m2); 
        extern std::tuple<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_COMPLEX_FLOAT>
		MatrixMultiplyV4Node(WeightedCFLOBDDComplexFloatBoostMulNodeHandle c1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle c2, int zero_exit_1, int zero_exit_2);
        
        // Initialization routine that needs to be called before any call to MatrixProjectVoc23Node
        extern void Matrix1234InitializerNode();  // Empty for now

        // extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4Node(
        // 	std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash>& hashMap,
        // 	WeightedCFLOBDDComplexFloatBoostMulNodeHandle c1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle c2);
        // extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4WithInfoNode(
        // 	std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash>& hashMap,
        // 	WeightedCFLOBDDComplexFloatBoostMulNodeHandle c1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle c2, int c1_zero_index, int c2_zero_index);
    }
 }

#endif

