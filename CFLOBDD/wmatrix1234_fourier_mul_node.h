#ifndef W_MATRIX1234_FOURIER_MUL_NODE_GUARD
#define W_MATRIX1234_FOURIER_MUL_NODE_GUARD

#include <map>
#include <unordered_map>
#include "weighted_cflobdd_node_t.h"
#include "return_map_T.h"
#include "weighted_matmult_map.h"
#include "fourier_semiring.h"

namespace CFL_OBDD {

    namespace WeightedMatrix1234FourierMul {

        typedef WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierMulNodeHandle;
        typedef WeightedCFLOBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierInternalNode;
        typedef WeightedCFLOBDDDontCareNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierDontCareNode;
        typedef WeightedCFLOBDDForkNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierForkNode;
        typedef WeightedCFLOBDDLeafNode<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDFourierLeafNode;
        typedef WConnection<fourierSemiring, std::multiplies<fourierSemiring>> Connection;
	    typedef ReturnMapHandle<WeightedMatMultMapHandle<fourierSemiring>> CFLOBDDMatMultMapHandle;
        typedef std::tuple<WeightedCFLOBDDFourierMulNodeHandle, CFLOBDDMatMultMapHandle, fourierSemiring> MatMultReturnT;

        class WeightedMatMultPair{
        public:
            WeightedCFLOBDDFourierMulNodeHandle m1;
            WeightedCFLOBDDFourierMulNodeHandle m2;
            WeightedMatMultPair(WeightedCFLOBDDFourierMulNodeHandle p1, WeightedCFLOBDDFourierMulNodeHandle p2)
            {
                m1 = p1;
                m2 = p2;
            }

            struct MatMultPairHash {
                size_t operator()(const WeightedMatMultPair& p) const
                {
                    WeightedCFLOBDDFourierMulNodeHandle t1 = p.m1;
                    WeightedCFLOBDDFourierMulNodeHandle t2 = p.m2;
                    auto hash1 = t1.Hash(997);
                    auto hash2 = t2.Hash(997);
                    return 117 * (hash1 + 1) + hash2;
                }
            };

            bool operator==(const WeightedMatMultPair& p) const
            {
                WeightedCFLOBDDFourierMulNodeHandle m11 = m1;
                WeightedCFLOBDDFourierMulNodeHandle m12 = m2;
                WeightedCFLOBDDFourierMulNodeHandle m21 = p.m1;
                WeightedCFLOBDDFourierMulNodeHandle m22 = p.m2;
                return (m11 == m21) && (m12 == m22);
            }
        };

        extern WeightedCFLOBDDFourierMulNodeHandle MkIdRelationInterleavedNode(unsigned int level);
        extern WeightedCFLOBDDFourierMulNodeHandle MkWalshInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFourierMulNodeHandle MkInverseReedMullerInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFourierMulNodeHandle MkNegationMatrixInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFourierMulNodeHandle MkCNOTInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFourierMulNodeHandle MkExchangeInterleavedNode(unsigned int i);
        extern WeightedCFLOBDDFourierMulNodeHandle MkCNOTNode(unsigned int level, unsigned int n, long int controller, long int controlled);
        extern WeightedCFLOBDDFourierMulNodeHandle MkCPGateNode(std::unordered_map<std::string, WeightedCFLOBDDFourierMulNodeHandle>& cp_hashMap, unsigned int level, long int controller, long int controlled, fourierSemiring theta_val);
        extern WeightedCFLOBDDFourierMulNodeHandle MkSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num);
        extern WeightedCFLOBDDFourierMulNodeHandle MkCSwapGateNode(unsigned int level, long int controller, long int i, long int j, int case_num);
        extern WeightedCFLOBDDFourierMulNodeHandle MkRZGateNode(unsigned int level, fourierSemiring theta);
        extern WeightedCFLOBDDFourierMulNodeHandle MkCADDGateNode(unsigned int level, long int c, long int x, WeightedCFLOBDDFourierMulNodeHandle f, int flag);
        extern WeightedCFLOBDDFourierMulNodeHandle MkCADDGate2Node(unsigned int level, long int x, WeightedCFLOBDDFourierMulNodeHandle f, int flag);
        extern bool CheckIfIndexIsNonZeroNode(unsigned int level, int index, WeightedCFLOBDDFourierMulNodeHandle f, int flag);
        extern WeightedCFLOBDDFourierMulNodeHandle MkSetBToZeroNode(unsigned int level, WeightedCFLOBDDFourierMulNodeHandle f, int flag);
        extern std::pair<WeightedCFLOBDDFourierMulNodeHandle, int> ComputeIQFTNode(unsigned int level, WeightedCFLOBDDFourierMulNodeHandle f, BIG_INT N, int n, int flag, int b_vars, int actual_level);
        extern std::pair<WeightedCFLOBDDFourierMulNodeHandle, int> MeasureAndResetNode(unsigned int level, long int n, WeightedCFLOBDDFourierMulNodeHandle f, fourierSemiring R);
        extern WeightedCFLOBDDFourierMulNodeHandle ResetStateNode(unsigned int level, WeightedCFLOBDDFourierMulNodeHandle f);

        extern WeightedCFLOBDDFourierMulNodeHandle KroneckerProduct2VocsNode(WeightedCFLOBDDFourierMulNodeHandle m1, WeightedCFLOBDDFourierMulNodeHandle m2, int zero_index_m1, int zero_index_m2); 
        extern std::tuple<WeightedCFLOBDDFourierMulNodeHandle, CFLOBDDMatMultMapHandle, fourierSemiring>
		MatrixMultiplyV4Node(WeightedCFLOBDDFourierMulNodeHandle c1, WeightedCFLOBDDFourierMulNodeHandle c2, int zero_exit_1, int zero_exit_2);
        
        // Initialization routine that needs to be called before any call to MatrixProjectVoc23Node
        extern void Matrix1234InitializerNode();  // Empty for now

        // extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4Node(
        // 	std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash>& hashMap,
        // 	WeightedCFLOBDDFourierMulNodeHandle c1, WeightedCFLOBDDFourierMulNodeHandle c2);
        // extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4WithInfoNode(
        // 	std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash>& hashMap,
        // 	WeightedCFLOBDDFourierMulNodeHandle c1, WeightedCFLOBDDFourierMulNodeHandle c2, int c1_zero_index, int c2_zero_index);
    }
 }

#endif

