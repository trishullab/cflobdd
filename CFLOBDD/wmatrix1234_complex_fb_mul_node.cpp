#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <vector>
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"
#include "wmatrix1234_complex_fb_mul_node.h"
#include "weighted_cross_product.h"
#include "return_map_T.h"
#include "weighted_matmult_map.h"
#include "wmatrix1234_complex_fb_mul_bdd_node.h"

namespace CFL_OBDD {

    namespace WeightedMatrix1234ComplexFloatBoostMul
    {
        std::vector<ReturnMapHandle<int>> commonly_used_return_maps;// m0, m1, m01, m10
        std::unordered_map<WeightedMatMultPair, MatMultReturnT, WeightedMatMultPair::MatMultPairHash> matmult_hash;
        std::unordered_map<std::string, WeightedCFLOBDDComplexFloatBoostMulNodeHandle> cnot_hashMap;
	    std::unordered_map<std::string, WeightedCFLOBDDComplexFloatBoostMulNodeHandle> swap_hashMap;

        void InitReturnMapHandles(){
            ReturnMapHandle<int> m0, m1, m01, m10;
            m0.AddToEnd(0);
            m0.Canonicalize();
            m1.AddToEnd(1);
            m1.Canonicalize();
            m01.AddToEnd(0);
            m01.AddToEnd(1);
            m01.Canonicalize();
            m10.AddToEnd(1);
            m10.AddToEnd(0);
            m10.Canonicalize();
            commonly_used_return_maps.push_back(m0);
            commonly_used_return_maps.push_back(m1);
            commonly_used_return_maps.push_back(m01);
            commonly_used_return_maps.push_back(m10);
        }

        void Matrix1234InitializerNode()
        {
            // Empty for now
            //projectMemoTable = new CFLOBDDLinearMapMemoTable;
            InitReturnMapHandles(); 
            return;
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkIdRelationInterleavedNode(unsigned int level, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkIdRelationInterleavedNode(bn->numberOfVars);
                bn->numberOfVars = level;
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            }

            WeightedCFLOBDDComplexFloatBoostInternalNode *n;

            if (level == 0) {
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle;
            }
            return WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level];
    //         else if (level == 1) {
    //             n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);

    //             WeightedCFLOBDDComplexFloatBoostMulNodeHandle temp = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle;
    //             n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
    //             n->numBConnections = 2;
    //             n->BConnection = new Connection[n->numBConnections];
    //             WeightedCFLOBDDComplexFloatBoostMulNodeHandle b0 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(1, 0)); 
    //             WeightedCFLOBDDComplexFloatBoostMulNodeHandle b1 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(0, 1));
    //             n->BConnection[0] = Connection(b0, commonly_used_return_maps[2]);//m01
    //             n->BConnection[1] = Connection(b1, commonly_used_return_maps[3]);//m10
    //         }
    //         else {  // Create an appropriate CFLOBDDInternalNode
    //             n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);

    //             WeightedCFLOBDDComplexFloatBoostMulNodeHandle temp = MkIdRelationInterleavedNode(level - 1);
    //             n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
    //             n->numBConnections = 2;
    //             n->BConnection = new Connection[n->numBConnections];
    //             n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
    //             n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps[1]);//m1
    //         }
    //         n->numExits = 2;
    // #ifdef PATH_COUNTING_ENABLED
    //         n->InstallPathCounts();
    // #endif
    //         return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
        } // MkIdRelationInterleavedNode
    

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkWalshInterleavedNode(unsigned int i, int cflobdd_kind)
        {
            assert(i >= 1);
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(i);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkWalshInterleavedNode(bn->numberOfVars);
                bn->numberOfVars = i;
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            }
            WeightedCFLOBDDComplexFloatBoostInternalNode *n = new WeightedCFLOBDDComplexFloatBoostInternalNode(i);
            if (i == 1) {  // Base case

                WeightedCFLOBDDComplexFloatBoostMulNodeHandle temp = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(1, 1));
                CFLOBDDReturnMapHandle m01;
                m01.AddToEnd(0);
                m01.AddToEnd(1);
                m01.Canonicalize();
                n->AConnection = Connection(temp, m01);
                n->numBConnections = 2;
                n->BConnection = new Connection[n->numBConnections];
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle b0 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostDontCareNode(1, 1));
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle b1 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostDontCareNode(1, -1));
                CFLOBDDReturnMapHandle m0;
                m0.AddToEnd(0);
                m0.Canonicalize();
                n->BConnection[0] = Connection(b0, m0);
                n->BConnection[1] = Connection(b1, m0);
            }
            else {
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle temp = MkWalshInterleavedNode(i - 1, cflobdd_kind);
                CFLOBDDReturnMapHandle m0;
                m0.AddToEnd(0);
                m0.Canonicalize();
                n->AConnection = Connection(temp, m0);

                n->numBConnections = 1;
                n->BConnection = new Connection[n->numBConnections];
                n->BConnection[0] = Connection(temp, m0);
            }
            n->numExits = 1;
    #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
    #endif
            return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
        } // MkWalshInterleavedNode

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle KroneckerProduct2VocsNode(WeightedCFLOBDDComplexFloatBoostMulNodeHandle m1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle m2, 
            int zero_index_m1, int zero_index_m2)
        {
            if (m1.handleContents->NodeKind() == W_BDD_TOPNODE)
            {
                WeightedBDDComplexFloatBoostTopNode* c1 = (WeightedBDDComplexFloatBoostTopNode *)m1.handleContents;
                WeightedBDDComplexFloatBoostTopNode* c2 = (WeightedBDDComplexFloatBoostTopNode *)m2.handleContents;
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(c1->numberOfVars + c2->numberOfVars);
                std::unordered_map<WeightedMatrix1234BDDComplexFloatBoostMul::WeightedBDDComplexFloatBoostMulNodeHandle, WeightedMatrix1234BDDComplexFloatBoostMul::WeightedBDDComplexFloatBoostMulNodeHandle, WeightedMatrix1234BDDComplexFloatBoostMul::WeightedBDDComplexFloatBoostMulNodeHandle::WeightedBDDNodeHandle_Hash> hashMap;
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::KroneckerProduct2VocsNode(hashMap, c1->bddContents, c2->bddContents, c1->numberOfVars);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn); 
            }
            int level = m1.handleContents->level;
            WeightedCFLOBDDComplexFloatBoostInternalNode* n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level + 1);
            CFLOBDDReturnMapHandle m;
            for (int i = 0; i < m1.handleContents->numExits; i++)
                m.AddToEnd(i);
            m.Canonicalize();
            n->AConnection =  Connection(m1, m);
            n->numBConnections = m1.handleContents->numExits;
            n->BConnection = new Connection[n->numBConnections];
            if (n->numBConnections == 1)
            {
                CFLOBDDReturnMapHandle mB;
                for (int i = 0; i < m2.handleContents->numExits; i++)
                    mB.AddToEnd(i);
                mB.Canonicalize();
                n->BConnection[0] = Connection(m2, mB);
                n->numExits = m2.handleContents->numExits;
            }
            else if (zero_index_m1 == 0)
            {
                CFLOBDDReturnMapHandle c0;
                c0.AddToEnd(0);
                c0.Canonicalize();
                n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level], c0);
                CFLOBDDReturnMapHandle c;
                if (zero_index_m2 == -1)
                    c.AddToEnd(1);
                else if (zero_index_m2 == 0)
                {
                    c.AddToEnd(0); c.AddToEnd(1);
                }
                else
                {
                    c.AddToEnd(1); c.AddToEnd(0);
                }
                c.Canonicalize();
                n->BConnection[1] = Connection(m2, c);
                n->numExits = 2;
            }
            else if (zero_index_m1 == 1)
            {
                CFLOBDDReturnMapHandle mB;
                CFLOBDDReturnMapHandle c;
                for (int i = 0; i < m2.handleContents->numExits; i++)
                    mB.AddToEnd(i);
                mB.Canonicalize();
                n->BConnection[0] = Connection(m2, mB);
                if (zero_index_m2 == -1)
                    c.AddToEnd(1);
                else if (zero_index_m2 == 0)
                    c.AddToEnd(0);
                else
                    c.AddToEnd(1);
                c.Canonicalize();
                n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level], c);
                n->numExits = 2; 
            }
        #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
        #endif
            return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
        }
    

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkExchangeInterleavedNode(unsigned int level)
        {
            WeightedCFLOBDDComplexFloatBoostInternalNode *n;

            if (level == 0) {
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01;
            }
            else if (level == 1) {
                n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);

                n->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
                n->numBConnections = 2;
                n->BConnection = new Connection[n->numBConnections];
                n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, commonly_used_return_maps[2]);//m01
                n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, commonly_used_return_maps[3]);//m10
            }
            else {  // Create an appropriate CFLOBDDInternalNode
                n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);

                WeightedCFLOBDDComplexFloatBoostMulNodeHandle temp = MkExchangeInterleavedNode(level - 1);
                n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
                n->numBConnections = 2;
                n->BConnection = new Connection[n->numBConnections];
                CFLOBDDReturnMapHandle m0, m01;
                m0.AddToEnd(0); m0.Canonicalize();
                m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], m0);//m01
                n->BConnection[1] = Connection(temp, m01);//m1
            }
            n->numExits = 2;
    #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
    #endif
            return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
        } // MkIdRelationInterleavedNode

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkNegationMatrixInterleavedNode(unsigned int level, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkNegationMatrixInterleavedNode(bn->numberOfVars);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            }
            // ADD code for Weighted CFLOBDD 
        } // MkIdRelationInterleavedNode
    


        std::tuple<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_COMPLEX_FLOAT>
        Add(WeightedCFLOBDDComplexFloatBoostMulNodeHandle b1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle b2,
            CFLOBDDMatMultMapHandle b1_m, CFLOBDDMatMultMapHandle b2_m, BIG_COMPLEX_FLOAT b_f1, BIG_COMPLEX_FLOAT b_f2)
        {   
            WeightedPairProductMapHandle<BIG_COMPLEX_FLOAT> MapHandle;
            auto tmp =  PairProduct2(b1, b2, b_f1, b_f2, MapHandle);

            CFLOBDDMatMultMapHandle returnMapHandle;
        
            boost::unordered_map<size_t, unsigned int> reductionMap;
            ReductionMapHandle reductionMapHandle;
            unsigned int iterator = 0;
            WeightedValuesListHandle<BIG_COMPLEX_FLOAT> valList;
            while (iterator < MapHandle.Size()){
                WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> c1, c2;
                int first, second;
                first = MapHandle[iterator].first.First();
                second = MapHandle[iterator].first.Second();
                WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> val;
                if (first == -1 && second == -1){
                    val = WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>();
                    BIG_COMPLEX_FLOAT v(0.0);
                    val.Add(std::make_pair(-1,-1), v);
                    val.mapContents->contains_zero_val = true;
                }
                else {
                    c1 = b1_m.Lookup(first);
                    c2 = b2_m.Lookup(second);
                    BIG_COMPLEX_FLOAT v1, v2;
                    v1 = MapHandle[iterator].second.First();
                    v2 = MapHandle[iterator].second.Second();
                    val = v1 * c1 + v2 * c2;
                }
                val.Canonicalize();
                if (reductionMap.find(val.getHashCheck()) == reductionMap.end()){
                    returnMapHandle.AddToEnd(val);
                    reductionMap.insert(std::make_pair(val.getHashCheck(), returnMapHandle.Size() - 1)); 
                    reductionMapHandle.AddToEnd(returnMapHandle.Size() - 1);
                    if (val.mapContents->contains_zero_val == false)
                        valList.AddToEnd(1.0);
                    else
                        valList.AddToEnd(0.0);
                }
                else{
                    reductionMapHandle.AddToEnd(reductionMap[val.getHashCheck()]);
                    if (val.mapContents->contains_zero_val == false)
                        valList.AddToEnd(1.0);
                    else
                        valList.AddToEnd(0.0);
                }
                iterator++;
            }

            returnMapHandle.Canonicalize();
            reductionMapHandle.Canonicalize();
            valList.Canonicalize();
            auto reduced_n = tmp.Reduce(reductionMapHandle, returnMapHandle.Size(), valList, false);
            BIG_COMPLEX_FLOAT factor = reduced_n.second;

            return std::make_tuple(reduced_n.first, returnMapHandle, factor);
        }

        MatMultReturnT
		MatrixMultiplyV4Node(WeightedCFLOBDDComplexFloatBoostMulNodeHandle c1, WeightedCFLOBDDComplexFloatBoostMulNodeHandle c2, int zero_exit_1, int zero_exit_2)
        {
            if (c1.handleContents->NodeKind() == W_BDD_TOPNODE)
            {
                WeightedBDDComplexFloatBoostTopNode* m1 = (WeightedBDDComplexFloatBoostTopNode *)c1.handleContents;
                WeightedBDDComplexFloatBoostTopNode* m2 = (WeightedBDDComplexFloatBoostTopNode *)c2.handleContents;
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(m1->numberOfVars);
                auto tmp = WeightedMatrix1234BDDComplexFloatBoostMul::MatrixMultiplyV4Node(m1->bddContents, m2->bddContents, m1->numberOfVars);
                bn->bddContents = std::get<0>(tmp);
                CFLOBDDMatMultMapHandle mt;
                WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> m_tmp;
                BIG_COMPLEX_FLOAT zero_0 = 0, one_1 = 1;
                for (auto i : std::get<1>(tmp)){
                    if (i == 0)
                        m_tmp.Add(std::make_pair(0, 0), zero_0);
                    else
                        m_tmp.Add(std::make_pair(1, 1), one_1);
                }
                m_tmp.Canonicalize();
                mt.AddToEnd(m_tmp);
                mt.Canonicalize();
                return std::make_tuple(WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn), mt, std::get<2>(tmp));
            }
            auto mmp = WeightedMatMultPair(c1, c2);
            if (matmult_hash.find(mmp) != matmult_hash.end())
                return matmult_hash[mmp];
            
            BIG_COMPLEX_FLOAT zero(0.0, 0.0);
            if (c1 == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1.handleContents->level] || c2 == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c2.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> tmp;
                tmp.Add(std::make_pair(-1,-1), zero);
                tmp.mapContents->contains_zero_val = true;
                tmp.Canonicalize();
                m.AddToEnd(tmp);
                m.Canonicalize();
                return std::make_tuple(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1.handleContents->level], m, zero);
            }

            if (c1 == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[c1.handleContents->level] && c2 == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[c2.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                BIG_COMPLEX_FLOAT one = 1.0;
                WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> m0, m1;
                m0.Add(std::make_pair(0,0), one); m0.Canonicalize();
                m1.Add(std::make_pair(-1,-1), zero); m1.mapContents->contains_zero_val = true; m1.Canonicalize();
                m.AddToEnd(m0); m.AddToEnd(m1); m.Canonicalize();
                return std::make_tuple(c1, m, one);
            }

            if (c1 == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[c1.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                BIG_COMPLEX_FLOAT one = 1.0;
                for (int i = 0; i < c2.handleContents->numExits; i++)
                {
                    WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> tmp;
                    if (i != zero_exit_2)
                        tmp.Add(std::make_pair(0, i), one);
                    else{
                        tmp.Add(std::make_pair(-1, -1), zero);
                        tmp.mapContents->contains_zero_val = true;
                    }
                    tmp.Canonicalize();
                    m.AddToEnd(tmp);
                }
                m.Canonicalize();
                return std::make_tuple(c2, m, one);
            }
            if (c2 == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[c1.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                BIG_COMPLEX_FLOAT one = 1.0;
                for (int i = 0; i < c1.handleContents->numExits; i++)
                {
                    WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> tmp;
                    if (i != zero_exit_1)
                        tmp.Add(std::make_pair(i, 0), one);
                    else{
                        tmp.Add(std::make_pair(-1, -1), zero);
                        tmp.mapContents->contains_zero_val = true;
                    }
                    tmp.Canonicalize();
                    m.AddToEnd(tmp);
                }
                m.Canonicalize();
                return std::make_tuple(c1, m, one);
            }

            CFLOBDDMatMultMapHandle g_return_map;
            WeightedValuesListHandle<BIG_COMPLEX_FLOAT> valList;
            WeightedCFLOBDDComplexFloatBoostInternalNode* g = new WeightedCFLOBDDComplexFloatBoostInternalNode(c1.handleContents->level);
            ReductionMapHandle reductionMapHandle;
            BIG_COMPLEX_FLOAT top_factor = 1.0;
            
            if (c1.handleContents->level == 1){
                WeightedCFLOBDDComplexFloatBoostInternalNode* c1_internal = (WeightedCFLOBDDComplexFloatBoostInternalNode *)c1.handleContents;
                WeightedCFLOBDDComplexFloatBoostInternalNode* c2_internal = (WeightedCFLOBDDComplexFloatBoostInternalNode *)c2.handleContents;

                /*
                    [[a_0, b_0],    [[a_1, b_1],
                     [c_0, d_0]]     [c_1, d_1]]
                */
                BIG_COMPLEX_FLOAT a0, b0, c0, d0, a1, b1, c1, d1;
                WeightedCFLOBDDComplexFloatBoostLeafNode* M1_A = (WeightedCFLOBDDComplexFloatBoostLeafNode *)c1_internal->AConnection.entryPointHandle->handleContents;
                WeightedCFLOBDDComplexFloatBoostLeafNode* M2_A = (WeightedCFLOBDDComplexFloatBoostLeafNode *)c2_internal->AConnection.entryPointHandle->handleContents;
                a0 = b0 = M1_A->lweight;
                c0 = d0 = M1_A->rweight;
                a1 = b1 = M2_A->lweight;
                c1 = d1 = M2_A->rweight;
                WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> v1, v2, v3, v4;
                int M1_numB = c1_internal->numBConnections;
                int M2_numB = c2_internal->numBConnections;
                WeightedCFLOBDDComplexFloatBoostLeafNode* M1_b0 = (WeightedCFLOBDDComplexFloatBoostLeafNode *)c1_internal->BConnection[0].entryPointHandle->handleContents;
                WeightedCFLOBDDComplexFloatBoostLeafNode* M1_b1 = (WeightedCFLOBDDComplexFloatBoostLeafNode *)c1_internal->BConnection[M1_numB-1].entryPointHandle->handleContents;
                a0 = a0 * M1_b0->lweight;
                b0 = b0 * M1_b0->rweight;
                c0 = c0 * M1_b1->lweight;
                d0 = d0 * M1_b1->rweight;
                WeightedCFLOBDDComplexFloatBoostLeafNode* M2_b0 = (WeightedCFLOBDDComplexFloatBoostLeafNode *)c2_internal->BConnection[0].entryPointHandle->handleContents;
                WeightedCFLOBDDComplexFloatBoostLeafNode* M2_b1 = (WeightedCFLOBDDComplexFloatBoostLeafNode *)c2_internal->BConnection[M2_numB-1].entryPointHandle->handleContents;
                a1 = a1 * M2_b0->lweight;
                b1 = b1 * M2_b0->rweight;
                c1 = c1 * M2_b1->lweight;
                d1 = d1 * M2_b1->rweight;

                int M1_b0_numE = M1_b0->numExits;
                int M1_b1_numE = M1_b1->numExits;
                int M2_b0_numE = M2_b0->numExits;
                int M2_b1_numE = M2_b1->numExits;

                // v1
                BIG_COMPLEX_FLOAT a0a1 = a0 * a1;
                BIG_COMPLEX_FLOAT b0c1 = b0 * c1;
                if (a0a1 == 0 && b0c1 == 0){
                    v1.Add(std::make_pair(-1,-1), zero);
                    v1.mapContents->contains_zero_val = true;
                }
                if (b0c1 != 0)
                    v1.Add(std::make_pair(c1_internal->BConnection[0].returnMapHandle[M1_b0_numE-1], c2_internal->BConnection[M2_numB-1].returnMapHandle[0]), b0c1);
                if (a0a1 != 0)
                    v1.Add(std::make_pair(c1_internal->BConnection[0].returnMapHandle[0], c2_internal->BConnection[0].returnMapHandle[0]), a0a1);
                v1.Canonicalize();

                // v2
                BIG_COMPLEX_FLOAT a0b1 = a0 * b1;
                BIG_COMPLEX_FLOAT b0d1 = b0 * d1;
                if (a0b1 == 0 && b0d1 == 0){
                    v2.Add(std::make_pair(-1,-1), zero);
                    v2.mapContents->contains_zero_val = true;
                }
                if (b0d1 != 0)
                    v2.Add(std::make_pair(c1_internal->BConnection[0].returnMapHandle[M1_b0_numE-1], c2_internal->BConnection[M2_numB-1].returnMapHandle[M2_b1_numE-1]), b0d1);
                if (a0b1 != 0)
                    v2.Add(std::make_pair(c1_internal->BConnection[0].returnMapHandle[0], c2_internal->BConnection[0].returnMapHandle[M2_b0_numE-1]), a0b1);
                v2.Canonicalize();

                // v3
                BIG_COMPLEX_FLOAT c0a1 = c0 * a1;
                BIG_COMPLEX_FLOAT d0c1 = d0 * c1;
                if (c0a1 == 0 && d0c1 == 0){
                    v3.Add(std::make_pair(-1,-1), zero);
                    v3.mapContents->contains_zero_val = true;
                }
                if (d0c1 != 0)
                    v3.Add(std::make_pair(c1_internal->BConnection[M1_numB-1].returnMapHandle[M1_b1_numE-1], c2_internal->BConnection[M2_numB-1].returnMapHandle[0]), d0c1);
                if (c0a1 != 0)
                    v3.Add(std::make_pair(c1_internal->BConnection[M1_numB-1].returnMapHandle[0], c2_internal->BConnection[0].returnMapHandle[0]), c0a1);
                v3.Canonicalize();

                // v4
                BIG_COMPLEX_FLOAT c0b1 = c0 * b1;
                BIG_COMPLEX_FLOAT d0d1 = d0 * d1;
                if (c0b1 == 0 && d0d1 == 0){
                    v4.Add(std::make_pair(-1,-1), zero);
                    v4.mapContents->contains_zero_val = true;
                }
                if (d0d1 != 0)
                    v4.Add(std::make_pair(c1_internal->BConnection[M1_numB-1].returnMapHandle[M1_b1_numE-1], c2_internal->BConnection[M2_numB-1].returnMapHandle[M2_b1_numE-1]), d0d1);
                if (c0b1 != 0)
                    v4.Add(std::make_pair(c1_internal->BConnection[M1_numB-1].returnMapHandle[0], c2_internal->BConnection[0].returnMapHandle[M2_b0_numE-1]), c0b1);
                v4.Canonicalize();

                if (v1 == v3 && v2 == v4){
                    CFLOBDDReturnMapHandle m0;
                    m0.AddToEnd(0); m0.Canonicalize();
                    if (v1 == v2){
                        g->numBConnections = 1;
                        g->BConnection = new Connection[g->numBConnections];
                        // if (v1.mapContents->contains_zero_val == true){
                        //     g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                        // }
                        // else {
                            g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0);
                            g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0); 
                        // }
                        reductionMapHandle.AddToEnd(0);
                        valList.AddToEnd(v1.mapContents->contains_zero_val ? 0 : 1);
                        g_return_map.AddToEnd(v1);
                    }
                    else{
                        g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0);
                        g->numBConnections = 1;
                        g->BConnection = new Connection[g->numBConnections];
                        CFLOBDDReturnMapHandle m01;
                        m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                        // if (v1.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                        // else if (v2.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                        // else
                            g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01); 
                        reductionMapHandle.AddToEnd(0);
                        reductionMapHandle.AddToEnd(1);
                        valList.AddToEnd(v1.mapContents->contains_zero_val ? 0 : 1);
                        valList.AddToEnd(v2.mapContents->contains_zero_val ? 0 : 1);
                        g_return_map.AddToEnd(v1); 
                        g_return_map.AddToEnd(v2); 
                    }
                }
                else
                {
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    int la = 1, ra = 1;
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    if (v1 == v2)
                    {
                        CFLOBDDReturnMapHandle m0;
                        m0.AddToEnd(0); m0.Canonicalize();
                        // if (v1.mapContents->contains_zero_val == true){
                        //     la = 0;
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                        // }
                        // else{
                        //     la = 1;
                            g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0);
                        // }
                        g_return_map.AddToEnd(v1);
                        reductionMapHandle.AddToEnd(0);
                        valList.AddToEnd(v1.mapContents->contains_zero_val ? 0 : 1);
                    }
                    else{
                        la = 1;
                        // if (v1.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                        // else if (v2.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                        // else
                            g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                        g_return_map.AddToEnd(v1);
                        g_return_map.AddToEnd(v2);
                        reductionMapHandle.AddToEnd(0);
                        reductionMapHandle.AddToEnd(1);
                        valList.AddToEnd(v1.mapContents->contains_zero_val ? 0 : 1);
                        valList.AddToEnd(v2.mapContents->contains_zero_val ? 0 : 1);
                    }
                    if (v3 == v4)
                    {
                        int k = 0;
                        for (k = 0; k < g_return_map.Size(); k++)
                        {
                            if (g_return_map[k] == v3)
                                break;
                        }
                        CFLOBDDReturnMapHandle m;
                        m.AddToEnd(k); m.Canonicalize();
                        // if (v3.mapContents->contains_zero_val == true)
                        // {
                        //     ra = 0;
                        //     g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m);
                        // }
                        // else {
                        //     ra = 1;
                            g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m);
                        // }
                        if (k >= g_return_map.Size()){
                            g_return_map.AddToEnd(v3);
                            reductionMapHandle.AddToEnd(k);
                            valList.AddToEnd(v3.mapContents->contains_zero_val ? 0 : 1);
                        }
                    }
                    else
                    {
                        ra = 1;
                        int k1 = g_return_map.Size(), k2 = -1, k = 0;
                        for (k = 0; k < g_return_map.Size(); k++)
                        {
                            if (g_return_map[k] == v3){
                                k1 = k;
                            }
                            else if (g_return_map[k] == v4){
                                k2 = k;
                            }
                        }
                        if (k2 == -1){
                            k2 = std::max(k1, k-1) + 1;
                        }
                        CFLOBDDReturnMapHandle m1;
                        m1.AddToEnd(k1); m1.AddToEnd(k2); m1.Canonicalize();
                        // if (v3.mapContents->contains_zero_val == true)
                        //     g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m1);
                        // else if (v4.mapContents->contains_zero_val == true)
                        //     g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m1);
                        // else
                            g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m1);
                        if (k1 >= k){
                            g_return_map.AddToEnd(v3);
                            reductionMapHandle.AddToEnd(k1);
                            valList.AddToEnd(v3.mapContents->contains_zero_val ? 0 : 1);
                        }
                        if (k2 >= k){
                            g_return_map.AddToEnd(v4);
                            reductionMapHandle.AddToEnd(k2);
                            valList.AddToEnd(v4.mapContents->contains_zero_val ? 0 : 1);
                        }
                    }
                    // if (la == 1 && ra == 1)
                    //     g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    // else if (la == 1 && ra == 0)
                    //     g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    // else if (la == 0 && ra == 1)
                    //     g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                    // assert(!(la == 0 && ra == 0));
                }

            }
            else{
                
                WeightedCFLOBDDComplexFloatBoostInternalNode* c1_internal = (WeightedCFLOBDDComplexFloatBoostInternalNode *)c1.handleContents;
                WeightedCFLOBDDComplexFloatBoostInternalNode* c2_internal = (WeightedCFLOBDDComplexFloatBoostInternalNode *)c2.handleContents;
                int level = c1.handleContents->level;
                // populate zero exit information
                int a_zero_exit_1 = -1, a_zero_exit_2 = -1;
                if (zero_exit_1 != -1){
                    for (int i = 0; i < c1_internal->numBConnections; i++)
                    {
                        if (*(c1_internal->BConnection[i].entryPointHandle) == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1])
                        {
                            a_zero_exit_1 = i;
                            break;
                        }
                    }
                }
                if (zero_exit_2 != -1)
                {
                    for (int i = 0; i < c2_internal->numBConnections; i++)
                    {
                        if (*(c2_internal->BConnection[i].entryPointHandle) == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1])
                        {
                            a_zero_exit_2 = i;
                            break;
                        }
                    }
                }
                auto aa = MatrixMultiplyV4Node(*(c1_internal->AConnection.entryPointHandle),
                    *(c2_internal->AConnection.entryPointHandle), a_zero_exit_1, a_zero_exit_2);
                CFLOBDDReturnMapHandle mI;
                for (unsigned int i = 0; i < std::get<1>(aa).Size(); i++)
                    mI.AddToEnd(i);
                mI.Canonicalize();
                top_factor = std::get<2>(aa);
                g->AConnection = Connection(std::get<0>(aa), mI);
                g->numBConnections = mI.Size();
                g->BConnection = new Connection[g->numBConnections];
                g->numExits = 0;
                std::unordered_map<unsigned int, unsigned int> mapFromHandleToIndex;
                for (unsigned int i = 0; i < g->numBConnections; i++){
                    WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> matmult_returnmap = std::get<1>(aa)[i];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle ans;
                    CFLOBDDMatMultMapHandle ans_matmult_map;
                    BIG_COMPLEX_FLOAT ans_factor = 1.0;
                    bool first = true;
                    // Consider Multiplication of M1 and M2
                    for (auto &v : matmult_returnmap.mapContents->map){
                        unsigned int M1_index = v.first.first;
                        unsigned int M2_index = v.first.second;
                        std::tuple<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_COMPLEX_FLOAT> bb_old;
                        CFLOBDDMatMultMapHandle new_bb_return;
                        WeightedCFLOBDDComplexFloatBoostMulNodeHandle n1, n2;
                        CFLOBDDReturnMapHandle n1_return, n2_return;
                        int b_zero_exit_1 = -1, b_zero_exit_2 = -1;
                        if (M1_index == -1 || v.second == 0)
                        {
                            n1 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1_internal->level-1];
                            n1_return.AddToEnd(0);
                            n1_return.Canonicalize();
                            b_zero_exit_1 = 0;
                        }
                        else{
                            n1 = *(c1_internal->BConnection[M1_index].entryPointHandle);
                            n1_return = c1_internal->BConnection[M1_index].returnMapHandle;
                            b_zero_exit_1 = c1_internal->BConnection[M1_index].returnMapHandle.LookupInv(zero_exit_1);
                        }
                        if (M2_index == -1 || v.second == 0)
                        {
                            n2 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1_internal->level-1];
                            n2_return.AddToEnd(0);
                            n2_return.Canonicalize();
                            b_zero_exit_2 = 0;
                        }
                        else{
                            n2 = *(c2_internal->BConnection[M2_index].entryPointHandle);
                            n2_return = c2_internal->BConnection[M2_index].returnMapHandle;
                            b_zero_exit_2 = c2_internal->BConnection[M2_index].returnMapHandle.LookupInv(zero_exit_2);
                        }
                        bb_old = MatrixMultiplyV4Node(n1,n2, b_zero_exit_1, b_zero_exit_2);
                        CFLOBDDMatMultMapHandle old_bb_return = std::get<1>(bb_old);
                        for (unsigned int j = 0; j < old_bb_return.Size(); j++)
                        {
                            WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> tmp;
                            for (auto& it : old_bb_return[j].mapContents->map){
                                if (it.first.first == -1 && it.first.second == -1)
                                    tmp.Add(std::make_pair(-1,-1), it.second);
                                else
                                    tmp.Add(std::make_pair(n1_return[it.first.first],
                                        n2_return[it.first.second]), it.second);
                            }
                            tmp.Canonicalize();
                            new_bb_return.AddToEnd(tmp);
                        }
                        new_bb_return.Canonicalize();
                        if (first){
                            ans = std::get<0>(bb_old);
                            ans_matmult_map = new_bb_return;
                            ans_factor = v.second * std::get<2>(bb_old);
                            first = false;
                        }
                        else{
                            // TODO: change this
                            auto t = Add(ans, std::get<0>(bb_old), ans_matmult_map, new_bb_return, ans_factor, v.second * std::get<2>(bb_old));
                            ans = std::get<0>(t);
                            ans_matmult_map = std::get<1>(t);
                            ans_factor = std::get<2>(t);
                        }
                    }

                    CFLOBDDReturnMapHandle ans_return_map;
                    for (unsigned int j = 0; j < ans_matmult_map.Size(); j++){
                        unsigned int map_hash_check = ans_matmult_map[j].getHashCheck();
                        if (mapFromHandleToIndex.find(map_hash_check) == mapFromHandleToIndex.end()){
                            ans_return_map.AddToEnd(g->numExits++);
                            g_return_map.AddToEnd(ans_matmult_map[j]);
                            reductionMapHandle.AddToEnd(g_return_map.Size() - 1);
                            mapFromHandleToIndex[map_hash_check] = g_return_map.Size() - 1;
                            if (ans_matmult_map[j].mapContents->contains_zero_val == true)
                                valList.AddToEnd(0.0);
                            else
                                valList.AddToEnd(ans_factor);
                        }
                        else{
                            unsigned int index = mapFromHandleToIndex[map_hash_check];
                            ans_return_map.AddToEnd(g->numExits++);
                            // ans_return_map.AddToEnd(index);
                            reductionMapHandle.AddToEnd(index);
                            if (ans_matmult_map[j].mapContents->contains_zero_val == true)
                                valList.AddToEnd(0.0);
                            else
                                valList.AddToEnd(ans_factor); 
                        }
                    }
                    ans_return_map.Canonicalize();
                    g->BConnection[i] = Connection(ans, ans_return_map);
                }

            }
            
            g_return_map.Canonicalize();
            reductionMapHandle.Canonicalize();
            valList.Canonicalize();
            g->numExits = g_return_map.Size();

    #ifdef PATH_COUNTING_ENABLED
            g->InstallPathCounts();
    #endif

            WeightedCFLOBDDComplexFloatBoostMulNodeHandle gHandle(g);
            MatMultReturnT ret;

            if (c1.handleContents->level == 1)
            {
                auto tmp = gHandle.Reduce(reductionMapHandle, g_return_map.Size(), valList, false);
                ret = std::make_tuple(tmp.first, g_return_map, tmp.second * top_factor);
                // ret = std::make_tuple(gHandle, g_return_map, top_factor);
            }
            else
            {
                auto tmp = gHandle.Reduce(reductionMapHandle, g_return_map.Size(), valList, true);
                ret = std::make_tuple(tmp.first, g_return_map, tmp.second * top_factor);
                // ret = std::make_tuple(gHandle, g_return_map, top_factor);
            }
            matmult_hash.insert(std::make_pair(mmp, ret));
            return ret;
        }


        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCNOTNode(unsigned int level, unsigned int n, long int controller, long int controlled, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkCNOTNode(bn->numberOfVars, n, controller, controlled);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            }
            std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled);
            if (cnot_hashMap.find(p) != cnot_hashMap.end()){
                return cnot_hashMap[p];
            }
            WeightedCFLOBDDComplexFloatBoostInternalNode *g = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);

            if (level == 1)
            {
                if (controller == 0)
                {
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    auto b0 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(1.0,0.0));
                    auto b1 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(0.0,1.0));
                    g->BConnection[0] = Connection(b0, m01);
                    CFLOBDDReturnMapHandle m12;
                    m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
                    g->BConnection[1] = Connection(b1, m12);
                    g->numExits = 3;
                }
                else if (controlled == 0)
                {
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    auto b0 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(1.0,0.0));
                    auto b1 = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(0.0,1.0));
                    g->BConnection[0] = Connection(b1, m01);
                    CFLOBDDReturnMapHandle m10;
                    m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                    g->BConnection[1] = Connection(b0, m10);
                    g->numExits = 2;
                }
            }
            else
            {
                if (controller < n/2 && controlled < n/2 && controlled >= 0 && controller >= 0)
                {
                    // Case 1: Both in A Connection
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    auto aa = MkCNOTNode(level-1, n/2, controller, controlled);
                    g->AConnection = Connection(aa, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 2: Both in B Connection
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    auto bb = MkCNOTNode(level-1, n/2, controller - n/2, controlled - n/2);
                    g->BConnection[0] = Connection(bb, m01);
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller < n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 3: controller in A and controlled in B
                    CFLOBDDReturnMapHandle m012;
                    m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                    auto aa = MkCNOTNode(level-1, n/2, controller, -1);
                    g->AConnection = Connection(aa, m012);
                    g->numBConnections = 3;
                    g->BConnection = new Connection[3];
                    CFLOBDDReturnMapHandle m01, m1, m10;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                    auto bb = MkCNOTNode(level-1, n/2, -1, controlled - n/2);
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->BConnection[2] = Connection(bb, m10);
                    g->numExits = 2;
                }
                else if (controlled == -1 && controller < n/2 && controller >= 0)
                {
                    // Case 4: controller in A and controlled == -1
                    CFLOBDDReturnMapHandle m012;
                    m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                    auto aa = MkCNOTNode(level-1, n/2, controller, -1);
                    g->AConnection = Connection(aa, m012);
                    g->numBConnections = 3;
                    g->BConnection = new Connection[g->numBConnections];
                    CFLOBDDReturnMapHandle m01, m1, m21;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->BConnection[2] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m21);
                    g->numExits = 3;
                }
                else if (controlled == -1 && controller >= n/2)
                {
                    // Case 5: controller in B and controlled == -1
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    auto bb = MkCNOTNode(level-1, n/2, controller - n/2, -1);
                    CFLOBDDReturnMapHandle m012, m1;
                    m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[0] = Connection(bb, m012);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 3;
                }
                else if (controller == -1 && controlled >= 0 && controlled < n/2)
                {
                    // Case 6: controller == -1 and controlled in A
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    auto aa = MkCNOTNode(level-1, n/2, -1, controlled);
                    g->AConnection = Connection(aa, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    CFLOBDDReturnMapHandle m0, m10;
                    m0.AddToEnd(0); m0.Canonicalize();
                    m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m0);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m10);
                    g->numExits = 2;
                }
                else if (controller == -1 && controlled >= n/2)
                {
                    // Case 7: controller == -1 and controlled in B
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    auto bb = MkCNOTNode(level-1, n/2, -1, controlled - n/2);
                    CFLOBDDReturnMapHandle m0;
                    m0.AddToEnd(0); m0.Canonicalize();
                    g->BConnection[0] = Connection(bb, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m0);
                    g->numExits = 2;
                }
            }

    #ifdef PATH_COUNTING_ENABLED
            g->InstallPathCounts();
    #endif
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle gHandle = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(g);
            cnot_hashMap.insert(std::make_pair(p, gHandle));
            return gHandle;
        }

    
        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkSwapGateNode(bn->numberOfVars, controller, controlled, -1);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            }

            std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled) + ";" + std::to_string(case_num);
            if (swap_hashMap.find(p) != swap_hashMap.end()){
                return swap_hashMap[p];
            }
            WeightedCFLOBDDComplexFloatBoostInternalNode *g = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);
            if (level == 1)
            {
                // std::cout << "CaseBase" << std::endl;
                if (case_num == -1)
                {
                    assert(controlled == -1);
                    CFLOBDDReturnMapHandle m01, m23;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m23.AddToEnd(2); m23.AddToEnd(3); m23.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m23);
                    g->numExits = 4;
                }
                else if (case_num == 0)
                {
                    // [[1 0] [0 0]]
                    assert(controller == -1);
                    CFLOBDDReturnMapHandle m01, m1;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m1);
                    g->numExits = 2;
                }
                else if (case_num == 1)
                {
                    // [[0 1] [0 0]]
                    assert(controller == -1);
                    CFLOBDDReturnMapHandle m01, m0;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m0.AddToEnd(0); m0.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                    g->numExits = 2;	
                }
                else if (case_num == 2)
                {
                    // [[0 0] [1 0]]
                    assert(controller == -1);
                    CFLOBDDReturnMapHandle m01, m0, m10;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m0.AddToEnd(0); m0.Canonicalize();
                    m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m10);
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                    g->numExits = 2;	
                }
                else if (case_num == 3)
                {
                    // [[0 0] [0 1]]
                    assert(controller == -1);
                    CFLOBDDReturnMapHandle m01, m0;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m0.AddToEnd(0); m0.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                    g->numExits = 2;	
                }
            }
            else if (level == 2 && controller == 0 && controlled == 1)
            {
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkSwapGateNode(level - 1, controller, -1, case_num);
                CFLOBDDReturnMapHandle m0123;
                m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
                g->AConnection = Connection(atmp, m0123);
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle b0 = MkSwapGateNode(level-1, -1, controlled, 0);
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle b1 = MkSwapGateNode(level-1, -1, controlled, 1);
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle b2 = MkSwapGateNode(level-1, -1, controlled, 2);
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle b3 = MkSwapGateNode(level-1, -1, controlled, 3);
                CFLOBDDReturnMapHandle m01, m10;
                m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                g->numBConnections = 4;
                g->BConnection = new Connection[4];
                g->BConnection[0] = Connection(b0, m01);
                g->BConnection[1] = Connection(b2, m10);
                g->BConnection[2] = Connection(b1, m10);
                g->BConnection[3] = Connection(b3, m10);
                g->numExits = 2;
            }
            else if (level == 2 && controller == 0 && controlled == -1)
            {
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkSwapGateNode(level - 1, controller, -1, case_num);
                CFLOBDDReturnMapHandle m0123;
                m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
                g->AConnection = Connection(atmp, m0123);
                CFLOBDDReturnMapHandle m01, m21, m31, m41;
                m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
                m31.AddToEnd(3); m31.AddToEnd(1); m31.Canonicalize();
                m41.AddToEnd(4); m41.AddToEnd(1); m41.Canonicalize();
                g->numBConnections = 4;
                g->BConnection = new Connection[4];
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                g->BConnection[0] = Connection(Id, m01);
                g->BConnection[1] = Connection(Id, m21);
                g->BConnection[2] = Connection(Id, m31);
                g->BConnection[3] = Connection(Id, m41);
                g->numExits = 5;	
            }
            else if (level == 2 && controller == 1 && controlled == -1)
            {
                CFLOBDDReturnMapHandle m01;
                m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                g->AConnection = Connection(Id, m01);
                g->numBConnections = 2;
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkSwapGateNode(level-1, 0, -1, case_num);
                CFLOBDDReturnMapHandle m0123, m4;
                m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
                m4.AddToEnd(4); m4.Canonicalize();
                g->BConnection = new Connection[2];
                g->BConnection[0] = Connection(btmp, m0123);
                g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m4);
                g->numExits = 5;

            }
            else 
            {
                unsigned int n = pow(2, level - 1);
                if (controller < n/2 && controlled < n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 1: Both fall in A
                    // std::cout << "Case1" << std::endl;
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle aTmp = MkSwapGateNode(level-1, controller, controlled, case_num);
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0);
                    m01.AddToEnd(1);
                    m01.Canonicalize();
                    g->AConnection = Connection(aTmp, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->BConnection[0] = Connection(Id, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 2: Both fall in B region
                    // std::cout << "Case2" << std::endl;
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0);
                    m01.AddToEnd(1);
                    m01.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->AConnection = Connection(Id, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkSwapGateNode(level-1, controller - n/2, controlled - n/2, case_num);
                    g->BConnection[0] = Connection(btmp, m01);
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller < n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 3: controller in A and controlled in B
                    // std::cout << "Case3" << std::endl;
                    CFLOBDDReturnMapHandle m01234;
                    m01234.AddToEnd(0);
                    m01234.AddToEnd(1);
                    m01234.AddToEnd(2);
                    m01234.AddToEnd(3);
                    m01234.AddToEnd(4);
                    m01234.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkSwapGateNode(level-1, controller, -1, case_num);
                    g->AConnection = Connection(atmp, m01234);
                    g->numBConnections = 5;
                    g->BConnection = new Connection[5];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle b0 = MkSwapGateNode(level-1, -1, controlled - n/2, 0);
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->BConnection[0] = Connection(b0, m01);
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle b1 = MkSwapGateNode(level-1, -1, controlled - n/2, 2);
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle b2 = MkSwapGateNode(level-1, -1, controlled - n/2, 1);
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle b3 = MkSwapGateNode(level-1, -1, controlled - n/2, 3);
                    CFLOBDDReturnMapHandle m10, m1;
                    m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    if (controller == n/2 - 1)
                    {
                        // std::cout << b1 << std::endl;
                        // std::cout << b2 << std::endl;
                        // std::cout << b3 << std::endl;
                        g->BConnection[1] = Connection(b1, m10);
                        g->BConnection[2] = Connection(b2, m10);
                        g->BConnection[3] = Connection(b3, m10);
                        g->BConnection[4] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    }
                    else
                    {
                        g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                        g->BConnection[2] = Connection(b1, m10);
                        g->BConnection[3] = Connection(b2, m10);
                        g->BConnection[4] = Connection(b3, m10);
                    }
                    g->numExits = 2;
                }
                else if (controller < n/2 && controlled == -1 && controller >= 0)
                {
                    // Case 4: controller in A and controlled == -1
                    // std::cout << "Case4" << std::endl;
                    CFLOBDDReturnMapHandle m01234;
                    m01234.AddToEnd(0);
                    m01234.AddToEnd(1);
                    m01234.AddToEnd(2);
                    m01234.AddToEnd(3);
                    m01234.AddToEnd(4);
                    m01234.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkSwapGateNode(level-1, controller, -1, case_num);
                    g->AConnection = Connection(atmp, m01234);
                    g->numBConnections = 5;
                    g->BConnection = new Connection[5];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    CFLOBDDReturnMapHandle m01, m21, m31, m41, m1;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
                    m31.AddToEnd(3); m31.AddToEnd(1); m31.Canonicalize();
                    m41.AddToEnd(4); m41.AddToEnd(1); m41.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    if (controller == n/2 - 1)
                    {
                        g->BConnection[0] = Connection(Id, m01);
                        g->BConnection[1] = Connection(Id, m21);
                        g->BConnection[2] = Connection(Id, m31);
                        g->BConnection[3] = Connection(Id, m41);
                        g->BConnection[4] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    }
                    else
                    {
                        g->BConnection[0] = Connection(Id, m01);
                        g->BConnection[2] = Connection(Id, m21);
                        g->BConnection[3] = Connection(Id, m31);
                        g->BConnection[4] = Connection(Id, m41);
                        g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);	
                    }
                    g->numExits = 5;
                }
                else if (controller >= n/2 && controlled == -1 && controller >= 0)
                {
                    // Case 5: controller in B and controlled == -1
                    // std::cout << "Case5" << std::endl;
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->AConnection = Connection(Id, m01);
                    g->numBConnections = 2;
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkSwapGateNode(level-1, controller - n/2, -1, case_num);
                    CFLOBDDReturnMapHandle m0123, m1, m2, m4, m01234;
                    m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
                    m01234.AddToEnd(0); m01234.AddToEnd(1); m01234.AddToEnd(2); m01234.AddToEnd(3); m01234.AddToEnd(4); m01234.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    m2.AddToEnd(2); m2.Canonicalize();
                    m4.AddToEnd(4); m4.Canonicalize();
                    g->BConnection = new Connection[2];
                    // if (controller == n - 1)
                    // {
                    // 	g->BConnection[0] = Connection(btmp, m0123);
                    // 	g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
                    // }
                    // else
                    // {
                        g->BConnection[0] = Connection(btmp, m01234);
                        if (controller == n - 1)
                            g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m4);
                        else
                            g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);

                    // }
                    g->numExits = 5;
                } 
                else if (controller == -1 && controlled < n/2 && controlled >= 0)
                {
                    // Case 6: controller == -1 && controlled in A
                    // std::cout << "Case6" << std::endl;
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkSwapGateNode(level-1, -1, controlled, case_num);
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(atmp, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    CFLOBDDReturnMapHandle m1, m0;
                    m1.AddToEnd(1); m1.Canonicalize();
                    m0.AddToEnd(0); m0.Canonicalize();
                    if (case_num == 0)
                    {
                        g->BConnection[0] = Connection(Id, m01);
                        g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    }
                    else{
                        CFLOBDDReturnMapHandle m10;
                        m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                        g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m0);
                        g->BConnection[1] = Connection(Id, m10);
                    }
                    g->numExits = 2;
                }
                else if (controller == -1 && controlled >= n/2 && controlled >= 0)
                {
                    // Case 7: controller == -1 && controlled in B
                    // std::cout << "Case7" << std::endl;
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->AConnection = Connection(Id, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    CFLOBDDReturnMapHandle m1, m0;
                    m1.AddToEnd(1); m1.Canonicalize();
                    m0.AddToEnd(0); m0.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkSwapGateNode(level-1, -1, controlled - n/2, case_num);
                    if (case_num == 0)
                    {
                        g->BConnection[0] = Connection(btmp, m01);
                        g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    }
                    else
                    {
                        g->BConnection[0] = Connection(btmp, m01);
                        g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m0);
                    }
                    g->numExits = 2;
                }
                
            }
            #ifdef PATH_COUNTING_ENABLED
            g->InstallPathCounts();
    #endif
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle gHandle = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(g);
            swap_hashMap.insert(std::make_pair(p, gHandle));
            return gHandle;
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCPGateNode(std::unordered_map<std::string, WeightedCFLOBDDComplexFloatBoostMulNodeHandle>& cp_hashMap, unsigned int level, long int controller, long int controlled, BIG_COMPLEX_FLOAT theta_val, int cflobdd_kind)
        {

            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkCPGateNode(bn->numberOfVars, controller, controlled, theta_val);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            }

            std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled);
            if (cp_hashMap.find(p) != cp_hashMap.end()){
                return cp_hashMap[p];
            }
            WeightedCFLOBDDComplexFloatBoostInternalNode *g = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);
            if (level == 1)
            {
                if (controlled == -1)
                {
                    CFLOBDDReturnMapHandle m01, m12;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m12);
                    g->numExits = 3;
                }
                else if (controller == -1)
                {
                    CFLOBDDReturnMapHandle m01, m10;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(1, theta_val));
                    g->AConnection = Connection(atmp, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m10);
                    g->numExits = 2;
                }
            }
            else 
            {
                unsigned int n = pow(2, level - 1);
                if (controller < n/2 && controlled < n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 1: Both fall in A
                    // std::cout << "Case1" << std::endl;
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle aTmp = MkCPGateNode(cp_hashMap, level-1, controller, controlled, theta_val);
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0);
                    m01.AddToEnd(1);
                    m01.Canonicalize();
                    g->AConnection = Connection(aTmp, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->BConnection[0] = Connection(Id, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 2: Both fall in B region
                    // std::cout << "Case2" << std::endl;
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0);
                    m01.AddToEnd(1);
                    m01.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->AConnection = Connection(Id, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkCPGateNode(cp_hashMap, level-1, controller - n/2, controlled - n/2, theta_val);
                    g->BConnection[0] = Connection(btmp, m01);
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller < n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 3: controller in A and controlled in B
                    // std::cout << "Case3" << std::endl;
                    CFLOBDDReturnMapHandle m012;
                    m012.AddToEnd(0);
                    m012.AddToEnd(1);
                    m012.AddToEnd(2);
                    m012.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkCPGateNode(cp_hashMap, level-1, controller, -1, theta_val);
                    g->AConnection = Connection(atmp, m012);
                    g->numBConnections = 3;
                    g->BConnection = new Connection[g->numBConnections];
                    CFLOBDDReturnMapHandle m01, m1;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->BConnection[0] = Connection(Id, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkCPGateNode(cp_hashMap, level-1, -1, controlled - n/2, theta_val);
                    g->BConnection[2] = Connection(btmp, m01);
                    g->numExits = 2;
                }
                else if (controller < n/2 && controlled == -1 && controller >= 0)
                {
                    // Case 4: controller in A and controlled == -1
                    // std::cout << "Case4" << std::endl;
                    CFLOBDDReturnMapHandle m012;
                    m012.AddToEnd(0);
                    m012.AddToEnd(1);
                    m012.AddToEnd(2);
                    m012.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkCPGateNode(cp_hashMap, level-1, controller, -1, theta_val);
                    g->AConnection = Connection(atmp, m012);
                    g->numBConnections = 3;
                    g->BConnection = new Connection[g->numBConnections];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    CFLOBDDReturnMapHandle m01, m21, m1;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[0] = Connection(Id, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->BConnection[2] = Connection(Id, m21);
                    g->numExits = 3;
                }
                else if (controller >= n/2 && controlled == -1 && controller >= 0)
                {
                    // Case 5: controller in B and controlled == -1
                    // std::cout << "Case5" << std::endl;
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->AConnection = Connection(Id, m01);
                    g->numBConnections = 2;
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkCPGateNode(cp_hashMap, level-1, controller - n/2, -1, theta_val);
                    CFLOBDDReturnMapHandle m012, m1;
                    m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection = new Connection[g->numBConnections];
                    g->BConnection[0] = Connection(btmp, m012);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 3;
                } 
                else if (controller == -1 && controlled < n/2 && controlled >= 0)
                {
                    // Case 6: controller == -1 && controlled in A
                    // std::cout << "Case6" << std::endl;
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkCPGateNode(cp_hashMap, level-1, -1, controlled, theta_val);
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(atmp, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[0] = Connection(Id, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller == -1 && controlled >= n/2 && controlled >= 0)
                {
                    // Case 7: controller == -1 && controlled in B
                    // std::cout << "Case7" << std::endl;
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    g->AConnection = Connection(Id, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkCPGateNode(cp_hashMap, level-1, -1, controlled - n/2, theta_val);
                    g->BConnection[0] = Connection(btmp, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                
            }
            #ifdef PATH_COUNTING_ENABLED
            g->InstallPathCounts();
    #endif
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle gHandle = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(g);
            cp_hashMap.insert(std::make_pair(p, gHandle));
            return gHandle;
        }

        // i and j are always in the second half of the variables at top node
        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCSwapGateNode(unsigned int level, long int controller, long int i, long int j, int flag)
        {
            // std::cout << level << " " << controller << " " << i << " " << j << " " << flag << std::endl;
            // top level
            WeightedCFLOBDDComplexFloatBoostInternalNode *g = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);
            if (flag == 1)
            {
                assert(level >= 3);
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkCSwapGateNode(level-1, controller, i, j, 0);
                CFLOBDDReturnMapHandle m012;
                m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                g->AConnection = Connection(atmp, m012);
                g->numBConnections = 3;
                g->BConnection = new Connection[3];
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                CFLOBDDReturnMapHandle m01;
                m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                g->BConnection[0] = Connection(Id, m01);
                CFLOBDDReturnMapHandle m1;
                m1.AddToEnd(1); m1.Canonicalize();
                g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                unsigned int n = pow(2, level - 1);
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkSwapGateNode(level-1, i-n/2,j-n/2, -1);
                g->BConnection[2] = Connection(btmp, m01);	
                g->numExits = 2;

            }
            else
            {
                unsigned int n = pow(2, level - 1);
                if (level == 1)
                {
                    CFLOBDDReturnMapHandle m01, m12;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    g->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m12);
                    g->numExits = 3;
                }
                else if (controller < n/2)
                {
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle atmp = MkCSwapGateNode(level-1,controller,i,j,flag);
                    CFLOBDDReturnMapHandle m012;
                    m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                    g->AConnection = Connection(atmp, m012);
                    g->numBConnections = 3;
                    g->BConnection = new Connection[3];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    CFLOBDDReturnMapHandle m01,m1,m21;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
                    g->BConnection[0] = Connection(Id, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1],m1);
                    g->BConnection[2] = Connection(Id, m21);
                    g->numExits = 3;
                }
                else
                {
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle Id = WeightedCFLOBDDComplexFloatBoostMulNodeHandle::IdentityNode[level-1];
                    CFLOBDDReturnMapHandle m01,m1;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(Id, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    WeightedCFLOBDDComplexFloatBoostMulNodeHandle btmp = MkCSwapGateNode(level-1, controller - n/2, i,j,flag);
                    CFLOBDDReturnMapHandle m012;
                    m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[0] = Connection(btmp, m012);
                    g->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 3;
                }
                
            }
    #ifdef PATH_COUNTING_ENABLED
            g->InstallPathCounts();
    #endif
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle gHandle = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(g);
            return gHandle;
        }


        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkPauliYGateNode(unsigned int level, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkPauliYGateNode(bn->numberOfVars);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            } 
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkPauliZGateNode(unsigned int level, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkPauliZGateNode(bn->numberOfVars);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            } 
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkSGateNode(unsigned int level, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkSGateNode(bn->numberOfVars);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            } 
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkiSwapGateNode(unsigned int level, long int index1, long int index2, int case_num, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkiSwapGateNode(bn->numberOfVars, index1, index2, case_num);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            } 
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkPhaseShiftGateNode(unsigned int level, BIG_COMPLEX_FLOAT theta_val, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkPhaseShiftGateNode(bn->numberOfVars, theta_val);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            } 
        }


        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCCNOTNode(unsigned int level, long int controller1, long int controller2, long int controlled, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkCCNOTNode(bn->numberOfVars, controller1, controller2, controlled);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            } 
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkCSwapGate2Node(unsigned int level, long int controller, long int index1, long int index2, int case_num, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkCSwapGateNode(bn->numberOfVars, controller, index1, index2, -1);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn);
            } 
        }
        
	    std::pair<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, int> MkRestrictNode(unsigned int level, std::string s, int cflobdd_kind)
        {
            if (cflobdd_kind == 0)
            {
                WeightedBDDComplexFloatBoostTopNode *bn = new WeightedBDDComplexFloatBoostTopNode(level);
                bn->bddContents = WeightedMatrix1234BDDComplexFloatBoostMul::MkRestrictNode(bn->numberOfVars, s);
                int val = 1;
                if (s.find('1') == std::string::npos)
                    val = 0;
                return std::make_pair(WeightedCFLOBDDComplexFloatBoostMulNodeHandle(bn), val);
            } 
        }



    }

}