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
#include "wmatrix1234_fb_mul_node.h"
#include "weighted_cross_product.h"
#include "return_map_T.h"
#include "weighted_matmult_map.h"

namespace CFL_OBDD {

    namespace WeightedMatrix1234FloatBoostMul
    {
        std::vector<ReturnMapHandle<int>> commonly_used_return_maps;// m0, m1, m01, m10
        std::unordered_map<WeightedMatMultPair, MatMultReturnT, WeightedMatMultPair::MatMultPairHash> matmult_hash;
        std::unordered_map<std::string, WeightedCFLOBDDFloatBoostMulNodeHandle> cnot_hashMap;

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

        WeightedCFLOBDDFloatBoostMulNodeHandle MkIdRelationInterleavedNode(unsigned int level)
        {
            WeightedCFLOBDDFloatBoostInternalNode *n;

            if (level == 0) {
                return WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle;
            }
            return WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level];
    //         else if (level == 1) {
    //             n = new WeightedCFLOBDDFloatBoostInternalNode(level);

    //             WeightedCFLOBDDFloatBoostMulNodeHandle temp = WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle;
    //             n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
    //             n->numBConnections = 2;
    //             n->BConnection = new Connection[n->numBConnections];
    //             WeightedCFLOBDDFloatBoostMulNodeHandle b0 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostForkNode(1, 0)); 
    //             WeightedCFLOBDDFloatBoostMulNodeHandle b1 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostForkNode(0, 1));
    //             n->BConnection[0] = Connection(b0, commonly_used_return_maps[2]);//m01
    //             n->BConnection[1] = Connection(b1, commonly_used_return_maps[3]);//m10
    //         }
    //         else {  // Create an appropriate CFLOBDDInternalNode
    //             n = new WeightedCFLOBDDFloatBoostInternalNode(level);

    //             WeightedCFLOBDDFloatBoostMulNodeHandle temp = MkIdRelationInterleavedNode(level - 1);
    //             n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
    //             n->numBConnections = 2;
    //             n->BConnection = new Connection[n->numBConnections];
    //             n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
    //             n->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps[1]);//m1
    //         }
    //         n->numExits = 2;
    // #ifdef PATH_COUNTING_ENABLED
    //         n->InstallPathCounts();
    // #endif
    //         return WeightedCFLOBDDFloatBoostMulNodeHandle(n);
        } // MkIdRelationInterleavedNode
    

        WeightedCFLOBDDFloatBoostMulNodeHandle MkWalshInterleavedNode(unsigned int i)
        {
            assert(i >= 1);
            WeightedCFLOBDDFloatBoostInternalNode *n = new WeightedCFLOBDDFloatBoostInternalNode(i);
            if (i == 1) {  // Base case

                WeightedCFLOBDDFloatBoostMulNodeHandle temp = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostForkNode(1, 1));
                CFLOBDDReturnMapHandle m01;
                m01.AddToEnd(0);
                m01.AddToEnd(1);
                m01.Canonicalize();
                n->AConnection = Connection(temp, m01);
                n->numBConnections = 2;
                n->BConnection = new Connection[n->numBConnections];
                WeightedCFLOBDDFloatBoostMulNodeHandle b0 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostDontCareNode(1, 1));
                WeightedCFLOBDDFloatBoostMulNodeHandle b1 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostDontCareNode(1, -1));
                CFLOBDDReturnMapHandle m0;
                m0.AddToEnd(0);
                m0.Canonicalize();
                n->BConnection[0] = Connection(b0, m0);
                n->BConnection[1] = Connection(b1, m0);
            }
            else {
                WeightedCFLOBDDFloatBoostMulNodeHandle temp = MkWalshInterleavedNode(i - 1);
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
            return WeightedCFLOBDDFloatBoostMulNodeHandle(n);
        } // MkWalshInterleavedNode

        WeightedCFLOBDDFloatBoostMulNodeHandle KroneckerProduct2VocsNode(WeightedCFLOBDDFloatBoostMulNodeHandle m1, WeightedCFLOBDDFloatBoostMulNodeHandle m2, 
            int zero_index_m1, int zero_index_m2)
        {
            int level = m1.handleContents->level;
            WeightedCFLOBDDFloatBoostInternalNode* n = new WeightedCFLOBDDFloatBoostInternalNode(level + 1);
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
                n->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level], c0);
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
                n->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level], c);
                n->numExits = 2; 
            }
        #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
        #endif
            return WeightedCFLOBDDFloatBoostMulNodeHandle(n);
        }
    

        WeightedCFLOBDDFloatBoostMulNodeHandle MkExchangeInterleavedNode(unsigned int level)
        {
            WeightedCFLOBDDFloatBoostInternalNode *n;

            if (level == 0) {
                return WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01;
            }
            else if (level == 1) {
                n = new WeightedCFLOBDDFloatBoostInternalNode(level);

                n->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
                n->numBConnections = 2;
                n->BConnection = new Connection[n->numBConnections];
                n->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, commonly_used_return_maps[2]);//m01
                n->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, commonly_used_return_maps[3]);//m10
            }
            else {  // Create an appropriate CFLOBDDInternalNode
                n = new WeightedCFLOBDDFloatBoostInternalNode(level);

                WeightedCFLOBDDFloatBoostMulNodeHandle temp = MkIdRelationInterleavedNode(level - 1);
                n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
                n->numBConnections = 2;
                n->BConnection = new Connection[n->numBConnections];
                CFLOBDDReturnMapHandle m0, m01;
                m0.AddToEnd(0); m0.Canonicalize();
                m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                n->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], m0);//m01
                n->BConnection[1] = Connection(temp, m01);//m1
            }
            n->numExits = 2;
    #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
    #endif
            return WeightedCFLOBDDFloatBoostMulNodeHandle(n);
        } // MkIdRelationInterleavedNode
    


        std::tuple<WeightedCFLOBDDFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_FLOAT>
        Add(WeightedCFLOBDDFloatBoostMulNodeHandle b1, WeightedCFLOBDDFloatBoostMulNodeHandle b2,
            CFLOBDDMatMultMapHandle b1_m, CFLOBDDMatMultMapHandle b2_m, BIG_FLOAT b_f1, BIG_FLOAT b_f2)
        {   
            WeightedPairProductMapHandle<BIG_FLOAT> MapHandle;
            auto tmp =  PairProduct2(b1, b2, b_f1, b_f2, MapHandle);

            CFLOBDDMatMultMapHandle returnMapHandle;
        
            boost::unordered_map<size_t, unsigned int> reductionMap;
            ReductionMapHandle reductionMapHandle;
            unsigned int iterator = 0;
            WeightedValuesListHandle<BIG_FLOAT> valList;
            while (iterator < MapHandle.Size()){
                WeightedMatMultMapHandle<BIG_FLOAT> c1, c2;
                int first, second;
                first = MapHandle[iterator].first.First();
                second = MapHandle[iterator].first.Second();
                WeightedMatMultMapHandle<BIG_FLOAT> val;
                if (first == -1 && second == -1){
                    val = WeightedMatMultMapHandle<BIG_FLOAT>();
                    BIG_FLOAT v = 0.0;
                    val.Add(std::make_pair(-1,-1), v);
                    val.mapContents->contains_zero_val = true;
                }
                else {
                    c1 = b1_m.Lookup(first);
                    c2 = b2_m.Lookup(second);
                    BIG_FLOAT v1, v2;
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
            BIG_FLOAT factor = reduced_n.second;

            return std::make_tuple(reduced_n.first, returnMapHandle, factor);
        }

        MatMultReturnT
		MatrixMultiplyV4Node(WeightedCFLOBDDFloatBoostMulNodeHandle c1, WeightedCFLOBDDFloatBoostMulNodeHandle c2, int zero_exit_1, int zero_exit_2)
        {
            auto mmp = WeightedMatMultPair(c1, c2);
            if (matmult_hash.find(mmp) != matmult_hash.end())
                return matmult_hash[mmp];
            
            BIG_FLOAT zero = 0.0;
            if (c1 == WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1.handleContents->level] || c2 == WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c2.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                WeightedMatMultMapHandle<BIG_FLOAT> tmp;
                tmp.Add(std::make_pair(-1,-1), zero);
                tmp.mapContents->contains_zero_val = true;
                tmp.Canonicalize();
                m.AddToEnd(tmp);
                m.Canonicalize();
                return std::make_tuple(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1.handleContents->level], m, zero);
            }

            if (c1 == WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[c1.handleContents->level] && c2 == WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[c2.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                BIG_FLOAT one = 1.0;
                WeightedMatMultMapHandle<BIG_FLOAT> m0, m1;
                m0.Add(std::make_pair(0,0), one); m0.Canonicalize();
                m1.Add(std::make_pair(-1,-1), zero); m1.mapContents->contains_zero_val = true; m1.Canonicalize();
                m.AddToEnd(m0); m.AddToEnd(m1); m.Canonicalize();
                return std::make_tuple(c1, m, one);
            }

            if (c1 == WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[c1.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                BIG_FLOAT one = 1.0;
                for (int i = 0; i < c2.handleContents->numExits; i++)
                {
                    WeightedMatMultMapHandle<BIG_FLOAT> tmp;
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
            if (c2 == WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[c1.handleContents->level])
            {
                CFLOBDDMatMultMapHandle m;
                BIG_FLOAT one = 1.0;
                for (int i = 0; i < c1.handleContents->numExits; i++)
                {
                    WeightedMatMultMapHandle<BIG_FLOAT> tmp;
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
            WeightedValuesListHandle<BIG_FLOAT> valList;
            WeightedCFLOBDDFloatBoostInternalNode* g = new WeightedCFLOBDDFloatBoostInternalNode(c1.handleContents->level);
            ReductionMapHandle reductionMapHandle;
            BIG_FLOAT top_factor = 1.0;
            
            if (c1.handleContents->level == 1){
                WeightedCFLOBDDFloatBoostInternalNode* c1_internal = (WeightedCFLOBDDFloatBoostInternalNode *)c1.handleContents;
                WeightedCFLOBDDFloatBoostInternalNode* c2_internal = (WeightedCFLOBDDFloatBoostInternalNode *)c2.handleContents;

                /*
                    [[a_0, b_0],    [[a_1, b_1],
                     [c_0, d_0]]     [c_1, d_1]]
                */
                BIG_FLOAT a0, b0, c0, d0, a1, b1, c1, d1;
                WeightedCFLOBDDFloatBoostLeafNode* M1_A = (WeightedCFLOBDDFloatBoostLeafNode *)c1_internal->AConnection.entryPointHandle->handleContents;
                WeightedCFLOBDDFloatBoostLeafNode* M2_A = (WeightedCFLOBDDFloatBoostLeafNode *)c2_internal->AConnection.entryPointHandle->handleContents;
                a0 = b0 = M1_A->lweight;
                c0 = d0 = M1_A->rweight;
                a1 = b1 = M2_A->lweight;
                c1 = d1 = M2_A->rweight;
                WeightedMatMultMapHandle<BIG_FLOAT> v1, v2, v3, v4;
                int M1_numB = c1_internal->numBConnections;
                int M2_numB = c2_internal->numBConnections;
                WeightedCFLOBDDFloatBoostLeafNode* M1_b0 = (WeightedCFLOBDDFloatBoostLeafNode *)c1_internal->BConnection[0].entryPointHandle->handleContents;
                WeightedCFLOBDDFloatBoostLeafNode* M1_b1 = (WeightedCFLOBDDFloatBoostLeafNode *)c1_internal->BConnection[M1_numB-1].entryPointHandle->handleContents;
                a0 = a0 * M1_b0->lweight;
                b0 = b0 * M1_b0->rweight;
                c0 = c0 * M1_b1->lweight;
                d0 = d0 * M1_b1->rweight;
                WeightedCFLOBDDFloatBoostLeafNode* M2_b0 = (WeightedCFLOBDDFloatBoostLeafNode *)c2_internal->BConnection[0].entryPointHandle->handleContents;
                WeightedCFLOBDDFloatBoostLeafNode* M2_b1 = (WeightedCFLOBDDFloatBoostLeafNode *)c2_internal->BConnection[M2_numB-1].entryPointHandle->handleContents;
                a1 = a1 * M2_b0->lweight;
                b1 = b1 * M2_b0->rweight;
                c1 = c1 * M2_b1->lweight;
                d1 = d1 * M2_b1->rweight;

                int M1_b0_numE = M1_b0->numExits;
                int M1_b1_numE = M1_b1->numExits;
                int M2_b0_numE = M2_b0->numExits;
                int M2_b1_numE = M2_b1->numExits;

                // v1
                BIG_FLOAT a0a1 = a0 * a1;
                BIG_FLOAT b0c1 = b0 * c1;
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
                BIG_FLOAT a0b1 = a0 * b1;
                BIG_FLOAT b0d1 = b0 * d1;
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
                BIG_FLOAT c0a1 = c0 * a1;
                BIG_FLOAT d0c1 = d0 * c1;
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
                BIG_FLOAT c0b1 = c0 * b1;
                BIG_FLOAT d0d1 = d0 * d1;
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
                        //     g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                        // }
                        // else {
                            g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0);
                            g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0); 
                        // }
                        reductionMapHandle.AddToEnd(0);
                        valList.AddToEnd(v1.mapContents->contains_zero_val ? 0 : 1);
                        g_return_map.AddToEnd(v1);
                    }
                    else{
                        g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0);
                        g->numBConnections = 1;
                        g->BConnection = new Connection[g->numBConnections];
                        CFLOBDDReturnMapHandle m01;
                        m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                        // if (v1.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                        // else if (v2.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                        // else
                            g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01); 
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
                    g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    int la = 1, ra = 1;
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    if (v1 == v2)
                    {
                        CFLOBDDReturnMapHandle m0;
                        m0.AddToEnd(0); m0.Canonicalize();
                        // if (v1.mapContents->contains_zero_val == true){
                        //     la = 0;
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                        // }
                        // else{
                        //     la = 1;
                            g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m0);
                        // }
                        g_return_map.AddToEnd(v1);
                        reductionMapHandle.AddToEnd(0);
                        valList.AddToEnd(v1.mapContents->contains_zero_val ? 0 : 1);
                    }
                    else{
                        la = 1;
                        // if (v1.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                        // else if (v2.mapContents->contains_zero_val == true)
                        //     g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                        // else
                            g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
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
                        //     g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m);
                        // }
                        // else {
                        //     ra = 1;
                            g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m);
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
                        //     g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m1);
                        // else if (v4.mapContents->contains_zero_val == true)
                        //     g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m1);
                        // else
                            g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m1);
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
                    //     g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    // else if (la == 1 && ra == 0)
                    //     g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle10, m01);
                    // else if (la == 0 && ra == 1)
                    //     g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle01, m01);
                    // assert(!(la == 0 && ra == 0));
                }

            }
            else{
                
                WeightedCFLOBDDFloatBoostInternalNode* c1_internal = (WeightedCFLOBDDFloatBoostInternalNode *)c1.handleContents;
                WeightedCFLOBDDFloatBoostInternalNode* c2_internal = (WeightedCFLOBDDFloatBoostInternalNode *)c2.handleContents;
                int level = c1.handleContents->level;
                // populate zero exit information
                int a_zero_exit_1 = -1, a_zero_exit_2 = -1;
                if (zero_exit_1 != -1){
                    for (int i = 0; i < c1_internal->numBConnections; i++)
                    {
                        if (*(c1_internal->BConnection[i].entryPointHandle) == WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1])
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
                        if (*(c2_internal->BConnection[i].entryPointHandle) == WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1])
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
                    WeightedMatMultMapHandle<BIG_FLOAT> matmult_returnmap = std::get<1>(aa)[i];
                    WeightedCFLOBDDFloatBoostMulNodeHandle ans;
                    CFLOBDDMatMultMapHandle ans_matmult_map;
                    BIG_FLOAT ans_factor = 1.0;
                    bool first = true;
                    // Consider Multiplication of M1 and M2
                    for (auto &v : matmult_returnmap.mapContents->map){
                        unsigned int M1_index = v.first.first;
                        unsigned int M2_index = v.first.second;
                        std::tuple<WeightedCFLOBDDFloatBoostMulNodeHandle, CFLOBDDMatMultMapHandle, BIG_FLOAT> bb_old;
                        CFLOBDDMatMultMapHandle new_bb_return;
                        WeightedCFLOBDDFloatBoostMulNodeHandle n1, n2;
                        CFLOBDDReturnMapHandle n1_return, n2_return;
                        int b_zero_exit_1 = -1, b_zero_exit_2 = -1;
                        if (M1_index == -1 || v.second == 0)
                        {
                            n1 = WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1_internal->level-1];
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
                            n2 = WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[c1_internal->level-1];
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
                            WeightedMatMultMapHandle<BIG_FLOAT> tmp;
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

            WeightedCFLOBDDFloatBoostMulNodeHandle gHandle(g);
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


        WeightedCFLOBDDFloatBoostMulNodeHandle MkCNOTNode(unsigned int level, unsigned int n, long int controller, long int controlled)
        {
            std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled);
            if (cnot_hashMap.find(p) != cnot_hashMap.end()){
                return cnot_hashMap[p];
            }
            WeightedCFLOBDDFloatBoostInternalNode *g = new WeightedCFLOBDDFloatBoostInternalNode(level);

            if (level == 1)
            {
                if (controller == 0)
                {
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    auto b0 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostForkNode(1.0,0.0));
                    auto b1 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostForkNode(0.0,1.0));
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
                    g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::CFLOBDDForkNodeHandle, m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[g->numBConnections];
                    auto b0 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostForkNode(1.0,0.0));
                    auto b1 = WeightedCFLOBDDFloatBoostMulNodeHandle(new WeightedCFLOBDDFloatBoostForkNode(0.0,1.0));
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
                    g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->numExits = 2;
                }
                else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
                {
                    // Case 2: Both in B Connection
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    auto bb = MkCNOTNode(level-1, n/2, controller - n/2, controlled - n/2);
                    g->BConnection[0] = Connection(bb, m01);
                    CFLOBDDReturnMapHandle m1;
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
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
                    g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
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
                    g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
                    g->BConnection[2] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m21);
                    g->numExits = 3;
                }
                else if (controlled == -1 && controller >= n/2)
                {
                    // Case 5: controller in B and controlled == -1
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    auto bb = MkCNOTNode(level-1, n/2, controller - n/2, -1);
                    CFLOBDDReturnMapHandle m012, m1;
                    m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
                    m1.AddToEnd(1); m1.Canonicalize();
                    g->BConnection[0] = Connection(bb, m012);
                    g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m1);
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
                    g->BConnection[0] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m0);
                    g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m10);
                    g->numExits = 2;
                }
                else if (controller == -1 && controlled >= n/2)
                {
                    // Case 7: controller == -1 and controlled in B
                    CFLOBDDReturnMapHandle m01;
                    m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
                    g->AConnection = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::IdentityNode[level-1], m01);
                    g->numBConnections = 2;
                    g->BConnection = new Connection[2];
                    auto bb = MkCNOTNode(level-1, n/2, -1, controlled - n/2);
                    CFLOBDDReturnMapHandle m0;
                    m0.AddToEnd(0); m0.Canonicalize();
                    g->BConnection[0] = Connection(bb, m01);
                    g->BConnection[1] = Connection(WeightedCFLOBDDFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1], m0);
                    g->numExits = 2;
                }
            }

    #ifdef PATH_COUNTING_ENABLED
            g->InstallPathCounts();
    #endif
            WeightedCFLOBDDFloatBoostMulNodeHandle gHandle = WeightedCFLOBDDFloatBoostMulNodeHandle(g);
            cnot_hashMap.insert(std::make_pair(p, gHandle));
            return gHandle;
        }

    }

}