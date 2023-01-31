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
#include "wmatrix1234_complex_fb_mul_bdd_node.h"
#include "weighted_cross_product_bdd.h"
#include "return_map_T.h"
#include "weighted_matmult_map.h"
#include "weighted_values.h"
#include "pair_t.h"
#include "bool_op.h"

namespace CFL_OBDD {

    namespace WeightedMatrix1234BDDComplexFloatBoostMul
    {

        std::unordered_map<RenamePair, WeightedBDDComplexFloatBoostMulNodeHandle, RenamePair::RenamePairHash> rename_map;
        std::unordered_map<WeightedBDDMatMultPair, BDDMatMultReturnT, WeightedBDDMatMultPair::MatMultPairHash> matmult_hash;

        WeightedBDDComplexFloatBoostMulNodeHandle MkIdRelationInterleavedNode(unsigned int numVars, unsigned int index)
        {
            if (numVars == 0)
            {
                return WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode; 
            }
            if (numVars == 2)
            {

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else
            {
                auto tmp = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkWalshInterleavedNode(unsigned int numVars, unsigned int index)
        {
            assert(numVars >= 2);
            if (numVars == 2)
            {

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y2->rweight = -1;
                y2->lweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);

            }
            else
            {
                auto tmp = MkWalshInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = tmp;
                y2->rweight = -1;
                y2->lweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = tmp;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            } 
        }

        WeightedBDDComplexFloatBoostMulNodeHandle RenameNodes(WeightedBDDComplexFloatBoostMulNodeHandle m, long int offset)
        {
            RenamePair p (m, offset);
            if (rename_map.find(p) != rename_map.end())
                return rename_map[p];
            if (m.handleContents->NodeKind() == LEAF)
                    return m;
            WeightedBDDComplexFloatBoostMulInternalNode* m2_i = (WeightedBDDComplexFloatBoostMulInternalNode *)m.handleContents;
            WeightedBDDComplexFloatBoostMulInternalNode* m2_new = new WeightedBDDComplexFloatBoostMulInternalNode(m2_i->GetIndex() + offset);
            m2_new->lweight = m2_i->lweight;
            m2_new->rweight = m2_i->rweight;
            m2_new->leftNode = RenameNodes(m2_i->leftNode, offset);
            m2_new->rightNode = RenameNodes(m2_i->rightNode, offset);
            WeightedBDDComplexFloatBoostMulNodeHandle t(m2_new);
            rename_map[p] = t;
            return t;
        }

        WeightedBDDComplexFloatBoostMulNodeHandle KroneckerProduct2VocsNode(std::unordered_map<WeightedBDDComplexFloatBoostMulNodeHandle, WeightedBDDComplexFloatBoostMulNodeHandle, WeightedBDDComplexFloatBoostMulNodeHandle::WeightedBDDNodeHandle_Hash>& hash_map, WeightedBDDComplexFloatBoostMulNodeHandle m1, WeightedBDDComplexFloatBoostMulNodeHandle m2, long int numVars)
        {
            if (hash_map.find(m1) != hash_map.end())
                return hash_map[m1];
            if (m1.handleContents->NodeKind() == LEAF){
                if (m1 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode)
                    return m1;
                else
                    return m2;
            }
            WeightedBDDComplexFloatBoostMulInternalNode* m1_i = (WeightedBDDComplexFloatBoostMulInternalNode *)m1.handleContents;
            WeightedBDDComplexFloatBoostMulInternalNode* m1_new = new WeightedBDDComplexFloatBoostMulInternalNode(m1_i->GetIndex());
            auto lNode = m1_i->leftNode;
            auto rNode = m1_i->rightNode;
            WeightedBDDComplexFloatBoostMulNodeHandle lNode_new;
            if (lNode.handleContents->NodeKind() == LEAF)
            {
                WeightedBDDComplexFloatBoostMulLeafNode* lNode_leaf = (WeightedBDDComplexFloatBoostMulLeafNode *)lNode.handleContents;
                if (lNode_leaf->value == getAnnhilatorValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>())
                    lNode_new = lNode;
                else {
                    lNode_new = RenameNodes(m2, numVars);
                }
            } else {
                lNode_new = KroneckerProduct2VocsNode(hash_map, lNode, m2, numVars);
            }

            WeightedBDDComplexFloatBoostMulNodeHandle rNode_new;
            if (rNode.handleContents->NodeKind() == LEAF)
            {
                WeightedBDDComplexFloatBoostMulLeafNode* rNode_leaf = (WeightedBDDComplexFloatBoostMulLeafNode *)rNode.handleContents;
                if (rNode_leaf->value == getAnnhilatorValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>())
                    rNode_new = rNode;
                else {
                    rNode_new = RenameNodes(m2, numVars);
                }
            } else {
                rNode_new = KroneckerProduct2VocsNode(hash_map, rNode, m2, numVars);
            }

            m1_new->leftNode = lNode_new;
            m1_new->rightNode = rNode_new;
            m1_new->lweight = m1_i->lweight;
            m1_new->rweight = m1_i->rweight;
            if (m1_new->leftNode == m1_new->rightNode && m1_new->lweight == m1_new->rweight)
                return m1_new->leftNode;
            WeightedBDDComplexFloatBoostMulNodeHandle m1h_new(m1_new);
            hash_map[m1] = m1h_new;
            return m1h_new;
        }

        std::pair<std::vector<WeightedBDDComplexFloatBoostMulNodeHandle>, std::vector<BIG_COMPLEX_FLOAT>> 
            FillEntries(WeightedBDDComplexFloatBoostMulNodeHandle m) {
                WeightedBDDComplexFloatBoostMulNodeHandle m11, m12, m13, m14;
                BIG_COMPLEX_FLOAT lf1, lf2, lf3, lf4;
                if (m.handleContents->GetIndex() % 2 == 0) {
                    WeightedBDDComplexFloatBoostMulInternalNode* m1_x = (WeightedBDDComplexFloatBoostMulInternalNode *)m.handleContents;
                    if (m1_x->leftNode.handleContents->NodeKind() == INTERNAL)
                    {
                        if (m1_x->leftNode.handleContents->GetIndex() == m1_x->GetIndex() + 1)
                        {
                            WeightedBDDComplexFloatBoostMulInternalNode* y = (WeightedBDDComplexFloatBoostMulInternalNode *)m1_x->leftNode.handleContents;
                            m11 = y->leftNode;
                            m12 = y->rightNode;
                            lf1 = m1_x->lweight * y->lweight;
                            lf2 = m1_x->lweight * y->rweight;
                        }
                        else
                        {
                            m11 = m1_x->leftNode;
                            m12 = m1_x->leftNode;
                            lf1 = m1_x->lweight;
                            lf2 = m1_x->lweight;
                        }
                    }
                    else 
                    {
                        WeightedBDDComplexFloatBoostMulLeafNode* y = (WeightedBDDComplexFloatBoostMulLeafNode *)m1_x->leftNode.handleContents; 
                        m11 = m12 = m1_x->leftNode;
                        lf1 = lf2 = m1_x->lweight * y->value;
                    }
                    if (m1_x->rightNode.handleContents->NodeKind() == INTERNAL)
                    {
                        if (m1_x->rightNode.handleContents->GetIndex() == m1_x->GetIndex() + 1)
                        {
                            WeightedBDDComplexFloatBoostMulInternalNode* y = (WeightedBDDComplexFloatBoostMulInternalNode *)m1_x->rightNode.handleContents;
                            m13 = y->leftNode;
                            m14 = y->rightNode;
                            lf3 = m1_x->rweight * y->lweight;
                            lf4 = m1_x->rweight * y->rweight;
                        }
                        else
                        {
                            m13 = m1_x->rightNode;
                            m14 = m1_x->rightNode;
                            lf3 = m1_x->rweight;
                            lf4 = m1_x->rweight;
                        }
                    }
                    else 
                    {
                        WeightedBDDComplexFloatBoostMulLeafNode* y = (WeightedBDDComplexFloatBoostMulLeafNode *)m1_x->rightNode.handleContents; 
                        m13 = m14 = m1_x->rightNode;
                        lf3 = lf4 = m1_x->rweight * y->value;
                    }
                }
                else {
                    WeightedBDDComplexFloatBoostMulInternalNode* m1_y = (WeightedBDDComplexFloatBoostMulInternalNode *)m.handleContents;
                    m11 = m13 = m1_y->leftNode;
                    m12 = m14 = m1_y->rightNode;
                    lf1 = lf3 = m1_y->lweight;
                    lf2 = lf4 = m1_y->rweight; 
                }

                std::vector<WeightedBDDComplexFloatBoostMulNodeHandle> ret_1 = {m11, m12, m13, m14};
                std::vector<BIG_COMPLEX_FLOAT> ret_2 = {lf1, lf2, lf3, lf4};

            return std::make_pair(ret_1, ret_2);
        }

        BDDMatMultReturnT
            Add(WeightedBDDComplexFloatBoostMulNodeHandle m1, WeightedBDDComplexFloatBoostMulNodeHandle m2, BIG_COMPLEX_FLOAT f1, BIG_COMPLEX_FLOAT f2, unsigned int numVars)
        {
            if (m1 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode && m2 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode)
            {
                std::vector<BIG_COMPLEX_FLOAT> v_ret;
                return std::make_tuple(m2, v_ret, 0);
            }
            if (m1 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode) {
                // Incorrect. Need to fill this. But won't affect algorithm.
                std::vector<BIG_COMPLEX_FLOAT> v_ret;
                return std::make_tuple(m2, v_ret, f2); 
            }
            else if (m2 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode) {
                std::vector<BIG_COMPLEX_FLOAT> v_ret;
                return std::make_tuple(m1, v_ret, f1);  
            }
                WeightedBDDPairProductMapHandle<BIG_COMPLEX_FLOAT> pairProductMapHandle;
                auto ans = BDDPairProduct2(m1, 
                            m2, 
                            numVars, 
                            f1,
                            f2,
                            pairProductMapHandle,
                            PlusFunc);
                std::vector<BIG_COMPLEX_FLOAT> v_ret = pairProductMapHandle.mapContents->mapArray;
                return std::make_tuple(ans, v_ret, pairProductMapHandle.mapContents->factor);
        } 
        
        BDDMatMultReturnT
              MatrixMultiplyV4Node(WeightedBDDComplexFloatBoostMulNodeHandle m1, WeightedBDDComplexFloatBoostMulNodeHandle m2, long int numVars, long int count)
        {
            WeightedBDDMatMultPair mp(m1, m2);
            if (matmult_hash.find(mp) != matmult_hash.end())
                return matmult_hash[mp];
            WeightedBDDComplexFloatBoostMulNodeHandle m11, m12, m13, m14;
            WeightedBDDComplexFloatBoostMulNodeHandle m21, m22, m23, m24;
            BIG_COMPLEX_FLOAT top_factor = 1.0;
            BIG_COMPLEX_FLOAT lf1, lf2, lf3, lf4;
            BIG_COMPLEX_FLOAT rf1, rf2, rf3, rf4;

            if (m1 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode || m2 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode)
            {
                std::vector<BIG_COMPLEX_FLOAT> tmp = {0.0};
                return std::make_tuple(WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode, 
                            tmp,
                            0.0); 
            }

            if (m1.handleContents->NodeKind() == LEAF && m2.handleContents->NodeKind() == LEAF)
            {
                if (m1 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode || m2 == WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode)
                {
                    std::vector<BIG_COMPLEX_FLOAT> tmp = {0.0};
                    return std::make_tuple(WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode, 
                                tmp,
                                0.0);
                }
                else {
                    std::vector<BIG_COMPLEX_FLOAT> tmp = {1.0};
                    return std::make_tuple(WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode, 
                                tmp,
                                1.0); 
                }
            }

            int prev_count = count;

            if (m1.handleContents->GetIndex() == m2.handleContents->GetIndex() && m1.handleContents->GetIndex() % 2 == 0)
            {
                // simplest case
                if (m1.handleContents->GetIndex() > count && m1.handleContents->GetIndex() % 2 == 0){
                    top_factor = 2 * (m1.handleContents->GetIndex() - count);
                    count = m1.handleContents->GetIndex() + 2;
                }
                else {
                    count += 2;
                }
                auto t1 = FillEntries(m1);
                m11 = t1.first[0]; m12 = t1.first[1]; m13 = t1.first[2]; m14 = t1.first[3];
                lf1 = t1.second[0]; lf2 = t1.second[1]; lf3 = t1.second[2]; lf4 = t1.second[3];

                auto t2 = FillEntries(m2);
                m21 = t2.first[0]; m22 = t2.first[1]; m23 = t2.first[2]; m24 = t2.first[3];
                rf1 = t2.second[0]; rf2 = t2.second[1]; rf3 = t2.second[2]; rf4 = t2.second[3];
            }

            else if (m1.handleContents->GetIndex() > m2.handleContents->GetIndex() && m2.handleContents->NodeKind() == LEAF)
            {
                if (m1.handleContents->GetIndex() > count && m1.handleContents->GetIndex() % 2 == 0){
                    top_factor = 2 * (m1.handleContents->GetIndex() - count);
                    count = m1.handleContents->GetIndex() + 2;
                }
                else {
                    count += 2;
                }

                auto t1 = FillEntries(m1);
                m11 = t1.first[0]; m12 = t1.first[1]; m13 = t1.first[2]; m14 = t1.first[3];
                lf1 = t1.second[0]; lf2 = t1.second[1]; lf3 = t1.second[2]; lf4 = t1.second[3]; 
                WeightedBDDComplexFloatBoostMulLeafNode* m2_x = (WeightedBDDComplexFloatBoostMulLeafNode *)m2.handleContents;
                m21 = m22 = m23 = m24 = m2;
                rf1 = rf2 = rf3 = rf4 = m2_x->value;
            }

            else if (m1.handleContents->GetIndex() < m2.handleContents->GetIndex() && m1.handleContents->NodeKind() == LEAF)
            {
               if (m2.handleContents->GetIndex() > count && m2.handleContents->GetIndex() % 2 == 0){
                    top_factor = 2 * (m2.handleContents->GetIndex() - count);
                    count = m2.handleContents->GetIndex() + 2;
                }
                else {
                    count += 2;
                }

                auto t1 = FillEntries(m2);
                m21 = t1.first[0]; m22 = t1.first[1]; m23 = t1.first[2]; m24 = t1.first[3];
                lf1 = t1.second[0]; lf2 = t1.second[1]; lf3 = t1.second[2]; lf4 = t1.second[3]; 
                WeightedBDDComplexFloatBoostMulLeafNode* m1_x = (WeightedBDDComplexFloatBoostMulLeafNode *)m1.handleContents;
                m11 = m12 = m13 = m14 = m1;
                lf1 = lf2 = lf3 = lf4 = m1_x->value; 
            }
        
            else if (m1.handleContents->GetIndex() > m2.handleContents->GetIndex())
            {
                if (m1.handleContents->GetIndex() > count && m1.handleContents->GetIndex() % 2 == 0){
                    top_factor = 2 * (m1.handleContents->GetIndex() - count);
                    count = m1.handleContents->GetIndex() + 2;
                }
                else {
                    count += 2;
                }
                auto t1 = FillEntries(m1);
                m11 = t1.first[0]; m12 = t1.first[1]; m13 = t1.first[2]; m14 = t1.first[3];
                lf1 = t1.second[0]; lf2 = t1.second[1]; lf3 = t1.second[2]; lf4 = t1.second[3];
                if (m2.handleContents->GetIndex() + 1 == m1.handleContents->GetIndex()){
                    auto t2 = FillEntries(m2);
                    m21 = t2.first[0]; m22 = t2.first[1]; m23 = t2.first[2]; m24 = t2.first[3];
                    rf1 = t2.second[0]; rf2 = t2.second[1]; rf3 = t2.second[2]; rf4 = t2.second[3]; 
                }
                else {
                    m21 = m22 = m23 = m24 = m2;
                    rf1 = rf2 = rf3 = rf4 = 1;
                }
            }

            else if (m1.handleContents->GetIndex() < m2.handleContents->GetIndex())
            {
                if (m2.handleContents->GetIndex() > count && m2.handleContents->GetIndex() % 2 == 0){
                    top_factor = 2 * (m2.handleContents->GetIndex() - count);
                    count = m2.handleContents->GetIndex() + 2;
                }
                else {
                    count += 2;
                }

                auto t1 = FillEntries(m2);
                m21 = t1.first[0]; m22 = t1.first[1]; m23 = t1.first[2]; m24 = t1.first[3];
                rf1 = t1.second[0]; rf2 = t1.second[1]; rf3 = t1.second[2]; rf4 = t1.second[3];
                if (m1.handleContents->GetIndex() + 1 == m2.handleContents->GetIndex()){
                    auto t2 = FillEntries(m1);
                    m11 = t2.first[0]; m12 = t2.first[1]; m13 = t2.first[2]; m14 = t2.first[3];
                    lf1 = t2.second[0]; lf2 = t2.second[1]; lf3 = t2.second[2]; lf4 = t2.second[3]; 
                }
                else {
                    m11 = m12 = m13 = m14 = m1;
                    lf1 = lf2 = lf3 = lf4 = 1;
                }
            }
            // std::cout << "count: " << prev_count << std::endl;
            // [[p1 p2] [p3 p4]] = [[m11 m12] [m13 m14]] * [[m21 m22] [m23 m24]]
            auto p1_1 = MatrixMultiplyV4Node(m11, m21, numVars, count);
            auto p1_2 = MatrixMultiplyV4Node(m12, m23, numVars, count);

            auto p1_t = Add(std::get<0>(p1_1), std::get<0>(p1_2), lf1 * rf1 * std::get<2>(p1_1), lf2 * rf3 * std::get<2>(p1_2), numVars);
            auto p1 = std::get<0>(p1_t);
            // std::cout << "p1" << std::endl;
            // p1.print(std::cout);
            auto f_p1 = std::get<2>(p1_t);
            // std::cout << "f_p1 " << f_p1 << std::endl;
            // std::cout << std::endl;
            std::vector<BIG_COMPLEX_FLOAT> v_p1 = std::get<1>(p1_t);

            auto p2_1 = MatrixMultiplyV4Node(m11, m22, numVars, count);
            auto p2_2 = MatrixMultiplyV4Node(m12, m24, numVars, count);

            auto p2_t = Add(std::get<0>(p2_1), std::get<0>(p2_2), lf1 * rf2 * std::get<2>(p2_1), lf2 * rf4 * std::get<2>(p2_2), numVars);
            auto p2 = std::get<0>(p2_t);
            // std::cout << "p2" << std::endl;
            // p2.print(std::cout);
            auto f_p2 = std::get<2>(p2_t);
            // std::cout << "f_p2 " << f_p2 << std::endl;
            // std::cout << std::endl;
            std::vector<BIG_COMPLEX_FLOAT> v_p2 = std::get<1>(p2_t);

            auto p3_1 = MatrixMultiplyV4Node(m13, m21, numVars, count);
            auto p3_2 = MatrixMultiplyV4Node(m14, m23, numVars, count);

            auto p3_t = Add(std::get<0>(p3_1), std::get<0>(p3_2), lf3 * rf1 * std::get<2>(p3_1), lf4 * rf3 * std::get<2>(p3_2), numVars);
            auto p3 = std::get<0>(p3_t);
            // std::cout << "p3" << std::endl;
            // p3.print(std::cout);
            auto f_p3 = std::get<2>(p3_t);
            // std::cout << "f_p3 " << f_p3 << std::endl;
            // std::cout << std::endl;
            std::vector<BIG_COMPLEX_FLOAT> v_p3 = std::get<1>(p3_t);

            auto p4_1 = MatrixMultiplyV4Node(m13, m22, numVars, count);
            auto p4_2 = MatrixMultiplyV4Node(m14, m24, numVars, count);

            auto p4_t = Add(std::get<0>(p4_1), std::get<0>(p4_2), lf3 * rf2 * std::get<2>(p4_1), lf4 * rf4 * std::get<2>(p4_2), numVars);
            auto p4 = std::get<0>(p4_t);
            // std::cout << "p4" << std::endl;
            // p4.print(std::cout);
            auto f_p4 = std::get<2>(p4_t);
            // std::cout << "f_p4 " << f_p4 << std::endl;
            // std::cout << std::endl;
            std::vector<BIG_COMPLEX_FLOAT> v_p4 = std::get<1>(p4_t);

            WeightedBDDComplexFloatBoostMulNodeHandle yh1;
            WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(prev_count + 1);
            y1->leftNode = p1;
            y1->rightNode = p2;
            auto w1 = computeInverseValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(f_p1, f_p2);
            y1->lweight = std::get<1>(w1);
            y1->rweight = std::get<2>(w1);

            if (y1->leftNode == y1->rightNode && y1->lweight == y1->rweight)
            {
                yh1 = y1->leftNode;
            }
            else {
                yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);
            }


            WeightedBDDComplexFloatBoostMulNodeHandle yh2;
            WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(prev_count + 1);
            y2->leftNode = p3;
            y2->rightNode = p4;
            auto w2 = computeInverseValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(f_p3, f_p4);
            y2->lweight = std::get<1>(w2);
            y2->rweight = std::get<2>(w2);

            if (y2->leftNode == y2->rightNode && y2->lweight == y2->rweight)
            {
                yh2 = y2->leftNode;
            }
            else {
                yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);
            }


            WeightedBDDComplexFloatBoostMulNodeHandle xh;
            WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(prev_count);

            x->leftNode = yh1;
            x->rightNode = yh2;
            auto w = computeInverseValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(std::get<0>(w1), std::get<0>(w2));
            x->lweight = std::get<1>(w);
            x->rweight = std::get<2>(w);

            if (x->leftNode == x->rightNode && x->lweight == x->rweight)
            {
                xh = x->leftNode;
            }
            else {
                xh = WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }

            std::vector<BIG_COMPLEX_FLOAT> ret_v;
            for (auto i : v_p1)
                ret_v.push_back(i);

            for (auto i : v_p2){
                int j = 0;
                for (j = 0; j < ret_v.size(); j++){
                    if (ret_v[j] == i)
                        break;
                }
                if (j == ret_v.size())
                    ret_v.push_back(i);
            }

            for (auto i : v_p3){
                int j = 0;
                for (j = 0; j < ret_v.size(); j++){
                    if (ret_v[j] == i)
                        break;
                }
                if (j == ret_v.size())
                    ret_v.push_back(i);
            }

            for (auto i : v_p4){
                int j = 0;
                for (j = 0; j < ret_v.size(); j++){
                    if (ret_v[j] == i)
                        break;
                }
                if (j == ret_v.size())
                    ret_v.push_back(i);
            }
            // std::cout << "xh" << std::endl;
            // xh.print(std::cout);
            auto ret = std::make_tuple(xh, ret_v, std::get<0>(w));
            matmult_hash.insert(std::make_pair(mp, ret));
            return ret;
        }
    
    
        WeightedBDDComplexFloatBoostMulNodeHandle MkCNOTNode(unsigned int numVars, unsigned int n, long int controller, long int controlled, unsigned int index)
        {
            // assuming controller < controlled
            if (numVars == 4 && controller == 0 && controlled == 1)
            {

                WeightedBDDComplexFloatBoostMulNodeHandle z1 = MkIdRelationInterleavedNode(2, index + 2);
                WeightedBDDComplexFloatBoostMulNodeHandle z2 = MkNegationMatrixInterleavedNode(2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = z1;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = z2;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else if (index/2 == controlled && numVars == 2)
            {
                return MkNegationMatrixInterleavedNode(numVars, index);
            }
            else if (index/2 == controller) {
                auto tmp = MkCNOTNode(numVars - 2, n - 1, -1, controlled, index + 2);
                auto I = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = MkIdRelationInterleavedNode(numVars - 2, index + 2);
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);  
            }
            else if (index/2 == controlled)
            {
                auto tmp = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->rightNode = tmp;
                y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->rweight = 1;
                y1->lweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->leftNode = tmp;
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->lweight = 1;
                y2->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
            else 
            {
                auto tmp = MkCNOTNode(numVars - 2, n - 1, controller, controlled, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkNegationMatrixInterleavedNode(unsigned int numVars, unsigned int index)
        {
            assert(numVars >= 2);
            if (numVars == 2)
            {

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh2;
                x->rightNode = yh1;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else
            {
                auto tmp = MkNegationMatrixInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh2;
                x->rightNode = yh1;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            } 
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkPauliYGateNode(unsigned int numVars, unsigned int index)
        {
            assert(numVars >= 2);
            if (numVars == 2)
            {

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = -1;
                x->leftNode = yh2;
                x->rightNode = yh1;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else
            {
                auto tmp = MkPauliYGateNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = -1;
                x->leftNode = yh2;
                x->rightNode = yh1;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }  
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkPauliZGateNode(unsigned int numVars, unsigned int index)
        {
            assert(numVars >= 2);
            if (numVars == 2)
            {

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = -1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else
            {
                auto tmp = MkPauliZGateNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = -1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }  
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkSGateNode(unsigned int numVars, unsigned int index)
        {
            assert(numVars >= 2);
            if (numVars == 2)
            {

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = BIG_COMPLEX_FLOAT(0, 1);
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else
            {
                auto tmp = MkPauliYGateNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = BIG_COMPLEX_FLOAT(0, 1);
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }  
        }
    
    
        WeightedBDDComplexFloatBoostMulNodeHandle MkSwapGateNode(unsigned int numVars, long int index1, long int index2, int case_num, unsigned int index)
        {
            if (numVars == 4 && index1 == 0 && index2 == 1)
            {
                WeightedBDDComplexFloatBoostMulInternalNode* m10 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 3);
                m10->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                m10->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m10->lweight = 1;
                m10->rweight = 0;
                auto m10_h = WeightedBDDComplexFloatBoostMulNodeHandle(m10);

                WeightedBDDComplexFloatBoostMulInternalNode* m01 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 3);
                m01->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m01->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                m01->lweight = 0;
                m01->rweight = 1;
                auto m01_h = WeightedBDDComplexFloatBoostMulNodeHandle(m01);

                WeightedBDDComplexFloatBoostMulInternalNode* m0010 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m0010->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m0010->rightNode = m10_h;
                m0010->lweight = 0;
                m0010->rweight = 1;
                auto m0010_h = WeightedBDDComplexFloatBoostMulNodeHandle(m0010); 
                
                WeightedBDDComplexFloatBoostMulInternalNode* m0100 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m0100->leftNode = m01_h;
                m0100->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m0100->lweight = 1;
                m0100->rweight = 0;
                auto m0100_h = WeightedBDDComplexFloatBoostMulNodeHandle(m0100);

                WeightedBDDComplexFloatBoostMulInternalNode* m0001 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m0001->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m0001->rightNode = m01_h;
                m0001->lweight = 0;
                m0001->rweight = 1;
                auto m0001_h = WeightedBDDComplexFloatBoostMulNodeHandle(m0001);
                
                WeightedBDDComplexFloatBoostMulInternalNode* m1000 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m1000->leftNode = m10_h;
                m1000->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m1000->lweight = 1;
                m1000->rweight = 0;
                auto m1000_h = WeightedBDDComplexFloatBoostMulNodeHandle(m1000);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = m1000;
                y1->rightNode = m0010;
                y1->lweight = 1;
                y1->rweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = m0001;
                y2->leftNode = m0100;
                y2->rweight = 1;
                y2->lweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else if (index/2 == index1) {
                auto tmp1 = MkSwapGateNode(numVars - 2, index1, index2, 0, index + 2);
                auto tmp2 = MkSwapGateNode(numVars - 2, index1, index2, 1, index + 2);
                auto tmp3 = MkSwapGateNode(numVars - 2, index1, index2, 2, index + 2);
                auto tmp4 = MkSwapGateNode(numVars - 2, index1, index2, 3, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp1;
                y1->rightNode = tmp2;
                y1->lweight = 1;
                y1->rweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp4;
                y2->leftNode = tmp3;
                y2->rweight = 1;
                y2->lweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);  
            }
            else if (index/2 == index2)
            {
                auto tmp = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                if (case_num == 0)
                {
                    // 1000
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->leftNode = tmp;
                    y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->rweight = 0;
                    y1->lweight = 1;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 1;
                    x->rweight = 0;
                    x->leftNode = yh1;
                    x->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

                if (case_num == 1)
                {
                    // 0010
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->leftNode = tmp;
                    y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->lweight = 1;
                    y1->rweight = 0;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 0;
                    x->rweight = 1;
                    x->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    x->rightNode = yh1;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

                if (case_num == 2)
                {
                    // 0100
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->rightNode = tmp;
                    y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->rweight = 1;
                    y1->lweight = 0;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 1;
                    x->rweight = 0;
                    x->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    x->leftNode = yh1;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

                if (case_num == 3)
                {
                    // 0001
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->rightNode = tmp;
                    y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->lweight = 0;
                    y1->rweight = 1;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 0;
                    x->rweight = 1;
                    x->rightNode = yh1;
                    x->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

            }
            else 
            {
                auto tmp = MkSwapGateNode(numVars - 2, index1, index2, case_num, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            } 
        }


        WeightedBDDComplexFloatBoostMulNodeHandle MkiSwapGateNode(unsigned int numVars, long int index1, long int index2, int case_num, unsigned int index)
        {
            if (numVars == 4 && index1 == 0 && index2 == 1)
            {
                WeightedBDDComplexFloatBoostMulInternalNode* m10 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 3);
                m10->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                m10->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m10->lweight = 1;
                m10->rweight = 0;
                auto m10_h = WeightedBDDComplexFloatBoostMulNodeHandle(m10);

                WeightedBDDComplexFloatBoostMulInternalNode* m01 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 3);
                m01->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m01->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                m01->lweight = 0;
                m01->rweight = 1;
                auto m01_h = WeightedBDDComplexFloatBoostMulNodeHandle(m01);

                WeightedBDDComplexFloatBoostMulInternalNode* m0010 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m0010->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m0010->rightNode = m10_h;
                m0010->lweight = 0;
                m0010->rweight = 1;
                auto m0010_h = WeightedBDDComplexFloatBoostMulNodeHandle(m0010); 
                
                WeightedBDDComplexFloatBoostMulInternalNode* m0100 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m0100->leftNode = m01_h;
                m0100->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m0100->lweight = 1;
                m0100->rweight = 0;
                auto m0100_h = WeightedBDDComplexFloatBoostMulNodeHandle(m0100);

                WeightedBDDComplexFloatBoostMulInternalNode* m0001 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m0001->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m0001->rightNode = m01_h;
                m0001->lweight = 0;
                m0001->rweight = 1;
                auto m0001_h = WeightedBDDComplexFloatBoostMulNodeHandle(m0001);
                
                WeightedBDDComplexFloatBoostMulInternalNode* m1000 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                m1000->leftNode = m10_h;
                m1000->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                m1000->lweight = 1;
                m1000->rweight = 0;
                auto m1000_h = WeightedBDDComplexFloatBoostMulNodeHandle(m1000);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = m1000;
                y1->rightNode = m0010;
                y1->lweight = 1;
                y1->rweight = BIG_COMPLEX_FLOAT(0, -1);
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = m0001;
                y2->leftNode = m0100;
                y2->rweight = BIG_COMPLEX_FLOAT(0, 1);
                y2->lweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = BIG_COMPLEX_FLOAT(0, 1);
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else if (index/2 == index1) {
                auto tmp1 = MkSwapGateNode(numVars - 2, index1, index2, 0, index + 2);
                auto tmp2 = MkSwapGateNode(numVars - 2, index1, index2, 1, index + 2);
                auto tmp3 = MkSwapGateNode(numVars - 2, index1, index2, 2, index + 2);
                auto tmp4 = MkSwapGateNode(numVars - 2, index1, index2, 3, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp1;
                y1->rightNode = tmp2;
                y1->lweight = 1;
                y1->rweight = BIG_COMPLEX_FLOAT(0, 1);
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp4;
                y2->leftNode = tmp3;
                y2->rweight = BIG_COMPLEX_FLOAT(0, -1);
                y2->lweight = 1;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = BIG_COMPLEX_FLOAT(0, 1);
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);  
            }
            else if (index/2 == index2)
            {
                auto tmp = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                if (case_num == 0)
                {
                    // 1000
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->leftNode = tmp;
                    y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->rweight = 0;
                    y1->lweight = 1;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 1;
                    x->rweight = 0;
                    x->leftNode = yh1;
                    x->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

                if (case_num == 1)
                {
                    // 0010
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->leftNode = tmp;
                    y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->lweight = 1;
                    y1->rweight = 0;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 0;
                    x->rweight = 1;
                    x->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    x->rightNode = yh1;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

                if (case_num == 2)
                {
                    // 0100
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->rightNode = tmp;
                    y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->rweight = 1;
                    y1->lweight = 0;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 1;
                    x->rweight = 0;
                    x->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    x->leftNode = yh1;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

                if (case_num == 3)
                {
                    // 0001
                    WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                    y1->rightNode = tmp;
                    y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    y1->lweight = 0;
                    y1->rweight = 1;
                    WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                    WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                    x->lweight = 0;
                    x->rweight = 1;
                    x->rightNode = yh1;
                    x->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;

                    return WeightedBDDComplexFloatBoostMulNodeHandle(x);    
                }

            }
            else 
            {
                auto tmp = MkSwapGateNode(numVars - 2, index1, index2, case_num, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            } 
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkCPGateNode(unsigned int numVars, long int controller, long int controlled, BIG_COMPLEX_FLOAT theta_val, unsigned int index)
        {
            // assuming controller < controlled
            if (numVars == 4 && controller == 0 && controlled == 1)
            {

                WeightedBDDComplexFloatBoostMulNodeHandle z1 = MkIdRelationInterleavedNode(2, index + 2);
                WeightedBDDComplexFloatBoostMulInternalNode* t1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 3);
                t1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                t1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                t1->lweight = 1;
                t1->rweight = 0;
                auto t1_h = WeightedBDDComplexFloatBoostMulNodeHandle(t1);
                WeightedBDDComplexFloatBoostMulInternalNode* t2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 3);
                t2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                t2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                t2->lweight = 0;
                t2->rweight = 1;
                auto t2_h = WeightedBDDComplexFloatBoostMulNodeHandle(t2);

                WeightedBDDComplexFloatBoostMulInternalNode* z2_i = new WeightedBDDComplexFloatBoostMulInternalNode(index + 2);
                z2_i->leftNode = t1_h;
                z2_i->rightNode = t2_h;
                z2_i->lweight = 1;
                z2_i->rweight = theta_val;
                WeightedBDDComplexFloatBoostMulNodeHandle z2 = WeightedBDDComplexFloatBoostMulNodeHandle(z2_i);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = z1;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = z2;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 

            }
            else if (index/2 == controller) {
                auto tmp = MkCPGateNode(numVars - 2, -1, controlled, theta_val, index + 2);
                auto I = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = MkIdRelationInterleavedNode(numVars - 2, index + 2);
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);  
            }
            else if (index/2 == controlled)
            {
                auto tmp = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = theta_val;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
            else 
            {
                auto tmp = MkCPGateNode(numVars - 2, controller, controlled, theta_val, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkPhaseShiftGateNode(unsigned int numVars, BIG_COMPLEX_FLOAT theta_val, unsigned int index)
        {
            assert (numVars == 2);


            WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
            y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
            y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
            y1->lweight = 1;
            y1->rweight = 0;

            WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

            WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
            y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
            y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
            y2->rweight = 1;
            y2->lweight = 0;
            WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

            WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
            x->lweight = 1;
            x->rweight = theta_val;
            x->leftNode = yh1;
            x->rightNode = yh2;

            return WeightedBDDComplexFloatBoostMulNodeHandle(x);

        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkCCNOTNode(unsigned int numVars, long int controller1, long int controller2, long int controlled, unsigned int index)
        {
            // assuming controller < controlled
            if (index/2 == controller1 || index/2 == controller2) {
                auto tmp = MkCCNOTNode(numVars - 2, controller1, controller2, controlled, index + 2);
                auto I = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = I;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);  
            }
            else if (index/2 == controlled)
            {
                auto tmp = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->rightNode = tmp;
                y1->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->rweight = 1;
                y1->lweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->leftNode = tmp;
                y2->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->lweight = 1;
                y2->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
            else 
            {
                auto tmp = MkCCNOTNode(numVars - 2, controller1, controller2, controlled, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkCSwapGateNode(unsigned int numVars, long int controller, long int index1, long int index2, int case_num, unsigned int index)
        {
            
            if (index/2 == controller)
            {
                auto tmp = MkSwapGateNode(numVars - 2, index1, index2, case_num, index + 2);
                auto I = MkIdRelationInterleavedNode(numVars - 2, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = I;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2; 

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);   
            }
            else 
            {
                auto tmp = MkCSwapGateNode(numVars - 2, controller, index1, index2, case_num, index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* y2 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y2->rightNode = tmp;
                y2->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y2->rweight = 1;
                y2->lweight = 0;
                WeightedBDDComplexFloatBoostMulNodeHandle yh2 = WeightedBDDComplexFloatBoostMulNodeHandle(y2);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh2;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            } 
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkRestrictNode(unsigned int numVars, std::string s, unsigned int index)
        {
            if (numVars == 0)
                return WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
            if (s[0] == 'X')
            {
                auto tmp = MkRestrictNode(numVars - 2, s.substr(1), index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 1;
                x->leftNode = yh1;
                x->rightNode = yh1;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
            else if (s[0] == '0')
            {
                auto tmp = MkRestrictNode(numVars - 2, s.substr(1), index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 1;
                x->rweight = 0;
                x->leftNode = yh1;
                x->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x);
            }
            else if (s[0] == '1')
            {
                auto tmp = MkRestrictNode(numVars - 2, s.substr(1), index + 2);

                WeightedBDDComplexFloatBoostMulInternalNode* y1 = new WeightedBDDComplexFloatBoostMulInternalNode(index + 1);
                y1->leftNode = tmp;
                y1->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                y1->lweight = 1;
                y1->rweight = 0;

                WeightedBDDComplexFloatBoostMulNodeHandle yh1 = WeightedBDDComplexFloatBoostMulNodeHandle(y1);

                WeightedBDDComplexFloatBoostMulInternalNode* x = new WeightedBDDComplexFloatBoostMulInternalNode(index);
                x->lweight = 0;
                x->rweight = 1;
                x->rightNode = yh1;
                x->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;

                return WeightedBDDComplexFloatBoostMulNodeHandle(x); 
            }
        }
    
    }
}