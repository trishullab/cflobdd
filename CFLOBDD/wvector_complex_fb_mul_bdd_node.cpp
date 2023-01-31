#include <cassert>
#include <cstdio>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include "wvector_complex_fb_mul_bdd_node.h"
#include "weighted_bdd_node_t.h"
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"

namespace CFL_OBDD {

    namespace WeightedVectorBDDComplexFloatBoostMul {

        WeightedBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int numVars, unsigned int index, long int id)
        {
            if (numVars == 1)
            {
                assert(index < 2);
                if (index == 0){
                    WeightedBDDComplexFloatBoostInternalNode* g = new WeightedBDDComplexFloatBoostInternalNode(id);
                    g->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                    g->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    g->lweight = 1;
                    g->rweight = 0;
                    return WeightedBDDComplexFloatBoostMulNodeHandle(g);
                }
                else {
                    WeightedBDDComplexFloatBoostInternalNode* g = new WeightedBDDComplexFloatBoostInternalNode(id);
                    g->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                    g->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    g->rweight = 1;
                    g->lweight = 0;
                    return WeightedBDDComplexFloatBoostMulNodeHandle(g);
                }
            }

            WeightedBDDComplexFloatBoostInternalNode *n = new WeightedBDDComplexFloatBoostInternalNode(id);

            unsigned int mid = index & (1UL << (numVars-1));
            if (mid == 0)
            {
                auto tmp = MkBasisVectorNode(numVars-1, index, id + 1);
                n->leftNode = tmp;
                n->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                n->lweight = 1;
                n->rweight = 0;
            }
            else 
            {
                auto tmp = MkBasisVectorNode(numVars-1, index - mid, id + 1);
                n->rightNode = tmp;
                n->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                n->lweight = 0;
                n->rweight = 1;
            }
            
            return WeightedBDDComplexFloatBoostMulNodeHandle(n);
        }

        WeightedBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int numVars, std::string s, long int id)
        {
            if (numVars == 1)
            {
                assert(s.length() == 1);
                if (s[0] =='0'){
                    WeightedBDDComplexFloatBoostInternalNode* g = new WeightedBDDComplexFloatBoostInternalNode(id);
                    g->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                    g->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    g->lweight = 1;
                    g->rweight = 0;
                    return WeightedBDDComplexFloatBoostMulNodeHandle(g);
                }
                else {
                    WeightedBDDComplexFloatBoostInternalNode* g = new WeightedBDDComplexFloatBoostInternalNode(id);
                    g->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::IdentityLeafNode;
                    g->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                    g->rweight = 1;
                    g->lweight = 0;
                    return WeightedBDDComplexFloatBoostMulNodeHandle(g);
                }
            }

            WeightedBDDComplexFloatBoostInternalNode *n = new WeightedBDDComplexFloatBoostInternalNode(id);

            if (s[0] == '0')
            {
                auto tmp = MkBasisVectorNode(numVars-1, s.substr(1), id + 1);
                n->leftNode = tmp;
                n->rightNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                n->lweight = 1;
                n->rweight = 0;
            }
            else 
            {
                auto tmp = MkBasisVectorNode(numVars-1, s.substr(1), id + 1);
                n->rightNode = tmp;
                n->leftNode = WeightedBDDComplexFloatBoostMulNodeHandle::AnnhilatorLeafNode;
                n->rweight = 1;
                n->lweight = 0;
            }
            
            return WeightedBDDComplexFloatBoostMulNodeHandle(n);
        }


        int chooseIndexRandomly(std::vector<std::pair<long double, unsigned int>>& weights, double random_value)
        {
            long double val = 0;
            for (int i = 0; i < weights.size(); i++)
            {
                if (val + weights[i].first >= random_value)
                    return weights[i].second;
                val += weights[i].first;
                
            }
            return weights[weights.size()-1].second;
        }

        std::pair<std::string,std::string> SamplingNode(WeightedBDDComplexFloatBoostMulNodeHandle nh, unsigned int numVars, std::mt19937 mt,            std::uniform_real_distribution<double> dis, long int count)
        {
            if (nh.handleContents->NodeKind() == LEAF){
                std::string add_string1 = "", add_string2 = "";
                while (numVars > 0)
                {
                    // auto r = rand()/RAND_MAX;
                    auto r = mt() % 2;
                    if (r == 0){
                        if (count % 2 == 0)
                            add_string1 += "0";
                        else
                            add_string2 += "0";
                    }
                    else{
                        if (count % 2 == 0)
                            add_string1 += "1";
                        else
                            add_string2 += "1";
                    }
                    numVars--;
                }
                return std::make_pair(add_string1, add_string2);
            }
            WeightedBDDComplexFloatBoostInternalNode * nNode = (WeightedBDDComplexFloatBoostInternalNode *)nh.handleContents;
            std::string add_string1 = "", add_string2 = "";
            while (count < nNode->GetIndex())
            {
                auto r = mt() % 2;
                if (r == 0){
                    if (count % 2 == 0)
                        add_string1 += "0";
                    else
                        add_string2 += "0";
                }
                else{
                    if (count % 2 == 0)
                        add_string1 += "1";
                    else
                        add_string2 += "1";
                }
                count++;
                numVars--;
            }
            long double leftProb = nNode->leftWeightOfPathsAsAmpsToExit / nNode->weightOfPathsAsAmpsToExit;
            long double rightProb = nNode->rightWeightOfPathsAsAmpsToExit / nNode->weightOfPathsAsAmpsToExit;

            std::vector<std::pair<long double, unsigned int>> weights = {std::make_pair(leftProb, 0), std::make_pair(rightProb, 1)};
            // double random_value = ((double)rand() / (RAND_MAX));
            double random_value = dis(mt);
            int chosen_index = -1;
            chosen_index = chooseIndexRandomly(weights, random_value);
            if (chosen_index == 0)
            {
                auto s = SamplingNode(nNode->leftNode, numVars - 1, mt, dis, count + 1);
                if (nNode->GetIndex() % 2 == 0)
                {
                    return std::make_pair(add_string1 + "0" + s.first, add_string2 +  s.second);
                }
                else
                {
                    return std::make_pair(add_string1 + s.first, add_string2 + "0" + s.second);
                }
            }
            else 
            {
                auto s = SamplingNode(nNode->rightNode, numVars - 1, mt, dis, count + 1);
                if (nNode->GetIndex() % 2 == 0)
                {
                    return std::make_pair(add_string1 + "1" + s.first, add_string2 + s.second);
                }
                else
                {
                    return std::make_pair(add_string1 + s.first, add_string2 + "1" + s.second);
                }
            }
        }
        

    //#ifdef PATH_COUNTING_ENABLED
        
    }
//#endif
}  // namespace CFL_OBDD
