#include <cassert>
#include <cstdio>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include "wvector_complex_fb_mul_node.h"
#include "weighted_cflobdd_node_t.h"
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"

namespace CFL_OBDD {

    namespace WeightedVectorComplexFloatBoostMul {

        std::vector<ReturnMapHandle<int>> commonly_used_return_maps_vector;// m0, m1, m01, m10

        void InitReturnMapVectorHandles(){
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
            commonly_used_return_maps_vector.push_back(m0);
            commonly_used_return_maps_vector.push_back(m1);
            commonly_used_return_maps_vector.push_back(m01);
            commonly_used_return_maps_vector.push_back(m10);
        }

        void VectorInitializerNode()
        {
            // Empty for now
            InitReturnMapVectorHandles();
            return;
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, unsigned int index, int cflobdd_kind)
        {
            if (cflobdd_kind == 0){
                WeightedBDDComplexFloatBoostTopNode *n = new WeightedBDDComplexFloatBoostTopNode(level); 
                n->bddContents = WeightedVectorBDDComplexFloatBoostMul::MkBasisVectorNode(level, index);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
            }
            if (level == 0)
            {
                assert(index < 2);
                if (index == 0)
                    return WeightedCFLOBDDComplexFloatBoostMulNodeHandle( new WeightedCFLOBDDComplexFloatBoostForkNode(1, 0));
                else
                    return WeightedCFLOBDDComplexFloatBoostMulNodeHandle( new WeightedCFLOBDDComplexFloatBoostForkNode(0, 1));
            }

            WeightedCFLOBDDComplexFloatBoostInternalNode *n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);

            unsigned int higherOrderIndex = index >> (1 << (level - 1));
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle tempANodeHandle = MkBasisVectorNode(level-1, higherOrderIndex);
            n->AConnection = Connection(tempANodeHandle, commonly_used_return_maps_vector[2]);//m01

            n->numBConnections = 2;
            n->BConnection = new Connection[n->numBConnections];
            unsigned int lowerOrderIndex = index & ((1 << (1 << (level - 1))) - 1);
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle tempBNodeHandle = MkBasisVectorNode(level - 1, lowerOrderIndex);
            if (higherOrderIndex == 0)
            {
                if (lowerOrderIndex == 0)
                    n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps_vector[1]);//m1
                else
                    n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps_vector[0]);//m0
                n->BConnection[0] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
            }
            else{
                n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps_vector[0]);//m0
                if (lowerOrderIndex == 0)
                    n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[3]);//m10
                else
                    n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
            }
            n->numExits = 2;
    #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
    #endif
            return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkBasisVectorNode(unsigned int level, std::string s, int cflobdd_kind)
        {

            if (cflobdd_kind == 0){
                WeightedBDDComplexFloatBoostTopNode *n = new WeightedBDDComplexFloatBoostTopNode(level); 
                n->bddContents = WeightedVectorBDDComplexFloatBoostMul::MkBasisVectorNode(level, s);
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
            }

            if (level == 0)
            {
                assert(s.length() == 1);
                if (s[0] == '0')
                    return WeightedCFLOBDDComplexFloatBoostMulNodeHandle( new WeightedCFLOBDDComplexFloatBoostForkNode(1, 0));
                else
                    return WeightedCFLOBDDComplexFloatBoostMulNodeHandle( new WeightedCFLOBDDComplexFloatBoostForkNode(0, 1));
            }

            WeightedCFLOBDDComplexFloatBoostInternalNode *n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);
            std::string first_half_s = s.substr(0, s.length() / 2);
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle tempANodeHandle = MkBasisVectorNode(level - 1, first_half_s);
            n->AConnection = Connection(tempANodeHandle, commonly_used_return_maps_vector[2]);//m01

            
            n->numBConnections = 2;
            n->BConnection = new Connection[n->numBConnections];
            std::string half_string = s.substr(s.length() / 2);
            WeightedCFLOBDDComplexFloatBoostMulNodeHandle tempBNodeHandle = MkBasisVectorNode(level - 1, half_string);
            if (first_half_s.find('1') == std::string::npos)
            {
                if (half_string.find('1') == std::string::npos)
                    n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps_vector[1]);//m1
                else
                    n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps_vector[0]);//m0
                n->BConnection[0] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
            }
            else{
                n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level - 1], commonly_used_return_maps_vector[0]);//m0
                if (half_string.find('1') == std::string::npos)
                    n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[3]);//m10
                else
                    n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
            }
            n->numExits = 2;
    #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
    #endif
            return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle VectorToMatrixInterleavedNode(std::unordered_map<WeightedCFLOBDDComplexFloatBoostMulNodeHandle, WeightedCFLOBDDComplexFloatBoostMulNodeHandle, WeightedCFLOBDDComplexFloatBoostMulNodeHandle::WeightedCFLOBDDNodeHandleT_Hash> hashMap, 
                                WeightedCFLOBDDComplexFloatBoostMulNodeHandle& nh)
        {
            // if (hashMap)
            // {
            //     WeightedCFLOBDDComplexFloatBoostMulNodeHandle lookupResult;
            //     hashMap->Fetch(nh, lookupResult);
            //     return lookupResult;
            // }
            if (hashMap.find(nh) != hashMap.end())
                return hashMap[nh];
            WeightedCFLOBDDComplexFloatBoostInternalNode *nhNode = (WeightedCFLOBDDComplexFloatBoostInternalNode *)nh.handleContents;
            
            if (nh == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode[nh.handleContents->level])
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode[nh.handleContents->level + 1];
            if (nh == WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[nh.handleContents->level])
                return WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[nh.handleContents->level + 1];

            WeightedCFLOBDDComplexFloatBoostInternalNode *n = new WeightedCFLOBDDComplexFloatBoostInternalNode(nh.handleContents->level + 1);

            if (nh.handleContents->level == 0)
            {
                CFLOBDDReturnMapHandle mI;
                for (int i = 0; i < nhNode->numExits; i++)
                    mI.AddToEnd(i);
                mI.Canonicalize();
                n->AConnection = Connection(nh, mI);
                n->numBConnections = mI.Size();
                n->BConnection = new Connection[n->numBConnections];
                WeightedCFLOBDDComplexFloatBoostLeafNode* nhL = (WeightedCFLOBDDComplexFloatBoostLeafNode *)nhNode;
                if (n->numBConnections == 1)
                {
                    if ((nhL->lweight != 0) && (nhL->rweight != 0))
                    {
                        CFLOBDDReturnMapHandle m0;
                        m0.AddToEnd(0); m0.AddToEnd(1); m0.Canonicalize();
                        n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode[0], m0);
                        n->numExits = 1;
                    }
                }
                else if (n->numBConnections == 2)
                {
                    if ((nhL->lweight != 0) && (nhL->rweight != 0))
                    {
                        CFLOBDDReturnMapHandle m0, m1;
                        m0.AddToEnd(0); m0.Canonicalize();
                        m1.AddToEnd(1); m1.Canonicalize();
                        n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode[0], m0);
                        n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode[0], m1);
                        n->numExits = 2; 
                    }
                    else if ((nhL->lweight == 0))
                    {
                        CFLOBDDReturnMapHandle m0, m1;
                        m0.AddToEnd(0); m0.Canonicalize();
                        m1.AddToEnd(1); m1.Canonicalize();
                        n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m0);
                        n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode[0], m1);
                        n->numExits = 2;
                    }
                    else if (nhL->rweight == 0)
                    {
                        CFLOBDDReturnMapHandle m0, m1;
                        m0.AddToEnd(0); m0.Canonicalize();
                        m1.AddToEnd(1); m1.Canonicalize();
                        n->BConnection[0] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode[0], m0);
                        n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[0], m1);
                        n->numExits = 2; 
                    }
                }
            }
            else
            {
                int level = nh.handleContents->level;
                auto tempHandle = VectorToMatrixInterleavedNode(hashMap, *(nhNode->AConnection.entryPointHandle));
                CFLOBDDReturnMapHandle mA;
                for (unsigned int i = 0; i < tempHandle.handleContents->numExits; i++)
                    mA.AddToEnd(i);
                mA.Canonicalize();
                n->AConnection = Connection(tempHandle, mA);
                n->numBConnections = mA.Size();
                n->BConnection = new Connection[n->numBConnections];
                for (int i = 0; i < n->numBConnections; i++)
                {
                    auto b = VectorToMatrixInterleavedNode(hashMap, *(nhNode->BConnection[i].entryPointHandle));
                    CFLOBDDReturnMapHandle mX;
                    for (int j = 0; j < nhNode->BConnection[i].returnMapHandle.Size(); j++)
                    {
                        int v = nhNode->BConnection[i].returnMapHandle[j];
                        mX.AddToEnd(v);
                    }
                    mX.Canonicalize();
                    n->BConnection[i] = Connection(b, mX);
                }
                n->numExits = nhNode->numExits;
            }
    #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
    #endif

            WeightedCFLOBDDComplexFloatBoostMulNodeHandle nHandle(n);
            // hashMap->Insert(nh, nHandle);
            hashMap[nh] = nHandle;
            return nHandle;
        }

        WeightedCFLOBDDComplexFloatBoostMulNodeHandle MkColumn1MatrixNode(unsigned int level)
        {
            WeightedCFLOBDDComplexFloatBoostInternalNode *n = new WeightedCFLOBDDComplexFloatBoostInternalNode(level);
            assert(level > 0);
            if (level == 1)
            {
                CFLOBDDReturnMapHandle m1, m2;
                m1.AddToEnd(0);
                m1.Canonicalize();

                m2.AddToEnd(0);
                m2.AddToEnd(1);
                m2.Canonicalize();

                n->AConnection = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::CFLOBDDDontCareNodeHandle, m1);
                n->numBConnections = 1;
                n->BConnection = new Connection[1];
                WeightedCFLOBDDComplexFloatBoostMulNodeHandle t = WeightedCFLOBDDComplexFloatBoostMulNodeHandle(new WeightedCFLOBDDComplexFloatBoostForkNode(1, 0));
                n->BConnection[0] = Connection(t, m2);
            }
            else
            {
                CFLOBDDReturnMapHandle m1, m2;
                m1.AddToEnd(0);
                m1.AddToEnd(1);
                m1.Canonicalize();

                m2.AddToEnd(1);
                m2.Canonicalize();

                WeightedCFLOBDDComplexFloatBoostMulNodeHandle tempHandle = MkColumn1MatrixNode(level - 1);
                n->AConnection = Connection(tempHandle, m1);
                n->numBConnections = 2;
                n->BConnection = new Connection[n->numBConnections];
                n->BConnection[0] = Connection(tempHandle, m1);
                n->BConnection[1] = Connection(WeightedCFLOBDDComplexFloatBoostMulNodeHandle::NoDistinctionNode_Ann[level-1],m2);
            }
            n->numExits = 2;
    #ifdef PATH_COUNTING_ENABLED
            n->InstallPathCounts();
    #endif
            return WeightedCFLOBDDComplexFloatBoostMulNodeHandle(n);
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
        

    //#ifdef PATH_COUNTING_ENABLED
        std::pair<std::string,std::string> SamplingNode(WeightedCFLOBDDComplexFloatBoostMulNodeHandle nh, unsigned int index, std::mt19937 mt, std::uniform_real_distribution<double> dis, bool VocTwo,  std::string func)
        {
            if (nh.handleContents->NodeKind() == W_BDD_TOPNODE)
            {
                WeightedBDDComplexFloatBoostTopNode* node = (WeightedBDDComplexFloatBoostTopNode *)nh.handleContents;
                return WeightedVectorBDDComplexFloatBoostMul::SamplingNode(node->bddContents, node->numberOfVars, mt, dis);
            }
            WeightedCFLOBDDComplexFloatBoostInternalNode *nhNode = (WeightedCFLOBDDComplexFloatBoostInternalNode *)nh.handleContents;

            if (nhNode->level == 0)
            {
                if (nhNode->numExits == 2)
                {
                    return std::make_pair(std::to_string(index),"");
                }
                else
                {
                    WeightedCFLOBDDComplexFloatBoostDontCareNode* nhL = (WeightedCFLOBDDComplexFloatBoostDontCareNode *)nh.handleContents;
                    assert(index == 0);
                    long double lw = ((nhL->lweight * nhL->lweight).convert_to<long double>() / nhL->numWeightsOfPathsAsAmpsToExit[0]);
                    long double rw = ((nhL->rweight * nhL->rweight).convert_to<long double>() / nhL->numWeightsOfPathsAsAmpsToExit[0]);
                    std::vector<std::pair<long double, unsigned int>> weights = {std::make_pair(lw, 0), std::make_pair(rw, 1)};
                    // sort(weights.begin(), weights.end(), sortNumPathPairs<long double>);
                    double random_value = ((double)rand() / (RAND_MAX));
                    int chosen_index = -1;
                    if (func != "Grovers")
                        chosen_index = chooseIndexRandomly(weights, random_value);
                    if (func == "Grovers")
                    {
                        chosen_index = lw > rw ? 0 : 1;
                    }
                    return std::make_pair(std::to_string(chosen_index), "");
                }
            }
            std::vector<std::pair<long double, unsigned int>> numBPaths;
            long double numBTotalPaths = 0.0;
            for (unsigned int i = 0; i < nhNode->numBConnections; i++)
            {
                int BIndex = nhNode->BConnection[i].returnMapHandle.LookupInv(index);
                // if (BIndex == -1){
                //     numBPaths.push_back(std::make_pair(0, i));
                // }
                // else{
                    if (BIndex != -1){
                        long double w = ((nhNode->BConnection[i].entryPointHandle->handleContents->numWeightsOfPathsAsAmpsToExit[BIndex] * 
                            nhNode->AConnection.entryPointHandle->handleContents->numWeightsOfPathsAsAmpsToExit[i]) / nhNode->numWeightsOfPathsAsAmpsToExit[index]);
                        numBPaths.push_back(std::make_pair(w, i));
                    }
                // }
            }
            // sort(numBPaths.begin(), numBPaths.end(), sortNumPathPairs<long double>);
            /*std::cout << numBPaths.size() << std::endl;
            for (unsigned int i = 0; i < numBPaths.size(); i++)
            {
                std::cout << i << " " << numBPaths[i].first << " " << numBPaths[i].second << std::endl;
            }*/
            int BConnectionIndex = -1;
            if (func != "Grovers"){
                long double random_value = ((double)rand() / RAND_MAX);
                
                BConnectionIndex = chooseIndexRandomly(numBPaths, random_value);
            }
            //if (BConnectionIndex == -1){
            //	std::cout << val << " " << random_value << " " << index << " " << nhNode->level << std::endl;
            //	for (unsigned int i = 0; i < numBPaths.size(); i++)
            //	{
            //		std::cout << numBPaths[i].first << " " << numBPaths[i].second << std::endl;
            //	}
            //	//BConnectionIndex = numBPaths.back().second;
            //}

            if (func == "Grovers")
            {
                const auto p = max_element(numBPaths.begin(), numBPaths.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
                BConnectionIndex = p->second;
            }
            assert(BConnectionIndex != -1);
            assert(nhNode->BConnection[BConnectionIndex].returnMapHandle.LookupInv(index) != -1);
            assert(BConnectionIndex < nhNode->numBConnections);
            assert(nhNode->BConnection[BConnectionIndex].returnMapHandle.LookupInv(index) < nhNode->BConnection[BConnectionIndex].returnMapHandle.mapContents->mapArray.size());
            std::pair<std::string,std::string> AString = SamplingNode(*(nhNode->AConnection.entryPointHandle), BConnectionIndex, mt, dis, VocTwo, func);
            std::pair<std::string, std::string> BString = SamplingNode(*(nhNode->BConnection[BConnectionIndex].entryPointHandle), nhNode->BConnection[BConnectionIndex].returnMapHandle.LookupInv(index), mt, dis, VocTwo, func);
            if (nhNode->level == 1)
                return std::make_pair(AString.first, BString.first);
            if (nhNode->level == 2 && !VocTwo)
                return std::make_pair(AString.first, BString.first);
            return std::make_pair(AString.first + BString.first, AString.second + BString.second);
            //return std::make_pair(AString.first + BString.first + AString.second + BString.second, "");
        }

        long double getNonZeroProbabilityNode(WeightedCFLOBDDComplexFloatBoostMulNodeHandle nh)
        {
            if (nh.handleContents->NodeKind() == W_BDD_TOPNODE)
            {
                WeightedBDDComplexFloatBoostTopNode* nNode = (WeightedBDDComplexFloatBoostTopNode *)nh.handleContents;
                return nNode->bddContents.handleContents->weightOfPathsAsAmpsToExit;
            }
            return 0;
        }
    }
//#endif
}  // namespace CFL_OBDD
