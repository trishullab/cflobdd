
#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include "weighted_bdd_node_t.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "inttriple.h"
#include "hash.h"
#include "hashset.h"
#include "weighted_cross_product_bdd.h"
#include "weighted_values.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/functional/hash.hpp>
#include "fourier_semiring.h"

// using namespace CFL_OBDD;

// ********************************************************************
// 2-Way Cross Product
// ********************************************************************

//***************************************************************
// WeightedBDDPairProductMapBody
//***************************************************************
namespace CFL_OBDD {

    // Initializations of static members ---------------------------------
    template <typename T>
    Hashset<WeightedBDDPairProductMapBody<T>> *WeightedBDDPairProductMapBody<T>::canonicalWeightedBDDPairProductMapBodySet = new Hashset<WeightedBDDPairProductMapBody<T>>(HASHSET_NUM_BUCKETS);

    // Constructor
    template <typename T>
    WeightedBDDPairProductMapBody<T>::WeightedBDDPairProductMapBody()
    : refCount(0), isCanonical(false)
    {
    }

    template <typename T>
    void WeightedBDDPairProductMapBody<T>::IncrRef()
    {
        refCount++;    // Warning: Saturation not checked
    }

    template <typename T>
    void WeightedBDDPairProductMapBody<T>::DecrRef()
    {
        if (--refCount == 0) {    // Warning: Saturation not checked
            if (isCanonical) {
            WeightedBDDPairProductMapBody<T>::canonicalWeightedBDDPairProductMapBodySet->DeleteEq(this);
            }
            delete this;
        }
    }

    template <typename T>
    unsigned int WeightedBDDPairProductMapBody<T>::Hash(unsigned int modsize)
    {
        unsigned int hvalue = 0;
        boost::hash<T> boost_hash;
        for (unsigned int i = 0; i < mapArray.size(); i++){
            hvalue = (997*hvalue + (unsigned int)97*boost_hash(mapArray[i])) % modsize;
        }
        hvalue = (117 * hvalue) + boost_hash(factor);

        return hvalue;
    }

    template <typename T>
    void WeightedBDDPairProductMapBody<T>::setHashCheck()
    {
        unsigned int hvalue = 0;
        boost::hash<T> boost_hash;
        for (auto &i : mapArray) {
            hvalue = (117 * (hvalue + 1) + (int)(97 * boost_hash(i)));
        }
        hvalue = (117 * hvalue) + boost_hash(factor); 
        hashCheck = hvalue;
    }

    template <typename T>
    void WeightedBDDPairProductMapBody<T>::AddToEnd(const std::vector<T>& y, const T v)
    {
        for (auto i : y)
            mapArray.push_back(i);
        
        factor = v;
    }

    template <typename T>
    void WeightedBDDPairProductMapBody<T>::AddToEnd(const std::vector<T>& y)
    {
        for (auto i : y)
            mapArray.push_back(i);
        // change this
        factor = 1.0;
    }

    template <>
    void WeightedBDDPairProductMapBody<fourierSemiring>::AddToEnd(const std::vector<fourierSemiring>& y)
    {
        for (auto i : y)
            mapArray.push_back(i);
        fourierSemiring one(1, 1);
        factor = one;
    }

    template <typename T>
    bool WeightedBDDPairProductMapBody<T>::operator==(const WeightedBDDPairProductMapBody<T> &o) const
    {
        if (mapArray.size() != o.mapArray.size())
            return false;

        for (unsigned int i = 0; i < mapArray.size(); i++){
            if (mapArray[i] != o.mapArray[i])
                return false;
        }
        if (factor != o.factor)
            return false;
        return true;
    }
    
    // template <typename T>
    // std::pair<intpair, Pair_T<T,T>> WeightedBDDPairProductMapBody<T>::operator[](unsigned int i){                       // Overloaded []
    //     return std::make_pair(mapArray[i],valueArray[i]);
    // }

    template <typename T>
    unsigned int WeightedBDDPairProductMapBody<T>::Size(){
        return (unsigned int)mapArray.size();
    }


    // namespace CFL_OBDD {
    template <typename T>
    std::ostream& operator<< (std::ostream & out, const WeightedBDDPairProductMapBody<T> &r)
    {
    //out << (List<int>&)r;
        for (unsigned int i = 0; i < r.mapArray.size(); i++)
        {
            out << r.mapArray[i] << ", ";
        }
        out << std::endl;
        out << "factor: " << r.factor << std::endl;
        return(out);
    }
    // }
    //***************************************************************
    // WeightedBDDPairProductMapHandle
    //***************************************************************


    // Default constructor
    template <typename T>
    WeightedBDDPairProductMapHandle<T>::WeightedBDDPairProductMapHandle()
    :  mapContents(new WeightedBDDPairProductMapBody<T>())
    {
        mapContents->IncrRef();
    }

    // Destructor
    template <typename T>
    WeightedBDDPairProductMapHandle<T>::~WeightedBDDPairProductMapHandle()
    {
        mapContents->DecrRef();
    }

    // Copy constructor
    template <typename T>
    WeightedBDDPairProductMapHandle<T>::WeightedBDDPairProductMapHandle(const WeightedBDDPairProductMapHandle<T> &r)
    :  mapContents(r.mapContents)
    {
        mapContents->IncrRef();
    }

    // Overloaded assignment
    template <typename T>
    WeightedBDDPairProductMapHandle<T>& WeightedBDDPairProductMapHandle<T>::operator= (const WeightedBDDPairProductMapHandle<T> &r)
    {
        if (this != &r)      // don't assign to self!
        {
            WeightedBDDPairProductMapBody<T> *temp = mapContents;
            mapContents = r.mapContents;
            mapContents->IncrRef();
            temp->DecrRef();
        }
        return *this;        
    }

    // Overloaded !=
    template <typename T>
    bool WeightedBDDPairProductMapHandle<T>::operator!=(const WeightedBDDPairProductMapHandle<T> &r)
    {
        return (mapContents != r.mapContents);
    }

    // Overloaded ==
    template <typename T>
    bool WeightedBDDPairProductMapHandle<T>::operator==(const WeightedBDDPairProductMapHandle<T> &r)
    {
        return (mapContents == r.mapContents);
    }

    template <typename T>
    std::ostream& operator<< (std::ostream & out, const WeightedBDDPairProductMapHandle<T> &r)
    {
        out << "[" << *r.mapContents << "]";
        return(out);
    }

    template <typename T>
    unsigned int WeightedBDDPairProductMapHandle<T>::Hash(unsigned int modsize)
    {
        return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
    }

    template <typename T>
    unsigned int WeightedBDDPairProductMapHandle<T>::Size()
    {
        return mapContents->Size();
    }

    // template <typename T>
    // std::pair<intpair, Pair_T<T,T>> WeightedBDDPairProductMapHandle<T>::operator[](unsigned int i)
    // {
    //     return std::make_pair(mapContents->mapArray[i], mapContents->valueArray[i]);
    // }

    template <typename T>
    void WeightedBDDPairProductMapHandle<T>::AddToEnd(const std::vector<T>& p, const T v)
    {
        assert(mapContents->refCount <= 1);
        mapContents->AddToEnd(p, v);
    }

    template <typename T>
    void WeightedBDDPairProductMapHandle<T>::AddToEnd(const std::vector<T>& p)
    {
        assert(mapContents->refCount <= 1);
        mapContents->AddToEnd(p);
    }

    // template <typename T>
    // bool WeightedBDDPairProductMapHandle<T>::Member(intpair& p)
    // {
    //     for (auto& i : mapContents->mapArray){
    //         if (i == p)
    //             return true;
    //     }
    //     return false;
    // }

    // template <typename T>
    // int WeightedBDDPairProductMapHandle<T>::Lookup(intpair& p, Pair_T<T,T>& v)
    // {
    //     for (unsigned int i = 0; i < mapContents->mapArray.size(); i++){
    //         if (mapContents->mapArray[i] == p && mapContents->valueArray[i] == v)
    //             return i;
    //     }
    //     return -1;
    // }

    template <typename T>
    void WeightedBDDPairProductMapHandle<T>::Canonicalize()
    {
        WeightedBDDPairProductMapBody<T> *answerContents;
        unsigned int hash = WeightedBDDPairProductMapBody<T>::canonicalWeightedBDDPairProductMapBodySet->GetHash(mapContents);
        answerContents = WeightedBDDPairProductMapBody<T>::canonicalWeightedBDDPairProductMapBodySet->Lookup(mapContents, hash);
        if (answerContents == NULL) {
            WeightedBDDPairProductMapBody<T>::canonicalWeightedBDDPairProductMapBodySet->Insert(mapContents, hash);
            mapContents->isCanonical = true;
        }
        else {
            answerContents->IncrRef();
            mapContents->DecrRef();
            mapContents = answerContents;
        }
    }

    // Create map with reversed entries
    // template <typename T>
    // WeightedBDDPairProductMapHandle<T> WeightedBDDPairProductMapHandle<T>::Flip()
    // {
    //     WeightedBDDPairProductMapHandle<T> answer;
    //     for (int i = 0; i < mapContents->Size(); i++){
    //         intpair p(mapContents->mapArray[i].Second(), mapContents->mapArray[i].First());
    //         Pair_T<T,T> v(mapContents->valueArray[i].Second(), mapContents->valueArray[i].First());
    //         answer.AddToEnd(p, v);
    //     }
    //     return answer;
    // }

    //***************************************************************
    // WeightedBDDPairProductKey
    //***************************************************************

    // Constructor
    template <typename T, typename Op>
    WeightedBDDPairProductKey<T,Op>::WeightedBDDPairProductKey(WeightedBDDNodeHandle<T,Op> nodeHandle1, WeightedBDDNodeHandle<T,Op> nodeHandle2)
    :  nodeHandle1(nodeHandle1), nodeHandle2(nodeHandle2)
    {
        factor1 = getIdentityValue<T,Op>();
        factor2 = getIdentityValue<T,Op>();
    }

    template <typename T, typename Op>
    WeightedBDDPairProductKey<T,Op>::WeightedBDDPairProductKey(WeightedBDDNodeHandle<T,Op> nodeHandle1, WeightedBDDNodeHandle<T,Op> nodeHandle2, T factor1, T factor2)
    :  nodeHandle1(nodeHandle1), nodeHandle2(nodeHandle2), factor1(factor1), factor2(factor2)
    {
    }

    // Hash
    template <typename T, typename Op>
    unsigned int WeightedBDDPairProductKey<T,Op>::Hash(unsigned int modsize)
    {
        unsigned int hvalue = 0;
        boost::hash<T> boost_hash;
        hvalue = (997 * nodeHandle1.Hash(modsize) + nodeHandle2.Hash(modsize)) % modsize;
        hvalue = (hvalue + (117 * boost_hash(factor1) % modsize + 17 * boost_hash(factor2) % modsize)) % modsize;
        return hvalue;
    }

    // print
    template <typename T, typename Op>
    std::ostream& WeightedBDDPairProductKey<T,Op>::print(std::ostream & out) const
    {
        out << "(" << nodeHandle1 << ", " << nodeHandle2 << "), (" << factor1 << ", " << factor2 << ")";
        return out;
    }

    template <typename T, typename Op>
    std::ostream& operator<< (std::ostream & out, const WeightedBDDPairProductKey<T,Op> &p)
    {
        p.print(out);
        return(out);
    }

    template <typename T, typename Op>
    WeightedBDDPairProductKey<T,Op>& WeightedBDDPairProductKey<T,Op>::operator= (const WeightedBDDPairProductKey<T,Op>& i)
    {
        if (this != &i)      // don't assign to self!
        {
            nodeHandle1 = i.nodeHandle1;
            nodeHandle2 = i.nodeHandle2;
            factor1 = i.factor1;
            factor2 = i.factor2;
        }
        return *this;        
    }

    // Overloaded !=
    template <typename T, typename Op>
    bool WeightedBDDPairProductKey<T,Op>::operator!=(const WeightedBDDPairProductKey<T,Op>& p)
    {
        return (factor1 != p.factor1) || (factor2 != p.factor2) || (nodeHandle1 != p.nodeHandle1) || (nodeHandle2 != p.nodeHandle2);
    }

    // Overloaded ==
    template <typename T, typename Op>
    bool WeightedBDDPairProductKey<T,Op>::operator==(const WeightedBDDPairProductKey<T,Op>& p)
    {
        return (factor1 == p.factor1) && (factor2 == p.factor2) && (nodeHandle1 == p.nodeHandle1) && (nodeHandle2 == p.nodeHandle2);
    }

    //***************************************************************
    // WeightedBDDPairProductMemo
    //***************************************************************

    // Default constructor
    template <typename T, typename Op>
    WeightedBDDPairProductMemo<T,Op>::WeightedBDDPairProductMemo()
    :  nodeHandle(WeightedBDDNodeHandle<T,Op>()), pairProductMapHandle(WeightedBDDPairProductMapHandle<T>())
    {
    }

    // Constructor
    template <typename T, typename Op>
    WeightedBDDPairProductMemo<T,Op>::WeightedBDDPairProductMemo(WeightedBDDNodeHandle<T,Op> nodeHandle, WeightedBDDPairProductMapHandle<T> pairProductMapHandle)
    :  nodeHandle(nodeHandle), pairProductMapHandle(pairProductMapHandle)
    {
    }

    template <typename T, typename Op>
    std::ostream& operator<< (std::ostream & out, const WeightedBDDPairProductMemo<T,Op> &p)
    {
        out << "(" << p.nodeHandle << ", " << p.pairProductMapHandle << ")";
        return(out);
    }

    template <typename T, typename Op>
    WeightedBDDPairProductMemo<T,Op>& WeightedBDDPairProductMemo<T,Op>::operator= (const WeightedBDDPairProductMemo<T,Op>& i)
    {
        if (this != &i)      // don't assign to self!
        {
            nodeHandle = i.nodeHandle;
            pairProductMapHandle = i.pairProductMapHandle;
        }
        return *this;        
    }

    // Overloaded !=
    template <typename T, typename Op>
    bool WeightedBDDPairProductMemo<T,Op>::operator!=(const WeightedBDDPairProductMemo<T,Op>& p)
    {
        return (nodeHandle != p.nodeHandle) || (pairProductMapHandle != p.pairProductMapHandle);
    }

    // Overloaded ==
    template <typename T, typename Op>
    bool WeightedBDDPairProductMemo<T,Op>::operator==(const WeightedBDDPairProductMemo<T,Op>& p)
    {
        return (nodeHandle == p.nodeHandle) && (pairProductMapHandle == p.pairProductMapHandle);
    }

    // --------------------------------------------------------------------
    // PairProduct
    //
    // Returns a new CFLOBDDNodeHandle, and (in WeightedBDDPairProductMapHandle) a descriptor of the
    // node's exits
    // --------------------------------------------------------------------

    // TODO: May change to something better
    template <typename T>
    class FourTuple {
        public:
            FourTuple(int c1, int c2, T f1, T f2) : c1(c1), c2(c2), f1(f1), f2(f2) {};
            FourTuple& operator= (const FourTuple& p)
            {
                if (this != &p)
                {
                    c1 = p.c1;
                    c2 = p.c2;
                    f1 = p.f1;
                    f2 = p.f2;
                }
                return *this;
            }
            bool operator== (const FourTuple& p) const
            {
                return (c1 == p.c1) && (c2 == p.c2) && (f1 == p.f1) && (f2 == p.f2);
            }
            struct FourTuple_Hash {
                size_t operator() (const FourTuple<T>& p) const
                {
                    boost::hash<T> boost_hash;
                    return (117 * p.c1 + p.c2) + (997 * boost_hash(p.f1) + 97 * boost_hash(p.f2)) % 2019;
                }
            };
            int c1;
            int c2;
            T f1;
            T f2;
    };
    

    template <typename T, typename Op>
    static Hashtable<WeightedBDDPairProductKey<T,Op>, WeightedBDDPairProductMemo<T,Op>> *pairProductCache = NULL;
    template <typename T, typename Op>
    static Hashtable<WeightedBDDPairProductKey<T,Op>, WeightedBDDPairProductMemo<T,Op>> *pairProduct2Cache = NULL;

    // namespace CFL_OBDD {

    template <typename T, typename Op>
    WeightedBDDNodeHandle<T,Op> BDDPairProduct(WeightedBDDInternalNode<T,Op> *n1,
                                WeightedBDDInternalNode<T,Op> *n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<T> &pairProductMapHandle,
                                T(*func)(T, T)
                                )
    {
        WeightedBDDNodeHandle<T,Op> leftNode1, rightNode1, leftNode2, rightNode2;
        T lweight, rweight;
        long int index;
        if (n1->GetIndex() == n2->GetIndex()){
            leftNode1 = n1->leftNode;
            rightNode1 = n1->rightNode;
            leftNode2 = n2->leftNode;
            rightNode2 = n2->rightNode;
            lweight = computeComposition<T,Op>(n1->lweight, n2->lweight);
            rweight = computeComposition<T,Op>(n1->rweight, n2->rweight);
            index = n1->GetIndex();
        } else if (n1->GetIndex() > n2->GetIndex()) {
            leftNode1 = n1;
            rightNode1 = n1;
            leftNode2 = n2->leftNode;
            rightNode2 = n2->rightNode;
            lweight = n2->lweight;
            rweight = n2->rweight;
            index = n2->GetIndex();
        } else if (n1->GetIndex() < n2->GetIndex()) {
            leftNode1 = n1->leftNode;
            rightNode1 = n1->rightNode;
            leftNode2 = n2;
            rightNode2 = n2;
            lweight = n1->lweight;
            rweight = n1->rweight;
            index = n1->GetIndex();
        }
        WeightedBDDPairProductMapHandle<T> lpairProductMapHandle;
        WeightedBDDPairProductMapHandle<T> rpairProductMapHandle;
        auto lvalue = BDDPairProduct<T,Op>(leftNode1,
                                    leftNode2,
                                    numVars,
                                    lpairProductMapHandle,
                                    func);
        auto rvalue = BDDPairProduct<T,Op>(rightNode1,
                                    rightNode2,
                                    numVars,
                                    rpairProductMapHandle,
                                    func); 
        WeightedBDDInternalNode<T,Op>* g = new WeightedBDDInternalNode<T,Op>(index);
        lweight = computeComposition<T,Op>(lweight, lpairProductMapHandle.mapContents->factor);
        rweight = computeComposition<T,Op>(rweight, rpairProductMapHandle.mapContents->factor);
        g->lweight = lweight;
        g->rweight = rweight;
        g->leftNode = lvalue;
        g->rightNode = rvalue;
        if (g->leftNode == g->rightNode && g->lweight == g->rweight){
            pairProductMapHandle = lpairProductMapHandle;
            return g->leftNode;
        }
        std::vector<T> ret_ans;
        auto lvector = lpairProductMapHandle.mapContents->mapArray;
        auto rvector = rpairProductMapHandle.mapContents->mapArray;
        for (auto i : lvector)
            ret_ans.push_back(i);
        
        if (lvector.size() == 1 && rvector.size() == 1)
        {
            if (ret_ans[0] != rvector[0])
                ret_ans.push_back(rvector[0]);
        }
        else if (lvector.size() == 1 && rvector.size() == 2)
        {
            for (auto i : rvector)
            {
                if (ret_ans[0] != i)
                    ret_ans.push_back(i);
            }
        }
        pairProductMapHandle.AddToEnd(ret_ans);
        return WeightedBDDNodeHandle<T,Op>(g);

    }

    template <typename T, typename Op>
    WeightedBDDNodeHandle<T,Op> BDDPairProduct(WeightedBDDNodeHandle<T,Op> n1,
                                WeightedBDDNodeHandle<T,Op> n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<T> &pairProductMapHandle,
                                T(*func)(T, T)
                                )
    {
        WeightedBDDPairProductMemo<T,Op> cachedWeightedBDDPairProductMemo;
        bool isCached = pairProductCache<T,Op>->Fetch(WeightedBDDPairProductKey(n1,n2), cachedWeightedBDDPairProductMemo);
        if (isCached) {
            pairProductMapHandle = cachedWeightedBDDPairProductMemo.pairProductMapHandle;
            return cachedWeightedBDDPairProductMemo.nodeHandle;
        }
        else if (pairProductCache<T,Op>->Fetch(WeightedBDDPairProductKey(n2,n1), cachedWeightedBDDPairProductMemo)) {
            pairProductMapHandle = cachedWeightedBDDPairProductMemo.pairProductMapHandle;
            // TODO: problem with vector
            return cachedWeightedBDDPairProductMemo.nodeHandle;
        }
        else {
            WeightedBDDNodeHandle<T,Op> answer;
        
            if (n1.handleContents->NodeKind() == INTERNAL && n2.handleContents->NodeKind() == INTERNAL) {
                answer = BDDPairProduct<T,Op>((WeightedBDDInternalNode<T,Op> *)n1.handleContents,
                                (WeightedBDDInternalNode<T,Op> *)n2.handleContents,
                                numVars,
                                pairProductMapHandle,
                                func
                                );
            }
            else if (n1.handleContents->NodeKind() == INTERNAL && n2.handleContents->NodeKind() == LEAF) {
                if (n2 == WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode) {
                    answer = n2;
                    std::vector<T> ret_ans;
                    ret_ans.push_back(getAnnhilatorValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
                else if (n2 == WeightedBDDNodeHandle<T,Op>::IdentityLeafNode) {
                    answer = n1;
                    std::vector<T> ret_ans;
                    // TODO: Incorrect but doesn't affect the algorithm
                    ret_ans.push_back(getIdentityValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
                else {
                    WeightedBDDInternalNode<T,Op>* n1_i = (WeightedBDDInternalNode<T,Op> *)n1.handleContents;
                    T lweight, rweight;
                    WeightedBDDPairProductMapHandle<T> lpairProductMapHandle;
                    WeightedBDDPairProductMapHandle<T> rpairProductMapHandle;
                    auto lvalue = BDDPairProduct<T,Op>(n1_i->leftNode, n2, numVars, lpairProductMapHandle, func);
                    auto rvalue = BDDPairProduct<T,Op>(n1_i->rightNode, n2, numVars, rpairProductMapHandle, func); 
                    WeightedBDDInternalNode<T,Op>* g = new WeightedBDDInternalNode<T,Op>(n1_i->GetIndex());
                    lweight = computeComposition<T,Op>(lweight, lpairProductMapHandle.mapContents->factor);
                    rweight = computeComposition<T,Op>(rweight, rpairProductMapHandle.mapContents->factor);
                    g->lweight = lweight;
                    g->rweight = rweight;
                    g->leftNode = lvalue;
                    g->rightNode = rvalue;
                    if (g->leftNode == g->rightNode && g->lweight == g->rweight){
                        answer = lvalue;
                        pairProductMapHandle = lpairProductMapHandle;
                    }
                    std::vector<T> ret_ans;
                    auto lvector = lpairProductMapHandle.mapContents->mapArray;
                    auto rvector = rpairProductMapHandle.mapContents->mapArray;
                    for (auto i : lvector)
                        ret_ans.push_back(i);
                    
                    if (lvector.size() == 1 && rvector.size() == 1)
                    {
                        if (ret_ans[0] != rvector[0])
                            ret_ans.push_back(rvector[0]);
                    }
                    else if (lvector.size() == 1 && rvector.size() == 2)
                    {
                        for (auto i : rvector)
                        {
                            if (ret_ans[0] != i)
                                ret_ans.push_back(i);
                        }
                    }
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
            }
            else if (n1.handleContents->NodeKind() == LEAF && n2.handleContents->NodeKind() == INTERNAL) { 
                /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // W_CFLOBDD_FORK, CFLOBDD_DONTCARE
                if (n1 == WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode) {
                    answer = n1;
                    std::vector<T> ret_ans;
                    ret_ans.push_back(getAnnhilatorValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
                else if (n1 == WeightedBDDNodeHandle<T,Op>::IdentityLeafNode) {
                    answer = n2;
                    std::vector<T> ret_ans;
                    // TODO: Incorrect but doesn't affect the algorithm
                    ret_ans.push_back(getIdentityValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
                else {
                    WeightedBDDInternalNode<T,Op>* n2_i = (WeightedBDDInternalNode<T,Op> *)n2.handleContents;
                    T lweight, rweight;
                    WeightedBDDPairProductMapHandle<T> lpairProductMapHandle;
                    WeightedBDDPairProductMapHandle<T> rpairProductMapHandle;
                    auto lvalue = BDDPairProduct<T,Op>(n1, n2_i->leftNode, numVars, lpairProductMapHandle, func);
                    auto rvalue = BDDPairProduct<T,Op>(n1, n2_i->rightNode, numVars, rpairProductMapHandle, func); 
                    WeightedBDDInternalNode<T,Op>* g = new WeightedBDDInternalNode<T,Op>(n2_i->GetIndex());
                    lweight = computeComposition<T,Op>(lweight, lpairProductMapHandle.mapContents->factor);
                    rweight = computeComposition<T,Op>(rweight, rpairProductMapHandle.mapContents->factor);
                    g->lweight = lweight;
                    g->rweight = rweight;
                    g->leftNode = lvalue;
                    g->rightNode = rvalue;
                    if (g->leftNode == g->rightNode && g->lweight == g->rweight){
                        answer = lvalue;
                        pairProductMapHandle = lpairProductMapHandle;
                    }
                    std::vector<T> ret_ans;
                    auto lvector = lpairProductMapHandle.mapContents->mapArray;
                    auto rvector = rpairProductMapHandle.mapContents->mapArray;
                    for (auto i : lvector)
                        ret_ans.push_back(i);
                    
                    if (lvector.size() == 1 && rvector.size() == 1)
                    {
                        if (ret_ans[0] != rvector[0])
                            ret_ans.push_back(rvector[0]);
                    }
                    else if (lvector.size() == 1 && rvector.size() == 2)
                    {
                        for (auto i : rvector)
                        {
                            if (ret_ans[0] != i)
                                ret_ans.push_back(i);
                        }
                    }
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
            }
            else { /* n1.handleContents->NodeKind() == CFLOBDD_DONTCARE */
                if (n1 == WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode || n2 == WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode)
                {
                    answer = WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode;
                    std::vector<T> ret_ans;
                    ret_ans.push_back(getAnnhilatorValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
                else {
                    answer = WeightedBDDNodeHandle<T,Op>::IdentityLeafNode;   
                    std::vector<T> ret_ans;
                    ret_ans.push_back(getIdentityValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans);
                }
            }
            
            pairProductCache<T,Op>->Insert(WeightedBDDPairProductKey(n1, n2),
                WeightedBDDPairProductMemo(answer, pairProductMapHandle));
            return answer;
        }    
    }


    template <typename T, typename Op>
    WeightedBDDNodeHandle<T,Op> BDDPairProduct2(WeightedBDDInternalNode<T,Op> *n1,
                                WeightedBDDInternalNode<T,Op> *n2,
                                unsigned int numVars,
                                T factor1,
                                T factor2,
                                WeightedBDDPairProductMapHandle<T> &pairProductMapHandle,
                                T(*func)(T, T)
                                )
    {   
        WeightedBDDNodeHandle<T,Op> leftNode1, rightNode1, leftNode2, rightNode2;
        T l1, r1, l2, r2;
        long int index;
        if (n1->GetIndex() == n2->GetIndex()){
            leftNode1 = n1->leftNode;
            rightNode1 = n1->rightNode;
            leftNode2 = n2->leftNode;
            rightNode2 = n2->rightNode;
            l1 = computeComposition<T,Op>(factor1, n1->lweight);
            r1 = computeComposition<T,Op>(factor2, n2->lweight);
            l2 = computeComposition<T,Op>(factor1, n1->rweight);
            r2 = computeComposition<T,Op>(factor2, n2->rweight);
            index = n1->GetIndex();
        } else if (n1->GetIndex() > n2->GetIndex()) {
            leftNode1 = n1;
            rightNode1 = n1;
            leftNode2 = n2->leftNode;
            rightNode2 = n2->rightNode;
            l1 = factor1;
            r1 = computeComposition<T,Op>(factor2, n2->lweight);
            l2 = factor1;
            r2 = computeComposition<T,Op>(factor2, n2->rweight);
            index = n2->GetIndex();
        } else if (n1->GetIndex() < n2->GetIndex()) {
            leftNode1 = n1->leftNode;
            rightNode1 = n1->rightNode;
            leftNode2 = n2;
            rightNode2 = n2;
            l1 = computeComposition<T,Op>(factor1, n1->lweight);
            r1 = factor2;
            l2 = computeComposition<T,Op>(factor1, n1->rweight);
            r2 = factor2;
            index = n1->GetIndex();
        }
        WeightedBDDPairProductMapHandle<T> lpairProductMapHandle;
        WeightedBDDPairProductMapHandle<T> rpairProductMapHandle;
        auto lvalue = BDDPairProduct2<T,Op>(leftNode1,
                                    leftNode2,
                                    numVars,
                                    l1,
                                    r1,
                                    lpairProductMapHandle,
                                    func);
        auto rvalue = BDDPairProduct2<T,Op>(rightNode1,
                                    rightNode2,
                                    numVars,
                                    l2,
                                    r2,
                                    rpairProductMapHandle,
                                    func); 
        WeightedBDDInternalNode<T,Op>* g = new WeightedBDDInternalNode<T,Op>(index);
        T lweight = lpairProductMapHandle.mapContents->factor;
        T rweight = rpairProductMapHandle.mapContents->factor;
        auto w = computeInverseValue<T, Op>(lweight, rweight);
        g->lweight = std::get<1>(w);
        g->rweight = std::get<2>(w);
        g->leftNode = lvalue;
        g->rightNode = rvalue;
        if (g->leftNode == g->rightNode && g->lweight == g->rweight){
            pairProductMapHandle = lpairProductMapHandle;
            return g->leftNode;
        }
        std::vector<T> ret_ans;
        auto lvector = lpairProductMapHandle.mapContents->mapArray;
        auto rvector = rpairProductMapHandle.mapContents->mapArray;
        for (auto i : lvector)
            ret_ans.push_back(i);
        
        if (lvector.size() == 1 && rvector.size() == 1)
        {
            if (ret_ans[0] != rvector[0])
                ret_ans.push_back(rvector[0]);
        }
        else if (lvector.size() == 1 && rvector.size() == 2)
        {
            for (auto i : rvector)
            {
                if (ret_ans[0] != i)
                    ret_ans.push_back(i);
            }
        }
        pairProductMapHandle.AddToEnd(ret_ans, std::get<0>(w));
        // std::cout << "hey " << g->GetIndex() << std::endl;
        // g->print(std::cout);
        return WeightedBDDNodeHandle<T,Op>(g);
    }


    template <typename T, typename Op>
    WeightedBDDNodeHandle<T,Op> BDDPairProduct2(WeightedBDDNodeHandle<T,Op> n1,
                                WeightedBDDNodeHandle<T,Op> n2,
                                unsigned int numVars,
                                T factor1,
                                T factor2,
                                WeightedBDDPairProductMapHandle<T> &pairProductMapHandle,
                                T(*func)(T, T)
                                )
    {
        WeightedBDDPairProductMemo<T,Op> cachedWeightedBDDPairProductMemo;
        bool isCached = pairProduct2Cache<T,Op>->Fetch(WeightedBDDPairProductKey(n1,n2,factor1,factor2), cachedWeightedBDDPairProductMemo);
        if (isCached) {
            pairProductMapHandle = cachedWeightedBDDPairProductMemo.pairProductMapHandle;
            return cachedWeightedBDDPairProductMemo.nodeHandle;
        }
        else if (pairProduct2Cache<T,Op>->Fetch(WeightedBDDPairProductKey(n2,n1,factor2,factor1), cachedWeightedBDDPairProductMemo)) {
            pairProductMapHandle = cachedWeightedBDDPairProductMemo.pairProductMapHandle;
            // TODO: problem with the vector
            return cachedWeightedBDDPairProductMemo.nodeHandle;
        }
        else {
            WeightedBDDNodeHandle<T,Op> answer;
        
            if (n1.handleContents->NodeKind() == INTERNAL && n2.handleContents->NodeKind() == INTERNAL) {
                answer = BDDPairProduct2<T,Op>((WeightedBDDInternalNode<T,Op> *)n1.handleContents,
                                (WeightedBDDInternalNode<T,Op> *)n2.handleContents,
                                numVars,
                                factor1,
                                factor2,
                                pairProductMapHandle,
                                func
                                );
            }
            else if (n1.handleContents->NodeKind() == INTERNAL && n2.handleContents->NodeKind() == LEAF) {
                if (n2 == WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode) {
                    answer = n1;
                    std::vector<T> ret_ans;
                    // TODO: Incorrect but doesn't affect the algorithm
                    ret_ans.push_back(getIdentityValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans, factor1);
                }
                else {
                    WeightedBDDInternalNode<T,Op>* n1_i = (WeightedBDDInternalNode<T,Op> *)n1.handleContents;
                    WeightedBDDLeafNode<T,Op>* n2_i = (WeightedBDDLeafNode<T,Op> *)n2.handleContents;
                    T lweight, rweight;
                    WeightedBDDPairProductMapHandle<T> lpairProductMapHandle;
                    WeightedBDDPairProductMapHandle<T> rpairProductMapHandle;

                    auto lvalue = BDDPairProduct2<T,Op>(n1_i->leftNode, n2, numVars, 
                                        computeComposition<T,Op>(n1_i->lweight, factor1),
                                        computeComposition<T,Op>(n2_i->value, factor2),
                                        lpairProductMapHandle,
                                        func);
                    auto rvalue = BDDPairProduct2<T,Op>(n1_i->rightNode, n2, numVars,
                                        computeComposition<T,Op>(n1_i->rweight, factor1),
                                        computeComposition<T,Op>(n2_i->value, factor2),
                                        rpairProductMapHandle,
                                        func); 
                    WeightedBDDInternalNode<T,Op>* g = new WeightedBDDInternalNode<T,Op>(n1_i->GetIndex());
                    lweight = lpairProductMapHandle.mapContents->factor;
                    rweight = rpairProductMapHandle.mapContents->factor;
                    auto w = computeInverseValue<T, Op>(lweight, rweight);
                    g->lweight = std::get<1>(w);
                    g->rweight = std::get<2>(w);
                    g->leftNode = lvalue;
                    g->rightNode = rvalue;
                    if (g->leftNode == g->rightNode && g->lweight == g->rweight){
                        answer = g->leftNode;
                        pairProductMapHandle = lpairProductMapHandle;
                    }
                    std::vector<T> ret_ans;
                    auto lvector = lpairProductMapHandle.mapContents->mapArray;
                    auto rvector = rpairProductMapHandle.mapContents->mapArray;
                    for (auto i : lvector)
                        ret_ans.push_back(i);
                    
                    if (lvector.size() == 1 && rvector.size() == 1)
                    {
                        if (ret_ans[0] != rvector[0])
                            ret_ans.push_back(rvector[0]);
                    }
                    else if (lvector.size() == 1 && rvector.size() == 2)
                    {
                        for (auto i : rvector)
                        {
                            if (ret_ans[0] != i)
                                ret_ans.push_back(i);
                        }
                    }
                    answer = WeightedBDDNodeHandle<T,Op>(g);
                    pairProductMapHandle.AddToEnd(ret_ans, std::get<0>(w));
                }
            }
            else if (n1.handleContents->NodeKind() == LEAF && n2.handleContents->NodeKind() == INTERNAL) { 
                if (n1 == WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode) {
                    answer = n2;
                    std::vector<T> ret_ans;
                    // TODO: Incorrect but doesn't affect the algorithm
                    ret_ans.push_back(getIdentityValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans, factor2);
                }
                else {
                    /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // W_CFLOBDD_FORK, CFLOBDD_DONTCARE
                    WeightedBDDInternalNode<T,Op>* n2_i = (WeightedBDDInternalNode<T,Op> *)n2.handleContents;
                    WeightedBDDLeafNode<T,Op>* n1_i = (WeightedBDDLeafNode<T,Op> *)n1.handleContents;
                    T lweight, rweight;
                    WeightedBDDPairProductMapHandle<T> lpairProductMapHandle;
                    WeightedBDDPairProductMapHandle<T> rpairProductMapHandle;

                    auto lvalue = BDDPairProduct2<T,Op>(n1, n2_i->leftNode, numVars, 
                                        computeComposition<T,Op>(n1_i->value, factor1),
                                        computeComposition<T,Op>(n2_i->lweight, factor2),
                                        lpairProductMapHandle,
                                        func);
                    auto rvalue = BDDPairProduct2<T,Op>(n1, n2_i->rightNode, numVars,
                                        computeComposition<T,Op>(n1_i->value, factor1),
                                        computeComposition<T,Op>(n2_i->rweight, factor2),
                                        rpairProductMapHandle,
                                        func); 
                    WeightedBDDInternalNode<T,Op>* g = new WeightedBDDInternalNode<T,Op>(n2_i->GetIndex());
                    lweight = lpairProductMapHandle.mapContents->factor;
                    rweight = rpairProductMapHandle.mapContents->factor;
                    auto w = computeInverseValue<T, Op>(lweight, rweight);
                    g->lweight = std::get<1>(w);
                    g->rweight = std::get<2>(w);
                    g->leftNode = lvalue;
                    g->rightNode = rvalue;
                    if (g->leftNode == g->rightNode && g->lweight == g->rweight){
                        answer = lvalue;
                        pairProductMapHandle = lpairProductMapHandle;
                    }
                    std::vector<T> ret_ans;
                    auto lvector = lpairProductMapHandle.mapContents->mapArray;
                    auto rvector = rpairProductMapHandle.mapContents->mapArray;
                    for (auto i : lvector)
                        ret_ans.push_back(i);
                    
                    if (lvector.size() == 1 && rvector.size() == 1)
                    {
                        if (ret_ans[0] != rvector[0])
                            ret_ans.push_back(rvector[0]);
                    }
                    else if (lvector.size() == 1 && rvector.size() == 2)
                    {
                        for (auto i : rvector)
                        {
                            if (ret_ans[0] != i)
                                ret_ans.push_back(i);
                        }
                    }
                    answer = WeightedBDDNodeHandle<T,Op>(g);
                    pairProductMapHandle.AddToEnd(ret_ans, std::get<0>(w));
                }
            }
            else { /* n1.handleContents->NodeKind() == CFLOBDD_DONTCARE */
                WeightedBDDLeafNode<T,Op>* n1_i = (WeightedBDDLeafNode<T,Op> *)n1.handleContents;
                WeightedBDDLeafNode<T,Op>* n2_i = (WeightedBDDLeafNode<T,Op> *)n2.handleContents;
                T lweight = computeComposition<T,Op>(factor1, n1_i->value);
                T rweight = computeComposition<T,Op>(factor2, n2_i->value);
                T weight = func(lweight, rweight);
                if (weight == getAnnhilatorValue<T,Op>())
                {
                    answer = WeightedBDDNodeHandle<T,Op>::AnnhilatorLeafNode;
                    std::vector<T> ret_ans;
                    ret_ans.push_back(weight);
                    pairProductMapHandle.AddToEnd(ret_ans, weight);
                }
                else {
                    answer = WeightedBDDNodeHandle<T,Op>::IdentityLeafNode;
                    std::vector<T> ret_ans;
                    ret_ans.push_back(getIdentityValue<T,Op>());
                    pairProductMapHandle.AddToEnd(ret_ans, weight); 
                }
            }
            
            pairProductCache<T,Op>->Insert(WeightedBDDPairProductKey(n1, n2),
                WeightedBDDPairProductMemo(answer, pairProductMapHandle));
            // std::cout << "hey " << answer.handleContents->GetIndex() << std::endl;
            // answer.print(std::cout);
            return answer;
        }
    }

    template <typename T, typename Op>
    void InitWeightedBDDPairProductCache()
    {
        pairProductCache<T,Op> = new Hashtable<WeightedBDDPairProductKey<T,Op>, WeightedBDDPairProductMemo<T,Op>>(HASH_NUM_BUCKETS);
        pairProduct2Cache<T,Op> = new Hashtable<WeightedBDDPairProductKey<T,Op>, WeightedBDDPairProductMemo<T,Op>>(HASH_NUM_BUCKETS);
    }

    template <typename T, typename Op>
    void DisposeOfWeightedBDDPairProductCache()
    {
        //std::cout << "PairProductCache Size: " << pairProductCache->Size() << std::endl;
        delete pairProductCache<T,Op>;
        pairProductCache<T,Op> = NULL;
        delete pairProduct2Cache<T,Op>;
        pairProduct2Cache<T,Op> = NULL;
    }


    typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
	typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
    // typedef double BIG_FLOAT;
    template class WeightedBDDPairProductMapHandle<BIG_FLOAT>;
    template WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                BDDPairProduct<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n1,
                                WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<BIG_FLOAT> &pairProductMap,
                                BIG_FLOAT(*func)(BIG_FLOAT, BIG_FLOAT)
                                );
    template WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                BDDPairProduct<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n1,
                                WeightedBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<BIG_FLOAT> &pairProductMap,
                                BIG_FLOAT(*func)(BIG_FLOAT, BIG_FLOAT)
                            );
    template WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                BDDPairProduct2<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n1,
                                WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n2,
                                unsigned int numVars,
                                BIG_FLOAT factor1,
                                BIG_FLOAT factor2,
                                WeightedBDDPairProductMapHandle<BIG_FLOAT> &pairProductMap,
                                BIG_FLOAT(*func)(BIG_FLOAT, BIG_FLOAT)
                                );
    template WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                BDDPairProduct2<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n1,
                                WeightedBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n2,
                                unsigned int numVars,
                                BIG_FLOAT factor1,
                                BIG_FLOAT factor2,
                                WeightedBDDPairProductMapHandle<BIG_FLOAT> &pairProductMap,
                                BIG_FLOAT(*func)(BIG_FLOAT, BIG_FLOAT)
                            );

    template void InitWeightedBDDPairProductCache<BIG_FLOAT, std::multiplies<BIG_FLOAT>>();
    template void DisposeOfWeightedBDDPairProductCache<BIG_FLOAT, std::multiplies<BIG_FLOAT>>();

    template class WeightedBDDPairProductMapHandle<BIG_COMPLEX_FLOAT>;
    template WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                BDDPairProduct<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n1,
                                WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap,
                                BIG_COMPLEX_FLOAT(*func)(BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT)
                                );
    template WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                BDDPairProduct<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n1,
                                WeightedBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap,
                                BIG_COMPLEX_FLOAT(*func)(BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT)
                            );
    template WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                BDDPairProduct2<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n1,
                                WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n2,
                                unsigned int numVars,
                                BIG_COMPLEX_FLOAT factor1,
                                BIG_COMPLEX_FLOAT factor2,
                                WeightedBDDPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap,
                                BIG_COMPLEX_FLOAT(*func)(BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT)
                                );
    template WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                BDDPairProduct2<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n1,
                                WeightedBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n2,
                                unsigned int numVars,
                                BIG_COMPLEX_FLOAT factor1,
                                BIG_COMPLEX_FLOAT factor2,
                                WeightedBDDPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap,
                                BIG_COMPLEX_FLOAT(*func)(BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT)
                            );

    template void InitWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    template void DisposeOfWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();

    template class WeightedBDDPairProductMapHandle<fourierSemiring>;
    template WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> 
                BDDPairProduct<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> n1,
                                WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<fourierSemiring> &pairProductMap,
                                fourierSemiring(*func)(fourierSemiring, fourierSemiring)
                                );
    template WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> 
                BDDPairProduct<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n1,
                                WeightedBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n2,
                                unsigned int numVars,
                                WeightedBDDPairProductMapHandle<fourierSemiring> &pairProductMap,
                                fourierSemiring(*func)(fourierSemiring, fourierSemiring)
                            );

    template WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> 
                BDDPairProduct2<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> n1,
                                WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> n2,
                                unsigned int numVars,
                                fourierSemiring factor1,
                                fourierSemiring factor2,
                                WeightedBDDPairProductMapHandle<fourierSemiring> &pairProductMap,
                                fourierSemiring(*func)(fourierSemiring, fourierSemiring)
                                );
    template WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> 
                BDDPairProduct2<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n1,
                                WeightedBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n2,
                                unsigned int numVars,
                                fourierSemiring factor1,
                                fourierSemiring factor2,
                                WeightedBDDPairProductMapHandle<fourierSemiring> &pairProductMap,
                                fourierSemiring(*func)(fourierSemiring, fourierSemiring)
                            );

    template void InitWeightedBDDPairProductCache<fourierSemiring, std::multiplies<fourierSemiring>>();
    template void DisposeOfWeightedBDDPairProductCache<fourierSemiring, std::multiplies<fourierSemiring>>();


}