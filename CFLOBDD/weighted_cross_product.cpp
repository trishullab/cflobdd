
#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "weighted_cflobdd_node_t.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "inttriple.h"
#include "weighted_cross_product.h"
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
// WeightedPairProductMapBody
//***************************************************************
namespace CFL_OBDD {

    // Initializations of static members ---------------------------------
    template <typename T>
    Hashset<WeightedPairProductMapBody<T>> *WeightedPairProductMapBody<T>::canonicalWeightedPairProductMapBodySet = new Hashset<WeightedPairProductMapBody<T>>(HASHSET_NUM_BUCKETS);

    // Constructor
    template <typename T>
    WeightedPairProductMapBody<T>::WeightedPairProductMapBody()
    : refCount(0), isCanonical(false)
    {
    }

    template <typename T>
    void WeightedPairProductMapBody<T>::IncrRef()
    {
        refCount++;    // Warning: Saturation not checked
    }

    template <typename T>
    void WeightedPairProductMapBody<T>::DecrRef()
    {
        if (--refCount == 0) {    // Warning: Saturation not checked
            if (isCanonical) {
            WeightedPairProductMapBody<T>::canonicalWeightedPairProductMapBodySet->DeleteEq(this);
            }
            delete this;
        }
    }

    template <typename T>
    unsigned int WeightedPairProductMapBody<T>::Hash(unsigned int modsize)
    {
        unsigned int hvalue = 0;

        for (unsigned int i = 0; i < mapArray.size(); i++){
            hvalue = (997*hvalue + (unsigned int)97*mapArray[i].First() + (unsigned int)mapArray[i].Second()) % modsize;
        }
        boost::hash<T> boost_hash;
        for (unsigned int i = 0; i < mapArray.size(); i++){
            hvalue = (997*hvalue + 97*boost_hash(valueArray[i].First()) + boost_hash(valueArray[i].Second())) % modsize;
        }

        return hvalue;
    }

    template <typename T>
    void WeightedPairProductMapBody<T>::setHashCheck()
    {
        unsigned int hvalue = 0;

        for (auto &i : mapArray) {
            hvalue = (117 * (hvalue + 1) + (int)(97 * i.First()) + i.Second());
        }
        boost::hash<T> boost_hash;
        for (unsigned int i = 0; i < valueArray.size(); i++){
            hvalue = (117* (hvalue + 1) + 97*boost_hash(valueArray[i].First()) + boost_hash(valueArray[i].Second())) % 2019;
        }
        hashCheck = hvalue;
    }

    template <typename T>
    void WeightedPairProductMapBody<T>::AddToEnd(const intpair& y, const Pair_T<T,T>& v)
    {
        mapArray.push_back(y);
        valueArray.push_back(v);
    }

    template <typename T>
    void WeightedPairProductMapBody<T>::AddToEnd(const intpair& y)
    {
        mapArray.push_back(y);
        // TODO: change this
        valueArray.push_back(Pair_T<T,T>(1.0, 1.0));
    }

    template <>
    void WeightedPairProductMapBody<fourierSemiring>::AddToEnd(const intpair& y)
    {
        mapArray.push_back(y);
        fourierSemiring one(1, 1);
        valueArray.push_back(Pair_T<fourierSemiring,fourierSemiring>(one, one)); 
    }

    template <typename T>
    bool WeightedPairProductMapBody<T>::operator==(const WeightedPairProductMapBody<T> &o) const
    {
        if (mapArray.size() != o.mapArray.size())
            return false;

        for (unsigned int i = 0; i < mapArray.size(); i++){
            if (mapArray[i] != o.mapArray[i] || valueArray[i] != o.valueArray[i])
                return false;
        }
        return true;
    }
    
    template <typename T>
    std::pair<intpair, Pair_T<T,T>> WeightedPairProductMapBody<T>::operator[](unsigned int i){                       // Overloaded []
        return std::make_pair(mapArray[i],valueArray[i]);
    }

    template <typename T>
    unsigned int WeightedPairProductMapBody<T>::Size(){
        return (unsigned int)mapArray.size();
    }


    // namespace CFL_OBDD {
    template <typename T>
    std::ostream& operator<< (std::ostream & out, const WeightedPairProductMapBody<T> &r)
    {
    //out << (List<int>&)r;
        for (unsigned int i = 0; i < r.mapArray.size(); i++)
        {
            out << "(" << r.mapArray[i] << ", " << r.valueArray[i] << ") ";
        }
        return(out);
    }
    // }
    //***************************************************************
    // WeightedPairProductMapHandle
    //***************************************************************


    // Default constructor
    template <typename T>
    WeightedPairProductMapHandle<T>::WeightedPairProductMapHandle()
    :  mapContents(new WeightedPairProductMapBody<T>())
    {
        mapContents->IncrRef();
    }

    // Destructor
    template <typename T>
    WeightedPairProductMapHandle<T>::~WeightedPairProductMapHandle()
    {
        mapContents->DecrRef();
    }

    // Copy constructor
    template <typename T>
    WeightedPairProductMapHandle<T>::WeightedPairProductMapHandle(const WeightedPairProductMapHandle<T> &r)
    :  mapContents(r.mapContents)
    {
        mapContents->IncrRef();
    }

    // Overloaded assignment
    template <typename T>
    WeightedPairProductMapHandle<T>& WeightedPairProductMapHandle<T>::operator= (const WeightedPairProductMapHandle<T> &r)
    {
        if (this != &r)      // don't assign to self!
        {
            WeightedPairProductMapBody<T> *temp = mapContents;
            mapContents = r.mapContents;
            mapContents->IncrRef();
            temp->DecrRef();
        }
        return *this;        
    }

    // Overloaded !=
    template <typename T>
    bool WeightedPairProductMapHandle<T>::operator!=(const WeightedPairProductMapHandle<T> &r)
    {
        return (mapContents != r.mapContents);
    }

    // Overloaded ==
    template <typename T>
    bool WeightedPairProductMapHandle<T>::operator==(const WeightedPairProductMapHandle<T> &r)
    {
        return (mapContents == r.mapContents);
    }

    template <typename T>
    std::ostream& operator<< (std::ostream & out, const WeightedPairProductMapHandle<T> &r)
    {
        out << "[" << *r.mapContents << "]";
        return(out);
    }

    template <typename T>
    unsigned int WeightedPairProductMapHandle<T>::Hash(unsigned int modsize)
    {
        return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
    }

    template <typename T>
    unsigned int WeightedPairProductMapHandle<T>::Size()
    {
        return mapContents->Size();
    }

    template <typename T>
    std::pair<intpair, Pair_T<T,T>> WeightedPairProductMapHandle<T>::operator[](unsigned int i)
    {
        return std::make_pair(mapContents->mapArray[i], mapContents->valueArray[i]);
    }

    template <typename T>
    void WeightedPairProductMapHandle<T>::AddToEnd(const intpair& p, const Pair_T<T,T>& v)
    {
        assert(mapContents->refCount <= 1);
        mapContents->AddToEnd(p, v);
    }

    template <typename T>
    void WeightedPairProductMapHandle<T>::AddToEnd(const intpair& p)
    {
        assert(mapContents->refCount <= 1);
        mapContents->AddToEnd(p);
    }

    template <typename T>
    bool WeightedPairProductMapHandle<T>::Member(intpair& p)
    {
        for (auto& i : mapContents->mapArray){
            if (i == p)
                return true;
        }
        return false;
    }

    template <typename T>
    int WeightedPairProductMapHandle<T>::Lookup(intpair& p, Pair_T<T,T>& v)
    {
        for (unsigned int i = 0; i < mapContents->mapArray.size(); i++){
            if (mapContents->mapArray[i] == p && mapContents->valueArray[i] == v)
                return i;
        }
        return -1;
    }

    template <typename T>
    void WeightedPairProductMapHandle<T>::Canonicalize()
    {
        WeightedPairProductMapBody<T> *answerContents;
        unsigned int hash = WeightedPairProductMapBody<T>::canonicalWeightedPairProductMapBodySet->GetHash(mapContents);
        answerContents = WeightedPairProductMapBody<T>::canonicalWeightedPairProductMapBodySet->Lookup(mapContents, hash);
        if (answerContents == NULL) {
            WeightedPairProductMapBody<T>::canonicalWeightedPairProductMapBodySet->Insert(mapContents, hash);
            mapContents->isCanonical = true;
        }
        else {
            answerContents->IncrRef();
            mapContents->DecrRef();
            mapContents = answerContents;
        }
    }

    // Create map with reversed entries
    template <typename T>
    WeightedPairProductMapHandle<T> WeightedPairProductMapHandle<T>::Flip()
    {
        WeightedPairProductMapHandle<T> answer;
        for (int i = 0; i < mapContents->Size(); i++){
            intpair p(mapContents->mapArray[i].Second(), mapContents->mapArray[i].First());
            Pair_T<T,T> v(mapContents->valueArray[i].Second(), mapContents->valueArray[i].First());
            answer.AddToEnd(p, v);
        }
        return answer;
    }

    //***************************************************************
    // WeightedPairProductKey
    //***************************************************************

    // Constructor
    template <typename T, typename Op>
    WeightedPairProductKey<T,Op>::WeightedPairProductKey(WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle1, WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle2)
    :  nodeHandle1(nodeHandle1), nodeHandle2(nodeHandle2)
    {
        factor1 = getIdentityValue<T,Op>();
        factor2 = getIdentityValue<T,Op>();
    }

    template <typename T, typename Op>
    WeightedPairProductKey<T,Op>::WeightedPairProductKey(WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle1, WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle2, T factor1, T factor2)
    :  nodeHandle1(nodeHandle1), nodeHandle2(nodeHandle2), factor1(factor1), factor2(factor2)
    {
    }

    // Hash
    template <typename T, typename Op>
    unsigned int WeightedPairProductKey<T,Op>::Hash(unsigned int modsize)
    {
        unsigned int hvalue = 0;
        boost::hash<T> boost_hash;
        hvalue = (997 * nodeHandle1.Hash(modsize) + nodeHandle2.Hash(modsize)) % modsize;
        hvalue = (hvalue + (117 * boost_hash(factor1) % modsize + 17 * boost_hash(factor2) % modsize)) % modsize;
        return hvalue;
    }

    // print
    template <typename T, typename Op>
    std::ostream& WeightedPairProductKey<T,Op>::print(std::ostream & out) const
    {
        out << "(" << nodeHandle1 << ", " << nodeHandle2 << "), (" << factor1 << ", " << factor2 << ")";
        return out;
    }

    template <typename T, typename Op>
    std::ostream& operator<< (std::ostream & out, const WeightedPairProductKey<T,Op> &p)
    {
        p.print(out);
        return(out);
    }

    template <typename T, typename Op>
    WeightedPairProductKey<T,Op>& WeightedPairProductKey<T,Op>::operator= (const WeightedPairProductKey<T,Op>& i)
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
    bool WeightedPairProductKey<T,Op>::operator!=(const WeightedPairProductKey<T,Op>& p)
    {
        return (factor1 != p.factor1) || (factor2 != p.factor2) || (nodeHandle1 != p.nodeHandle1) || (nodeHandle2 != p.nodeHandle2);
    }

    // Overloaded ==
    template <typename T, typename Op>
    bool WeightedPairProductKey<T,Op>::operator==(const WeightedPairProductKey<T,Op>& p)
    {
        return (factor1 == p.factor1) && (factor2 == p.factor2) && (nodeHandle1 == p.nodeHandle1) && (nodeHandle2 == p.nodeHandle2);
    }

    //***************************************************************
    // WeightedPairProductMemo
    //***************************************************************

    // Default constructor
    template <typename T, typename Op>
    WeightedPairProductMemo<T,Op>::WeightedPairProductMemo()
    :  nodeHandle(WeightedCFLOBDDNodeHandleT<T,Op>()), pairProductMapHandle(WeightedPairProductMapHandle<T>())
    {
    }

    // Constructor
    template <typename T, typename Op>
    WeightedPairProductMemo<T,Op>::WeightedPairProductMemo(WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle, WeightedPairProductMapHandle<T> pairProductMapHandle)
    :  nodeHandle(nodeHandle), pairProductMapHandle(pairProductMapHandle)
    {
    }

    template <typename T, typename Op>
    std::ostream& operator<< (std::ostream & out, const WeightedPairProductMemo<T,Op> &p)
    {
        out << "(" << p.nodeHandle << ", " << p.pairProductMapHandle << ")";
        return(out);
    }

    template <typename T, typename Op>
    WeightedPairProductMemo<T,Op>& WeightedPairProductMemo<T,Op>::operator= (const WeightedPairProductMemo<T,Op>& i)
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
    bool WeightedPairProductMemo<T,Op>::operator!=(const WeightedPairProductMemo<T,Op>& p)
    {
        return (nodeHandle != p.nodeHandle) || (pairProductMapHandle != p.pairProductMapHandle);
    }

    // Overloaded ==
    template <typename T, typename Op>
    bool WeightedPairProductMemo<T,Op>::operator==(const WeightedPairProductMemo<T,Op>& p)
    {
        return (nodeHandle == p.nodeHandle) && (pairProductMapHandle == p.pairProductMapHandle);
    }

    // --------------------------------------------------------------------
    // PairProduct
    //
    // Returns a new CFLOBDDNodeHandle, and (in WeightedPairProductMapHandle) a descriptor of the
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
    static Hashtable<WeightedPairProductKey<T,Op>, WeightedPairProductMemo<T,Op>> *pairProductCache = NULL;
    template <typename T, typename Op>
    static Hashtable<WeightedPairProductKey<T,Op>, WeightedPairProductMemo<T,Op>> *pairProduct2Cache = NULL;

    // namespace CFL_OBDD {

    template <typename T, typename Op>
    WeightedCFLOBDDNodeHandleT<T,Op> PairProduct(WeightedCFLOBDDInternalNode<T,Op> *n1,
                                WeightedCFLOBDDInternalNode<T,Op> *n2,
                                WeightedPairProductMapHandle<T> &pairProductMapHandle
                                )
    {
        if (n1 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[n1->level].handleContents) {
            if (n2 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[n2->level].handleContents) {   // ND, ND
            pairProductMapHandle.AddToEnd(intpair(0,0));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n1);
            }
            else if (n2 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level].handleContents) {
            pairProductMapHandle.AddToEnd(intpair(-1,-1));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n2);  
            }
            else {                                                                        // ND, XX
            for (unsigned int kk = 0; kk < n2->numExits; kk++) {
                pairProductMapHandle.AddToEnd(intpair(0,kk));
            }
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n2);
            }
        }
        else if (n1 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level].handleContents) {
            pairProductMapHandle.AddToEnd(intpair(-1,-1));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n1);
        }
        else if (n2 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level].handleContents) {
            pairProductMapHandle.AddToEnd(intpair(-1,-1));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n2);
        }
        else {
            if (n2 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[n2->level].handleContents) {   // XX, ND
            for (unsigned int kk = 0; kk < n1->numExits; kk++) {
                pairProductMapHandle.AddToEnd(intpair(kk,0));
            }
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n1);
            }
            else {                                                                        // XX, XX
            WeightedPairProductMapHandle<T> AMap;
            WeightedCFLOBDDInternalNode<T,Op> *n;
            unsigned int j;
            unsigned int curExit;
            int b1, b2;
            
            n = new WeightedCFLOBDDInternalNode<T,Op>(n1->level);
            
            // Perform the cross product of the AConnections
            WeightedCFLOBDDNodeHandleT<T,Op> aHandle =
                        PairProduct<T,Op>(*(n1->AConnection.entryPointHandle),
                                    *(n2->AConnection.entryPointHandle),
                                    AMap
                                    );
            // Fill in n->AConnection.returnMapHandle
            // Correctness relies on AMap having no duplicates
            CFLOBDDReturnMapHandle aReturnHandle;
                for (unsigned int k = 0; k < AMap.Size(); k++) {
                    aReturnHandle.AddToEnd(k);
                }
                aReturnHandle.Canonicalize();
                n->AConnection = WConnection<T,Op>(aHandle, aReturnHandle);
            // Perform the appropriate cross products of the BConnections
                j = 0;
                curExit = 0;
                n->numBConnections = AMap.Size();
                n->BConnection = new WConnection<T,Op>[n->numBConnections];
                unsigned int Aiterator = 0;
                std::unordered_map<intpair, unsigned int, intpair::intpair_hash> pair_to_index;
                //while (!AMapIterator.AtEnd()) {
                while (Aiterator < AMap.Size()) {
                WeightedPairProductMapHandle<T> BMap;
                auto AVal = AMap[Aiterator];
                b1 = AVal.first.First();
                b2 = AVal.first.Second();
                WeightedCFLOBDDNodeHandleT<T,Op> bHandle;
                if (b1 == -1 && b2 == -1){
                    bHandle = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level-1];
                    BMap.AddToEnd(intpair(-1,-1));
                }
                else{
                    bHandle = PairProduct(*(n1->BConnection[b1].entryPointHandle),
                                        *(n2->BConnection[b2].entryPointHandle),
                                        BMap
                                    );
                }
                CFLOBDDReturnMapHandle bReturnHandle;
                unsigned int Biterator = 0;
                while (Biterator < BMap.Size()){
                        int c1, c2;
                        if (BMap[Biterator].first.First() == -1 && BMap[Biterator].first.Second()){
                            c1 = -1;
                            c2 = -1;
                        }
                        else {
                            c1 = n1->BConnection[b1].returnMapHandle.Lookup(BMap[Biterator].first.First());
                            c2 = n2->BConnection[b2].returnMapHandle.Lookup(BMap[Biterator].first.Second());
                        }
                        // Test whether the pair (c1,c2) occurs in WeightedPairProductMapHandle
                        intpair p = intpair(c1, c2);
                        auto it = pair_to_index.find(p);
                        if (it != pair_to_index.end()){
                            bReturnHandle.AddToEnd(it->second);
                        }
                        else {   // New pair found (i.e., new exit node found)
                            pairProductMapHandle.AddToEnd(p);
                            bReturnHandle.AddToEnd(curExit);
                            pair_to_index.emplace(p, curExit);
                            curExit++;
                        }
                        Biterator++;
                    }
                bReturnHandle.Canonicalize();
                n->BConnection[j] = WConnection<T,Op>(bHandle, bReturnHandle);
                Aiterator++;
                j++;
                }
                n->numExits = curExit;
        #ifdef PATH_COUNTING_ENABLED
                n->InstallPathCounts();
        #endif
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n);
            }
        }      
    }

    template <typename T, typename Op>
    WeightedCFLOBDDNodeHandleT<T,Op> PairProduct(WeightedCFLOBDDNodeHandleT<T,Op> n1,
                                WeightedCFLOBDDNodeHandleT<T,Op> n2,
                                WeightedPairProductMapHandle<T> &pairProductMapHandle
                                )
    {
        WeightedPairProductMemo<T,Op> cachedWeightedPairProductMemo;
        bool isCached = pairProductCache<T,Op>->Fetch(WeightedPairProductKey(n1,n2), cachedWeightedPairProductMemo);
        if (isCached) {
            pairProductMapHandle = cachedWeightedPairProductMemo.pairProductMapHandle;
            return cachedWeightedPairProductMemo.nodeHandle;
        }
        else if (pairProductCache<T,Op>->Fetch(WeightedPairProductKey(n2,n1), cachedWeightedPairProductMemo)) {
            pairProductMapHandle = cachedWeightedPairProductMemo.pairProductMapHandle.Flip();
            return cachedWeightedPairProductMemo.nodeHandle;
        }
        else {
            WeightedCFLOBDDNodeHandleT<T,Op> answer;
        
            if (n1.handleContents->NodeKind() == W_CFLOBDD_INTERNAL) {
            answer = PairProduct<T,Op>((WeightedCFLOBDDInternalNode<T,Op> *)n1.handleContents,
                                (WeightedCFLOBDDInternalNode<T,Op> *)n2.handleContents,
                                pairProductMapHandle
                                );
            }
            else if (n1.handleContents->NodeKind() == W_CFLOBDD_FORK) {
            if (n2.handleContents->NodeKind() == W_CFLOBDD_FORK) {                 // W_CFLOBDD_FORK, W_CFLOBDD_FORK
                WeightedCFLOBDDForkNode<T,Op>* c1 = (WeightedCFLOBDDForkNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDForkNode<T,Op>* c2 = (WeightedCFLOBDDForkNode<T,Op> *)(n2.handleContents);

                T lw = computeComposition<T,Op>(c1->lweight, c2->lweight);
                T rw = computeComposition<T,Op>(c1->rweight, c2->rweight);

                if (lw == getAnnhilatorValue<T,Op>() && rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                }
                else if (lw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1, -1));
                    pairProductMapHandle.AddToEnd(intpair(1, 1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                }
                else if (rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0));
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10;
                }
                else
                {
                    pairProductMapHandle.AddToEnd(intpair(0,0));
                    pairProductMapHandle.AddToEnd(intpair(1,1));
                    pairProductMapHandle.Canonicalize();   
                    answer = WeightedCFLOBDDNodeHandleT(new WeightedCFLOBDDForkNode<T,Op>(lw, rw));
                }
            }
            else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // W_CFLOBDD_FORK, CFLOBDD_DONTCARE
                WeightedCFLOBDDForkNode<T,Op>* c1 = (WeightedCFLOBDDForkNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDDontCareNode<T,Op>* c2 = (WeightedCFLOBDDDontCareNode<T,Op> *)(n2.handleContents);

                T lw = computeComposition<T,Op>(c1->lweight, c2->lweight);
                T rw = computeComposition<T,Op>(c1->rweight, c2->rweight);

                if (lw == getAnnhilatorValue<T,Op>() && rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                }
                else if (lw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1, -1));
                    pairProductMapHandle.AddToEnd(intpair(1, 0));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                }
                else if (rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0));
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10;
                }
                else
                {
                    pairProductMapHandle.AddToEnd(intpair(0,0));
                    pairProductMapHandle.AddToEnd(intpair(1,0));
                    pairProductMapHandle.Canonicalize();   
                    answer = WeightedCFLOBDDNodeHandleT(new WeightedCFLOBDDForkNode<T,Op>(lw, rw));
                }
            }
            }
            else { /* n1.handleContents->NodeKind() == CFLOBDD_DONTCARE */
            if (n2.handleContents->NodeKind() == W_CFLOBDD_FORK) {                 // CFLOBDD_DONTCARE, W_CFLOBDD_FORK
                WeightedCFLOBDDDontCareNode<T,Op>* c1 = (WeightedCFLOBDDDontCareNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDForkNode<T,Op>* c2 = (WeightedCFLOBDDForkNode<T,Op> *)(n2.handleContents);

                T lw = computeComposition<T,Op>(c1->lweight, c2->lweight);
                T rw = computeComposition<T,Op>(c1->rweight, c2->rweight);

                if (lw == getAnnhilatorValue<T,Op>() && rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                }
                else if (lw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1, -1));
                    pairProductMapHandle.AddToEnd(intpair(0, 1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                }
                else if (rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0));
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10;
                }
                else
                {
                    pairProductMapHandle.AddToEnd(intpair(0,0));
                    pairProductMapHandle.AddToEnd(intpair(0,1));
                    pairProductMapHandle.Canonicalize();   
                    answer = WeightedCFLOBDDNodeHandleT(new WeightedCFLOBDDForkNode<T,Op>(lw, rw));
                }
            }
            else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // CFLOBDD_DONTCARE, CFLOBDD_DONTCARE
                WeightedCFLOBDDDontCareNode<T,Op>* c1 = (WeightedCFLOBDDDontCareNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDDontCareNode<T,Op>* c2 = (WeightedCFLOBDDDontCareNode<T,Op> *)(n2.handleContents);

                T lw = computeComposition<T,Op>(c1->lweight, c2->lweight);
                T rw = computeComposition<T,Op>(c1->rweight, c2->rweight);

                if (lw == getAnnhilatorValue<T,Op>() && rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                }
                else if (lw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(-1, -1));
                    pairProductMapHandle.AddToEnd(intpair(1, 1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                }
                else if (rw == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0));
                    pairProductMapHandle.AddToEnd(intpair(-1,-1));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10;
                }
                else
                {
                    pairProductMapHandle.AddToEnd(intpair(0,0));
                    pairProductMapHandle.Canonicalize();   
                    answer = WeightedCFLOBDDNodeHandleT(new WeightedCFLOBDDDontCareNode<T,Op>(lw, rw));
                }
            }
            }
                pairProductCache<T,Op>->Insert(WeightedPairProductKey(n1, n2),
                WeightedPairProductMemo(answer, pairProductMapHandle));
            return answer;
        }    
    }


    template <typename T, typename Op>
    WeightedCFLOBDDNodeHandleT<T,Op> PairProduct2(WeightedCFLOBDDInternalNode<T,Op> *n1,
                                WeightedCFLOBDDInternalNode<T,Op> *n2,
                                T factor1,
                                T factor2,
                                WeightedPairProductMapHandle<T> &pairProductMapHandle
                                )
    {
        if (n1 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level] && n2 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n2->level])
        {
            pairProductMapHandle.AddToEnd(intpair(-1,-1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n1);
        }
        if (n1 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level].handleContents) {
            for (int i = 0; i < n2->numExits; i++)
                pairProductMapHandle.AddToEnd(intpair(0,i), Pair_T<T,T>(factor1, factor2));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n2);
        }
        else if (n2 == WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level].handleContents) {
            for (int i = 0; i < n1->numExits; i++)
                pairProductMapHandle.AddToEnd(intpair(i,0), Pair_T<T,T>(factor1, factor2));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n1);
        }
        else if (n1 == n2) {
            for (int i = 0; i < n1->numExits; i++)
                pairProductMapHandle.AddToEnd(intpair(i,i), Pair_T<T,T>(factor1, factor2));
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n1);
        }
        else {
            WeightedPairProductMapHandle<T> AMap;
            WeightedCFLOBDDInternalNode<T,Op> *n;
            unsigned int j;
            unsigned int curExit;
            int b1, b2;
            
            n = new WeightedCFLOBDDInternalNode<T,Op>(n1->level);
            
            // Perform the cross product of the AConnections
            WeightedCFLOBDDNodeHandleT<T,Op> aHandle =
                        PairProduct2<T,Op>(*(n1->AConnection.entryPointHandle),
                                    *(n2->AConnection.entryPointHandle),
                                    factor1,
                                    factor2,
                                    AMap
                                    );
            // Fill in n->AConnection.returnMapHandle
            // Correctness relies on AMap having no duplicates
            CFLOBDDReturnMapHandle aReturnHandle;
                for (unsigned int k = 0; k < AMap.Size(); k++) {
                    aReturnHandle.AddToEnd(k);
                }
                aReturnHandle.Canonicalize();
                n->AConnection = WConnection<T,Op>(aHandle, aReturnHandle);
            // Perform the appropriate cross products of the BConnections
                j = 0;
                curExit = 0;
                n->numBConnections = AMap.Size();
                n->BConnection = new WConnection<T,Op>[n->numBConnections];
                unsigned int Aiterator = 0;
                std::unordered_map<FourTuple<T>, unsigned int, typename FourTuple<T>::FourTuple_Hash> pair_to_index;
                //while (!AMapIterator.AtEnd()) {
                while (Aiterator < AMap.Size()) {
                WeightedPairProductMapHandle<T> BMap;
                auto AVal = AMap[Aiterator];
                b1 = AVal.first.First();
                b2 = AVal.first.Second();
                WeightedCFLOBDDNodeHandleT<T,Op> bHandle;
                if (b1 == -1 && b2 == -1){
                    bHandle = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[n1->level-1];
                    BMap.AddToEnd(intpair(-1,-1));
                }
                else{
                    T v1 = AVal.second.First();
                    T v2 = AVal.second.Second();
                    bHandle = PairProduct2(*(n1->BConnection[b1].entryPointHandle),
                                        *(n2->BConnection[b2].entryPointHandle),
                                        v1,
                                        v2,
                                        BMap
                                    );
                }
                CFLOBDDReturnMapHandle bReturnHandle;
                unsigned int Biterator = 0;
                while (Biterator < BMap.Size()){
                        int c1, c2;
                        T v1 = getAnnhilatorValue<T,Op>(), v2 = getAnnhilatorValue<T,Op>();
                        if (BMap[Biterator].first.First() == -1 && BMap[Biterator].first.Second()){
                            c1 = -1;
                            c2 = -1;
                        }
                        else {
                            c1 = n1->BConnection[b1].returnMapHandle.Lookup(BMap[Biterator].first.First());
                            c2 = n2->BConnection[b2].returnMapHandle.Lookup(BMap[Biterator].first.Second());
                            v1 = BMap[Biterator].second.First();
                            v2 = BMap[Biterator].second.Second();
                        }
                        // Test whether the pair (c1,c2) occurs in WeightedPairProductMapHandle
                        auto p = FourTuple(c1, c2, v1, v2);
                        auto it = pair_to_index.find(p);
                        if (it != pair_to_index.end()){
                            bReturnHandle.AddToEnd(it->second);
                        }
                        else {   // New pair found (i.e., new exit node found)
                            pairProductMapHandle.AddToEnd(intpair(c1, c2), Pair_T<T,T>(v1,v2));
                            bReturnHandle.AddToEnd(curExit);
                            pair_to_index.emplace(p, curExit);
                            curExit++;
                        }
                        Biterator++;
                    }
                bReturnHandle.Canonicalize();
                n->BConnection[j] = WConnection<T,Op>(bHandle, bReturnHandle);
                Aiterator++;
                j++;
                }
                n->numExits = curExit;
        #ifdef PATH_COUNTING_ENABLED
                n->InstallPathCounts();
        #endif
            pairProductMapHandle.Canonicalize();
            return WeightedCFLOBDDNodeHandleT(n);
        }      
    }


    template <typename T, typename Op>
    WeightedCFLOBDDNodeHandleT<T,Op> PairProduct2(WeightedCFLOBDDNodeHandleT<T,Op> n1,
                                WeightedCFLOBDDNodeHandleT<T,Op> n2,
                                T factor1,
                                T factor2,
                                WeightedPairProductMapHandle<T> &pairProductMapHandle
                                )
    {
        WeightedPairProductMemo<T,Op> cachedWeightedPairProductMemo;
        bool isCached = pairProduct2Cache<T,Op>->Fetch(WeightedPairProductKey(n1,n2,factor1,factor2), cachedWeightedPairProductMemo);
        if (isCached) {
            pairProductMapHandle = cachedWeightedPairProductMemo.pairProductMapHandle;
            return cachedWeightedPairProductMemo.nodeHandle;
        }
        else if (pairProduct2Cache<T,Op>->Fetch(WeightedPairProductKey(n2,n1,factor2,factor1), cachedWeightedPairProductMemo)) {
            pairProductMapHandle = cachedWeightedPairProductMemo.pairProductMapHandle.Flip();
            return cachedWeightedPairProductMemo.nodeHandle;
        }
        else {
            WeightedCFLOBDDNodeHandleT<T,Op> answer;
        
            if (n1.handleContents->NodeKind() == W_CFLOBDD_INTERNAL) {
            answer = PairProduct2<T,Op>((WeightedCFLOBDDInternalNode<T,Op> *)n1.handleContents,
                                (WeightedCFLOBDDInternalNode<T,Op> *)n2.handleContents,
                                factor1,
                                factor2,
                                pairProductMapHandle
                                );
            }
            else if (n1.handleContents->NodeKind() == W_CFLOBDD_FORK) {
            if (n2.handleContents->NodeKind() == W_CFLOBDD_FORK) {                 // W_CFLOBDD_FORK, W_CFLOBDD_FORK
                WeightedCFLOBDDForkNode<T,Op>* c1 = (WeightedCFLOBDDForkNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDForkNode<T,Op>* c2 = (WeightedCFLOBDDForkNode<T,Op> *)(n2.handleContents);
                T lw1 = computeComposition<T,Op>(factor1, c1->lweight);
                T lw2 = computeComposition<T,Op>(factor2, c2->lweight);
                T rw1 = computeComposition<T,Op>(factor1, c1->rweight);
                T rw2 = computeComposition<T,Op>(factor2, c2->rweight);

                if (lw1 == getAnnhilatorValue<T,Op>() && lw2 == getAnnhilatorValue<T,Op>())
                {
                    if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                    {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                    }
                    else {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.AddToEnd(intpair(1, 1), Pair_T<T,T>(rw1, rw2));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                    }
                }
                else if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10; 
                }
                else {
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(1, 1), Pair_T<T,T>(rw1, rw2));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle;  
                }
            }
            else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // W_CFLOBDD_FORK, CFLOBDD_DONTCARE
                WeightedCFLOBDDForkNode<T,Op>* c1 = (WeightedCFLOBDDForkNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDForkNode<T,Op>* c2 = (WeightedCFLOBDDForkNode<T,Op> *)(n2.handleContents);
                T lw1 = computeComposition<T,Op>(factor1, c1->lweight);
                T lw2 = computeComposition<T,Op>(factor2, c2->lweight);
                T rw1 = computeComposition<T,Op>(factor1, c1->rweight);
                T rw2 = computeComposition<T,Op>(factor2, c2->rweight);

                if (lw1 == getAnnhilatorValue<T,Op>() && lw2 == getAnnhilatorValue<T,Op>())
                {
                    if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                    {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                    }
                    else {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.AddToEnd(intpair(1, 0), Pair_T<T,T>(rw1, rw2));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                    }
                }
                else if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10; 
                }
                else{
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(1, 0), Pair_T<T,T>(rw1, rw2));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle; 
                }
            }
            }
            else { /* n1.handleContents->NodeKind() == CFLOBDD_DONTCARE */
            if (n2.handleContents->NodeKind() == W_CFLOBDD_FORK) {                 // CFLOBDD_DONTCARE, W_CFLOBDD_FORK
                WeightedCFLOBDDForkNode<T,Op>* c1 = (WeightedCFLOBDDForkNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDForkNode<T,Op>* c2 = (WeightedCFLOBDDForkNode<T,Op> *)(n2.handleContents);
                T lw1 = computeComposition<T,Op>(factor1, c1->lweight);
                T lw2 = computeComposition<T,Op>(factor2, c2->lweight);
                T rw1 = computeComposition<T,Op>(factor1, c1->rweight);
                T rw2 = computeComposition<T,Op>(factor2, c2->rweight);

                if (lw1 == getAnnhilatorValue<T,Op>() && lw2 == getAnnhilatorValue<T,Op>())
                {
                    if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                    {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                    }
                    else {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.AddToEnd(intpair(0, 1), Pair_T<T,T>(rw1, rw2));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                    }
                }
                else if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10; 
                }
                else {
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(0, 1), Pair_T<T,T>(rw1, rw2));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle; 
                }
            }
            else { /* n2.handleContents->NodeKind() == CFLOBDD_DONTCARE */       // CFLOBDD_DONTCARE, CFLOBDD_DONTCARE
                // TODO: To optimize
                WeightedCFLOBDDForkNode<T,Op>* c1 = (WeightedCFLOBDDForkNode<T,Op> *)(n1.handleContents);
                WeightedCFLOBDDForkNode<T,Op>* c2 = (WeightedCFLOBDDForkNode<T,Op> *)(n2.handleContents);
                T lw1 = computeComposition<T,Op>(factor1, c1->lweight);
                T lw2 = computeComposition<T,Op>(factor2, c2->lweight);
                T rw1 = computeComposition<T,Op>(factor1, c1->rweight);
                T rw2 = computeComposition<T,Op>(factor2, c2->rweight);

                if (lw1 == getAnnhilatorValue<T,Op>() && lw2 == getAnnhilatorValue<T,Op>())
                {
                    if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                    {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[0];
                    }
                    else {
                        pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                        pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(rw1, rw2));
                        pairProductMapHandle.Canonicalize();
                        answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle01;
                    }
                }
                else if (rw1 == getAnnhilatorValue<T,Op>() && rw2 == getAnnhilatorValue<T,Op>())
                {
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(-1, -1), Pair_T<T,T>(getAnnhilatorValue<T,Op>(), getAnnhilatorValue<T,Op>()));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle10; 
                }
                else {
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(lw1, lw2));
                    pairProductMapHandle.AddToEnd(intpair(0, 0), Pair_T<T,T>(rw1, rw2));
                    pairProductMapHandle.Canonicalize();
                    answer = WeightedCFLOBDDNodeHandleT<T,Op>::CFLOBDDForkNodeHandle; 
                }
            }
            }
                pairProduct2Cache<T,Op>->Insert(WeightedPairProductKey(n1, n2, factor1, factor2),
                WeightedPairProductMemo(answer, pairProductMapHandle));
            return answer;
        }    
    }

    template <typename T, typename Op>
    void InitWeightedPairProductCache()
    {
        pairProductCache<T,Op> = new Hashtable<WeightedPairProductKey<T,Op>, WeightedPairProductMemo<T,Op>>(HASH_NUM_BUCKETS);
        pairProduct2Cache<T,Op> = new Hashtable<WeightedPairProductKey<T,Op>, WeightedPairProductMemo<T,Op>>(HASH_NUM_BUCKETS);
    }

    template <typename T, typename Op>
    void DisposeOfWeightedPairProductCache()
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
    template class WeightedPairProductMapHandle<BIG_FLOAT>;
    template WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                PairProduct<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n1,
                                WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n2,
                                WeightedPairProductMapHandle<BIG_FLOAT> &pairProductMap
                                );
    template WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                PairProduct(WeightedCFLOBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n1,
                                WeightedCFLOBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n2,
                                WeightedPairProductMapHandle<BIG_FLOAT> &pairProductMap
                            );
    template WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                PairProduct2<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n1,
                                WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> n2,
                                BIG_FLOAT factor1,
                                BIG_FLOAT factor2,
                                WeightedPairProductMapHandle<BIG_FLOAT> &pairProductMap
                                );
    template WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> 
                PairProduct2(WeightedCFLOBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n1,
                                WeightedCFLOBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>> *n2,
                                BIG_FLOAT factor1,
                                BIG_FLOAT factor2,
                                WeightedPairProductMapHandle<BIG_FLOAT> &pairProductMap
                            );

    template void InitWeightedPairProductCache<BIG_FLOAT, std::multiplies<BIG_FLOAT>>();
    template void DisposeOfWeightedPairProductCache<BIG_FLOAT, std::multiplies<BIG_FLOAT>>();

    template class WeightedPairProductMapHandle<BIG_COMPLEX_FLOAT>;
    template WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                PairProduct<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n1,
                                WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n2,
                                WeightedPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap
                                );
    template WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                PairProduct(WeightedCFLOBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n1,
                                WeightedCFLOBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n2,
                                WeightedPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap
                            );
    template WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                PairProduct2<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n1,
                                WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> n2,
                                BIG_COMPLEX_FLOAT factor1,
                                BIG_COMPLEX_FLOAT factor2,
                                WeightedPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap
                                );
    template WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> 
                PairProduct2(WeightedCFLOBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n1,
                                WeightedCFLOBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> *n2,
                                BIG_COMPLEX_FLOAT factor1,
                                BIG_COMPLEX_FLOAT factor2,
                                WeightedPairProductMapHandle<BIG_COMPLEX_FLOAT> &pairProductMap
                            );

    template void InitWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    template void DisposeOfWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();

    template class WeightedPairProductMapHandle<fourierSemiring>;
    template WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> 
                PairProduct<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> n1,
                                WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> n2,
                                WeightedPairProductMapHandle<fourierSemiring> &pairProductMap
                                );
    template WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> 
                PairProduct(WeightedCFLOBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n1,
                                WeightedCFLOBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n2,
                                WeightedPairProductMapHandle<fourierSemiring> &pairProductMap
                            );

    template WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> 
                PairProduct2<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> n1,
                                WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> n2,
                                fourierSemiring factor1,
                                fourierSemiring factor2,
                                WeightedPairProductMapHandle<fourierSemiring> &pairProductMap
                                );
    template WeightedCFLOBDDNodeHandleT<fourierSemiring, std::multiplies<fourierSemiring>> 
                PairProduct2(WeightedCFLOBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n1,
                                WeightedCFLOBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>> *n2,
                                fourierSemiring factor1,
                                fourierSemiring factor2,
                                WeightedPairProductMapHandle<fourierSemiring> &pairProductMap
                            );

    template void InitWeightedPairProductCache<fourierSemiring, std::multiplies<fourierSemiring>>();
    template void DisposeOfWeightedPairProductCache<fourierSemiring, std::multiplies<fourierSemiring>>();


}