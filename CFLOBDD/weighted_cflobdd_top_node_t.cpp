#include "weighted_cflobdd_top_node_t.h"
#include "weighted_values.h"
#include "ref_ptr.h"
#include "weighted_values_list.h"
#include "weighted_cross_product.h"
#include "weighted_cross_product_bdd.h"


namespace CFL_OBDD {

    // Initializations of static members ---------------------------------

    template <typename T, typename Op>
    unsigned int const WeightedCFLOBDDTopNodeT<T, Op>::maxLevel = CFLOBDDMaxlevel;

    template <typename T, typename Op>
    Hashset<WeightedCFLOBDDTopNodeT<T, Op>> *WeightedCFLOBDDTopNodeT<T, Op>::computedCache = new Hashset<WeightedCFLOBDDTopNodeT<T, Op>>(10000);

    // Constructors/Destructor -------------------------------------------

    template<typename T, typename Op>
    WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeT(WeightedCFLOBDDNode<T, Op> *n, ReturnMapHandle<T> &mapHandle, T factor)
    {
    #ifdef WeightedCFLOBDDTopNodeTDebug
        if (n->level >= 1) {
            WeightedCFLOBDDInternalNode *in = (WeightedCFLOBDDInternalNode *)n;
            // Check for inconsistencies between the entries in the return maps of in's BConnections and mapHandle.Size()
            unsigned int bound = mapHandle.Size();
            for (unsigned int i = 0; i < in->numBConnections; i++) {
                for (unsigned int j = 0; j < in->BConnection[i].returnMapHandle.mapContents->mapArray.size(); j++) {
                    if (in->BConnection[i].returnMapHandle.Lookup(j) >= bound) {
                        std::cout << "Inconsistent WeightedCFLOBDDTopNodeT construction" << std::endl;
                        std::cout << "numExits = " << n->numExits << std::endl;
                        std::cout << "level = " << n->level << std::endl;
                        std::cout << "bound = " << bound << std::endl;
                        std::cout << "index j = " << j << std::endl;
                        std::cout << "map at j = " << in->BConnection[i].returnMapHandle.Lookup(j) << std::endl;
                        std::cout << "map = " << in->BConnection[i].returnMapHandle << std::endl;
                        abort();
                    }
                }
            }
        }
    #endif

        rootConnection = WRootConnection<ReturnMapHandle<T>, T, Op>(n, mapHandle, factor);
        level = n->level;
    }

    template<typename T, typename Op>
    WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeT(WeightedCFLOBDDNodeHandleT<T, Op> &nodeHandle, ReturnMapHandle<T> &mapHandle, T factor)
    {
    #ifdef WeightedCFLOBDDTopNodeTDebug
        if (nodeHandle.handleContents->level >= 1) {
            WeightedCFLOBDDInternalNode *in = (WeightedCFLOBDDInternalNode *)nodeHandle.handleContents;
            // Check for inconsistencies between the entries in the return maps of in's BConnections and mapHandle.Size()
            unsigned int bound = mapHandle.Size();
            for (unsigned int i = 0; i < in->numBConnections; i++) {
                for (unsigned int j = 0; j < in->BConnection[i].returnMapHandle.mapContents->mapArray.size(); j++) {
                    if (in->BConnection[i].returnMapHandle.Lookup(j) >= bound) {
                        std::cout << "Inconsistent WeightedCFLOBDDTopNodeT construction" << std::endl;
                        std::cout << "numExits = " << nodeHandle.handleContents->numExits << std::endl;
                        std::cout << "level = " << nodeHandle.handleContents->level << std::endl;
                        std::cout << "bound = " << bound << std::endl;
                        std::cout << "index j = " << j << std::endl;
                        std::cout << "map at j = " << in->BConnection[i].returnMapHandle.Lookup(j) << std::endl;
                        std::cout << "map = " << in->BConnection[i].returnMapHandle << std::endl;
                        abort();
                    }
                }
            }
        }
    #endif

        rootConnection = WRootConnection<ReturnMapHandle<T>, T, Op>(nodeHandle, mapHandle, factor);
        level = nodeHandle.handleContents->level;
    }

    template<typename T, typename Op>
    WeightedCFLOBDDTopNodeT<T, Op>::~WeightedCFLOBDDTopNodeT()
    {
    }

    // Evaluate
    //    Return the value of the Boolean function under the given assignment
    template<typename T, typename Op>
    T WeightedCFLOBDDTopNodeT<T, Op>::Evaluate(SH_OBDD::Assignment &assignment)
    {
        SH_OBDD::AssignmentIterator ai(assignment);
        int i = rootConnection.entryPointHandle->handleContents->Traverse(ai);
        T ans = rootConnection.returnMapHandle.Lookup(i);
        return ans;
    }


    // EvaluateIteratively
    //    Return the value of the Boolean function under the given assignment
    // template <typename T, typename Op>
    // T WeightedCFLOBDDTopNodeT<T, Op>::EvaluateIteratively(SH_OBDD::Assignment &assignment)
    // {
    //     SH_OBDD::AssignmentIterator ai(assignment);

    //     unsigned int exitIndex = 0;
    //     ConsCell<TraverseState<WeightedCFLOBDDNode>> *S = NULL;
    //     TraverseState<WeightedCFLOBDDNode> ts;
    //     T ans;

    //     S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(rootConnection.entryPointHandle->handleContents, FirstVisit), S);
    //     while (S != NULL) {
    //         ts = S->Item();
    //         S = S->Next();
    //         if (ts.node->NodeKind() == CFLOBDD_DONTCARE) {
    //             ai.Next();
    //             exitIndex = 0;
    //         }
    //         else if (ts.node->NodeKind() == CFLOBDD_FORK) {
    //             bool val = ai.Current();
    //             ai.Next();
    //             exitIndex = (int)val;
    //         }
    //         else {  // Must be a WeightedCFLOBDDInternalNode
    //             WeightedCFLOBDDInternalNode *n = (WeightedCFLOBDDInternalNode *)ts.node;

    //             if (ts.visitState == FirstVisit) {
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n, SecondVisit), S);
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n->AConnection.entryPointHandle->handleContents, FirstVisit), S);
    //             }
    //             else if (ts.visitState == SecondVisit) {
    //                 int i = n->AConnection.returnMapHandle.Lookup(exitIndex);
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n, ThirdVisit, i), S);
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n->BConnection[i].entryPointHandle->handleContents, FirstVisit), S);
    //             }
    //             else {  // if (ts.visitState == ThirdVisit)
    //                 exitIndex = n->BConnection[ts.index].returnMapHandle.Lookup(exitIndex);
    //             }
    //         }
    //     }
    //     ans = rootConnection.returnMapHandle.Lookup(exitIndex);
    //     return ans;
    // }

    // PrintYield -----------------------------------------------------

    // // PrintYieldAux
    // template <typename T, typename Op>
    // void WeightedCFLOBDDTopNodeT<T, Op>::PrintYieldAux(std::ostream * out, List<ConsCell<TraverseState<WeightedCFLOBDDNode<T,Op>>> *> &L, ConsCell<TraverseState<WeightedCFLOBDDNode<T,Op>>> *S)
    // {
    //     unsigned int exitIndex = 0;
    //     TraverseState<WeightedCFLOBDDNode> ts;
    //     T ans;

    //     while (S != NULL) {
    //         ts = S->Item();
    //         S = S->Next();
    //         if (ts.node->NodeKind() == CFLOBDD_DONTCARE) {
    //             if (ts.visitState == FirstVisit) {
    //                 L.AddToFront(new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(ts.node, Restart), S));
    //                 exitIndex = 0;
    //             }
    //             else {     // ts.visitState == Restart
    //                 exitIndex = 0;
    //             }
    //         }
    //         else if (ts.node->NodeKind() == CFLOBDD_FORK) {
    //             if (ts.visitState == FirstVisit) {
    //                 L.AddToFront(new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(ts.node, Restart), S));
    //                 exitIndex = 0;
    //             }
    //             else {     // ts.visitState == Restart
    //                 exitIndex = 1;
    //             }
    //         }
    //         else {  // Must be a WeightedCFLOBDDInternalNode
    //             WeightedCFLOBDDInternalNode *n = (WeightedCFLOBDDInternalNode *)ts.node;

    //             if (ts.visitState == FirstVisit) {
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n, SecondVisit), S);
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n->AConnection.entryPointHandle->handleContents, FirstVisit), S);
    //             }
    //             else if (ts.visitState == SecondVisit) {
    //                 int i = n->AConnection.returnMapHandle.Lookup(exitIndex);
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n, ThirdVisit, i), S);
    //                 S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(n->BConnection[i].entryPointHandle->handleContents, FirstVisit), S);
    //             }
    //             else {  // if (ts.visitState == ThirdVisit)
    //                 exitIndex = n->BConnection[ts.index].returnMapHandle.Lookup(exitIndex);
    //             }
    //         }
    //     }
    //     ans = rootConnection.returnMapHandle.Lookup(exitIndex);
    //     if (out != NULL) *out << ans;
    // }

    // // PrintYield
    // //
    // // print the yield of the CFLOBDDTopNode (i.e., the leaves of 0's and 1's
    // // in "left-to-right order").
    // //
    // template <typename T, typename Op>
    // void WeightedCFLOBDDTopNodeT<T, Op>::PrintYield(std::ostream * out)
    // {
    //     ConsCell<TraverseState<WeightedCFLOBDDNode>> *S = NULL;   // Traversal stack
    //     List<ConsCell<TraverseState<WeightedCFLOBDDNode>> *> L;   // Snapshot stack

    //     S = new ConsCell<TraverseState<WeightedCFLOBDDNode>>(TraverseState(rootConnection.entryPointHandle->handleContents, FirstVisit), S);
    //     PrintYieldAux(out, L, S);
    //     while (!L.IsEmpty()) {
    //         S = L.RemoveFirst();
    //         PrintYieldAux(out, L, S);
    //     }
    // }

    // PrintYieldSemantic
    //
    // print the yield of the CFLOBDD (i.e., the leaves of values
    // in "left-to-right order").
    //
    // template<typename T, typename Op>
    // void WeightedCFLOBDDTopNodeT<T, Op>::PrintYieldSemantic(std::ostream & out)
    // {
    //     if (level >= 1 && level <= 4) {
    //         unsigned int size = 1 << level;
    //         SH_OBDD::Assignment a(size);
    //         T b;
    //         unsigned long int range = 1UL << size;
    //         for (unsigned long int i = 0UL; i < range; i++) {
    //             unsigned long int mask = 1UL;
    //             for (int j = size - 1; j >= 0; j--) {
    //                 a[j] = (i & mask);
    //                 mask = mask << 1;
    //             }
    //             b = EvaluateIteratively(a);
    //             out << b;
    //         }
    //         out << std::endl;
    //     }
    //     else {
    //         std::cerr << "Cannot test all assignments: level must be in [1 .. 4]" << std::endl;
    //     }
    // }

    // template<typename T, typename Op>
    // bool WeightedCFLOBDDTopNodeT<T, Op>::IsValid()
    // {
    //     return rootConnection.entryPointHandle->handleContents->IsValid();
    // }


    // Satisfaction Operations ------------------------------------




    template <typename T, typename Op>
    void WeightedCFLOBDDTopNodeT<T, Op>::DumpConnections(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visited, std::ostream & out /* = std::cout */)
    {
        rootConnection.entryPointHandle->handleContents->DumpConnections(visited, out);
        out << rootConnection << std::endl;
    }

    template <typename T, typename Op>
    void WeightedCFLOBDDTopNodeT<T, Op>::CountNodes(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes, unsigned int &nodeCount)
    {
        rootConnection.entryPointHandle->handleContents->CountNodes(visitedNodes, nodeCount);
    }

    template <typename T, typename Op>
    void WeightedCFLOBDDTopNodeT<T, Op>::CountPaths(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes)
    {
        rootConnection.entryPointHandle->handleContents->CountPaths(visitedNodes);
    }


    template <typename T, typename Op>
    void WeightedCFLOBDDTopNodeT<T, Op>::ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes)
    {
        rootConnection.entryPointHandle->handleContents->ComputeWeightOfPathsAsAmpsToExits(visitedNodes);
    }


    template <typename T, typename Op>
    void WeightedCFLOBDDTopNodeT<T, Op>::CountNodesAndEdges(Hashset<WeightedCFLOBDDNodeHandleT<T, Op>> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
        unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount, unsigned int& returnEdgesObjCount)
    {
        rootConnection.entryPointHandle->handleContents->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount, returnEdgesCount);
        edgeCount += rootConnection.returnMapHandle.Size();
    }


    template <typename T, typename Op>
    void WeightedCFLOBDDTopNodeT<T, Op>::DeallocateMemory()
    {
        // WeightedCFLOBDDTopNodeT<T, Op>::~WeightedCFLOBDDTopNodeT();
    }

    // Hash
    template <typename T, typename Op>
    unsigned int WeightedCFLOBDDTopNodeT<T, Op>::Hash(unsigned int modsize)
    {
        return rootConnection.Hash(modsize);
    }

    // Overloaded !=
    template <typename T, typename Op>
    bool WeightedCFLOBDDTopNodeT<T, Op>::operator!= (const WeightedCFLOBDDTopNodeT<T, Op> & C)
    {
        return rootConnection != C.rootConnection;
    }

    // Overloaded ==
    template <typename T, typename Op>
    bool WeightedCFLOBDDTopNodeT<T, Op>::operator== (const WeightedCFLOBDDTopNodeT<T, Op> & C)
    {
        return rootConnection == C.rootConnection;
    }

    // print
    template <typename T, typename Op>
    std::ostream& WeightedCFLOBDDTopNodeT<T, Op>::print(std::ostream & out) const
    {
        out << *(rootConnection.entryPointHandle) << std::endl;
        out << rootConnection.returnMapHandle << std::endl;
        out << "{ factor: " << rootConnection.factor << " }" << std::endl;
        return out;
    }

    template <typename T, typename Op>
    std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNodeT<T, Op> &d)
    {
        d.print(out);
        return(out);
    }

    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    // ApplyAndReduce -----------------------------------------------------
    ApplyAndReduce(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n1,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n2,
                BoolOp op
                )
    {
        // Perform 2-way cross product of n1 and n2
        // PairProductMapHandle MapHandle;
        // WeightedCFLOBDDNodeHandleT<T,Op> n = PairProduct(*(n1->rootConnection.entryPointHandle),
        //     *(n2->rootConnection.entryPointHandle),
        //     MapHandle);

        // // Create returnMapHandle from MapHandle: Fold the pairs in MapHandle by applying
        // // [n1->rootConnection.returnMapHandle, n2->rootConnection.returnMapHandle]
        // // (component-wise) to each pair.
        // ReturnMapHandle<T> returnMapHandle;
        // //PairProductMapBodyIterator MapIterator(*MapHandle.mapContents);
        // //MapIterator.Reset();
        // std::unordered_map<T, unsigned int> reduction_map;
        // ReductionMapHandle reductionMapHandle;
        // unsigned int iterator = 0;
        // //while (!MapIterator.AtEnd()) {
        // while (iterator < MapHandle.Size()){
        //     T c1, c2;
        //     int first, second;
        //     //first = MapIterator.Current().First();
        //     //second = MapIterator.Current().Second();
        //     first = MapHandle[iterator].First();
        //     second = MapHandle[iterator].Second();
        //     c1 = n1->rootConnection.returnMapHandle.Lookup(first);
        //     c2 = n2->rootConnection.returnMapHandle.Lookup(second);
        //     T val = op[c1][c2];
        //     if (reduction_map.find(val) == reduction_map.end()){
        //         returnMapHandle.AddToEnd(val);
        //         reduction_map.insert(std::make_pair(val, returnMapHandle.Size() - 1));
        //         reductionMapHandle.AddToEnd(returnMapHandle.Size() - 1);
        //     }
        //     else{
        //         reductionMapHandle.AddToEnd(reduction_map[val]);
        //     }
        //     //MapIterator.Next();
        //     iterator++;
        // }
        // returnMapHandle.Canonicalize();
        // reductionMapHandle.Canonicalize();

        // // Perform reduction on n, with respect to the common elements that returnMapHandle maps together
        // //ReductionMapHandle inducedReductionMapHandle;
        // //ReturnMapHandle<T> inducedReturnMap;
        // //returnMapHandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
        // //WeightedCFLOBDDNodeHandleT<T,Op>::InitReduceCache();
        // //WeightedCFLOBDDNodeHandleT<T,Op> reduced_n = n.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
        // WeightedCFLOBDDNodeHandleT<T,Op> reduced_n = n.Reduce(reductionMapHandle, returnMapHandle.Size());
        // //WeightedCFLOBDDNodeHandleT<T,Op>::DisposeOfReduceCache();

        // // Create and return CFLOBDDTopNode
        // //return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_n, inducedReturnMap));
        // return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_n, returnMapHandle));
        return n1;
    }


    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n1,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n2,
                T(*func)(T, T),
                bool flag // False - Plus, True - Times
                )
    {
        if (n1->rootConnection.entryPointHandle->handleContents->NodeKind() == W_BDD_TOPNODE)
        {
            WeightedBDDTopNode<T,Op>* n1h = (WeightedBDDTopNode<T,Op> *) n1->rootConnection.entryPointHandle->handleContents;
            WeightedBDDTopNode<T,Op>* n2h = (WeightedBDDTopNode<T,Op> *) n2->rootConnection.entryPointHandle->handleContents;
            WeightedBDDTopNode<T,Op>* g = new WeightedBDDTopNode<T,Op>(n1h->numberOfVars);
            if (flag)
            {
                WeightedBDDPairProductMapHandle<T> pairProductMapHandle;
                auto ans = BDDPairProduct<T,Op>(n1h->bddContents, 
                            n2h->bddContents, 
                            n1h->numberOfVars,
                            pairProductMapHandle,
                            func);
                g->bddContents = ans;
                std::vector<T> v_ret = pairProductMapHandle.mapContents->mapArray;
                ReturnMapHandle<T> ret_ans;
                for (auto i : v_ret)
                    ret_ans.AddToEnd(i);
                return new WeightedCFLOBDDTopNodeT<T,Op>(g, ret_ans, pairProductMapHandle.mapContents->factor);
            }
            else {
                WeightedBDDPairProductMapHandle<T> pairProductMapHandle;
                auto ans = BDDPairProduct2<T,Op>(n1h->bddContents, 
                            n2h->bddContents, 
                            n1h->numberOfVars, 
                            n1->rootConnection.factor,
                            n2->rootConnection.factor,
                            pairProductMapHandle,
                            func);
                g->bddContents = ans;
                std::vector<T> v_ret = pairProductMapHandle.mapContents->mapArray;
                ReturnMapHandle<T> ret_ans;
                for (auto i : v_ret)
                    ret_ans.AddToEnd(i);
                return new WeightedCFLOBDDTopNodeT<T,Op>(g, ret_ans, pairProductMapHandle.mapContents->factor); 
            }
        }
        // Perform 2-way cross product of n1 and n2
        WeightedPairProductMapHandle<T> MapHandle;
        WeightedCFLOBDDNodeHandleT<T,Op> n;
        if (flag)
            n = PairProduct<T,Op>(*(n1->rootConnection.entryPointHandle),
                *(n2->rootConnection.entryPointHandle),
                MapHandle);
        else
            n = PairProduct2<T,Op>(*(n1->rootConnection.entryPointHandle),
                *(n2->rootConnection.entryPointHandle),
                n1->rootConnection.factor,
                n2->rootConnection.factor,
                MapHandle);

        // Create returnMapHandle from MapHandle: Fold the pairs in MapHandle by applying
        // [n1->rootConnection.returnMapHandle, n2->rootConnection.returnMapHandle]
        // (component-wise) to each pair.
        
        ReturnMapHandle<T> returnMapHandle;
        
        boost::unordered_map<T, unsigned int> reductionMap;
        ReductionMapHandle reductionMapHandle;
        WeightedValuesListHandle<T> valList;
        unsigned int iterator = 0;
        while (iterator < MapHandle.Size()){
            T c1, c2;
            int first, second;
            first = MapHandle[iterator].first.First();
            second = MapHandle[iterator].first.Second();
            T val;
            if (first == -1 && second == -1)
                val = getAnnhilatorValue<T,Op>();
            else {
                c1 = n1->rootConnection.returnMapHandle.Lookup(first);
                c2 = n2->rootConnection.returnMapHandle.Lookup(second);
                if (flag) 
                    val = (*func)(c1, c2);
                else{
                    T v1, v2;
                    v1 = MapHandle[iterator].second.First();
                    v2 = MapHandle[iterator].second.Second();
                    val = (*func)(computeComposition<T,Op>(v1, c1), computeComposition<T,Op>(v2, c2));
                }
            }
            T value_to_check = val;
            if (!flag)
                value_to_check = (val == getAnnhilatorValue<T,Op>()) ? val : getIdentityValue<T,Op>();
            if (reductionMap.find(value_to_check) == reductionMap.end()){
                returnMapHandle.AddToEnd(value_to_check);
                reductionMap.insert(std::make_pair(value_to_check, returnMapHandle.Size() - 1)); 
                valList.AddToEnd(val);
                reductionMapHandle.AddToEnd(returnMapHandle.Size() - 1);
            }
            else{
                reductionMapHandle.AddToEnd(reductionMap[val]);
                valList.AddToEnd(val);
            }
            iterator++;
        }

            returnMapHandle.Canonicalize();
            reductionMapHandle.Canonicalize();
            valList.Canonicalize();
            std::pair<WeightedCFLOBDDNodeHandleT<T,Op>, T> reduced_n = n.Reduce(reductionMapHandle, returnMapHandle.Size(), valList, true);
            T factor = reduced_n.second;
            if (flag)
                factor = reduced_n.second * n1->rootConnection.factor * n2->rootConnection.factor;
            return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_n.first, returnMapHandle, factor));
    }


    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n1,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n2,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr n3,
                BoolOp3 op
                )
    {
        // Perform 3-way cross product of n1, n2, and n3
        // TripleProductMapHandle MapHandle;
        // WeightedCFLOBDDNodeHandleT<T,Op> n = TripleProduct(*(n1->rootConnection.entryPointHandle),
        //     *(n2->rootConnection.entryPointHandle),
        //     *(n3->rootConnection.entryPointHandle),
        //     MapHandle);

        // // Create returnMapHandle from MapHandle: Fold the pairs in MapHandle by applying
        // // [n1->rootConnection.returnMapHandle, n2->rootConnection.returnMapHandle, n3->rootConnection.returnMapHandle]
        // // (component-wise) to each triple.
        // ReturnMapHandle<T> returnMapHandle;
        // TripleProductMapBodyIterator MapIterator(*MapHandle.mapContents);
        // MapIterator.Reset();
        // while (!MapIterator.AtEnd()) {
        //     T c1, c2, c3;
        //     int first, second, third;
        //     first = MapIterator.Current().First();
        //     second = MapIterator.Current().Second();
        //     third = MapIterator.Current().Third();
        //     c1 = n1->rootConnection.returnMapHandle.Lookup(first);
        //     c2 = n2->rootConnection.returnMapHandle.Lookup(second);
        //     c3 = n3->rootConnection.returnMapHandle.Lookup(third);
        //     returnMapHandle.AddToEnd(op[c1][c2][c3]);
        //     MapIterator.Next();
        // }
        // returnMapHandle.Canonicalize();

        // // Perform reduction on n, with respect to the common elements that returnMapHandle maps together
        // ReductionMapHandle inducedReductionMapHandle;
        // ReturnMapHandle<T> inducedReturnMap;
        // returnMapHandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
        // //WeightedCFLOBDDNodeHandleT<T,Op>::InitReduceCache();
        // WeightedCFLOBDDNodeHandleT<T,Op> reduced_n = n.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
        // //WeightedCFLOBDDNodeHandleT<T,Op>::DisposeOfReduceCache();
        // // Create and return CFLOBDDTopNode
        // return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_n, inducedReturnMap));
        return n1;
    }

    // \f.\g.(f + g)
    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkPlusTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g
                )
    {
        return ApplyAndReduce<T,Op>(f, g, PlusFunc, false);
    }

    // Pointwise addition: \f.\g.(f + g)
    template<typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator+(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g)
    {
        try{
            return MkPlusTopNode<T,Op>(f, g);
        }
        catch (std::exception e){
            std::cout << e.what() << std::endl;
            std::cout << "Addition" << std::endl;
            abort();
        }
    }


    // \f.\g.(f ^ g)
    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkExorTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f,
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g
    )
    {
        return ApplyAndReduce<T,Op>(f, g, exclusiveOrOp, false);
    }

    // Pointwise addition: \f.\g.(f ^ g)
    template<typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator^(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g)
    {
        return MkExorTopNode<T>(f, g);
    }


    //// Left scalar-multiplication: \c:unsigned int.\g.(c * g)
    //template <typename T, typename Op>
    //typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>>
    //MkLeftScalarTimesTopNode(unsigned int c, typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> g)
    //{
    //	if (c == 1) return g;  // Nothing need be performed
    //
    //	WeightedCFLOBDDNodeHandleT<T,Op> eph;
    //	typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> ans;
    //	if (c == 0) {  // Special case
    //		eph = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[g->level];
    //	}
    //	else {
    //		eph = *(g->rootConnection.entryPointHandle);
    //	}
    //
    //	ReturnMapHandle<T> rmh = c * g->rootConnection.returnMapHandle;
    //
    //	// Perform reduction on eph, with respect to the common elements that rmh maps together
    //	ReductionMapHandle inducedReductionMapHandle;
    //	ReturnMapHandle<T> inducedReturnMap;
    //	rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
    //	//     WeightedCFLOBDDNodeHandleT<T,Op>::InitReduceCache();
    //	WeightedCFLOBDDNodeHandleT<T,Op> reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
    //	//     WeightedCFLOBDDNodeHandleT<T,Op>::DisposeOfReduceCache();
    //
    //	// Create and return CFLOBDDTopNode
    //	return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_eph, inducedReturnMap));
    //}
    //
    //// Left scalar-multiplication: \c:int.\g.(c * g)
    //template <typename T, typename Op>
    //typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>>
    //MkLeftScalarTimesTopNode(int c, typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> g)
    //{
    //	if (c == 1) return g;  // Nothing need be performed
    //
    //	WeightedCFLOBDDNodeHandleT<T,Op> eph;
    //	typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> ans;
    //	if (c == 0) {  // Special case
    //		eph = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[g->level];
    //	}
    //	else {
    //		eph = *(g->rootConnection.entryPointHandle);
    //	}
    //
    //	ReturnMapHandle<T> rmh = c * g->rootConnection.returnMapHandle;
    //
    //	// Perform reduction on eph, with respect to the common elements that rmh maps together
    //	ReductionMapHandle inducedReductionMapHandle;
    //	ReturnMapHandle<T> inducedReturnMap;
    //	rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
    //	//     WeightedCFLOBDDNodeHandleT<T,Op>::InitReduceCache();
    //	WeightedCFLOBDDNodeHandleT<T,Op> reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
    //	//     WeightedCFLOBDDNodeHandleT<T,Op>::DisposeOfReduceCache();
    //
    //	// Create and return CFLOBDDTopNode
    //	return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_eph, inducedReturnMap));
    //}

    //// Left scalar-multiplication: \c:int.\g.(c * g)
    //template <typename T, typename Op>
    //typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>>
    //MkLeftScalarTimesTopNode(double c, typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> g)
    //{
    //	if (c == 1) return g;  // Nothing need be performed
    //
    //	WeightedCFLOBDDNodeHandleT<T,Op> eph;
    //	typename ref_ptr<WeightedCFLOBDDTopNodeT<T, Op>> ans;
    //	if (c == 0) {  // Special case
    //		eph = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[g->level];
    //	}
    //	else {
    //		eph = *(g->rootConnection.entryPointHandle);
    //	}
    //
    //	ReturnMapHandle<T> rmh = c * g->rootConnection.returnMapHandle;
    //
    //	// Perform reduction on eph, with respect to the common elements that rmh maps together
    //	ReductionMapHandle inducedReductionMapHandle;
    //	ReturnMapHandle<T> inducedReturnMap;
    //	rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
    //	//     WeightedCFLOBDDNodeHandleT<T,Op>::InitReduceCache();
    //	WeightedCFLOBDDNodeHandleT<T,Op> reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
    //	//     WeightedCFLOBDDNodeHandleT<T,Op>::DisposeOfReduceCache();
    //
    //	// Create and return CFLOBDDTopNode
    //	return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_eph, inducedReturnMap));
    //}

    // Left scalar-multiplication: \c:int.\g.(c * g)
    template <typename T, typename T1, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkLeftScalarTimesTopNode(T1 c, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g)
    {
        if (c == getIdentityValue<T, Op>()) return g;  // Nothing need be performed

        WeightedCFLOBDDNodeHandleT<T,Op> eph;
        typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr ans;
        if (c == getAnnhilatorValue<T,Op>()) {  // Special case
            eph = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode_Ann[g->level];
            // ReturnMapHandle<T> rmh = c * g->rootConnection.returnMapHandle;
            // // // Perform reduction on eph, with respect to the common elements that rmh maps tssogether
            // ReductionMapHandle inducedReductionMapHandle;
            // ReturnMapHandle<T> inducedReturnMap;
            // rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
            // WeightedCFLOBDDNodeHandleT<T,Op> reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
            ReturnMapHandle<T> rmh;
            rmh.AddToEnd(c);
            rmh.Canonicalize();
            // Create and return CFLOBDDTopNode
            return(new WeightedCFLOBDDTopNodeT<T, Op>(eph, rmh, getIdentityValue<T,Op>()));
        }
        else {
            eph = *(g->rootConnection.entryPointHandle);
            ReturnMapHandle<T> rmh = g->rootConnection.returnMapHandle;
            T factor = c * g->rootConnection.factor;
            return(new WeightedCFLOBDDTopNodeT<T, Op>(eph, rmh, factor));
        }
    }



    // Left scalar-multiplication: \c:int.\g.(c * g)
    template<typename T, typename T1, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator*(T1 c, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g)
    {
        try{
            return MkLeftScalarTimesTopNode<T, T1, Op>(c, g);
        }
        catch (std::exception e){
            std::cout << e.what() << std::endl;
            std::cout << "Multiply " << c << std::endl;
            throw e;
        }
    }

    // Right scalar-multiplication: \f.\c:int.(f * c)
    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkRightScalarTimesTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, int c)
    {
        if (c == 1) return f;  // Nothing need be performed

        WeightedCFLOBDDNodeHandleT<T,Op> eph;
        typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr ans;
        if (c == 0) {  // Special case
            eph = WeightedCFLOBDDNodeHandleT<T,Op>::NoDistinctionNode[f->level];
            ReturnMapHandle<T> rmh = f->rootConnection.returnMapHandle * c;

            // Perform reduction on eph, with respect to the common elements that rmh maps together
            ReductionMapHandle inducedReductionMapHandle;
            ReturnMapHandle<T> inducedReturnMap;
            rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
            //WeightedCFLOBDDNodeHandleT<T,Op>::InitReduceCache();
            WeightedCFLOBDDNodeHandleT<T,Op> reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
            //WeightedCFLOBDDNodeHandleT<T,Op>::DisposeOfReduceCache();

            // Create and return CFLOBDDTopNode
            return(new WeightedCFLOBDDTopNodeT<T, Op>(reduced_eph, inducedReturnMap));
        }
        else {
            eph = *(f->rootConnection.entryPointHandle);
            ReturnMapHandle<T> rmh = f->rootConnection.returnMapHandle * c;
            return(new WeightedCFLOBDDTopNodeT<T, Op>(eph, rmh));
        }
    }

    // Right scalar-multiplication: \f.\c:int.(f * c)
    template<typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator*(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, int c)
    {
        return MkRightScalarTimesTopNode<T>(f, c);
    }

    // Pointwise multiplication: \f.\g.(f * g) ---------------------------------------
    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    MkTimesTopNode(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f,
                typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g
                )
    {
        return ApplyAndReduce<T, Op>(f, g, TimesFunc, true);
    }

    // Pointwise multiplication: \f.\g.(f * g) 
    template <typename T, typename Op>
    typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr
    operator*(typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr f, typename WeightedCFLOBDDTopNodeT<T, Op>::WeightedCFLOBDDTopNodeTRefPtr g)
    {
        return MkTimesTopNode<T, Op>(f, g);
    }


}