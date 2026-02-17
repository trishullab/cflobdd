#include "cflobdd_top_node_t.h"
#include "ref_ptr.h"


namespace CFL_OBDD{

    // Initializations of static members ---------------------------------

    template <typename T>
    unsigned int const CFLOBDDTopNodeT<T>::maxLevel = CFLOBDDMaxLevel;

    template <typename T>
    Hashset<CFLOBDDTopNodeT<T>> *CFLOBDDTopNodeT<T>::computedCache = new Hashset<CFLOBDDTopNodeT<T>>(10000);

    // Constructors/Destructor -------------------------------------------

    template<typename T>
    CFLOBDDTopNodeT<T>::CFLOBDDTopNodeT(CFLOBDDNode *n, ReturnMapHandle<T> &mapHandle)
    {
    #ifdef CFLOBDDTopNodeTDebug
        if (n->level >= 1) {
            CFLOBDDInternalNode *in = (CFLOBDDInternalNode *)n;
            // Check for inconsistencies between the entries in the return maps of in's BConnections and mapHandle.Size()
            unsigned int bound = mapHandle.Size();
            for (unsigned int i = 0; i < in->numBConnections; i++) {
                for (unsigned int j = 0; j < in->BConnection[i].returnMapHandle.mapContents->mapArray.size(); j++) {
                    if (in->BConnection[i].returnMapHandle.Lookup(j) >= bound) {
                        std::cout << "Inconsistent CFLOBDDTopNodeT construction" << std::endl;
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

        rootConnection = ConnectionT<ReturnMapHandle<T>>(n, mapHandle);
        level = n->level;
    }

    template<typename T>
    CFLOBDDTopNodeT<T>::CFLOBDDTopNodeT(CFLOBDDNodeHandle &nodeHandle, ReturnMapHandle<T> &mapHandle)
    {
    #ifdef CFLOBDDTopNodeTDebug
        if (nodeHandle.handleContents->level >= 1) {
            CFLOBDDInternalNode *in = (CFLOBDDInternalNode *)nodeHandle.handleContents;
            // Check for inconsistencies between the entries in the return maps of in's BConnections and mapHandle.Size()
            unsigned int bound = mapHandle.Size();
            for (unsigned int i = 0; i < in->numBConnections; i++) {
                for (unsigned int j = 0; j < in->BConnection[i].returnMapHandle.mapContents->mapArray.size(); j++) {
                    if (in->BConnection[i].returnMapHandle.Lookup(j) >= bound) {
                        std::cout << "Inconsistent CFLOBDDTopNodeT construction" << std::endl;
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

        rootConnection = ConnectionT<ReturnMapHandle<T>>(nodeHandle, mapHandle);
        level = nodeHandle.handleContents->level;
    }

    template<typename T>
    CFLOBDDTopNodeT<T>::~CFLOBDDTopNodeT()
    {
    }

    // Evaluate
    //    Return the value of the Boolean function under the given assignment
    template<typename T>
    T CFLOBDDTopNodeT<T>::Evaluate(SH_OBDD::Assignment &assignment)
    {
        SH_OBDD::AssignmentIterator ai(assignment);
        int i = rootConnection.entryPointHandle->handleContents->Traverse(ai);
        T ans = rootConnection.returnMapHandle.Lookup(i);
        return ans;
    }


    // EvaluateIteratively
    //    Return the value of the Boolean function under the given assignment
    template <typename T>
    T CFLOBDDTopNodeT<T>::EvaluateIteratively(SH_OBDD::Assignment &assignment)
    {
        SH_OBDD::AssignmentIterator ai(assignment);

        unsigned int exitIndex = 0;
        ConsCell<TraverseState> *S = NULL;
        TraverseState ts;
        T ans;

        S = new ConsCell<TraverseState>(TraverseState(rootConnection.entryPointHandle->handleContents, FirstVisit), S);
        while (S != NULL) {
            ts = S->Item();
            S = S->Next();
            if (ts.node->NodeKind() == CFLOBDD_DONTCARE) {
                ai.Next();
                exitIndex = 0;
            }
            else if (ts.node->NodeKind() == CFLOBDD_FORK) {
                bool val = ai.Current();
                ai.Next();
                exitIndex = (int)val;
            }
            else {  // Must be a CFLOBDDInternalNode
                CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)ts.node;

                if (ts.visitState == FirstVisit) {
                    S = new ConsCell<TraverseState>(TraverseState(n, SecondVisit), S);
                    S = new ConsCell<TraverseState>(TraverseState(n->AConnection.entryPointHandle->handleContents, FirstVisit), S);
                }
                else if (ts.visitState == SecondVisit) {
                    int i = n->AConnection.returnMapHandle.Lookup(exitIndex);
                    S = new ConsCell<TraverseState>(TraverseState(n, ThirdVisit, i), S);
                    S = new ConsCell<TraverseState>(TraverseState(n->BConnection[i].entryPointHandle->handleContents, FirstVisit), S);
                }
                else {  // if (ts.visitState == ThirdVisit)
                    exitIndex = n->BConnection[ts.index].returnMapHandle.Lookup(exitIndex);
                }
            }
        }
        ans = rootConnection.returnMapHandle.Lookup(exitIndex);
        return ans;
    }

    // PrintYield -----------------------------------------------------

    // PrintYieldAux
    template <typename T>
    void CFLOBDDTopNodeT<T>::PrintYieldAux(std::ostream * out, List<ConsCell<TraverseState> *> &L, ConsCell<TraverseState> *S)
    {
        unsigned int exitIndex = 0;
        TraverseState ts;
        T ans;

        while (S != NULL) {
            ts = S->Item();
            S = S->Next();
            if (ts.node->NodeKind() == CFLOBDD_DONTCARE) {
                if (ts.visitState == FirstVisit) {
                    L.AddToFront(new ConsCell<TraverseState>(TraverseState(ts.node, Restart), S));
                    exitIndex = 0;
                }
                else {     // ts.visitState == Restart
                    exitIndex = 0;
                }
            }
            else if (ts.node->NodeKind() == CFLOBDD_FORK) {
                if (ts.visitState == FirstVisit) {
                    L.AddToFront(new ConsCell<TraverseState>(TraverseState(ts.node, Restart), S));
                    exitIndex = 0;
                }
                else {     // ts.visitState == Restart
                    exitIndex = 1;
                }
            }
            else {  // Must be a CFLOBDDInternalNode
                CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)ts.node;

                if (ts.visitState == FirstVisit) {
                    S = new ConsCell<TraverseState>(TraverseState(n, SecondVisit), S);
                    S = new ConsCell<TraverseState>(TraverseState(n->AConnection.entryPointHandle->handleContents, FirstVisit), S);
                }
                else if (ts.visitState == SecondVisit) {
                    int i = n->AConnection.returnMapHandle.Lookup(exitIndex);
                    S = new ConsCell<TraverseState>(TraverseState(n, ThirdVisit, i), S);
                    S = new ConsCell<TraverseState>(TraverseState(n->BConnection[i].entryPointHandle->handleContents, FirstVisit), S);
                }
                else {  // if (ts.visitState == ThirdVisit)
                    exitIndex = n->BConnection[ts.index].returnMapHandle.Lookup(exitIndex);
                }
            }
        }
        ans = rootConnection.returnMapHandle.Lookup(exitIndex);
        if (out != NULL) *out << ans;
    }

    // PrintYield
    //
    // print the yield of the CFLOBDDTopNode (i.e., the leaves of 0's and 1's
    // in "left-to-right order").
    //
    template <typename T>
    void CFLOBDDTopNodeT<T>::PrintYield(std::ostream * out)
    {
        ConsCell<TraverseState> *S = NULL;   // Traversal stack
        List<ConsCell<TraverseState> *> L;   // Snapshot stack

        S = new ConsCell<TraverseState>(TraverseState(rootConnection.entryPointHandle->handleContents, FirstVisit), S);
        PrintYieldAux(out, L, S);
        while (!L.IsEmpty()) {
            S = L.RemoveFirst();
            PrintYieldAux(out, L, S);
        }
    }

    // PrintYieldSemantic
    //
    // print the yield of the CFLOBDD (i.e., the leaves of values
    // in "left-to-right order").
    //
    template<typename T>
    void CFLOBDDTopNodeT<T>::PrintYieldSemantic(std::ostream & out)
    {
        if (level >= 1 && level <= 4) {
            unsigned int size = 1 << level;
            SH_OBDD::Assignment a(size);
            T b;
            unsigned long int range = 1UL << size;
            for (unsigned long int i = 0UL; i < range; i++) {
                unsigned long int mask = 1UL;
                for (int j = size - 1; j >= 0; j--) {
                    a[j] = (i & mask);
                    mask = mask << 1;
                }
                b = EvaluateIteratively(a);
                out << b;
            }
            out << std::endl;
        }
        else {
            std::cerr << "Cannot test all assignments: level must be in [1 .. 4]" << std::endl;
        }
    }

    template<typename T>
    bool CFLOBDDTopNodeT<T>::IsValid()
    {
        return rootConnection.entryPointHandle->handleContents->IsValid();
    }


    // Satisfaction Operations ------------------------------------




    template <typename T>
    void CFLOBDDTopNodeT<T>::DumpConnections(Hashset<CFLOBDDNodeHandle> *visited, std::ostream & out /* = std::cout */)
    {
        rootConnection.entryPointHandle->handleContents->DumpConnections(visited, out);
        out << rootConnection << std::endl;
    }

    template <typename T>
    void CFLOBDDTopNodeT<T>::CountNodes(Hashset<CFLOBDDNodeHandle> *visitedNodes, unsigned int &nodeCount)
    {
        rootConnection.entryPointHandle->handleContents->CountNodes(visitedNodes, nodeCount);
    }

    template <typename T>
    void CFLOBDDTopNodeT<T>::CountPaths(Hashset<CFLOBDDNodeHandle> *visitedNodes)
    {
        rootConnection.entryPointHandle->handleContents->CountPaths(visitedNodes);
    }


    template <typename T>
    void CFLOBDDTopNodeT<T>::CountNodesAndEdges(Hashset<CFLOBDDNodeHandle> *visitedNodes, Hashset<CFLOBDDReturnMapBody> *visitedEdges, 
        unsigned int &nodeCount, unsigned int &edgeCount, unsigned int& returnEdgesCount, unsigned int& returnEdgesObjCount)
    {
        rootConnection.entryPointHandle->handleContents->CountNodesAndEdges(visitedNodes, visitedEdges, nodeCount, edgeCount, returnEdgesCount);
        edgeCount += rootConnection.returnMapHandle.Size();
    }


    template <typename T>
    void CFLOBDDTopNodeT<T>::DeallocateMemory()
    {
        // CFLOBDDTopNodeT<T>::~CFLOBDDTopNodeT();
    }

    // Hash
    template <typename T>
    size_t CFLOBDDTopNodeT<T>::Hash()
    {
        return rootConnection.Hash();
    }

    // Overloaded !=
    template <typename T>
    bool CFLOBDDTopNodeT<T>::operator!= (const CFLOBDDTopNodeT<T> & C)
    {
        return rootConnection != C.rootConnection;
    }

    // Overloaded ==
    template <typename T>
    bool CFLOBDDTopNodeT<T>::operator== (const CFLOBDDTopNodeT<T> & C)
    {
        return rootConnection == C.rootConnection;
    }

    // print
    template <typename T>
    std::ostream& CFLOBDDTopNodeT<T>::print(std::ostream & out) const
    {
        out << *(rootConnection.entryPointHandle) << std::endl;
        out << rootConnection.returnMapHandle << std::endl;
        return out;
    }

    template <typename T>
    std::ostream& operator<< (std::ostream & out, const CFLOBDDTopNodeT<T> &d)
    {
        d.print(out);
        return(out);
    }

    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    // ApplyAndReduce -----------------------------------------------------
    ApplyAndReduce(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n1,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n2,
                BoolOp op
                )
    {
        // Perform 2-way cross product of n1 and n2
        PairProductMapHandle MapHandle;
        CFLOBDDNodeHandle n = PairProduct(*(n1->rootConnection.entryPointHandle),
            *(n2->rootConnection.entryPointHandle),
            MapHandle);

        // Create returnMapHandle from MapHandle: Fold the pairs in MapHandle by applying
        // [n1->rootConnection.returnMapHandle, n2->rootConnection.returnMapHandle]
        // (component-wise) to each pair.
        ReturnMapHandle<T> returnMapHandle;
        //PairProductMapBodyIterator MapIterator(*MapHandle.mapContents);
        //MapIterator.Reset();
        std::unordered_map<T, unsigned int> reduction_map;
        ReductionMapHandle reductionMapHandle;
        unsigned int iterator = 0;
        //while (!MapIterator.AtEnd()) {
        while (iterator < MapHandle.Size()){
            T c1, c2;
            int first, second;
            //first = MapIterator.Current().First();
            //second = MapIterator.Current().Second();
            first = MapHandle[iterator].First();
            second = MapHandle[iterator].Second();
            c1 = n1->rootConnection.returnMapHandle.Lookup(first);
            c2 = n2->rootConnection.returnMapHandle.Lookup(second);
            T val = op[c1][c2];
            if (reduction_map.find(val) == reduction_map.end()){
                returnMapHandle.AddToEnd(val);
                reduction_map.insert(std::make_pair(val, returnMapHandle.Size() - 1));
                reductionMapHandle.AddToEnd(returnMapHandle.Size() - 1);
            }
            else{
                reductionMapHandle.AddToEnd(reduction_map[val]);
            }
            //MapIterator.Next();
            iterator++;
        }
        returnMapHandle.Canonicalize();
        reductionMapHandle.Canonicalize();

        // Perform reduction on n, with respect to the common elements that returnMapHandle maps together
        //ReductionMapHandle inducedReductionMapHandle;
        //ReturnMapHandle<T> inducedReturnMap;
        //returnMapHandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
        //CFLOBDDNodeHandle::InitReduceCache();
        //CFLOBDDNodeHandle reduced_n = n.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
        CFLOBDDNodeHandle reduced_n = n.Reduce(reductionMapHandle, returnMapHandle.Size());
        //CFLOBDDNodeHandle::DisposeOfReduceCache();

        // Create and return CFLOBDDTopNode
        //return(new CFLOBDDTopNodeT<T>(reduced_n, inducedReturnMap));
        return(new CFLOBDDTopNodeT<T>(reduced_n, returnMapHandle));
    }


    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n1,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n2,
                T(*func)(T, T)
                )
    {
        // Perform 2-way cross product of n1 and n2
        PairProductMapHandle MapHandle;
        /*DisposeOfPairProductCache();
        InitPairProductCache();*/
        CFLOBDDNodeHandle n = PairProduct(*(n1->rootConnection.entryPointHandle),
                *(n2->rootConnection.entryPointHandle),
                MapHandle);

        // Create returnMapHandle from MapHandle: Fold the pairs in MapHandle by applying
        // [n1->rootConnection.returnMapHandle, n2->rootConnection.returnMapHandle]
        // (component-wise) to each pair.
        
        ReturnMapHandle<T> returnMapHandle;
        //PairProductMapBodyIterator MapIterator(*MapHandle.mapContents);
        //MapIterator.Reset();
        
        boost::unordered_map<T, unsigned int> reductionMap;
        ReductionMapHandle reductionMapHandle;
        unsigned int iterator = 0;
        //while (!MapIterator.AtEnd()) {
        while (iterator < MapHandle.Size()){
            T c1, c2;
            int first, second;
            //first = MapIterator.Current().First();
            //second = MapIterator.Current().Second();
            first = MapHandle[iterator].First();
            second = MapHandle[iterator].Second();
            c1 = n1->rootConnection.returnMapHandle.Lookup(first);
            c2 = n2->rootConnection.returnMapHandle.Lookup(second);
            T val = (*func)(c1, c2);
            unsigned int k = returnMapHandle.Size();
            for ( k = 0; k < returnMapHandle.Size(); k++)
            {
                if (returnMapHandle[k] == val)
                {
                    break;
                }
            }
            if (k < returnMapHandle.Size())
            {
                reductionMapHandle.AddToEnd(k);
            }
            else
            {
                returnMapHandle.AddToEnd(val);
                reductionMapHandle.AddToEnd(returnMapHandle.Size() - 1);
            }
            // if (reductionMap.find(val) == reductionMap.end()){
            //     returnMapHandle.AddToEnd(val);
            //     reductionMap.insert(std::make_pair(val, returnMapHandle.Size() - 1));
            //     reductionMapHandle.AddToEnd(returnMapHandle.Size() - 1);
            // }
            // else{
            //     reductionMapHandle.AddToEnd(reductionMap[val]);
            // }
            //MapIterator.Next();
            iterator++;
        }

            returnMapHandle.Canonicalize();
            reductionMapHandle.Canonicalize();

            // Perform reduction on n, with respect to the common elements that returnMapHandle maps together
            //ReductionMapHandle inducedReductionMapHandle;
            //ReturnMapHandle<T> inducedReturnMap;
            //returnMapHandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
            //CFLOBDDNodeHandle::InitReduceCache();
            //CFLOBDDNodeHandle reduced_n = n.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
            CFLOBDDNodeHandle reduced_n = n.Reduce(reductionMapHandle, returnMapHandle.Size());
            //CFLOBDDNodeHandle::DisposeOfReduceCache();
            //std::cout << returnMapHandle << std::endl;
                // Create and return CFLOBDDTopNode
            //return(new CFLOBDDTopNodeT<T>(reduced_n, inducedReturnMap));
            return(new CFLOBDDTopNodeT<T>(reduced_n, returnMapHandle));
    }


    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    ApplyAndReduce(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n1,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n2,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr n3,
                BoolOp3 op
                )
    {
        // Perform 3-way cross product of n1, n2, and n3
        TripleProductMapHandle MapHandle;
        CFLOBDDNodeHandle n = TripleProduct(*(n1->rootConnection.entryPointHandle),
            *(n2->rootConnection.entryPointHandle),
            *(n3->rootConnection.entryPointHandle),
            MapHandle);

        // Create returnMapHandle from MapHandle: Fold the pairs in MapHandle by applying
        // [n1->rootConnection.returnMapHandle, n2->rootConnection.returnMapHandle, n3->rootConnection.returnMapHandle]
        // (component-wise) to each triple.
        ReturnMapHandle<T> returnMapHandle;
        TripleProductMapBodyIterator MapIterator(*MapHandle.mapContents);
        MapIterator.Reset();
        while (!MapIterator.AtEnd()) {
            T c1, c2, c3;
            int first, second, third;
            first = MapIterator.Current().First();
            second = MapIterator.Current().Second();
            third = MapIterator.Current().Third();
            c1 = n1->rootConnection.returnMapHandle.Lookup(first);
            c2 = n2->rootConnection.returnMapHandle.Lookup(second);
            c3 = n3->rootConnection.returnMapHandle.Lookup(third);
            returnMapHandle.AddToEnd(op[c1][c2][c3]);
            MapIterator.Next();
        }
        returnMapHandle.Canonicalize();

        // Perform reduction on n, with respect to the common elements that returnMapHandle maps together
        ReductionMapHandle inducedReductionMapHandle;
        ReturnMapHandle<T> inducedReturnMap;
        returnMapHandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
        //CFLOBDDNodeHandle::InitReduceCache();
        CFLOBDDNodeHandle reduced_n = n.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
        //CFLOBDDNodeHandle::DisposeOfReduceCache();
        // Create and return CFLOBDDTopNode
        return(new CFLOBDDTopNodeT<T>(reduced_n, inducedReturnMap));
    }

    // \f.\g.(f + g)
    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkPlusTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g
                )
    {
        return ApplyAndReduce<T>(f, g, PlusFunc);
    }

    // Pointwise addition: \f.\g.(f + g)
    template<typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator+(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g)
    {
        try{
            return MkPlusTopNode<T>(f, g);
        }
        catch (std::exception e){
            std::cout << e.what() << std::endl;
            std::cout << "Addition" << std::endl;
            abort();
        }
    }


    // \f.\g.(f ^ g)
    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkExorTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f,
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g
    )
    {
        return ApplyAndReduce<T>(f, g, exclusiveOrOp);
    }

    // Pointwise addition: \f.\g.(f ^ g)
    template<typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator^(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g)
    {
        return MkExorTopNode<T>(f, g);
    }


    //// Left scalar-multiplication: \c:unsigned int.\g.(c * g)
    //template <typename T>
    //typename ref_ptr<CFLOBDDTopNodeT<T>>
    //MkLeftScalarTimesTopNode(unsigned int c, typename ref_ptr<CFLOBDDTopNodeT<T>> g)
    //{
    //	if (c == 1) return g;  // Nothing need be performed
    //
    //	CFLOBDDNodeHandle eph;
    //	typename ref_ptr<CFLOBDDTopNodeT<T>> ans;
    //	if (c == 0) {  // Special case
    //		eph = CFLOBDDNodeHandle::NoDistinctionNode[g->level];
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
    //	//     CFLOBDDNodeHandle::InitReduceCache();
    //	CFLOBDDNodeHandle reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
    //	//     CFLOBDDNodeHandle::DisposeOfReduceCache();
    //
    //	// Create and return CFLOBDDTopNode
    //	return(new CFLOBDDTopNodeT<T>(reduced_eph, inducedReturnMap));
    //}
    //
    //// Left scalar-multiplication: \c:int.\g.(c * g)
    //template <typename T>
    //typename ref_ptr<CFLOBDDTopNodeT<T>>
    //MkLeftScalarTimesTopNode(int c, typename ref_ptr<CFLOBDDTopNodeT<T>> g)
    //{
    //	if (c == 1) return g;  // Nothing need be performed
    //
    //	CFLOBDDNodeHandle eph;
    //	typename ref_ptr<CFLOBDDTopNodeT<T>> ans;
    //	if (c == 0) {  // Special case
    //		eph = CFLOBDDNodeHandle::NoDistinctionNode[g->level];
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
    //	//     CFLOBDDNodeHandle::InitReduceCache();
    //	CFLOBDDNodeHandle reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
    //	//     CFLOBDDNodeHandle::DisposeOfReduceCache();
    //
    //	// Create and return CFLOBDDTopNode
    //	return(new CFLOBDDTopNodeT<T>(reduced_eph, inducedReturnMap));
    //}

    //// Left scalar-multiplication: \c:int.\g.(c * g)
    //template <typename T>
    //typename ref_ptr<CFLOBDDTopNodeT<T>>
    //MkLeftScalarTimesTopNode(double c, typename ref_ptr<CFLOBDDTopNodeT<T>> g)
    //{
    //	if (c == 1) return g;  // Nothing need be performed
    //
    //	CFLOBDDNodeHandle eph;
    //	typename ref_ptr<CFLOBDDTopNodeT<T>> ans;
    //	if (c == 0) {  // Special case
    //		eph = CFLOBDDNodeHandle::NoDistinctionNode[g->level];
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
    //	//     CFLOBDDNodeHandle::InitReduceCache();
    //	CFLOBDDNodeHandle reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
    //	//     CFLOBDDNodeHandle::DisposeOfReduceCache();
    //
    //	// Create and return CFLOBDDTopNode
    //	return(new CFLOBDDTopNodeT<T>(reduced_eph, inducedReturnMap));
    //}

    // Left scalar-multiplication: \c:int.\g.(c * g)
    template <typename T, typename T1>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkLeftScalarTimesTopNode(T1 c, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g)
    {
        if (c == 1) return g;  // Nothing need be performed

        CFLOBDDNodeHandle eph;
        typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr ans;
        if (c == 0) {  // Special case
            eph = CFLOBDDNodeHandle::NoDistinctionNode[g->level];
            ReturnMapHandle<T> rmh = c * g->rootConnection.returnMapHandle;
            // Perform reduction on eph, with respect to the common elements that rmh maps tssogether
            ReductionMapHandle inducedReductionMapHandle;
            ReturnMapHandle<T> inducedReturnMap;
            rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
            //CFLOBDDNodeHandle::InitReduceCache();
            CFLOBDDNodeHandle reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
            //CFLOBDDNodeHandle::DisposeOfReduceCache();

            // Create and return CFLOBDDTopNode
            return(new CFLOBDDTopNodeT<T>(reduced_eph, inducedReturnMap));
        }
        else {
            eph = *(g->rootConnection.entryPointHandle);
            ReturnMapHandle<T> rmh = c * g->rootConnection.returnMapHandle;
            return(new CFLOBDDTopNodeT<T>(eph, rmh));
        }
    }



    // Left scalar-multiplication: \c:int.\g.(c * g)
    template<typename T, typename T1>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator*(T1 c, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g)
    {
        try{
            return MkLeftScalarTimesTopNode<T, T1>(c, g);
        }
        catch (std::exception e){
            std::cout << e.what() << std::endl;
            std::cout << "Multiply " << c << std::endl;
            throw e;
        }
    }

    // Right scalar-multiplication: \f.\c:int.(f * c)
    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkRightScalarTimesTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, int c)
    {
        if (c == 1) return f;  // Nothing need be performed

        CFLOBDDNodeHandle eph;
        typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr ans;
        if (c == 0) {  // Special case
            eph = CFLOBDDNodeHandle::NoDistinctionNode[f->level];
            ReturnMapHandle<T> rmh = f->rootConnection.returnMapHandle * c;

            // Perform reduction on eph, with respect to the common elements that rmh maps together
            ReductionMapHandle inducedReductionMapHandle;
            ReturnMapHandle<T> inducedReturnMap;
            rmh.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
            //CFLOBDDNodeHandle::InitReduceCache();
            CFLOBDDNodeHandle reduced_eph = eph.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
            //CFLOBDDNodeHandle::DisposeOfReduceCache();

            // Create and return CFLOBDDTopNode
            return(new CFLOBDDTopNodeT<T>(reduced_eph, inducedReturnMap));
        }
        else {
            eph = *(f->rootConnection.entryPointHandle);
            ReturnMapHandle<T> rmh = f->rootConnection.returnMapHandle * c;
            return(new CFLOBDDTopNodeT<T>(eph, rmh));
        }
    }

    // Right scalar-multiplication: \f.\c:int.(f * c)
    template<typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator*(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, int c)
    {
        return MkRightScalarTimesTopNode<T>(f, c);
    }

    // Pointwise multiplication: \f.\g.(f * g) ---------------------------------------
    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    MkTimesTopNode(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f,
                typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g
                )
    {
        return ApplyAndReduce<T>(f, g, TimesFunc);
    }

    // Pointwise multiplication: \f.\g.(f * g) 
    template <typename T>
    typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr
    operator*(typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr f, typename CFLOBDDTopNodeT<T>::CFLOBDDTopNodeTRefPtr g)
    {
        return MkTimesTopNode<T>(f, g);
    }


}