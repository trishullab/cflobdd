#ifndef WEIGHTED_BDD_NODE_T_GUARD
#define WEIGHTED_BDD_NODE_T_GUARD

#include <iostream>
#include "hashset.h"

namespace CFL_OBDD {
    template <typename, typename>
    class WeightedBDDNode;
    template <typename, typename>
    class WeightedBDDNodeHandle;
}

namespace CFL_OBDD {
    template <typename T, typename Op>
    class WeightedBDDNodeHandle {
        public:
            WeightedBDDNodeHandle();
            WeightedBDDNodeHandle(WeightedBDDNode<T,Op> *n);
            WeightedBDDNodeHandle(const WeightedBDDNodeHandle<T,Op> &n);
            ~WeightedBDDNodeHandle();
            unsigned int Hash(unsigned int modsize);
            bool operator!= (const WeightedBDDNodeHandle &nh);
            bool operator== (const WeightedBDDNodeHandle &nh) const;
            WeightedBDDNodeHandle<T,Op>& operator= (const WeightedBDDNodeHandle &nh);
            std::ostream& print(std::ostream &out = std::cout) const;

            public:

            static WeightedBDDNodeHandle AnnhilatorLeafNode;
            static WeightedBDDNodeHandle IdentityLeafNode;

            static void InitLeafNodes(); 

            WeightedBDDNode<T,Op> *handleContents;
            static Hashset<WeightedBDDNode<T,Op>> *canonicalBDDNodeTable;
            void Canonicalize();
            void ComputeWeightofPathsAsAmpsToExits();

            struct WeightedBDDNodeHandle_Hash {
            public:
                size_t operator()(const WeightedBDDNodeHandle<T,Op>& c) const {
                    return ((reinterpret_cast<std::uintptr_t>(c.handleContents) >> 2) % 997);
                }
            };

    };

    template <typename T, typename Op>
    std::ostream& operator<< (std::ostream & out, const WeightedBDDNodeHandle<T, Op> &d);

    enum BDD_NODEKIND { INTERNAL, LEAF };

    template <typename T, typename Op>
    class WeightedBDDNode {
        friend void WeightedBDDNodeHandle<T,Op>::InitLeafNodes();
        public:
            WeightedBDDNode();
            WeightedBDDNode(long int i);
            virtual long int GetIndex() = 0;
            virtual ~WeightedBDDNode();
            virtual BDD_NODEKIND NodeKind() const = 0;
            virtual bool operator!= (const WeightedBDDNode<T,Op> &n) = 0;
            virtual bool operator== (const WeightedBDDNode<T,Op> &n) = 0;
            virtual void IncrRef() = 0;
            virtual void DecrRef() = 0;
            unsigned int GetRefCount(){ return refCount; }
            const bool IsCanonical() const { return isCanonical; }
            void SetCanonical() { isCanonical = true;  }
            virtual std::ostream& print(std::ostream &out = std::cout) const = 0;
            virtual unsigned int Hash(unsigned int modsize) = 0;
            virtual void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedBDDNodeHandle<T,Op>>* visitedNodes);
            long double weightOfPathsAsAmpsToExit;
            protected:
                unsigned int refCount;
                bool isCanonical;
    };

    template <typename T, typename Op>
    class WeightedBDDInternalNode : public WeightedBDDNode<T,Op> {
        friend void WeightedBDDNodeHandle<T,Op>::InitLeafNodes();
        public:
            WeightedBDDInternalNode(long int index);
            ~WeightedBDDInternalNode();
            long int GetIndex() {return index;}
            bool operator!= (const WeightedBDDNode<T,Op> &n);
            bool operator== (const WeightedBDDNode<T,Op> &n);
            void IncrRef();
            void DecrRef();
            BDD_NODEKIND NodeKind() const { return INTERNAL; }

            std::ostream& print(std::ostream &out = std::cout) const;
            WeightedBDDNodeHandle<T,Op> leftNode;
            WeightedBDDNodeHandle<T,Op> rightNode;
            T lweight;
            T rweight;
            long double leftWeightOfPathsAsAmpsToExit;
            long double rightWeightOfPathsAsAmpsToExit;

            unsigned int Hash(unsigned int modsize);
            void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedBDDNodeHandle<T,Op>>* visitedNodes);

            private:
                long int index;
    };

    template <typename T, typename Op>
    class WeightedBDDLeafNode : public WeightedBDDNode<T,Op> {
        friend void WeightedBDDNodeHandle<T,Op>::InitLeafNodes();
        public:
            WeightedBDDLeafNode(T value);
            ~WeightedBDDLeafNode();
            long int GetIndex() {return -1;}
            bool operator!= (const WeightedBDDNode<T,Op> &n);
            bool operator== (const WeightedBDDNode<T,Op> &n);
            void IncrRef();
            void DecrRef();
            BDD_NODEKIND NodeKind() const { return LEAF; }
            std::ostream& print(std::ostream &out = std::cout) const;
            T value;
            unsigned int Hash(unsigned int modsize);
            void ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedBDDNodeHandle<T,Op>>* visitedNodes);
    };
}

#endif

