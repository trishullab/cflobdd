#include "weighted_bdd_node_t.h"
#include "hashset.h"
#include "weighted_values.h"
#include <boost/container_hash/hash.hpp>

namespace CFL_OBDD {

template <typename T, typename Op>
WeightedBDDNodeHandle<T, Op> WeightedBDDNodeHandle<T, Op>::AnnhilatorLeafNode;
template <typename T, typename Op>
WeightedBDDNodeHandle<T, Op> WeightedBDDNodeHandle<T, Op>::IdentityLeafNode;

template <typename T, typename Op>
void WeightedBDDNodeHandle<T,Op>::InitLeafNodes()
{
    AnnhilatorLeafNode = WeightedBDDNodeHandle<T,Op>(new WeightedBDDLeafNode<T,Op>(getAnnhilatorValue<T,Op>()));
    IdentityLeafNode = WeightedBDDNodeHandle<T,Op>(new WeightedBDDLeafNode<T,Op>(getIdentityValue<T,Op>()));
}

template <typename T, typename Op>
Hashset<WeightedBDDNode<T,Op>> *WeightedBDDNodeHandle<T,Op>::canonicalBDDNodeTable = new Hashset<WeightedBDDNode<T,Op>>(HASHSET_NUM_BUCKETS);

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op>::WeightedBDDNodeHandle() : handleContents (NULL)
{}

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op>::WeightedBDDNodeHandle(WeightedBDDNode<T,Op> *n)
: handleContents(n)
{
    assert(n != NULL);
    handleContents->IncrRef();
    Canonicalize();
}

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op>::WeightedBDDNodeHandle(const WeightedBDDNodeHandle<T,Op>& n)
{
    handleContents = n.handleContents;
    if (handleContents != NULL)
        handleContents->IncrRef();
}

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op>::~WeightedBDDNodeHandle()
{
    if (handleContents != NULL) {
        handleContents->DecrRef();
    }
}

template <typename T, typename Op>
unsigned int WeightedBDDNodeHandle<T,Op>::Hash(unsigned int modsize)
{
    return ((unsigned int) reinterpret_cast<uintptr_t>(handleContents) >> 2) % modsize;
}

// Overloaded !=
template <typename T, typename Op>
bool WeightedBDDNodeHandle<T,Op>::operator!= (const WeightedBDDNodeHandle<T,Op> & C)
{
  return handleContents != C.handleContents;
}

// Overloaded ==
template <typename T, typename Op>
bool WeightedBDDNodeHandle<T,Op>::operator== (const WeightedBDDNodeHandle & C) const
{
  return handleContents == C.handleContents;
}

// Overloaded assignment
template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op> & WeightedBDDNodeHandle<T,Op>::operator= (const WeightedBDDNodeHandle<T,Op> &c)
{
  if (this != &c)      // don't assign to self!
  {
    WeightedBDDNode<T,Op> *temp = handleContents;
    handleContents = c.handleContents;
    if (handleContents != NULL) {
      handleContents->IncrRef();
    }
    if (temp != NULL) {
      temp->DecrRef();
    }
  }
  return *this;        
}

// Canonicalization --------------------------------------------
template <typename T, typename Op> 
void WeightedBDDNodeHandle<T,Op>::Canonicalize()
{
  WeightedBDDNode<T,Op> *answerContents;

  if (!handleContents->IsCanonical()) {
	unsigned int hash = canonicalBDDNodeTable->GetHash(handleContents);
    answerContents = canonicalBDDNodeTable->Lookup(handleContents, hash);
    if (answerContents == NULL) {
      canonicalBDDNodeTable->Insert(handleContents, hash);
      handleContents->SetCanonical();
    }
    else {
      answerContents->IncrRef();
      handleContents->DecrRef();
      handleContents = answerContents;
    }
  }
}

template <typename T, typename Op> 
void WeightedBDDNodeHandle<T,Op>::ComputeWeightofPathsAsAmpsToExits()
{
    Hashset<WeightedBDDNodeHandle<T,Op>> *visitedNodes = new Hashset<WeightedBDDNodeHandle<T,Op>>;
    handleContents->ComputeWeightOfPathsAsAmpsToExits(visitedNodes);
    delete visitedNodes;
}


// print
template <typename T, typename Op> 
std::ostream& WeightedBDDNodeHandle<T,Op>::print(std::ostream & out) const
{
    handleContents->print(out);
    return (out);
}

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedBDDNodeHandle<T, Op> &d)
{
    d.print(out);
    return (out);
}

// ***********
// WeightedBDDNode
// ***********

template <typename T, typename Op>
WeightedBDDNode<T,Op>::WeightedBDDNode() : refCount(0), isCanonical(false)
{}

template <typename T, typename Op>
WeightedBDDNode<T,Op>::~WeightedBDDNode() {}

template <typename T, typename Op>
void WeightedBDDNode<T,Op>::ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedBDDNodeHandle<T,Op>>* visitedNodes) {}

// ************
// WeightedBDDInternalNode
// ************

template <typename T, typename Op>
WeightedBDDInternalNode<T,Op>::WeightedBDDInternalNode(long int i) : WeightedBDDNode<T,Op>()
{
    index = i;
}

template <typename T, typename Op>
WeightedBDDInternalNode<T,Op>::~WeightedBDDInternalNode()
{
}

template <typename T, typename Op>
bool WeightedBDDInternalNode<T,Op>::operator!= (const WeightedBDDNode<T,Op> &n)
{
    return !(*this == n);
}

template <typename T, typename Op>
bool WeightedBDDInternalNode<T,Op>::operator== (const WeightedBDDNode<T,Op> &nh)
{
    if (nh.NodeKind() != INTERNAL)
        return false;
    WeightedBDDInternalNode& n = (WeightedBDDInternalNode&)nh;
    return (lweight == n.lweight) && (rweight == n.rweight) && (leftNode == n.leftNode) && (rightNode == n.rightNode) && (GetIndex() == n.GetIndex());
}

template <typename T, typename Op>
void WeightedBDDInternalNode<T,Op>::IncrRef()
{
    this->refCount++;
}

template <typename T, typename Op>
void WeightedBDDInternalNode<T,Op>::DecrRef()
{
    if (--this->refCount == 0) {    // Warning: Saturation not checked
        if (this->isCanonical) {
            WeightedBDDNodeHandle<T,Op>::canonicalBDDNodeTable->DeleteEq(this);
        }
        delete this;
    }
}

template <typename T, typename Op>
unsigned int WeightedBDDInternalNode<T,Op>::Hash(unsigned int modsize)
{
    unsigned int l_hvalue = leftNode.Hash(modsize);
    unsigned int r_hvalue = rightNode.Hash(modsize);

    boost::hash<T> boost_hash;
    return (997 * boost_hash(lweight) + 97 * boost_hash(rweight) + 117 * (l_hvalue) + r_hvalue) % modsize;
}

template <typename T, typename Op>
void WeightedBDDInternalNode<T,Op>::ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedBDDNodeHandle<T,Op>>* visitedNodes)
{
    WeightedBDDNodeHandle<T,Op>* handle = new WeightedBDDNodeHandle<T,Op>(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
        leftNode.handleContents->ComputeWeightOfPathsAsAmpsToExits(visitedNodes);
        long double lp_weight = leftNode.handleContents->weightOfPathsAsAmpsToExit;
        if (leftNode.handleContents->NodeKind() == INTERNAL)
        {
            long double skip_l = std::pow(2, leftNode.handleContents->GetIndex() - GetIndex() - 1);
            lp_weight = computeComposition<long double,std::multiplies<long double>>(lp_weight, skip_l);
        }
        rightNode.handleContents->ComputeWeightOfPathsAsAmpsToExits(visitedNodes);
        long double rp_weight = rightNode.handleContents->weightOfPathsAsAmpsToExit;
        if (rightNode.handleContents->NodeKind() == INTERNAL)
        {
            long double skip_l = std::pow(2, rightNode.handleContents->GetIndex() - GetIndex() - 1);
            rp_weight = computeComposition<long double,std::multiplies<long double>>(rp_weight, skip_l);
        }
        long double l_weight = computeComposition<long double,std::multiplies<long double>>(computeProbabilityFromAmplitude<T,Op>(lweight), lp_weight);
        long double r_weight = computeComposition<long double,std::multiplies<long double>>(computeProbabilityFromAmplitude<T,Op>(rweight), rp_weight);
		long double weight = l_weight + r_weight;
        this->weightOfPathsAsAmpsToExit = weight;
        leftWeightOfPathsAsAmpsToExit = l_weight;
        rightWeightOfPathsAsAmpsToExit = r_weight;
	}
}

template <typename T, typename Op>
std::ostream& WeightedBDDInternalNode<T,Op>::print(std::ostream & out) const
{
    if (leftNode.handleContents->NodeKind() == INTERNAL && rightNode.handleContents->NodeKind() == INTERNAL){
        WeightedBDDInternalNode<T,Op>* lNode = (WeightedBDDInternalNode<T,Op> *)leftNode.handleContents;
        WeightedBDDInternalNode<T,Op>* rNode = (WeightedBDDInternalNode<T,Op> *)rightNode.handleContents;
        out << "index: " << index << " (lindex, rindex): " << lNode->GetIndex() << " " << rNode->GetIndex() << " (" << lweight << ", " << rweight << ") " << std::endl;
        leftNode.print(out);
        rightNode.print(out);
    }
    else if (leftNode.handleContents->NodeKind() == INTERNAL && rightNode.handleContents->NodeKind() == LEAF)
    {
        WeightedBDDInternalNode<T,Op>* lNode = (WeightedBDDInternalNode<T,Op> *)leftNode.handleContents;
        WeightedBDDLeafNode<T,Op>* rNode = (WeightedBDDLeafNode<T,Op> *)rightNode.handleContents;
        out << "index: " << index << " (lindex, rindex): " << lNode->GetIndex() << " " << rNode->GetIndex() << " (" << lweight << ", " << rweight << ") " << std::endl;
        leftNode.print(out);
        rightNode.print(out); 
    }
    else if (leftNode.handleContents->NodeKind() == LEAF && rightNode.handleContents->NodeKind() == INTERNAL)
    {
        WeightedBDDLeafNode<T,Op>* lNode = (WeightedBDDLeafNode<T,Op> *)leftNode.handleContents;
        WeightedBDDInternalNode<T,Op>* rNode = (WeightedBDDInternalNode<T,Op> *)rightNode.handleContents;
        out << "index: " << index << " (lindex, rindex): " << lNode->GetIndex() << " " << rNode->GetIndex() << " (" << lweight << ", " << rweight << ") " << std::endl;
        leftNode.print(out);
        rightNode.print(out);
    }
    else if (leftNode.handleContents->NodeKind() == LEAF && rightNode.handleContents->NodeKind() == LEAF)
    {
        WeightedBDDLeafNode<T,Op>* lNode = (WeightedBDDLeafNode<T,Op> *)leftNode.handleContents;
        WeightedBDDLeafNode<T,Op>* rNode = (WeightedBDDLeafNode<T,Op> *)rightNode.handleContents;
        out << "index: " << index << " (lindex, rindex): " << lNode->GetIndex() << " " << rNode->GetIndex() << " (" << lweight << ", " << rweight << ") " << std::endl;
        leftNode.print(out);
        rightNode.print(out);
    }

    return (out);
}

// ************
// WeightedBDDLeafNode
// ************

template <typename T, typename Op>
WeightedBDDLeafNode<T,Op>::WeightedBDDLeafNode(T v) : WeightedBDDNode<T,Op>()
{
    value = v;
    this->refCount = 1;
}

template <typename T, typename Op>
WeightedBDDLeafNode<T,Op>::~WeightedBDDLeafNode()
{
}

template <typename T, typename Op>
bool WeightedBDDLeafNode<T,Op>::operator!= (const WeightedBDDNode<T,Op> &n)
{
    return !(*this == n);
}

template <typename T, typename Op>
bool WeightedBDDLeafNode<T,Op>::operator== (const WeightedBDDNode<T,Op> &n)
{
    if (n.NodeKind() != LEAF)
        return false;
    WeightedBDDLeafNode& nh = (WeightedBDDLeafNode&)n;
    return (value == nh.value);
}

template <typename T, typename Op>
void WeightedBDDLeafNode<T,Op>::IncrRef()
{
}

template <typename T, typename Op>
void WeightedBDDLeafNode<T,Op>::DecrRef()
{
}

template <typename T, typename Op>
unsigned int WeightedBDDLeafNode<T,Op>::Hash(unsigned int modsize)
{
    boost::hash<T> boost_hash;
    return ((117 * boost_hash(value) + 1) % modsize);
}

template <typename T, typename Op>
std::ostream& WeightedBDDLeafNode<T,Op>::print(std::ostream & out) const
{
    out << "value: " << value << std::endl;
    return (out);
}

template <typename T, typename Op>
void WeightedBDDLeafNode<T,Op>::ComputeWeightOfPathsAsAmpsToExits(Hashset<WeightedBDDNodeHandle<T,Op>>* visitedNodes)
{
    WeightedBDDNodeHandle<T,Op>* handle = new WeightedBDDNodeHandle<T,Op>(this);
	if (visitedNodes->Lookup(handle) == NULL) {
		visitedNodes->Insert(handle);
        this->weightOfPathsAsAmpsToExit = computeProbabilityFromAmplitude<T,Op>(value);
	}
}



}