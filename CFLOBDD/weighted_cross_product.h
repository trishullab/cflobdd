#ifndef WEIGHTED_CROSS_PRODUCT_CFL_GUARD
#define WEIGHTED_CROSS_PRODUCT_CFL_GUARD

#include "intpair.h"
#include "hashset.h"
#include "pair_t.h"
#include "weighted_cflobdd_node_t.h"

// ********************************************************************
// 2-Way Cross Product
// ********************************************************************

// Classes and types declared in this file ---------------------
namespace CFL_OBDD {
template <typename>
class WeightedPairProductMapHandle;
template <typename>
class WeightedPairProductMapBody;
// class WeightedPairProductKey;
}


//***************************************************************
// WeightedPairProductMapHandle
//***************************************************************
namespace CFL_OBDD {

template <typename T>
class WeightedPairProductMapHandle {
 public:
  WeightedPairProductMapHandle();                               // Default constructor
  ~WeightedPairProductMapHandle();                              // Destructor
  WeightedPairProductMapHandle(const WeightedPairProductMapHandle<T> &r);             // Copy constructor
  WeightedPairProductMapHandle& operator= (const WeightedPairProductMapHandle<T> &r); // Overloaded assignment
  bool operator!= (const WeightedPairProductMapHandle<T> &r);      // Overloaded !=
  bool operator== (const WeightedPairProductMapHandle<T> &r);      // Overloaded ==
  unsigned int Hash(unsigned int modsize);
  unsigned int Size();
  std::pair<intpair, Pair_T<T,T>> operator[](unsigned int i);                       // Overloaded []
  void AddToEnd(const intpair& p, const Pair_T<T,T>& v);
  void AddToEnd(const intpair& p);
  bool Member(intpair& p);
  int Lookup(intpair& p, Pair_T<T,T>& v);
  void Canonicalize();
  WeightedPairProductMapHandle<T> Flip();                          // Create map with reversed entries
  WeightedPairProductMapBody<T> *mapContents;
};

//***************************************************************
// WeightedPairProductMapBody
//***************************************************************

template <typename T>
class WeightedPairProductMapBody {//: public List<intpair> {

  friend void WeightedPairProductMapHandle<T>::Canonicalize();

 public:
  WeightedPairProductMapBody();    // Constructor
  void IncrRef();
  void DecrRef();
  unsigned int Hash(unsigned int modsize);
  unsigned int refCount;         // reference-count value
  void setHashCheck();
  void AddToEnd(const intpair& y, const Pair_T<T,T>& v);          // Override AddToEnd
  void AddToEnd(const intpair& y);          // Override AddToEnd
  std::vector<intpair> mapArray;
  std::vector<Pair_T<T,T>> valueArray;
  bool operator==(const WeightedPairProductMapBody<T> &p) const;
  std::pair<intpair,Pair_T<T,T>> operator[](unsigned int i);                       // Overloaded []
  unsigned int Size();
  unsigned int hashCheck;
 public:
  bool isCanonical;              // Is this WeightedPairProductMapBody in *canonicalWeightedPairProductMapBodySet?
  static Hashset<WeightedPairProductMapBody> *canonicalWeightedPairProductMapBodySet;

};

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedPairProductMapBody<T> &r);
}
//***************************************************************
// WeightedPairProductKey
//***************************************************************

namespace CFL_OBDD{

template <typename T, typename Op>
class WeightedPairProductKey {

 public:
  WeightedPairProductKey(WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle1, WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle2); // Constructor
  WeightedPairProductKey(WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle1, WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle2, T factor1, T factor2); // Constructor
  unsigned int Hash(unsigned int modsize);
  WeightedPairProductKey& operator= (const WeightedPairProductKey& p);  // Overloaded assignment
  bool operator!= (const WeightedPairProductKey& p);        // Overloaded !=
  bool operator== (const WeightedPairProductKey& p);        // Overloaded ==
  WeightedCFLOBDDNodeHandleT<T,Op> NodeHandle1() const { return nodeHandle1; }      // Access function
  WeightedCFLOBDDNodeHandleT<T,Op> NodeHandle2() const { return nodeHandle2; }      // Access function
  T Factor1() const { return factor1; } // Access function
  T Factor2() const { return factor2; } // Access function
  std::ostream& print(std::ostream & out) const;
  ~WeightedPairProductKey(){}
 private:
  WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle1;
  WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle2;
  T factor1;
  T factor2;

  WeightedPairProductKey();                                 // Default constructor (hidden)
};

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedPairProductKey<T, Op> &p);

//***************************************************************
// WeightedPairProductMemo
//***************************************************************

template <typename T, typename Op>
class WeightedPairProductMemo {

 public:
  WeightedPairProductMemo();                                 // Default constructor
  WeightedPairProductMemo(WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle, WeightedPairProductMapHandle<T> pairProductMapHandle); // Constructor
  WeightedPairProductMemo& operator= (const WeightedPairProductMemo& p);  // Overloaded assignment
  bool operator!= (const WeightedPairProductMemo& p);        // Overloaded !=
  bool operator== (const WeightedPairProductMemo& p);        // Overloaded ==

  WeightedCFLOBDDNodeHandleT<T,Op> nodeHandle;
  WeightedPairProductMapHandle<T> pairProductMapHandle;
};

// Auxiliary functions -----------------------------------------------

template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op> PairProduct(WeightedCFLOBDDNodeHandleT<T,Op> n1,
                              WeightedCFLOBDDNodeHandleT<T,Op> n2,
                              WeightedPairProductMapHandle<T> &pairProductMap
                             );

template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op> PairProduct(WeightedCFLOBDDInternalNode<T,Op> *n1,
                              WeightedCFLOBDDInternalNode<T,Op> *n2,
                              WeightedPairProductMapHandle<T> &pairProductMap
                             );

template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op> PairProduct2(WeightedCFLOBDDNodeHandleT<T,Op> n1,
                              WeightedCFLOBDDNodeHandleT<T,Op> n2,
                              T factor1,
                              T factor2,
                              WeightedPairProductMapHandle<T> &pairProductMap
                             );

template <typename T, typename Op>
WeightedCFLOBDDNodeHandleT<T,Op> PairProduct2(WeightedCFLOBDDInternalNode<T,Op> *n1,
                              WeightedCFLOBDDInternalNode<T,Op> *n2,
                              T factor1,
                              T factor2,
                              WeightedPairProductMapHandle<T> &pairProductMap
                             );


template <typename T, typename Op>
void InitWeightedPairProductCache();
template <typename T, typename Op>
void DisposeOfWeightedPairProductCache();
}

#endif //CROSS_PRODUCT_CFL_GUARD
