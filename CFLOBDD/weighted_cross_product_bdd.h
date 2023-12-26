#ifndef WEIGHTED_CROSS_PRODUCT_BDD_GUARD
#define WEIGHTED_CROSS_PRODUCT_BDD_GUARD

#include "intpair.h"
#include "hashset.h"
#include "pair_t.h"
#include "weighted_bdd_node_t.h"

// ********************************************************************
// 2-Way Cross Product
// ********************************************************************

// Classes and types declared in this file ---------------------
namespace CFL_OBDD {
template <typename>
class WeightedBDDPairProductMapHandle;
template <typename>
class WeightedBDDPairProductMapBody;
// class WeightedBDDPairProductKey;
}


//***************************************************************
// WeightedBDDPairProductMapHandle
//***************************************************************
namespace CFL_OBDD {

template <typename T>
class WeightedBDDPairProductMapHandle {
 public:
  WeightedBDDPairProductMapHandle();                               // Default constructor
  ~WeightedBDDPairProductMapHandle();                              // Destructor
  WeightedBDDPairProductMapHandle(const WeightedBDDPairProductMapHandle<T> &r);             // Copy constructor
  WeightedBDDPairProductMapHandle& operator= (const WeightedBDDPairProductMapHandle<T> &r); // Overloaded assignment
  bool operator!= (const WeightedBDDPairProductMapHandle<T> &r);      // Overloaded !=
  bool operator== (const WeightedBDDPairProductMapHandle<T> &r);      // Overloaded ==
  unsigned int Hash(unsigned int modsize);
  unsigned int Size();
  // std::pair<intpair, Pair_T<T,T>> operator[](unsigned int i);                       // Overloaded []
  void AddToEnd(const std::vector<T>& p, const T v);
  void AddToEnd(const std::vector<T>& p);
  // bool Member(intpair& p);
  // int Lookup(intpair& p, Pair_T<T,T>& v);
  void Canonicalize();
  // WeightedBDDPairProductMapHandle<T> Flip();                          // Create map with reversed entries
  WeightedBDDPairProductMapBody<T> *mapContents;
};

//***************************************************************
// WeightedBDDPairProductMapBody
//***************************************************************

template <typename T>
class WeightedBDDPairProductMapBody {//: public List<intpair> {

  friend void WeightedBDDPairProductMapHandle<T>::Canonicalize();

 public:
  WeightedBDDPairProductMapBody();    // Constructor
  void IncrRef();
  void DecrRef();
  unsigned int Hash(unsigned int modsize);
  unsigned int refCount;         // reference-count value
  void setHashCheck();
  void AddToEnd(const std::vector<T>& y, const T v);          // Override AddToEnd
  void AddToEnd(const std::vector<T>& y);          // Override AddToEnd
  std::vector<T> mapArray;
  T factor;
  bool operator==(const WeightedBDDPairProductMapBody<T> &p) const;
  // std::pair<intpair,Pair_T<T,T>> operator[](unsigned int i);                       // Overloaded []
  unsigned int Size();
  unsigned int hashCheck;
 public:
  bool isCanonical;              // Is this WeightedBDDPairProductMapBody in *canonicalWeightedBDDPairProductMapBodySet?
  static Hashset<WeightedBDDPairProductMapBody> *canonicalWeightedBDDPairProductMapBodySet;

};

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedBDDPairProductMapBody<T> &r);
}
//***************************************************************
// WeightedBDDPairProductKey
//***************************************************************

namespace CFL_OBDD{

template <typename T, typename Op>
class WeightedBDDPairProductKey {

 public:
  WeightedBDDPairProductKey(WeightedBDDNodeHandle<T,Op> nodeHandle1, WeightedBDDNodeHandle<T,Op> nodeHandle2); // Constructor
  WeightedBDDPairProductKey(WeightedBDDNodeHandle<T,Op> nodeHandle1, WeightedBDDNodeHandle<T,Op> nodeHandle2, T factor1, T factor2); // Constructor
  unsigned int Hash(unsigned int modsize);
  WeightedBDDPairProductKey& operator= (const WeightedBDDPairProductKey& p);  // Overloaded assignment
  bool operator!= (const WeightedBDDPairProductKey& p);        // Overloaded !=
  bool operator== (const WeightedBDDPairProductKey& p);        // Overloaded ==
  WeightedBDDNodeHandle<T,Op> NodeHandle1() const { return nodeHandle1; }      // Access function
  WeightedBDDNodeHandle<T,Op> NodeHandle2() const { return nodeHandle2; }      // Access function
  T Factor1() const { return factor1; } // Access function
  T Factor2() const { return factor2; } // Access function
  std::ostream& print(std::ostream & out) const;
  ~WeightedBDDPairProductKey(){}
 private:
  WeightedBDDNodeHandle<T,Op> nodeHandle1;
  WeightedBDDNodeHandle<T,Op> nodeHandle2;
  T factor1;
  T factor2;

  WeightedBDDPairProductKey();                                 // Default constructor (hidden)
};

template <typename T, typename Op>
std::ostream& operator<< (std::ostream & out, const WeightedBDDPairProductKey<T, Op> &p);

//***************************************************************
// WeightedBDDPairProductMemo
//***************************************************************

template <typename T, typename Op>
class WeightedBDDPairProductMemo {

 public:
  WeightedBDDPairProductMemo();                                 // Default constructor
  WeightedBDDPairProductMemo(WeightedBDDNodeHandle<T,Op> nodeHandle, WeightedBDDPairProductMapHandle<T> pairProductMapHandle); // Constructor
  WeightedBDDPairProductMemo& operator= (const WeightedBDDPairProductMemo& p);  // Overloaded assignment
  bool operator!= (const WeightedBDDPairProductMemo& p);        // Overloaded !=
  bool operator== (const WeightedBDDPairProductMemo& p);        // Overloaded ==

  WeightedBDDNodeHandle<T,Op> nodeHandle;
  WeightedBDDPairProductMapHandle<T> pairProductMapHandle;
};

// Auxiliary functions -----------------------------------------------

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op> BDDPairProduct(WeightedBDDNodeHandle<T,Op> n1,
                              WeightedBDDNodeHandle<T,Op> n2,
                              unsigned int numVars,
                              WeightedBDDPairProductMapHandle<T>& pairProductMapHandle,
                              T(*func)(T, T)
                             );

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op> BDDPairProduct(WeightedBDDInternalNode<T,Op> *n1,
                              WeightedBDDInternalNode<T,Op> *n2,
                              unsigned int numVars,
                              WeightedBDDPairProductMapHandle<T>& pairProductMapHandle,
                              T(*func)(T, T)
                             );

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op> BDDPairProduct2(WeightedBDDNodeHandle<T,Op> n1,
                              WeightedBDDNodeHandle<T,Op> n2,
                              unsigned int numVars,
                              T factor1,
                              T factor2,
                              WeightedBDDPairProductMapHandle<T>& pairProductMapHandle,
                              T(*func)(T, T)
                             );

template <typename T, typename Op>
WeightedBDDNodeHandle<T,Op> BDDPairProduct2(WeightedBDDInternalNode<T,Op> *n1,
                              WeightedBDDInternalNode<T,Op> *n2,
                              unsigned int numVars,
                              T factor1,
                              T factor2,
                              WeightedBDDPairProductMapHandle<T>& pairProductMapHandle,
                              T(*func)(T, T)
                             );


template <typename T, typename Op>
void InitWeightedBDDPairProductCache();
template <typename T, typename Op>
void DisposeOfWeightedBDDPairProductCache();
}

#endif //CROSS_PRODUCT_CFL_GUARD
