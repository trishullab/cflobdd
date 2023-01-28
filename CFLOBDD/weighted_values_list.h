#ifndef WEIGHTED_VALUES_LIST_GUARD
#define WEIGHTED_VALUES_LIST_GUARD


#include <iostream>
#include <fstream>
#include <unordered_set>
#include "list_T.h"
#include "list_TPtr.h"
#include "hashset.h"
#include "intpair.h"

template <typename>
class WeightedValuesListHandle;
template <typename>
class WeightedValuesListBody;


//***************************************************************
// WeightedValuesListHandle
//***************************************************************

template <typename T>
class WeightedValuesListHandle {
 public:
  WeightedValuesListHandle();                               // Default constructor
  ~WeightedValuesListHandle();                              // Destructor
  WeightedValuesListHandle(const WeightedValuesListHandle &r);             // Copy constructor
  WeightedValuesListHandle(unsigned int capacity);
  WeightedValuesListHandle& operator= (const WeightedValuesListHandle &r); // Overloaded assignment
  bool operator!= (const WeightedValuesListHandle &r);      // Overloaded !=
  bool operator== (const WeightedValuesListHandle &r);      // Overloaded ==
  unsigned int Hash(unsigned int modsize);
  unsigned int Size();
  void AddToEnd(T y);
  T& operator[](unsigned int i);                       // Overloaded []
  T& Lookup(int x);
  int LookupInv(T y);
  void Canonicalize();
  WeightedValuesListBody<T> *mapContents;
  static Hashset<WeightedValuesListBody<T>> *canonicalWeightedValuesListBodySet;
  std::ostream& print(std::ostream & out = std::cout) const;
};

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<T> &r);
template <typename T>
extern std::size_t hash_value(const WeightedValuesListHandle<T>& val);

//***************************************************************
// WeightedValuesListBody
//***************************************************************

template <typename T>
class WeightedValuesListBody {//: public List<int> {

  friend void WeightedValuesListHandle<T>::Canonicalize();

 public:
  WeightedValuesListBody();    // Constructor
  WeightedValuesListBody(unsigned int capacity);
  void IncrRef();
  void DecrRef();
  unsigned int Hash(unsigned int modsize);
  void setHashCheck();
  void AddToEnd(T y);          // Override AddToEnd
  unsigned int refCount;         // reference-count value
  std::vector<T> mapArray;
  bool isAllSame;
  bool isOneOrZero;
  T value; // the value that is same in all
  bool operator==(const WeightedValuesListBody &o) const;
  T& operator[](unsigned int i);                       // Overloaded []
  unsigned int Size();

 //protected:
  unsigned int hashCheck;
  bool isCanonical;              // Is this WeightedValuesListBody in *canonicalWeightedValuesListBodySet?

  struct PointerHash {
  public:
	  size_t operator()(const WeightedValuesListBody* r){
		  return r->hashCheck;
	  }
  };

  struct PointerEqual {
  public:
	  bool operator()(const WeightedValuesListBody* r1, const WeightedValuesListBody* r2){
		  return (r1->hashCheck == r2->hashCheck);
	  };
  };

};

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedValuesListBody<T> &r);
//static std::unordered_set<WeightedValuesListBody *, WeightedValuesListBody::PointerHash, WeightedValuesListBody::PointerEqual> canonicalWeightedValuesListBodySet;
#endif