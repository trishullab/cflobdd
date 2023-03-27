

#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "return_map_T.h"
#include "weighted_values_list.h"
#include "fourier_semiring.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "intpair.h"
#include "hashset.h"
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>


//***************************************************************
// WeightedValuesListBody<T>
//***************************************************************

// Constructor
template <typename T>
WeightedValuesListBody<T>::WeightedValuesListBody()
  : refCount(0), isCanonical(false), isAllSame(true), isOneOrZero(true)
{
}

template <typename T>
WeightedValuesListBody<T>::WeightedValuesListBody(unsigned int capacity)
	: refCount(0), isCanonical(false), isAllSame(true), isOneOrZero(true)
{
	mapArray.reserve(capacity);
}

template <typename T>
void WeightedValuesListBody<T>::IncrRef()
{
  refCount++;    // Warning: Saturation not checked
}

template <typename T>
void WeightedValuesListBody<T>::DecrRef()
{
  if (--refCount == 0) {    // Warning: Saturation not checked
    if (isCanonical) {
      WeightedValuesListHandle<T>::canonicalWeightedValuesListBodySet->DeleteEq(this);
		//canonicalWeightedValuesListBody<T>Set.erase(this);
    }
    delete this;
  }
}

template <typename T>
unsigned int WeightedValuesListBody<T>::Hash(unsigned int modsize)
{
  boost::hash<T> boost_hash;
  unsigned int hvalue = 0;
  for (unsigned int i = 0; i < mapArray.size(); i++){
	  hvalue = (117 * (hvalue + 1) + boost_hash(mapArray[i])) % modsize;
  }

  return hvalue;
}

template <typename T>
void WeightedValuesListBody<T>::setHashCheck()
{
    boost::hash<T> boost_hash;
	unsigned int hvalue = 0;

	for (auto &i : mapArray) {
		hvalue = (117 * (hvalue + 1) + boost_hash((i)));
	}
	hashCheck = hvalue;
}

template <typename T>
void WeightedValuesListBody<T>::AddToEnd(T y)
{
    // TODO: To change this
    if (mapArray.size() == 0)
    {
        isAllSame = true;
        value = y;
    }
    else if (mapArray.size() > 0 && isAllSame == true)
    {
        isAllSame = isAllSame && (y == mapArray[mapArray.size()-1]);
    }
    if (isOneOrZero)
        isOneOrZero = isOneOrZero && ((y == 0) || (y == 1));
	mapArray.push_back(y);
}

template <>
void WeightedValuesListBody<fourierSemiring>::AddToEnd(fourierSemiring y)
{
    // TODO: To change this
    if (mapArray.size() == 0)
    {
        isAllSame = true;
        value = y;
    }
    else if (mapArray.size() > 0 && isAllSame == true)
    {
        isAllSame = isAllSame && (y == mapArray[mapArray.size()-1]);
    }
    if (isOneOrZero){
        if (!y.isComplexValueSet)
            isOneOrZero = isOneOrZero && ((y == fourierSemiring(0, 1)) || (y == fourierSemiring(1, 1)));
        else if (y.isComplexValueSet)
            isOneOrZero = isOneOrZero && (y.complex_value == BIG_COMPLEX(0) || y.complex_value == BIG_COMPLEX(1));
    }
	mapArray.push_back(y);
}

template <typename T>
bool WeightedValuesListBody<T>::operator==(const WeightedValuesListBody<T> &o) const
{
	if (hashCheck != o.hashCheck)
		return false;

	if (mapArray.size() != o.mapArray.size())
		return false;

	for (unsigned int i = 0; i < mapArray.size(); i++){
		if (mapArray[i] != o.mapArray[i])
			return false;
	}
	return true;
}

template <typename T>
T& WeightedValuesListBody<T>::operator[](unsigned int i){                       // Overloaded []
	return mapArray[i];
}

template <typename T>
unsigned int WeightedValuesListBody<T>::Size(){
	return (unsigned int)mapArray.size();
}

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedValuesListBody<T> &r)
{
  //out << (List<int>&)r;
	for (unsigned int i = 0; i < r.mapArray.size(); i++)
	{
		out << r.mapArray[i] << " ";
	}
  return(out);
}

//***************************************************************
// WeightedValuesListHandle<T>
//***************************************************************

// Initializations of static members ---------------------------------
template <typename T>
Hashset<WeightedValuesListBody<T>> *WeightedValuesListHandle<T>::canonicalWeightedValuesListBodySet = new Hashset<WeightedValuesListBody<T>>(HASHSET_NUM_BUCKETS);

// Default constructor
template <typename T>
WeightedValuesListHandle<T>::WeightedValuesListHandle()
  :  mapContents(new WeightedValuesListBody<T>)
{
  mapContents->IncrRef();
}

// Destructor
template <typename T>
WeightedValuesListHandle<T>::~WeightedValuesListHandle()
{
  mapContents->DecrRef();
}

// Copy constructor
template <typename T>
WeightedValuesListHandle<T>::WeightedValuesListHandle(const WeightedValuesListHandle<T> &r)
  :  mapContents(r.mapContents)
{
  mapContents->IncrRef();
}

template <typename T>
WeightedValuesListHandle<T>::WeightedValuesListHandle(unsigned int capacity)
	: mapContents(new WeightedValuesListBody<T>(capacity))
{
	mapContents->IncrRef();
}

// Overloaded assignment
template <typename T>
WeightedValuesListHandle<T>& WeightedValuesListHandle<T>::operator= (const WeightedValuesListHandle<T> &r)
{
  if (this != &r)      // don't assign to self!
  {
    WeightedValuesListBody<T> *temp = mapContents;
    mapContents = r.mapContents;
    mapContents->IncrRef();
    temp->DecrRef();
  }
  return *this;        
}

// Overloaded !=
template <typename T>
bool WeightedValuesListHandle<T>::operator!=(const WeightedValuesListHandle<T> &r)
{
  return (mapContents != r.mapContents);
}

// Overloaded ==
template <typename T>
bool WeightedValuesListHandle<T>::operator==(const WeightedValuesListHandle<T> &r)
{
  return (mapContents == r.mapContents);
}

// print
template <typename T>
std::ostream& WeightedValuesListHandle<T>::print(std::ostream & out) const
{
  out << *mapContents << std::endl;
  return out;
}

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<T> &r)
{
  r.print(out);
  return(out);
}

template <typename T>
unsigned int WeightedValuesListHandle<T>::Hash(unsigned int modsize)
{
  return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

template <typename T>
unsigned int WeightedValuesListHandle<T>::Size()
{
	return mapContents->Size();
}

template <typename T>
void WeightedValuesListHandle<T>::AddToEnd(T y)
{
  assert(mapContents->refCount <= 1);
  mapContents->AddToEnd(y);
}

template <typename T>
T& WeightedValuesListHandle<T>::operator[](unsigned int i){                       // Overloaded []
	return mapContents->mapArray[i];
}

template <typename T>
T& WeightedValuesListHandle<T>::Lookup(int x)
{
	return mapContents->mapArray[x];
}

template <typename T>
int WeightedValuesListHandle<T>::LookupInv(T y)
{
	for (unsigned int i = 0; i < Size(); i++){
		if (mapContents->mapArray[i] == y)
			return i;
	}
	return -1;
}

template <typename T>
void WeightedValuesListHandle<T>::Canonicalize()
{
  WeightedValuesListBody<T> *answerContents;
  mapContents->setHashCheck();

  if (!mapContents->isCanonical) {
	unsigned int hash = canonicalWeightedValuesListBodySet->GetHash(mapContents);
    answerContents = canonicalWeightedValuesListBodySet->Lookup(mapContents, hash);
    if (answerContents == NULL) {
      canonicalWeightedValuesListBodySet->Insert(mapContents, hash);
      mapContents->isCanonical = true;
    }
    else {
      answerContents->IncrRef();
      mapContents->DecrRef();
      mapContents = answerContents;
    }
  }
}

template <typename T>
std::size_t hash_value(const WeightedValuesListHandle<T>& val)
{
	return val.mapContents->hashCheck;
}

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
template class WeightedValuesListHandle<BIG_FLOAT>;
template std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<BIG_FLOAT> &r);


typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
template class WeightedValuesListHandle<BIG_COMPLEX_FLOAT>;
template std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<BIG_COMPLEX_FLOAT> &r);

template class WeightedValuesListHandle<fourierSemiring>;
template std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<fourierSemiring> &r);

