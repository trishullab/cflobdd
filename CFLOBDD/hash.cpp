#ifndef _HASH_CPP_GUARD
#define _HASH_CPP_GUARD

#include <cstdlib>
#include "hash.h"
#include "list_T.h"
#include "list_TPtr.h"

const bool DEBUG_HASH = false;

static const unsigned int HASHBASE = 8388593;
const int HASH_NUM_BUCKETS = 200000;

// Implementation of Hashtable template.
// See hash.h for documentation.

// **********************************************************************
// Hashtable constructor  1
// **********************************************************************
template<class KeyT, class ItemT> Hashtable<KeyT, ItemT>::Hashtable():
  numBuckets(1000),
  mySize(0),
  myItems(new apvector<List <Pair<KeyT, ItemT> *> > (numBuckets))
{ }

// **********************************************************************
// Hashtable constructor 2
// **********************************************************************
template<class KeyT, class ItemT> Hashtable<KeyT, ItemT>::Hashtable(int size):
  numBuckets(size),
  mySize(0),
  myItems(new apvector<List <Pair<KeyT, ItemT> *> > (numBuckets))
{ }

// **********************************************************************
// Hashtable destructor
//   free storage
// **********************************************************************
template<class KeyT, class ItemT> Hashtable<KeyT, ItemT>::~Hashtable()
{
  delete myItems;
}

// **********************************************************************
// Hashtable copy constructor
// **********************************************************************
template<class KeyT, class ItemT>
Hashtable<KeyT, ItemT>::Hashtable(const Hashtable<KeyT, ItemT> & H):
  numBuckets(H.numBuckets),
  mySize(H.mySize),
  myItems(new apvector< List<Pair<KeyT, ItemT> *> > (numBuckets))
// postcondition: table is a copy of H
{
  *myItems = *(H.myItems);
}

template<class KeyT, class ItemT>
void Hashtable<KeyT, ItemT>::DeallocateMemory()
{
	// Hashtable<KeyT, ItemT>::~Hashtable();
  delete myItems;
}

// **********************************************************************
// Hashtable assignment
// **********************************************************************
template<class KeyT, class ItemT>
Hashtable<KeyT, ItemT> & Hashtable<KeyT, ItemT>::operator = (const Hashtable & H)
// postcondition: hashtable is a copy of H
{
    if (this != & H)                      // watch aliasing
    {
      // Constant: numBuckets = H.numBuckets;
      mySize = H.mySize;
      myItems = new apvector< List <Pair<KeyT, ItemT> *> >;
      *myItems = *(H.myItems);
    }
    return * this;
}

// **********************************************************************
// Insert
//    Insert the pair <k,item>
// **********************************************************************
template<class KeyT, class ItemT>
void Hashtable<KeyT, ItemT>::Insert(KeyT k, ItemT item)
{
  int j = k.Hash(HASHBASE) % numBuckets;
  Pair<KeyT, ItemT> *p = new Pair<KeyT, ItemT>(k, item);

  if (DEBUG_HASH) {
    std::cerr << "hashcode: " << j << std::endl;
    std::cerr << "length of list: " << (*myItems)[j].Length() << std::endl;
  }
  (*myItems)[j].AddToEnd(p);
  mySize++;
}

// **********************************************************************
// Size
//   return the number of items currently in this table
// **********************************************************************
template<class KeyT, class ItemT> int Hashtable<KeyT, ItemT>::Size() const
{
  return(mySize);
}

// **********************************************************************
// Lookup
//   return true iff k is in this table
// **********************************************************************
template<class KeyT, class ItemT>
bool Hashtable<KeyT, ItemT>::Lookup(KeyT k) const
{
  int j = k.Hash(HASHBASE) % numBuckets;
  
  ListIterator<Pair<KeyT, ItemT> *> li((*myItems)[j]);
  li.Reset();
  while (!li.AtEnd()) {
    Pair<KeyT, ItemT> *p = li.Current();
    KeyT oneKey = p->key;
    if (oneKey == k) return true;
    li.Next();
  }
  return false;

  //  (*myItems)[j].Reset();
  //  while (!(*myItems)[j].AtEnd()) {
  //    Pair<KeyT, ItemT> *p = (*myItems)[j].Current();
  //    KeyT oneKey = p->key;
  //    if (oneKey == k) return true;
  //    (*myItems)[j].Next();
  //  }
  //  return false;
}

// **********************************************************************
// Fetch
//   return item = table(k)
//   otherwise, return NULL
// **********************************************************************
template<class KeyT, class ItemT>
bool Hashtable<KeyT, ItemT>::Fetch(KeyT k, ItemT &i) const
{
  int j = k.Hash(HASHBASE) % numBuckets;
  
  ListIterator<Pair<KeyT, ItemT> *> li((*myItems)[j]);
  li.Reset();
  while (!li.AtEnd()) {
    Pair<KeyT, ItemT> *p = li.Current();
    KeyT oneKey = p->key;
    if (oneKey == k) {
      i = p->item;
      return true;
    }
    li.Next();
  }
  return false;
}

//template<class KeyT, class ItemT>
//ItemT & Hashtable<KeyT, ItemT>::operator[](KeyT k){
//	
//}

// **********************************************************************
// print
// **********************************************************************
template<class KeyT, class ItemT>
std::ostream& Hashtable<KeyT, ItemT>::print(std::ostream & out)
{
  HashtableIterator<KeyT, ItemT> Iter(*this);
  
  Pair<KeyT, ItemT> *p;
  while ((p = Iter.Next()) != NULL) {
    out << p->key << ": " << p->item << std::endl;
  }
  return out;
}

// Implementation of HashtableIterator class

// **********************************************************************
// constructor
// **********************************************************************
template<class KeyT, class ItemT> HashtableIterator<KeyT, ItemT>::
HashtableIterator(const Hashtable<KeyT, ItemT> & H):
    myItems(H.myItems),
    mySize(H.mySize),
    numBuckets(H.numBuckets),
    numAccessed(0)
{
  if (mySize != 0) {
    k = FirstNonEmptyList();
    myCurrList = new ListIterator<Pair<KeyT, ItemT> *>((*myItems)[k]);
    myCurrList->Reset();
  }
}

// **********************************************************************
// Next
//   return the next <key, item> pair or NULL if at end
// **********************************************************************
template<class KeyT, class ItemT> Pair<KeyT, ItemT> *
HashtableIterator<KeyT, ItemT>::Next()
{
  if (numAccessed == mySize) return NULL;
  Pair<KeyT, ItemT> *p = myCurrList->Current();
  numAccessed++;
  if (numAccessed < mySize) {
    myCurrList->Next();
    if (myCurrList->AtEnd()) {
      k = NextNonEmptyList(k);
      myCurrList = new ListIterator<Pair<KeyT, ItemT> *>((*myItems)[k]);
      myCurrList->Reset();
    }
  }
  return p;
}

template<class KeyT, class ItemT> 
bool HashtableIterator<KeyT, ItemT>::hasNext()
{
	if (numAccessed == mySize) return 0;
	return 1;
}

// **********************************************************************
// FirstNonEmptyList
//   return the index of the first non-empty list in the given hashtable
// **********************************************************************
template<class KeyT, class ItemT>
int HashtableIterator<KeyT, ItemT>::FirstNonEmptyList()
{
  List<Pair<KeyT, ItemT> *> L;

  int k = 0;
  L = (*myItems)[k];
  while (L.Length() == 0) {
    k++;
    if (k >= numBuckets) {
      std::cerr << "HashtableIterator::FirstNonEmptyList: can't find first item" << std::endl;
      abort();
    }
    L = (*myItems)[k];
  }
  return k;
}

// **********************************************************************
// NextNonEmptyList
//   return the index of the next non-empty list in the hashtable after
//   the given one
// **********************************************************************
template<class KeyT, class ItemT>
int HashtableIterator<KeyT, ItemT>::NextNonEmptyList(int k)
{
  List<Pair<KeyT, ItemT> *> L;

  k++;
  if (k >= numBuckets) {
    std::cerr << "HashtableIterator::NextNonEmptylist can't find next list" << std::endl;
    abort();
  }
  L = (*myItems)[k];
  while (L.Length() == 0) {
    k++;
    if (k >= numBuckets) {
      std::cerr << "HashtableIterator::NextNonEmptyList can't find next list" << std::endl;
      abort();
    }
    L = (*myItems)[k];
  }
  // found a non-empty list
  return k;
}

#endif

