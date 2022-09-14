#ifndef _HASHSET_CPP_GUARD
#define _HASHSET_CPP_GUARD

#include <cstdlib>
#include "hashset.h"
#include "list_T.h"
#include "list_TPtr.h"

const bool DEBUG_HASHSET = false;

const unsigned int HASHSETBASE = 8388593;
const int HASHSET_NUM_BUCKETS = 200000;

// Implementation of Hashset template.
// See hashset.h for documentation.

// **********************************************************************
// Hashset constructor  1
// **********************************************************************
template<class ItemT> Hashset<ItemT>::Hashset():
  numBuckets(40000),
  mySize(0),
  myItems(new apvector<List <ItemT *> > (numBuckets))
{ }

// **********************************************************************
// Hashset constructor 2
// **********************************************************************
template<class ItemT> Hashset<ItemT>::Hashset(int size):
  numBuckets(size),
  mySize(0),
  myItems(new apvector<List <ItemT *> > (numBuckets))
{ }

// **********************************************************************
// Hashset destructor
//   free storage
// **********************************************************************
template<class ItemT> Hashset<ItemT>::~Hashset()
{
  delete myItems;
}

// **********************************************************************
// Hashset copy constructor
// **********************************************************************
template<class ItemT>
Hashset<ItemT>::Hashset(const Hashset<ItemT> & H):
  numBuckets(H.numBuckets),
  mySize(H.mySize),
  myItems(new apvector< List<ItemT *> > (numBuckets))
// postcondition: table is a copy of H
{
  *myItems = *(H.myItems);
}

// **********************************************************************
// Hashset assignment
// **********************************************************************
template<class ItemT>
Hashset<ItemT> & Hashset<ItemT>::operator = (const Hashset & H)
// postcondition: hashtable is a copy of H
{
    if (this != & H)                      // watch aliasing
    {
      // Constant: numBuckets = H.numBuckets;
      mySize = H.mySize;
      myItems = new apvector< List <ItemT *> >;
      *myItems = *(H.myItems);
    }
    return * this;
}

// **********************************************************************
// Insert
//    Insert item
// **********************************************************************
template<class ItemT>
void Hashset<ItemT>::Insert(ItemT *item)
{
  unsigned int j = item->Hash(HASHSETBASE) % numBuckets;

  if (DEBUG_HASHSET) {
    std::cerr << "hashcode: " << j << std::endl;
    std::cerr << "length of list: " << (*myItems)[j].Length() << std::endl;
  }
  (*myItems)[j].AddToFront(item);
  mySize++;
}

template<class ItemT>
void Hashset<ItemT>::Insert(ItemT *item, unsigned int hash)
{
	unsigned int j = hash % numBuckets;

	if (DEBUG_HASHSET) {
		std::cerr << "hashcode: " << j << std::endl;
		std::cerr << "length of list: " << (*myItems)[j].Length() << std::endl;
	}
	(*myItems)[j].AddToFront(item);
	mySize++;
}

// **********************************************************************
// Delete
//    Delete item
//    Return true on success; false if *item not found
// **********************************************************************
template<class ItemT>
bool Hashset<ItemT>::Delete(ItemT *item)
{
  unsigned int j = item->Hash(HASHSETBASE) % numBuckets;

  (*myItems)[j].Reset();
  while (!(*myItems)[j].AtEnd()) {
    if (*item == *(*myItems)[j].Current()) {
      (*myItems)[j].RemoveCurrent();
      mySize--;
      return true;
    }
    (*myItems)[j].Next();
  }
  return false;
}

// **********************************************************************
// DeleteEq
//    Delete item (based on pointer equality)
//    Return true on success; false if *item not found
// **********************************************************************
template<class ItemT>
bool Hashset<ItemT>::DeleteEq(ItemT *item)
{
  unsigned int j = item->Hash(HASHSETBASE) % numBuckets;

  (*myItems)[j].Reset();
  while (!(*myItems)[j].AtEnd()) {
    if (item == (*myItems)[j].Current()) {
      (*myItems)[j].RemoveCurrent();
      mySize--;
      return true;
    }
    (*myItems)[j].Next();
  }
  return false;
}

// **********************************************************************
// Size
//   return the number of items currently in this table
// **********************************************************************
template<class ItemT> int Hashset<ItemT>::Size() const
{
  return(mySize);
}

// **********************************************************************
// Lookup
//   if the item is not in the set, return NULL 
//   if the item is in the set, return a pointer to the existing item
//   Note: We cannot use (*myItems)[j].Member(item) because we wish
//   to compare object values, not pointer values
// **********************************************************************
template<class ItemT>
ItemT *Hashset<ItemT>::Lookup(ItemT *item) const
{
  unsigned int j = item->Hash(HASHSETBASE) % numBuckets;

  ListIterator<ItemT *> li((*myItems)[j]);
  li.Reset();
  while (!li.AtEnd()) {
    if (*item == *li.Current())
      return li.Current();
    li.Next();
  }
  return NULL;

  //  (*myItems)[j].Reset();
  //  while (!(*myItems)[j].AtEnd()) {
  //    if (*item == *(*myItems)[j].Current())
  //      return (*myItems)[j].Current();
  //    (*myItems)[j].Next();
  //  }
  //  return NULL;
}

template<class ItemT>
ItemT *Hashset<ItemT>::Lookup(ItemT *item, unsigned int hash) const
{
	unsigned int j = hash % numBuckets;

	ListIterator<ItemT *> li((*myItems)[j]);
	li.Reset();
	while (!li.AtEnd()) {
		if (*item == *li.Current())
			return li.Current();
		li.Next();
	}
	return NULL;
}

// **********************************************************************
// GetHash
// **********************************************************************
template<class ItemT>
unsigned int Hashset<ItemT>::GetHash(ItemT *item) const
{
	return item->Hash(HASHSETBASE);
}

// **********************************************************************
// print
// **********************************************************************
template<class ItemT>
std::ostream& Hashset<ItemT>::print(std::ostream & out) const
{
  HashsetIterator<ItemT> Iter(*this);
  
  ItemT *item;
  while ((item = Iter.Next()) != NULL) {
    out << *item << std::endl;
  }
  return out;
}

template<class ItemT>
void Hashset<ItemT>::PrintTable(std::ostream & out) const
{
	HashsetIterator<ItemT> Iter(*this);
	Iter.PrintTable(out);
}
// Implementation of HashsetIterator class

// **********************************************************************
// constructor
// **********************************************************************
template<class ItemT> HashsetIterator<ItemT>::
HashsetIterator(const Hashset<ItemT> & H):
    myItems(H.myItems),
    mySize(H.mySize),
    numBuckets(H.numBuckets),
    numAccessed(0)
{
  if (mySize != 0) {
    k = FirstNonEmptyList();
    myCurrList = new ListIterator<ItemT *>((*myItems)[k]);
    myCurrList->Reset();
  }
}

// **********************************************************************
// Next
//   return the next item or NULL if at end
// **********************************************************************
template<class ItemT> ItemT *HashsetIterator<ItemT>::Next()
{
  if (numAccessed == mySize) return NULL;
  ItemT *item = myCurrList->Current();
  numAccessed++;
  if (numAccessed < mySize) {
    myCurrList->Next();
    if (myCurrList->AtEnd()) {
      k = NextNonEmptyList(k);
      myCurrList = new ListIterator<ItemT *>((*myItems)[k]);
      myCurrList->Reset();
    }
  }
  return item;
}

template<class ItemT> void HashsetIterator<ItemT>::PrintTable(std::ostream & out)
{
  int tempSize = 0;
  while (numAccessed != mySize)
  {
	  ItemT *item = myCurrList->Current();
	  numAccessed++;
	  tempSize++;
	  if (numAccessed < mySize) {
	    myCurrList->Next();
	    if (myCurrList->AtEnd()) {
		      std::cout << k << ": " << tempSize << std::endl;
		      k = NextNonEmptyList(k);
			  tempSize = 0;
		      myCurrList = new ListIterator<ItemT *>((*myItems)[k]);
		      myCurrList->Reset();
		}
    }
  }
  std::cout << k << ": " << tempSize << std::endl;
  std::cout << "size: " << mySize << std::endl;
  std::cout << "buckets: " << numBuckets << std::endl;
}

// **********************************************************************
// FirstNonEmptyList
//   return the index of the first non-empty list in the given hashtable
// **********************************************************************
template<class ItemT>
int HashsetIterator<ItemT>::FirstNonEmptyList()
{
  List<ItemT *> L;

  int k = 0;
  L = (*myItems)[k];
  while (L.Length() == 0) {
    k++;
    if (k >= numBuckets) {
      std::cerr << "HashsetIterator::FirstNonEmptyList: can't find first item" << std::endl;
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
template<class ItemT>
int HashsetIterator<ItemT>::NextNonEmptyList(int k)
{
  List<ItemT *> L;

  k++;
  if (k >= numBuckets) {
    std::cerr << "HashsetIterator::NextNonEmptylist can't find next list" << std::endl;
    abort();
  }
  L = (*myItems)[k];
  while (L.Length() == 0) {
    k++;
    if (k >= numBuckets) {
      std::cerr << "HashsetIterator::NextNonEmptyList can't find next list" << std::endl;
      abort();
    }
    L = (*myItems)[k];
  }
  // found a non-empty list
  return k;
}

#endif

