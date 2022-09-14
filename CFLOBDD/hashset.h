#ifndef HASHSETGUARD
#define HASHSETGUARD

#include <iostream>
#include "apvector_T.h"
#include "list_T.h"
#include "list_TPtr.h"

// **********************************************************************
// The Hashset template implements a set of "pointers to ItemT", with fast
// lookup and insert operations.
//
// NOTE: The item type (ItemT) must provide:
//       1. a Hash member function that takes one parameter k and
//          returns an integer in the range 0..k-1
//       2. a definition of std::ostream &operator<<(std::ostream, const KeyT &)
//       3. a definition of bool operator==(const ItemT &)
//
// To declare a Hashset of ItemT* use:
// Hashset<ItemT>,
// e.g., Hashset<string> table;
// **********************************************************************
//
// Member functions
// ================
//
// constructors/destructor
// =======================
// Hashset()                    -- constructor: creates an empty hashtable
// Hashset(int size)            -- constructor: creates an empty hashtable
// Hashset(const Hashset & H) -- copy constructor: creates a copy of H
// ~Hashset()                   -- destructor: cleans up as necessary
//
// mutators/modifiers
// ==================
// void Insert(ItemT *item) -- insert item (no check whether this is a duplicate)
// bool Delete(ItemT *item) -- delete item
// bool DeleteEq(ItemT *item) -- delete item (based on pointer equality)
//
// other operations
// ================
// int Size()                -- number of items currently in hashtable
// ItemT *Lookup(ItemT *item)  -- NULL iff item is not in the hashtable
// Hashset & operator =    -- (assignment)
//                                               
// **********************************************************************

template <class ItemT> class HashsetIterator;  // forward declaration

template <class ItemT> class Hashset
{
  friend class HashsetIterator<ItemT>;

  public:
    Hashset();                       // constructor (default)
    Hashset(int size);               // constructor
    ~Hashset();                      // destructor
    Hashset(const Hashset<ItemT> & H); // copy constructor

    // mutator/modifier member functions

    void Insert(ItemT *item);
	void Insert(ItemT *item, unsigned int hash);
    bool Delete(ItemT *item);
    bool DeleteEq(ItemT *item);
        
    // other operations

    int Size() const;
	unsigned int GetHash(ItemT *item) const;
    ItemT *Lookup(ItemT *item) const;
	ItemT *Lookup(ItemT *item, unsigned int hash) const;
    Hashset & operator = (const Hashset<ItemT> & H); // assignment
  public:
	  std::ostream& print(std::ostream & out = std::cout) const;
	void PrintTable(std::ostream & out = std::cout) const;

  private:
    const int numBuckets;
    int mySize;                        // current size
    apvector<List<ItemT *> > *myItems; // pointer to array of list of ItemT*'s
};

template<class ItemT> std::ostream& operator<< (std::ostream & out, const Hashset<ItemT> &ht)
{
  ht.print(out);
  return(out);
}

// **********************************************************************
//
// HashsetIterator class template
//
// **********************************************************************
template <class ItemT> class HashsetIterator
{
  public:
  // constructor
    HashsetIterator(const Hashset<ItemT> & H);

  // next item  
    ItemT *Next();
	void PrintTable(std::ostream & out);

  private:
    apvector<List<ItemT *> > *myItems; // pointer to array of items
    int mySize;                        // current size
    int numBuckets;
    int k;
    int numAccessed;
    ListIterator<ItemT *> *myCurrList;

    // private operations
    int FirstNonEmptyList();
    int NextNonEmptyList(int k);
};

#include "hashset.cpp"

#endif
