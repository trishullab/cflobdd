#ifndef HASHTABGUARD
#define HASHTABGUARD

#include <iostream>
#include "apvector_T.h"
#include "list_T.h"
#include "list_TPtr.h"
#include "ref_ptr.h"


// **********************************************************************
// The Hashtable template implements a set of <key,item> pairs, with fast
// lookup and insert operations.
//
// NOTE: The key type (KeyT) must provide:
//       1. a Hash member function that takes one parameter k and
//          returns an integer in the range 0..k-1
//       2. a definition of std::ostream &operator<<(std::ostream, const KeyT &)
//       3. a definition of bool operator==(const KeyT &)
//       The item type (ItemT) must also provide operator<<.
//
// To declare a Hashtable of <key,item> pairs of types KeyT and ItemT use:
// Hashtable<KeyT, ItemT>,
// e.g., Hashtable<int,string> table;
// **********************************************************************
//
// Member functions
// ================
//
// constructors/destructor
// =======================
// Hashtable()                    -- constructor: creates an empty hashtable
// Hashtable(int size)            -- constructor: creates an empty hashtable
// Hashtable(const Hashtable & H) -- copy constructor: creates a copy of H
// ~Hashtable()                   -- destructor: cleans up as necessary
//
// mutators/modifiers
// ==================
// void Insert(keyT k, ItemT item) -- insert <k,item> (no check whether there
//                                     is already a key k)
//
// other operations
// ================
// int Size()                -- number of items currently in hashtable
// bool Lookup(KeyT k)       -- true iff a pair with key k is in the hashtable
// bool Fetch(KeyT k, ItemT &i)
//                           -- returns true iff hashtable has an entry <k,item>, and
//                              places item in i; otherwise returns false
// Hashtable & operator =    -- (assignment)
//                                               
// **********************************************************************

template <class KeyT, class ItemT> class HashtableIterator;  // forward declaration

// Pair type -- returned by the Hashtable iterator
template <class KeyT, class ItemT> class Pair {
 public: 
  // constructor
  Pair(KeyT k, ItemT i): key(k), item(i) { }
  ~Pair(){ }
  // data members
  KeyT key;
  ItemT item;
};

template <class KeyT, class ItemT> class Hashtable
{
  friend class HashtableIterator<KeyT, ItemT>;

  public:
    Hashtable();                       // constructor (default)
    Hashtable(int size);               // constructor
    ~Hashtable();                      // destructor
    Hashtable(const Hashtable<KeyT, ItemT> & H); // copy constructor
	RefCounter count;
	void DeallocateMemory();
    // mutator/modifier member functions

    void Insert(KeyT k, ItemT item);
        
    // other operations

    int Size() const;
    bool Lookup(KeyT k) const;
    bool Fetch(KeyT k, ItemT &i) const;
	//ItemT & operator [] (KeyT k);
    Hashtable & operator = (const Hashtable<KeyT, ItemT> & H); // assignment
	std::ostream& print(std::ostream & out);

  private:
    const int numBuckets;
    int mySize;                        // current size
    apvector<List<Pair<KeyT, ItemT> *> > *myItems; // pointer to array of items
};

// **********************************************************************
//
// HashtableIterator class template
//
// **********************************************************************
template <class KeyT, class ItemT> class HashtableIterator
{
  public:
  // constructor
    HashtableIterator(const Hashtable<KeyT, ItemT> & H);

  // next item  
    Pair<KeyT, ItemT> *Next();
	bool hasNext();

  private:
    apvector<List<Pair<KeyT, ItemT> *> > *myItems; // pointer to array of items
    int mySize;                        // current size
    int numBuckets;
    int k;
    int numAccessed;
    ListIterator<Pair<KeyT, ItemT> *> *myCurrList;

    // private operations
    int FirstNonEmptyList();
    int NextNonEmptyList(int k);
};

#include "hash.cpp"

#endif
