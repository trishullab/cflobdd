#ifndef _LIST_GUARD
#define _LIST_GUARD

// **********************************************************************
// The List template implements an ordered list of items with a "current" item.
// To declare a list of items of type T use: List<T>, e.g., List<int> intlist;
// **********************************************************************
//
// Member functions
// ================
//
// constructors/destructor
// =======================
// List()               -- constructor: creates a list containing no items
// List(const List & L) -- copy constructor: creates a copy of L
// ~List()              -- destructor: cleans up as necessary
//
// mutators/modifiers
// ==================
// void AddToEnd(T k) -- add k to the end of the list
// void AddToFront(T k) -- add k to the front of the list
// void Union(const List & L) -- add all items in L to the list if they aren't
//                               already there
// void Union(const T & k) -- add item k to the list if it isn't already there
// void Append(const List & L) -- add all items in L to the list
// void Intersect(const List & L) -- remove from the list any items not in L
// void RemoveCurrent() -- remove the current item (current becomes next item)
// T RemoveFirst() -- remove the first item
// void FreeList()    -- reinitialize the list to be empty, freeing storage
//
// other operations
// ================
// int Length()         -- length of list
// void Reset()         -- reset "current" item to be first item
// void Next()          -- advance current item
// bool AtEnd()         -- is current item past last item?
// bool IsEmpty()    -- is the list empty (i.e. has length 0)?
// T & Current()        -- return reference to current item (error if AtEnd)
// bool Member(T k)     -- true iff k is in the list
// List & operator =    -- Overloaded =
// bool operator!=      -- Overloaded !=
// bool operator==      -- Overloaded ==
//
// Note: to iterate through a List L use code like this:
//      L.Reset();
//    while (!L.AtEnd()) {
//      ...some code involving L.Current()...
//      L.Next();
//    }
//                                               
// **********************************************************************

#include <cassert>
#include <iostream>
#include <fstream>
#include "ref_ptr.h"

template <typename T> class ListIterator;

template <typename T> class List
{
    friend class ListIterator<T>;
  public:
    List();                  // constructor (default)
    ~List();                 // destructor
    List(const List<T> & L); // copy constructor

    // mutator/modifier member functions

    void AddToEnd(T k);
    void AddToFront(T k);
    void Union(const List<T> & L);
    void Union(const T & k);
    void Append(const List<T> & L);
    void Intersect(const List<T> &L);
    void RemoveCurrent();
    T RemoveFirst();
    void FreeList();
        
    // other operations

    unsigned int Length() const;
    void Reset();
    void Next();
    bool AtEnd() const;
    bool IsEmpty() const;
    T & Current() const;
    bool Member(T k) const;
    List & operator = (const List<T> & dl);  // Overloaded =
    bool operator!= (const List & L);        // Overloaded !=
    bool operator== (const List & L);        // Overloaded ==
    RefCounter count;

  private:
   struct ListNode {
      T item;
      ListNode *next;
    };

    unsigned int mySize;                    // current length of the list
    ListNode * myItems;            // ptr to header node
    ListNode * myCurrent;          // ptr to "current" item
    ListNode * myLast;             // ptr to last node in list
  public:
	  std::ostream& print(std::ostream & out = std::cout) const;
};

template<typename T> std::ostream& operator<< (std::ostream & out, List<T> const &l)
{
  l.print(out);
  return(out);
}

// **********************************************************************
// ListIterator class template
//
// Note: There are two ways to initialize a list iterator:
//       1. From a List (share list of items; iterator is initially AtEnd)
//       2. From another ListIterator (share list of items, new iterator's
//          current item is initialized to be the same as given iterator's
//          current item).
// **********************************************************************
template <typename T> class ListIterator
{
  public:
    ListIterator(const List<T> & L);           // constructor 1
    ListIterator(const ListIterator<T> & LI);  // constructor 2

    void Reset();
    void Next();
    bool AtEnd() const;
    T & Current() const;

  private:
    typename List<T>::ListNode * myItems;      // ptr to header node
    typename List<T>::ListNode * myCurrent;    // ptr to current item in list being iterated
};

// Implementation of List template using a linked list with a header node.

// **********************************************************************
// List constructor
//   initialize the list to be empty
// **********************************************************************
template<class T> List<T>::List():
    mySize(0),                 // contains no items, just a header node
    myItems(new ListNode),
    myCurrent(NULL),
    myLast(myItems)
{
  myItems->next = NULL;
}

// **********************************************************************
// List destructor
//   free storage
// **********************************************************************
template<class T> List<T>::~List()
{
  ListNode *tmp1 = myItems;
  ListNode *tmp2;
  
	do {
		tmp2 = tmp1->next;
		delete tmp1;
		tmp1 = tmp2;
	} while (tmp1 != NULL);
  
}

// **********************************************************************
// List copy constructor
// **********************************************************************
template<class T> List<T>::List(const List<T> & L):
  mySize(0),
  myItems(new ListNode),
  myCurrent(NULL),
  myLast(myItems)
// postcondition: list is a copy of L
{

  myItems->next = NULL;
  // add all items in list L to the end of this list
  Append(L);
}

// **********************************************************************
// List assignment
// **********************************************************************
template<class T> List<T> & List<T>::operator = (const List & L)
// postcondition: list is a copy of L
{
    if (this != & L)                      // watch aliasing
    {
      FreeList();
      Append(L);
    }
    return * this;
}

// **********************************************************************
// List inequality
// **********************************************************************
template<class T> bool List<T>::operator!= (const List & L)
{
  return !(*this == L);
}

// **********************************************************************
// List equality
// **********************************************************************
template<class T> bool List<T>::operator== (const List & L)
{
  if (Length() != L.Length())
    return false;

  bool ans = true;
  ListIterator<T> LI1(*this);
  ListIterator<T> LI2(L);
  LI1.Reset();
  LI2.Reset();
  while (!LI1.AtEnd()) {
    if (LI1.Current() != LI2.Current()) {
      return false;
    }
    LI1.Next();
    LI2.Next();
  }
  return ans;
}

// **********************************************************************
// FreeList
//    reinitialize the list to be empty, freeing storage
// **********************************************************************
template<class T> void List<T>::FreeList()
{
ListNode *tmp1 = myItems->next;
ListNode *tmp2;

  while (tmp1 != NULL) {
    tmp2 = tmp1->next;
    delete tmp1;
    tmp1 = tmp2;
  }
  myItems->next = NULL;
  myCurrent = NULL;
  myLast = myItems;
  mySize = 0;
}

// **********************************************************************
// AddToEnd
//    Add the given item k to the end of the list
// **********************************************************************
template<class T> void List<T>::AddToEnd(T k)
{
ListNode *tmp;

    tmp = new ListNode;
    tmp->item = k;
    tmp->next = NULL;
    myLast->next = tmp;
    myLast = tmp;
    mySize++;
}

// **********************************************************************
// AddToFront
//    Add the given item k to the front of the list
// **********************************************************************
template<class T> void List<T>::AddToFront(T k)
{
ListNode *tmp;

    tmp = new ListNode;
    tmp->item = k;
    tmp->next = myItems->next;
    myItems->next = tmp;
    if (mySize == 0) myLast = tmp;
    mySize++;
}

// **********************************************************************
// Union
//    Union this list with the given list
// **********************************************************************
template<class T> void List<T>::Union(const List<T> & L)
{
  T item;
  ListIterator<T> LI(L);

  LI.Reset();
  while (!LI.AtEnd()) {
    item = LI.Current();
    if (!Member(item)) {
      AddToEnd(item);
    }
    LI.Next();
  }
}

// **********************************************************************
// Union
//    Add the given item to this list if it isn't already there
// **********************************************************************
template<class T> void List<T>::Union(const T & k)
{
  if (!Member(k)) {
    AddToEnd(k);
  }
}

// **********************************************************************
// Append
//    Add all items in the given list to this list
// **********************************************************************
template<class T> void List<T>::Append(const List<T> & L)
{
  T item;
  ListIterator<T> LI(L);

  LI.Reset();
  while (!LI.AtEnd()) {
    item = LI.Current();
    AddToEnd(item);
    LI.Next();
  }
}

// **********************************************************************
// Intersect
//    Intersect this list with the given list
// **********************************************************************
template<class T> void List<T>::Intersect(const List<T> & L)
{
  T item;

  Reset();
  while (!AtEnd()) {
    item = Current();
    if (!L.Member(item)) {
      RemoveCurrent();
    }
    else Next();
  }
}

// **********************************************************************
// RemoveCurrent
//   remove the current item; new current item is old next item
//   note: if current item is the LAST item on the list, must set myLast
// **********************************************************************
template<class T> void List<T>::RemoveCurrent()
{
  ListNode *tmp = myItems;

  assert(myCurrent != NULL);
  while (tmp->next != myCurrent) {
    tmp = tmp->next;
    assert(tmp != NULL);
  }
  tmp->next = myCurrent->next;
  delete myCurrent;
  myCurrent = tmp->next;
  mySize--;
  if (myCurrent == NULL) {
    myLast = myItems;
    while (myLast->next != NULL) myLast = myLast->next;
  }
}

// **********************************************************************
// RemoveFirst
//   remove the first item;
//   if the current item is the first item on the list, the current item
//      becomes the old next item
//   if the first item is the only item on the list, must set myLast
// **********************************************************************
template<class T> T List<T>::RemoveFirst()
{
  ListNode *tmp = myItems->next;
  T answer;

  assert(tmp != NULL);
  answer = tmp->item;
  
  if (myCurrent == tmp) {
    myCurrent = tmp->next;
  }
  if (myLast == tmp) {
    myLast = tmp->next;
  }
  myItems->next = tmp->next;
  delete tmp;
  mySize--;
  return answer;
}

// **********************************************************************
// Length
//   return the length of the list
// **********************************************************************
template<class T> unsigned int List<T>::Length() const
{
  return(mySize);
}

// **********************************************************************
// Reset
//   reset the "current" item to be the first one
// **********************************************************************
template<class T> void List<T>::Reset()
{
  myCurrent = myItems->next;
}

// **********************************************************************
// Next
//   advance the "current" item
// **********************************************************************
template<class T> void List<T>::Next()
{
    if (myCurrent != NULL) myCurrent = myCurrent->next;
}

// **********************************************************************
// AtEnd
//   return true iff current item is past last item
// **********************************************************************
template<class T> bool List<T>::AtEnd() const
{
  return (NULL == myCurrent);
}

// **********************************************************************
// IsEmpty
//   return true iff list has no items
// **********************************************************************
template<class T> bool List<T>::IsEmpty() const
{
  return (0 == mySize);
}

// **********************************************************************
// Current
//   return reference to current item (error if AtEnd)
// **********************************************************************
template<class T> T & List<T>::Current() const
{
    assert(myCurrent !=  NULL);
    return( myCurrent->item );
}

// **********************************************************************
// Member
//   return true iff k is in the list
// **********************************************************************
template<class T> bool List<T>::Member(T k) const
{
  ListIterator<T> LI(*this);

  LI.Reset();
  while (!LI.AtEnd()) {
    if (LI.Current() == k) return true;
    LI.Next();
  }
  return false;
}

// print --------------------------------------------------------------
template<class T> std::ostream& List<T>::print(std::ostream & out) const
{
  out << "[";
  for (ListNode *p = myItems->next; p != NULL; p = p->next)
  {
    out << p->item;
    if (p->next != NULL)
      out << ", ";
  }
  out << "]" << std::endl;
  return out;
}

// Implementation of ListIterator class

// **********************************************************************
// ListIterator constructor 1
// **********************************************************************
template<class T> ListIterator<T>::ListIterator(const List<T> & L):

    myItems(L.myItems),
    myCurrent(NULL)
{
}

// **********************************************************************
// ListIterator constructor 2
// **********************************************************************
template<class T> ListIterator<T>::ListIterator(const ListIterator<T> & LI):

    myItems(LI.myItems),
    myCurrent(LI.myCurrent)
{
}

// **********************************************************************
// Reset
//   reset the "current" item to be the first one
// **********************************************************************
template<class T> void ListIterator<T>::Reset()
{
  myCurrent = myItems->next;
}

// **********************************************************************
// Next
//   advance the "current" item
// **********************************************************************
template<class T> void ListIterator<T>::Next()
{
    if (myCurrent != NULL) myCurrent = myCurrent->next;
}

// **********************************************************************
// AtEnd
//   return true iff current item is past last item
// **********************************************************************
template<class T> bool ListIterator<T>::AtEnd() const
{
  return (NULL == myCurrent);
}

// **********************************************************************
// Current
//   return reference to current item (error if AtEnd)
// **********************************************************************
template<class T> T & ListIterator<T>::Current() const
{
    assert(myCurrent !=  NULL);
    return( myCurrent->item );
}

#endif
