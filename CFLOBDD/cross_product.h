#ifndef CROSS_PRODUCT_CFL_GUARD
#define CROSS_PRODUCT_CFL_GUARD

//
//    Copyright (c) 1999 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "intpair.h"
#include "inttriple.h"

// ********************************************************************
// 2-Way Cross Product
// ********************************************************************

// Classes and types declared in this file ---------------------
namespace CFL_OBDD {
class PairProductMapHandle;
class PairProductMapBody;
// class PairProductKey;
}

//***************************************************************
// PairProductMapBodyIterator
//***************************************************************

typedef ListIterator<intpair> PairProductMapBodyIterator;

//***************************************************************
// PairProductMapHandle
//***************************************************************
namespace CFL_OBDD {
class PairProductMapHandle {
 public:
  PairProductMapHandle();                               // Default constructor
  ~PairProductMapHandle();                              // Destructor
  PairProductMapHandle(const PairProductMapHandle &r);             // Copy constructor
  PairProductMapHandle& operator= (const PairProductMapHandle &r); // Overloaded assignment
  bool operator!= (const PairProductMapHandle &r);      // Overloaded !=
  bool operator== (const PairProductMapHandle &r);      // Overloaded ==
  unsigned int Hash(unsigned int modsize);
  unsigned int Size();
  intpair& operator[](unsigned int i);                       // Overloaded []
  void AddToEnd(const intpair& p);
  bool Member(intpair& p);
  int Lookup(intpair& p);
  void Canonicalize();
  PairProductMapHandle Flip();                          // Create map with reversed entries
  PairProductMapBody *mapContents;
};

//***************************************************************
// PairProductMapBody
//***************************************************************

class PairProductMapBody {//: public List<intpair> {

  friend void PairProductMapHandle::Canonicalize();

 public:
  PairProductMapBody();    // Constructor
  void IncrRef();
  void DecrRef();
  unsigned int Hash(unsigned int modsize);
  unsigned int refCount;         // reference-count value
  void setHashCheck();
  void AddToEnd(const intpair& y);          // Override AddToEnd
  std::vector<intpair> mapArray;
  bool operator==(const PairProductMapBody &p) const;
  intpair& operator[](unsigned int i);                       // Overloaded []
  unsigned int Size();
  unsigned int hashCheck;
 public:
  bool isCanonical;              // Is this PairProductMapBody in *canonicalPairProductMapBodySet?
  static Hashset<PairProductMapBody> *canonicalPairProductMapBodySet;

};

std::ostream& operator<< (std::ostream & out, const PairProductMapBody &r);
}
//***************************************************************
// PairProductKey
//***************************************************************

namespace CFL_OBDD{

class PairProductKey {

 public:
  PairProductKey(CFLOBDDNodeHandle nodeHandle1, CFLOBDDNodeHandle nodeHandle2); // Constructor
  unsigned int Hash(unsigned int modsize);
  PairProductKey& operator= (const PairProductKey& p);  // Overloaded assignment
  bool operator!= (const PairProductKey& p);        // Overloaded !=
  bool operator== (const PairProductKey& p);        // Overloaded ==
  CFLOBDDNodeHandle NodeHandle1() const { return nodeHandle1; }      // Access function
  CFLOBDDNodeHandle NodeHandle2() const { return nodeHandle2; }      // Access function
  std::ostream& print(std::ostream & out) const;
  ~PairProductKey(){}
 private:
  CFLOBDDNodeHandle nodeHandle1;
  CFLOBDDNodeHandle nodeHandle2;
  PairProductKey();                                 // Default constructor (hidden)
};

std::ostream& operator<< (std::ostream & out, const PairProductKey &p);

//***************************************************************
// PairProductMemo
//***************************************************************

class PairProductMemo {

 public:
  PairProductMemo();                                 // Default constructor
  PairProductMemo(CFLOBDDNodeHandle nodeHandle, PairProductMapHandle pairProductMapHandle); // Constructor
  PairProductMemo& operator= (const PairProductMemo& p);  // Overloaded assignment
  bool operator!= (const PairProductMemo& p);        // Overloaded !=
  bool operator== (const PairProductMemo& p);        // Overloaded ==

  CFLOBDDNodeHandle nodeHandle;
  PairProductMapHandle pairProductMapHandle;
};

// Auxiliary functions -----------------------------------------------
CFLOBDDNodeHandle PairProduct(CFLOBDDNodeHandle n1,
                              CFLOBDDNodeHandle n2,
                              PairProductMapHandle &pairProductMap
                             );
CFLOBDDNodeHandle PairProduct(CFLOBDDInternalNode *n1,
                              CFLOBDDInternalNode *n2,
                              PairProductMapHandle &pairProductMap
                             );


void InitPairProductCache();
void DisposeOfPairProductCache();
}
// ********************************************************************
// 3-Way Cross Product
// ********************************************************************

// Classes and types declared in this file ---------------------
namespace CFL_OBDD {
class TripleProductMapHandle;
class TripleProductMapBody;
class TripleProductKey;
}
//***************************************************************
// TripleProductMapBodyIterator
//***************************************************************

typedef ListIterator<inttriple> TripleProductMapBodyIterator;

//***************************************************************
// TripleProductMapHandle
//***************************************************************
namespace CFL_OBDD {
class TripleProductMapHandle {
 public:
  TripleProductMapHandle();                               // Default constructor
  ~TripleProductMapHandle();                              // Destructor
  TripleProductMapHandle(const TripleProductMapHandle &r);             // Copy constructor
  TripleProductMapHandle& operator= (const TripleProductMapHandle &r); // Overloaded assignment
  bool operator!= (const TripleProductMapHandle &r);      // Overloaded !=
  bool operator== (const TripleProductMapHandle &r);      // Overloaded ==
  unsigned int Hash(unsigned int modsize);
  unsigned int Size();
  void AddToEnd(inttriple t);
  bool Member(inttriple t);
  int Lookup(inttriple t);
  void Canonicalize();
  TripleProductMapBody *mapContents;
};

//***************************************************************
// TripleProductMapBody
//***************************************************************

class TripleProductMapBody : public List<inttriple> {

  friend void TripleProductMapHandle::Canonicalize();

 public:
  TripleProductMapBody();    // Constructor
  void IncrRef();
  void DecrRef();
  unsigned int Hash(unsigned int modsize);
  unsigned int refCount;         // reference-count value

 public:
  bool isCanonical;              // Is this TripleProductMapBody in *canonicalTripleProductMapBodySet?
  static Hashset<TripleProductMapBody> *canonicalTripleProductMapBodySet;

};
}
std::ostream& operator<< (std::ostream & out, const CFL_OBDD::TripleProductMapBody &r);



//***************************************************************
// TripleProductKey
//***************************************************************
namespace CFL_OBDD {
class TripleProductKey {

 public:
  TripleProductKey(CFLOBDDNodeHandle nodeHandle1, CFLOBDDNodeHandle nodeHandle2, CFLOBDDNodeHandle nodeHandle3); // Constructor
  unsigned int Hash(unsigned int modsize);
  TripleProductKey& operator= (const TripleProductKey& p);  // Overloaded assignment
  bool operator!= (const TripleProductKey& p);        // Overloaded !=
  bool operator== (const TripleProductKey& p);        // Overloaded ==
  CFLOBDDNodeHandle NodeHandle1() const { return nodeHandle1; }      // Access function
  CFLOBDDNodeHandle NodeHandle2() const { return nodeHandle2; }      // Access function
  CFLOBDDNodeHandle NodeHandle3() const { return nodeHandle3; }      // Access function
  std::ostream& print(std::ostream & out) const;

 private:
  CFLOBDDNodeHandle nodeHandle1;
  CFLOBDDNodeHandle nodeHandle2;
  CFLOBDDNodeHandle nodeHandle3;
  TripleProductKey();                                 // Default constructor (hidden)
};

std::ostream& operator<< (std::ostream & out, const TripleProductKey &p);

//***************************************************************
// TripleProductMemo
//***************************************************************

class TripleProductMemo {

 public:
  TripleProductMemo();                                 // Default constructor
  TripleProductMemo(CFLOBDDNodeHandle nodeHandle, TripleProductMapHandle tripleProductMapHandle); // Constructor
  TripleProductMemo& operator= (const TripleProductMemo& p);  // Overloaded assignment
  bool operator!= (const TripleProductMemo& p);        // Overloaded !=
  bool operator== (const TripleProductMemo& p);        // Overloaded ==

  CFLOBDDNodeHandle nodeHandle;
  TripleProductMapHandle tripleProductMapHandle;
};


// Auxiliary functions -----------------------------------------------
CFLOBDDNodeHandle TripleProduct(CFLOBDDNodeHandle n1,
                                CFLOBDDNodeHandle n2,
                                CFLOBDDNodeHandle n3,
                                TripleProductMapHandle &tripleProductMap
                               );
CFLOBDDNodeHandle TripleProduct(CFLOBDDInternalNode *n1,
                                CFLOBDDInternalNode *n2,
                                CFLOBDDInternalNode *n3,
                                TripleProductMapHandle &tripleProductMap
                               );


void InitTripleProductCache();
void DisposeOfTripleProductCache();

} // namespace CFL_OBDD

#endif //CROSS_PRODUCT_CFL_GUARD
