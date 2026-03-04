#ifndef REDUCTION_MAP_GUARD
#define REDUCTION_MAP_GUARD

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

#include <iostream>
#include <fstream>
#include <unordered_set>
#include "list_T.h"
#include "list_TPtr.h"
#include "hashset.h"
#include "intpair.h"

class ReductionMapHandle;
class ReductionMapBody;


//***************************************************************
// ReductionMapBody
//***************************************************************

class ReductionMapBody {

  friend class ReductionMapHandle;

 public:
  ReductionMapBody();    // Constructor
  ReductionMapBody(unsigned int capacity);
  void IncrRef();
  void DecrRef();
  size_t Hash();
  void setHashCheck();
  void AddToEnd(int y);          // Override AddToEnd
  unsigned int refCount;         // reference-count value
  bool isIdentityMap;            // Is this ReductionMapBody an identity map?
  std::vector<int> mapArray;
  bool operator==(const ReductionMapBody &o) const;
  int& operator[](unsigned int i);                       // Overloaded []
  unsigned int Size();

 //protected:
  unsigned int hashCheck;
  bool isCanonical;              // Is this ReductionMapBody in *canonicalReductionMapBodySet?

};

std::ostream& operator<< (std::ostream & out, const ReductionMapBody &r);

// Functors for std::unordered_set<ReductionMapBody*, ...>
struct ReductionMapBodyPtrHash {
    size_t operator()(ReductionMapBody* p) const { return p->Hash(); }
};

struct ReductionMapBodyPtrEq {
    bool operator()(ReductionMapBody* a, ReductionMapBody* b) const { return *a == *b; }
};


//***************************************************************
// ReductionMapHandle
//***************************************************************

class ReductionMapHandle {
 public:
  ReductionMapHandle();                               // Default constructor
  ~ReductionMapHandle();                              // Destructor
  ReductionMapHandle(const ReductionMapHandle &r);             // Copy constructor
  ReductionMapHandle(unsigned int capacity);
  ReductionMapHandle& operator= (const ReductionMapHandle &r); // Overloaded assignment
  bool operator!= (const ReductionMapHandle &r);      // Overloaded !=
  bool operator== (const ReductionMapHandle &r);      // Overloaded ==
  size_t Hash();
  unsigned int Size();
  void AddToEnd(int y);
  int& operator[](unsigned int i);                       // Overloaded []
  int Lookup(int x);
  intpair Lookup(intpair& x);
  int LookupInv(int y);
  void Canonicalize();
  ReductionMapBody *mapContents;

  using CanonicalReductionMapBodySet = std::unordered_set<ReductionMapBody*,
                                                           ReductionMapBodyPtrHash,
                                                           ReductionMapBodyPtrEq>;
  static CanonicalReductionMapBodySet *canonicalReductionMapBodySet;

  std::ostream& print(std::ostream & out = std::cout) const;

 private:
  static CanonicalReductionMapBodySet *initCanonicalSet();
};

std::ostream& operator<< (std::ostream & out, const ReductionMapHandle &r);

#endif
