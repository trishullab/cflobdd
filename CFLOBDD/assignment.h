#ifndef _ASSIGNMENT_GUARD
#define _ASSIGNMENT_GUARD

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

#include <cassert>
#include <iostream>
#include <fstream>


namespace SH_OBDD {

class Assignment;
class AssignmentIterator;

// **********************************************************************
// Assignment class
// **********************************************************************
class Assignment
{
    friend class AssignmentIterator;
  public:
    Assignment(unsigned int i);          // constructor
    ~Assignment();                       // destructor
    Assignment(const Assignment & A);    // copy constructor
    Assignment(double fake, bool x0, bool x1, bool x2)
  {
    const int sz = 3;
    size = sz; binding = new bool[sz];
    (*this)[0] = x0;
    (*this)[1] = x1;
    (*this)[2] = x2;
  }
	Assignment(double fake, bool x0, bool x1, bool x2, bool x3)
	{
		const int sz = 4;
		size = sz; binding = new bool[sz];
		(*this)[0] = x0;
		(*this)[1] = x1;
		(*this)[2] = x2;
		(*this)[3] = x3;
	}
  Assignment(unsigned int i, unsigned int bits)
    // least significant bit of bits is x0 (*this)[0], next bit is x1, etc.
    // for i bits.  i must be strictly positive.
  {
    assert(i > 0);
    size = i;
    binding = new bool[i];
    for (unsigned int j = 0; j < i; j++) {
      (*this)[j] = bits & 1;
      bits >>= 1;
    }
  }

    Assignment & operator = (const Assignment & dl); // assignment
    bool& operator[] (unsigned int i);    // indexing
  private:
    unsigned int size;                   // number of entries
    bool *binding;                       // bool binding[size];
  public:
	std::ostream& print(std::ostream & out = std::cout) const;
};

std::ostream& operator<< (std::ostream & out, const Assignment &A);

// **********************************************************************
// AssignmentIterator class
// **********************************************************************
class AssignmentIterator
{
  public:
    AssignmentIterator(const Assignment & A);
    void Reset();
    void Next();
    void Advance(unsigned int i);
    bool MoreToDo() const;
    bool & Current() const;
  private:
    bool *binding;
    unsigned int size;             // number of entries
    unsigned int current;          // current item in the iteration
};
}
#endif
