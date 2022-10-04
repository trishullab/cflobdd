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
#include "assignment.h"
#include <fstream> // for definition of NULL

using namespace SH_OBDD;
// **********************************************************************
// Implementation of Assignment class
// **********************************************************************

// Constructor
Assignment::Assignment(unsigned int i):
    size(i),
    binding(new bool[i])
{
}

// Destructor
Assignment::~Assignment()
{
  delete [] binding;
}

// Copy constructor
Assignment::Assignment(const Assignment & A):
  size(A.size),
  binding(new bool[A.size])
{
  for (unsigned int i = 0; i < size; i++) {
    binding[i] = A.binding[i];
  }
}

// Overloaded assignment
Assignment & Assignment::operator = (const Assignment & A)
{
    if (this != & A)                      // watch aliasing
    {
      delete [] binding;
      size = A.size;
      binding = new bool[size];
      for (unsigned int i = 0; i < size; i++) {
        binding[i] = A.binding[i];
      }
    }
    return * this;
}

// Overloaded indexing
bool& Assignment::operator[] (unsigned int i)
{
  assert(i < size);
  return binding[i];
}


// print --------------------------------------------------------------
std::ostream& Assignment::print(std::ostream & out) const
{
  out << "[";
  for (unsigned int i = 0; i < size; i++)
  {
    out << "<" << "x" << i << " |-> " << binding[i] << ">";
    if (i < size-1)
      out << ", ";
  }
  out << "]";
  return out;
}

std::ostream& operator<< (std::ostream & out, const Assignment &A)
{
  A.print(out);
  return(out);
}

// **********************************************************************
// Implementation of AssignmentIterator class
// **********************************************************************

// Constructor
AssignmentIterator::AssignmentIterator(const Assignment & A)
  :  binding(A.binding), size(A.size), current(0)
{
}

// Reset
//   reset the "current" item to be the first one
void AssignmentIterator::Reset()
{
  current = 0;
}

// Next
//   advance the "current" item
void AssignmentIterator::Next()
{
  if (current < size)
    current++;
}

// Advance
//   advance the "current" item by a given amount
void AssignmentIterator::Advance(unsigned int i)
{
  if (current + i < size + 1)
    current += i;
}

// MoreToDo
//   return true iff current item is in bounds
bool AssignmentIterator::MoreToDo() const
{
  return (current < size);
}

// Current
//   return reference to current item (error if out of bounds)
bool & AssignmentIterator::Current() const
{
    assert(current < size);
    return(binding[current]);
}
