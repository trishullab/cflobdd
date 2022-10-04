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
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "inttriple.h"

inttriple::inttriple()
  :  first(0), second(0), third(0)
{
}

inttriple::inttriple(const int i1, const int i2, const int i3)
  :  first(i1), second(i2), third(i3)
{
}

std::ostream& operator<< (std::ostream & out, const inttriple &p)
{
  out << "(" << p.First() << ", " << p.Second() << ", " << p.Third() << ")";
  return(out);
}

inttriple& inttriple::operator= (const inttriple& i)
{
  if (this != &i)      // don't assign to self!
  {
    first = i.first;
    second = i.second;
    third = i.third;
  }
  return *this;        
}

// Overloaded !=
bool inttriple::operator!=(const inttriple& p)
{
  return (first != p.first) || (second != p.second) || (third != p.third);
}

// Overloaded ==
bool inttriple::operator==(const inttriple& p)
{
  return (first == p.first) && (second == p.second) && (third == p.third);
}

