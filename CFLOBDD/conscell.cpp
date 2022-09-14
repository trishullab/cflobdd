#ifndef _CONSCELL_CPP_GUARD
#define _CONSCELL_CPP_GUARD

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
#include "conscell.h"

template<typename T>
ConsCell<T>::ConsCell()
  :  next(NULL)
{
}

template<typename T>
ConsCell<T>::ConsCell(T i, ConsCell<T> *n)
  :  item(i), next(n)
{
}

template<typename T>
ConsCell<T>& ConsCell<T>::operator= (const ConsCell<T>& i)
{
  if (this != &i)      // don't assign to self!
  {
    item = i.item;
    next = i.next;
  }
  return *this;        
}

// Overloaded !=
template<typename T>
bool ConsCell<T>::operator!=(const ConsCell<T>& p)
{
  return !(*this == p);
}

// Overloaded ==
template<typename T>
bool ConsCell<T>::operator==(const ConsCell<T>& p)
{
  return (item == p.item) && (next == p.next);
}

// print
template<typename T>
std::ostream& ConsCell<T>::print(std::ostream & out) const
{
  out << Item() << " " << *Next();
  return out;
}

template<typename T>
std::ostream& operator<< (std::ostream & out, const ConsCell<T> &d)
{
  d.print(out);
  return(out);
}

#endif
