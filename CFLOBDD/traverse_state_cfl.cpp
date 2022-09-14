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
#include "cflobdd_node.h"
#include "traverse_state_cfl.h"

using namespace CFL_OBDD;

//**********************************************************************
// TraverseState
//**********************************************************************

TraverseState::TraverseState()
  :  node(NULL), visitState(FirstVisit)
{
}

TraverseState::TraverseState(CFLOBDDNode *n, VisitState vs)
  :  node(n), visitState(vs)
{
}

TraverseState::TraverseState(CFLOBDDNode *n, VisitState vs, int i)
  :  node(n), visitState(vs), index(i)
{
}

std::ostream& operator<< (std::ostream & out, const TraverseState &p)
{
  p.print(out);
  return(out);
}

TraverseState& TraverseState::operator= (const TraverseState& i)
{
  if (this != &i)      // don't assign to self!
  {
    node = i.node;
    visitState = i.visitState;
    if (i.visitState == ThirdVisit)
      index = i.index;
  }
  return *this;        
}

// Overloaded !=
bool TraverseState::operator!=(const TraverseState& p)
{
  return !(*this == p);
}

// Overloaded ==
bool TraverseState::operator==(const TraverseState& p)
{
  if (visitState == ThirdVisit) {
    return (node == p.node) && (visitState == p.visitState) && (index == p.index);
  }
  return (node == p.node) && (visitState == p.visitState);
}

// print
std::ostream& TraverseState::print(std::ostream & out) const
{
  if (visitState == ThirdVisit) {
    out << "(" << node << ", " << visitState << ", " << index << ")";
  }
  else {
    out << "(" << node << ", " << visitState << ")";
  } 
  return out;
}
