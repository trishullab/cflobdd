#ifndef _CONNECTION_T_H
#define _CONNECTION_T_H

//
//    Copyright (c) 1999, 2017 Thomas W. Reps
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
#include <cstdlib>
#include <iostream>
#include "cflobdd_node.h"

namespace CFL_OBDD {

	
	class CFLOBDDNode;
	class CFLOBDDInternalNode;   //  : public CFLOBDDNode
	class CFLOBDDLeafNode;       //  : public CFLOBDDNode
	class CFLOBDDForkNode;       //  : public CFLOBDDLeafNode
	class CFLOBDDDontCareNode;   //  : public CFLOBDDLeafNode
	class CFLOBDDNodeHandle;
}


namespace CFL_OBDD {

	template <typename Handle>
	class ConnectionT
	{
	public:
		ConnectionT();                                  // Default constructor
		ConnectionT(CFLOBDDNode *entryPoint, Handle &returnMapHandle);
		ConnectionT(CFLOBDDNodeHandle &entryPointHandle, Handle &returnMapHandle);
		~ConnectionT();                                 // Destructor

		unsigned int Hash(unsigned int modsize);
		ConnectionT& operator= (const ConnectionT &C);   // Overloaded =
		bool operator!= (const ConnectionT & C);        // Overloaded !=
		bool operator== (const ConnectionT & C);        // Overloaded ==

		CFLOBDDNodeHandle* entryPointHandle;
		Handle returnMapHandle;

	public:
		std::ostream& print(std::ostream & out = std::cout) const;
	};


	template <typename Handle>
	std::ostream& operator<< (std::ostream & out, const ConnectionT<Handle> &c);

}

#endif

