#ifndef _WEIGHTED_ROOT_CONNECTION_T_H
#define _WEIGHTED_ROOT_CONNECTION_T_H

#include <cassert>
#include <cstdlib>
#include <iostream>
#include "weighted_cflobdd_node_t.h"

namespace CFL_OBDD {

	
	template <typename, typename>
	class WeightedCFLOBDDNode;
	template <typename, typename>
	class WeightedCFLOBDDInternalNode;   //  : public CFLOBDDNode
	template <typename, typename>
	class WeightedCFLOBDDLeafNode;       //  : public CFLOBDDNode
	template <typename, typename>
	class WeightedCFLOBDDForkNode;       //  : public CFLOBDDLeafNode
	template <typename, typename>
	class WeightedCFLOBDDDontCareNode;   //  : public CFLOBDDLeafNode
	template <typename, typename>
	class WeightedCFLOBDDNodeHandleT;
}


namespace CFL_OBDD {

	template <typename Handle, typename T, typename Op>
	class WRootConnection
	{
	public:
		WRootConnection();                                  // Default constructor
		WRootConnection(WeightedCFLOBDDNode<T, Op> *entryPoint, Handle &returnMapHandle, T factor = 1);
		WRootConnection(WeightedCFLOBDDNodeHandleT<T, Op> &entryPointHandle, Handle &returnMapHandle, T factor = 1);
		~WRootConnection();                                 // Destructor

		unsigned int Hash(unsigned int modsize);
		WRootConnection& operator= (const WRootConnection &C);   // Overloaded =
		bool operator!= (const WRootConnection & C);        // Overloaded !=
		bool operator== (const WRootConnection & C);        // Overloaded ==

		WeightedCFLOBDDNodeHandleT<T, Op>* entryPointHandle;
		Handle returnMapHandle;
		T factor;

	public:
		std::ostream& print(std::ostream & out = std::cout) const;
	};


	template <typename Handle, typename T, typename Op>
	std::ostream& operator<< (std::ostream & out, const WRootConnection<Handle, T, Op> &c);

}

#endif

