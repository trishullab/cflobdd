#ifndef _WEIGHTED_CONNECTION_T_H
#define _WEIGHTED_CONNECTION_T_H

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

	template <typename T, typename Op>
	class WConnection
	{
	public:
		WConnection();                                  // Default constructor
		WConnection(WeightedCFLOBDDNode<T, Op> *entryPoint, CFLOBDDReturnMapHandle &returnMapHandle);
		WConnection(WeightedCFLOBDDNodeHandleT<T, Op> &entryPointHandle, CFLOBDDReturnMapHandle &returnMapHandle);
		~WConnection();                                 // Destructor

		unsigned int Hash(unsigned int modsize);
		WConnection& operator= (const WConnection &C);   // Overloaded =
		bool operator!= (const WConnection & C);        // Overloaded !=
		bool operator== (const WConnection & C);        // Overloaded ==

		WeightedCFLOBDDNodeHandleT<T, Op>* entryPointHandle;
		CFLOBDDReturnMapHandle returnMapHandle;

	public:
		std::ostream& print(std::ostream & out = std::cout) const;
	};


	template <typename T, typename Op>
	std::ostream& operator<< (std::ostream & out, const WConnection<T, Op> &c);

}

#endif

