
#include "weighted_cflobdd_node_t.h"
#include "weighted_connectionT.h"
// #include "return_map_specializations.h"

namespace CFL_OBDD {


// Default constructor
template <typename T, typename Op>
WConnection<T, Op>::WConnection()
{
}

// Constructor
template <typename T, typename Op>
WConnection<T, Op>::WConnection(WeightedCFLOBDDNodeHandleT<T,Op> &entryPointHandle, CFLOBDDReturnMapHandle &returnMapHandle)
	: entryPointHandle(new WeightedCFLOBDDNodeHandleT<T,Op>(entryPointHandle.handleContents)), returnMapHandle(returnMapHandle)
{
}

// Constructor
template <typename T, typename Op>
WConnection<T, Op>::WConnection(WeightedCFLOBDDNode<T,Op> *entryPoint, CFLOBDDReturnMapHandle &returnMapHandle)
	: entryPointHandle(new WeightedCFLOBDDNodeHandleT<T,Op>(entryPoint)), returnMapHandle(returnMapHandle)
{
}

template <typename T, typename Op>
WConnection<T, Op>::~WConnection()
{
}

// Hash
template <typename T, typename Op>
unsigned int WConnection<T, Op>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;
	hvalue = (997 * returnMapHandle.Hash(modsize) + entryPointHandle->Hash(modsize)) % modsize;
	return hvalue;
}

// Overloaded =
template <typename T, typename Op>
WConnection<T, Op>& WConnection<T, Op>::operator= (const WConnection<T, Op>& C)
{
	if (this != &C)      // don't assign to self!
	{
		entryPointHandle = C.entryPointHandle;
		returnMapHandle = C.returnMapHandle;
	}
	return *this;
}

// Overloaded !=
template <typename T, typename Op>
bool WConnection<T, Op>::operator!= (const WConnection<T, Op> & C)
{
	return (returnMapHandle != C.returnMapHandle) || ((*entryPointHandle) != (*C.entryPointHandle));
}

// Overloaded ==
template <typename T, typename Op>
bool WConnection<T, Op>::operator== (const WConnection<T, Op> & C)
{
	return (returnMapHandle == C.returnMapHandle) && ((*entryPointHandle) == (*C.entryPointHandle));
}

// print
template <typename T, typename Op>
std::ostream& WConnection<T, Op>::print(std::ostream & out) const
{
	out << (*entryPointHandle);
	out << returnMapHandle;
	return out;
}

template <typename T, typename Op>	
std::ostream& operator<< (std::ostream & out, const WConnection<T, Op> &c)
	{
		c.print(out);
		return(out);
	}
} // namespace CFL_OBDD
