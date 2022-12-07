
#include "weighted_cflobdd_node_t.h"
#include "weighted_rootConnectionT.h"
// #include "return_map_specializations.h"

namespace CFL_OBDD {


// Default constructor
template <typename Handle, typename T, typename Op>
WRootConnection<Handle, T, Op>::WRootConnection()
{
}

// Constructor
template <typename Handle, typename T, typename Op>
WRootConnection<Handle, T, Op>::WRootConnection(WeightedCFLOBDDNodeHandleT<T,Op> &entryPointHandle, Handle &returnMapHandle, T factor)
	: entryPointHandle(new WeightedCFLOBDDNodeHandleT<T,Op>(entryPointHandle.handleContents)), returnMapHandle(returnMapHandle), factor(factor)
{
}

// Constructor
template <typename Handle, typename T, typename Op>
WRootConnection<Handle, T, Op>::WRootConnection(WeightedCFLOBDDNode<T,Op> *entryPoint, Handle &returnMapHandle, T factor)
	: entryPointHandle(new WeightedCFLOBDDNodeHandleT<T,Op>(entryPoint)), returnMapHandle(returnMapHandle), factor(factor)
{
}

template <typename Handle, typename T, typename Op>
WRootConnection<Handle, T, Op>::~WRootConnection()
{
}

// Hash
template <typename Handle, typename T, typename Op>
unsigned int WRootConnection<Handle, T, Op>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;
    boost::hash<T> boost_hash;
	hvalue = (997 * returnMapHandle.Hash(modsize) + 97 * entryPointHandle->Hash(modsize) + boost_hash(factor) % modsize) % modsize;
	return hvalue;
}

// Overloaded =
template <typename Handle, typename T, typename Op>
WRootConnection<Handle, T, Op>& WRootConnection<Handle, T, Op>::operator= (const WRootConnection<Handle, T, Op>& C)
{
	if (this != &C)      // don't assign to self!
	{
		entryPointHandle = C.entryPointHandle;
		returnMapHandle = C.returnMapHandle;
        factor = C.factor;
	}
	return *this;
}

// Overloaded !=
template <typename Handle, typename T, typename Op>
bool WRootConnection<Handle, T, Op>::operator!= (const WRootConnection<Handle, T, Op> & C)
{
	return (returnMapHandle != C.returnMapHandle) || ((*entryPointHandle) != (*C.entryPointHandle)) || (factor != C.factor);
}

// Overloaded ==
template <typename Handle, typename T, typename Op>
bool WRootConnection<Handle, T, Op>::operator== (const WRootConnection<Handle, T, Op> & C)
{
	return (returnMapHandle == C.returnMapHandle) && ((*entryPointHandle) == (*C.entryPointHandle)) && (factor == C.factor);
}

// print
template <typename Handle, typename T, typename Op>
std::ostream& WRootConnection<Handle, T, Op>::print(std::ostream & out) const
{
	out << (*entryPointHandle);
	out << returnMapHandle;
    out << factor;
	return out;
}

template <typename Handle, typename T, typename Op>	
std::ostream& operator<< (std::ostream & out, const WRootConnection<Handle, T, Op> &c)
	{
		c.print(out);
		return(out);
	}
} // namespace CFL_OBDD
