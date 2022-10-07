#include "cflobdd_node.h"
#include "connectionT.h"
// #include "return_map_specializations.h"

namespace CFL_OBDD {


// Default constructor
template <typename Handle>
ConnectionT<Handle>::ConnectionT()
{
}

// Constructor
template <typename Handle>
ConnectionT<Handle>::ConnectionT(CFLOBDDNodeHandle &entryPointHandle, Handle &returnMapHandle)
	: entryPointHandle(new CFLOBDDNodeHandle(entryPointHandle.handleContents)), returnMapHandle(returnMapHandle)
{
}

// Constructor
template <typename Handle>
ConnectionT<Handle>::ConnectionT(CFLOBDDNode *entryPoint, Handle &returnMapHandle)
	: entryPointHandle(new CFLOBDDNodeHandle(entryPoint)), returnMapHandle(returnMapHandle)
{
}

template <typename Handle>
ConnectionT<Handle>::~ConnectionT()
{
}

// Hash
template <typename Handle>
unsigned int ConnectionT<Handle>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;
	hvalue = (997 * returnMapHandle.Hash(modsize) + entryPointHandle->Hash(modsize)) % modsize;
	return hvalue;
}

// Overloaded =
template <typename Handle>
ConnectionT<Handle>& ConnectionT<Handle>::operator= (const ConnectionT<Handle>& C)
{
	if (this != &C)      // don't assign to self!
	{
		entryPointHandle = C.entryPointHandle;
		returnMapHandle = C.returnMapHandle;
	}
	return *this;
}

// Overloaded !=
template <typename Handle>
bool ConnectionT<Handle>::operator!= (const ConnectionT<Handle> & C)
{
	return (returnMapHandle != C.returnMapHandle) || ((*entryPointHandle) != (*C.entryPointHandle));
}

// Overloaded ==
template <typename Handle>
bool ConnectionT<Handle>::operator== (const ConnectionT<Handle> & C)
{
	return (returnMapHandle == C.returnMapHandle) && ((*entryPointHandle) == (*C.entryPointHandle));
}

// print
template <typename Handle>
std::ostream& ConnectionT<Handle>::print(std::ostream & out) const
{
	out << (*entryPointHandle);
	out << returnMapHandle;
	return out;
}

template <typename Handle>	
std::ostream& operator<< (std::ostream & out, const ConnectionT<Handle> &c)
	{
		c.print(out);
		return(out);
	}
} // namespace CFL_OBDD
