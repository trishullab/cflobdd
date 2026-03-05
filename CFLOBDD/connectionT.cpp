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
	: entryPointHandle(new CFLOBDDNodeHandle(entryPointHandle)), returnMapHandle(returnMapHandle)
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
	delete entryPointHandle;
	entryPointHandle = NULL;
}

// Hash
template <typename Handle>
size_t ConnectionT<Handle>::Hash()
{
	size_t hvalue = 0;
	hvalue = (997 * returnMapHandle.Hash() + entryPointHandle->Hash());
	return hvalue;
}

// Overloaded =
template <typename Handle>
ConnectionT<Handle>& ConnectionT<Handle>::operator= (const ConnectionT<Handle>& C)
{
	if (this != &C)      // don't assign to self!
	{
		entryPointHandle = new CFLOBDDNodeHandle(*(C.entryPointHandle));
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
