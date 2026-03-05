#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "vector_node.h"
#include "vector_int_top_node.h"

namespace CFL_OBDD {

	namespace VectorInt {

		void VectorInitializerTop()
		{
			VectorInitializerNode();
			return;
		}

		CFLOBDDTopNodeIntRefPtr MkBasisVectorTop(unsigned int level, unsigned int index)
		{
			CFLOBDDTopNodeIntRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle rhandle;

			tempHandle = MkBasisVectorNode(level, index);
			if (index == 0)
			{
				rhandle.AddToEnd(1);
				rhandle.AddToEnd(0);
			}
			else
			{
				rhandle.AddToEnd(0);
				rhandle.AddToEnd(1);
			}
			rhandle.Canonicalize();

			ptr = new CFLOBDDTopNode(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeIntRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeIntRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorToMatrixInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeIntRefPtr MatrixToVectorTop(CFLOBDDTopNodeIntRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixToVectorNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeIntRefPtr NoDistinctionNodeTop(unsigned int level)
		{
			CFLOBDDReturnMapHandle m1;
			m1.AddToEnd(1);
			m1.Canonicalize();

			return new CFLOBDDTopNode(CFLOBDDNodeHandle::NoDistinctionNode[level],m1);
		}

		CFLOBDDTopNodeIntRefPtr MkColumn1MatrixTop(unsigned int level)
		{
			CFLOBDDTopNodeIntRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			CFLOBDDReturnMapHandle rhandle;

			tempHandle = MkColumn1MatrixNode(level);
			
			rhandle.AddToEnd(1);
			rhandle.AddToEnd(0);
			rhandle.Canonicalize();

			ptr = new CFLOBDDTopNode(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeIntRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeIntRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;
			CFLOBDDNodeHandle tempHandle = MkVectorWithVoc12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeIntRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeIntRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 1);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorShiftVocs1To2Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeIntRefPtr v = new CFLOBDDTopNode(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}
	}
}

