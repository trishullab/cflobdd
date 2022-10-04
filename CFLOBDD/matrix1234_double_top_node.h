#ifndef MATRIX1234_DOUBLE_TOP_NODE_GUARD
#define MATRIX1234_DOUBLE_TOP_NODE_GUARD

#include <iostream>
#include <fstream>
#include "matrix1234_double.h"
#include "return_map_T.h"

namespace CFL_OBDD {
	typedef ReturnMapBody<double> DoubleReturnMapBody;
	typedef ReturnMapHandle<double> DoubleReturnMapHandle;
	typedef ReturnMapBody<int> CFLOBDDReturnMapBody;
	typedef ReturnMapHandle<int> CFLOBDDReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<CFLOBDDReturnMapHandle> Connection;
}

namespace CFL_OBDD {

	typedef CFLOBDDTopNodeT<double> CFLOBDDTopNodeDouble;
	typedef CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeDoubleRefPtr;
	typedef CFLOBDDTopNodeT<int> CFLOBDDTopNode;
	typedef CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeIntRefPtr;

	namespace Matrix1234Double {

		// Initialization routine
		extern void Matrix1234InitializerTop();


		extern CFLOBDDTopNodeDoubleRefPtr MkIdRelationInterleavedTop(unsigned int i); // Representation of identity relation
		extern CFLOBDDTopNodeDoubleRefPtr MkWalshInterleavedTop(unsigned int i); // Representation of Walsh matrix
		extern CFLOBDDTopNodeDoubleRefPtr MkNegationMatrixInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeDoubleRefPtr MkCNOTInterleavedTop(unsigned int i);

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDDTopNodeDoubleRefPtr MkWalshVoc12Top(unsigned int i); // Representation of Walsh matrix with room for two additional vocabularies
		extern CFLOBDDTopNodeDoubleRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeDoubleRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeDoubleRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeDoubleRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeDoubleRefPtr MkDetensorConstraintInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeDoubleRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeDoubleRefPtr n); // Vocabulary projection
		extern CFLOBDDTopNodeDoubleRefPtr ReverseColumnsTop(CFLOBDDTopNodeDoubleRefPtr n);
		extern CFLOBDDTopNodeDoubleRefPtr MatrixTransposeTop(CFLOBDDTopNodeDoubleRefPtr n);

		extern CFLOBDDTopNodeDoubleRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeDoubleRefPtr c);
		extern CFLOBDDTopNodeDoubleRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeDoubleRefPtr c);

		extern CFLOBDDTopNodeDoubleRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeDoubleRefPtr c);
		extern CFLOBDDTopNodeDoubleRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeDoubleRefPtr c);
		extern CFLOBDDTopNodeDoubleRefPtr ConvertIntToDoubleTop(CFLOBDDTopNodeIntRefPtr c);

		extern void MatrixPrintRowMajorTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out);
		extern void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out);

		extern CFLOBDDTopNodeDoubleRefPtr SMatrixTop(std::string c);
		extern CFLOBDDTopNodeDoubleRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeDoubleRefPtr c1, CFLOBDDTopNodeDoubleRefPtr c2);


	}
}

#endif

