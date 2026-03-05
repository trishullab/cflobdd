#ifndef MATRIX1234_FLOAT_BOOST_TOP_NODE_GUARD
#define MATRIX1234_FLOAT_BOOST_TOP_NODE_GUARD

#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "matrix1234_float_boost.h"
#include "return_map_T.h"
// #include "general_map.h"

 namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_FLOAT;
	typedef ReturnMapBody<BIG_FLOAT> FloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_FLOAT> FloatBoostReturnMapHandle;
	typedef ReturnMapBody<int> CFLOBDDReturnMapBody;
	typedef ReturnMapHandle<int> CFLOBDDReturnMapHandle;
	// typedef ReturnMapHandle<GeneralMapHandle> GeneralReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<CFLOBDDReturnMapHandle> Connection;
}

namespace CFL_OBDD {

	typedef CFLOBDDTopNodeT<BIG_FLOAT> CFLOBDDTopNodeFloatBoost;
	typedef CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeFloatBoostRefPtr;
	typedef CFLOBDDTopNodeT<int> CFLOBDDTopNode;
	typedef CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeIntRefPtr;
	// typedef CFLOBDDTopNodeT<GeneralMapHandle> CFLOBDDTopNodeGeneralMap;
	// typedef CFLOBDDTopNodeT<GeneralMapHandle>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeGeneralMapRefPtr;
	// typedef CFLOBDD_T<GeneralMapHandle> CFLOBDD_GENERAL;
	// typedef Hashtable<CFLOBDDNodeHandle, CFLOBDD_GENERAL> CFLOBDDGeneralMapNodeMemoTable;
	// typedef ref_ptr<CFLOBDDGeneralMapNodeMemoTable> CFLOBDDGeneralMapNodeMemoTableRefPtr;


	namespace Matrix1234FloatBoost {

		// Initialization routine
		extern void Matrix1234InitializerTop();

		extern CFLOBDDTopNodeFloatBoostRefPtr MkIdRelationInterleavedTop(unsigned int i); // Representation of identity relation
		extern CFLOBDDTopNodeFloatBoostRefPtr MkWalshInterleavedTop(unsigned int i); // Representation of Walsh matrix
		extern CFLOBDDTopNodeFloatBoostRefPtr MkNegationMatrixInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeFloatBoostRefPtr MkCNOTInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeFloatBoostRefPtr MkExchangeInterleavedTop(unsigned int i); // Representation of exchange matrix
		extern CFLOBDDTopNodeFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled);

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDDTopNodeFloatBoostRefPtr MkWalshVoc12Top(unsigned int i); // Representation of Walsh matrix with room for two additional vocabularies
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftVocs12To34Top(CFLOBDDTopNodeFloatBoostRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftVoc42Top(CFLOBDDTopNodeFloatBoostRefPtr n); // Vocabulary shift
		extern CFLOBDDTopNodeFloatBoostRefPtr MkDetensorConstraintInterleavedTop(unsigned int i);
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixProjectVoc23Top(CFLOBDDTopNodeFloatBoostRefPtr n); // Vocabulary projection
		extern CFLOBDDTopNodeFloatBoostRefPtr ReverseColumnsTop(CFLOBDDTopNodeFloatBoostRefPtr n);
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixTransposeTop(CFLOBDDTopNodeFloatBoostRefPtr n);

		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftToAConnectionTop(CFLOBDDTopNodeFloatBoostRefPtr c);
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixShiftToBConnectionTop(CFLOBDDTopNodeFloatBoostRefPtr c);

		extern CFLOBDDTopNodeFloatBoostRefPtr PromoteInterleavedTo12Top(CFLOBDDTopNodeFloatBoostRefPtr c);
		extern CFLOBDDTopNodeFloatBoostRefPtr Demote12ToInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr c);
		extern CFLOBDDTopNodeFloatBoostRefPtr ConvertIntToFloatBoostTop(CFLOBDDTopNodeIntRefPtr c);
		extern CFLOBDDTopNodeFloatBoostRefPtr PromoteInterleavedTo13Top(CFLOBDDTopNodeFloatBoostRefPtr c);
		extern CFLOBDDTopNodeFloatBoostRefPtr AddMatrixRowsTopNode(CFLOBDDTopNodeFloatBoostRefPtr c);

		extern void MatrixPrintRowMajorTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);
		extern void MatrixPrintRowMajorInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);

		extern CFLOBDDTopNodeFloatBoostRefPtr SMatrixTop(std::string c);
		extern CFLOBDDTopNodeFloatBoostRefPtr MkDistinctionTwoVarsTopNode(int x, int y, unsigned int var_level, unsigned int total_level);
		extern CFLOBDDTopNodeFloatBoostRefPtr MultiplyOperationTopNode(CFLOBDDTopNodeFloatBoostRefPtr c);
		extern CFLOBDDTopNodeFloatBoostRefPtr MkMatrixMultiplyConstraintTopNode(unsigned int level);
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixMultiplyV4TopNode(CFLOBDDTopNodeFloatBoostRefPtr c1, CFLOBDDTopNodeFloatBoostRefPtr c2);
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixMultiplyV4WithInfoTopNode(CFLOBDDTopNodeFloatBoostRefPtr c1, CFLOBDDTopNodeFloatBoostRefPtr c2);
		extern CFLOBDDTopNodeFloatBoostRefPtr ComputeShortestPathTop(CFLOBDDTopNodeFloatBoostRefPtr c1, CFLOBDDTopNodeFloatBoostRefPtr c2);
	}
}

#endif

