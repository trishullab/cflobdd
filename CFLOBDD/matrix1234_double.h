#ifndef MATRIX1234_DOUBLE_GUARD
#define MATRIX1234_DOUBLE_GUARD


#include <iostream>
#include <fstream>
#include "cflobdd_t.h"

namespace CFL_OBDD {

	typedef CFLOBDD_T<double> CFLOBDD_DOUBLE;
	typedef CFLOBDD_T<int> CFLOBDD;

	namespace Matrix1234Double {
		// Initialization routine
		extern void Matrix1234Initializer();

		extern CFLOBDD_DOUBLE MkIdRelationInterleaved(unsigned int i); // Representation of identity relation
		extern CFLOBDD_DOUBLE MkWalshInterleaved(unsigned int i);              // Representation of Walsh matrix
		extern CFLOBDD_DOUBLE MkNegationMatrixInterleaved(unsigned int i);
		extern CFLOBDD_DOUBLE MkCNOTInterleaved(unsigned int i);

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDD_DOUBLE MkWalshVoc12(unsigned int i);                    // Create representation of Walsh matrix with room for two extra vocabularies
		extern CFLOBDD_DOUBLE MatrixShiftVocs12To34(CFLOBDD_DOUBLE c);                // Vocabulary shift in a matrix
		extern CFLOBDD_DOUBLE MatrixShiftVoc42(CFLOBDD_DOUBLE c);                     // Vocabulary shift in a matrix
		extern CFLOBDD_DOUBLE KroneckerProduct(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2);        // Kronecker product (on matrices with room for two extra vocabularies)
		extern CFLOBDD_DOUBLE MkDetensorConstraintInterleaved(unsigned int i); // Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal: (W,X,Y,Z) s.t. X==Y with interleaved variables
		extern CFLOBDD_DOUBLE MatrixProjectVoc23(CFLOBDD_DOUBLE c);                   // Vocabulary projection
		extern CFLOBDD_DOUBLE MatrixDetensor(CFLOBDD_DOUBLE k);                       // Detensor of a 4-vocabulary matrix
		extern CFLOBDD_DOUBLE MatrixMultiply(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2);          // Matrix multiplication
		extern CFLOBDD_DOUBLE MatrixPadWithZeros(CFLOBDD_DOUBLE c, unsigned int level);
		extern CFLOBDD_DOUBLE ReverseColumns(CFLOBDD_DOUBLE c);
		extern CFLOBDD_DOUBLE MatrixTranspose(CFLOBDD_DOUBLE c);
		extern CFLOBDD_DOUBLE MatrixMultiplyV4(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2);          // Matrix multiplication

		extern CFLOBDD_DOUBLE MatrixShiftToAConnection(CFLOBDD_DOUBLE c);
		extern CFLOBDD_DOUBLE MatrixShiftToBConnection(CFLOBDD_DOUBLE c);
		extern CFLOBDD_DOUBLE KroneckerProduct2Vocs(CFLOBDD_DOUBLE m1, CFLOBDD_DOUBLE m2);

		extern CFLOBDD_DOUBLE PromoteInterleavedTo12(CFLOBDD_DOUBLE c);
		extern CFLOBDD_DOUBLE Demote12ToInterleaved(CFLOBDD_DOUBLE c);
		extern CFLOBDD_DOUBLE ConvertIntToDouble(CFLOBDD c);

		extern void MatrixPrintRowMajor(CFLOBDD_DOUBLE c, std::ostream & out);
		extern void MatrixPrintRowMajorInterleaved(CFLOBDD_DOUBLE c, std::ostream & out);

		// Simon's Algo Test Matrix Construction
		extern CFLOBDD_DOUBLE SMatrix(std::string s);
		extern CFLOBDD_DOUBLE Func2To1CFLOBDDMatrix(std::string s, std::string t);

		extern CFLOBDD_DOUBLE NormalizeOutputTo1(CFLOBDD_DOUBLE c);
	}
}

#endif

