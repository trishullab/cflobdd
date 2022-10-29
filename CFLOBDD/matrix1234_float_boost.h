#ifndef MATRIX1234_FLOAT_BOOST_GUARD
#define MATRIX1234_FLOAT_BOOST_GUARD


#include <iostream>
#include <fstream>
#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "cflobdd_t.h"


namespace CFL_OBDD {


	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
	//typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000> > BIG_FLOAT;
	typedef CFLOBDD_T<BIG_FLOAT> CFLOBDD_FLOAT_BOOST;
	typedef CFLOBDD_T<int> CFLOBDD;

	namespace Matrix1234FloatBoost {
		// Initialization routine
		extern void Matrix1234Initializer();

		extern CFLOBDD_FLOAT_BOOST MkIdRelationInterleaved(unsigned int i); // Representation of identity relation
		extern CFLOBDD_FLOAT_BOOST MkWalshInterleaved(unsigned int i);              // Representation of Walsh matrix
		extern CFLOBDD_FLOAT_BOOST MkNegationMatrixInterleaved(unsigned int i);
		extern CFLOBDD_FLOAT_BOOST MkCNOTInterleaved(unsigned int i);
		extern CFLOBDD_FLOAT_BOOST MkExchangeInterleaved(unsigned int i); // Representation of exchange matrix
		extern CFLOBDD_FLOAT_BOOST MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled); // Representation of CNOT matrix with index1 as controller and index2 as controlled bits

		// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
		extern CFLOBDD_FLOAT_BOOST MkWalshVoc12(unsigned int i);                    // Create representation of Walsh matrix with room for two extra vocabularies
		extern CFLOBDD_FLOAT_BOOST MatrixShiftVocs12To34(CFLOBDD_FLOAT_BOOST c);                // Vocabulary shift in a matrix
		extern CFLOBDD_FLOAT_BOOST MatrixShiftVoc42(CFLOBDD_FLOAT_BOOST c);                     // Vocabulary shift in a matrix
		extern CFLOBDD_FLOAT_BOOST KroneckerProduct(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);        // Kronecker product (on matrices with room for two extra vocabularies)
		extern CFLOBDD_FLOAT_BOOST MkDetensorConstraintInterleaved(unsigned int i); // Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal: (W,X,Y,Z) s.t. X==Y with interleaved variables
		extern CFLOBDD_FLOAT_BOOST MatrixProjectVoc23(CFLOBDD_FLOAT_BOOST c);                   // Vocabulary projection
		extern CFLOBDD_FLOAT_BOOST MatrixDetensor(CFLOBDD_FLOAT_BOOST k);                       // Detensor of a 4-vocabulary matrix
		extern CFLOBDD_FLOAT_BOOST MatrixMultiply(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);          // Matrix multiplication
		extern CFLOBDD_FLOAT_BOOST MatrixPadWithZeros(CFLOBDD_FLOAT_BOOST c, unsigned int level);
		extern CFLOBDD_FLOAT_BOOST ReverseColumns(CFLOBDD_FLOAT_BOOST c);
		extern CFLOBDD_FLOAT_BOOST MatrixTranspose(CFLOBDD_FLOAT_BOOST c);

		extern CFLOBDD_FLOAT_BOOST MatrixShiftToAConnection(CFLOBDD_FLOAT_BOOST c);
		extern CFLOBDD_FLOAT_BOOST MatrixShiftToBConnection(CFLOBDD_FLOAT_BOOST c);
		extern CFLOBDD_FLOAT_BOOST KroneckerProduct2Vocs(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);

		extern CFLOBDD_FLOAT_BOOST PromoteInterleavedTo12(CFLOBDD_FLOAT_BOOST c);
		extern CFLOBDD_FLOAT_BOOST Demote12ToInterleaved(CFLOBDD_FLOAT_BOOST c);
		extern CFLOBDD_FLOAT_BOOST ConvertIntToFloatBoost(CFLOBDD c);
		extern CFLOBDD_FLOAT_BOOST PromoteInterleavedTo13(CFLOBDD_FLOAT_BOOST c);
		extern CFLOBDD_FLOAT_BOOST AddMatrixRows(CFLOBDD_FLOAT_BOOST m);

		extern void MatrixPrintRowMajor(CFLOBDD_FLOAT_BOOST c, std::ostream & out);
		extern void MatrixPrintRowMajorInterleaved(CFLOBDD_FLOAT_BOOST c, std::ostream & out);

		// Simon's Algo Test Matrix Construction
		extern CFLOBDD_FLOAT_BOOST SMatrix(std::string s);
		extern CFLOBDD_FLOAT_BOOST Func2To1CFLOBDDMatrix_DivideAndConquer(std::string s, std::mt19937 mt);
		extern CFLOBDD_FLOAT_BOOST Func2To1CFLOBDDMatrix(std::string s, std::vector<int> v);
		extern CFLOBDD_FLOAT_BOOST NormalizeOutputTo1(CFLOBDD_FLOAT_BOOST c);
		CFLOBDD_FLOAT_BOOST MultiplyOperation(CFLOBDD_FLOAT_BOOST m);
		extern CFLOBDD_FLOAT_BOOST MatrixMultiplyV3(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);
		extern CFLOBDD_FLOAT_BOOST MatrixMultiplyV4(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);
		extern CFLOBDD_FLOAT_BOOST MatrixMultiplyV4WithInfo(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);
		extern CFLOBDD_FLOAT_BOOST MkMatrixMultiplyConstraint(unsigned int level);
		extern CFLOBDD_FLOAT_BOOST CreateBalancedFn(int n, std::mt19937 mt);
		extern CFLOBDD_FLOAT_BOOST ApplyExchangeAndIdentity(std::string s);
		extern CFLOBDD_FLOAT_BOOST ComputeShortestPath(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);
	}
}

#endif

