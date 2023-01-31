#ifndef W_MATRIX1234_COMPLEX_FB_MUL_GUARD
#define W_MATRIX1234_COMPLEX_FB_MUL_GUARD


#include <iostream>
#include <fstream>
#include <random>
#include <boost/multiprecision/cpp_complex.hpp>
#include "weighted_cflobdd_t.h"


namespace CFL_OBDD {


	typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
    // typedef double BIG_COMPLEX_FLOAT;
	typedef WEIGHTED_CFLOBDD_T<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL;
	// cflobdd_kind = 1 -> CFLOBDD, 0 -> BDD
	namespace WeightedMatrix1234ComplexFloatBoostMul {
		// Initialization routine
		extern void Matrix1234Initializer();

		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkIdRelationInterleaved(unsigned int i, int cflobdd_kind = 1); // Representation of identity relation
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkWalshInterleaved(unsigned int i, int cflobdd_kind = 1);              // Representation of Walsh matrix
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkNegationMatrixInterleaved(unsigned int i, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkPauliYGate(unsigned int i, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkPauliZGate(unsigned int i, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkSGate(unsigned int i, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCNOTInterleaved(unsigned int i);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkExchangeInterleaved(unsigned int i); // Representation of exchange matrix
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCNOT(unsigned int level, unsigned int n, long int controller, long int controlled, int cflobdd_kind = 1); // Representation of CNOT matrix with index1 as controller and index2 as controlled bits
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCCNOT(unsigned int level, long int controller1, long int controller2, long int controlled, int cflobdd_kind = 1); // Representation of CNOT matrix with index1 as controller and index2 as controlled bits
		
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkRestrictMatrix(unsigned int level, std::string s, int cflobdd_kind = 1); 

		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL KroneckerProduct2Vocs(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m1, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m2);

		extern void MatrixPrintRowMajor(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c, std::ostream & out);
		extern void MatrixPrintRowMajorInterleaved(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL c, std::ostream & out);

		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MatrixMultiplyV4(WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m1, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL m2);
        extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL CreateBalancedFn(int n, std::mt19937 mt);

		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkSwapGate(unsigned int i, long int c1, long int c2, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkiSwapGate(unsigned int i, long int c1, long int c2, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCPGate(unsigned int i, long int c1, long int c2, double theta, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkPhaseShiftGate(unsigned int i, double theta, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCZGate(unsigned int i, long int c1, long int c2, double theta, int cflobdd_kind = 1);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL MkCSwapGate(unsigned int i, long int c1, long int x1, long int x2, int cflobdd_kind = 1);
	}
}

#endif

