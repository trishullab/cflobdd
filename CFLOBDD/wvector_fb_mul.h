#ifndef W_VECTOR_FB_MUL_GUARD
#define W_VECTOR_FB_MUL_GUARD


#include <iostream>
#include <fstream>
#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "weighted_cflobdd_t.h"


namespace CFL_OBDD {


	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
    // typedef double BIG_FLOAT;
	typedef WEIGHTED_CFLOBDD_T<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL;

	namespace WeightedVectorFloatBoostMul {
		// Initialization routine
		extern void VectorInitializer();

		extern WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL MkBasisVector(unsigned int level, unsigned int index);         // Representation of basis vector with index i
		extern WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL MkBasisVector(unsigned int level, std::string s); // index i represented as string
		extern WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL VectorToMatrixInterleaved(WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL c); // Convert vector to matrix with variables interleaved
		extern WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL MkColumn1Matrix(unsigned int level);
		extern WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL KroneckerProduct(WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL m1, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL m2);
		extern WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL NoDistinctionNode(unsigned int level, BIG_FLOAT val);
		extern std::string Sampling(WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL c, bool voctwo, std::string = "");
		extern WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL VectorWithAmplitude(WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL c);
		extern void VectorPrintColumnMajor(WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL c, std::ostream & out);
		extern void VectorPrintColumnMajorInterleaved(WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL c, std::ostream & out);
	}
}

#endif

