#ifndef VECTOR_COMPLEX_FLOAT_BOOST_GUARD
#define VECTOR_COMPLEX_FLOAT_BOOST_GUARD


#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <boost/multiprecision/cpp_complex.hpp>
#include "cflobdd_t.h"

namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef mp::cpp_complex_100 BIG_COMPLEX_FLOAT;
	typedef CFLOBDD_T<BIG_COMPLEX_FLOAT> CFLOBDD_COMPLEX_BIG;
	//typedef CFLOBDD_T<BIG_FLOAT> CFLOBDD_COMPLEX_BIG;

	namespace VectorComplexFloatBoost {
		// Initialization routine
		extern void VectorInitializer();

		extern CFLOBDD_COMPLEX_BIG MkBasisVector(unsigned int level, unsigned int index);         // Representation of basis vector with index i
		extern CFLOBDD_COMPLEX_BIG MkBasisVector(unsigned int level, std::string s); // index i represented as string
		extern CFLOBDD_COMPLEX_BIG VectorToMatrixInterleaved(CFLOBDD_COMPLEX_BIG c); // Convert vector to matrix with variables interleaved
		//extern CFLOBDD_COMPLEX_BIG MatrixToVector(CFLOBDD_COMPLEX_BIG c); // Convert vector to matrix with variables interleaved
		extern CFLOBDD_COMPLEX_BIG MkColumn1Matrix(unsigned int level);
		extern CFLOBDD_COMPLEX_BIG MkVectorWithVoc12(CFLOBDD_COMPLEX_BIG c); // convert the vector into another vector with extra volcabularies
		extern CFLOBDD_COMPLEX_BIG KroneckerProduct(CFLOBDD_COMPLEX_BIG m1, CFLOBDD_COMPLEX_BIG m2);
		extern CFLOBDD_COMPLEX_BIG VectorShiftVocs1To2(CFLOBDD_COMPLEX_BIG m1);
		extern CFLOBDD_COMPLEX_BIG VectorPadWithZeros(CFLOBDD_COMPLEX_BIG c, unsigned int level);
		extern CFLOBDD_COMPLEX_BIG NoDistinctionNode(unsigned int level, BIG_COMPLEX_FLOAT val);
		//#ifdef PATH_COUNTING_ENABLED
		extern std::string Sampling(CFLOBDD_COMPLEX_BIG c, bool voctwo);
		extern std::string SamplingV2(CFLOBDD_COMPLEX_BIG c);
		//#endif
		extern CFLOBDD_COMPLEX_BIG VectorWithAmplitude(CFLOBDD_COMPLEX_BIG c);
		extern void VectorPrintColumnMajor(CFLOBDD_COMPLEX_BIG c, std::ostream & out);
		extern void VectorPrintColumnMajorInterleaved(CFLOBDD_COMPLEX_BIG c, std::ostream & out);
		extern long double getNonZeroProbability(CFLOBDD_COMPLEX_BIG n);
		extern unsigned long long int GetPathCount(CFLOBDD_COMPLEX_BIG n, long double p);
	}
}

#endif

