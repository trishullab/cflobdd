#ifndef VECTOR_FLOAT_BOOST_GUARD
#define VECTOR_FLOAT_BOOST_GUARD


#include <iostream>
#include <string>
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "cflobdd_t.h"

 namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
	//typedef boost::multiprecision::number<mp::cpp_dec_float<1000> > BIG_FLOAT;
	typedef CFLOBDD_T<BIG_FLOAT> CFLOBDD_FLOAT_BOOST;

	namespace VectorFloatBoost {
		// Initialization routine
		extern void VectorInitializer();

		extern CFLOBDD_FLOAT_BOOST MkBasisVector(unsigned int level, unsigned int index);         // Representation of basis vector with index i
		extern CFLOBDD_FLOAT_BOOST MkBasisVector(unsigned int level, std::string s); // index i represented as string
		extern CFLOBDD_FLOAT_BOOST VectorToMatrixInterleaved(CFLOBDD_FLOAT_BOOST c); // Convert vector to matrix with variables interleaved
		extern CFLOBDD_FLOAT_BOOST MatrixToVector(CFLOBDD_FLOAT_BOOST c); // Convert vector to matrix with variables interleaved
		extern CFLOBDD_FLOAT_BOOST MkColumn1Matrix(unsigned int level);
		extern CFLOBDD_FLOAT_BOOST MkVectorWithVoc12(CFLOBDD_FLOAT_BOOST c); // convert the vector into another vector with extra volcabularies
		extern CFLOBDD_FLOAT_BOOST KroneckerProduct(CFLOBDD_FLOAT_BOOST m1, CFLOBDD_FLOAT_BOOST m2);
		extern CFLOBDD_FLOAT_BOOST VectorShiftVocs1To2(CFLOBDD_FLOAT_BOOST m1);
		extern CFLOBDD_FLOAT_BOOST VectorPadWithZeros(CFLOBDD_FLOAT_BOOST c, unsigned int level);
		extern CFLOBDD_FLOAT_BOOST NoDistinctionNode(unsigned int level, int val);
		extern CFLOBDD_FLOAT_BOOST ConvertToDouble(CFLOBDD_FLOAT_BOOST n);
//#ifdef PATH_COUNTING_ENABLED
		extern std::string Sampling(CFLOBDD_FLOAT_BOOST c, bool voctwo, std::string = "");
		extern std::string SamplingV2(CFLOBDD_FLOAT_BOOST c);
//#endif
		extern CFLOBDD_FLOAT_BOOST VectorWithAmplitude(CFLOBDD_FLOAT_BOOST c);
		extern void VectorPrintColumnMajor(CFLOBDD_FLOAT_BOOST c, std::ostream & out);
		extern void VectorPrintColumnMajorInterleaved(CFLOBDD_FLOAT_BOOST c, std::ostream & out);
	}
}

#endif

