#ifndef QUANTUMALGOS_INT_GUARD
#define QUANTUMALGOS_INT_GUARD


#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include "cflobdd_int.h"
#include "matrix1234_float_boost.h"
#include "matrix1234_complex_float_boost.h"
#include "matrix1234_double.h"

namespace CFL_OBDD {

	namespace QuantumAlgos {

		// Initialization routine
		extern void QuantumAlgosInitializer();

		extern CFLOBDD_DOUBLE DeutschJozsaAlgo(unsigned int n, CFLOBDD_DOUBLE F);
		extern std::pair<std::string, CFLOBDD_FLOAT_BOOST> DeutschJozsaAlgo(unsigned int n, CFLOBDD_FLOAT_BOOST F);
		extern std::vector<std::string> SimonsAlgo(int n, CFLOBDD_DOUBLE F);
		extern std::vector<std::string> SimonsAlgo(int n, CFLOBDD_FLOAT_BOOST F);
		extern std::vector<std::string> SimonsAlgoV2(int n, CFLOBDD_FLOAT_BOOST F);
		extern std::vector<std::string> SimonsAlgoV3(int n, CFLOBDD_FLOAT_BOOST F);
		extern std::pair<CFLOBDD_FLOAT_BOOST, std::vector<std::string>> SimonsAlgoV4(int n, CFLOBDD_FLOAT_BOOST F);
		extern std::pair<CFLOBDD_FLOAT_BOOST, std::vector<std::string>> SimonsAlgoV4New(int n, std::string s);
		extern std::vector<std::string> SimonsAlgoV4_Voc2(int n, CFLOBDD_FLOAT_BOOST F);
		extern std::string GroversAlgo(int n, std::string s);
		extern std::pair<std::string, CFLOBDD_FLOAT_BOOST> GroversAlgoWithV4(int n, std::string s);
		extern std::string GroversAlgoWithV4_double(int n, std::string s);
		extern std::pair<std::string, CFLOBDD_FLOAT_BOOST> BV(long long int n, CFLOBDD_FLOAT_BOOST F);
		extern std::pair<std::string, CFLOBDD_FLOAT_BOOST> GHZ(unsigned long long int n);
		extern CFLOBDD_COMPLEX_BIG QFT(long long int n, std::string m);
		extern CFLOBDD_COMPLEX_BIG ShorsAlgoNew(int a, int n);
		std::pair<CFLOBDD_COMPLEX_BIG, std::vector<std::string>> ShorsAlgo(int n, CFLOBDD_COMPLEX_BIG F);
	}
}

#endif

