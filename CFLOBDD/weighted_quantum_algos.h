#ifndef WEIGHTED_QUANTUMALGOS_INT_GUARD
#define WEIGHTED_QUANTUMALGOS_INT_GUARD


#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include "weighted_cflobdd_node_t.h"
#include "wmatrix1234_fb_mul.h"
#include "wmatrix1234_complex_fb_mul.h"
#include "wmatrix1234_fourier_mul.h"

namespace CFL_OBDD {

	namespace WeightedQuantumAlgos {

		// Initialization routine
		extern void QuantumAlgosInitializer();

		extern std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> DeutschJozsaAlgo(unsigned int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F);
		extern std::vector<std::string> SimonsAlgo(int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F);
		extern std::vector<std::string> SimonsAlgoV2(int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F);
		extern std::vector<std::string> SimonsAlgoV3(int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F);
		extern std::pair<WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL, std::vector<std::string>> SimonsAlgoV4(int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F);
		extern std::pair<WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL, std::vector<std::string>> SimonsAlgoV4New(int n, std::string s);
		extern std::vector<std::string> SimonsAlgoV4_Voc2(int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F);
		extern std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> GroversAlgoWithV4(int n, std::string s);
		extern std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> BV(long long int n, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL F);
		extern std::pair<std::string, WEIGHTED_CFLOBDD_FLOAT_BOOST_MUL> GHZ(unsigned long long int n);
		extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL QFT(long long int n, std::string m);
        extern WEIGHTED_CFLOBDD_FOURIER_MUL QFT_fourier(long long int n, std::string s);
        extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ShorsAlgoNew(int a, int N);
        extern WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ShorsAlgo(int a, int N, int bits);
        extern std::tuple<WEIGHTED_CFLOBDD_FOURIER_MUL, std::string, BIG_INT> ShorsFourier(int a, int N);
	}
}

#endif

