#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <chrono>
#include <random>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "vector_node.h"
#include "vector_complex_float_boost_top_node.h"

//using namespace boost::random;
//mt19937 gen;

namespace CFL_OBDD {

	namespace VectorComplexFloatBoost {

		void VectorInitializerTop()
		{
			VectorInitializerNode();
			return;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkBasisVectorTop(unsigned int level, unsigned int index)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle rhandle;

			tempHandle = MkBasisVectorNode(level, index);
			if (index == 0)
			{
				rhandle.AddToEnd(1);
				rhandle.AddToEnd(0);
			}
			else
			{
				rhandle.AddToEnd(0);
				rhandle.AddToEnd(1);
			}
			rhandle.Canonicalize();

			ptr = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkBasisVectorTop(unsigned int level, std::string s)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle rhandle;

			tempHandle = MkBasisVectorNode(level, s);
			if (s.find('1') == std::string::npos)
			{
				rhandle.AddToEnd(1);
				rhandle.AddToEnd(0);
			}
			else
			{
				rhandle.AddToEnd(0);
				rhandle.AddToEnd(1);
			}
			rhandle.Canonicalize();

			ptr = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorToMatrixInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));

			ComplexFloatBoostReturnMapHandle rhandle;
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
				rhandle.AddToEnd(n->rootConnection.returnMapHandle[i]);
			rhandle.Canonicalize();
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, rhandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MatrixToVectorTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixToVectorNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr NoDistinctionNodeTop(unsigned int level, BIG_COMPLEX_FLOAT val)
		{
			ComplexFloatBoostReturnMapHandle m1;
			m1.AddToEnd(val);
			m1.Canonicalize();

			return new CFLOBDDTopNodeComplexFloatBoost(CFLOBDDNodeHandle::NoDistinctionNode[level], m1);
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkColumn1MatrixTop(unsigned int level)
		{
			CFLOBDDTopNodeComplexFloatBoostRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			ComplexFloatBoostReturnMapHandle rhandle;

			tempHandle = MkColumn1MatrixNode(level);

			rhandle.AddToEnd(1);
			rhandle.AddToEnd(0);
			rhandle.Canonicalize();

			ptr = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;
			CFLOBDDNodeHandle tempHandle = MkVectorWithVoc12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 1);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorShiftVocs1To2Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeComplexFloatBoostRefPtr v = new CFLOBDDTopNodeComplexFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeComplexFloatBoostRefPtr VectorWithAmplitudeTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			ComplexFloatBoostReturnMapHandle rhandle;
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				BIG_COMPLEX_FLOAT val = n->rootConnection.returnMapHandle[i];
				auto real_val = val.real(); auto imag_val = val.imag();
				auto abs_val = real_val * real_val + imag_val * imag_val;
				rhandle.AddToEnd(abs_val);
			}

			ReductionMapHandle inducedReductionMapHandle;
			ReturnMapHandle<BIG_COMPLEX_FLOAT> inducedReturnMap;
			rhandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//     CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_n = n->rootConnection.entryPointHandle->Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			return (new CFLOBDDTopNodeComplexFloatBoost(reduced_n, inducedReturnMap));
		}

		long double getNonZeroProbabilityTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			long double prob = 0;

			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				if (n->rootConnection.returnMapHandle.Lookup(i) != 0){
					long double v = n->rootConnection.returnMapHandle.Lookup(i).real().convert_to<long double>();
					long double logNumPaths = n->rootConnection.entryPointHandle->handleContents->numPathsToExit[i];
					prob += v * std::pow(2, logNumPaths);
				}
			}
			return prob;
		}

		//#ifdef PATH_COUNTING_ENABLED
		std::string SamplingTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, bool VocTwo)
		{
			std::vector<std::pair<BIG_FLOAT, unsigned int>> values;
			long double prob = -1 * std::numeric_limits<long double>::infinity();

			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				if (n->rootConnection.returnMapHandle.Lookup(i) == 0){
					values.push_back(std::make_pair(-1 * std::numeric_limits<BIG_FLOAT>::infinity(), i));
				}
				else{
					BIG_FLOAT v = n->rootConnection.returnMapHandle.Lookup(i).real().convert_to<BIG_FLOAT>();
					BIG_FLOAT amplitude = mp::log2(v);
					long double logNumPaths = n->rootConnection.entryPointHandle->handleContents->numPathsToExit[i];
					values.push_back(std::make_pair(amplitude + logNumPaths, i));
				}
			}


			sort(values.begin(), values.end(), sortNumPathPairs<BIG_FLOAT>);
			prob = getLogSumNumPaths(values, values.size()).convert_to<long double>();

			BIG_FLOAT val = -1 * std::numeric_limits<BIG_FLOAT>::infinity();
			long double random_value = 0.0;
			if (prob >= 64){
				std::random_device rd;
				std::default_random_engine generator(rd());
				std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF);
				random_value = log2l(distribution(generator)) + prob - 64;
			}
			else{
				auto rand_val = rand();
				random_value = log2l((((double)rand_val) / RAND_MAX)*pow(2, prob));
			}

			unsigned int index = 0;
			for (unsigned int i = 0; i < values.size(); i++)
			{
				val = getLogSumNumPaths(values, i + 1);
				if (val >= random_value)
				{
					index = values[i].second;
					break;
				}
			}
			std::pair<std::string, std::string> stringPair = SamplingNode(*(n->rootConnection.entryPointHandle), index, VocTwo);
			return stringPair.first + stringPair.second;
		}

		// V2 has binary values
		std::string SamplingV2Top(CFLOBDDTopNodeComplexFloatBoostRefPtr n)
		{
			unsigned int index = 0;
			if (n->rootConnection.returnMapHandle.Size() == 2)
				index = 1;
			std::pair<std::string, std::string> stringPair = SamplingNode(*(n->rootConnection.entryPointHandle), index);
			return stringPair.first + stringPair.second;
		}
		//#endif

		void VectorPrintColumnMajorTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 2 && level <= 4 || true) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << indexBits;
				unsigned long int cols = rows;
				SH_OBDD::Assignment a(totalBits);
				for (unsigned long int i = 0UL; i < rows; i++) { // Vocs 1
					// Fill in the even positions of a 
					unsigned long int maskVoc1 = 1UL;
					for (int k = (indexBits - 1); k >= 0; k--) {
						a[2 * k] = (i & maskVoc1);
						maskVoc1 = maskVoc1 << 1;
					}
					for (unsigned long int j = 0UL; j < cols; j++) {  // Vocs 2
						// Fill in the odd positions of a
						unsigned long int maskVoc2 = 1UL;
						for (int k = (indexBits - 1); k >= 0; k--) {
							a[2 * k + 1] = (j & maskVoc2);
							maskVoc2 = maskVoc2 << 1;
						}
						//out << "i = " << i << ", j = " << j << ", a = ";
						//a.print(out);
						//out << std::endl;

						double b1 = n->EvaluateIteratively(a).convert_to<double>();
						if (b1 != 0)
							out << i << " " << j << std::endl;
						/*std::string b;
						b = n->EvaluateIteratively(a).str(5);
						out << b << " ";*/
					}
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [2 .. 4]" << std::endl;
			}
		}

		void VectorPrintColumnMajorInterleavedTop(CFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 1 && level <= 4) {
				unsigned int indexBits = 1 << (level - 1);
				unsigned int totalBits = 2 * indexBits;
				unsigned long int rows = 1UL << totalBits;
				SH_OBDD::Assignment a(totalBits);
				for (unsigned long int i = 0UL; i < rows; i++) {
					// Fill in the even positions of a 
					unsigned long int mask1 = 1UL;
					for (int k = totalBits - 1; k >= 0; k--) {
						a[k] = (i & mask1);
						mask1 = mask1 << 1;
					}
					std::string b;
					b = n->EvaluateIteratively(a).str(5);
					out << b << " ";
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [1 .. 4]" << std::endl;
			}
		}
	}
}

