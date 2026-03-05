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
#include "vector_double_top_node.h"
#include "vector_float_boost_top_node.h"

//using namespace boost::random;
//mt19937 gen;

namespace CFL_OBDD {

	namespace VectorFloatBoost {

		void VectorInitializerTop()
		{
			VectorInitializerNode();
			return;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkBasisVectorTop(unsigned int level, unsigned int index)
		{
			CFLOBDDTopNodeFloatBoostRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle rhandle;

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

			ptr = new CFLOBDDTopNodeFloatBoost(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkBasisVectorTop(unsigned int level, std::string s)
		{
			CFLOBDDTopNodeFloatBoostRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle rhandle;

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

			ptr = new CFLOBDDTopNodeFloatBoost(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeFloatBoostRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorToMatrixInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));

			FloatBoostReturnMapHandle rhandle;
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
				rhandle.AddToEnd(n->rootConnection.returnMapHandle[i]);
			rhandle.Canonicalize();
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, rhandle);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MatrixToVectorTop(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixToVectorNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr NoDistinctionNodeTop(unsigned int level, int val)
		{
			FloatBoostReturnMapHandle m1;
			m1.AddToEnd(val);
			m1.Canonicalize();

			return new CFLOBDDTopNodeFloatBoost(CFLOBDDNodeHandle::NoDistinctionNode[level], m1);
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkColumn1MatrixTop(unsigned int level)
		{
			CFLOBDDTopNodeFloatBoostRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			FloatBoostReturnMapHandle rhandle;

			tempHandle = MkColumn1MatrixNode(level);

			rhandle.AddToEnd(1);
			rhandle.AddToEnd(0);
			rhandle.Canonicalize();

			ptr = new CFLOBDDTopNodeFloatBoost(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeFloatBoostRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;
			CFLOBDDNodeHandle tempHandle = MkVectorWithVoc12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 1);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorShiftVocs1To2Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeFloatBoostRefPtr v = new CFLOBDDTopNodeFloatBoost(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeFloatBoostRefPtr ConvertToDoubleTop(CFLOBDDTopNodeFloatBoostRefPtr c)
		{
			FloatBoostReturnMapHandle tmp;
			for (int i = 0; i < c->rootConnection.returnMapHandle.Size(); i++)
			{
				if (!tmp.Member(c->rootConnection.returnMapHandle[i]))
					tmp.AddToEnd(c->rootConnection.returnMapHandle[i].convert_to<double>());
			}
			tmp.Canonicalize();
			ReductionMapHandle inducedReductionMapHandle;
			ReturnMapHandle<BIG_FLOAT> inducedReturnMap;
			tmp.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//     CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_n = c->rootConnection.entryPointHandle->Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			return (new CFLOBDDTopNodeFloatBoost(reduced_n, inducedReturnMap));
		}

		CFLOBDDTopNodeFloatBoostRefPtr VectorWithAmplitudeTop(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			FloatBoostReturnMapHandle rhandle;
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				BIG_FLOAT val = n->rootConnection.returnMapHandle[i];
				rhandle.AddToEnd(val * val);
			}

			ReductionMapHandle inducedReductionMapHandle;
			ReturnMapHandle<BIG_FLOAT> inducedReturnMap;
			rhandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//     CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_n = n->rootConnection.entryPointHandle->Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			return (new CFLOBDDTopNodeFloatBoost(reduced_n, inducedReturnMap));
		}

//#ifdef PATH_COUNTING_ENABLED
		std::string SamplingTop(CFLOBDDTopNodeFloatBoostRefPtr n, bool VocTwo, std::string func)
		{
			std::vector<std::pair<BIG_FLOAT, unsigned int>> values;
			long double prob = -1 * std::numeric_limits<long double>::infinity();
			
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				if (n->rootConnection.returnMapHandle.Lookup(i) == 0){
					values.push_back(std::make_pair(-1 * std::numeric_limits<BIG_FLOAT>::infinity(), i));
				}
				else{
					BIG_FLOAT amplitude = boost::multiprecision::log2(n->rootConnection.returnMapHandle.Lookup(i));
					long double logNumPaths = n->rootConnection.entryPointHandle->handleContents->numPathsToExit[i];
					values.push_back(std::make_pair(amplitude + logNumPaths, i));
				}
			}


			sort(values.begin(), values.end(), sortNumPathPairs<BIG_FLOAT>);
			prob = getLogSumNumPaths(values, values.size()).convert_to<long double>();
			
			// for (int j = 0; j < values.size(); j++)
			// 	std::cout << values[j].first << " " << values[j].second << std::endl;
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
				// std::cout << rand_val << " " << (((double)rand_val) / RAND_MAX) << " " << pow(2, prob) << " " << random_value;
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
			if (n->rootConnection.returnMapHandle.Size() == 2){
				if (n->rootConnection.returnMapHandle.Lookup(0) == 0)
					index = 1;
				else if (n->rootConnection.returnMapHandle.Lookup(1) == 0)
					index = 0;
			}
			if (n->rootConnection.returnMapHandle.Size() == 2 && func == "Grovers"){
				BIG_FLOAT a = n->rootConnection.returnMapHandle[0];
				BIG_FLOAT b = n->rootConnection.returnMapHandle[1];
				if (abs(a - (a - b)) < 0.01)
					index = 0;
				else if (abs(b - (b - a)) < 0.01)
					index = 1;
			}
			std::pair<std::string, std::string> stringPair = SamplingNode(*(n->rootConnection.entryPointHandle), index, VocTwo);
			//std::cout << stringPair.first << " " << stringPair.second << std::endl;
			return stringPair.first + stringPair.second;
		}

		// V2 has binary values
		std::string SamplingV2Top(CFLOBDDTopNodeFloatBoostRefPtr n)
		{
			unsigned int index = 0;
			if (n->rootConnection.returnMapHandle.Size() == 2)
				index = 1;
			std::pair<std::string, std::string> stringPair = SamplingNode(*(n->rootConnection.entryPointHandle), index);
			return stringPair.first + stringPair.second;
		}
//#endif

		void VectorPrintColumnMajorTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out)
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

		void VectorPrintColumnMajorInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out)
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

