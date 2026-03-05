#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <chrono>
#include <random>
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/random.hpp>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "vector_node.h"
#include "vector_double_top_node.h"

//using namespace boost::random;
//mt19937 gen;

namespace CFL_OBDD {

	namespace VectorDouble {

		void VectorInitializerTop()
		{
			VectorInitializerNode();
			return;
		}

		CFLOBDDTopNodeDoubleRefPtr MkBasisVectorTop(unsigned int level, unsigned int index)
		{
			CFLOBDDTopNodeDoubleRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle rhandle;

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

			ptr = new CFLOBDDTopNodeDouble(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeDoubleRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeDoubleRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorToMatrixInterleavedNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr MatrixToVectorTop(CFLOBDDTopNodeDoubleRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = MatrixToVectorNode(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr NoDistinctionNodeTop(unsigned int level)
		{
			DoubleReturnMapHandle m1;
			m1.AddToEnd(1);
			m1.Canonicalize();

			return new CFLOBDDTopNodeDouble(CFLOBDDNodeHandle::NoDistinctionNode[level], m1);
		}

		CFLOBDDTopNodeDoubleRefPtr MkColumn1MatrixTop(unsigned int level)
		{
			CFLOBDDTopNodeDoubleRefPtr ptr;
			CFLOBDDNodeHandle tempHandle;
			DoubleReturnMapHandle rhandle;

			tempHandle = MkColumn1MatrixNode(level);

			rhandle.AddToEnd(1);
			rhandle.AddToEnd(0);
			rhandle.Canonicalize();

			ptr = new CFLOBDDTopNodeDouble(tempHandle, rhandle);
			return ptr;
		}

		CFLOBDDTopNodeDoubleRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeDoubleRefPtr n)
		{
			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;
			CFLOBDDNodeHandle tempHandle = MkVectorWithVoc12Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeDoubleRefPtr n)
		{
			assert(n->rootConnection.entryPointHandle->handleContents->level >= 1);

			CFLOBDDNodeMemoTableRefPtr memoTable = new CFLOBDDNodeMemoTable;

			CFLOBDDNodeHandle tempHandle = VectorShiftVocs1To2Node(memoTable, *(n->rootConnection.entryPointHandle));
			CFLOBDDTopNodeDoubleRefPtr v = new CFLOBDDTopNodeDouble(tempHandle, n->rootConnection.returnMapHandle);
			return v;
		}

		CFLOBDDTopNodeDoubleRefPtr VectorWithAmplitudeTop(CFLOBDDTopNodeDoubleRefPtr n)
		{
			DoubleReturnMapHandle rhandle;
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				double val = n->rootConnection.returnMapHandle[i];
				//rhandle.AddToEnd(val * val);
				rhandle.AddToEnd(fabs(val));
			}

			ReductionMapHandle inducedReductionMapHandle;
			ReturnMapHandle<double> inducedReturnMap;
			rhandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			//     CFLOBDDNodeHandle::InitReduceCache();
			CFLOBDDNodeHandle reduced_n = n->rootConnection.entryPointHandle->Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
			return (new CFLOBDDTopNodeDouble(reduced_n, inducedReturnMap));
		}

#ifdef PATH_COUNTING_ENABLED
		std::string SamplingTop(CFLOBDDTopNodeDoubleRefPtr n)
		{
			//std::vector<cpp_bin_float_50> values;
			//std::vector<int> values;
			std::vector<std::pair<long double, unsigned int>> values;
			//assert(n->rootConnection.returnMapHandle.Size() == n->rootConnection.entryPointHandle.handleContents->numExits);
			//cpp_bin_float_50 prob = 0.0;
			//int prob = 0.0;
			long double prob = -1 * std::numeric_limits<long double>::infinity();
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				//cpp_bin_float_50 amplitude = n->rootConnection.returnMapHandle.Lookup(i);
				//cpp_int numPaths = n->rootConnection.entryPointHandle.handleContents->numPathsToExit[i];
				//values.push_back((amplitude * numPaths.convert_to<cpp_bin_float_50>()));
				if (n->rootConnection.returnMapHandle.Lookup(i) == 0){
					values.push_back(std::make_pair(-1 * std::numeric_limits<long double>::infinity(), i));
				}
				else{
					long double amplitude = log2l(n->rootConnection.returnMapHandle.Lookup(i));
					long double logNumPaths = n->rootConnection.entryPointHandle->handleContents->numPathsToExit[i];
					//std::cout << amplitude << " " << logNumPaths << " " << (amplitude + logNumPaths) << std::endl;
					values.push_back(std::make_pair(amplitude + logNumPaths, i));
				}
				//prob += (amplitude * numPaths.convert_to<cpp_bin_float_50>());
			}
			sort(values.begin(), values.end(), sortNumPathPairs<long double>);
			prob = getLogSumNumPaths(values, values.size());
			//std::cout << prob << std::endl;
			//cpp_bin_float_50 val = 0.0;
			//cpp_bin_float_50 random_val = generate_canonical<cpp_bin_float_50, std::numeric_limits<cpp_bin_float_50>::digits>(gen);
			long double val = -1 * std::numeric_limits<long double>::infinity();
			long double random_value = 0.0;
			if (prob >= 64){
				std::random_device rd;
				std::default_random_engine generator(rd());
				std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF);
				random_value = log2l(distribution(generator)) + prob - 64;
			}
			else{
				random_value = log2l((((double)rand()) / RAND_MAX)*pow(2, prob));
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
			/*if (n->rootConnection.returnMapHandle.Lookup(0) == 0)
				index = 1;
			else
				index = 0;*/
			std::pair<std::string,std::string> stringPair = SamplingNode(*(n->rootConnection.entryPointHandle), index);
			return stringPair.first + stringPair.second;
		}
#endif

		void VectorPrintColumnMajorTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			if (level >= 2 && level <= 4) {
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
						double b;
						b = n->EvaluateIteratively(a);
						if (b != 0)
							out << b << " " << i << std::endl;
					}
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [2 .. 4]" << std::endl;
			}
		}

		void VectorPrintColumnMajorInterleavedTop(CFLOBDDTopNodeDoubleRefPtr n, std::ostream & out)
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
					double b;
					b = n->EvaluateIteratively(a);
					/*if (b != 0)
					{
						out << i << std::endl;
					}*/
					//out << i << " " << b << std::endl;
					if (b != 0)
						out << b << " " << i << std::endl;
					//out << std::endl;
				}
				out << std::endl;
			}
			else {
				std::cerr << "Cannot print matrix: level must be in [1 .. 4]" << std::endl;
			}
		}
	}
}

