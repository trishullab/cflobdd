#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <chrono>
#include <random>
#include <unordered_map>
#include "wvector_top_node_fourier_mul.h"
#include "wvector_fourier_mul_node.h"


namespace CFL_OBDD {

	namespace WeightedVectorFourierMul {

		void VectorInitializerTop()
		{
			VectorInitializerNode();
			return;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkBasisVectorTop(unsigned int level, unsigned int index)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr ptr;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle rhandle;

			tempHandle = MkBasisVectorNode(level, index);
			if (index == 0)
			{
				rhandle.AddToEnd(fourierSemiring(1, 1));
				rhandle.AddToEnd(fourierSemiring(0, 1));
			}
			else
			{
				rhandle.AddToEnd(fourierSemiring(0, 1));
				rhandle.AddToEnd(fourierSemiring(1, 1));
			}
			rhandle.Canonicalize();

			ptr = new WeightedCFLOBDDTopNodeFourier(tempHandle, rhandle, fourierSemiring(1, 1));
			return ptr;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkBasisVectorTop(unsigned int level, std::string s)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr ptr;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle rhandle;

			tempHandle = MkBasisVectorNode(level, s);
			if (s.find('1') == std::string::npos)
			{
				rhandle.AddToEnd(fourierSemiring(1, 1));
				rhandle.AddToEnd(fourierSemiring(0, 1));
			}
			else
			{
				rhandle.AddToEnd(fourierSemiring(0, 1));
				rhandle.AddToEnd(fourierSemiring(1, 1));
			}
			rhandle.Canonicalize();

			ptr = new WeightedCFLOBDDTopNodeFourier(tempHandle, rhandle, fourierSemiring(1, 1));
			return ptr;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr VectorToMatrixInterleavedTop(WeightedCFLOBDDTopNodeFourierRefPtr n)
		{
			// std::unordered_map<WeightedCFLOBDDFourierMulNodeHandle, WeightedCFLOBDDFourierMulNodeHandle, WeightedCFLOBDDFourierMulNodeHandle::WeightedCFLOBDDNodeHandleT_Hash> hashMap;
			
			// WeightedCFLOBDDFloatBoostMulNodeMemoTableRefPtr memoTable = new WeightedCFLOBDDFloatBoostMulNodeMemoTable;
			std::unordered_map<WeightedCFLOBDDFourierMulNodeHandle, WeightedCFLOBDDFourierMulNodeHandle, WeightedCFLOBDDFourierMulNodeHandle::WeightedCFLOBDDNodeHandleT_Hash> hashMap;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle = VectorToMatrixInterleavedNode(hashMap, *(n->rootConnection.entryPointHandle));
			WeightedCFLOBDDTopNodeFourierRefPtr v = new WeightedCFLOBDDTopNodeFourier(tempHandle, n->rootConnection.returnMapHandle, n->rootConnection.factor);
			return v;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr NoDistinctionNodeTop(unsigned int level, fourierSemiring val)
		{
			if (val == fourierSemiring(0, 1))
			{
				FourierReturnMapHandle m0;
				m0.AddToEnd(fourierSemiring(0, 1));
				m0.Canonicalize();
				return new WeightedCFLOBDDTopNodeFourier(WeightedCFLOBDDFourierMulNodeHandle::NoDistinctionNode_Ann[level], m0, fourierSemiring(0, 1));
			}
			FourierReturnMapHandle m1;
			m1.AddToEnd(fourierSemiring(1, 1));
			m1.Canonicalize();

			return new WeightedCFLOBDDTopNodeFourier(WeightedCFLOBDDFourierMulNodeHandle::NoDistinctionNode[level], m1, val);
		}

		WeightedCFLOBDDTopNodeFourierRefPtr MkColumn1MatrixTop(unsigned int level)
		{
			WeightedCFLOBDDTopNodeFourierRefPtr ptr;
			WeightedCFLOBDDFourierMulNodeHandle tempHandle;
			FourierReturnMapHandle rhandle;

			tempHandle = MkColumn1MatrixNode(level);

			rhandle.AddToEnd(fourierSemiring(1, 1));
			rhandle.AddToEnd(fourierSemiring(0, 1));
			rhandle.Canonicalize();

			ptr = new WeightedCFLOBDDTopNodeFourier(tempHandle, rhandle, fourierSemiring(1, 1));
			return ptr;
		}

		WeightedCFLOBDDTopNodeFourierRefPtr VectorWithAmplitudeTop(WeightedCFLOBDDTopNodeFourierRefPtr n)
		{
			FourierReturnMapHandle rhandle;
			for (unsigned int i = 0; i < n->rootConnection.returnMapHandle.Size(); i++)
			{
				fourierSemiring val = n->rootConnection.returnMapHandle[i];
				rhandle.AddToEnd(val * val);
			}

			ReductionMapHandle inducedReductionMapHandle;
			ReturnMapHandle<fourierSemiring> inducedReturnMap;
			rhandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
			return n;
		}

//#ifdef PATH_COUNTING_ENABLED
		std::string SamplingTop(WeightedCFLOBDDTopNodeFourierRefPtr n, bool VocTwo, std::string func)
		{
			if (n->rootConnection.factor == fourierSemiring(0, 1))
				return "";
			int index = n->rootConnection.returnMapHandle.LookupInv(fourierSemiring(1, 1));
			std::pair<std::string, std::string> stringPair = SamplingNode(*(n->rootConnection.entryPointHandle), index, VocTwo, func);
			//std::cout << stringPair.first << " " << stringPair.second << std::endl;
			return stringPair.first + stringPair.second;
		}

		void VectorPrintColumnMajorTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out)
		{
			// unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			// if (level >= 2 && level <= 4 || true) {
			// 	unsigned int indexBits = 1 << (level - 1);
			// 	unsigned int totalBits = 2 * indexBits;
			// 	unsigned long int rows = 1UL << indexBits;
			// 	unsigned long int cols = rows;
			// 	SH_OBDD::Assignment a(totalBits);
			// 	for (unsigned long int i = 0UL; i < rows; i++) { // Vocs 1
			// 		// Fill in the even positions of a 
			// 		unsigned long int maskVoc1 = 1UL;
			// 		for (int k = (indexBits - 1); k >= 0; k--) {
			// 			a[2 * k] = (i & maskVoc1);
			// 			maskVoc1 = maskVoc1 << 1;
			// 		}
			// 		for (unsigned long int j = 0UL; j < cols; j++) {  // Vocs 2
			// 			// Fill in the odd positions of a
			// 			unsigned long int maskVoc2 = 1UL;
			// 			for (int k = (indexBits - 1); k >= 0; k--) {
			// 				a[2 * k + 1] = (j & maskVoc2);
			// 				maskVoc2 = maskVoc2 << 1;
			// 			}
			// 			//out << "i = " << i << ", j = " << j << ", a = ";
			// 			//a.print(out);
			// 			//out << std::endl;

			// 			double b1 = n->EvaluateIteratively(a).convert_to<double>();
			// 			if (b1 != 0)
			// 				out << i << " " << j << std::endl;
			// 			/*std::string b;
			// 			b = n->EvaluateIteratively(a).str(5);
			// 			out << b << " ";*/
			// 		}
			// 	}
			// 	out << std::endl;
			// }
			// else {
			// 	std::cerr << "Cannot print matrix: level must be in [2 .. 4]" << std::endl;
			// }
		}

		void VectorPrintColumnMajorInterleavedTop(WeightedCFLOBDDTopNodeFourierRefPtr n, std::ostream & out)
		{
			unsigned int level = n->rootConnection.entryPointHandle->handleContents->level;
			// if (level >= 1 && level <= 4) {
			// 	unsigned int indexBits = 1 << (level - 1);
			// 	unsigned int totalBits = 2 * indexBits;
			// 	unsigned long int rows = 1UL << totalBits;
			// 	SH_OBDD::Assignment a(totalBits);
			// 	for (unsigned long int i = 0UL; i < rows; i++) {
			// 		// Fill in the even positions of a 
			// 		unsigned long int mask1 = 1UL;
			// 		for (int k = totalBits - 1; k >= 0; k--) {
			// 			a[k] = (i & mask1);
			// 			mask1 = mask1 << 1;
			// 		}
			// 		std::string b;
			// 		b = n->EvaluateIteratively(a).str(5);
			// 		out << b << " ";
			// 	}
			// 	out << std::endl;
			// }
			// else {
			// 	std::cerr << "Cannot print matrix: level must be in [1 .. 4]" << std::endl;
			// }
		}
	}
}

