#include <cassert>
#include <cstdio>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
//#include <boost/multiprecision/cpp_int.hpp>
//#include <boost/random.hpp>
#include "vector_node.h"
// #include "linear_map_subst.h"
#include "cflobdd_node.h"
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"
// #include "cflobdd_top_node_linear_map.h"

//using namespace boost::multiprecision;
//using namespace boost::random;

//typedef independent_bits_engine<mt19937, 256, cpp_int> generator_type;
//generator_type gen;

namespace CFL_OBDD {

	std::vector<ReturnMapHandle<int>> commonly_used_return_maps_vector;// m0, m1, m01, m10

	void InitReturnMapVectorHandles(){
		ReturnMapHandle<int> m0, m1, m01, m10;
		m0.AddToEnd(0);
		m0.Canonicalize();
		m1.AddToEnd(1);
		m1.Canonicalize();
		m01.AddToEnd(0);
		m01.AddToEnd(1);
		m01.Canonicalize();
		m10.AddToEnd(1);
		m10.AddToEnd(0);
		m10.Canonicalize();
		commonly_used_return_maps_vector.push_back(m0);
		commonly_used_return_maps_vector.push_back(m1);
		commonly_used_return_maps_vector.push_back(m01);
		commonly_used_return_maps_vector.push_back(m10);
	}

	void VectorInitializerNode()
	{
		// Empty for now
		InitReturnMapVectorHandles();
		return;
	}

	CFLOBDDNodeHandle MkBasisVectorNode(unsigned int level, unsigned int index)
	{
		if (level == 0)
		{
			assert(index < 2);
			return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
		}

		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			assert(index < 4);

			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[2]);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			

			if ((index & 2 ) == 0)
			{
				if ((index & 1) == 0){
					n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[2]);//m01
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps_vector[1]);//m1
				}
				else{
					n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[2]);//m01
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps_vector[0]);//m0
				}
			}
			else
			{
				if ((index & 1) == 0)
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[3]);//m10
				else
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[2]);//m01
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps_vector[0]);//m0
			}
		}
		else
		{

			// unsigned int higherOrderIndex = index & ((1 << (1 << level)) - (1 << (1 << (level - 1))));
			unsigned int higherOrderIndex = index >> (1 << (level - 1));
			CFLOBDDNodeHandle tempANodeHandle = MkBasisVectorNode(level-1, higherOrderIndex);
			n->AConnection = Connection(tempANodeHandle, commonly_used_return_maps_vector[2]);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			unsigned int lowerOrderIndex = index & ((1 << (1 << (level - 1))) - 1);
			CFLOBDDNodeHandle tempBNodeHandle = MkBasisVectorNode(level - 1, lowerOrderIndex);
			if (higherOrderIndex == 0)
			{
				if (lowerOrderIndex == 0)
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps_vector[1]);//m1
				else
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps_vector[0]);//m0
				n->BConnection[0] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
			}
			else{
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps_vector[0]);//m0
				if (lowerOrderIndex == 0)
					n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[3]);//m10
				else
					n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
			}
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkBasisVectorNode(unsigned int level, std::string s)
	{
		if (level == 0)
		{
			assert(s.length() == 1);
			return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
		}

		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			assert(s.length() == 2);
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];


			if (s == "00" || s == "01")
			{
				if (s == "00"){
					n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[2]);//m01
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps_vector[1]);//m1
				}
				else{
					n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[2]);//m01
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps_vector[0]);//m0
				}
			}
			else
			{
				if (s == "10")
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[3]);//m10
				else
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps_vector[2]);//m01
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps_vector[0]);//m0
			}
		}
		else
		{

			// unsigned int higherOrderIndex = index & ((1 << (1 << level)) - (1 << (1 << (level - 1))));
			std::string first_half_s = s.substr(0, s.length() / 2);
			CFLOBDDNodeHandle tempANodeHandle = MkBasisVectorNode(level - 1, first_half_s);
			n->AConnection = Connection(tempANodeHandle, commonly_used_return_maps_vector[2]);//m01

			
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			std::string half_string = s.substr(s.length() / 2);
			CFLOBDDNodeHandle tempBNodeHandle = MkBasisVectorNode(level - 1, half_string);
			if (first_half_s.find('1') == std::string::npos)
			{
				if (half_string.find('1') == std::string::npos)
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps_vector[1]);//m1
				else
					n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps_vector[0]);//m0
				n->BConnection[0] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
			}
			else{
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps_vector[0]);//m0
				if (half_string.find('1') == std::string::npos)
					n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[3]);//m10
				else
					n->BConnection[1] = Connection(tempBNodeHandle, commonly_used_return_maps_vector[2]);//m01
			}
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle VectorToMatrixInterleavedNode(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;

		if (nh == CFLOBDDNodeHandle::NoDistinctionNode[nh.handleContents->level])
			return CFLOBDDNodeHandle::NoDistinctionNode[nh.handleContents->level + 1];

		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nh.handleContents->level + 1);

		if (nh.handleContents->level == 0)
		{
			if (nhNode->numExits == 1)
				return CFLOBDDNodeHandle::NoDistinctionNode[1];
			assert(nhNode->numExits == 2);

			n->AConnection = Connection(nh, commonly_used_return_maps_vector[2]);//m01
			n->numBConnections = nhNode->numExits;
			n->BConnection = new Connection[n->numBConnections];
			for (int i = 0; i < nhNode->numExits; i++)
			{
				CFLOBDDReturnMapHandle m2;
				m2.AddToEnd(i);
				m2.Canonicalize();
				n->BConnection[i] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[0], m2);
			}
		}
		else
		{
			CFLOBDDNodeHandle tempHandle = VectorToMatrixInterleavedNode(memoTable, *(nhNode->AConnection.entryPointHandle));
			CFLOBDDReturnMapHandle mA;
			for (unsigned int i = 0; i < nhNode->AConnection.returnMapHandle.Size(); i++)
				mA.AddToEnd(nhNode->AConnection.returnMapHandle[i]);
			mA.Canonicalize();
			n->AConnection = Connection(tempHandle,mA);
			assert(nhNode->numBConnections == tempHandle.handleContents->numExits);
			n->numBConnections = tempHandle.handleContents->numExits;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDReturnMapHandle mB;
				for (unsigned int j = 0; j < nhNode->BConnection[i].returnMapHandle.Size(); j++)
					mB.AddToEnd(nhNode->BConnection[i].returnMapHandle[j]);
				mB.Canonicalize();
				CFLOBDDNodeHandle tempB = VectorToMatrixInterleavedNode(memoTable, *(nhNode->BConnection[i].entryPointHandle));
				n->BConnection[i] = Connection(tempB, mB);
			}
		}
		n->numExits = nhNode->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif

		CFLOBDDNodeHandle nHandle(n);
		memoTable->Insert(nh, nHandle);
		return nHandle;
	}

	CFLOBDDNodeHandle MatrixToVectorNode(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		assert(nh.handleContents->level > 0);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nh.handleContents->level - 1);
		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;

		/*for (unsigned int i = 0; i < nhNode->numBConnections; i++)
		{
			assert(nhNode->BConnection[i].entryPointHandle == CFLOBDDNodeHandle::NoDistinctionNode[nh.handleContents->level - 1]);
		}*/

		if (nh.handleContents->level == 1)
		{
			return *(nhNode->AConnection.entryPointHandle);
			
		}
		else
		{
			CFLOBDDNodeHandle tempHandle = MatrixToVectorNode(memoTable, *(nhNode->AConnection.entryPointHandle));
			n->AConnection = Connection(tempHandle, nhNode->AConnection.returnMapHandle);
			//assert(nhNode->numExits == tempHandle.handleContents->numExits);
			n->numBConnections = tempHandle.handleContents->numExits;
			n->BConnection = new Connection[n->numBConnections];
			for (int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDNodeHandle tmpB = MatrixToVectorNode(memoTable, *(nhNode->BConnection[i].entryPointHandle));
				n->BConnection[i] = Connection(tmpB, nhNode->BConnection[i].returnMapHandle);
			}
		}
		n->numExits = nhNode->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif

		CFLOBDDNodeHandle nHandle(n);
		memoTable->Insert(nh, nHandle);
		return nHandle;
	}

	CFLOBDDNodeHandle MkColumn1MatrixNode(unsigned int level)
	{
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(level);
		assert(level > 0);
		if (level == 1)
		{
			CFLOBDDReturnMapHandle m1, m2;
			m1.AddToEnd(0);
			m1.Canonicalize();

			m2.AddToEnd(0);
			m2.AddToEnd(1);
			m2.Canonicalize();

			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
			n->numBConnections = 1;
			n->BConnection = new Connection[1];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m2);
		}
		else
		{
			CFLOBDDReturnMapHandle m1, m2;
			m1.AddToEnd(0);
			m1.AddToEnd(1);
			m1.Canonicalize();

			m2.AddToEnd(1);
			m2.Canonicalize();

			CFLOBDDNodeHandle tempHandle = MkColumn1MatrixNode(level - 1);
			n->AConnection = Connection(tempHandle, m1);
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(tempHandle, m1);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1],m2);
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkVectorWithVoc12Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nhNode->level + 1);
		
		if (nhNode->level == 0)
		{
			CFLOBDDNodeHandle tempHandle(nh);
			CFLOBDDReturnMapHandle m1;
			for (int i = 0; i < tempHandle.handleContents->numExits; i++)
			{
				m1.AddToEnd(i);
			}
			m1.Canonicalize();
			n->AConnection = Connection(nh, m1);
			n->numBConnections = tempHandle.handleContents->numExits;
			n->BConnection = new Connection[n->numBConnections];
			for (int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDReturnMapHandle m2;
				m2.AddToEnd(i);
				m2.Canonicalize();
				n->BConnection[i] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[0],m2);
			}
		}
		else
		{
			CFLOBDDNodeHandle temp = MkVectorWithVoc12Node(memoTable, *(nhNode->AConnection.entryPointHandle));
			n->AConnection = Connection(temp, nhNode->AConnection.returnMapHandle);
			n->numBConnections = nhNode->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (int i = 0; i < n->numBConnections; i++){
				CFLOBDDNodeHandle tmpB = MkVectorWithVoc12Node(memoTable, *(nhNode->BConnection[i].entryPointHandle));
				n->BConnection[i] = Connection(tmpB, nhNode->BConnection[i].returnMapHandle);
			}
		}
		n->numExits = nhNode->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif

		CFLOBDDNodeHandle nHandle(n);
		memoTable->Insert(nh, nHandle);
		return nHandle;
	}

	CFLOBDDNodeHandle VectorShiftVocs1To2Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		//assert(nh.handleContents->level >= 1);
		//assert(nh.handleContents->NodeKind() == CFLOBDD_INTERNAL);

		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nhNode->level);
		

		if (nhNode->level == 1)
		{
			if (nh == CFLOBDDNodeHandle::NoDistinctionNode[1])
				return nhNode;

			for (int i = 0; i < nhNode->numBConnections; i++)
				assert(*(nhNode->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[0]);

			CFLOBDDReturnMapHandle m0;
			m0.AddToEnd(0);
			m0.Canonicalize();

			n->AConnection = Connection(CFLOBDDNodeHandle::NoDistinctionNode[0], m0);
			n->numBConnections = 1;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = nhNode->AConnection;
			n->numExits = nhNode->numBConnections;
		}
		else
		{
			CFLOBDDNodeHandle tmpA = VectorShiftVocs1To2Node(memoTable, *(nhNode->AConnection.entryPointHandle));
			n->AConnection = Connection(tmpA, nhNode->AConnection.returnMapHandle);
			n->numBConnections = nhNode->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int i = 0; i < n->numBConnections; i++) {
				CFLOBDDNodeHandle temp = VectorShiftVocs1To2Node(memoTable, *(nhNode->BConnection[i].entryPointHandle));
				n->BConnection[i] = Connection(temp, nhNode->BConnection[i].returnMapHandle);
			}
			n->numExits = nhNode->numExits;
		}
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		CFLOBDDNodeHandle nHandle(n);
		memoTable->Insert(nh, nHandle);
		return nHandle;
	}

	long double addNumPathsToExit(long double path1, long double path2){
		long double maxPath = std::max(path1, path2);
		long double minPath = std::min(path1, path2);
		if (minPath == -1 * std::numeric_limits<double>::infinity())
			return maxPath;
		long double totalPaths = minPath;
		if ((maxPath - minPath) > 15)
			totalPaths = maxPath;
		else{
			unsigned int powerVal = pow(2, maxPath - minPath) + 1;
			totalPaths += log2l(powerVal);
		}
		return totalPaths;
	}


	long double addNumPathsToExit(std::vector<long double>& logOfPaths){
		if (logOfPaths.size() == 1)
			return logOfPaths.back();
		long double sum = 0.0;
		for (int i = 0; i < logOfPaths.size() - 1; i++){
			if (logOfPaths[i] != -1.0 * std::numeric_limits<double>::infinity())
				sum += pow(2, logOfPaths[i] - logOfPaths.back());
		}
		long double logOfSum = logOfPaths.back() + log1p(sum) / ((double)log(2));
		return logOfSum;
	}

	BIG_FLOAT addNumPathsToExit(std::vector<BIG_FLOAT>& logOfPaths){
		if (logOfPaths.size() == 1)
			return logOfPaths.back();
		BIG_FLOAT sum = 0.0;
		for (int i = 0; i < logOfPaths.size() - 1; i++){
			if (logOfPaths[i] != -1.0 * std::numeric_limits<BIG_FLOAT>::infinity())
				sum += boost::multiprecision::pow(2, logOfPaths[i] - logOfPaths.back());
		}
		BIG_FLOAT logOfSum = logOfPaths.back() + (boost::multiprecision::log1p(sum) / ((double)log(2)));
		return logOfSum;
	}

	long double getLogSumNumPaths(std::vector<std::pair<long double, unsigned int>>& numBPaths, unsigned int size){
		std::vector<long double> paths;
		assert(numBPaths.size() >= size);
		for (int i = 0; i < size; i++)
			paths.push_back(numBPaths[i].first);
		return addNumPathsToExit(paths);
	}

	BIG_FLOAT getLogSumNumPaths(std::vector<std::pair<BIG_FLOAT, unsigned int>>& numBPaths, unsigned int size){
		std::vector<BIG_FLOAT> paths;
		assert(numBPaths.size() >= size);
		for (int i = 0; i < size; i++)
			paths.push_back(numBPaths[i].first);
		return addNumPathsToExit(paths);
	}

//#ifdef PATH_COUNTING_ENABLED
	std::pair<std::string,std::string> SamplingNode(CFLOBDDNodeHandle nh, unsigned int index, bool VocTwo)
	{
		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;

		if (nhNode->level == 0)
		{
			if (nh == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle)
			{
				return std::make_pair(std::to_string(index),"");
			}
			else
			{
				assert(index == 0);
				double random_value = ((double)rand() / (RAND_MAX));
				if (random_value < 0.5)
					return std::make_pair("0","");
				else
					return std::make_pair("1","");
			}
		}
		std::vector<std::pair<long double, unsigned int>> numBPaths;
		//std::vector<cpp_int> numBPaths;
		//std::vector<unsigned long long int> numBPaths;
		//cpp_int numBTotalPaths = 0;
		long double numBTotalPaths = 0.0;
		//unsigned long long int numBTotalPaths = 0;
		for (unsigned int i = 0; i < nhNode->numBConnections; i++)
		{
			int BIndex = nhNode->BConnection[i].returnMapHandle.LookupInv(index);
			if (BIndex == -1){
				numBPaths.push_back(std::make_pair(-1 * std::numeric_limits<long double>::infinity(), i));
			}
			else{
				numBPaths.push_back(std::make_pair(nhNode->BConnection[i].entryPointHandle->handleContents->numPathsToExit[BIndex] + 
					nhNode->AConnection.entryPointHandle->handleContents->numPathsToExit[i], i));
			}
				//numBPaths.push_back(nhNode->BConnection[i].entryPointHandle.handleContents->numPathsToExit[BIndex]);
		}
		sort(numBPaths.begin(), numBPaths.end(), sortNumPathPairs<long double>);
		/*std::cout << numBPaths.size() << std::endl;
		for (unsigned int i = 0; i < numBPaths.size(); i++)
		{
			std::cout << i << " " << numBPaths[i].first << " " << numBPaths[i].second << std::endl;
		}*/
		numBTotalPaths = getLogSumNumPaths(numBPaths, numBPaths.size());
		long double random_value = 0.0;
		if (numBTotalPaths >= 64){
			std::random_device rd;
			std::default_random_engine generator(rd());
			std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF);
			random_value = log2l(distribution(generator)) + numBTotalPaths - 64;
		}
		else{
			random_value = log2l((((double)rand()) / RAND_MAX)*pow(2, numBTotalPaths));
		}
		long double val = -1 * std::numeric_limits<long double>::infinity();
		/*cpp_int random_value = gen() % numBTotalPaths;
		cpp_int val = 0;*/
		/*unsigned long long int random_value = rand() % numBTotalPaths;
		unsigned long long int val = 0;*/
		int BConnectionIndex = -1;
		for (unsigned int i = 0; i < numBPaths.size(); i++)
		{
			if (numBPaths[i].first == -1 * std::numeric_limits<long double>::infinity())
				continue;
			val = getLogSumNumPaths(numBPaths, i+1);
			if (val >= random_value)
			{
				BConnectionIndex = numBPaths[i].second;
				break;
			}
		}
		//if (BConnectionIndex == -1){
		//	std::cout << val << " " << random_value << " " << index << " " << nhNode->level << std::endl;
		//	for (unsigned int i = 0; i < numBPaths.size(); i++)
		//	{
		//		std::cout << numBPaths[i].first << " " << numBPaths[i].second << std::endl;
		//	}
		//	//BConnectionIndex = numBPaths.back().second;
		//}
		assert(BConnectionIndex != -1);
		assert(nhNode->BConnection[BConnectionIndex].returnMapHandle.LookupInv(index) != -1);
		assert(BConnectionIndex < nhNode->numBConnections);
		assert(nhNode->BConnection[BConnectionIndex].returnMapHandle.LookupInv(index) < nhNode->BConnection[BConnectionIndex].returnMapHandle.mapContents->mapArray.size());
		std::pair<std::string,std::string> AString = SamplingNode(*(nhNode->AConnection.entryPointHandle), BConnectionIndex, VocTwo);
		std::pair<std::string, std::string> BString = SamplingNode(*(nhNode->BConnection[BConnectionIndex].entryPointHandle), nhNode->BConnection[BConnectionIndex].returnMapHandle.LookupInv(index), VocTwo);
		if (nhNode->level == 1)
			return std::make_pair(AString.first, BString.first);
		if (nhNode->level == 2 && !VocTwo)
			return std::make_pair(AString.first, BString.first);
		return std::make_pair(AString.first + BString.first, AString.second + BString.second);
		//return std::make_pair(AString.first + BString.first + AString.second + BString.second, "");
	}
//#endif
}  // namespace CFL_OBDD
