#ifndef VECTOR_NODE_GUARD
#define VECTOR_NODE_GUARD

#include <string.h>
#include "cflobdd_node.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
// #include "cflobdd_top_node_linear_map.h"

namespace CFL_OBDD {

	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
	//typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000> > BIG_FLOAT;

	extern CFLOBDDNodeHandle MkBasisVectorNode(unsigned int level, unsigned int index);
	extern CFLOBDDNodeHandle MkBasisVectorNode(unsigned int level, std::string s);
	
	extern void VectorInitializerNode();  // Empty for now

	extern CFLOBDDNodeHandle VectorToMatrixInterleavedNode(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle MatrixToVectorNode(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle MkColumn1MatrixNode(unsigned int level);
	extern CFLOBDDNodeHandle MkVectorWithVoc12Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle VectorShiftVocs1To2Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
//#ifdef PATH_COUNTING_ENABLED
	extern std::pair<std::string,std::string> SamplingNode(CFLOBDDNodeHandle nh, unsigned int index, bool voctwo = false);
//#endif
	// needs to be removed and linked to the one in cflobdd_node.cpp
	extern long double addNumPathsToExit(long double path1, long double path2);
	extern long double addNumPathsToExit(std::vector<long double>& paths);
	extern BIG_FLOAT addNumPathsToExit(std::vector<BIG_FLOAT>& paths);
	long double getLogSumNumPaths(std::vector<std::pair<long double, unsigned int>>& numBPaths, unsigned int size);
	BIG_FLOAT getLogSumNumPaths(std::vector<std::pair<BIG_FLOAT, unsigned int>>& numBPaths, unsigned int size);
	
	template <typename T>
	bool sortNumPathPairs(const std::pair<T, unsigned int>& p1, const std::pair<T, unsigned int> &p2){
		if (p1.first < p2.first)
			return true;
		else if (p1.first > p2.first)
			return false;
		return p1.second < p2.second;
	}
}

#endif

