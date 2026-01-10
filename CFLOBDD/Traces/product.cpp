#include "../cflobdd_int.h"
#include "../cflobdd_top_node_int.h"

#include <ctime>
#include <time.h>
#include <string>
#include <chrono>

using namespace std;
using namespace CFL_OBDD;
using namespace std::chrono;
using namespace SH_OBDD;
static void iscasProduct(){

  auto start = high_resolution_clock::now();

  int max_level = 4;
// product -----------------------------------------------------
// CFLOBDD x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
//         x11, x12, x13, x14, x15;
CFLOBDD x0 = MkProjection(0, max_level);
CFLOBDD x1 = MkProjection(1, max_level);
CFLOBDD x2 = MkProjection(2, max_level);
CFLOBDD x3 = MkProjection(3, max_level);
CFLOBDD x4 = MkProjection(4, max_level);
CFLOBDD x5 = MkProjection(5, max_level);
CFLOBDD x6 = MkProjection(6, max_level);
CFLOBDD x7 = MkProjection(7, max_level);
CFLOBDD x8 = MkProjection(8, max_level);
CFLOBDD x9 = MkProjection(9, max_level);
CFLOBDD x10 = MkProjection(10, max_level);
CFLOBDD x11 = MkProjection(11, max_level);
CFLOBDD x12 = MkProjection(12, max_level);
CFLOBDD x13 = MkProjection(13, max_level);
CFLOBDD x14 = MkProjection(14, max_level);
CFLOBDD x15 = MkProjection(15, max_level);
CFLOBDD x0_1 = MkAnd(x0, x1);
CFLOBDD x2_3 = MkAnd(x2, x3);
CFLOBDD x4_5 = MkAnd(x4, x5);
CFLOBDD x6_7 = MkAnd(x6, x7);
CFLOBDD x8_9 = MkAnd(x8, x9);
CFLOBDD x10_11 = MkAnd(x10, x11);
CFLOBDD x12_13 = MkAnd(x12, x13);
CFLOBDD x14_15 = MkAnd(x14, x15);

CFLOBDD x0_3 = MkOr(x0_1, x2_3);
CFLOBDD x4_7 = MkOr(x4_5, x6_7);
CFLOBDD x8_11 = MkOr(x8_9, x10_11);
CFLOBDD x12_15 = MkOr(x12_13, x14_15);

CFLOBDD x0_7 = MkOr(x0_3, x4_7);
CFLOBDD x8_15 = MkOr(x8_11, x12_15);

CFLOBDD x0_15 = MkOr(x0_7, x8_15);

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(end - start); 

unsigned int nodeCount, edgeCount;
cout << "CFLOBDD sizes" << endl;
GroupCountNodesAndEdgesStart(nodeCount, edgeCount);
x0_1.GroupCountNodesAndEdges(nodeCount, edgeCount);
x2_3.GroupCountNodesAndEdges(nodeCount, edgeCount);
x4_5.GroupCountNodesAndEdges(nodeCount, edgeCount);
x6_7.GroupCountNodesAndEdges(nodeCount, edgeCount);
x8_9.GroupCountNodesAndEdges(nodeCount, edgeCount);
x10_11.GroupCountNodesAndEdges(nodeCount, edgeCount);
x12_13.GroupCountNodesAndEdges(nodeCount, edgeCount);
x14_15.GroupCountNodesAndEdges(nodeCount, edgeCount);

x0_3.GroupCountNodesAndEdges(nodeCount, edgeCount);
x4_7.GroupCountNodesAndEdges(nodeCount, edgeCount);
x8_11.GroupCountNodesAndEdges(nodeCount, edgeCount);
x12_15.GroupCountNodesAndEdges(nodeCount, edgeCount);

x0_7.GroupCountNodesAndEdges(nodeCount, edgeCount);
x8_15.GroupCountNodesAndEdges(nodeCount, edgeCount);

x0_15.GroupCountNodesAndEdges(nodeCount, edgeCount);
GroupCountNodesAndEdgesEnd();
// cout << nodeCount << ", " << edgeCount << endl;

// unsigned int returnEdgesCount, returnEdgesObjCount;
// x0_15.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
  std::cout << "Duration: " << duration.count() << " Memory: " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;

// if MaxLevel is 4, test all assignments
// if (CFLOBDD::maxLevel == 4) {

//   cout << "Testing all assignments" << endl;

//   unsigned int size = 1 << CFLOBDD::maxLevel;
//   Assignment a(size);
//   bool val0_15;
//   unsigned long int range = 1UL << size;
//   for (unsigned long int i = 0UL; i < range; i++) {
//     unsigned long int mask = 1UL;
//     for (int j = size - 1; j >= 0; j--) {
//       a[j] = (i & mask);
//       mask = mask << 1;
//     }
//     val0_15 = x0_15.Evaluate(a);

//     bool temp =    ((bool)a[0] && (bool)a[1])
//                 || ((bool)a[2] && (bool)a[3])
//                 || ((bool)a[4] && (bool)a[5])
//                 || ((bool)a[6] && (bool)a[7])
//                 || ((bool)a[8] && (bool)a[9])
//                 || ((bool)a[10] && (bool)a[11])
//                 || ((bool)a[12] && (bool)a[13])
//                 || ((bool)a[14] && (bool)a[15]);
//     if (val0_15 != temp) {
//       cerr << "val0_15: " << val0_15 << " != " << temp << endl;
//     }
//   }
// }
// else {
//   cout << "Cannot test all assignments: maxLevel must be 4" << endl;
// }
}

static void InitModules()
{

	CFLOBDDNodeHandle::InitNoDistinctionTable();
	// CFLOBDDNodeHandle::InitAdditionInterleavedTable();
	CFLOBDDNodeHandle::InitReduceCache();
	InitPairProductCache();
	InitTripleProductCache();
}

// static void ClearModules()
// {
// 	DisposeOfTripleProductCache();
// 	DisposeOfPairProductCache();
// 	CFLOBDDNodeHandle::DisposeOfReduceCache();
// }

int main () {

  InitModules();
  iscasProduct();
//   ClearModules();
  return 0;
}

