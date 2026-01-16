#include "../cflobdd_int.h"
#include "../cflobdd_top_node_int.h"
#include <ctime>
#include <time.h>
#include <string>
#include <chrono>

using namespace std;
using namespace CFL_OBDD;
using namespace std::chrono;

static void c17(){
// iscas85: c17 -----------------------------------------------------
// CFLOBDD gat1, gat2, gat3, gat6, gat7, gat10, gat11,
//         gat16, gat19, gat22, gat23;

auto start = high_resolution_clock::now();

int max_level = 3;

CFLOBDD gat1 = MkProjection(0, max_level);
CFLOBDD gat2 = MkProjection(1, max_level);
CFLOBDD gat3 = MkProjection(2, max_level);
CFLOBDD gat6 = MkProjection(3, max_level);
CFLOBDD gat7 = MkProjection(4, max_level);
CFLOBDD gat10 = MkNand(gat1, gat3);
CFLOBDD gat11 = MkNand(gat3, gat6);
CFLOBDD gat16 = MkNand(gat2, gat11);
CFLOBDD gat19 = MkNand(gat11, gat7);
CFLOBDD gat22 = MkNand(gat10, gat16);
CFLOBDD gat23 = MkNand(gat16, gat19);

auto end = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(end - start);

unsigned int nodeCount, edgeCount;
// cout << "CFLOBDD sizes" << endl;
GroupCountNodesAndEdgesStart(nodeCount, edgeCount);
gat22.GroupCountNodesAndEdges(nodeCount, edgeCount);
// cout << nodeCount << ", " << edgeCount << ", " << (nodeCount + edgeCount) << endl;
gat23.GroupCountNodesAndEdges(nodeCount, edgeCount);
// cout << nodeCount << ", " << edgeCount << ", " << (nodeCount + edgeCount) << endl;
GroupCountNodesAndEdgesEnd();
// cout << nodeCount << ", " << edgeCount << ", " << (nodeCount + edgeCount) << endl;

std::cout << "Duration: " << duration.count()/1000.0 << " Memory: " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;

// gat23.print(std::cout);

// nodeCount = 0;
// edgeCount = 0;
// unsigned int returnEdgesCount = 0, returnEdgesObjCount = 0;
// gat1.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat1| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat2.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat2| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat3.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat3| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat6.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat6| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat7.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat7| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat10.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat10| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat11.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat11| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat16.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat16| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat19.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat19| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat22.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat22| = " << nodeCount << ", " << edgeCount << std::endl;
// nodeCount = 0;
// edgeCount = 0;
// returnEdgesCount = 0;
// returnEdgesObjCount = 0;
// gat23.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|gat23| = " << nodeCount << ", " << edgeCount << std::endl;

// if MaxLevel is 3 or 4, test all assignments
// if (max_level == 3 || max_level == 4) {

//   cout << "Testing all assignments" << endl;

//   unsigned int size = 1 << max_level;
//   SH_OBDD::Assignment a(size);
//   bool val22, val23;
//   unsigned long int range = 1UL << size;
//   for (unsigned long int i = 0UL; i < range; i++) {
//     unsigned long int mask = 1UL;
//     for (int j = size - 1; j >= 0; j--) {
//       a[j] = (i & mask);
//       mask = mask << 1;
//     }
//     val22 = gat22.Evaluate(a);
//     val23 = gat23.Evaluate(a);
  
//     bool temp1, temp2, temp3, temp6, temp7, temp10, temp11, temp16, temp19, temp22, temp23;
  
//     temp1 = (bool)a[0];
//     temp2 = (bool)a[1];
//     temp3 = (bool)a[2];
//     temp6 = (bool)a[3];
//     temp7 = (bool)a[4];
//     temp10 = !(temp1 && temp3);
//     temp11 = !(temp3 && temp6);
//     temp16 = !(temp2 && temp11);
//     temp19 = !(temp11 && temp7);
//     temp22 = !(temp10 && temp16);
//     temp23 = !(temp16 && temp19);
  
//     if (val22 != temp22) {
//       cerr << "val22: " << val22 << " != " << temp22 << endl;
//     }
//     if (val23 != temp23) {
//       cerr << "val23: " << val23 << " != " << temp23 << endl;
//     }
//   }
// }
// else {
//   cout << "Cannot test all assignments: maxLevel must be 3 or 4" << endl;
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

static void ClearModules()
{
	DisposeOfTripleProductCache();
	DisposeOfPairProductCache();
	CFLOBDDNodeHandle::DisposeOfReduceCache();
}

int main () {

  InitModules();
  c17();
  // ClearModules();
  return 0;
}
