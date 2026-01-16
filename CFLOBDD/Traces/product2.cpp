#include "../cflobdd_int.h"
#include "../cflobdd_top_node_int.h"

#include <ctime>
#include <time.h>
#include <string>
#include <chrono>

using namespace std;
using namespace CFL_OBDD;
using namespace std::chrono;

static void iscasProduct2(){
// product -----------------------------------------------------
// if (CFLOBDD::maxLevel >= 27) {
//   cerr << "Cannot create product CFLOBDD: maxLevel must be < 27" << endl;
// }
CFLOBDD **x = new CFLOBDD* [CFLOBDD::maxLevel+1];
for (int i = CFLOBDD::maxLevel; i >= 0; i--) {
  x[i] = new CFLOBDD[1 << i];
  for (int j = 0; j < (1 << i); j++) {
    if (i == (int)CFLOBDD::maxLevel) {
      x[i][j] = MkProjection(j);
    }
    else if (i == (int)(CFLOBDD::maxLevel-1)) {
      x[i][j] = MkAnd(x[i+1][2*j],
                      x[i+1][2*j+1]);
    }
    else {
      x[i][j] = MkOr(x[i+1][2*j],
                     x[i+1][2*j+1]);
    }
  }
}

unsigned int nodeCount, edgeCount;
unsigned int returnEdgesCount, returnEdgesObjCount;
// cout << "CFLOBDD sizes" << endl;
GroupCountNodesAndEdgesStart(nodeCount, edgeCount);
for (int i = CFLOBDD::maxLevel; i >= 0; i--) {
  for (int j = 0; j < (1 << i); j++) {
    x[i][j].GroupCountNodesAndEdges(nodeCount, edgeCount);
  }
}
GroupCountNodesAndEdgesEnd();
// cout << nodeCount << ", " << edgeCount << endl;

x[0][0].CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
cout << nodeCount << ", " << edgeCount << endl;

if (CFLOBDD::maxLevel == 4) {
  CFLOBDD x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
          x11, x12, x13, x14, x15;
  x0 = MkProjection(0);
  x1 = MkProjection(1);
  x2 = MkProjection(2);
  x3 = MkProjection(3);
  x4 = MkProjection(4);
  x5 = MkProjection(5);
  x6 = MkProjection(6);
  x7 = MkProjection(7);
  x8 = MkProjection(8);
  x9 = MkProjection(9);
  x10 = MkProjection(10);
  x11 = MkProjection(11);
  x12 = MkProjection(12);
  x13 = MkProjection(13);
  x14 = MkProjection(14);
  x15 = MkProjection(15);
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

  if (x[0][0] == x0_15) {
    cout << "Equal" << endl;
  }
  else {
    cout << "Not equal" << endl;
  }
}
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
  iscasProduct2();
//   ClearModules();
  return 0;
}
