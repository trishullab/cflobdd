#include "../cflobdd_int.h"
#include "../cflobdd_top_node_int.h"

#include <ctime>
#include <time.h>
#include <string>
#include <chrono>

using namespace std;
using namespace CFL_OBDD;
using namespace std::chrono;

// iscas85: c6288-8
static void c6288_8()
{


auto start = high_resolution_clock::now();

int coeff = 1;
int offset = 0;
int max_level = 4;

// #define DESCENDING_NONINTERLEAVED 0
// #define DESCENDING_INTERLEAVED 1
// #define ASCENDING_NONINTERLEAVED 1
// #define ASCENDING_INTERLEAVED 1
// #define AMANO_8BIT_OPTIMAL 1
// #define NEW_ORDERING 1
#define NEW_ORDERING_REV 1

#ifdef DESCENDING_NONINTERLEAVED
CFLOBDD a7 = MkProjection(coeff*0 + offset, max_level);
CFLOBDD a6 = MkProjection(coeff*1 + offset, max_level);
CFLOBDD a5 = MkProjection(coeff*2 + offset, max_level);
CFLOBDD a4 = MkProjection(coeff*3 + offset, max_level);
CFLOBDD a3 = MkProjection(coeff*4 + offset, max_level);
CFLOBDD a2 = MkProjection(coeff*5 + offset, max_level);
CFLOBDD a1 = MkProjection(coeff*6 + offset, max_level);
CFLOBDD a0 = MkProjection(coeff*7 + offset, max_level);
CFLOBDD b7 = MkProjection(coeff*8 + offset, max_level);
CFLOBDD b6 = MkProjection(coeff*9 + offset, max_level);
CFLOBDD b5 = MkProjection(coeff*10 + offset, max_level);
CFLOBDD b4 = MkProjection(coeff*11 + offset, max_level);
CFLOBDD b3 = MkProjection(coeff*12 + offset, max_level);
CFLOBDD b2 = MkProjection(coeff*13 + offset, max_level);
CFLOBDD b1 = MkProjection(coeff*14 + offset, max_level);
CFLOBDD b0 = MkProjection(coeff*15 + offset, max_level);
#endif

#ifdef DESCENDING_INTERLEAVED
CFLOBDD a7 = MkProjection(coeff*0 + offset, max_level);
CFLOBDD b7 = MkProjection(coeff*1 + offset, max_level);
CFLOBDD a6 = MkProjection(coeff*2 + offset, max_level);
CFLOBDD b6 = MkProjection(coeff*3 + offset, max_level);
CFLOBDD a5 = MkProjection(coeff*4 + offset, max_level);
CFLOBDD b5 = MkProjection(coeff*5 + offset, max_level);
CFLOBDD a4 = MkProjection(coeff*6 + offset, max_level);
CFLOBDD b4 = MkProjection(coeff*7 + offset, max_level);
CFLOBDD a3 = MkProjection(coeff*8 + offset, max_level);
CFLOBDD b3 = MkProjection(coeff*9 + offset, max_level);
CFLOBDD a2 = MkProjection(coeff*10 + offset, max_level);
CFLOBDD b2 = MkProjection(coeff*11 + offset, max_level);
CFLOBDD a1 = MkProjection(coeff*12 + offset, max_level);
CFLOBDD b1 = MkProjection(coeff*13 + offset, max_level);
CFLOBDD a0 = MkProjection(coeff*14 + offset, max_level);
CFLOBDD b0 = MkProjection(coeff*15 + offset, max_level);
#endif

#ifdef ASCENDING_NONINTERLEAVED
CFLOBDD a0 = MkProjection(coeff*0 + offset, max_level);
CFLOBDD a1 = MkProjection(coeff*1 + offset, max_level);
CFLOBDD a2 = MkProjection(coeff*2 + offset, max_level);
CFLOBDD a3 = MkProjection(coeff*3 + offset, max_level);
CFLOBDD a4 = MkProjection(coeff*4 + offset, max_level);
CFLOBDD a5 = MkProjection(coeff*5 + offset, max_level);
CFLOBDD a6 = MkProjection(coeff*6 + offset, max_level);
CFLOBDD a7 = MkProjection(coeff*7 + offset, max_level);
CFLOBDD b0 = MkProjection(coeff*8 + offset, max_level);
CFLOBDD b1 = MkProjection(coeff*9 + offset, max_level);
CFLOBDD b2 = MkProjection(coeff*10 + offset, max_level);
CFLOBDD b3 = MkProjection(coeff*11 + offset, max_level);
CFLOBDD b4 = MkProjection(coeff*12 + offset, max_level);
CFLOBDD b5 = MkProjection(coeff*13 + offset, max_level);
CFLOBDD b6 = MkProjection(coeff*14 + offset, max_level);
CFLOBDD b7 = MkProjection(coeff*15 + offset, max_level);
#endif

#ifdef ASCENDING_INTERLEAVED
CFLOBDD a0 = MkProjection(coeff*0 + offset, max_level);
CFLOBDD b0 = MkProjection(coeff*1 + offset, max_level);
CFLOBDD a1 = MkProjection(coeff*2 + offset, max_level);
CFLOBDD b1 = MkProjection(coeff*3 + offset, max_level);
CFLOBDD a2 = MkProjection(coeff*4 + offset, max_level);
CFLOBDD b2 = MkProjection(coeff*5 + offset, max_level);
CFLOBDD a3 = MkProjection(coeff*6 + offset, max_level);
CFLOBDD b3 = MkProjection(coeff*7 + offset, max_level);
CFLOBDD a4 = MkProjection(coeff*8 + offset, max_level);
CFLOBDD b4 = MkProjection(coeff*9 + offset, max_level);
CFLOBDD a5 = MkProjection(coeff*10 + offset, max_level);
CFLOBDD b5 = MkProjection(coeff*11 + offset, max_level);
CFLOBDD a6 = MkProjection(coeff*12 + offset, max_level);
CFLOBDD b6 = MkProjection(coeff*13 + offset, max_level);
CFLOBDD a7 = MkProjection(coeff*14 + offset, max_level);
CFLOBDD b7 = MkProjection(coeff*15 + offset, max_level);
#endif

#ifdef DESCENDING_ASCENDING
CFLOBDD a7 = MkProjection(coeff*0 + offset, max_level);
CFLOBDD a6 = MkProjection(coeff*1 + offset, max_level);
CFLOBDD a5 = MkProjection(coeff*2 + offset, max_level);
CFLOBDD a4 = MkProjection(coeff*3 + offset, max_level);
CFLOBDD a3 = MkProjection(coeff*4 + offset, max_level);
CFLOBDD a2 = MkProjection(coeff*5 + offset, max_level);
CFLOBDD a1 = MkProjection(coeff*6 + offset, max_level);
CFLOBDD a0 = MkProjection(coeff*7 + offset, max_level);
CFLOBDD b0 = MkProjection(coeff*8 + offset, max_level);
CFLOBDD b1 = MkProjection(coeff*9 + offset, max_level);
CFLOBDD b2 = MkProjection(coeff*10 + offset, max_level);
CFLOBDD b3 = MkProjection(coeff*11 + offset, max_level);
CFLOBDD b4 = MkProjection(coeff*12 + offset, max_level);
CFLOBDD b5 = MkProjection(coeff*13 + offset, max_level);
CFLOBDD b6 = MkProjection(coeff*14 + offset, max_level);
CFLOBDD b7 = MkProjection(coeff*15 + offset, max_level);
#endif

#ifdef ASCENDING_DESCENDING
CFLOBDD a0 = MkProjection(coeff*0 + offset, max_level);
CFLOBDD a1 = MkProjection(coeff*1 + offset, max_level);
CFLOBDD a2 = MkProjection(coeff*2 + offset, max_level);
CFLOBDD a3 = MkProjection(coeff*3 + offset, max_level);
CFLOBDD a4 = MkProjection(coeff*4 + offset, max_level);
CFLOBDD a5 = MkProjection(coeff*5 + offset, max_level);
CFLOBDD a6 = MkProjection(coeff*6 + offset, max_level);
CFLOBDD a7 = MkProjection(coeff*7 + offset, max_level);
CFLOBDD b7 = MkProjection(coeff*8 + offset, max_level);
CFLOBDD b6 = MkProjection(coeff*9 + offset, max_level);
CFLOBDD b5 = MkProjection(coeff*10 + offset, max_level);
CFLOBDD b4 = MkProjection(coeff*11 + offset, max_level);
CFLOBDD b3 = MkProjection(coeff*12 + offset, max_level);
CFLOBDD b2 = MkProjection(coeff*13 + offset, max_level);
CFLOBDD b1 = MkProjection(coeff*14 + offset, max_level);
CFLOBDD b0 = MkProjection(coeff*15 + offset, max_level);
#endif

// Amano & Maruoka, DAC 2007: the variable order
// that gives the smallest 8-bit multiplication BDD for MUL_7,8 is
// (a1, a2, a3, a4, b3, b4, b2, a5, b5, b1, a6, b6, a0, b7, a7, b0)
#ifdef AMANO_8BIT_OPTIMAL
CFLOBDD a1 = MkProjection(coeff*0 + offset, max_level);
CFLOBDD a2 = MkProjection(coeff*1 + offset, max_level);
CFLOBDD a3 = MkProjection(coeff*2 + offset, max_level);
CFLOBDD a4 = MkProjection(coeff*3 + offset, max_level);
CFLOBDD b3 = MkProjection(coeff*4 + offset, max_level);
CFLOBDD b4 = MkProjection(coeff*5 + offset, max_level);
CFLOBDD b2 = MkProjection(coeff*6 + offset, max_level);
CFLOBDD a5 = MkProjection(coeff*7 + offset, max_level);
CFLOBDD b5 = MkProjection(coeff*8 + offset, max_level);
CFLOBDD b1 = MkProjection(coeff*9 + offset, max_level);
CFLOBDD a6 = MkProjection(coeff*10 + offset, max_level);
CFLOBDD b6 = MkProjection(coeff*11 + offset, max_level);
CFLOBDD a0 = MkProjection(coeff*12 + offset, max_level);
CFLOBDD b7 = MkProjection(coeff*13 + offset, max_level);
CFLOBDD a7 = MkProjection(coeff*14 + offset, max_level);
CFLOBDD b0 = MkProjection(coeff*15 + offset, max_level);
#endif

#ifdef NEW_ORDERING
CFLOBDD a0 = MkProjection(0, 16);//32768
CFLOBDD b0 = MkProjection(32768, 16);//16384
CFLOBDD a1 = MkProjection(49152, 16);//8192
CFLOBDD b1 = MkProjection(57344, 16);//4096
CFLOBDD a2 = MkProjection(61440, 16);//2048
CFLOBDD b2 = MkProjection(63488, 16);//1024
CFLOBDD a3 = MkProjection(64512, 16);//512
CFLOBDD b3 = MkProjection(65024, 16);//256
CFLOBDD a4 = MkProjection(65280, 16);//128
CFLOBDD b4 = MkProjection(65408, 16);//64
CFLOBDD a5 = MkProjection(65472, 16);//32
CFLOBDD b5 = MkProjection(65504, 16);//16
CFLOBDD a6 = MkProjection(65520, 16);//8
CFLOBDD b6 = MkProjection(65528, 16);//4
CFLOBDD a7 = MkProjection(65532, 16);//2
CFLOBDD b7 = MkProjection(65534, 16);
#endif

// #ifdef NEW_ORDERING_REV
// CFLOBDD a0 = MkProjection(0, 16);//32768
// CFLOBDD b0 = MkProjection(2, 16);//16384
// CFLOBDD a1 = MkProjection(4, 16);//8192
// CFLOBDD b1 = MkProjection(8, 16);//4096
// CFLOBDD a2 = MkProjection(16, 16);//2048
// CFLOBDD b2 = MkProjection(32, 16);//1024
// CFLOBDD a3 = MkProjection(64, 16);//512
// CFLOBDD b3 = MkProjection(128, 16);//256
// CFLOBDD a4 = MkProjection(256, 16);//128
// CFLOBDD b4 = MkProjection(512, 16);//64
// CFLOBDD a5 = MkProjection(1024, 16);//32
// CFLOBDD b5 = MkProjection(2048, 16);//16
// CFLOBDD a6 = MkProjection(4096, 16);//8
// CFLOBDD b6 = MkProjection(8192, 16);//4
// CFLOBDD a7 = MkProjection(16384, 16);//2
// CFLOBDD b7 = MkProjection(32768, 16);
// #endif

#ifdef NEW_ORDERING_REV
CFLOBDD a0 = MkProjection(0, 3);//32768
CFLOBDD b0 = MkProjection(1, 3);//16384
CFLOBDD a1 = MkProjection(2, 3);//8192
CFLOBDD b1 = MkProjection(4, 3);//4096
#endif


std::cout << "	pp0_0	 " << std::endl;
CFLOBDD	pp0_0	 = MkAnd(a0, b0); // a0 & b0
CFLOBDD	pp1_0	 = MkAnd(a0, b1);
CFLOBDD	sum0	= pp0_0;
CFLOBDD	pp0_1	 = MkAnd(a1, b0);
std::cout << "	pp1_1	 " << std::endl;
CFLOBDD	pp1_1	 = MkAnd(a1, b1);
CFLOBDD	gate0_1_0	=  MkNor(pp1_0, pp0_1);
CFLOBDD	gate0_1_1	=  MkNor(pp1_0, gate0_1_0);
CFLOBDD	gate0_1_2	=  MkNor(pp0_1, gate0_1_0);
CFLOBDD	gate0_1_3	=  MkNor(gate0_1_1, gate0_1_2);
CFLOBDD	gate0_1_4	= MkNot(gate0_1_3); 
CFLOBDD	gate0_1_5	= MkNor(gate0_1_4, gate0_1_3); 
CFLOBDD	gate0_1_6	= MkNot(gate0_1_4); 
CFLOBDD	sum1	= MkNor(gate0_1_5,gate0_1_6); 
CFLOBDD	c0_1	= MkNor(gate0_1_0,gate0_1_4); 
CFLOBDD	gate1_1_0	=  MkNor(c0_1, pp1_1);
CFLOBDD	gate1_1_1	=  MkNor(c0_1, gate1_1_0);
CFLOBDD	gate1_1_2	=  MkNor(pp1_1, gate1_1_0);
CFLOBDD	gate1_1_3	=  MkNor(gate1_1_1, gate1_1_2);
std::cout << "	gate1_1_4	" << std::endl;
CFLOBDD	gate1_1_4	= MkNot(gate1_1_3); 
CFLOBDD	gate1_1_5	= MkNor(gate1_1_4, gate1_1_3); 
CFLOBDD	gate1_1_6	= MkNot(gate1_1_4); 
CFLOBDD	sum2	= MkNor(gate1_1_5,gate1_1_6); 
CFLOBDD	c0_2	= MkNor(gate1_1_0,gate1_1_4); 

auto end = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(end - start);   

std::cout << sum0 << std::endl;
std::cout << sum1 << std::endl;
std::cout << sum2 << std::endl;


unsigned int nodeCount, edgeCount;
std::cout << "CFLOBDD sizes" << std::endl;
GroupCountNodesAndEdgesStart(nodeCount, edgeCount);
sum0.GroupCountNodesAndEdges(nodeCount, edgeCount);
sum1.GroupCountNodesAndEdges(nodeCount, edgeCount);
sum2.GroupCountNodesAndEdges(nodeCount, edgeCount);
GroupCountNodesAndEdgesEnd();
std::cout << "Duration: " << duration.count() << " Memory: " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;

// nodeCount = 0;
// edgeCount = 0;
// unsigned int returnEdgesCount = 0, returnEdgesObjCount = 0;
// std::cout << "sum6 size -------------------------" << std::endl;
// sum6.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
// std::cout << "|sum6| = " << nodeCount << ", " << edgeCount << std::endl;

/*
// Test a specific multiplication problem
// Only for the order a7 .. a0 b7 .. b0 in positions 0 .. 15
unsigned int a = 238; //  0356: 0 <= a <= 255
unsigned int b = 165; //  0245: 0 <= b <= 255
Assignment assign(16);
unsigned long int mask = 1UL;
for (int j = 15; j >= 8; j--) {
  assign[j] = (b & mask);
  mask = mask << 1;
}
mask = 1UL;
for (int j = 7; j >= 0; j--) {
  assign[j] = (a & mask);
  mask = mask << 1;
}
bool val[16];

val[0] = sum0.Evaluate(assign);
val[1] = sum1.Evaluate(assign);
val[2] = sum2.Evaluate(assign);
val[3] = sum3.Evaluate(assign);
val[4] = sum4.Evaluate(assign);
val[5] = sum5.Evaluate(assign);
val[6] = sum6.Evaluate(assign);
val[7] = sum7.Evaluate(assign);
val[8] = sum8.Evaluate(assign);
val[9] = sum9.Evaluate(assign);
val[10] = sum10.Evaluate(assign);
val[11] = sum11.Evaluate(assign);
val[12] = sum12.Evaluate(assign);
val[13] = sum13.Evaluate(assign);
val[14] = sum14.Evaluate(assign);
val[15] = sum15.Evaluate(assign);

// 238 * 165 = 0356 * 0245 = 0x9966 = 1001 1001 0110 0110
for (int k = 15; k >= 0; k--) {
  std::cout << (val[k]?1:0);
}
std::cout << std::endl;
std::cout << a*b << std::endl;
*/
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
  c6288_8();
//   ClearModules();
  return 0;
}
