#include "../../cudd-3.0.0/cplusplus/cuddObj.hh"
#include <iostream>

#include <ctime>
#include <time.h>
#include <string>
#include <chrono>


using namespace std;
using namespace std::chrono;

BDD MkAnd(BDD a, BDD b)
{
	return a & b;
}

BDD MkNor(BDD a, BDD b)
{
	return ~(a | b);
}

BDD MkNand(BDD a, BDD b)
{
	return ~(a & b);
}

BDD MkExclusiveOr(BDD a, BDD b)
{
	return a ^ b;
}

BDD MkNot(BDD a)
{
	return !a;
}

static void c17(){

auto start = high_resolution_clock::now();
	Cudd mgr(0,0);
int coeff = 1;
int offset = 0;
// iscas85: c17 -----------------------------------------------------
// BDD gat1, gat2, gat3, gat6, gat7, gat10, gat11,
//         gat16, gat19, gat22, gat23;

int max_level = 3;

BDD gat1 = mgr.bddVar(0);
BDD gat2 = mgr.bddVar(1);
BDD gat3 = mgr.bddVar(2);
BDD gat6 = mgr.bddVar(3);
BDD gat7 = mgr.bddVar(4);
BDD gat10 = MkNand(gat1, gat3);
BDD gat11 = MkNand(gat3, gat6);
BDD gat16 = MkNand(gat2, gat11);
BDD gat19 = MkNand(gat11, gat7);
BDD gat22 = MkNand(gat10, gat16);
BDD gat23 = MkNand(gat16, gat19);

auto end = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(end - start);

unsigned int nodeCount = 0;
cout << "BDD sizes" << endl;
std::vector<BDD> bdds = {gat22, gat23};
nodeCount = mgr.nodeCount(bdds);
std::cout << "Duration: " << duration.count() << " Memory: " << nodeCount << std::endl;

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


int main () {

  c17();
  return 0;
}
