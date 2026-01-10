#include "../../cudd-3.0.0/cplusplus/cuddObj.hh"
#include <iostream>

#include <ctime>
#include <time.h>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

BDD MkOr(BDD a, BDD b)
{
	return a | b;
}

BDD MkAnd(BDD a, BDD b)
{
	return a & b;
}

BDD MkNor(BDD a, BDD b)
{
	return !(a | b);
}

BDD MkNand(BDD a, BDD b)
{
	return !(a & b);
}

BDD MkExclusiveOr(BDD a, BDD b)
{
	return a ^ b;
}

BDD MkNot(BDD a)
{
	return !a;
}

BDD MkIfThenElse(BDD f, BDD g, BDD h)
{
   return f.Ite(g, h);
}

static void iscasProduct(){

  auto start = high_resolution_clock::now();

  Cudd mgr(0,0);

  int max_level = 4;
// product -----------------------------------------------------
// BDD x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
//         x11, x12, x13, x14, x15;
BDD x0 = mgr.bddVar(0);
BDD x1 = mgr.bddVar(1);
BDD x2 = mgr.bddVar(2);
BDD x3 = mgr.bddVar(3);
BDD x4 = mgr.bddVar(4);
BDD x5 = mgr.bddVar(5);
BDD x6 = mgr.bddVar(6);
BDD x7 = mgr.bddVar(7);
BDD x8 = mgr.bddVar(8);
BDD x9 = mgr.bddVar(9);
BDD x10 = mgr.bddVar(10);
BDD x11 = mgr.bddVar(11);
BDD x12 = mgr.bddVar(12);
BDD x13 = mgr.bddVar(13);
BDD x14 = mgr.bddVar(14);
BDD x15 = mgr.bddVar(15);
BDD x0_1 = MkAnd(x0, x1);
BDD x2_3 = MkAnd(x2, x3);
BDD x4_5 = MkAnd(x4, x5);
BDD x6_7 = MkAnd(x6, x7);
BDD x8_9 = MkAnd(x8, x9);
BDD x10_11 = MkAnd(x10, x11);
BDD x12_13 = MkAnd(x12, x13);
BDD x14_15 = MkAnd(x14, x15);

BDD x0_3 = MkOr(x0_1, x2_3);
BDD x4_7 = MkOr(x4_5, x6_7);
BDD x8_11 = MkOr(x8_9, x10_11);
BDD x12_15 = MkOr(x12_13, x14_15);

BDD x0_7 = MkOr(x0_3, x4_7);
BDD x8_15 = MkOr(x8_11, x12_15);

BDD x0_15 = MkOr(x0_7, x8_15);

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(end - start); 

unsigned int nodeCount = 0;
cout << "BDD sizes" << endl;
x0_1.nodeCount();
x2_3.nodeCount();
x4_5.nodeCount();
x6_7.nodeCount();
x8_9.nodeCount();
x10_11.nodeCount();
x12_13.nodeCount();
x14_15.nodeCount();

x0_3.nodeCount();
x4_7.nodeCount();
x8_11.nodeCount();
x12_15.nodeCount();

x0_7.nodeCount();
x8_15.nodeCount();

x0_15.nodeCount();

std::cout << "Duration: " << duration.count() << " Memory: " << nodeCount << std::endl;

// if MaxLevel is 4, test all assignments
// if (BDD::maxLevel == 4) {

//   cout << "Testing all assignments" << endl;

//   unsigned int size = 1 << BDD::maxLevel;
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

int main () {

  iscasProduct();
  return 0;
}

