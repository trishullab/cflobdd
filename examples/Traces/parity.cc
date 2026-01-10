#include "../../cudd-3.0.0/cplusplus/cuddObj.hh"
#include <iostream>

#include <ctime>
#include <time.h>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

ADD MkExclusiveOr(ADD a, ADD b)
{
	return a.Xor(b);
}

ADD MkIfThenElse(ADD f, ADD g, ADD h)
{
   return f.Ite(g, h);
}

static void iscasParity(){

// if MaxLevel is 4, run test
// if (ADD::maxLevel == 4) {

auto start = high_resolution_clock::now();

Cudd mgr(0,0);
  int max_level = 4;
  // ADD x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
  //         x11, x12, x13, x14, x15;

  ADD x0 = mgr.addVar(0);
  ADD x1 = mgr.addVar(1);
  ADD x2 = mgr.addVar(2);
  ADD x3 = mgr.addVar(3);
  ADD x4 = mgr.addVar(4);
  ADD x5 = mgr.addVar(5);
  ADD x6 = mgr.addVar(6);
  ADD x7 = mgr.addVar(7);
  ADD x8 = mgr.addVar(8);
  ADD x9 = mgr.addVar(9);
  ADD x10 = mgr.addVar(10);
  ADD x11 = mgr.addVar(11);
  ADD x12 = mgr.addVar(12);
  ADD x13 = mgr.addVar(13);
  ADD x14 = mgr.addVar(14);
  ADD x15 = mgr.addVar(15);

  ADD parity16;
  std::vector<ADD> x_s, y_s;
  x_s = {x0, x1, x2, x3, x4, x5, x6, x7};
  y_s = {x8, x9, x10, x11, x12, x13, x14, x15};

  parity16 = mgr.Xneqy(x_s, y_s);
  
  ADD x0_1 = MkExclusiveOr(x0, x1);
  ADD x2_3 = MkExclusiveOr(x2, x3);
  ADD x4_5 = MkExclusiveOr(x4, x5);
  ADD x6_7 = MkExclusiveOr(x6, x7);
  ADD x8_9 = MkExclusiveOr(x8, x9);
  ADD x10_11 = MkExclusiveOr(x10, x11);
  ADD x12_13 = MkExclusiveOr(x12, x13);
  ADD x14_15 = MkExclusiveOr(x14, x15);
  
  ADD x0_3 = MkExclusiveOr(x0_1, x2_3);
  ADD x4_7 = MkExclusiveOr(x4_5, x6_7);
  ADD x8_11 = MkExclusiveOr(x8_9, x10_11);
  ADD x12_15 = MkExclusiveOr(x12_13, x14_15);
  
  ADD x0_7 = MkExclusiveOr(x0_3, x4_7);
  ADD x8_15 = MkExclusiveOr(x8_11, x12_15);
  
  ADD x0_15 = MkExclusiveOr(x0_7, x8_15);
  
  if (x0_15 == parity16) {
    cout << "x0_15 == parity16" << endl;
  }
  else {
    cout << "x0_15 != parity16" << endl;
  }

  auto end = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(end - start); 
  
  unsigned int nodeCount = 0;
  nodeCount = x0_15.nodeCount();
  std::cout << "Duration: " << duration.count() << " Memory: " << nodeCount << std::endl;
  
  // cout << "Testing all assignments" << endl;

  // unsigned int size = 1 << ADD::maxLevel;
  // Assignment a(size);
  // bool val0_15;
  // unsigned long int range = 1UL << size;
  // for (unsigned long int i = 0UL; i < range; i++) {
  //   unsigned long int mask = 1UL;
  //   bool bb = false;
  //   for (int j = size - 1; j >= 0; j--) {
  //     a[j] = (i & mask);
  //     bb ^= a[j];
  //     mask = mask << 1;
  //   }
  //   val0_15 = x0_15.Evaluate(a);
  //   if (val0_15 != bb) {
  //     cerr << a << ": " << val0_15 << endl;
  //   }
  // }
// }
// else {
//   cout << "Cannot run test: maxLevel must be 4" << endl;
// }
}

int main () {
  iscasParity();
  return 0;
}
