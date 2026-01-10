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

// iscas85: c6288-8
static void c6288_8()
{


auto start = high_resolution_clock::now();

Cudd mgr(0,0);
int coeff = 1;
int offset = 0;
int max_level = 4;

// #define DESCENDING_NONINTERLEAVED 1
// #define DESCENDING_INTERLEAVED 1
// #define ASCENDING_NONINTERLEAVED 0
#define ASCENDING_INTERLEAVED 1
// #define AMANO_8BIT_OPTIMAL 1
// #define NEW_ORDERING 1
// #define NEW_ORDERING_REV 1

#ifdef DESCENDING_NONINTERLEAVED
BDD a7 = mgr.bddVar(coeff*0 + offset);
BDD a6 = mgr.bddVar(coeff*1 + offset);
BDD a5 = mgr.bddVar(coeff*2 + offset);
BDD a4 = mgr.bddVar(coeff*3 + offset);
BDD a3 = mgr.bddVar(coeff*4 + offset);
BDD a2 = mgr.bddVar(coeff*5 + offset);
BDD a1 = mgr.bddVar(coeff*6 + offset);
BDD a0 = mgr.bddVar(coeff*7 + offset);
BDD b7 = mgr.bddVar(coeff*8 + offset);
BDD b6 = mgr.bddVar(coeff*9 + offset);
BDD b5 = mgr.bddVar(coeff*10 + offset);
BDD b4 = mgr.bddVar(coeff*11 + offset);
BDD b3 = mgr.bddVar(coeff*12 + offset);
BDD b2 = mgr.bddVar(coeff*13 + offset);
BDD b1 = mgr.bddVar(coeff*14 + offset);
BDD b0 = mgr.bddVar(coeff*15 + offset);
#endif

#ifdef DESCENDING_INTERLEAVED
BDD a7 = mgr.bddVar(coeff*0 + offset);
BDD b7 = mgr.bddVar(coeff*1 + offset);
BDD a6 = mgr.bddVar(coeff*2 + offset);
BDD b6 = mgr.bddVar(coeff*3 + offset);
BDD a5 = mgr.bddVar(coeff*4 + offset);
BDD b5 = mgr.bddVar(coeff*5 + offset);
BDD a4 = mgr.bddVar(coeff*6 + offset);
BDD b4 = mgr.bddVar(coeff*7 + offset);
BDD a3 = mgr.bddVar(coeff*8 + offset);
BDD b3 = mgr.bddVar(coeff*9 + offset);
BDD a2 = mgr.bddVar(coeff*10 + offset);
BDD b2 = mgr.bddVar(coeff*11 + offset);
BDD a1 = mgr.bddVar(coeff*12 + offset);
BDD b1 = mgr.bddVar(coeff*13 + offset);
BDD a0 = mgr.bddVar(coeff*14 + offset);
BDD b0 = mgr.bddVar(coeff*15 + offset);
#endif

#ifdef ASCENDING_NONINTERLEAVED
BDD a0 = mgr.bddVar(coeff*0 + offset);
BDD a1 = mgr.bddVar(coeff*1 + offset);
BDD a2 = mgr.bddVar(coeff*2 + offset);
BDD a3 = mgr.bddVar(coeff*3 + offset);
BDD a4 = mgr.bddVar(coeff*4 + offset);
BDD a5 = mgr.bddVar(coeff*5 + offset);
BDD a6 = mgr.bddVar(coeff*6 + offset);
BDD a7 = mgr.bddVar(coeff*7 + offset);
BDD b0 = mgr.bddVar(coeff*8 + offset);
BDD b1 = mgr.bddVar(coeff*9 + offset);
BDD b2 = mgr.bddVar(coeff*10 + offset);
BDD b3 = mgr.bddVar(coeff*11 + offset);
BDD b4 = mgr.bddVar(coeff*12 + offset);
BDD b5 = mgr.bddVar(coeff*13 + offset);
BDD b6 = mgr.bddVar(coeff*14 + offset);
BDD b7 = mgr.bddVar(coeff*15 + offset);
#endif

#ifdef ASCENDING_INTERLEAVED
BDD a0 = mgr.bddVar(coeff*0 + offset);
BDD b0 = mgr.bddVar(coeff*1 + offset);
BDD a1 = mgr.bddVar(coeff*2 + offset);
BDD b1 = mgr.bddVar(coeff*3 + offset);
BDD a2 = mgr.bddVar(coeff*4 + offset);
BDD b2 = mgr.bddVar(coeff*5 + offset);
BDD a3 = mgr.bddVar(coeff*6 + offset);
BDD b3 = mgr.bddVar(coeff*7 + offset);
BDD a4 = mgr.bddVar(coeff*8 + offset);
BDD b4 = mgr.bddVar(coeff*9 + offset);
BDD a5 = mgr.bddVar(coeff*10 + offset);
BDD b5 = mgr.bddVar(coeff*11 + offset);
BDD a6 = mgr.bddVar(coeff*12 + offset);
BDD b6 = mgr.bddVar(coeff*13 + offset);
BDD a7 = mgr.bddVar(coeff*14 + offset);
BDD b7 = mgr.bddVar(coeff*15 + offset);
#endif

#ifdef DESCENDING_ASCENDING
BDD a7 = mgr.bddVar(coeff*0 + offset);
BDD a6 = mgr.bddVar(coeff*1 + offset);
BDD a5 = mgr.bddVar(coeff*2 + offset);
BDD a4 = mgr.bddVar(coeff*3 + offset);
BDD a3 = mgr.bddVar(coeff*4 + offset);
BDD a2 = mgr.bddVar(coeff*5 + offset);
BDD a1 = mgr.bddVar(coeff*6 + offset);
BDD a0 = mgr.bddVar(coeff*7 + offset);
BDD b0 = mgr.bddVar(coeff*8 + offset);
BDD b1 = mgr.bddVar(coeff*9 + offset);
BDD b2 = mgr.bddVar(coeff*10 + offset);
BDD b3 = mgr.bddVar(coeff*11 + offset);
BDD b4 = mgr.bddVar(coeff*12 + offset);
BDD b5 = mgr.bddVar(coeff*13 + offset);
BDD b6 = mgr.bddVar(coeff*14 + offset);
BDD b7 = mgr.bddVar(coeff*15 + offset);
#endif

#ifdef ASCENDING_DESCENDING
BDD a0 = mgr.bddVar(coeff*0 + offset);
BDD a1 = mgr.bddVar(coeff*1 + offset);
BDD a2 = mgr.bddVar(coeff*2 + offset);
BDD a3 = mgr.bddVar(coeff*3 + offset);
BDD a4 = mgr.bddVar(coeff*4 + offset);
BDD a5 = mgr.bddVar(coeff*5 + offset);
BDD a6 = mgr.bddVar(coeff*6 + offset);
BDD a7 = mgr.bddVar(coeff*7 + offset);
BDD b7 = mgr.bddVar(coeff*8 + offset);
BDD b6 = mgr.bddVar(coeff*9 + offset);
BDD b5 = mgr.bddVar(coeff*10 + offset);
BDD b4 = mgr.bddVar(coeff*11 + offset);
BDD b3 = mgr.bddVar(coeff*12 + offset);
BDD b2 = mgr.bddVar(coeff*13 + offset);
BDD b1 = mgr.bddVar(coeff*14 + offset);
BDD b0 = mgr.bddVar(coeff*15 + offset);
#endif

// Amano & Maruoka, DAC 2007: the variable order
// that gives the smallest 8-bit multiplication BDD for MUL_7,8 is
// (a1, a2, a3, a4, b3, b4, b2, a5, b5, b1, a6, b6, a0, b7, a7, b0)
#ifdef AMANO_8BIT_OPTIMAL
BDD a1 = mgr.bddVar(coeff*0 + offset);
BDD a2 = mgr.bddVar(coeff*1 + offset);
BDD a3 = mgr.bddVar(coeff*2 + offset);
BDD a4 = mgr.bddVar(coeff*3 + offset);
BDD b3 = mgr.bddVar(coeff*4 + offset);
BDD b4 = mgr.bddVar(coeff*5 + offset);
BDD b2 = mgr.bddVar(coeff*6 + offset);
BDD a5 = mgr.bddVar(coeff*7 + offset);
BDD b5 = mgr.bddVar(coeff*8 + offset);
BDD b1 = mgr.bddVar(coeff*9 + offset);
BDD a6 = mgr.bddVar(coeff*10 + offset);
BDD b6 = mgr.bddVar(coeff*11 + offset);
BDD a0 = mgr.bddVar(coeff*12 + offset);
BDD b7 = mgr.bddVar(coeff*13 + offset);
BDD a7 = mgr.bddVar(coeff*14 + offset);
BDD b0 = mgr.bddVar(coeff*15 + offset);
#endif

#ifdef NEW_ORDERING
BDD a0 = mgr.bddVar(0);
BDD b0 = mgr.bddVar(32768);
BDD a1 = mgr.bddVar(49152);
BDD b1 = mgr.bddVar(57344);
BDD a2 = mgr.bddVar(61440);
BDD b2 = mgr.bddVar(63488);
BDD a3 = mgr.bddVar(64512);
BDD b3 = mgr.bddVar(65024);
BDD a4 = mgr.bddVar(65280);
BDD b4 = mgr.bddVar(65408);
BDD a5 = mgr.bddVar(65472);
BDD b5 = mgr.bddVar(65504);
BDD a6 = mgr.bddVar(65520);
BDD b6 = mgr.bddVar(65528);
BDD a7 = mgr.bddVar(65532);
BDD b7 = mgr.bddVar(65534);
#endif

// #ifdef NEW_ORDERING_REV
// BDD a0 = mgr.bddVar(0);
// BDD b0 = mgr.bddVar(2);
// BDD a1 = mgr.bddVar(4);
// BDD b1 = mgr.bddVar(8);
// BDD a2 = mgr.bddVar(16);
// BDD b2 = mgr.bddVar(32);
// BDD a3 = mgr.bddVar(64);
// BDD b3 = mgr.bddVar(128);
// BDD a4 = mgr.bddVar(256);
// BDD b4 = mgr.bddVar(512);
// BDD a5 = mgr.bddVar(1024);
// BDD b5 = mgr.bddVar(2048);
// BDD a6 = mgr.bddVar(4096);
// BDD b6 = mgr.bddVar(8192);
// BDD a7 = mgr.bddVar(16384);
// BDD b7 = mgr.bddVar(32768);
// #endif

#ifdef NEW_ORDERING_REV
BDD a0 = mgr.bddVar(0);
BDD b0 = mgr.bddVar(1);
BDD a1 = mgr.bddVar(2);
BDD b1 = mgr.bddVar(4);
BDD a2 = mgr.bddVar(8);
BDD b2 = mgr.bddVar(16);
BDD a3 = mgr.bddVar(32);
BDD b3 = mgr.bddVar(64);
BDD a4 = mgr.bddVar(128);
BDD b4 = mgr.bddVar(256);
BDD a5 = mgr.bddVar(512);
BDD b5 = mgr.bddVar(1024);
BDD a6 = mgr.bddVar(2048);
BDD b6 = mgr.bddVar(4096);
BDD a7 = mgr.bddVar(8192);
BDD b7 = mgr.bddVar(16384);
#endif

std::cout << "	pp0_0	 " << std::endl;
BDD	pp0_0	 = MkAnd(a0, b0);
BDD	pp1_0	 = MkAnd(a0, b1);
BDD	pp2_0	 = MkAnd(a0, b2);
BDD	pp3_0	 = MkAnd(a0, b3);
BDD	pp4_0	 = MkAnd(a0, b4);
BDD	pp5_0	 = MkAnd(a0, b5);
BDD	pp6_0	 = MkAnd(a0, b6);
BDD	pp7_0	 = MkAnd(a0, b7);
BDD	sum0	= pp0_0;
BDD	pp0_1	 = MkAnd(a1, b0);
std::cout << "	pp1_1	 " << std::endl;
BDD	pp1_1	 = MkAnd(a1, b1);
BDD	pp2_1	 = MkAnd(a1, b2);
BDD	pp3_1	 = MkAnd(a1, b3);
BDD	pp4_1	 = MkAnd(a1, b4);
BDD	pp5_1	 = MkAnd(a1, b5);
BDD	pp6_1	 = MkAnd(a1, b6);
BDD	pp7_1	 = MkAnd(a1, b7);
BDD	gate0_1_0	=  MkNor(pp1_0, pp0_1);
BDD	gate0_1_1	=  MkNor(pp1_0, gate0_1_0);
BDD	gate0_1_2	=  MkNor(pp0_1, gate0_1_0);
std::cout << "	gate0_1_3	" << std::endl;
BDD	gate0_1_3	=  MkNor(gate0_1_1, gate0_1_2);
BDD	gate0_1_4	= MkNot(gate0_1_3); 
BDD	gate0_1_5	= MkNor(gate0_1_4, gate0_1_3); 
BDD	gate0_1_6	= MkNot(gate0_1_4); 
BDD	sum1	= MkNor(gate0_1_5,gate0_1_6); 
BDD	c0_1	= MkNor(gate0_1_0,gate0_1_4); 
BDD	gate1_1_0	=  MkNor(pp2_0, pp1_1);
BDD	gate1_1_1	=  MkNor(pp2_0, gate1_1_0);
BDD	gate1_1_2	=  MkNor(pp1_1, gate1_1_0);
BDD	gate1_1_3	=  MkNor(gate1_1_1, gate1_1_2);
std::cout << "	gate1_1_4	" << std::endl;
BDD	gate1_1_4	= MkNot(gate1_1_3); 
BDD	gate1_1_5	= MkNor(gate1_1_4, gate1_1_3); 
BDD	gate1_1_6	= MkNot(gate1_1_4); 
BDD	s1_1	= MkNor(gate1_1_5,gate1_1_6); 
BDD	c1_1	= MkNor(gate1_1_0,gate1_1_4); 
BDD	gate2_1_0	=  MkNor(pp3_0, pp2_1);
BDD	gate2_1_1	=  MkNor(pp3_0, gate2_1_0);
BDD	gate2_1_2	=  MkNor(pp2_1, gate2_1_0);
BDD	gate2_1_3	=  MkNor(gate2_1_1, gate2_1_2);
BDD	gate2_1_4	= MkNot(gate2_1_3); 
std::cout << "	gate2_1_5	" << std::endl;
BDD	gate2_1_5	= MkNor(gate2_1_4, gate2_1_3); 
BDD	gate2_1_6	= MkNot(gate2_1_4); 
BDD	s2_1	= MkNor(gate2_1_5,gate2_1_6); 
BDD	c2_1	= MkNor(gate2_1_0,gate2_1_4); 
BDD	gate3_1_0	=  MkNor(pp4_0, pp3_1);
BDD	gate3_1_1	=  MkNor(pp4_0, gate3_1_0);
BDD	gate3_1_2	=  MkNor(pp3_1, gate3_1_0);
BDD	gate3_1_3	=  MkNor(gate3_1_1, gate3_1_2);
BDD	gate3_1_4	= MkNot(gate3_1_3); 
BDD	gate3_1_5	= MkNor(gate3_1_4, gate3_1_3); 
std::cout << "	gate3_1_6	" << std::endl;
BDD	gate3_1_6	= MkNot(gate3_1_4); 
BDD	s3_1	= MkNor(gate3_1_5,gate3_1_6); 
BDD	c3_1	= MkNor(gate3_1_0,gate3_1_4); 
BDD	gate4_1_0	=  MkNor(pp5_0, pp4_1);
BDD	gate4_1_1	=  MkNor(pp5_0, gate4_1_0);
BDD	gate4_1_2	=  MkNor(pp4_1, gate4_1_0);
BDD	gate4_1_3	=  MkNor(gate4_1_1, gate4_1_2);
BDD	gate4_1_4	= MkNot(gate4_1_3); 
BDD	gate4_1_5	= MkNor(gate4_1_4, gate4_1_3); 
BDD	gate4_1_6	= MkNot(gate4_1_4); 
std::cout << "	s4_1	" << std::endl;
BDD	s4_1	= MkNor(gate4_1_5,gate4_1_6); 
BDD	c4_1	= MkNor(gate4_1_0,gate4_1_4); 
BDD	gate5_1_0	=  MkNor(pp6_0, pp5_1);
BDD	gate5_1_1	=  MkNor(pp6_0, gate5_1_0);
BDD	gate5_1_2	=  MkNor(pp5_1, gate5_1_0);
BDD	gate5_1_3	=  MkNor(gate5_1_1, gate5_1_2);
BDD	gate5_1_4	= MkNot(gate5_1_3); 
BDD	gate5_1_5	= MkNor(gate5_1_4, gate5_1_3); 
BDD	gate5_1_6	= MkNot(gate5_1_4); 
BDD	s5_1	= MkNor(gate5_1_5,gate5_1_6); 
std::cout << "	c5_1	" << std::endl;
BDD	c5_1	= MkNor(gate5_1_0,gate5_1_4); 
BDD	gate6_1_0	=  MkNor(pp7_0, pp6_1);
BDD	gate6_1_1	=  MkNor(pp7_0, gate6_1_0);
BDD	gate6_1_2	=  MkNor(pp6_1, gate6_1_0);
BDD	gate6_1_3	=  MkNor(gate6_1_1, gate6_1_2);
BDD	gate6_1_4	= MkNot(gate6_1_3); 
BDD	gate6_1_5	= MkNor(gate6_1_4, gate6_1_3); 
BDD	gate6_1_6	= MkNot(gate6_1_4); 
BDD	s6_1	= MkNor(gate6_1_5,gate6_1_6); 
BDD	c6_1	= MkNor(gate6_1_0,gate6_1_4); 
std::cout << "	s7_1	" << std::endl;
BDD	s7_1	= pp7_1;
BDD	pp0_2	 = MkAnd(a2, b0);
BDD	pp1_2	 = MkAnd(a2, b1);
BDD	pp2_2	 = MkAnd(a2, b2);
BDD	pp3_2	 = MkAnd(a2, b3);
BDD	pp4_2	 = MkAnd(a2, b4);
BDD	pp5_2	 = MkAnd(a2, b5);
BDD	pp6_2	 = MkAnd(a2, b6);
BDD	pp7_2	 = MkAnd(a2, b7);
BDD	gate0_2_0	=  MkNor(s1_1, c0_1);
std::cout << "	gate0_2_1	" << std::endl;
BDD	gate0_2_1	=  MkNor(s1_1, gate0_2_0);
BDD	gate0_2_2	=  MkNor(c0_1, gate0_2_0);
BDD	gate0_2_3	=  MkNor(gate0_2_1, gate0_2_2);
BDD	gate0_2_4	=  MkNor(gate0_2_3, pp0_2);
BDD	gate0_2_5	=  MkNor(gate0_2_3, gate0_2_4);
BDD	gate0_2_6	=  MkNor(pp0_2, gate0_2_4);
BDD	sum2	=  MkNor(gate0_2_5, gate0_2_6);
BDD	c0_2	= MkNor(gate0_2_0,gate0_2_4); 
BDD	gate1_2_0	=  MkNor(s2_1, c1_1);
BDD	gate1_2_1	=  MkNor(s2_1, gate1_2_0);
std::cout << "	gate1_2_2	" << std::endl;
BDD	gate1_2_2	=  MkNor(c1_1, gate1_2_0);
BDD	gate1_2_3	=  MkNor(gate1_2_1, gate1_2_2);
BDD	gate1_2_4	=  MkNor(gate1_2_3, pp1_2);
BDD	gate1_2_5	=  MkNor(gate1_2_3, gate1_2_4);
BDD	gate1_2_6	=  MkNor(pp1_2, gate1_2_4);
BDD	s1_2	=  MkNor(gate1_2_5, gate1_2_6);
BDD	c1_2	= MkNor(gate1_2_0,gate1_2_4); 
BDD	gate2_2_0	=  MkNor(s3_1, c2_1);
BDD	gate2_2_1	=  MkNor(s3_1, gate2_2_0);
BDD	gate2_2_2	=  MkNor(c2_1, gate2_2_0);
std::cout << "	gate2_2_3	" << std::endl;
BDD	gate2_2_3	=  MkNor(gate2_2_1, gate2_2_2);
BDD	gate2_2_4	=  MkNor(gate2_2_3, pp2_2);
BDD	gate2_2_5	=  MkNor(gate2_2_3, gate2_2_4);
BDD	gate2_2_6	=  MkNor(pp2_2, gate2_2_4);
BDD	s2_2	=  MkNor(gate2_2_5, gate2_2_6);
BDD	c2_2	= MkNor(gate2_2_0,gate2_2_4); 
BDD	gate3_2_0	=  MkNor(s4_1, c3_1);
BDD	gate3_2_1	=  MkNor(s4_1, gate3_2_0);
BDD	gate3_2_2	=  MkNor(c3_1, gate3_2_0);
BDD	gate3_2_3	=  MkNor(gate3_2_1, gate3_2_2);
std::cout << "	gate3_2_4	" << std::endl;
BDD	gate3_2_4	=  MkNor(gate3_2_3, pp3_2);
BDD	gate3_2_5	=  MkNor(gate3_2_3, gate3_2_4);
BDD	gate3_2_6	=  MkNor(pp3_2, gate3_2_4);
BDD	s3_2	=  MkNor(gate3_2_5, gate3_2_6);
BDD	c3_2	= MkNor(gate3_2_0,gate3_2_4); 
BDD	gate4_2_0	=  MkNor(s5_1, c4_1);
BDD	gate4_2_1	=  MkNor(s5_1, gate4_2_0);
BDD	gate4_2_2	=  MkNor(c4_1, gate4_2_0);
BDD	gate4_2_3	=  MkNor(gate4_2_1, gate4_2_2);
BDD	gate4_2_4	=  MkNor(gate4_2_3, pp4_2);
std::cout << "	gate4_2_5	" << std::endl;
BDD	gate4_2_5	=  MkNor(gate4_2_3, gate4_2_4);
BDD	gate4_2_6	=  MkNor(pp4_2, gate4_2_4);
BDD	s4_2	=  MkNor(gate4_2_5, gate4_2_6);
BDD	c4_2	= MkNor(gate4_2_0,gate4_2_4); 
BDD	gate5_2_0	=  MkNor(s6_1, c5_1);
BDD	gate5_2_1	=  MkNor(s6_1, gate5_2_0);
BDD	gate5_2_2	=  MkNor(c5_1, gate5_2_0);
BDD	gate5_2_3	=  MkNor(gate5_2_1, gate5_2_2);
BDD	gate5_2_4	=  MkNor(gate5_2_3, pp5_2);
BDD	gate5_2_5	=  MkNor(gate5_2_3, gate5_2_4);
std::cout << "	gate5_2_6	" << std::endl;
BDD	gate5_2_6	=  MkNor(pp5_2, gate5_2_4);
BDD	s5_2	=  MkNor(gate5_2_5, gate5_2_6);
BDD	c5_2	= MkNor(gate5_2_0,gate5_2_4); 
BDD	gate6_2_0	=  MkNor(s7_1, c6_1);
BDD	gate6_2_1	=  MkNor(s7_1, gate6_2_0);
BDD	gate6_2_2	=  MkNor(c6_1, gate6_2_0);
BDD	gate6_2_3	=  MkNor(gate6_2_1, gate6_2_2);
BDD	gate6_2_4	=  MkNor(gate6_2_3, pp6_2);
BDD	gate6_2_5	=  MkNor(gate6_2_3, gate6_2_4);
BDD	gate6_2_6	=  MkNor(pp6_2, gate6_2_4);
std::cout << "	s6_2	" << std::endl;
BDD	s6_2	=  MkNor(gate6_2_5, gate6_2_6);
BDD	c6_2	= MkNor(gate6_2_0,gate6_2_4); 
BDD	s7_2	= pp7_2;
BDD	pp0_3	 = MkAnd(a3, b0);
BDD	pp1_3	 = MkAnd(a3, b1);
BDD	pp2_3	 = MkAnd(a3, b2);
BDD	pp3_3	 = MkAnd(a3, b3);
BDD	pp4_3	 = MkAnd(a3, b4);
BDD	pp5_3	 = MkAnd(a3, b5);
BDD	pp6_3	 = MkAnd(a3, b6);
std::cout << "	pp7_3	 " << std::endl;
BDD	pp7_3	 = MkAnd(a3, b7);
BDD	gate0_3_0	=  MkNor(s1_2, c0_2);
BDD	gate0_3_1	=  MkNor(s1_2, gate0_3_0);
BDD	gate0_3_2	=  MkNor(c0_2, gate0_3_0);
BDD	gate0_3_3	=  MkNor(gate0_3_1, gate0_3_2);
BDD	gate0_3_4	=  MkNor(gate0_3_3, pp0_3);
BDD	gate0_3_5	=  MkNor(gate0_3_3, gate0_3_4);
BDD	gate0_3_6	=  MkNor(pp0_3, gate0_3_4);
BDD	sum3	=  MkNor(gate0_3_5, gate0_3_6);
BDD	c0_3	= MkNor(gate0_3_0,gate0_3_4); 
std::cout << "	gate1_3_0	" << std::endl;
BDD	gate1_3_0	=  MkNor(s2_2, c1_2);
BDD	gate1_3_1	=  MkNor(s2_2, gate1_3_0);
BDD	gate1_3_2	=  MkNor(c1_2, gate1_3_0);
BDD	gate1_3_3	=  MkNor(gate1_3_1, gate1_3_2);
BDD	gate1_3_4	=  MkNor(gate1_3_3, pp1_3);
BDD	gate1_3_5	=  MkNor(gate1_3_3, gate1_3_4);
BDD	gate1_3_6	=  MkNor(pp1_3, gate1_3_4);
BDD	s1_3	=  MkNor(gate1_3_5, gate1_3_6);
BDD	c1_3	= MkNor(gate1_3_0,gate1_3_4); 
BDD	gate2_3_0	=  MkNor(s3_2, c2_2);
std::cout << "	gate2_3_1	" << std::endl;
BDD	gate2_3_1	=  MkNor(s3_2, gate2_3_0);
BDD	gate2_3_2	=  MkNor(c2_2, gate2_3_0);
BDD	gate2_3_3	=  MkNor(gate2_3_1, gate2_3_2);
BDD	gate2_3_4	=  MkNor(gate2_3_3, pp2_3);
BDD	gate2_3_5	=  MkNor(gate2_3_3, gate2_3_4);
BDD	gate2_3_6	=  MkNor(pp2_3, gate2_3_4);
BDD	s2_3	=  MkNor(gate2_3_5, gate2_3_6);
BDD	c2_3	= MkNor(gate2_3_0,gate2_3_4); 
BDD	gate3_3_0	=  MkNor(s4_2, c3_2);
BDD	gate3_3_1	=  MkNor(s4_2, gate3_3_0);
std::cout << "	gate3_3_2	" << std::endl;
BDD	gate3_3_2	=  MkNor(c3_2, gate3_3_0);
BDD	gate3_3_3	=  MkNor(gate3_3_1, gate3_3_2);
BDD	gate3_3_4	=  MkNor(gate3_3_3, pp3_3);
BDD	gate3_3_5	=  MkNor(gate3_3_3, gate3_3_4);
BDD	gate3_3_6	=  MkNor(pp3_3, gate3_3_4);
BDD	s3_3	=  MkNor(gate3_3_5, gate3_3_6);
BDD	c3_3	= MkNor(gate3_3_0,gate3_3_4); 
BDD	gate4_3_0	=  MkNor(s5_2, c4_2);
BDD	gate4_3_1	=  MkNor(s5_2, gate4_3_0);
BDD	gate4_3_2	=  MkNor(c4_2, gate4_3_0);
std::cout << "	gate4_3_3	" << std::endl;
BDD	gate4_3_3	=  MkNor(gate4_3_1, gate4_3_2);
BDD	gate4_3_4	=  MkNor(gate4_3_3, pp4_3);
BDD	gate4_3_5	=  MkNor(gate4_3_3, gate4_3_4);
BDD	gate4_3_6	=  MkNor(pp4_3, gate4_3_4);
BDD	s4_3	=  MkNor(gate4_3_5, gate4_3_6);
BDD	c4_3	= MkNor(gate4_3_0,gate4_3_4); 
BDD	gate5_3_0	=  MkNor(s6_2, c5_2);
BDD	gate5_3_1	=  MkNor(s6_2, gate5_3_0);
BDD	gate5_3_2	=  MkNor(c5_2, gate5_3_0);
BDD	gate5_3_3	=  MkNor(gate5_3_1, gate5_3_2);
std::cout << "	gate5_3_4	" << std::endl;
BDD	gate5_3_4	=  MkNor(gate5_3_3, pp5_3);
BDD	gate5_3_5	=  MkNor(gate5_3_3, gate5_3_4);
BDD	gate5_3_6	=  MkNor(pp5_3, gate5_3_4);
BDD	s5_3	=  MkNor(gate5_3_5, gate5_3_6);
BDD	c5_3	= MkNor(gate5_3_0,gate5_3_4); 
BDD	gate6_3_0	=  MkNor(s7_2, c6_2);
BDD	gate6_3_1	=  MkNor(s7_2, gate6_3_0);
BDD	gate6_3_2	=  MkNor(c6_2, gate6_3_0);
BDD	gate6_3_3	=  MkNor(gate6_3_1, gate6_3_2);
BDD	gate6_3_4	=  MkNor(gate6_3_3, pp6_3);
std::cout << "	gate6_3_5	" << std::endl;
BDD	gate6_3_5	=  MkNor(gate6_3_3, gate6_3_4);
BDD	gate6_3_6	=  MkNor(pp6_3, gate6_3_4);
BDD	s6_3	=  MkNor(gate6_3_5, gate6_3_6);
BDD	c6_3	= MkNor(gate6_3_0,gate6_3_4); 
BDD	s7_3	= pp7_3;
BDD	pp0_4	 = MkAnd(a4, b0);
BDD	pp1_4	 = MkAnd(a4, b1);
BDD	pp2_4	 = MkAnd(a4, b2);
BDD	pp3_4	 = MkAnd(a4, b3);
BDD	pp4_4	 = MkAnd(a4, b4);
std::cout << "	pp5_4	 " << std::endl;
BDD	pp5_4	 = MkAnd(a4, b5);
BDD	pp6_4	 = MkAnd(a4, b6);
BDD	pp7_4	 = MkAnd(a4, b7);
BDD	gate0_4_0	=  MkNor(s1_3, c0_3);
BDD	gate0_4_1	=  MkNor(s1_3, gate0_4_0);
BDD	gate0_4_2	=  MkNor(c0_3, gate0_4_0);
BDD	gate0_4_3	=  MkNor(gate0_4_1, gate0_4_2);
BDD	gate0_4_4	=  MkNor(gate0_4_3, pp0_4);
BDD	gate0_4_5	=  MkNor(gate0_4_3, gate0_4_4);
BDD	gate0_4_6	=  MkNor(pp0_4, gate0_4_4);
std::cout << "	sum4	" << std::endl;
BDD	sum4	=  MkNor(gate0_4_5, gate0_4_6);
BDD	c0_4	= MkNor(gate0_4_0,gate0_4_4); 
BDD	gate1_4_0	=  MkNor(s2_3, c1_3);
BDD	gate1_4_1	=  MkNor(s2_3, gate1_4_0);
BDD	gate1_4_2	=  MkNor(c1_3, gate1_4_0);
BDD	gate1_4_3	=  MkNor(gate1_4_1, gate1_4_2);
BDD	gate1_4_4	=  MkNor(gate1_4_3, pp1_4);
BDD	gate1_4_5	=  MkNor(gate1_4_3, gate1_4_4);
BDD	gate1_4_6	=  MkNor(pp1_4, gate1_4_4);
BDD	s1_4	=  MkNor(gate1_4_5, gate1_4_6);
std::cout << "	c1_4	" << std::endl;
BDD	c1_4	= MkNor(gate1_4_0,gate1_4_4); 
BDD	gate2_4_0	=  MkNor(s3_3, c2_3);
BDD	gate2_4_1	=  MkNor(s3_3, gate2_4_0);
BDD	gate2_4_2	=  MkNor(c2_3, gate2_4_0);
BDD	gate2_4_3	=  MkNor(gate2_4_1, gate2_4_2);
BDD	gate2_4_4	=  MkNor(gate2_4_3, pp2_4);
BDD	gate2_4_5	=  MkNor(gate2_4_3, gate2_4_4);
BDD	gate2_4_6	=  MkNor(pp2_4, gate2_4_4);
BDD	s2_4	=  MkNor(gate2_4_5, gate2_4_6);
BDD	c2_4	= MkNor(gate2_4_0,gate2_4_4); 
std::cout << "	gate3_4_0	" << std::endl;
BDD	gate3_4_0	=  MkNor(s4_3, c3_3);
BDD	gate3_4_1	=  MkNor(s4_3, gate3_4_0);
BDD	gate3_4_2	=  MkNor(c3_3, gate3_4_0);
BDD	gate3_4_3	=  MkNor(gate3_4_1, gate3_4_2);
BDD	gate3_4_4	=  MkNor(gate3_4_3, pp3_4);
BDD	gate3_4_5	=  MkNor(gate3_4_3, gate3_4_4);
BDD	gate3_4_6	=  MkNor(pp3_4, gate3_4_4);
BDD	s3_4	=  MkNor(gate3_4_5, gate3_4_6);
BDD	c3_4	= MkNor(gate3_4_0,gate3_4_4); 
BDD	gate4_4_0	=  MkNor(s5_3, c4_3);
std::cout << "	gate4_4_1	" << std::endl;
BDD	gate4_4_1	=  MkNor(s5_3, gate4_4_0);
BDD	gate4_4_2	=  MkNor(c4_3, gate4_4_0);
BDD	gate4_4_3	=  MkNor(gate4_4_1, gate4_4_2);
BDD	gate4_4_4	=  MkNor(gate4_4_3, pp4_4);
BDD	gate4_4_5	=  MkNor(gate4_4_3, gate4_4_4);
BDD	gate4_4_6	=  MkNor(pp4_4, gate4_4_4);
BDD	s4_4	=  MkNor(gate4_4_5, gate4_4_6);
BDD	c4_4	= MkNor(gate4_4_0,gate4_4_4); 
BDD	gate5_4_0	=  MkNor(s6_3, c5_3);
BDD	gate5_4_1	=  MkNor(s6_3, gate5_4_0);
std::cout << "	gate5_4_2	" << std::endl;
BDD	gate5_4_2	=  MkNor(c5_3, gate5_4_0);
BDD	gate5_4_3	=  MkNor(gate5_4_1, gate5_4_2);
BDD	gate5_4_4	=  MkNor(gate5_4_3, pp5_4);
BDD	gate5_4_5	=  MkNor(gate5_4_3, gate5_4_4);
BDD	gate5_4_6	=  MkNor(pp5_4, gate5_4_4);
BDD	s5_4	=  MkNor(gate5_4_5, gate5_4_6);
BDD	c5_4	= MkNor(gate5_4_0,gate5_4_4); 
BDD	gate6_4_0	=  MkNor(s7_3, c6_3);
BDD	gate6_4_1	=  MkNor(s7_3, gate6_4_0);
BDD	gate6_4_2	=  MkNor(c6_3, gate6_4_0);
std::cout << "	gate6_4_3	" << std::endl;
BDD	gate6_4_3	=  MkNor(gate6_4_1, gate6_4_2);
BDD	gate6_4_4	=  MkNor(gate6_4_3, pp6_4);
BDD	gate6_4_5	=  MkNor(gate6_4_3, gate6_4_4);
BDD	gate6_4_6	=  MkNor(pp6_4, gate6_4_4);
BDD	s6_4	=  MkNor(gate6_4_5, gate6_4_6);
BDD	c6_4	= MkNor(gate6_4_0,gate6_4_4); 
BDD	s7_4	= pp7_4;
BDD	pp0_5	 = MkAnd(a5, b0);
BDD	pp1_5	 = MkAnd(a5, b1);
BDD	pp2_5	 = MkAnd(a5, b2);
std::cout << "	pp3_5	 " << std::endl;
BDD	pp3_5	 = MkAnd(a5, b3);
BDD	pp4_5	 = MkAnd(a5, b4);
BDD	pp5_5	 = MkAnd(a5, b5);
BDD	pp6_5	 = MkAnd(a5, b6);
BDD	pp7_5	 = MkAnd(a5, b7);
BDD	gate0_5_0	=  MkNor(s1_4, c0_4);
BDD	gate0_5_1	=  MkNor(s1_4, gate0_5_0);
BDD	gate0_5_2	=  MkNor(c0_4, gate0_5_0);
BDD	gate0_5_3	=  MkNor(gate0_5_1, gate0_5_2);
BDD	gate0_5_4	=  MkNor(gate0_5_3, pp0_5);
std::cout << "	gate0_5_5	" << std::endl;
BDD	gate0_5_5	=  MkNor(gate0_5_3, gate0_5_4);
BDD	gate0_5_6	=  MkNor(pp0_5, gate0_5_4);
BDD	sum5	=  MkNor(gate0_5_5, gate0_5_6);
BDD	c0_5	= MkNor(gate0_5_0,gate0_5_4); 
BDD	gate1_5_0	=  MkNor(s2_4, c1_4);
BDD	gate1_5_1	=  MkNor(s2_4, gate1_5_0);
BDD	gate1_5_2	=  MkNor(c1_4, gate1_5_0);
BDD	gate1_5_3	=  MkNor(gate1_5_1, gate1_5_2);
BDD	gate1_5_4	=  MkNor(gate1_5_3, pp1_5);
BDD	gate1_5_5	=  MkNor(gate1_5_3, gate1_5_4);
std::cout << "	gate1_5_6	" << std::endl;
BDD	gate1_5_6	=  MkNor(pp1_5, gate1_5_4);
BDD	s1_5	=  MkNor(gate1_5_5, gate1_5_6);
BDD	c1_5	= MkNor(gate1_5_0,gate1_5_4); 
BDD	gate2_5_0	=  MkNor(s3_4, c2_4);
BDD	gate2_5_1	=  MkNor(s3_4, gate2_5_0);
BDD	gate2_5_2	=  MkNor(c2_4, gate2_5_0);
BDD	gate2_5_3	=  MkNor(gate2_5_1, gate2_5_2);
BDD	gate2_5_4	=  MkNor(gate2_5_3, pp2_5);
BDD	gate2_5_5	=  MkNor(gate2_5_3, gate2_5_4);
BDD	gate2_5_6	=  MkNor(pp2_5, gate2_5_4);
std::cout << "	s2_5	" << std::endl;
BDD	s2_5	=  MkNor(gate2_5_5, gate2_5_6);
BDD	c2_5	= MkNor(gate2_5_0,gate2_5_4); 
BDD	gate3_5_0	=  MkNor(s4_4, c3_4);
BDD	gate3_5_1	=  MkNor(s4_4, gate3_5_0);
BDD	gate3_5_2	=  MkNor(c3_4, gate3_5_0);
BDD	gate3_5_3	=  MkNor(gate3_5_1, gate3_5_2);
BDD	gate3_5_4	=  MkNor(gate3_5_3, pp3_5);
BDD	gate3_5_5	=  MkNor(gate3_5_3, gate3_5_4);
BDD	gate3_5_6	=  MkNor(pp3_5, gate3_5_4);
BDD	s3_5	=  MkNor(gate3_5_5, gate3_5_6);
std::cout << "	c3_5	" << std::endl;
BDD	c3_5	= MkNor(gate3_5_0,gate3_5_4); 
BDD	gate4_5_0	=  MkNor(s5_4, c4_4);
BDD	gate4_5_1	=  MkNor(s5_4, gate4_5_0);
BDD	gate4_5_2	=  MkNor(c4_4, gate4_5_0);
BDD	gate4_5_3	=  MkNor(gate4_5_1, gate4_5_2);
BDD	gate4_5_4	=  MkNor(gate4_5_3, pp4_5);
BDD	gate4_5_5	=  MkNor(gate4_5_3, gate4_5_4);
BDD	gate4_5_6	=  MkNor(pp4_5, gate4_5_4);
BDD	s4_5	=  MkNor(gate4_5_5, gate4_5_6);
BDD	c4_5	= MkNor(gate4_5_0,gate4_5_4); 
std::cout << "	gate5_5_0	" << std::endl;
BDD	gate5_5_0	=  MkNor(s6_4, c5_4);
BDD	gate5_5_1	=  MkNor(s6_4, gate5_5_0);
BDD	gate5_5_2	=  MkNor(c5_4, gate5_5_0);
BDD	gate5_5_3	=  MkNor(gate5_5_1, gate5_5_2);
BDD	gate5_5_4	=  MkNor(gate5_5_3, pp5_5);
BDD	gate5_5_5	=  MkNor(gate5_5_3, gate5_5_4);
BDD	gate5_5_6	=  MkNor(pp5_5, gate5_5_4);
BDD	s5_5	=  MkNor(gate5_5_5, gate5_5_6);
BDD	c5_5	= MkNor(gate5_5_0,gate5_5_4); 
BDD	gate6_5_0	=  MkNor(s7_4, c6_4);
std::cout << "	gate6_5_1	" << std::endl;
BDD	gate6_5_1	=  MkNor(s7_4, gate6_5_0);
BDD	gate6_5_2	=  MkNor(c6_4, gate6_5_0);
BDD	gate6_5_3	=  MkNor(gate6_5_1, gate6_5_2);
BDD	gate6_5_4	=  MkNor(gate6_5_3, pp6_5);
BDD	gate6_5_5	=  MkNor(gate6_5_3, gate6_5_4);
BDD	gate6_5_6	=  MkNor(pp6_5, gate6_5_4);
BDD	s6_5	=  MkNor(gate6_5_5, gate6_5_6);
BDD	c6_5	= MkNor(gate6_5_0,gate6_5_4); 
BDD	s7_5	= pp7_5;
BDD	pp0_6	 = MkAnd(a6, b0);
std::cout << "	pp1_6	 " << std::endl;
BDD	pp1_6	 = MkAnd(a6, b1);
BDD	pp2_6	 = MkAnd(a6, b2);
BDD	pp3_6	 = MkAnd(a6, b3);
BDD	pp4_6	 = MkAnd(a6, b4);
BDD	pp5_6	 = MkAnd(a6, b5);
BDD	pp6_6	 = MkAnd(a6, b6);
BDD	pp7_6	 = MkAnd(a6, b7);
BDD	gate0_6_0	=  MkNor(s1_5, c0_5);
BDD	gate0_6_1	=  MkNor(s1_5, gate0_6_0);
BDD	gate0_6_2	=  MkNor(c0_5, gate0_6_0);
std::cout << "	gate0_6_3	" << std::endl;
BDD	gate0_6_3	=  MkNor(gate0_6_1, gate0_6_2);
BDD	gate0_6_4	=  MkNor(gate0_6_3, pp0_6);
BDD	gate0_6_5	=  MkNor(gate0_6_3, gate0_6_4);
BDD	gate0_6_6	=  MkNor(pp0_6, gate0_6_4);
BDD	sum6	=  MkNor(gate0_6_5, gate0_6_6);
BDD	c0_6	= MkNor(gate0_6_0,gate0_6_4); 
BDD	gate1_6_0	=  MkNor(s2_5, c1_5);
BDD	gate1_6_1	=  MkNor(s2_5, gate1_6_0);
BDD	gate1_6_2	=  MkNor(c1_5, gate1_6_0);
BDD	gate1_6_3	=  MkNor(gate1_6_1, gate1_6_2);
std::cout << "	gate1_6_4	" << std::endl;
BDD	gate1_6_4	=  MkNor(gate1_6_3, pp1_6);
BDD	gate1_6_5	=  MkNor(gate1_6_3, gate1_6_4);
BDD	gate1_6_6	=  MkNor(pp1_6, gate1_6_4);
BDD	s1_6	=  MkNor(gate1_6_5, gate1_6_6);
BDD	c1_6	= MkNor(gate1_6_0,gate1_6_4); 
BDD	gate2_6_0	=  MkNor(s3_5, c2_5);
BDD	gate2_6_1	=  MkNor(s3_5, gate2_6_0);
BDD	gate2_6_2	=  MkNor(c2_5, gate2_6_0);
BDD	gate2_6_3	=  MkNor(gate2_6_1, gate2_6_2);
BDD	gate2_6_4	=  MkNor(gate2_6_3, pp2_6);
std::cout << "	gate2_6_5	" << std::endl;
BDD	gate2_6_5	=  MkNor(gate2_6_3, gate2_6_4);
BDD	gate2_6_6	=  MkNor(pp2_6, gate2_6_4);
BDD	s2_6	=  MkNor(gate2_6_5, gate2_6_6);
BDD	c2_6	= MkNor(gate2_6_0,gate2_6_4); 
BDD	gate3_6_0	=  MkNor(s4_5, c3_5);
BDD	gate3_6_1	=  MkNor(s4_5, gate3_6_0);
BDD	gate3_6_2	=  MkNor(c3_5, gate3_6_0);
BDD	gate3_6_3	=  MkNor(gate3_6_1, gate3_6_2);
BDD	gate3_6_4	=  MkNor(gate3_6_3, pp3_6);
BDD	gate3_6_5	=  MkNor(gate3_6_3, gate3_6_4);
std::cout << "	gate3_6_6	" << std::endl;
BDD	gate3_6_6	=  MkNor(pp3_6, gate3_6_4);
BDD	s3_6	=  MkNor(gate3_6_5, gate3_6_6);
BDD	c3_6	= MkNor(gate3_6_0,gate3_6_4); 
BDD	gate4_6_0	=  MkNor(s5_5, c4_5);
BDD	gate4_6_1	=  MkNor(s5_5, gate4_6_0);
BDD	gate4_6_2	=  MkNor(c4_5, gate4_6_0);
BDD	gate4_6_3	=  MkNor(gate4_6_1, gate4_6_2);
BDD	gate4_6_4	=  MkNor(gate4_6_3, pp4_6);
BDD	gate4_6_5	=  MkNor(gate4_6_3, gate4_6_4);
BDD	gate4_6_6	=  MkNor(pp4_6, gate4_6_4);
std::cout << "	s4_6	" << std::endl;
BDD	s4_6	=  MkNor(gate4_6_5, gate4_6_6);
BDD	c4_6	= MkNor(gate4_6_0,gate4_6_4); 
BDD	gate5_6_0	=  MkNor(s6_5, c5_5);
BDD	gate5_6_1	=  MkNor(s6_5, gate5_6_0);
BDD	gate5_6_2	=  MkNor(c5_5, gate5_6_0);
BDD	gate5_6_3	=  MkNor(gate5_6_1, gate5_6_2);
BDD	gate5_6_4	=  MkNor(gate5_6_3, pp5_6);
BDD	gate5_6_5	=  MkNor(gate5_6_3, gate5_6_4);
BDD	gate5_6_6	=  MkNor(pp5_6, gate5_6_4);
BDD	s5_6	=  MkNor(gate5_6_5, gate5_6_6);
std::cout << "	c5_6	" << std::endl;
BDD	c5_6	= MkNor(gate5_6_0,gate5_6_4); 
BDD	gate6_6_0	=  MkNor(s7_5, c6_5);
BDD	gate6_6_1	=  MkNor(s7_5, gate6_6_0);
BDD	gate6_6_2	=  MkNor(c6_5, gate6_6_0);
BDD	gate6_6_3	=  MkNor(gate6_6_1, gate6_6_2);
BDD	gate6_6_4	=  MkNor(gate6_6_3, pp6_6);
BDD	gate6_6_5	=  MkNor(gate6_6_3, gate6_6_4);
BDD	gate6_6_6	=  MkNor(pp6_6, gate6_6_4);
BDD	s6_6	=  MkNor(gate6_6_5, gate6_6_6);
BDD	c6_6	= MkNor(gate6_6_0,gate6_6_4); 
std::cout << "	s7_6	" << std::endl;
BDD	s7_6	= pp7_6;
BDD	pp0_7	 = MkAnd(a7, b0);
BDD	pp1_7	 = MkAnd(a7, b1);
BDD	pp2_7	 = MkAnd(a7, b2);
BDD	pp3_7	 = MkAnd(a7, b3);
BDD	pp4_7	 = MkAnd(a7, b4);
BDD	pp5_7	 = MkAnd(a7, b5);
BDD	pp6_7	 = MkAnd(a7, b6);
BDD	pp7_7	 = MkAnd(a7, b7);
BDD	gate0_7_0	=  MkNor(s1_6, c0_6);
std::cout << "	gate0_7_1	" << std::endl;
BDD	gate0_7_1	=  MkNor(s1_6, gate0_7_0);
BDD	gate0_7_2	=  MkNor(c0_6, gate0_7_0);
BDD	gate0_7_3	=  MkNor(gate0_7_1, gate0_7_2);
BDD	gate0_7_4	=  MkNor(gate0_7_3, pp0_7);
BDD	gate0_7_5	=  MkNor(gate0_7_3, gate0_7_4);
BDD	gate0_7_6	=  MkNor(pp0_7, gate0_7_4);
BDD	sum7	=  MkNor(gate0_7_5, gate0_7_6);
BDD	c0_7	= MkNor(gate0_7_0,gate0_7_4); 
BDD	gate1_7_0	=  MkNor(s2_6, c1_6);
BDD	gate1_7_1	=  MkNor(s2_6, gate1_7_0);
std::cout << "	gate1_7_2	" << std::endl;
BDD	gate1_7_2	=  MkNor(c1_6, gate1_7_0);
BDD	gate1_7_3	=  MkNor(gate1_7_1, gate1_7_2);
BDD	gate1_7_4	=  MkNor(gate1_7_3, pp1_7);
BDD	gate1_7_5	=  MkNor(gate1_7_3, gate1_7_4);
BDD	gate1_7_6	=  MkNor(pp1_7, gate1_7_4);
BDD	s1_7	=  MkNor(gate1_7_5, gate1_7_6);
BDD	c1_7	= MkNor(gate1_7_0,gate1_7_4); 
BDD	gate2_7_0	=  MkNor(s3_6, c2_6);
BDD	gate2_7_1	=  MkNor(s3_6, gate2_7_0);
BDD	gate2_7_2	=  MkNor(c2_6, gate2_7_0);
std::cout << "	gate2_7_3	" << std::endl;
BDD	gate2_7_3	=  MkNor(gate2_7_1, gate2_7_2);
BDD	gate2_7_4	=  MkNor(gate2_7_3, pp2_7);
BDD	gate2_7_5	=  MkNor(gate2_7_3, gate2_7_4);
BDD	gate2_7_6	=  MkNor(pp2_7, gate2_7_4);
BDD	s2_7	=  MkNor(gate2_7_5, gate2_7_6);
BDD	c2_7	= MkNor(gate2_7_0,gate2_7_4); 
BDD	gate3_7_0	=  MkNor(s4_6, c3_6);
BDD	gate3_7_1	=  MkNor(s4_6, gate3_7_0);
BDD	gate3_7_2	=  MkNor(c3_6, gate3_7_0);
BDD	gate3_7_3	=  MkNor(gate3_7_1, gate3_7_2);
std::cout << "	gate3_7_4	" << std::endl;
BDD	gate3_7_4	=  MkNor(gate3_7_3, pp3_7);
BDD	gate3_7_5	=  MkNor(gate3_7_3, gate3_7_4);
BDD	gate3_7_6	=  MkNor(pp3_7, gate3_7_4);
BDD	s3_7	=  MkNor(gate3_7_5, gate3_7_6);
BDD	c3_7	= MkNor(gate3_7_0,gate3_7_4); 
BDD	gate4_7_0	=  MkNor(s5_6, c4_6);
BDD	gate4_7_1	=  MkNor(s5_6, gate4_7_0);
BDD	gate4_7_2	=  MkNor(c4_6, gate4_7_0);
BDD	gate4_7_3	=  MkNor(gate4_7_1, gate4_7_2);
BDD	gate4_7_4	=  MkNor(gate4_7_3, pp4_7);
std::cout << "	gate4_7_5	" << std::endl;
BDD	gate4_7_5	=  MkNor(gate4_7_3, gate4_7_4);
BDD	gate4_7_6	=  MkNor(pp4_7, gate4_7_4);
BDD	s4_7	=  MkNor(gate4_7_5, gate4_7_6);
BDD	c4_7	= MkNor(gate4_7_0,gate4_7_4); 
BDD	gate5_7_0	=  MkNor(s6_6, c5_6);
BDD	gate5_7_1	=  MkNor(s6_6, gate5_7_0);
BDD	gate5_7_2	=  MkNor(c5_6, gate5_7_0);
BDD	gate5_7_3	=  MkNor(gate5_7_1, gate5_7_2);
BDD	gate5_7_4	=  MkNor(gate5_7_3, pp5_7);
BDD	gate5_7_5	=  MkNor(gate5_7_3, gate5_7_4);
std::cout << "	gate5_7_6	" << std::endl;
BDD	gate5_7_6	=  MkNor(pp5_7, gate5_7_4);
BDD	s5_7	=  MkNor(gate5_7_5, gate5_7_6);
BDD	c5_7	= MkNor(gate5_7_0,gate5_7_4); 
BDD	gate6_7_0	=  MkNor(s7_6, c6_6);
BDD	gate6_7_1	=  MkNor(s7_6, gate6_7_0);
BDD	gate6_7_2	=  MkNor(c6_6, gate6_7_0);
BDD	gate6_7_3	=  MkNor(gate6_7_1, gate6_7_2);
BDD	gate6_7_4	=  MkNor(gate6_7_3, pp6_7);
BDD	gate6_7_5	=  MkNor(gate6_7_3, gate6_7_4);
BDD	gate6_7_6	=  MkNor(pp6_7, gate6_7_4);
std::cout << "	s6_7	" << std::endl;
BDD	s6_7	=  MkNor(gate6_7_5, gate6_7_6);
BDD	c6_7	= MkNor(gate6_7_0,gate6_7_4); 
BDD	s7_7	= pp7_7;
BDD	gate0_8_0	= MkNot(s1_7); 
BDD	gate0_8_1	= MkNor(gate0_8_0, s1_7); 
BDD	gate0_8_2	= MkNot(gate0_8_0); 
BDD	gate0_8_3	= MkNor(gate0_8_1,gate0_8_2); 
BDD	gate0_8_4	=  MkNor(gate0_8_3, c0_7);
BDD	gate0_8_5	=  MkNor(gate0_8_3, gate0_8_4);
BDD	gate0_8_6	=  MkNor(c0_7, gate0_8_4);
std::cout << "	sum8	" << std::endl;
BDD	sum8	=  MkNor(gate0_8_5, gate0_8_6);
BDD	c0_8	= MkNor(gate0_8_0,gate0_8_4); 
BDD	gate1_8_0	=  MkNor(c1_7, c0_8);
BDD	gate1_8_1	=  MkNor(c1_7, gate1_8_0);
BDD	gate1_8_2	=  MkNor(c0_8, gate1_8_0);
BDD	gate1_8_3	=  MkNor(gate1_8_1, gate1_8_2);
BDD	gate1_8_4	=  MkNor(gate1_8_3, s2_7);
BDD	gate1_8_5	=  MkNor(gate1_8_3, gate1_8_4);
BDD	gate1_8_6	=  MkNor(s2_7, gate1_8_4);
BDD	sum9	=  MkNor(gate1_8_5, gate1_8_6);
std::cout << "	c1_8	" << std::endl;
BDD	c1_8	= MkNor(gate1_8_0,gate1_8_4); 
BDD	gate2_8_0	=  MkNor(c2_7, c1_8);
BDD	gate2_8_1	=  MkNor(c2_7, gate2_8_0);
BDD	gate2_8_2	=  MkNor(c1_8, gate2_8_0);
BDD	gate2_8_3	=  MkNor(gate2_8_1, gate2_8_2);
BDD	gate2_8_4	=  MkNor(gate2_8_3, s3_7);
BDD	gate2_8_5	=  MkNor(gate2_8_3, gate2_8_4);
BDD	gate2_8_6	=  MkNor(s3_7, gate2_8_4);
BDD	sum10	=  MkNor(gate2_8_5, gate2_8_6);
BDD	c2_8	= MkNor(gate2_8_0,gate2_8_4); 
std::cout << "	gate3_8_0	" << std::endl;
BDD	gate3_8_0	=  MkNor(c3_7, c2_8);
BDD	gate3_8_1	=  MkNor(c3_7, gate3_8_0);
BDD	gate3_8_2	=  MkNor(c2_8, gate3_8_0);
BDD	gate3_8_3	=  MkNor(gate3_8_1, gate3_8_2);
BDD	gate3_8_4	=  MkNor(gate3_8_3, s4_7);
BDD	gate3_8_5	=  MkNor(gate3_8_3, gate3_8_4);
BDD	gate3_8_6	=  MkNor(s4_7, gate3_8_4);
BDD	sum11	=  MkNor(gate3_8_5, gate3_8_6);
BDD	c3_8	= MkNor(gate3_8_0,gate3_8_4); 
BDD	gate4_8_0	=  MkNor(c4_7, c3_8);
std::cout << "	gate4_8_1	" << std::endl;
BDD	gate4_8_1	=  MkNor(c4_7, gate4_8_0);
BDD	gate4_8_2	=  MkNor(c3_8, gate4_8_0);
BDD	gate4_8_3	=  MkNor(gate4_8_1, gate4_8_2);
BDD	gate4_8_4	=  MkNor(gate4_8_3, s5_7);
BDD	gate4_8_5	=  MkNor(gate4_8_3, gate4_8_4);
BDD	gate4_8_6	=  MkNor(s5_7, gate4_8_4);
BDD	sum12	=  MkNor(gate4_8_5, gate4_8_6);
BDD	c4_8	= MkNor(gate4_8_0,gate4_8_4); 
BDD	gate5_8_0	=  MkNor(c5_7, c4_8);
BDD	gate5_8_1	=  MkNor(c5_7, gate5_8_0);
std::cout << "	gate5_8_2	" << std::endl;
BDD	gate5_8_2	=  MkNor(c4_8, gate5_8_0);
BDD	gate5_8_3	=  MkNor(gate5_8_1, gate5_8_2);
BDD	gate5_8_4	=  MkNor(gate5_8_3, s6_7);
BDD	gate5_8_5	=  MkNor(gate5_8_3, gate5_8_4);
BDD	gate5_8_6	=  MkNor(s6_7, gate5_8_4);
BDD	sum13	=  MkNor(gate5_8_5, gate5_8_6);
BDD	c5_8	= MkNor(gate5_8_0,gate5_8_4); 
BDD	gate6_8_0	=  MkNor(c6_7, c5_8);
BDD	gate6_8_1	=  MkNor(c6_7, gate6_8_0);
BDD	gate6_8_2	=  MkNor(c5_8, gate6_8_0);
std::cout << "	gate6_8_3	" << std::endl;
BDD	gate6_8_3	=  MkNor(gate6_8_1, gate6_8_2);
BDD	gate6_8_4	=  MkNor(gate6_8_3, s7_7);
BDD	gate6_8_5	=  MkNor(gate6_8_3, gate6_8_4);
BDD	gate6_8_6	=  MkNor(s7_7, gate6_8_4);
BDD	sum14	=  MkNor(gate6_8_5, gate6_8_6);
BDD	sum15	= MkNor(gate6_8_0,gate6_8_4); 

auto end = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(end - start);   

std::vector<BDD> sums = {sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7,
                        sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15};

unsigned int nodeCount = 0;
std::cout << "BDD sizes" << std::endl;
nodeCount = mgr.nodeCount(sums);
std::cout << "Duration: " << duration.count() << " Memory: " << nodeCount << std::endl;

for (unsigned int i = 0; i < 16; i++) {
  std::cout << "i: " << i << " size: " << sums[i].nodeCount() << std::endl;
}

// for (unsigned int i = 0; i < 16; i++) {
//   std::cout << "Printing sum " << i << std::endl;
//   sums[i].print(16, 2);
// }


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

int main () {

  c6288_8();
  return 0;
}
