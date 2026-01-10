// iscas85: c432 ------------------------------------------

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

static void c432(){
auto start = high_resolution_clock::now();
	Cudd mgr(0,0);
int coeff = 1;
int offset = 0;

// Inputs
// BDD gat1 = mgr.bddVar(coeff*0 + offset);
// BDD gat4 = mgr.bddVar(coeff*1 + offset);
// BDD gat11 = mgr.bddVar(coeff*2 + offset);
// BDD gat17 = mgr.bddVar(coeff*3 + offset);
// BDD gat24 = mgr.bddVar(coeff*4 + offset);
// BDD gat30 = mgr.bddVar(coeff*5 + offset);
// BDD gat37 = mgr.bddVar(coeff*6 + offset);
// BDD gat43 = mgr.bddVar(coeff*7 + offset);
// BDD gat50 = mgr.bddVar(coeff*8 + offset);
// BDD gat56 = mgr.bddVar(coeff*9 + offset);
// BDD gat63 = mgr.bddVar(coeff*10 + offset);
// BDD gat69 = mgr.bddVar(coeff*11 + offset);
// BDD gat76 = mgr.bddVar(coeff*12 + offset);
// BDD gat82 = mgr.bddVar(coeff*13 + offset);
// BDD gat89 = mgr.bddVar(coeff*14 + offset);
// BDD gat95 = mgr.bddVar(coeff*15 + offset);
// BDD gat102 = mgr.bddVar(coeff*16 + offset);
// BDD gat108 = mgr.bddVar(coeff*17 + offset);
// BDD gat8 = mgr.bddVar(coeff*18 + offset);
// BDD gat21 = mgr.bddVar(coeff*19 + offset);
// BDD gat34 = mgr.bddVar(coeff*20 + offset);
// BDD gat47 = mgr.bddVar(coeff*21 + offset);
// BDD gat60 = mgr.bddVar(coeff*22 + offset);
// BDD gat73 = mgr.bddVar(coeff*23 + offset);
// BDD gat86 = mgr.bddVar(coeff*24 + offset);
// BDD gat99 = mgr.bddVar(coeff*25 + offset);
// BDD gat112 = mgr.bddVar(coeff*26 + offset);
// BDD gat14 = mgr.bddVar(coeff*27 + offset);
// BDD gat27 = mgr.bddVar(coeff*28 + offset);
// BDD gat40 = mgr.bddVar(coeff*29 + offset);
// BDD gat53 = mgr.bddVar(coeff*30 + offset);
// BDD gat66 = mgr.bddVar(coeff*31 + offset);
// BDD gat79 = mgr.bddVar(coeff*32 + offset);
// BDD gat92 = mgr.bddVar(coeff*33 + offset);
// BDD gat105 = mgr.bddVar(coeff*34 + offset);
// BDD gat115 = mgr.bddVar(coeff*35 + offset);

// std::vector<int> var_order = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
// 							  16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
// 							  26, 27, 28, 29, 30, 31, 32, 33, 34, 35};

std::vector<int> var_order = {0, 1, 2, 4, 3, 5, 6, 8, 7, 9, 10, 12, 11, 13, 14, 16,
							  15, 17, 18, 20, 19, 21, 22, 24, 23, 25,
							  26, 28, 27, 29, 30, 32, 31, 33, 34, 35};

// 1 4 8 11 14 17 21 24 27 30 34 37 40 43 47 50 53 56 60 63 66 69 73 76 79 82 86 89 92 95 99 102 105 108 112 115
// 0 1 18 2 27 3 19 4 28 5 20 6 29 7 21 8 30 9 22 10 31 11 23 12 32 13 24 14 33 15 25 16 34 17 26
// (1, 0) (4, 1) (11, 2) (17, 3) (24, 4) (30, 5) (37, 6) (43, 7) (50, 8) (56, 9) (63, 10) (69, 11) (76, 12) (82, 13) (89, 14) (95, 15) (102, 16) (108, 17) (8, 18)
// (21, 19) (34, 20) (47, 21) (60, 22) (73, 23) (86, 24) (99, 25) (112, 26) (14, 27) (27, 28) (40, 29) (53, 30) (66, 31) (79, 32) (92, 33) (105, 34) (115, 35)
// std::vector<int> var_order = {0, 1, 18, 2, 27, 3, 19, 4, 28, 5, 20, 6, 29, 7, 21, 8, 30, 9, 22, 10, 31, 11, 23, 12, 32, 13, 24, 14, 33, 15, 25, 16, 34, 17, 26};

BDD gat1 = mgr.bddVar(var_order[0]);
BDD gat4 = mgr.bddVar(var_order[1]);
BDD gat8 = mgr.bddVar(var_order[2]);
BDD gat11 = mgr.bddVar(var_order[3]);
BDD gat14 = mgr.bddVar(var_order[4]);
BDD gat17 = mgr.bddVar(var_order[5]);
BDD gat21 = mgr.bddVar(var_order[6]);
BDD gat24 = mgr.bddVar(var_order[7]);
BDD gat27 = mgr.bddVar(var_order[8]);
BDD gat30 = mgr.bddVar(var_order[9]);
BDD gat34 = mgr.bddVar(var_order[10]);
BDD gat37 = mgr.bddVar(var_order[11]);
BDD gat40 = mgr.bddVar(var_order[12]);
BDD gat43 = mgr.bddVar(var_order[13]);
BDD gat47 = mgr.bddVar(var_order[14]);
BDD gat50 = mgr.bddVar(var_order[15]);
BDD gat53 = mgr.bddVar(var_order[16]);
BDD gat56 = mgr.bddVar(var_order[17]);
BDD gat60 = mgr.bddVar(var_order[18]);
BDD gat63 = mgr.bddVar(var_order[19]);
BDD gat66 = mgr.bddVar(var_order[20]);
BDD gat69 = mgr.bddVar(var_order[21]);
BDD gat73 = mgr.bddVar(var_order[22]);
BDD gat76 = mgr.bddVar(var_order[23]);
BDD gat79 = mgr.bddVar(var_order[24]);
BDD gat82 = mgr.bddVar(var_order[25]);
BDD gat86 = mgr.bddVar(var_order[26]);
BDD gat89 = mgr.bddVar(var_order[27]);
BDD gat92 = mgr.bddVar(var_order[28]);
BDD gat95 = mgr.bddVar(var_order[29]);
BDD gat99 = mgr.bddVar(var_order[30]);
BDD gat102 = mgr.bddVar(var_order[31]);
BDD gat105 = mgr.bddVar(var_order[32]);
BDD gat108 = mgr.bddVar(var_order[33]);
BDD gat112 = mgr.bddVar(var_order[34]);
BDD gat115 = mgr.bddVar(var_order[35]);

// Outputs: gat223,gat329,gat370,gat421,gat430,gat431,gat432;
   
// The circuit
cout << "Starting circuit" << endl;
cout << "gat118" << endl;
BDD gat118 = MkNot(gat1); // done
cout << "gat119" << endl;
BDD gat119 = MkNot(gat4); // done
cout << "gat122" << endl;
BDD gat122 = MkNot(gat11); // done
cout << "gat123" << endl;
BDD gat123 = MkNot(gat17); // done
cout << "gat126" << endl;
BDD gat126 = MkNot(gat24); // done
cout << "gat127" << endl;
BDD gat127 = MkNot(gat30); // done
cout << "gat130" << endl;
BDD gat130 = MkNot(gat37); // done
cout << "gat131" << endl;
BDD gat131 = MkNot(gat43); // done
cout << "gat134" << endl;
BDD gat134 = MkNot(gat50); // done
cout << "gat135" << endl;
BDD gat135 = MkNot(gat56); // done 
cout << "gat138" << endl;
BDD gat138 = MkNot(gat63); // done
cout << "gat139" << endl;
BDD gat139 = MkNot(gat69); // done
cout << "gat142" << endl;
BDD gat142 = MkNot(gat76); // done
cout << "gat143" << endl;
BDD gat143 = MkNot(gat82); // done
cout << "gat146" << endl;
BDD gat146 = MkNot(gat89); // done
cout << "gat147" << endl;
BDD gat147 = MkNot(gat95); // done
cout << "gat150" << endl;
BDD gat150 = MkNot(gat102); // done
cout << "gat151" << endl;
BDD gat151 = MkNot(gat108); // done
cout << "gat154" << endl;
BDD gat154 = MkNand(gat118, gat4);
cout << "gat157" << endl;
BDD gat157 = MkNor(gat8, gat119);
cout << "gat158" << endl;
BDD gat158 = MkNor(gat14, gat119);
cout << "gat159" << endl;
BDD gat159 = MkNand(gat122, gat17);
cout << "gat162" << endl;
BDD gat162 = MkNand(gat126, gat30);
cout << "gat165" << endl;
BDD gat165 = MkNand(gat130, gat43);
cout << "gat168" << endl;
BDD gat168 = MkNand(gat134, gat56);
cout << "gat171" << endl;
BDD gat171 = MkNand(gat138, gat69);
cout << "gat174" << endl;
BDD gat174 = MkNand(gat142, gat82);
cout << "gat177" << endl;
BDD gat177 = MkNand(gat146, gat95);
cout << "gat180" << endl;
BDD gat180 = MkNand(gat150, gat108);
cout << "gat183" << endl;
BDD gat183 = MkNor(gat21, gat123);
cout << "gat184" << endl;
BDD gat184 = MkNor(gat27, gat123);
cout << "gat185" << endl;
BDD gat185 = MkNor(gat34, gat127);
cout << "gat186" << endl;
BDD gat186 = MkNor(gat40, gat127);
cout << "gat187" << endl;
BDD gat187 = MkNor(gat47, gat131);
cout << "gat188" << endl;
BDD gat188 = MkNor(gat53, gat131);
cout << "gat189" << endl;
BDD gat189 = MkNor(gat60, gat135);
cout << "gat190" << endl;
BDD gat190 = MkNor(gat66, gat135);
cout << "gat191" << endl;
BDD gat191 = MkNor(gat73, gat139);
cout << "gat192" << endl;
BDD gat192 = MkNor(gat79, gat139);
cout << "gat193" << endl;
BDD gat193 = MkNor(gat86, gat143);
cout << "gat194" << endl;
BDD gat194 = MkNor(gat92, gat143);
cout << "gat195" << endl;
BDD gat195 = MkNor(gat99, gat147);
cout << "gat196" << endl;
BDD gat196 = MkNor(gat105, gat147);
cout << "gat197" << endl;
BDD gat197 = MkNor(gat112, gat151);
cout << "gat198" << endl;
BDD gat198 = MkNor(gat115, gat151); // done
cout << "gat199" << endl;
// BDD gat199 = MkAnd(9, gat154.root, gat159.root, gat162.root, gat165.root, gat168.root, gat171.root, gat174.root, gat177.root, gat180.root);
BDD gat199 = MkAnd(MkAnd(MkAnd(gat154, gat159), MkAnd(gat162, MkAnd(gat165, gat168))), MkAnd(MkAnd(gat171, gat174), MkAnd(gat177, gat180)));
cout << "gat203" << endl;
BDD gat203 = MkNot(gat199);
cout << "gat213" << endl;
BDD gat213 = MkNot(gat199);
cout << "gat223" << endl;
BDD gat223 = MkNot(gat199);
cout << "gat224" << endl;
BDD gat224 = MkExclusiveOr(gat203, gat154);
cout << "gat227" << endl;
BDD gat227 = MkExclusiveOr(gat203, gat159);
cout << "gat230" << endl;
BDD gat230 = MkExclusiveOr(gat203, gat162);
cout << "gat233" << endl;
BDD gat233 = MkExclusiveOr(gat203, gat165);
cout << "gat236" << endl;
BDD gat236 = MkExclusiveOr(gat203, gat168);
cout << "gat239" << endl;
BDD gat239 = MkExclusiveOr(gat203, gat171);
cout << "gat242" << endl;
BDD gat242 = MkNand(gat1, gat213);
cout << "gat243" << endl;
BDD gat243 = MkExclusiveOr(gat203, gat174);
cout << "gat246" << endl;
BDD gat246 = MkNand(gat213, gat11);
cout << "gat247" << endl;
BDD gat247 = MkExclusiveOr(gat203, gat177);
cout << "gat250" << endl;
BDD gat250 = MkNand(gat213, gat24);
cout << "gat251" << endl;
BDD gat251 = MkExclusiveOr(gat203, gat180);
cout << "gat254" << endl;
BDD gat254 = MkNand(gat213, gat37);
cout << "gat255" << endl;
BDD gat255 = MkNand(gat213, gat50);
cout << "gat256" << endl;
BDD gat256 = MkNand(gat213, gat63);
cout << "gat257" << endl;
BDD gat257 = MkNand(gat213, gat76);
cout << "gat258" << endl;
BDD gat258 = MkNand(gat213, gat89);
cout << "gat259" << endl;
BDD gat259 = MkNand(gat213, gat102);
cout << "gat260" << endl;
BDD gat260 = MkNand(gat224, gat157);
cout << "gat263" << endl;
BDD gat263 = MkNand(gat224, gat158);
cout << "gat264" << endl;
BDD gat264 = MkNand(gat227, gat183);
cout << "gat267" << endl;
BDD gat267 = MkNand(gat230, gat185);
cout << "gat270" << endl;
BDD gat270 = MkNand(gat233, gat187);
cout << "gat273" << endl;
BDD gat273 = MkNand(gat236, gat189);
cout << "gat276" << endl;
BDD gat276 = MkNand(gat239, gat191);
cout << "gat279" << endl;
BDD gat279 = MkNand(gat243, gat193);
cout << "gat282" << endl;
BDD gat282 = MkNand(gat247, gat195);
cout << "gat285" << endl;
BDD gat285 = MkNand(gat251, gat197);
cout << "gat288" << endl;
BDD gat288 = MkNand(gat227, gat184);
cout << "gat289" << endl;
BDD gat289 = MkNand(gat230, gat186);
cout << "gat290" << endl;
BDD gat290 = MkNand(gat233, gat188);
cout << "gat291" << endl;
BDD gat291 = MkNand(gat236, gat190);
cout << "gat292" << endl;
BDD gat292 = MkNand(gat239, gat192);
cout << "gat293" << endl;
BDD gat293 = MkNand(gat243, gat194);
cout << "gat294" << endl;
BDD gat294 = MkNand(gat247, gat196);
cout << "gat295" << endl;
BDD gat295 = MkNand(gat251, gat198);
cout << "gat296" << endl;
// BDD gat296 = MkAnd(9, gat260.root, gat264.root, gat267.root, gat270.root, gat273.root, gat276.root, gat279.root, gat282.root, gat285.root);
BDD gat296 = MkAnd(MkAnd(MkAnd(gat260, gat264), MkAnd(gat267, gat270)), MkAnd(MkAnd(gat273, gat276), MkAnd(MkAnd(gat279, gat282), gat285)));
cout << "gat300" << endl;
BDD gat300 = MkNot(gat263);
cout << "gat301" << endl;
BDD gat301 = MkNot(gat288);
cout << "gat302" << endl;
BDD gat302 = MkNot(gat289);
cout << "gat303" << endl;
BDD gat303 = MkNot(gat290);
cout << "gat304" << endl;
BDD gat304 = MkNot(gat291);
cout << "gat305" << endl;
BDD gat305 = MkNot(gat292);
cout << "gat306" << endl;
BDD gat306 = MkNot(gat293);
cout << "gat307" << endl;
BDD gat307 = MkNot(gat294);
cout << "gat308" << endl;
BDD gat308 = MkNot(gat295);
cout << "gat309" << endl;
BDD gat309 = MkNot(gat296);
cout << "gat319" << endl;
BDD gat319 = MkNot(gat296);
cout << "gat329" << endl;
BDD gat329 = MkNot(gat296);
cout << "gat330" << endl;
BDD gat330 = MkExclusiveOr(gat309, gat260);
cout << "gat331" << endl;
BDD gat331 = MkExclusiveOr(gat309, gat264);
cout << "gat332" << endl;
BDD gat332 = MkExclusiveOr(gat309, gat267);
cout << "gat333" << endl;
BDD gat333 = MkExclusiveOr(gat309, gat270);
cout << "gat334" << endl;
BDD gat334 = MkNand(gat8, gat319);
cout << "gat335" << endl;
BDD gat335 = MkExclusiveOr(gat309, gat273);
cout << "gat336" << endl;
BDD gat336 = MkNand(gat319, gat21);
cout << "gat337" << endl;
BDD gat337 = MkExclusiveOr(gat309, gat276);
cout << "gat338" << endl;
BDD gat338 = MkNand(gat319, gat34);
cout << "gat339" << endl;
BDD gat339 = MkExclusiveOr(gat309, gat279);
cout << "gat340" << endl;
BDD gat340 = MkNand(gat319, gat47);
cout << "gat341" << endl;
BDD gat341 = MkExclusiveOr(gat309, gat282);
cout << "gat342" << endl;
BDD gat342 = MkNand(gat319, gat60);
cout << "gat343" << endl;
BDD gat343 = MkExclusiveOr(gat309, gat285);
cout << "gat344" << endl;
BDD gat344 = MkNand(gat319, gat73);
cout << "gat345" << endl;
BDD gat345 = MkNand(gat319, gat86);
cout << "gat346" << endl;
BDD gat346 = MkNand(gat319, gat99);
cout << "gat347" << endl;
BDD gat347 = MkNand(gat319, gat112);
cout << "gat348" << endl;
BDD gat348 = MkNand(gat330, gat300);
cout << "gat349" << endl;
BDD gat349 = MkNand(gat331, gat301);
cout << "gat350" << endl;
BDD gat350 = MkNand(gat332, gat302);
cout << "gat351" << endl;
BDD gat351 = MkNand(gat333, gat303);
cout << "gat352" << endl;
BDD gat352 = MkNand(gat335, gat304);
cout << "gat353" << endl;
BDD gat353 = MkNand(gat337, gat305);
cout << "gat354" << endl;
BDD gat354 = MkNand(gat339, gat306);
cout << "gat355" << endl;
BDD gat355 = MkNand(gat341, gat307);
cout << "gat356" << endl;
BDD gat356 = MkNand(gat343, gat308);
cout << "gat357" << endl;
// BDD gat357 = MkAnd(9, gat348.root, gat349.root, gat350.root, gat351.root, gat352.root, gat353.root, gat354.root, gat355.root, gat356.root);
BDD gat357 = MkAnd(MkAnd(MkAnd(gat348, gat349), MkAnd(gat350, gat351)), MkAnd(MkAnd(gat352, gat353), MkAnd(MkAnd(gat354, gat355), gat356)));
cout << "gat360" << endl;
BDD gat360 = MkNot(gat357);
cout << "gat370" << endl;
BDD gat370 = MkNot(gat357);
cout << "gat371" << endl;
BDD gat371 = MkNand(gat14, gat360);
cout << "gat372" << endl;
BDD gat372 = MkNand(gat360, gat27);
cout << "gat373" << endl;
BDD gat373 = MkNand(gat360, gat40);
cout << "gat374" << endl;
BDD gat374 = MkNand(gat360, gat53);
cout << "gat375" << endl;
BDD gat375 = MkNand(gat360, gat66);
cout << "gat376" << endl;
BDD gat376 = MkNand(gat360, gat79);
cout << "gat377" << endl;
BDD gat377 = MkNand(gat360, gat92);
cout << "gat378" << endl;
BDD gat378 = MkNand(gat360, gat105);
cout << "gat379" << endl;
BDD gat379 = MkNand(gat360, gat115);
cout << "gat380" << endl;
BDD gat380 = MkNand(MkNand(gat4, gat242), MkNand(gat334, gat371));
cout << "gat381" << endl;
BDD gat381 = MkNand(MkNand(gat246, gat336), MkNand(gat372, gat17));
cout << "gat386" << endl;
BDD gat386 = MkNand(MkNand(gat250, gat338), MkNand(gat373, gat30));
cout << "gat393" << endl;
BDD gat393 = MkNand(MkNand(gat254, gat340), MkNand(gat374, gat43));
cout << "gat399" << endl;
BDD gat399 = MkNand(MkNand(gat255, gat342), MkNand(gat375, gat56));
cout << "gat404" << endl;
BDD gat404 = MkNand(MkNand(gat256, gat344), MkNand(gat376, gat69));
cout << "gat407" << endl;
BDD gat407 = MkNand(MkNand(gat257, gat345), MkNand(gat377, gat82));
cout << "gat411" << endl;
BDD gat411 = MkNand(MkNand(gat258, gat346), MkNand(gat378, gat95));
cout << "gat414" << endl;
BDD gat414 = MkNand(MkNand(gat259, gat347), MkNand(gat379, gat108));
cout << "gat415" << endl;
BDD gat415 = MkNot(gat380);
cout << "gat416" << endl;
// BDD gat416 = MkAnd(8, gat381.root, gat386.root, gat393.root, gat399.root, gat404.root, gat407.root, gat411.root, gat414.root);
BDD gat416 = MkAnd(MkAnd(MkAnd(gat381, gat386), MkAnd(gat393, gat399)), MkAnd(MkAnd(gat404, gat407), MkAnd(gat411, gat414)));
cout << "gat417" << endl;
BDD gat417 = MkNot(gat393);
cout << "gat418" << endl;
BDD gat418 = MkNot(gat404);
cout << "gat419" << endl;
BDD gat419 = MkNot(gat407);
cout << "gat420" << endl;
BDD gat420 = MkNot(gat411);
cout << "gat421" << endl;
BDD gat421 = MkNor(gat415, gat416);
cout << "gat422" << endl;
BDD gat422 = MkNand(gat386, gat417);
cout << "gat425" << endl;
BDD gat425 = MkNand(MkNand(gat386, gat393), MkNand(gat418, gat399));
cout << "gat428" << endl;
BDD gat428 = MkNand(gat399, MkNand(gat393, gat419));
cout << "gat429" << endl;
BDD gat429 = MkNand(MkNand(gat386, gat393), MkNand(gat407, gat420));
cout << "gat430" << endl;
BDD gat430 = MkNand(MkNand(gat381, gat386), MkNand(gat422, gat399));
cout << "gat431" << endl;
BDD gat431 = MkNand(MkNand(gat381, gat386), MkNand(gat425, gat428));
cout << "gat432" << endl;
BDD gat432 = MkNand(MkNand(gat381, gat422), MkNand(gat425, gat429));
auto end = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(end - start);

std::vector<BDD> bdds = {gat223, gat329, gat370, gat421, gat430, gat431, gat432};
   
unsigned int nodeCount = 0;
// cout << "BDD sizes" << endl;
// nodeCount += gat223.nodeCount();
// // cout << "nodeCount: " << nodeCount << endl;
// nodeCount += gat329.nodeCount();
// // cout << "nodeCount: " << nodeCount << endl;
// nodeCount += gat370.nodeCount();
// // cout << "nodeCount: " << nodeCount << endl;
// nodeCount += gat421.nodeCount();
// // cout << "nodeCount: " << nodeCount << endl;
// nodeCount += gat430.nodeCount();
// // cout << "nodeCount: " << nodeCount << endl;
// nodeCount += gat431.nodeCount();
// // cout << "nodeCount: " << nodeCount << endl;
// nodeCount += gat432.nodeCount();

nodeCount = mgr.nodeCount(bdds);

std::cout << "Duration: " << duration.count() << " Memory: " << nodeCount << std::endl;

for (unsigned int i = 0; i < bdds.size(); i++) {
	cout << "Output BDD " << i << " size: " << bdds[i].nodeCount() << endl;
}
}

int main () {

  c432();
//   ClearModules();
  return 0;
}

