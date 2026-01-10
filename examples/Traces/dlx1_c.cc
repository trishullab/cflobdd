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
// sss.1.0_trace: dlx1_c

static void dlx1_c()
{

auto start = high_resolution_clock::now();

Cudd mgr(0,0);
int coeff = 1;
int offset = 0;
int max_level = 6;

BDD ID_EX_RegWrite = mgr.bddVar(coeff*0 + offset);
BDD ID_EX_MemToReg = mgr.bddVar(coeff*1 + offset);
BDD _Taken_Branch_1_1 = mgr.bddVar(coeff*2 + offset);
BDD EX_MEM_Jump = mgr.bddVar(coeff*3 + offset);
BDD MEM_WB_RegWrite = mgr.bddVar(coeff*4 + offset);
BDD EX_MEM_RegWrite = mgr.bddVar(coeff*5 + offset);
BDD IF_ID_RegWrite = mgr.bddVar(coeff*6 + offset);
BDD ID_EX_Jump = mgr.bddVar(coeff*7 + offset);
BDD ID_EX_Branch = mgr.bddVar(coeff*8 + offset);
BDD TakeBranchALU_0 = mgr.bddVar(coeff*9 + offset);
BDD IF_ID_Flush = mgr.bddVar(coeff*10 + offset);
BDD e_1_1 = mgr.bddVar(coeff*11 + offset);
BDD IF_ID_UseData2 = mgr.bddVar(coeff*12 + offset);
BDD IF_ID_Branch = mgr.bddVar(coeff*13 + offset);
BDD IF_ID_MemWrite = mgr.bddVar(coeff*14 + offset);
BDD IF_ID_MemToReg = mgr.bddVar(coeff*15 + offset);
BDD e_2_1 = mgr.bddVar(coeff*16 + offset);
BDD e_2_2 = mgr.bddVar(coeff*17 + offset);
BDD e_3_2 = mgr.bddVar(coeff*18 + offset);
BDD EX_MEM_MemToReg = mgr.bddVar(coeff*19 + offset);
BDD MEM_WB_MemToReg = mgr.bddVar(coeff*20 + offset);
BDD e_3_1 = mgr.bddVar(coeff*21 + offset);
BDD e_1_3 = mgr.bddVar(coeff*22 + offset);
BDD e_1_4 = mgr.bddVar(coeff*23 + offset);
BDD UseData2_0 = mgr.bddVar(coeff*24 + offset);
BDD e_2_3 = mgr.bddVar(coeff*25 + offset);
BDD e_2_4 = mgr.bddVar(coeff*26 + offset);
BDD e_5_2 = mgr.bddVar(coeff*27 + offset);
BDD e_1_2 = mgr.bddVar(coeff*28 + offset);
BDD e_4_2 = mgr.bddVar(coeff*29 + offset);
BDD IF_ID_Jump = mgr.bddVar(coeff*30 + offset);
BDD TakeBranchALU_1 = mgr.bddVar(coeff*31 + offset);
BDD RegWrite_0 = mgr.bddVar(coeff*32 + offset);
BDD e_5_5 = mgr.bddVar(coeff*33 + offset);
BDD e_3_4 = mgr.bddVar(coeff*34 + offset);
BDD e_4_1 = mgr.bddVar(coeff*35 + offset);
BDD e_5_1 = mgr.bddVar(coeff*36 + offset);
BDD e_3_3 = mgr.bddVar(coeff*37 + offset);
BDD EX_MEM_MemWrite = mgr.bddVar(coeff*38 + offset);
BDD ID_EX_MemWrite = mgr.bddVar(coeff*39 + offset);
BDD e_4_3 = mgr.bddVar(coeff*40 + offset);
BDD e_4_4 = mgr.bddVar(coeff*41 + offset);
BDD e_5_4 = mgr.bddVar(coeff*42 + offset);
BDD MemWrite_0 = mgr.bddVar(coeff*43 + offset);
BDD e_5_3 = mgr.bddVar(coeff*44 + offset);
BDD Jump_0 = mgr.bddVar(coeff*45 + offset);
BDD Branch_0 = mgr.bddVar(coeff*46 + offset);
BDD TakeBranchALU_2 = mgr.bddVar(coeff*47 + offset);
BDD MemToReg_0 = mgr.bddVar(coeff*48 + offset);
BDD TakeBranchALU_3 = mgr.bddVar(coeff*49 + offset);
BDD TakeBranchALU_4 = mgr.bddVar(coeff*50 + offset);

// The circuit
BDD _squash_1_1 = MkOr(_Taken_Branch_1_1, EX_MEM_Jump);
BDD _squash_bar_1_1 = MkNot(_squash_1_1);
BDD _EX_Jump_1_1 = MkAnd(_squash_bar_1_1, ID_EX_Jump);
BDD _Taken_Branch_9_1 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_Branch), TakeBranchALU_0);
BDD _squash_9_1 = MkOr(_EX_Jump_1_1, _Taken_Branch_9_1);
BDD _squash_bar_9_1 = MkNot(_squash_9_1);
BDD _IF_ID_Flush_bar_1_1 = MkNot(IF_ID_Flush);
BDD _Reg2Used_1_1 = MkOr(MkOr(MkOr(IF_ID_UseData2, IF_ID_Branch), IF_ID_MemWrite), IF_ID_MemToReg);
BDD _temp_967 = MkAnd(_Reg2Used_1_1, e_2_1);
BDD _temp_968 = MkOr(e_1_1, _temp_967);
BDD _temp_969 = MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968);
BDD _temp_970 = MkOr(MkOr(MkOr(_Taken_Branch_1_1, IF_ID_Flush), EX_MEM_Jump), _temp_969);
BDD _temp_971 = MkNot(_temp_970);
BDD _temp_972 = MkAnd(MkAnd(IF_ID_Jump, _squash_bar_9_1), _temp_971);
BDD _temp_973 = MkAnd(MkAnd(MkAnd(IF_ID_Branch, _squash_bar_9_1), TakeBranchALU_1), _temp_971);
BDD _temp_974 = MkOr(_temp_972, _temp_973);
BDD _temp_975 = MkNot(_temp_974);
BDD _temp_976 = MkIfThenElse(_temp_969, IF_ID_Jump, Jump_0);
BDD _temp_977 = MkOr(MkOr(MkOr(MkOr(MkOr(_Taken_Branch_1_1, EX_MEM_Jump), _EX_Jump_1_1), _Taken_Branch_9_1), _temp_972), _temp_973);
BDD _temp_978 = MkNot(_temp_977);
BDD _temp_979 = MkAnd(_temp_976, _temp_978);
BDD _temp_980 = MkIfThenElse(_temp_969, IF_ID_Branch, Branch_0);
BDD _temp_368 = MkNot(ID_EX_MemToReg);
BDD _temp_981 = MkIfThenElse(_temp_969, e_2_1, e_3_1);
BDD _temp_982 = MkIfThenElse(_temp_969, e_2_2, e_3_2);
BDD _temp_987 = MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _squash_bar_9_1), _temp_971), _temp_982);
BDD _temp_988 = MkNot(_temp_987);
BDD _temp_989 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_2_1), _temp_981), _temp_988);
BDD _temp_268 = MkNot(MEM_WB_MemToReg);
BDD _temp_983 = MkIfThenElse(_temp_969, e_2_3, e_3_3);
BDD _temp_990 = MkAnd(EX_MEM_RegWrite, e_2_4);
BDD _temp_991 = MkNot(_temp_990);
BDD _temp_992 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), e_2_1);
BDD _temp_993 = MkNot(_temp_992);
BDD _temp_984 = MkIfThenElse(_temp_969, e_2_4, e_3_4);
BDD _temp_994 = MkAnd(EX_MEM_RegWrite, _temp_984);
BDD _temp_995 = MkNot(_temp_994);
BDD _temp_996 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_981);
BDD _temp_997 = MkNot(_temp_996);
BDD _temp_998 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_2_3), _temp_983), _temp_988), _temp_991), _temp_993), _temp_995), _temp_997);
BDD _temp_999 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_2_3), _temp_983), _temp_988), _temp_991), _temp_993), _temp_995), _temp_997);
BDD _temp_279 = MkNot(EX_MEM_MemToReg);
BDD _temp_1000 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_2_4), _temp_984), _temp_988), _temp_993), _temp_997);
BDD _temp_1001 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_2_4), _temp_984), _temp_988), _temp_993), _temp_997);
BDD _temp_1002 = MkAnd(MEM_WB_RegWrite, e_2_3);
BDD _temp_1003 = MkNot(_temp_1002);
BDD _temp_1004 = MkAnd(MEM_WB_RegWrite, _temp_983);
BDD _temp_1005 = MkNot(_temp_1004);
BDD _temp_1006 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_988), _temp_991), _temp_993), _temp_995), _temp_997), _temp_1003), _temp_1005);
BDD _temp_1007 = MkOr(MkOr(MkOr(MkOr(MkOr(_temp_989, _temp_998), _temp_999), _temp_1000), _temp_1001), _temp_1006);
BDD _temp_1008 = MkIfThenElse(_temp_969, e_1_1, e_4_1);
BDD _temp_1009 = MkIfThenElse(_temp_969, e_1_2, e_4_2);
BDD _temp_1014 = MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _squash_bar_9_1), _temp_971), _temp_1009);
BDD _temp_1015 = MkNot(_temp_1014);
BDD _temp_1016 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_1_1), _temp_1008), _temp_1015);
BDD _temp_1017 = MkAnd(MEM_WB_RegWrite, e_1_3);
BDD _temp_1018 = MkNot(_temp_1017);
BDD _temp_1019 = MkAnd(EX_MEM_RegWrite, e_1_4);
BDD _temp_1020 = MkNot(_temp_1019);
BDD _temp_1021 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), e_1_1);
BDD _temp_1022 = MkNot(_temp_1021);
BDD _temp_1010 = MkIfThenElse(_temp_969, e_1_3, e_4_3);
BDD _temp_1023 = MkAnd(MEM_WB_RegWrite, _temp_1010);
BDD _temp_1024 = MkNot(_temp_1023);
BDD _temp_1011 = MkIfThenElse(_temp_969, e_1_4, e_4_4);
BDD _temp_1025 = MkAnd(EX_MEM_RegWrite, _temp_1011);
BDD _temp_1026 = MkNot(_temp_1025);
BDD _temp_1027 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_1008);
BDD _temp_1028 = MkNot(_temp_1027);
BDD _temp_1029 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1015), _temp_1018), _temp_1020), _temp_1022), _temp_1024), _temp_1026), _temp_1028);
BDD _temp_1030 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_1_3), _temp_1010), _temp_1015), _temp_1020), _temp_1022), _temp_1026), _temp_1028);
BDD _temp_1031 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_1_3), _temp_1010), _temp_1015), _temp_1020), _temp_1022), _temp_1026), _temp_1028);
BDD _temp_1032 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_1_4), _temp_1011), _temp_1015), _temp_1022), _temp_1028);
BDD _temp_1033 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_1_4), _temp_1011), _temp_1015), _temp_1022), _temp_1028);
BDD _temp_1034 = MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1016, _temp_1029), _temp_1030), _temp_1031), _temp_1032), _temp_1033);
BDD _temp_1035 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1007), _temp_1034);
BDD _temp_1036 = MkIfThenElse(_temp_1035, TakeBranchALU_1, TakeBranchALU_2);
BDD _temp_1037 = MkAnd(MkAnd(_temp_978, _temp_980), _temp_1036);
BDD _temp_1038 = MkOr(_temp_979, _temp_1037);
BDD _temp_1039 = MkNot(_temp_1038);
BDD _ID_squash_57_1 = MkOr(MkOr(MkOr(MkOr(_Taken_Branch_1_1, IF_ID_Flush), EX_MEM_Jump), _EX_Jump_1_1), _Taken_Branch_9_1);
BDD _ID_squash_bar_57_1 = MkNot(_ID_squash_57_1);
BDD _ID_Jump_57_1 = MkAnd(IF_ID_Jump, _ID_squash_bar_57_1);
BDD _temp_1040 = MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_2_3), _temp_991), _temp_993);
BDD _temp_1041 = MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_2_3), _temp_991), _temp_993);
BDD _temp_1042 = MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_2_4), _temp_993);
BDD _temp_1043 = MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_2_4), _temp_993);
BDD _temp_1044 = MkAnd(MkAnd(_temp_991, _temp_993), _temp_1003);
BDD _temp_1045 = MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_2_1);
BDD _temp_1046 = MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1040, _temp_1041), _temp_1042), _temp_1043), _temp_1044), _temp_1045);
BDD _temp_1047 = MkAnd(MkAnd(_temp_1018, _temp_1020), _temp_1022);
BDD _temp_1048 = MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_1_3), _temp_1020), _temp_1022);
BDD _temp_1049 = MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_1_3), _temp_1020), _temp_1022);
BDD _temp_1050 = MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_1_4), _temp_1022);
BDD _temp_1051 = MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_1_4), _temp_1022);
BDD _temp_1052 = MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_1_1);
BDD _temp_1053 = MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1047, _temp_1048), _temp_1049), _temp_1050), _temp_1051), _temp_1052);
BDD _temp_1054 = MkAnd(_temp_1046, _temp_1053);
BDD _temp_1055 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), e_2_1), _temp_981), _temp_988);
BDD _temp_1056 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_989, _temp_998), _temp_999), _temp_1000), _temp_1001), _temp_1006), _temp_1055);
BDD _temp_1057 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), e_1_1), _temp_1008), _temp_1015);
BDD _temp_1058 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1016, _temp_1029), _temp_1030), _temp_1031), _temp_1032), _temp_1033), _temp_1057);
BDD _temp_1059 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1056), _temp_1058);
BDD _temp_1060 = MkIfThenElse(_temp_1059, TakeBranchALU_2, TakeBranchALU_3);
BDD _temp_1061 = MkIfThenElse(_temp_1054, TakeBranchALU_1, _temp_1060);
BDD _temp_1062 = MkAnd(MkAnd(IF_ID_Branch, _ID_squash_bar_57_1), _temp_1061);
BDD _temp_1063 = MkOr(_ID_Jump_57_1, _temp_1062);
BDD _temp_1064 = MkNot(_temp_1063);
BDD _temp_1065 = MkAnd(MkAnd(MkAnd(_squash_9_1, _temp_975), _temp_1039), _temp_1064);
BDD _temp_1066 = MkNot(_temp_969);
BDD _temp_1067 = MkOr(MkOr(_Taken_Branch_1_1, EX_MEM_Jump), _temp_1066);
BDD _temp_1068 = MkNot(_temp_1067);
BDD _temp_1069 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, _squash_bar_9_1), _temp_975), _temp_1039), _temp_1064), _temp_1068);
BDD _temp_1070 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_1_1, _squash_bar_9_1), _temp_975), _temp_1039), _temp_1064), _temp_1067);
BDD _temp_1071 = MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1034);
BDD _temp_1072 = MkNot(_temp_1071);
BDD _temp_1073 = MkNot(_temp_1053);
BDD _temp_1074 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1038), _temp_1058), _temp_1063), _temp_1072), _temp_1073);
BDD _temp_1075 = MkIfThenElse(_temp_1038, _temp_1071, _temp_974);
BDD _temp_1076 = MkAnd(MkAnd(_temp_1053, _temp_1063), _temp_1075);
BDD _temp_1077 = MkOr(MkOr(MkOr(MkOr(_temp_1065, _temp_1069), _temp_1070), _temp_1074), _temp_1076);
BDD _EX_MemWrite_1_1 = MkAnd(_squash_bar_1_1, ID_EX_MemWrite);
BDD _temp_595 = MkNot(_EX_MemWrite_1_1);
BDD _ID_MemWrite_57_1 = MkAnd(IF_ID_MemWrite, _ID_squash_bar_57_1);
BDD _temp_647 = MkNot(_ID_MemWrite_57_1);
BDD _temp_1078 = MkAnd(MkAnd(IF_ID_MemWrite, _squash_bar_9_1), _temp_971);
BDD _temp_1079 = MkNot(_temp_1078);
BDD _temp_1080 = MkIfThenElse(_temp_969, IF_ID_MemWrite, MemWrite_0);
BDD _temp_1081 = MkAnd(_temp_978, _temp_1080);
BDD _temp_1082 = MkNot(_temp_1081);
BDD _temp_1083 = MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_MemWrite, _temp_595), _temp_647), _temp_1079), _temp_1082);
BDD _temp_593 = MkNot(EX_MEM_MemWrite);
BDD _temp_1084 = MkAnd(MkAnd(MkAnd(MkAnd(_temp_593, _temp_595), _temp_647), _temp_1079), _temp_1082);
BDD _temp_1085 = MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_MemWrite), _temp_647), _temp_1079), _temp_1082);
BDD _temp_1086 = MkAnd(MkAnd(EX_MEM_MemWrite, _temp_595), _temp_1079);
BDD _temp_1087 = MkAnd(MkAnd(_temp_593, _temp_595), _temp_1079);
BDD _temp_1088 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_MemWrite), _temp_1079);
BDD _temp_1089 = MkOr(MkOr(_temp_1086, _temp_1087), _temp_1088);
BDD _temp_1090 = MkIfThenElse(_temp_969, IF_ID_UseData2, UseData2_0);
BDD _temp_1091 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, _squash_bar_1_1), ID_EX_RegWrite), _temp_368), e_2_1), _temp_981), _temp_988), _temp_1090);
BDD _temp_1092 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, MEM_WB_RegWrite), _temp_268), e_2_3), _temp_983), _temp_988), _temp_991), _temp_993), _temp_995), _temp_997), _temp_1090);
BDD _temp_1093 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, IF_ID_UseData2), MEM_WB_RegWrite), e_2_3), _temp_983), _temp_988), _temp_991), _temp_993), _temp_995), _temp_997), _temp_1090);
BDD _temp_1094 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, EX_MEM_RegWrite), _temp_279), e_2_4), _temp_984), _temp_988), _temp_993), _temp_997), _temp_1090);
BDD _temp_1095 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, EX_MEM_RegWrite), EX_MEM_MemToReg), e_2_4), _temp_984), _temp_988), _temp_993), _temp_997), _temp_1090);
BDD _temp_353 = MkNot(IF_ID_UseData2);
BDD _temp_1096 = MkNot(_temp_1090);
BDD _temp_1097 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_353), _temp_968), _temp_1096);
BDD _temp_1098 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, _IF_ID_Flush_bar_1_1), ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_988), _temp_991), _temp_993), _temp_995), _temp_997), _temp_1003), _temp_1005), _temp_1090);
BDD _temp_1099 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1091, _temp_1092), _temp_1093), _temp_1094), _temp_1095), _temp_1097), _temp_1098);
BDD _temp_1100 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1007), _temp_1034), _temp_1089), _temp_1099);
BDD _temp_1101 = MkNot(_temp_1100);
BDD _temp_1102 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, MEM_WB_RegWrite), _temp_268), e_2_3), _temp_991), _temp_993);
BDD _temp_1103 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, IF_ID_UseData2), MEM_WB_RegWrite), e_2_3), _temp_991), _temp_993);
BDD _temp_1104 = MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, EX_MEM_RegWrite), _temp_279), e_2_4), _temp_993);
BDD _temp_1105 = MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, EX_MEM_RegWrite), EX_MEM_MemToReg), e_2_4), _temp_993);
BDD _temp_1106 = MkAnd(MkAnd(MkAnd(IF_ID_UseData2, _temp_991), _temp_993), _temp_1003);
BDD _temp_1107 = MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, _squash_bar_1_1), ID_EX_RegWrite), _temp_368), e_2_1);
BDD _temp_1108 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_353, _temp_1102), _temp_1103), _temp_1104), _temp_1105), _temp_1106), _temp_1107);
BDD _temp_1109 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1034), _temp_1053), _temp_1099), _temp_1108);
BDD _temp_1110 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1034), _temp_1099);
BDD _temp_1111 = MkNot(_temp_1110);
BDD _temp_1112 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_UseData2, _squash_bar_1_1), ID_EX_RegWrite), ID_EX_MemToReg), e_2_1), _temp_981), _temp_988), _temp_1090);
BDD _temp_1113 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1091, _temp_1092), _temp_1093), _temp_1094), _temp_1095), _temp_1097), _temp_1098), _temp_1112);
BDD _temp_1114 = MkAnd(_temp_1053, _temp_1108);
BDD _temp_1115 = MkNot(_temp_1114);
BDD _temp_1116 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1058), _temp_1111), _temp_1113), _temp_1115);
BDD _temp_1117 = MkOr(_temp_1109, _temp_1116);
BDD _temp_1118 = MkAnd(MkAnd(_temp_1046, _temp_1053), _temp_1108);
BDD _temp_1119 = MkNot(_temp_1118);
BDD _temp_1120 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemWrite, _IF_ID_Flush_bar_1_1), ID_EX_RegWrite), ID_EX_MemToReg), _ID_squash_bar_57_1), _temp_968), _temp_978), _temp_1056), _temp_1080), _temp_1089), _temp_1101), _temp_1117), _temp_1119);
BDD _temp_1121 = MkIfThenElse(_temp_1081, _temp_1100, _temp_1078);
BDD _temp_1122 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemWrite, _ID_squash_bar_57_1), _temp_1046), _temp_1053), _temp_1108), _temp_1121);
BDD _temp_1123 = MkOr(MkOr(MkOr(MkOr(_temp_1083, _temp_1084), _temp_1085), _temp_1120), _temp_1122);
BDD _temp_1124 = MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _squash_bar_9_1), _temp_971), e_5_2);
BDD _temp_1125 = MkNot(_temp_1124);
BDD _temp_1126 = MkIfThenElse(_temp_969, IF_ID_RegWrite, RegWrite_0);
BDD _temp_1128 = MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), e_5_2);
BDD _temp_1129 = MkAnd(_temp_1066, e_5_5);
BDD _temp_1130 = MkOr(_temp_1128, _temp_1129);
BDD _temp_1131 = MkAnd(MkAnd(_temp_978, _temp_1126), _temp_1130);
BDD _temp_1132 = MkNot(_temp_1131);
BDD _temp_1133 = MkAnd(MkAnd(IF_ID_RegWrite, _ID_squash_bar_57_1), e_5_2);
BDD _temp_1134 = MkNot(_temp_1133);
BDD _temp_1135 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_5_1), _temp_1125), _temp_1132), _temp_1134);
BDD _temp_1136 = MkAnd(EX_MEM_RegWrite, e_5_4);
BDD _temp_1137 = MkNot(_temp_1136);
BDD _temp_1138 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), e_5_1);
BDD _temp_1139 = MkNot(_temp_1138);
BDD _temp_1140 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_5_3), _temp_1125), _temp_1132), _temp_1134), _temp_1137), _temp_1139);
BDD _temp_1141 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_5_3), _temp_1125), _temp_1132), _temp_1134), _temp_1137), _temp_1139);
BDD _temp_1142 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_5_4), _temp_1125), _temp_1132), _temp_1134), _temp_1139);
BDD _temp_1143 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_5_4), _temp_1125), _temp_1132), _temp_1134), _temp_1139);
BDD _temp_363 = MkNot(IF_ID_MemToReg);
BDD _temp_1144 = MkIfThenElse(_temp_969, IF_ID_MemToReg, MemToReg_0);
BDD _temp_1145 = MkNot(_temp_1144);
BDD _temp_1146 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1034), _temp_1099), _temp_1145);
BDD _temp_1147 = MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _squash_bar_9_1), _temp_363), _temp_971), e_5_2);
BDD _temp_1148 = MkIfThenElse(_temp_1131, _temp_1146, _temp_1147);
BDD _temp_1149 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _ID_squash_bar_57_1), _temp_363), _temp_1053), _temp_1108), e_5_2), _temp_1148);
BDD _temp_1150 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), _temp_968), _temp_1007), _temp_1034), _temp_1089), _temp_1099), _temp_1144);
BDD _temp_1151 = MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, IF_ID_RegWrite), _squash_bar_9_1), _temp_971), e_5_2);
BDD _temp_1152 = MkIfThenElse(_temp_1131, _temp_1150, _temp_1151);
BDD _temp_1153 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, IF_ID_RegWrite), _ID_squash_bar_57_1), _temp_1046), _temp_1053), _temp_1108), e_5_2), _temp_1152);
BDD _temp_1154 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), e_5_1), _temp_1125), _temp_1132), _temp_1134);
BDD _temp_1155 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), IF_ID_RegWrite), _ID_squash_bar_57_1), _temp_363), _temp_968), _temp_978), _temp_1058), _temp_1111), _temp_1113), _temp_1115), e_5_2), _temp_1126), _temp_1130), _temp_1145);
BDD _temp_1156 = MkAnd(MEM_WB_RegWrite, e_5_3);
BDD _temp_1157 = MkNot(_temp_1156);
BDD _temp_1158 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_1125, _temp_1132), _temp_1134), _temp_1137), _temp_1139), _temp_1157);
BDD _temp_1159 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, _IF_ID_Flush_bar_1_1), ID_EX_RegWrite), ID_EX_MemToReg), IF_ID_RegWrite), _ID_squash_bar_57_1), _temp_968), _temp_978), _temp_1056), _temp_1089), _temp_1101), _temp_1117), _temp_1119), e_5_2), _temp_1126), _temp_1130), _temp_1144);
BDD _temp_1160 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1135, _temp_1140), _temp_1141), _temp_1142), _temp_1143), _temp_1149), _temp_1153), _temp_1154), _temp_1155), _temp_1158), _temp_1159);
BDD _temp_1161 = MkAnd(MkAnd(_temp_1077, _temp_1123), _temp_1160);
BDD _temp_1162 = MkAnd(RegWrite_0, e_5_5);
BDD _temp_1163 = MkNot(_temp_1162);
BDD _temp_1164 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_5_1), _temp_1125), _temp_1132), _temp_1134), _temp_1163);
BDD _temp_1165 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_5_3), _temp_1125), _temp_1132), _temp_1134), _temp_1137), _temp_1139), _temp_1163);
BDD _temp_1166 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_5_3), _temp_1125), _temp_1132), _temp_1134), _temp_1137), _temp_1139), _temp_1163);
BDD _temp_1167 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_5_4), _temp_1125), _temp_1132), _temp_1134), _temp_1139), _temp_1163);
BDD _temp_1168 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_5_4), _temp_1125), _temp_1132), _temp_1134), _temp_1139), _temp_1163);
BDD _temp_1169 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _ID_squash_bar_57_1), _temp_363), _temp_1053), _temp_1108), e_5_2), _temp_1148), _temp_1163);
BDD _temp_1170 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, IF_ID_RegWrite), _ID_squash_bar_57_1), _temp_1046), _temp_1053), _temp_1108), e_5_2), _temp_1152), _temp_1163);
BDD _temp_1171 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), e_5_1), _temp_1125), _temp_1132), _temp_1134), _temp_1163);
BDD _temp_893 = MkNot(MemToReg_0);
BDD _temp_1172 = MkAnd(MkAnd(IF_ID_RegWrite, _ID_squash_bar_57_1), e_4_2);
BDD _temp_1173 = MkNot(_temp_1172);
BDD _temp_1174 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_4_1), _temp_1008), _temp_1015), _temp_1173);
BDD _temp_1175 = MkAnd(EX_MEM_RegWrite, e_4_4);
BDD _temp_1176 = MkNot(_temp_1175);
BDD _temp_1177 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), e_4_1);
BDD _temp_1178 = MkNot(_temp_1177);
BDD _temp_1179 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_4_3), _temp_1010), _temp_1015), _temp_1026), _temp_1028), _temp_1173), _temp_1176), _temp_1178);
BDD _temp_1180 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_4_3), _temp_1010), _temp_1015), _temp_1026), _temp_1028), _temp_1173), _temp_1176), _temp_1178);
BDD _temp_1181 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_4_4), _temp_1011), _temp_1015), _temp_1028), _temp_1173), _temp_1178);
BDD _temp_1182 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_4_4), _temp_1011), _temp_1015), _temp_1028), _temp_1173), _temp_1178);
BDD _temp_1183 = MkAnd(MEM_WB_RegWrite, e_4_3);
BDD _temp_1184 = MkNot(_temp_1183);
BDD _temp_1185 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_1015, _temp_1024), _temp_1026), _temp_1028), _temp_1066), _temp_1173), _temp_1176), _temp_1178), _temp_1184);
BDD _temp_1186 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _squash_bar_9_1), _ID_squash_bar_57_1), _temp_363), _temp_971), e_4_2), _temp_1009), _temp_1053), _temp_1108);
BDD _temp_1187 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, IF_ID_RegWrite), _squash_bar_9_1), _ID_squash_bar_57_1), _temp_971), e_4_2), _temp_1009), _temp_1046), _temp_1053), _temp_1108);
BDD _temp_1188 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), e_4_1), _temp_1008), _temp_1015), _temp_1173);
BDD _temp_1189 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1174, _temp_1179), _temp_1180), _temp_1181), _temp_1182), _temp_1185), _temp_1186), _temp_1187), _temp_1188);
BDD _temp_1190 = MkAnd(MkAnd(IF_ID_RegWrite, _ID_squash_bar_57_1), e_3_2);
BDD _temp_1191 = MkNot(_temp_1190);
BDD _temp_1192 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), UseData2_0), _temp_368), e_3_1), _temp_981), _temp_988), _temp_1090), _temp_1191);
BDD _temp_1193 = MkAnd(EX_MEM_RegWrite, e_3_4);
BDD _temp_1194 = MkNot(_temp_1193);
BDD _temp_1195 = MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), e_3_1);
BDD _temp_1196 = MkNot(_temp_1195);
BDD _temp_1197 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, UseData2_0), _temp_268), e_3_3), _temp_983), _temp_988), _temp_995), _temp_997), _temp_1090), _temp_1191), _temp_1194), _temp_1196);
BDD _temp_1198 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), UseData2_0), e_3_3), _temp_983), _temp_988), _temp_995), _temp_997), _temp_1090), _temp_1191), _temp_1194), _temp_1196);
BDD _temp_1199 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, UseData2_0), _temp_279), e_3_4), _temp_984), _temp_988), _temp_997), _temp_1090), _temp_1191), _temp_1196);
BDD _temp_1200 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, UseData2_0), EX_MEM_MemToReg), e_3_4), _temp_984), _temp_988), _temp_997), _temp_1090), _temp_1191), _temp_1196);
BDD _temp_1201 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, UseData2_0), _squash_bar_9_1), _ID_squash_bar_57_1), _temp_363), _temp_971), e_3_2), _temp_982), _temp_1053), _temp_1090), _temp_1108);
BDD _temp_1202 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, IF_ID_RegWrite), UseData2_0), _squash_bar_9_1), _ID_squash_bar_57_1), _temp_971), e_3_2), _temp_982), _temp_1046), _temp_1053), _temp_1090), _temp_1108);
BDD _temp_1203 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), UseData2_0), e_3_1), _temp_981), _temp_988), _temp_1090), _temp_1191);
BDD _temp_1204 = MkAnd(MEM_WB_RegWrite, e_3_3);
BDD _temp_1205 = MkNot(_temp_1204);
BDD _temp_1206 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(UseData2_0, _temp_988), _temp_995), _temp_997), _temp_1005), _temp_1066), _temp_1090), _temp_1191), _temp_1194), _temp_1196), _temp_1205);
BDD _temp_848 = MkNot(UseData2_0);
BDD _temp_1207 = MkAnd(MkAnd(_temp_848, _temp_1066), _temp_1096);
BDD _temp_1208 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1192, _temp_1197), _temp_1198), _temp_1199), _temp_1200), _temp_1201), _temp_1202), _temp_1203), _temp_1206), _temp_1207);
BDD _temp_1209 = MkAnd(MkAnd(MkAnd(_temp_893, _temp_1066), _temp_1189), _temp_1208);
BDD _temp_1210 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_IF_ID_Flush_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), IF_ID_RegWrite), _ID_squash_bar_57_1), _temp_363), _temp_968), _temp_1058), _temp_1113), _temp_1115), e_5_2);
BDD _temp_1211 = MkIfThenElse(_temp_1162, _temp_1209, _temp_1210);
BDD _temp_1212 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_978, _temp_1111), _temp_1126), _temp_1130), _temp_1145), _temp_1211);
BDD _temp_1213 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_1125, _temp_1132), _temp_1134), _temp_1137), _temp_1139), _temp_1157), _temp_1163);
BDD _temp_1214 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), _temp_368), e_3_1), _temp_981), _temp_988), _temp_1191);
BDD _temp_1215 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_RegWrite, _temp_268), e_3_3), _temp_983), _temp_988), _temp_995), _temp_997), _temp_1191), _temp_1194), _temp_1196);
BDD _temp_1216 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MEM_WB_MemToReg, MEM_WB_RegWrite), e_3_3), _temp_983), _temp_988), _temp_995), _temp_997), _temp_1191), _temp_1194), _temp_1196);
BDD _temp_1217 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, _temp_279), e_3_4), _temp_984), _temp_988), _temp_997), _temp_1191), _temp_1196);
BDD _temp_1218 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_RegWrite, EX_MEM_MemToReg), e_3_4), _temp_984), _temp_988), _temp_997), _temp_1191), _temp_1196);
BDD _temp_1219 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_RegWrite, _squash_bar_9_1), _ID_squash_bar_57_1), _temp_363), _temp_971), e_3_2), _temp_982), _temp_1053), _temp_1108);
BDD _temp_1220 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, IF_ID_RegWrite), _squash_bar_9_1), _ID_squash_bar_57_1), _temp_971), e_3_2), _temp_982), _temp_1046), _temp_1053), _temp_1108);
BDD _temp_1221 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_RegWrite), ID_EX_MemToReg), e_3_1), _temp_981), _temp_988), _temp_1191);
BDD _temp_1222 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_988, _temp_995), _temp_997), _temp_1005), _temp_1066), _temp_1191), _temp_1194), _temp_1196), _temp_1205);
BDD _temp_1223 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1214, _temp_1215), _temp_1216), _temp_1217), _temp_1218), _temp_1219), _temp_1220), _temp_1221), _temp_1222);
BDD _temp_1224 = MkAnd(MkAnd(MkAnd(EX_MEM_MemWrite, _temp_595), _temp_647), _temp_1079);
BDD _temp_1225 = MkAnd(MkAnd(MkAnd(_temp_593, _temp_595), _temp_647), _temp_1079);
BDD _temp_1226 = MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_MemWrite), _temp_647), _temp_1079);
BDD _temp_1227 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemWrite, _squash_bar_9_1), _ID_squash_bar_57_1), _temp_971), _temp_1046), _temp_1053), _temp_1108);
BDD _temp_1228 = MkOr(MkOr(MkOr(_temp_1224, _temp_1225), _temp_1226), _temp_1227);
BDD _temp_1229 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MemToReg_0, _temp_1066), _temp_1111), _temp_1189), _temp_1208), _temp_1223), _temp_1228);
BDD _temp_1230 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemToReg, _IF_ID_Flush_bar_1_1), ID_EX_RegWrite), ID_EX_MemToReg), IF_ID_RegWrite), _ID_squash_bar_57_1), _temp_968), _temp_1056), _temp_1089), _temp_1117), _temp_1119), e_5_2);
BDD _temp_1231 = MkIfThenElse(_temp_1162, _temp_1229, _temp_1230);
BDD _temp_1232 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_978, _temp_1101), _temp_1126), _temp_1130), _temp_1144), _temp_1231);
BDD _temp_1233 = MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(MkOr(_temp_1164, _temp_1165), _temp_1166), _temp_1167), _temp_1168), _temp_1169), _temp_1170), _temp_1171), _temp_1212), _temp_1213), _temp_1232);
BDD _temp_933 = MkNot(MemWrite_0);
BDD _temp_1234 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(EX_MEM_MemWrite, _temp_595), _temp_647), _temp_933), _temp_1079), _temp_1082);
BDD _temp_1235 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_593, _temp_595), _temp_647), _temp_933), _temp_1079), _temp_1082);
BDD _temp_1236 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, ID_EX_MemWrite), _temp_647), _temp_933), _temp_1079), _temp_1082);
BDD _temp_1237 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_temp_1066, _temp_1111), _temp_1189), _temp_1208), _temp_1223), _temp_1228);
BDD _temp_1238 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemWrite, _IF_ID_Flush_bar_1_1), ID_EX_RegWrite), ID_EX_MemToReg), _ID_squash_bar_57_1), _temp_968), _temp_1056), _temp_1089), _temp_1117), _temp_1119);
BDD _temp_1239 = MkIfThenElse(MemWrite_0, _temp_1237, _temp_1238);
BDD _temp_1240 = MkAnd(MkAnd(MkAnd(_temp_978, _temp_1080), _temp_1101), _temp_1239);
BDD _temp_1241 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(IF_ID_MemWrite, _ID_squash_bar_57_1), _temp_933), _temp_1046), _temp_1053), _temp_1108), _temp_1121);
BDD _temp_1242 = MkOr(MkOr(MkOr(MkOr(_temp_1234, _temp_1235), _temp_1236), _temp_1240), _temp_1241);
BDD _temp_1243 = MkAnd(MkAnd(_temp_1066, _temp_1189), _temp_1223);
BDD _temp_1244 = MkIfThenElse(_temp_1243, TakeBranchALU_2, TakeBranchALU_4);
BDD _temp_1245 = MkAnd(Branch_0, _temp_1244);
BDD _temp_1246 = MkOr(Jump_0, _temp_1245);
BDD _temp_1247 = MkNot(_temp_1246);
BDD _temp_1248 = MkAnd(MkAnd(MkAnd(MkAnd(MkAnd(_squash_bar_1_1, _squash_bar_9_1), _temp_975), _temp_1039), _temp_1067), _temp_1247);
BDD _temp_1249 = MkAnd(MkAnd(MkAnd(MkAnd(_temp_1038, _temp_1066), _temp_1072), _temp_1189), _temp_1246);
BDD _temp_1250 = MkOr(_temp_1248, _temp_1249);
BDD _temp_1251 = MkAnd(MkAnd(_temp_1233, _temp_1242), _temp_1250);
BDD _temp_1252 = MkOr(_temp_1161, _temp_1251);

BDD true_value = mgr.bddOne();
cout << ((_temp_1252 == true_value) ? "equal" : "not equal") << endl;

auto end = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>(end - start); 
std::cout << "Duration: " << duration.count() << std::endl;

}

int main () {

  dlx1_c();
  return 0;
}
