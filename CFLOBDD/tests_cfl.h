#ifndef TESTS_CFL_GUARD
#define TESTS_CFL_GUARD

// #include "cflobdd_int.h"
// #include "matrix1234_int.h"
// #include "vector_int.h"

namespace CFL_OBDD {
class CFLTests
{
    public:
        static void testInverseReedMull();
        static void testWalsh();
        static void testStepFunction();
        static void testPulseFunction();
        static void testIscas85();
        static void testPermute();
        static void testTopNodes();
        static void testCanonicalness();
        static void testAnd();
		static void testPlus();
		static void testTimes();
		static void testMatrixMultiplication();
        static void testSatisfyingAssignments();
        static void testAllAssignments();
        static void test3();
        static void testMkDetensorConstraintInterleaved();
		static void testMkCFLOBDDMatrixEqVoc14();
		static void testMkFourierMatrix(int p);
        static void testMkIdRelationInterleaved();
        static void testProbability();
        static void testParity();
        static void test1();
        static void test2();
        static void testMaxLevelGreaterThan10();
		static void testMkAdditionInterleaved();
        static int test_demorgans(void);
        static int test_irm4(void);
        static int test_arbitrary_step_functions(unsigned int mlev);
        static int test_restrict_exists_and_forall(void);
        // static void factoringTest(unsigned int testNo, CFL_OBDD::CFLOBDD rel, unsigned int product);
		static void testKaratsuba();
		static void testBasisVector();
		static void testBasicOperations(int size);
		static void testShorsAlgo();
		static void testShortestPath();
		static void testGHZAlgo(int size);
		static void testBVAlgo(int size, int seed);
		static void testDJAlgo(int size, int seed);
		static void testGroversAlgo(int size, int seed);
		static void testSimonsAlgo(int size, int seed);
        static void testSimonsAlgoNew(int size);
        static void testXOR(int size);
        static void testMatMul(int size);
        static void testQFT(int size, int seed);

        static void testWeightedOps(unsigned int size);
		static void testGHZAlgo_W(int size);
		static void testBVAlgo_W(int size, int seed);
		static void testDJAlgo_W(int size, int seed);
        static void testGroversAlgo_W(int size, int seed);
        static void testQFT_W(int size, int seed);
		static void testShorsAlgo_W(int N, int a);
	    
        static void InitModules();
	    static void ClearModules();
		static bool runTests(const char *arg, int size = 0, int seed = 0);
};
}

#endif
