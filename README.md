# CFLOBDDs: Context-Free-Language Ordered Binary Decision Diagrams
## Paper: https://arxiv.org/abs/2211.06818

CFLOBDDs are data structures that provide double-exponential compression of Boolean functions in the best case (and exponential compression over BDDs). They can be used as a plug-compatible replacement for BDDs. In this paper, we have used them in the context of simulating quantum algorithms.

### Installation
1. Git clone this repository (git clone https://github.com/trishullab/cflobdd.git .)
2. Install Boost C++ from https://www.boost.org/users/download/ (any version should work). The version used in the paper is 1.80.0.
3. Change the directory to the CFLOBDD subfolder (cd CFLOBDD/). The other folders contain code for testing the baselines mentioned in the paper.
4. Compile using the command (creates a.out executable): 
    g++ -g -std=c++17 -w -I. -I.<path_to_boost> -I.Solver/uwr/bit_vector/ -I.Solver/uwr/assert/ -I.Solver/uwr/matrix/ -I.Solver/uwr/parsing/ -lm *.cpp Solver/uwr/bit_vector/*.cpp Solver/uwr/parsing/*.cpp

### Code Structure:

The main CFLOBDD class is present in cflobdd_t.h with definitions of internal groupings and other classes in cflobdd_node.cpp.<br>
The matrix and vector operations are present in files matrix1234_* and vector_*.<br>
As all the classes are templatized, the classes are instantiated with different types (int, double, float_boost -> arbitary precision floating point).<br>
A CFLOBDD object (say 'int') generally follows :-> cflobdd_int.cpp, cflobdd_top_node_int.cpp (deals with terminal values) and then the core operations are in cflobdd_node.cpp. <br>
cflobdd_node.cpp (or matrix1234_node.cpp) is common to all data types are they handle proto-CFLOBDDs (see the paper for more details). <br>
All the tests are done via runTests function in tests_cfl.cpp which is the only function called in main.cpp.<br>
InitModules in tests_cfl.cpp should be run before creating any CFLOBDD object. InitModules initializes the cache. ClearModules should be invoked before exiting the program.
