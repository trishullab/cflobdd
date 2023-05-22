#ifndef MATRIX1234_NODE_GUARD
#define MATRIX1234_NODE_GUARD

//
//    Copyright (c) 2017, 2018 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "cflobdd_node.h"
#include "cflobdd_top_node_matmult_map.h"
#include "zero_indices_map.h"
#include <map>
#include <unordered_map>

namespace CFL_OBDD {

	typedef std::pair<long, long> INT_PAIR;
	enum VisitPosition { AVisit, BVisit, TopLevelVisit };
	extern int nodeNum;
	extern int cacheNodeNum;
	extern int notL1NodeNum;
	extern void clearMultMap();
	extern CFLOBDDNodeHandle MkIdRelationInterleavedNode(unsigned int level);
	extern CFLOBDDNodeHandle MkWalshInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkInverseReedMullerInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkNegationMatrixInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkPauliYInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkPauliZInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkSGateInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkCNOTInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkExchangeInterleavedNode(unsigned int i);
	extern CFLOBDDNodeHandle MkCNOTNode(unsigned int level, unsigned int n, long int controller, long int controlled);
	extern CFLOBDDNodeHandle MkCNOT2Node(unsigned int level, unsigned int n, long int controller, long int controlled);
	extern CFLOBDDNodeHandle MkCCNOTNode(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled);
	extern CFLOBDDNodeHandle MkCPGateNode(unsigned int level, long int controller, long int controlled);
	extern CFLOBDDNodeHandle MkSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num);
	extern CFLOBDDNodeHandle MkiSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num);
	extern CFLOBDDNodeHandle MkCSwapGateNode(unsigned int level, long int controller, long int i, long int j, int case_num);
	extern CFLOBDDNodeHandle MkCSwapGate2Node(unsigned int level, long int controller, long int i, long int j, int case_num);
	extern CFLOBDDNodeHandle MkCCPNode(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled);
	extern std::pair<CFLOBDDNodeHandle, int> MkRestrictNode(unsigned int level, std::string s);
	
	// Initialization routine that needs to be called before any call to MatrixProjectVoc23Node
	extern void Matrix1234InitializerNode();  // Empty for now

	// Matrix-related operations (on matrices with room for two extra vocabularies) ------------
	extern CFLOBDDNodeHandle MkWalshVoc13Node(unsigned int i);
	extern CFLOBDDNodeHandle MkWalshVoc12Node(unsigned int i);
	extern CFLOBDDNodeHandle MatrixShiftVoc43Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle MatrixShiftVoc42Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle MatrixShiftVocs13To24Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle MatrixShiftVocs12To34Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle MkDetensorConstraintInterleavedNode(unsigned int i);
	// extern CFLOBDDTopNodeLinearMapRefPtr MatrixProjectVoc23Node(CFLOBDDLinearMapMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh, VisitPosition position); // Vocabulary projection
	extern CFLOBDDNodeHandle ReverseColumnsNode(CFLOBDDNodeHandle nh);
	extern std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>
		MatrixTransposeNode(std::unordered_map<CFLOBDDNodeHandle, std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>, 
		CFLOBDDNodeHandle::CFLOBDDNodeHandle_Hash>& hashMap,
		CFLOBDDNodeHandle nh);

	extern CFLOBDDNodeHandle MatrixShiftToAConnectionNode(CFLOBDDNodeHandle c);
	extern CFLOBDDNodeHandle MatrixShiftToBConnectionNode(CFLOBDDNodeHandle c);

	// Subroutines for Discrete Fourier Transform
	extern CFLOBDDNodeHandle MkCFLOBDDMatrixEqVoc14Node(unsigned int i);
	extern CFLOBDDNodeHandle MkFourierDiagonalComponentNode(unsigned int i);
	extern CFLOBDDNodeHandle PromoteInterleavedTo12Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle Demote12ToInterleavedNode(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);
	extern CFLOBDDNodeHandle PromoteInterleavedTo13Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh);

	extern CFLOBDDNodeHandle SMatrixNode(std::string s);
	extern CFLOBDDNodeHandle MkDistinctionTwoVarsNode(int x, int y, unsigned int var_level, unsigned int matrix_level);

	// extern CFLOBDD_GENERAL AddMatrixRowsNode(CFLOBDDGeneralMapNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nhHandle);
	extern std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>> MultiplyOperationNode(CFLOBDDNodeHandle c);
	std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>> MultiplyOperationNodeInternal(CFLOBDDNodeHandle c, char position, unsigned int maxLevel, CFLOBDDReturnMapHandle cReturnMapHandle);
	std::vector<std::vector<std::pair<int, int>>> multiplyGeneralMap(int multiple, std::vector<std::vector<std::pair<int, int>>> &tmpReturnMap);
	std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>> addNodes(CFLOBDDNodeHandle n1, CFLOBDDNodeHandle n2, std::vector<std::vector<std::pair<int, int>>> n1GeneralMap, std::vector<std::vector<std::pair<int, int>>> n2GeneralMap);
	extern CFLOBDDNodeHandle MkMatrixMultiplyConstraintNode(unsigned int level);
	CFLOBDDNodeHandle MkMatrixMultiplyConstraintNodeInternal(unsigned int level);
	CFLOBDDNodeHandle MkRowWithOne(unsigned int bits, unsigned int rowNo, std::map<std::pair<int, int>, CFLOBDDNodeHandle>&  memoTable);
	void changeGeneralMap(std::vector<std::vector<std::pair<int, int>>> &generalMap, CFLOBDDReturnMapHandle returnMapHandle);
	extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4Node(
		std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash>& hashMap,
		CFLOBDDNodeHandle c1, CFLOBDDNodeHandle c2);
	extern CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4WithInfoNode(
		std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash>& hashMap,
		CFLOBDDNodeHandle c1, CFLOBDDNodeHandle c2, int c1_zero_index, int c2_zero_index);
 }

#endif

