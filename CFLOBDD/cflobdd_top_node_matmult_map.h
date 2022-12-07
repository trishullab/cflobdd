#ifndef CFLOBDD_TOP_NODE_MATMULT_MAP_GUARD
#define CFLOBDD_TOP_NODE_MATMULT_MAP_GUARD

//
//    Copyright (c) 2017 Thomas W. Reps
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


// #include <iostream>
// #include <fstream>
// #include <complex>
#include "matmult_map.h"
#include "connectionT.h"
#include "cflobdd_top_node_t.h"
// #include "matrix1234_complex_double_top_node.h"
// #include "matrix1234_double_top_node.h"
// #include "matrix1234_float_boost_top_node.h"

namespace CFL_OBDD {

	typedef ReturnMapHandle<MatMultMapHandle> CFLOBDDMatMultMapHandle;
	typedef ConnectionT<CFLOBDDMatMultMapHandle> MatMultMapConnection;

	typedef CFLOBDDTopNodeT<MatMultMapHandle> CFLOBDDTopNodeMatMultMap;
	typedef CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeMatMultMapRefPtr;

	class MatMultPair{
	public:
		CFLOBDDNodeHandle m1;
		CFLOBDDNodeHandle m2;
		MatMultPair(CFLOBDDNodeHandle p1, CFLOBDDNodeHandle p2)
		{
			m1 = p1;
			m2 = p2;
		}

		struct MatMultPairHash {
			size_t operator()(const MatMultPair& p) const
			{
				CFLOBDDNodeHandle t1 = p.m1;
				CFLOBDDNodeHandle t2 = p.m2;
				auto hash1 = t1.Hash(997);
				auto hash2 = t2.Hash(997);
				return 117 * (hash1 + 1) + hash2;
			}
		};

		bool operator==(const MatMultPair& p) const
		{
			CFLOBDDNodeHandle m11 = m1;
			CFLOBDDNodeHandle m12 = m2;
			CFLOBDDNodeHandle m21 = p.m1;
			CFLOBDDNodeHandle m22 = p.m2;
			return (m11 == m21) && (m12 == m22);
		}
	};

	class MatMultPairWithInfo {
	public:
		CFLOBDDNodeHandle m1;
		CFLOBDDNodeHandle m2;
		int c1_index;
		int c2_index;
		MatMultPairWithInfo(CFLOBDDNodeHandle p1, CFLOBDDNodeHandle p2, int c1, int c2)
		{
			m1 = p1;
			m2 = p2;
			c1_index = c1;
			c2_index = c2;
		}

		struct MatMultPairWithInfoHash {
			size_t operator()(const MatMultPairWithInfo& p) const
			{
				CFLOBDDNodeHandle t1 = p.m1;
				CFLOBDDNodeHandle t2 = p.m2;
				auto hash1 = t1.Hash(997);
				auto hash2 = t2.Hash(997);
				return 117 * (hash1 + 1) + 97 * 97 * hash2 + (size_t)(97 * p.c1_index + p.c2_index);
			}
		};

		bool operator==(const MatMultPairWithInfo& p) const
		{
			CFLOBDDNodeHandle m11 = m1;
			CFLOBDDNodeHandle m12 = m2;
			CFLOBDDNodeHandle m21 = p.m1;
			CFLOBDDNodeHandle m22 = p.m2;
			return (m11 == m21) && (m12 == m22) && (c1_index == p.c1_index) && (c2_index == p.c2_index);
		}
	};

	class MatMultAddPair{
	public:
		CFLOBDDTopNodeMatMultMapRefPtr m1;
		CFLOBDDTopNodeMatMultMapRefPtr m2;
		MatMultAddPair(CFLOBDDTopNodeMatMultMapRefPtr p1, CFLOBDDTopNodeMatMultMapRefPtr p2)
		{
			m1 = p1;
			m2 = p2;
		}

		struct MatMultAddPairHash {
			size_t operator()(const MatMultAddPair& p) const
			{
				CFLOBDDTopNodeMatMultMapRefPtr t1 = p.m1;
				CFLOBDDTopNodeMatMultMapRefPtr t2 = p.m2;
				auto hash1 = t1->Hash(997);
				auto hash2 = t2->Hash(997);
				return 117 * hash1 + hash2;
			}
		};

		bool operator==(const MatMultAddPair& p) const
		{
			CFLOBDDTopNodeMatMultMapRefPtr m11 = m1;
			CFLOBDDTopNodeMatMultMapRefPtr m12 = m2;
			CFLOBDDTopNodeMatMultMapRefPtr m21 = p.m1;
			CFLOBDDTopNodeMatMultMapRefPtr m22 = p.m2;
			return (m11 == m21) && (m12 == m22);
		}
	};

	class ZeroValNodeInfo {
	public:
		CFLOBDDNodeHandle m;
		unsigned int index;
		ZeroValNodeInfo(CFLOBDDNodeHandle p, unsigned int i) : m(p), index(i)
		{
		}

		struct ZeroValNodeInfoHash {
			size_t operator()(const ZeroValNodeInfo& p) const
			{
				CFLOBDDNodeHandle t1 = p.m;
				auto hash1 = t1.Hash(997);
				return 117 * (hash1 + 1) + 997 * p.index;
			}
		};

		bool operator==(const ZeroValNodeInfo& p) const
		{
			CFLOBDDNodeHandle m11 = m;
			CFLOBDDNodeHandle m21 = p.m;
			return (m11 == m21) && (index == p.index);
		}
	};


} // namespace CFL_OBDD

#endif
