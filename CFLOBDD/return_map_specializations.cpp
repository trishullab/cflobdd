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


#include <cassert>
#include <complex>
#include <functional>
#include <unordered_map>
#include <vector>
// #include "return_map_T.h"
#include "intpair.h"
#include "cflobdd_node.h"
#include "matmult_map.h"
#include "fourier_semiring.h"
#ifdef WCFLOBDD_SUPPORTED
#include "weighted_matmult_map.h"
#endif


#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
//typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000> > BIG_FLOAT;
typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
//typedef mp::number<mp::cpp_dec_float<200> > BIG_FLOAT;

// Instantiation and specialization of class ReturnMapHandle<LinearMapHandle> ----------

// template<>
// unsigned int ReturnMapBody<LinearMapHandle>::Hash()
// {
// 	unsigned int hvalue = 0;

// 	for (unsigned i = 0; i < mapArray.size(); i++)
// 	{
// 		hvalue = (997 * hvalue + mapArray[i].Hash(modsize)) % modsize;
// 	}
// 	return hvalue;
// }

template<>
void ReturnMapBody<MatMultMapHandle>::setHashCheck()
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		if (mapArray[i].mapContents->hashCheck == NULL) {
			mapArray[i].mapContents->setHashCheck();
		}
		hvalue = (117 * (hvalue + 1) + mapArray[i].mapContents->hashCheck);
	}
	hashCheck = hvalue;
}

template<>
size_t ReturnMapBody<MatMultMapHandle>::Hash()
{
	size_t hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash());
	}
	return hvalue;
}

#ifdef WCFLOBDD_SUPPORTED
template<>
void ReturnMapBody<WeightedMatMultMapHandle<BIG_FLOAT>>::setHashCheck()
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		if (mapArray[i].mapContents->hashCheck == NULL) {
			mapArray[i].mapContents->setHashCheck();
		}
		hvalue = (117 * (hvalue + 1) + mapArray[i].mapContents->hashCheck);
	}
	hashCheck = hvalue;
}

template<>
size_t ReturnMapBody<WeightedMatMultMapHandle<BIG_FLOAT>>::Hash()
{
	size_t hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash());
	}
	return hvalue;
}

template<>
void ReturnMapBody<WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>>::setHashCheck()
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		if (mapArray[i].mapContents->hashCheck == NULL) {
			mapArray[i].mapContents->setHashCheck();
		}
		hvalue = (117 * (hvalue + 1) + mapArray[i].mapContents->hashCheck);
	}
	hashCheck = hvalue;
}

template<>
size_t ReturnMapBody<WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>>::Hash()
{
	size_t hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash());
	}
	return hvalue;
}

template<>
void ReturnMapBody<WeightedMatMultMapHandle<fourierSemiring>>::setHashCheck()
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		if (mapArray[i].mapContents->hashCheck == NULL) {
			mapArray[i].mapContents->setHashCheck();
		}
		hvalue = (117 * (hvalue + 1) + mapArray[i].mapContents->hashCheck);
	}
	hashCheck = hvalue;
}

template<>
size_t ReturnMapBody<WeightedMatMultMapHandle<fourierSemiring>>::Hash()
{
	size_t hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash());
	}
	return hvalue;
}
#endif // WCFLOBDD_SUPPORTED


// template<>
// void ReturnMapBody<LinearMapHandle>::setHashCheck()
// {
// 	unsigned int hvalue = 0;

// 	for (unsigned i = 0; i < mapArray.size(); i++)
// 	{
// 		if (mapArray[i].mapContents->hashCheck == NULL) {
// 			mapArray[i].mapContents->setHashCheck();
// 		}
// 		hvalue = (117 * (hvalue + 1) + mapArray[i].mapContents->hashCheck);
// 	}
// 	hashCheck = hvalue;
// }

/*
template <>
std::ostream& operator<< (std::ostream & out, const ReturnMapBody<LinearMapHandle> &r)
{
	out << "{RMB<LinearMapHandle>: ";
	size_t last = r.mapArray.size() - 1;
	for (size_t i = 0; i <= last; ++i) {
		out << r.mapArray[i];
		if (i != last)
			out << ", ";
	}
	out << " RMB<LinearMapHandle>}";
	return out;
}
*/


// template<>
// ReturnMapHandle<LinearMapHandle> ReturnMapHandle<LinearMapHandle>::Complement()
// {
// 	assert(false);
// 	return *this;   // Should never be executed; included to supporess VS 2013 error report
// }

// template <>
// ReturnMapHandle<LinearMapHandle> ReturnMapHandle<LinearMapHandle>::Compose(ReductionMapHandle redMapHandle)
// {
// 	assert(false);
// 	return *this;   // Should never be executed; included to suppress VS 2013 error report
// }

// template class ReturnMapHandle<LinearMapHandle>;

// Instantiation and specialization of class ReturnMapHandle<int> ---------------------

// Murmur3 finalizer — ensures full avalanche (each output bit depends on all input bits)
static inline size_t fmix64(size_t h) {
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

template<>
size_t ReturnMapBody<int>::Hash()
{
  return fmix64(hashCheck);
}

template<>
void ReturnMapBody<int>::setHashCheck()
{
  unsigned int hvalue = 0;
  const auto* data = mapArray.data();
  unsigned int sz = mapArray.size();
  for (unsigned i = 0; i < sz; i++)
  {
	  hvalue = (131*(hvalue+1) + data[i]);
  }
  hashCheck = hvalue;
}

/*
template <>
std::ostream& operator<< (std::ostream & out, const ReturnMapBody<int> &r)
{
	out << "{RMB<int>: ";
	size_t last = r.mapArray.size() - 1;
	for (size_t i = 0; i <= last; ++i) {
		out << r.mapArray[i];
		if (i != last)
			out << ", ";
	}
	out << " RMB<int>}";
	return out;
}
*/

template<>
CFL_OBDD::CFLOBDDReturnMapHandle CFL_OBDD::CFLOBDDReturnMapHandle::Complement()
// template<> ReturnMapHandle<int> ReturnMapHandle<int>::Complement()
{
	ReturnMapHandle<int> answer;
	unsigned size = mapContents->mapArray.size();
	const auto* srcData = mapContents->mapArray.data();
	for (unsigned i = 0; i < size; i++)
	{
		answer.mapContents->mapArray.push_back(!srcData[i]);
	}
	answer.Canonicalize();
	return answer;
}


// Instantiation and specialization of class ReturnMapHandle<double> ---------------------

template<>
size_t ReturnMapBody<double>::Hash()
{
	size_t hvalue = 0;
	std::hash<double> double_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + double_hash(mapArray[i]));
	}
	return hvalue;
}

template<>
void ReturnMapBody<double>::setHashCheck()
{
	unsigned int hvalue = 0;
	std::hash<double> double_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (117 * (hvalue + 1) + double_hash(mapArray[i]));
	}
	hashCheck = hvalue;
}


template<>
ReturnMapHandle<double> ReturnMapHandle<double>::Complement()
{
	assert(false);
	return *this;   // Should never be executed; included to supporess VS 2013 error report
}


// Instantiation and specialization of class ReturnMapHandle<BIG_FLOAT> ---------------------
// double_hash needs to be changed

template<>
size_t ReturnMapBody<BIG_FLOAT>::Hash()
{
	size_t hvalue = 0;
	std::hash<BIG_FLOAT> big_float_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + big_float_hash(mapArray[i]));
	}
	return hvalue;
}

template<>
void ReturnMapBody<BIG_FLOAT>::setHashCheck()
{
	unsigned int hvalue = 0;
	std::hash<BIG_FLOAT> big_float_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (117 * (hvalue + 1) + big_float_hash(mapArray[i]));
	}
	hashCheck = hvalue;
}


template<>
ReturnMapHandle<BIG_FLOAT> ReturnMapHandle<BIG_FLOAT>::Complement()
{
	assert(false);
	return *this;   // Should never be executed; included to supporess VS 2013 error report
}


// Instantiation and specialization of class ReturnMapHandle<BIG_COMPLEX_FLOAT> ---------------------
// double_hash needs to be changed

template<>
size_t ReturnMapBody<BIG_COMPLEX_FLOAT>::Hash()
{
	size_t hvalue = 0;
	std::hash<BIG_COMPLEX_FLOAT> big_float_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + big_float_hash(mapArray[i]));
	}
	return hvalue;
}

template<>
void ReturnMapBody<BIG_COMPLEX_FLOAT>::setHashCheck()
{
	unsigned int hvalue = 0;
	std::hash<BIG_COMPLEX_FLOAT> big_float_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (117 * (hvalue + 1) + big_float_hash(mapArray[i]));
	}
	hashCheck = hvalue;
}


template<>
ReturnMapHandle<BIG_COMPLEX_FLOAT> ReturnMapHandle<BIG_COMPLEX_FLOAT>::Complement()
{
	assert(false);
	return *this;   // Should never be executed; included to supporess VS 2013 error report
}

// Instantiation and specialization of class ReturnMapHandle<GeneralMapHandle> ---------------------
// double_hash needs to be changed

// template<>
// unsigned int ReturnMapBody<GeneralMapHandle>::Hash()
// {
// 	unsigned int hvalue = 0;
// 	std::hash<double> double_hash;

// 	for (unsigned i = 0; i < mapArray.size(); i++)
// 	{
// 		hvalue = (997 * hvalue + 97*mapArray[i].Hash(modsize)) % modsize;
// 	}
// 	return hvalue;
// }


// template<>
// void ReturnMapBody<GeneralMapHandle>::setHashCheck()
// {
// 	unsigned int hvalue = 0;
// 	std::hash<double> double_hash;

// 	for (unsigned i = 0; i < mapArray.size(); i++)
// 	{
// 		hvalue = (117 * (hvalue + 1) + mapArray[i].Hash(997));
// 	}
// 	hashCheck = hvalue;
// }


// template<>
// ReturnMapHandle<GeneralMapHandle> ReturnMapHandle<GeneralMapHandle>::Complement()
// {
// 	assert(false);
// 	return *this;   // Should never be executed; included to supporess VS 2013 error report
// }

// Instantiation and specialization of class ReturnMapHandle<intPair> ---------------------

template<>
size_t ReturnMapBody<intpair>::Hash()
{
  return fmix64(hashCheck);
}

template<>
void ReturnMapBody<intpair>::setHashCheck()
{
  unsigned int hvalue = 0;
  const auto* data = mapArray.data();
  unsigned int sz = mapArray.size();
  for (unsigned i = 0; i < sz; i++)
  {
	  hvalue = (131*(hvalue+1) + (int)(97 * data[i].First()) + data[i].Second());
  }
  hashCheck = hvalue;
}

template<>
ReturnMapHandle<intpair> ReturnMapHandle<intpair>::Complement()
{
	ReturnMapHandle<intpair> answer;
	unsigned size = mapContents->mapArray.size();
	const auto* srcData = mapContents->mapArray.data();
	for (unsigned i = 0; i < size; i++)
	{
		int c0 = !srcData[i].First();
		int c1 = !srcData[i].Second();
		answer.mapContents->mapArray.push_back(intpair(c0, c1)); 
	}
	answer.Canonicalize();
	return answer;
}

// Instantiation and specialization of class ReturnMapHandle<fourierSemiring> ---------------------

template<>
size_t ReturnMapBody<fourierSemiring>::Hash()
{
	size_t hvalue = 0;
	boost::hash<BIG_INT> boost_hash;
	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		if (!mapArray[i].isComplexValueSet)
			hvalue = (997 * hvalue + 117 * boost_hash(mapArray[i].GetVal()) + boost_hash(mapArray[i].GetRingSize()));
		else
		{
			boost::hash<BIG_COMPLEX> h;
			hvalue = (997 * hvalue + 97 * h(mapArray[i].complex_value));
		}
	}
	return hvalue;
}

template<>
void ReturnMapBody<fourierSemiring>::setHashCheck()
{
	unsigned int hvalue = 0;
	boost::hash<BIG_INT> boost_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		if (!mapArray[i].isComplexValueSet)
			hvalue = (117 * (hvalue + 1) + 97 * boost_hash(mapArray[i].GetVal()) + boost_hash(mapArray[i].GetRingSize()));
		else{
			boost::hash<BIG_COMPLEX> h;
			hvalue = (117 * (hvalue + 1) + 97 * h(mapArray[i].complex_value));
		}
	}
	hashCheck = hvalue;
}

template<>
ReturnMapHandle<fourierSemiring> ReturnMapHandle<fourierSemiring>::Complement()
{
	assert(false);
	return *this;
}

// Instantiation and specialization of class ReturnMapHandle<std::complex<double>> ---------------------

size_t hash_complex_double(std::complex<double> c)
{
	std::hash<double> double_hash;
	return (997 * double_hash(real(c)) + double_hash(imag(c)));
}

template<>
size_t ReturnMapBody<std::complex<double>>::Hash()
{
	size_t hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + hash_complex_double(mapArray[i]));
	}
	return hvalue;
}

template<>
void ReturnMapBody<std::complex<double>>::setHashCheck()
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (117 * (hvalue + 1) + hash_complex_double(mapArray[i]));
	}
	hashCheck = hvalue;
}

/*
template <>
std::ostream& operator<< (std::ostream & out, const ReturnMapBody<std::complex<double>> &r)
{
    out << "{CMB<std::complex<double>>: ";
    size_t last = r.mapArray.size() - 1;
    for (size_t i = 0; i <= last; ++i) {
        out << r.mapArray[i];
        if (i != last)
            out << ", ";
    }
    out << " CMB<std::complex<double>>}";
    return out;
}
*/

template<>
ReturnMapHandle<std::complex<double>> ReturnMapHandle<std::complex<double>>::Complement()
{
	assert(false);
	return *this;   // Should never be executed; included to supporess VS 2013 error report
}

// ============================================================================
// MakeIdentityReturnMap -- cached construction of identity return maps
// ============================================================================
//
// Returns a canonicalized CFLOBDDReturnMapHandle [0, 1, 2, ..., k-1].
// Previously computed maps are cached:
//   - For k <= IDENTITY_MAP_ARRAY_THRESHOLD (1024): flat array with O(1) lookup.
//   - For k > threshold: std::unordered_map for sparse large-size cases.
//
// The flat array uses a parallel bool vector to distinguish "not yet computed"
// from "computed" (since a default-constructed CFLOBDDReturnMapHandle is a
// valid but empty map, not a usable sentinel).

namespace CFL_OBDD {

static constexpr unsigned int IDENTITY_MAP_ARRAY_THRESHOLD = 1024;

CFLOBDDReturnMapHandle MakeIdentityReturnMap(unsigned int k)
{
    // Function-local statics: guaranteed to be initialized on first call
    // (avoids static-initialization-order fiasco across translation units).

    // Flat array cache for small sizes (indices 0..THRESHOLD)
    static std::vector<CFLOBDDReturnMapHandle> identityMapArray(IDENTITY_MAP_ARRAY_THRESHOLD + 1);
    static std::vector<bool> identityMapValid(IDENTITY_MAP_ARRAY_THRESHOLD + 1, false);
    // Overflow cache for large sizes (> THRESHOLD)
    static std::unordered_map<unsigned int, CFLOBDDReturnMapHandle> identityMapOverflow;

    if (k <= IDENTITY_MAP_ARRAY_THRESHOLD) {
        if (identityMapValid[k]) {
            return identityMapArray[k];
        }
        // Build, canonicalize, and cache the identity map of size k
        CFLOBDDReturnMapHandle idMap;
        for (unsigned int i = 0; i < k; i++) {
            idMap.AddToEnd(i);
        }
        idMap.Canonicalize();
        identityMapArray[k] = idMap;
        identityMapValid[k] = true;
        return idMap;
    }
    else {
        // Large-size path: use unordered_map
        auto it = identityMapOverflow.find(k);
        if (it != identityMapOverflow.end()) {
            return it->second;
        }
        CFLOBDDReturnMapHandle idMap;
        for (unsigned int i = 0; i < k; i++) {
            idMap.AddToEnd(i);
        }
        idMap.Canonicalize();
        identityMapOverflow.emplace(k, idMap);
        return idMap;
    }
}

} // namespace CFL_OBDD

// ComposeAndReduce: compose a return map with a reduction map, producing
// a reduced return map and an induced reduction map.
// Placed here so the compiler can inline ReturnMapHandle<int> and
// ReductionMapHandle operations in the hot loop.
CFL_OBDD::CFLOBDDReturnMapHandle ComposeAndReduce(CFL_OBDD::CFLOBDDReturnMapHandle& mapHandle, ReductionMapHandle& redMapHandle, ReductionMapHandle& inducedRedMapHandle)
{
	using CFL_OBDD::CFLOBDDReturnMapHandle;
	int c2, c3;
	int size = mapHandle.mapContents->mapArray.size();
	CFLOBDDReturnMapHandle answer;// (size);
	if (redMapHandle.mapContents->isIdentityMap){
		inducedRedMapHandle = redMapHandle;
		return mapHandle;
	}
	unsigned int redSize = redMapHandle.Size();
	// Static flat array reused across calls to avoid repeated allocation.
	// Only indices actually written are tracked in dirtyIndices for O(size) cleanup.
	static std::vector<int> flatMap;
	static std::vector<unsigned int> dirtyIndices;
	constexpr unsigned int COMPOSE_FLAT_THRESHOLD = 33554432; // 2^25 (> 5783^2 for 4096-bit)
	if (redSize <= COMPOSE_FLAT_THRESHOLD) {
		if (flatMap.size() < redSize) {
			flatMap.resize(redSize, -1);
		}
		dirtyIndices.clear();
		const auto* mapArrayData = mapHandle.mapContents->mapArray.data();
		const auto* redArrayData = redMapHandle.mapContents->mapArray.data();
		const unsigned int redArraySize = redMapHandle.mapContents->mapArray.size();
		for (int i = 0; i < size; i++)
		{
			c2 = mapArrayData[i];
			assert((unsigned int)c2 < redArraySize);
			c3 = redArrayData[c2];
			if (flatMap[c3] == -1){
				answer.AddToEnd(c3);
				flatMap[c3] = answer.Size() - 1;
				dirtyIndices.push_back(c3);
				inducedRedMapHandle.AddToEnd(answer.Size() - 1);
			}
			else{
				inducedRedMapHandle.AddToEnd(flatMap[c3]);
			}
		}
		// Reset only the indices we touched
		for (unsigned int idx : dirtyIndices) {
			flatMap[idx] = -1;
		}
	} else {
		// Fallback: unordered_map for very large reduction maps
		std::unordered_map<int, unsigned int> reductionMap(size);
		const auto* mapArrayData2 = mapHandle.mapContents->mapArray.data();
		const auto* redArrayData2 = redMapHandle.mapContents->mapArray.data();
		const unsigned int redArraySize2 = redMapHandle.mapContents->mapArray.size();
		for (int i = 0; i < size; i++)
		{
			c2 = mapArrayData2[i];
			assert((unsigned int)c2 < redArraySize2);
			c3 = redArrayData2[c2];
			if (reductionMap.find(c3) == reductionMap.end()){
				answer.AddToEnd(c3);
				reductionMap.emplace(c3, answer.Size() - 1);
				inducedRedMapHandle.AddToEnd(answer.Size() - 1);
			}
			else{
				inducedRedMapHandle.AddToEnd(reductionMap[c3]);
			}
		}
	}
	inducedRedMapHandle.Canonicalize();
	answer.Canonicalize();
	return answer;
}
