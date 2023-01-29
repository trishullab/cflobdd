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
// #include "return_map_T.h"
#include "intpair.h"
#include "cflobdd_node.h"
#include "matmult_map.h"
#include "fourier_semiring.h"
#include "weighted_matmult_map.h"


#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
//typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<1000> > BIG_FLOAT;
typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
//typedef mp::number<mp::cpp_dec_float<200> > BIG_FLOAT;

// Instantiation and specialization of class ReturnMapHandle<LinearMapHandle> ----------

// template<>
// unsigned int ReturnMapBody<LinearMapHandle>::Hash(unsigned int modsize)
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
unsigned int ReturnMapBody<MatMultMapHandle>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash(modsize)) % modsize;
	}
	return hvalue;
}

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
unsigned int ReturnMapBody<WeightedMatMultMapHandle<BIG_FLOAT>>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash(modsize)) % modsize;
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
unsigned int ReturnMapBody<WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash(modsize)) % modsize;
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
unsigned int ReturnMapBody<WeightedMatMultMapHandle<fourierSemiring>>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + mapArray[i].Hash(modsize)) % modsize;
	}
	return hvalue;
}


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

template<>
unsigned int ReturnMapBody<int>::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;

  for (unsigned i = 0; i < mapArray.size(); i++)
  {
	  hvalue = (997* hvalue + mapArray[i]) % modsize;
  }
  return hvalue;
}

template<>
void ReturnMapBody<int>::setHashCheck()
{
  unsigned int hvalue = 0;

  for (unsigned i = 0; i < mapArray.size(); i++)
  {
	  hvalue = (117*(hvalue+1) + mapArray[i]);
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
	for (unsigned i = 0; i < size; i++)
	{
		answer.mapContents->mapArray.push_back(!mapContents->mapArray[i]);
	}
	answer.Canonicalize();
	return answer;
}


// Instantiation and specialization of class ReturnMapHandle<double> ---------------------

template<>
unsigned int ReturnMapBody<double>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;
	std::hash<double> double_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + double_hash(mapArray[i])) % modsize;
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
unsigned int ReturnMapBody<BIG_FLOAT>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;
	std::hash<BIG_FLOAT> big_float_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + big_float_hash(mapArray[i])) % modsize;
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
unsigned int ReturnMapBody<BIG_COMPLEX_FLOAT>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;
	std::hash<BIG_COMPLEX_FLOAT> big_float_hash;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + big_float_hash(mapArray[i])) % modsize;
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
// unsigned int ReturnMapBody<GeneralMapHandle>::Hash(unsigned int modsize)
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
unsigned int ReturnMapBody<intpair>::Hash(unsigned int modsize)
{
  unsigned int hvalue = 0;

  for (unsigned i = 0; i < mapArray.size(); i++)
  {
	  hvalue = (997 * hvalue + mapArray[i].First() + mapArray[i].Second()) % modsize;
  } 
  return hvalue;
}

template<>
void ReturnMapBody<intpair>::setHashCheck()
{
  unsigned int hvalue = 0;

  for (unsigned i = 0; i < mapArray.size(); i++)
  {
	  hvalue = (117*(hvalue+1) + mapArray[i].First() + mapArray[i].Second());
  }
  hashCheck = hvalue;
}

template<>
ReturnMapHandle<intpair> ReturnMapHandle<intpair>::Complement()
{
	ReturnMapHandle<intpair> answer;
	unsigned size = mapContents->mapArray.size();
	for (unsigned i = 0; i < size; i++)
	{
		int c0 = !mapContents->mapArray[i].First();
		int c1 = !mapContents->mapArray[i].Second();
		answer.mapContents->mapArray.push_back(intpair(c0, c1)); 
	}
	answer.Canonicalize();
	return answer;
}

// Instantiation and specialization of class ReturnMapHandle<fourierSemiring> ---------------------

template<>
unsigned int ReturnMapBody<fourierSemiring>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;
	boost::hash<BIG_INT> boost_hash;
	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		if (!mapArray[i].isComplexValueSet)
			hvalue = (997 * hvalue + 117 * boost_hash(mapArray[i].GetVal()) + boost_hash(mapArray[i].GetRingSize())) % modsize;
		else
		{
			boost::hash<BIG_COMPLEX> h;
			hvalue = (997 * hvalue + 97 * h(mapArray[i].complex_value)) % modsize;	
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

/*
std::size_t hash_value(std::complex<double> c)
{
	std::hash<double> double_hash;
	return (997 * double_hash(real(c)) + double_hash(imag(c)));
}
*/

template<>
unsigned int ReturnMapBody<std::complex<double>>::Hash(unsigned int modsize)
{
	unsigned int hvalue = 0;

	for (unsigned i = 0; i < mapArray.size(); i++)
	{
		hvalue = (997 * hvalue + hash_complex_double(mapArray[i])) % modsize;
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



