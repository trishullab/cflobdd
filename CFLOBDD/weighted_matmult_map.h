#ifndef WEIGHTED_MATMULT_MAP_GUARD
#define WEIGHTED_MATMULT_MAP_GUARD

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <map>

#include "hashset.h"

template <typename>
class WeightedMatMultMapHandle;
template <typename>
class WeightedMatMultMapBody;


typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
// typedef double BIG_FLOAT;

typedef std::pair<long int, long int> INT_PAIR;


//***************************************************************
// WeightedMatMultMapHandle
//***************************************************************

template <typename T>
class WeightedMatMultMapHandle {
public:
	WeightedMatMultMapHandle();                               // Default constructor
	~WeightedMatMultMapHandle();                              // Destructor
	WeightedMatMultMapHandle(const WeightedMatMultMapHandle<T> &r);             // Copy constructor
	WeightedMatMultMapHandle(const int size);             // constructor with size for vector
	WeightedMatMultMapHandle<T>& operator= (const WeightedMatMultMapHandle<T> &r); // Overloaded assignment
	bool operator!= (const WeightedMatMultMapHandle<T> &r) const;      // Overloaded !=
	bool operator== (const WeightedMatMultMapHandle<T> &r) const;      // Overloaded ==
	T& operator[](INT_PAIR& p);                        // Overloaded []
	unsigned int Hash(unsigned int modsize);
	void Add(const INT_PAIR& p, T& v);
	void ForceAdd(const INT_PAIR& p, T& v);
	bool Member(INT_PAIR& p);
	T Lookup(INT_PAIR& x);
	std::string ToString();
	unsigned int Size();
	void Canonicalize();
    size_t getHashCheck();
	WeightedMatMultMapBody<T> *mapContents;
	static Hashset<WeightedMatMultMapBody<T>> *canonicalWeightedMatMultMapBodySet;
	std::ostream& print(std::ostream & out = std::cout) const;
	WeightedMatMultMapHandle<T> operator+ (const WeightedMatMultMapHandle<T>&) const; // map concatenation with summation as merge operation
};

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedMatMultMapHandle<T> &r);

template <typename T>
extern WeightedMatMultMapHandle<T> operator* (const T&, const WeightedMatMultMapHandle<T>&);
template <typename T>
extern WeightedMatMultMapHandle<T> operator* (const WeightedMatMultMapHandle<T>&, const T&);
template <typename T>
extern std::size_t hash_value(const WeightedMatMultMapHandle<T>& val);

//***************************************************************
// WeightedMatMultMapBody
//***************************************************************

template <typename T>
class WeightedMatMultMapBody {

	friend void WeightedMatMultMapHandle<T>::Canonicalize();
	friend unsigned int WeightedMatMultMapHandle<T>::Hash(unsigned int modsize);

public:
	WeightedMatMultMapBody();    // Constructor
	void IncrRef();
	void DecrRef();
	unsigned int Hash(unsigned int modsize);
	void setHashCheck();
	unsigned int refCount;         // reference-count value
	std::map<INT_PAIR, T> map;
	bool operator==(const WeightedMatMultMapBody<T> &o) const;
	T& operator[](INT_PAIR& i);                        // Overloaded []
	long int hashCheck;
	bool contains_zero_val;
	bool isCanonical;              // Is this WeightedMatMultMapBody in *canonicalWeightedMatMultMapBodySet?

};

#endif
