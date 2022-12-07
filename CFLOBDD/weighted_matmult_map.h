#ifndef WEIGHTED_MATMULT_MAP_GUARD
#define WEIGHTED_MATMULT_MAP_GUARD

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <map>

#include "hashset.h"

class WeightedMatMultMapHandle;
class WeightedMatMultMapBody;


typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;

typedef std::pair<long int, long int> INT_PAIR;


//***************************************************************
// WeightedMatMultMapHandle
//***************************************************************

class WeightedMatMultMapHandle {
public:
	WeightedMatMultMapHandle();                               // Default constructor
	~WeightedMatMultMapHandle();                              // Destructor
	WeightedMatMultMapHandle(const WeightedMatMultMapHandle &r);             // Copy constructor
	WeightedMatMultMapHandle(const int size);             // constructor with size for vector
	WeightedMatMultMapHandle& operator= (const WeightedMatMultMapHandle &r); // Overloaded assignment
	bool operator!= (const WeightedMatMultMapHandle &r) const;      // Overloaded !=
	bool operator== (const WeightedMatMultMapHandle &r) const;      // Overloaded ==
	BIG_FLOAT& operator[](INT_PAIR& p);                        // Overloaded []
	unsigned int Hash(unsigned int modsize);
	void Add(const INT_PAIR& p, BIG_FLOAT& v);
	void ForceAdd(const INT_PAIR& p, BIG_FLOAT& v);
	bool Member(INT_PAIR& p);
	BIG_FLOAT Lookup(INT_PAIR& x);
	std::string ToString();
	unsigned int Size();
	void Canonicalize();
    size_t getHashCheck();
	WeightedMatMultMapBody *mapContents;
	static Hashset<WeightedMatMultMapBody> *canonicalWeightedMatMultMapBodySet;
	std::ostream& print(std::ostream & out = std::cout) const;
	WeightedMatMultMapHandle operator+ (const WeightedMatMultMapHandle&) const; // map concatenation with summation as merge operation
};

std::ostream& operator<< (std::ostream & out, const WeightedMatMultMapHandle &r);

extern WeightedMatMultMapHandle operator* (const BIG_FLOAT&, const WeightedMatMultMapHandle&);
extern WeightedMatMultMapHandle operator* (const WeightedMatMultMapHandle&, const BIG_FLOAT&);
extern std::size_t hash_value(const WeightedMatMultMapHandle& val);

//***************************************************************
// WeightedMatMultMapBody
//***************************************************************

class WeightedMatMultMapBody {

	friend void WeightedMatMultMapHandle::Canonicalize();
	friend unsigned int WeightedMatMultMapHandle::Hash(unsigned int modsize);

public:
	WeightedMatMultMapBody();    // Constructor
	void IncrRef();
	void DecrRef();
	unsigned int Hash(unsigned int modsize);
	void setHashCheck();
	unsigned int refCount;         // reference-count value
	std::map<INT_PAIR, BIG_FLOAT> map;
	bool operator==(const WeightedMatMultMapBody &o) const;
	BIG_FLOAT& operator[](INT_PAIR& i);                        // Overloaded []
	long int hashCheck;
	bool contains_zero_val;
	bool isCanonical;              // Is this WeightedMatMultMapBody in *canonicalWeightedMatMultMapBodySet?

};

#endif
