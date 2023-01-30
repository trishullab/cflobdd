
#include <iostream>
#include <fstream>
#include <cassert>
#include "weighted_matmult_map.h"
#include "fourier_semiring.h"
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

//***************************************************************
// WeightedMatMultMapBody
//***************************************************************

// Constructor
template <typename T>
WeightedMatMultMapBody<T>::WeightedMatMultMapBody()
	: refCount(0), isCanonical(false), contains_zero_val(false)
{
	hashCheck = NULL;
}

template <typename T>
void WeightedMatMultMapBody<T>::IncrRef()
{
	refCount++;    // Warning: Saturation not checked
}

template <typename T>
void WeightedMatMultMapBody<T>::DecrRef()
{
	if (--refCount == 0) {    // Warning: Saturation not checked
		if (isCanonical) {
			WeightedMatMultMapHandle<T>::canonicalWeightedMatMultMapBodySet->DeleteEq(this);
		}
		delete this;
	}
}

template <typename T>
unsigned int WeightedMatMultMapBody<T>::Hash(unsigned int modsize)
{
	/*if (modsize == HASHSETBASE)
		return hashCheck;*/
	unsigned int hvalue = 0;
	boost::hash<T> boost_hash;
	for (auto &i : map)
	{
		hvalue = (997 * hvalue + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second))) % modsize;
	}
	return hvalue;
}

template <typename T>
void WeightedMatMultMapBody<T>::setHashCheck()
{
	long int hvalue = 0;
	boost::hash<T> boost_hash;
	for (auto &i : map) {
		hvalue = (117 * (hvalue + 1) + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second)));
	}
	hashCheck = hvalue;
}

template <>
void WeightedMatMultMapBody<BIG_FLOAT>::setHashCheck()
{
	long int hvalue = 0;
	boost::hash<BIG_FLOAT> boost_hash;
	for (auto &i : map) {
		hvalue = (117 * (hvalue + 1) + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second)));
	}
	hashCheck = hvalue;
}

template <>
void WeightedMatMultMapBody<BIG_COMPLEX_FLOAT>::setHashCheck()
{
	long int hvalue = 0;
	boost::hash<BIG_COMPLEX_FLOAT> boost_hash;
	for (auto &i : map) {
		hvalue = (117 * (hvalue + 1) + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second)));
	}
	hashCheck = hvalue;
}

template <>
void WeightedMatMultMapBody<fourierSemiring>::setHashCheck()
{
	long int hvalue = 0;
	for (auto &i : map) {
        if (!i.second.isComplexValueSet)
		    hvalue = (117 * (hvalue + 1) + (int)(i.first.first + 97 * i.first.second + 97 * 97 * (i.second.GetVal() + 17 * i.second.GetRingSize())));
        else
        {
            boost::hash<BIG_COMPLEX> boost_hash;
            hvalue = (117 * (hvalue + 1) + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second.complex_value))); 
        }
	}
	hashCheck = hvalue;
}

template <typename T>
bool WeightedMatMultMapBody<T>::operator==(const WeightedMatMultMapBody<T> &o) const
{
	if (hashCheck != o.hashCheck) {
		return false;
	}
	else if (map.size() != o.map.size()) {
		return false;
	}
	else {
		for (auto m1_it = map.begin(), m2_it = o.map.begin(); m1_it != map.end() && m2_it != o.map.end(); m1_it++, m2_it++){
			if ((m1_it->first.first != m2_it->first.first) || (m1_it->first.second != m2_it->first.second) || 
                ((m1_it->second != m2_it->second)))
				return false;
		}
	}
	return true;
}

// Overloaded []
template <typename T>
T& WeightedMatMultMapBody<T>::operator[](INT_PAIR& p)
{
	return map[p];
}

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedMatMultMapBody<T> &r)
{
	out << "{WMMMB: ";
	for (auto &i : r.map) {
		out << "(" << i.first.first << "," << i.first.second << "," << i.second << ");";
	}
	out << " WMMMB}";
	return out;
}

//***************************************************************
// WeightedMatMultMapHandle
//***************************************************************

// Initializations of static members ---------------------------------
template <typename T>
Hashset<WeightedMatMultMapBody<T>> *WeightedMatMultMapHandle<T>::canonicalWeightedMatMultMapBodySet = new Hashset<WeightedMatMultMapBody<T>>(HASHSET_NUM_BUCKETS);

// Default constructor
template <typename T>
WeightedMatMultMapHandle<T>::WeightedMatMultMapHandle()
	: mapContents(new WeightedMatMultMapBody<T>)
{
	mapContents->IncrRef();
}

// Destructor
template <typename T>
WeightedMatMultMapHandle<T>::~WeightedMatMultMapHandle()
{
	mapContents->DecrRef();
}

// Copy constructor
template <typename T>
WeightedMatMultMapHandle<T>::WeightedMatMultMapHandle(const WeightedMatMultMapHandle<T> &r)
	: mapContents(r.mapContents)
{
	mapContents->IncrRef();
}

// Overloaded assignment
template <typename T>
WeightedMatMultMapHandle<T>& WeightedMatMultMapHandle<T>::operator= (const WeightedMatMultMapHandle<T> &r)
{
	if (this != &r)      // don't assign to self!
	{
		WeightedMatMultMapBody<T> *temp = mapContents;
		mapContents = r.mapContents;
		mapContents->IncrRef();
		temp->DecrRef();
	}
	return *this;
}

// Overloaded !=
template <typename T>
bool WeightedMatMultMapHandle<T>::operator!=(const WeightedMatMultMapHandle<T> &r) const
{
	return (mapContents != r.mapContents);
}

// Overloaded ==
template <typename T>
bool WeightedMatMultMapHandle<T>::operator==(const WeightedMatMultMapHandle<T> &r) const
{

	return (mapContents == r.mapContents);
}

// Overloaded []
template <typename T>
T& WeightedMatMultMapHandle<T>::operator[](INT_PAIR& i)
{
	return (*(this->mapContents))[i];
}

// print
template <typename T>
std::ostream& WeightedMatMultMapHandle<T>::print(std::ostream & out) const
{
	out << *mapContents;
	return out;
}

template <typename T>
std::ostream& operator<< (std::ostream & out, const WeightedMatMultMapHandle<T> &r)
{
	r.print(out);
	return(out);
}

template <typename T>
unsigned int WeightedMatMultMapHandle<T>::Hash(unsigned int modsize)
{
	if (!(mapContents->isCanonical)) {
		std::cout << "Hash of a non-canonical LinearMapHandle occurred" << std::endl;
		abort();
		this->Canonicalize();
	}
	assert(mapContents->isCanonical);
	return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

template <typename T>
void WeightedMatMultMapHandle<T>::Add(const INT_PAIR& p, T& v)
{
	assert(mapContents->refCount <= 1);
	auto it = mapContents->map.find(p);
    if (it == mapContents->map.end()){
        mapContents->map.emplace(p, v);
    }
    else{
        it->second += v;
        if (it->second == 0)
            mapContents->map.erase(it);
    }
}

template <>
void WeightedMatMultMapHandle<fourierSemiring>::Add(const INT_PAIR& p, fourierSemiring& v)
{
	assert(mapContents->refCount <= 1);
	auto it = mapContents->map.find(p);
    if (it == mapContents->map.end()){
        mapContents->map.emplace(p, v);
    }
    else{
        it->second = it->second + v;
        if (!it->second.isComplexValueSet){
            if (it->second == fourierSemiring(0, 1))
                mapContents->map.erase(it);
        }
        else
        {
            if (it->second.complex_value == BIG_COMPLEX(0))
                mapContents->map.erase(it);
        }
    }
}

template <typename T>
void WeightedMatMultMapHandle<T>::ForceAdd(const INT_PAIR& p, T& v)
{
	assert(mapContents->refCount <= 1);
	auto it = mapContents->map.find(p);
	if (it == mapContents->map.end()){
		mapContents->map.emplace(p, v);
		if (p.first == -1 || p.second == -1)
			mapContents->contains_zero_val = true;
	}
	else{
		if (p.first != -1 && p.second != -1){
			it->second = it->second + v;
        }
		else{
			assert(it->first.first == -1 && it->first.second == -1);
		}
	}
}

template <typename T>
bool WeightedMatMultMapHandle<T>::Member(INT_PAIR& y)
{
	return (mapContents->map.find(y) != mapContents->map.end());
}

template <typename T>
T WeightedMatMultMapHandle<T>::Lookup(INT_PAIR& p)
{
	return mapContents->map[p];
}

template <typename T>
unsigned int WeightedMatMultMapHandle<T>::Size(){
	return (unsigned int)mapContents->map.size();
}

template <typename T>
std::string WeightedMatMultMapHandle<T>::ToString(){
	std::string map_string;
	for (auto &i : mapContents->map){
		//map_string += std::to_string(i.first.first) + "," + std::to_string(i.first.second);
		//map_string += "," + std::to_string(i.second);
		map_string += boost::lexical_cast<std::string>(i.first.first) + "," + boost::lexical_cast<std::string>(i.first.second);
		map_string += "," + boost::lexical_cast<std::string>(i.second);
		map_string += ";";
	}
	return map_string;
}


template <typename T>
void WeightedMatMultMapHandle<T>::Canonicalize()
{

    // TODO: maybe change this to a function
    // for (auto it = mapContents->map.begin(); it != mapContents->map.end();)
    // {
    //     if (it->second == 0)
    //         mapContents->map.erase(it++);
    //     else
    //         ++it;
    // }

    if (mapContents->map.empty())
        mapContents->map.insert(std::make_pair(std::make_pair(-1, -1), 0));

	WeightedMatMultMapBody<T> *answerContents;
	mapContents->setHashCheck();

	if (!mapContents->isCanonical) {
		unsigned int hash = canonicalWeightedMatMultMapBodySet->GetHash(mapContents);
		answerContents = canonicalWeightedMatMultMapBodySet->Lookup(mapContents, hash);
		if (answerContents == NULL) {
			canonicalWeightedMatMultMapBodySet->Insert(mapContents, hash);
			mapContents->isCanonical = true;
			//mapContents->setHashCheck();
		}
		else {
			answerContents->IncrRef();
			mapContents->DecrRef();
			mapContents = answerContents;
		}
	}
}

template <>
void WeightedMatMultMapHandle<fourierSemiring>::Canonicalize()
{

    if (mapContents->map.empty())
        mapContents->map.insert(std::make_pair(std::make_pair(-1, -1), fourierSemiring(0, 1)));

	WeightedMatMultMapBody<fourierSemiring> *answerContents;
	mapContents->setHashCheck();

	if (!mapContents->isCanonical) {
		unsigned int hash = canonicalWeightedMatMultMapBodySet->GetHash(mapContents);
		answerContents = canonicalWeightedMatMultMapBodySet->Lookup(mapContents, hash);
		if (answerContents == NULL) {
			canonicalWeightedMatMultMapBodySet->Insert(mapContents, hash);
			mapContents->isCanonical = true;
		}
		else {
			answerContents->IncrRef();
			mapContents->DecrRef();
			mapContents = answerContents;
		}
	}
}

template <typename T>
size_t WeightedMatMultMapHandle<T>::getHashCheck()
{
    return mapContents->hashCheck;
}

// Linear operators on WeightedMatMultMapHandles ------------------------------------------------

// Binary addition
template <typename T>
WeightedMatMultMapHandle<T> WeightedMatMultMapHandle<T>::operator+ (const WeightedMatMultMapHandle<T>& mapHandle) const
{
	WeightedMatMultMapHandle<T> ans;
	for (auto &i : mapContents->map){
		ans.Add(i.first, i.second);
	}
	for (auto &i : mapHandle.mapContents->map){
		ans.Add(i.first, i.second);
	}
	if (ans.Size() == 0){
		T one(0);
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

template <>
WeightedMatMultMapHandle<fourierSemiring> WeightedMatMultMapHandle<fourierSemiring>::operator+ (const WeightedMatMultMapHandle<fourierSemiring>& mapHandle) const
{
	WeightedMatMultMapHandle<fourierSemiring> ans;
	for (auto &i : mapContents->map){
		ans.Add(i.first, i.second);
	}
	for (auto &i : mapHandle.mapContents->map){
		ans.Add(i.first, i.second);
	}
	if (ans.Size() == 0){
		fourierSemiring one(0, 1);
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

// Left scalar-multiplication
template <typename T>
WeightedMatMultMapHandle<T> operator* (const T& factor, const WeightedMatMultMapHandle<T>& mapHandle)
{
	if (factor == 1)
		return mapHandle;
	WeightedMatMultMapHandle<T> ans;
	for (auto &i : mapHandle.mapContents->map){
		T v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		T one(0);
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

template <>
WeightedMatMultMapHandle<fourierSemiring> operator* (const fourierSemiring& factor, const WeightedMatMultMapHandle<fourierSemiring>& mapHandle)
{
	if (!factor.isComplexValueSet && factor == fourierSemiring(1, 1))
		return mapHandle;
    if (factor.isComplexValueSet && factor.complex_value == BIG_COMPLEX(1))
        return mapHandle;
	WeightedMatMultMapHandle<fourierSemiring> ans;
	for (auto &i : mapHandle.mapContents->map){
		fourierSemiring v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		fourierSemiring one(0, 1);
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

// Right scalar-multiplication
template <typename T>
WeightedMatMultMapHandle<T> operator* (const WeightedMatMultMapHandle<T>& mapHandle, const T& factor)
{
	if (factor == 1)
		return mapHandle;
	WeightedMatMultMapHandle<T> ans;
	for (auto &i : mapHandle.mapContents->map){
		T v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		T one(0);
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

template <>
WeightedMatMultMapHandle<fourierSemiring> operator* (const WeightedMatMultMapHandle<fourierSemiring>& mapHandle, const fourierSemiring& factor)
{
	if (!factor.isComplexValueSet && factor == fourierSemiring(1, 1))
		return mapHandle;
    if (factor.isComplexValueSet && factor.complex_value == BIG_COMPLEX(1))
        return mapHandle;
	WeightedMatMultMapHandle<fourierSemiring> ans;
	for (auto &i : mapHandle.mapContents->map){
		fourierSemiring v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		fourierSemiring zero(0, 1);
		ans.ForceAdd(std::make_pair(-1, -1), zero);
	}
	ans.Canonicalize();
	return ans;
}

template <typename T>
std::size_t hash_value(const WeightedMatMultMapHandle<T>& val)
{
	return val.mapContents->hashCheck;
}

template class WeightedMatMultMapHandle<BIG_FLOAT>;
template WeightedMatMultMapHandle<BIG_FLOAT> operator*<BIG_FLOAT>(const BIG_FLOAT&, const WeightedMatMultMapHandle<BIG_FLOAT>&);
template WeightedMatMultMapHandle<BIG_FLOAT> operator*<BIG_FLOAT>(const WeightedMatMultMapHandle<BIG_FLOAT>&, const BIG_FLOAT&);

template class WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>;
template WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> operator*<BIG_COMPLEX_FLOAT>(const BIG_COMPLEX_FLOAT&, const WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>&);
template WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT> operator*<BIG_COMPLEX_FLOAT>(const WeightedMatMultMapHandle<BIG_COMPLEX_FLOAT>&, const BIG_COMPLEX_FLOAT&);


template class WeightedMatMultMapHandle<fourierSemiring>;
template WeightedMatMultMapHandle<fourierSemiring> operator*<fourierSemiring>(const fourierSemiring&, const WeightedMatMultMapHandle<fourierSemiring>&);
template WeightedMatMultMapHandle<fourierSemiring> operator*<fourierSemiring>(const WeightedMatMultMapHandle<fourierSemiring>&, const fourierSemiring&);