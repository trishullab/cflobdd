
#include <iostream>
#include <fstream>
#include <cassert>
#include "weighted_matmult_map.h"
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

//***************************************************************
// WeightedMatMultMapBody
//***************************************************************

// Constructor
WeightedMatMultMapBody::WeightedMatMultMapBody()
	: refCount(0), isCanonical(false), contains_zero_val(false)
{
	hashCheck = NULL;
}

void WeightedMatMultMapBody::IncrRef()
{
	refCount++;    // Warning: Saturation not checked
}

void WeightedMatMultMapBody::DecrRef()
{
	if (--refCount == 0) {    // Warning: Saturation not checked
		if (isCanonical) {
			WeightedMatMultMapHandle::canonicalWeightedMatMultMapBodySet->DeleteEq(this);
		}
		delete this;
	}
}

unsigned int WeightedMatMultMapBody::Hash(unsigned int modsize)
{
	/*if (modsize == HASHSETBASE)
		return hashCheck;*/
	unsigned int hvalue = 0;
	boost::hash<BIG_FLOAT> boost_hash;
	for (auto &i : map)
	{
		hvalue = (997 * hvalue + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second))) % modsize;
	}
	return hvalue;
}
void WeightedMatMultMapBody::setHashCheck()
{
	long int hvalue = 0;
	boost::hash<BIG_FLOAT> boost_hash;
	for (auto &i : map) {
		hvalue = (117 * (hvalue + 1) + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second)));
	}
	hashCheck = hvalue;
}

bool WeightedMatMultMapBody::operator==(const WeightedMatMultMapBody &o) const
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
BIG_FLOAT& WeightedMatMultMapBody::operator[](INT_PAIR& p)
{
	return map[p];
}

std::ostream& operator<< (std::ostream & out, const WeightedMatMultMapBody &r)
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
Hashset<WeightedMatMultMapBody> *WeightedMatMultMapHandle::canonicalWeightedMatMultMapBodySet = new Hashset<WeightedMatMultMapBody>(HASHSET_NUM_BUCKETS);

// Default constructor
WeightedMatMultMapHandle::WeightedMatMultMapHandle()
	: mapContents(new WeightedMatMultMapBody)
{
	mapContents->IncrRef();
}

// Destructor
WeightedMatMultMapHandle::~WeightedMatMultMapHandle()
{
	mapContents->DecrRef();
}

// Copy constructor
WeightedMatMultMapHandle::WeightedMatMultMapHandle(const WeightedMatMultMapHandle &r)
	: mapContents(r.mapContents)
{
	mapContents->IncrRef();
}

// Overloaded assignment
WeightedMatMultMapHandle& WeightedMatMultMapHandle::operator= (const WeightedMatMultMapHandle &r)
{
	if (this != &r)      // don't assign to self!
	{
		WeightedMatMultMapBody *temp = mapContents;
		mapContents = r.mapContents;
		mapContents->IncrRef();
		temp->DecrRef();
	}
	return *this;
}

// Overloaded !=
bool WeightedMatMultMapHandle::operator!=(const WeightedMatMultMapHandle &r) const
{
	return (mapContents != r.mapContents);
}

// Overloaded ==
bool WeightedMatMultMapHandle::operator==(const WeightedMatMultMapHandle &r) const
{

	return (mapContents == r.mapContents);
}

// Overloaded []
BIG_FLOAT& WeightedMatMultMapHandle::operator[](INT_PAIR& i)
{
	return (*(this->mapContents))[i];
}

// print
std::ostream& WeightedMatMultMapHandle::print(std::ostream & out) const
{
	out << *mapContents;
	return out;
}

std::ostream& operator<< (std::ostream & out, const WeightedMatMultMapHandle &r)
{
	r.print(out);
	return(out);
}

unsigned int WeightedMatMultMapHandle::Hash(unsigned int modsize)
{
	if (!(mapContents->isCanonical)) {
		std::cout << "Hash of a non-canonical LinearMapHandle occurred" << std::endl;
		abort();
		this->Canonicalize();
	}
	assert(mapContents->isCanonical);
	return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

void WeightedMatMultMapHandle::Add(const INT_PAIR& p, BIG_FLOAT& v)
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

void WeightedMatMultMapHandle::ForceAdd(const INT_PAIR& p, BIG_FLOAT& v)
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
			it->second += v;
        }
		else{
			assert(it->first.first == -1 && it->first.second == -1);
		}
	}
}

bool WeightedMatMultMapHandle::Member(INT_PAIR& y)
{
	return (mapContents->map.find(y) != mapContents->map.end());
}

BIG_FLOAT WeightedMatMultMapHandle::Lookup(INT_PAIR& p)
{
	return mapContents->map[p];
}

unsigned int WeightedMatMultMapHandle::Size(){
	return (unsigned int)mapContents->map.size();
}

std::string WeightedMatMultMapHandle::ToString(){
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


void WeightedMatMultMapHandle::Canonicalize()
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

	WeightedMatMultMapBody *answerContents;
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

size_t WeightedMatMultMapHandle::getHashCheck()
{
    return mapContents->hashCheck;
}

// Linear operators on WeightedMatMultMapHandles ------------------------------------------------

// Binary addition
WeightedMatMultMapHandle WeightedMatMultMapHandle::operator+ (const WeightedMatMultMapHandle& mapHandle) const
{
	WeightedMatMultMapHandle ans;
	for (auto &i : mapContents->map){
		ans.Add(i.first, i.second);
	}
	for (auto &i : mapHandle.mapContents->map){
		ans.Add(i.first, i.second);
	}
	if (ans.Size() == 0){
		BIG_FLOAT one = 1;
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

// Left scalar-multiplication
WeightedMatMultMapHandle operator* (const BIG_FLOAT& factor, const WeightedMatMultMapHandle& mapHandle)
{
	if (factor == 1)
		return mapHandle;
	WeightedMatMultMapHandle ans;
	for (auto &i : mapHandle.mapContents->map){
		BIG_FLOAT v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		BIG_FLOAT one = 1;
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

// Right scalar-multiplication
WeightedMatMultMapHandle operator* (const WeightedMatMultMapHandle& mapHandle, const BIG_FLOAT& factor)
{
	if (factor == 1)
		return mapHandle;
	WeightedMatMultMapHandle ans;
	for (auto &i : mapHandle.mapContents->map){
		BIG_FLOAT v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		BIG_FLOAT one = 1;
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

std::size_t hash_value(const WeightedMatMultMapHandle& val)
{
	return val.mapContents->hashCheck;
}