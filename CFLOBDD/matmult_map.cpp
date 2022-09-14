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

#include <iostream>
#include <fstream>
#include <cassert>
#include "matmult_map.h"
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

//***************************************************************
// MatMultMapBody
//***************************************************************

// Constructor
MatMultMapBody::MatMultMapBody()
	: refCount(0), isCanonical(false)
{
	hashCheck = NULL;
}

void MatMultMapBody::IncrRef()
{
	refCount++;    // Warning: Saturation not checked
}

void MatMultMapBody::DecrRef()
{
	if (--refCount == 0) {    // Warning: Saturation not checked
		if (isCanonical) {
			MatMultMapHandle::canonicalMatMultMapBodySet->DeleteEq(this);
		}
		delete this;
	}
}

unsigned int MatMultMapBody::Hash(unsigned int modsize)
{
	/*if (modsize == HASHSETBASE)
		return hashCheck;*/
	unsigned int hvalue = 0;
	boost::hash<VAL_TYPE> boost_hash;
	for (auto &i : map)
	{
		hvalue = (997 * hvalue + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second))) % modsize;
	}
	return hvalue;
}
void MatMultMapBody::setHashCheck()
{
	long int hvalue = 0;
	boost::hash<VAL_TYPE> boost_hash;
	for (auto &i : map) {
		hvalue = (117 * (hvalue + 1) + (int)(i.first.first + 97 * i.first.second + 97 * 97 * boost_hash(i.second)));
	}
	hashCheck = hvalue;
}

bool MatMultMapBody::operator==(const MatMultMapBody &o) const
{
	if (hashCheck != o.hashCheck) {
		return false;
	}
	else if (map.size() != o.map.size()) {
		return false;
	}
	else {
		for (auto m1_it = map.begin(), m2_it = o.map.begin(); m1_it != map.end() && m2_it != o.map.end(); m1_it++, m2_it++){
			if ((m1_it->first.first != m2_it->first.first) || (m1_it->first.second != m2_it->first.second) || ((m1_it->second != m2_it->second)))
				return false;
		}
	}
	return true;
}

// Overloaded []
VAL_TYPE& MatMultMapBody::operator[](INT_PAIR& p)
{
	return map[p];
}

std::ostream& operator<< (std::ostream & out, const MatMultMapBody &r)
{
	out << "{MMMB: ";
	for (auto &i : r.map) {
		out << "(" << i.first.first << "," << i.first.second << "," << i.second << ");";
	}
	out << " MMMB}";
	return out;
}

//***************************************************************
// MatMultMapHandle
//***************************************************************

// Initializations of static members ---------------------------------
Hashset<MatMultMapBody> *MatMultMapHandle::canonicalMatMultMapBodySet = new Hashset<MatMultMapBody>(HASHSET_NUM_BUCKETS);

// Default constructor
MatMultMapHandle::MatMultMapHandle()
	: mapContents(new MatMultMapBody)
{
	mapContents->IncrRef();
}

// Destructor
MatMultMapHandle::~MatMultMapHandle()
{
	mapContents->DecrRef();
}

// Copy constructor
MatMultMapHandle::MatMultMapHandle(const MatMultMapHandle &r)
	: mapContents(r.mapContents)
{
	mapContents->IncrRef();
}

// Overloaded assignment
MatMultMapHandle& MatMultMapHandle::operator= (const MatMultMapHandle &r)
{
	if (this != &r)      // don't assign to self!
	{
		MatMultMapBody *temp = mapContents;
		mapContents = r.mapContents;
		mapContents->IncrRef();
		temp->DecrRef();
	}
	return *this;
}

// Overloaded !=
bool MatMultMapHandle::operator!=(const MatMultMapHandle &r) const
{
	return (mapContents != r.mapContents);
}

// Overloaded ==
bool MatMultMapHandle::operator==(const MatMultMapHandle &r) const
{

	return (mapContents == r.mapContents);
}

// Overloaded []
VAL_TYPE& MatMultMapHandle::operator[](INT_PAIR& i)
{
	return (*(this->mapContents))[i];
}

// print
std::ostream& MatMultMapHandle::print(std::ostream & out) const
{
	out << *mapContents;
	return out;
}

std::ostream& operator<< (std::ostream & out, const MatMultMapHandle &r)
{
	r.print(out);
	return(out);
}

unsigned int MatMultMapHandle::Hash(unsigned int modsize)
{
	if (!(mapContents->isCanonical)) {
		std::cout << "Hash of a non-canonical LinearMapHandle occurred" << std::endl;
		abort();
		this->Canonicalize();
	}
	assert(mapContents->isCanonical);
	return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

void MatMultMapHandle::Add(const INT_PAIR& p, VAL_TYPE& v)
{
	assert(mapContents->refCount <= 1);
	if (p.first != -1 && p.second != -1){
		auto it = mapContents->map.find(p);
		if (it == mapContents->map.end()){
			mapContents->map.emplace(p, v);
		}
		else{
			it->second += v;
		}
	}
}

void MatMultMapHandle::ForceAdd(const INT_PAIR& p, VAL_TYPE& v)
{
	assert(mapContents->refCount <= 1);
	auto it = mapContents->map.find(p);
	if (it == mapContents->map.end()){
		mapContents->map.emplace(p, v);
		if (p.first == -1 || p.second == -1)
			mapContents->contains_zero_val = true;
	}
	else{
		if (p.first != -1 && p.second != -1)
			it->second += v;
		else{
			assert(it->first.first == -1 && it->first.second == -1);
		}
	}
}

bool MatMultMapHandle::Member(INT_PAIR& y)
{
	return (mapContents->map.find(y) != mapContents->map.end());
}

VAL_TYPE MatMultMapHandle::Lookup(INT_PAIR& p)
{
	return mapContents->map[p];
}

unsigned int MatMultMapHandle::Size(){
	return (unsigned int)mapContents->map.size();
}

std::string MatMultMapHandle::ToString(){
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

void MatMultMapHandle::Canonicalize()
{
	MatMultMapBody *answerContents;
	mapContents->setHashCheck();

	if (!mapContents->isCanonical) {
		unsigned int hash = canonicalMatMultMapBodySet->GetHash(mapContents);
		answerContents = canonicalMatMultMapBodySet->Lookup(mapContents, hash);
		if (answerContents == NULL) {
			canonicalMatMultMapBodySet->Insert(mapContents, hash);
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

// Linear operators on MatMultMapHandles ------------------------------------------------

// Binary addition
MatMultMapHandle MatMultMapHandle::operator+ (const MatMultMapHandle& mapHandle) const
{
	MatMultMapHandle ans;
	for (auto &i : mapContents->map){
		ans.Add(i.first, i.second);
	}
	for (auto &i : mapHandle.mapContents->map){
		ans.Add(i.first, i.second);
	}
	if (ans.Size() == 0){
		VAL_TYPE one = 1;
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

// Left scalar-multiplication
MatMultMapHandle operator* (const VAL_TYPE& factor, const MatMultMapHandle& mapHandle)
{
	if (factor == 1)
		return mapHandle;
	MatMultMapHandle ans;
	for (auto &i : mapHandle.mapContents->map){
		VAL_TYPE v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		VAL_TYPE one = 1;
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

// Right scalar-multiplication
MatMultMapHandle operator* (const MatMultMapHandle& mapHandle, const VAL_TYPE& factor)
{
	if (factor == 1)
		return mapHandle;
	MatMultMapHandle ans;
	for (auto &i : mapHandle.mapContents->map){
		VAL_TYPE v = i.second * factor;
		ans.Add(i.first, v);
	}
	if (ans.Size() == 0){
		VAL_TYPE one = 1;
		ans.ForceAdd(std::make_pair(-1, -1), one);
	}
	ans.Canonicalize();
	return ans;
}

std::size_t hash_value(const MatMultMapHandle& val)
{
	return val.mapContents->hashCheck;
}