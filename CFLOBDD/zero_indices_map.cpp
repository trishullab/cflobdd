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
#include "zero_indices_map.h"

//***************************************************************
// ZeroIndicesMapBody
//***************************************************************

// Constructor
ZeroIndicesMapBody::ZeroIndicesMapBody()
	: refCount(0), isCanonical(false)
{
	hashCheck = NULL;
}

void ZeroIndicesMapBody::IncrRef()
{
	refCount++;    // Warning: Saturation not checked
}

void ZeroIndicesMapBody::DecrRef()
{
	if (--refCount == 0) {    // Warning: Saturation not checked
		if (isCanonical) {
			ZeroIndicesMapHandle::canonicalZeroIndicesMapBodySet->DeleteEq(this);
		}
		delete this;
	}
}

unsigned int ZeroIndicesMapBody::Hash(unsigned int modsize)
{
	/*if (modsize == HASHSETBASE)
	return hashCheck;*/
	unsigned int hvalue = (997 * a_index) % modsize;
	for (auto &i : b_indices)
	{
		hvalue = (997 * hvalue + (int)(i)) % modsize;
	}
	return hvalue;
}
void ZeroIndicesMapBody::setHashCheck()
{
	unsigned int hvalue = 117 + 97 * a_index;
	for (auto &i : b_indices) {
		hvalue = (117 * (hvalue + 1) + (int)(i));
	}
	hashCheck = hvalue;
}

bool ZeroIndicesMapBody::operator==(const ZeroIndicesMapBody &o) const
{
	if (hashCheck != o.hashCheck) {
		return false;
	}
	else if (b_indices.size() != o.b_indices.size()) {
		return false;
	}
	else if (a_index != o.a_index){
		return false;
	}
	else {
		for (unsigned int i = 0; i < b_indices.size(); i++){
			if (b_indices[i] != o.b_indices[i])
				return false;
		}
	}
	return true;
}

// Overloaded []
int ZeroIndicesMapBody::operator[](int p)
{
	return b_indices[p];
}

std::ostream& operator<< (std::ostream & out, const ZeroIndicesMapBody &r)
{
	out << "{MMMB: ";
	out << r.a_index << ",";
	for (auto &i : r.b_indices) {
		out << i << ",";
	}
	out << " MMMB}";
	return out;
}

//***************************************************************
// ZeroIndicesMapHandle
//***************************************************************

// Initializations of static members ---------------------------------
Hashset<ZeroIndicesMapBody> *ZeroIndicesMapHandle::canonicalZeroIndicesMapBodySet = new Hashset<ZeroIndicesMapBody>(HASHSET_NUM_BUCKETS);

// Default constructor
ZeroIndicesMapHandle::ZeroIndicesMapHandle()
	: mapContents(new ZeroIndicesMapBody)
{
	mapContents->IncrRef();
}

// Destructor
ZeroIndicesMapHandle::~ZeroIndicesMapHandle()
{
	mapContents->DecrRef();
}

// Copy constructor
ZeroIndicesMapHandle::ZeroIndicesMapHandle(const ZeroIndicesMapHandle &r)
	: mapContents(r.mapContents)
{
	mapContents->IncrRef();
}

// Overloaded assignment
ZeroIndicesMapHandle& ZeroIndicesMapHandle::operator= (const ZeroIndicesMapHandle &r)
{
	if (this != &r)      // don't assign to self!
	{
		ZeroIndicesMapBody *temp = mapContents;
		mapContents = r.mapContents;
		mapContents->IncrRef();
		temp->DecrRef();
	}
	return *this;
}

// Overloaded !=
bool ZeroIndicesMapHandle::operator!=(const ZeroIndicesMapHandle &r) const
{
	return (mapContents != r.mapContents);
}

// Overloaded ==
bool ZeroIndicesMapHandle::operator==(const ZeroIndicesMapHandle &r) const
{

	return (mapContents == r.mapContents);
}

// Overloaded []
int ZeroIndicesMapHandle::operator[](int i)
{
	return (*(this->mapContents))[i];
}

int ZeroIndicesMapHandle::Get_AIndex(){
	return mapContents->a_index;
}

void ZeroIndicesMapHandle::Set_AIndex(int i){
	mapContents->a_index = i;
}

// print
std::ostream& ZeroIndicesMapHandle::print(std::ostream & out) const
{
	out << *mapContents;
	return out;
}

std::ostream& operator<< (std::ostream & out, const ZeroIndicesMapHandle &r)
{
	r.print(out);
	return(out);
}

unsigned int ZeroIndicesMapHandle::Hash(unsigned int modsize)
{
	if (!(mapContents->isCanonical)) {
		std::cout << "Hash of a non-canonical LinearMapHandle occurred" << std::endl;
		abort();
		this->Canonicalize();
	}
	assert(mapContents->isCanonical);
	return ((unsigned int) reinterpret_cast<uintptr_t>(mapContents) >> 2) % modsize;
}

void ZeroIndicesMapHandle::Add_BIndex(int v)
{
	assert(mapContents->refCount <= 1);
	mapContents->b_indices.push_back(v);
}

bool ZeroIndicesMapHandle::Member(int y)
{
	for (auto i : mapContents->b_indices){
		if (i == y)
			return true;
	}
	return false;
}

int ZeroIndicesMapHandle::Lookup(int p)
{
	return mapContents->b_indices[p];
}

unsigned int ZeroIndicesMapHandle::Size(){
	return (unsigned int)mapContents->b_indices.size();
}

void ZeroIndicesMapHandle::Canonicalize()
{
	ZeroIndicesMapBody *answerContents;
	mapContents->setHashCheck();

	if (!mapContents->isCanonical) {
		unsigned int hash = canonicalZeroIndicesMapBodySet->GetHash(mapContents);
		answerContents = canonicalZeroIndicesMapBodySet->Lookup(mapContents, hash);
		if (answerContents == NULL) {
			canonicalZeroIndicesMapBodySet->Insert(mapContents, hash);
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

std::size_t hash_value(const ZeroIndicesMapHandle& val)
{
	return val.mapContents->hashCheck;
}
