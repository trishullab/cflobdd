#include "../assert/uw_assert.h"
#include <algorithm> 
#include <iomanip>
#include <cstring>
#include "../bit_vector/bit_vector_1.h"
#include "../bit_vector/bit_vector_ops.h"
#include "ModularRow.h"

namespace HowellMatrix
{
	//-------------------------------------------------
	// TArray: Just a wrapper around T* that does array
	// bound checking.
	//-------------------------------------------------

	template <typename T>
	TArray<T>::TArray(size_t l)
		: len(l), m(new T[l])
	{ }

	template <typename T>
	TArray<T>::TArray(const TArray &ia)
		: len(ia.len), m(new T[ia.len])
	{
		memcpy(m, ia.m, len * sizeof(T));
	}

	template <typename T>
	TArray<T>::TArray(TArray* ia)
		: len(ia->len), m(new T[ia->len])
	{
		memcpy(m, ia->m, len * sizeof(T));
	}

	template <typename T>
	TArray<T>::TArray(T* ia, size_t length)
		: len(length), m(new T[length])
	{
		memcpy(m, ia, len * sizeof(T));
	}

	template <typename T>
	TArray<T> TArray<T>::operator=(const TArray &ia) {
		if (this != &ia) {
			len = ia.len;
			delete[] m;
			m = new T[len];
			memcpy(m, ia.m, len * sizeof(T));
		}
		return *this;
	}

	template <typename T>
	TArray<T>::~TArray() {
		delete[] m;
	}

	template <typename T>
	bool TArray<T>::operator==(const TArray &t) const {
		if (len != t.len)
		{
			return false;
		}

		const T * t_arr = t.getConstArray();

		for (size_t i = 0; i < len; i++)
		{
			if (m[i] != t_arr[i])
				return false;
		}

		return true;
	}

	template <typename T>
	inline void TArray<T>::set(size_t index, T val) {
#ifdef UW_ASSERTS_ENABLED
		if (index >= len)
			UWAssert::UWAssert::shouldNeverHappen();
#endif
		m[index] = val;
	}


	template <typename T>
	inline T TArray<T>::get(size_t index) const {
#ifdef UW_ASSERTS_ENABLED
		if (index >= len)
			UWAssert::UWAssert::shouldNeverHappen();
#endif
		return m[index];
	}

	template <typename T>
	size_t TArray<T>::getLeadingIndex() const {
		size_t i = 0;
		while (i < len && m[i] == 0) {
			++i;
		}
		return i;
	}

	template <typename T>
	inline T & TArray<T>::operator [](size_t j) const {
#ifdef UW_ASSERTS_ENABLED
		if (j >= len)
			UWAssert::UWAssert::shouldNeverHappen();
#endif
		return m[j];
	}

	//Set m[begin] through m[begin + length - 1] to value val.
	template <typename T>
	inline void TArray<T>::setRange(size_t begin, size_t length, T val) {
		for (size_t i = begin; i < begin + length; i++)
			set(i, val);
	}

	// setRange: Set m[begin] through m[begin + length - 1] to value
	// valArray[val_begin] to valArray[val_begin + length - 1].
	template <typename T>
	inline void TArray<T>::setRange(size_t begin, size_t length,
		const T* valArray, size_t val_begin) {
		for (size_t i = 0; i < length; i++)
			set(begin + i, valArray[val_begin + i]);
	}

	// insertVals: insert `length` copies of `val` starting at index `begin`.
	// This inserts values, where setRange overwrites them.
	// So, this increases `len`, and the size of the array `m`.
	template <typename T>
	inline void TArray<T>::insertVals(const size_t begin, const size_t length, const T val) {
		size_t new_len = len + length;
		T * new_m = new T[new_len];

		// Insert values before, inside, and after the new range of values.
		memcpy(new_m, m, begin*sizeof(T));
		for (size_t i = 0; i < length; ++i) { new_m[begin + i] = val; }
		memcpy(new_m + (length + begin), m + (begin), (len - begin)*sizeof(T));

		delete[] m;
		m = new_m;
		len = new_len;
	}

	template <typename T>
	inline void TArray<T>::insertRange(const size_t begin, const size_t length,
		const T* valArray, size_t val_begin) {
		size_t new_len = len + length;
		T * new_m = new T[new_len];

		// Insert values before, inside, and after the new range of values.
		memcpy(new_m, m, begin*sizeof(T));
		memcpy(new_m + begin, valArray + val_begin, length*sizeof(T));
		memcpy(new_m + (length + begin), m + begin, (len - begin)*sizeof(T));

		delete[] m;
		m = new_m;
		len = new_len;
	}

	template <typename T>
	inline void TArray<T>::insertRange(const size_t begin, const size_t length,
		const TArray & valArray, size_t val_begin) {
		insertRange(begin, length, valArray.m, val_begin);
	}

	template <typename T>
	TArray<T>* TArray<T>::projectOut(const std::vector<size_t> & projectOutIndices) const
	{
		TArray<T>* ret = new TArray(len - projectOutIndices.size());
		size_t curIndex = 0;
		for (size_t i = 0; i < len; i++)
		{
			if (std::find(projectOutIndices.begin(), projectOutIndices.end(), i) == projectOutIndices.end())
			{
				ret->set(curIndex, get(i));
				curIndex++;
			}
		}
		UWAssert::UWAssert::shouldNeverHappen(curIndex != (len - projectOutIndices.size()));
		return ret;
	}

	template <typename T>
	std::ostream & TArray<T>::print(std::ostream & out) const {
		for (size_t i = 0; i < getLength(); i++) {
			out << std::setw(4) << get(i);
			if (i>0) out << " ";
		}
		return out;
	}



	//----------------------------------
	// Multiply two matrices
	// TODO: DEPRECATE.
	//----------------------------------
	template <typename T>
	void multiplyMatrices(TArray<T>* multResult, const TArray<T>* op1, const TArray<T>* op2, unsigned N)
	{
		unsigned i, j, k;
		for (i = 0; i < N; ++i) {
			for (j = 0; j < N; ++j) {
				multResult->set(i * N + j, 0);
				for (k = 0; k < N; ++k) {
					T val = multResult->get(i * N + j) + op1->get(i * N + k) * op2->get(k * N + j);
					multResult->set(i * N + j, val);
				}
			}
		}
	}



	/*---------------------------------------------------------------
	*  ModularRow represents a row in HowellMatrix.
	*
	*  It wraps TArray and maintains leading index and leading rank.
	*--------------------------------------------------------------*/

	template <typename T>
	void ModularRow<T>::init(const T* _m, size_t n) {
		m = new TArray<T>(_m, n);
		set_leading();
	}

	// i, if given, is an index no less than the new leading index.
	template <typename T>
	void ModularRow<T>::set_leading(size_t i) {
		size_t n = m->getLength();
		while (i < n && m->get(i) == 0) {
			++i;
		}
		if (i == n) {
			leading_index = n;
			leading_rank = T::highest_power();
		}
		else {
			leading_index = i;
			leading_rank = BitVector::compute_rank(m->get(i));
		}
	}

	//Constructors
	template <typename T>
	ModularRow<T>::ModularRow(size_t l) {
		m = new TArray<T>(l);
		set_leading();
	}

	template <typename T>
	ModularRow<T>::ModularRow(const ModularRow &a) {
		m = a.m;
		leading_index = a.leading_index;
		leading_rank = a.leading_rank;
	}
	template <typename T>
	ModularRow<T>::ModularRow(const T* a, size_t len) {
		init(a, len);
	}

	template <typename T>
	ModularRow<T>::ModularRow(const std::vector<T> & vec)
		: m(new TArray<T>(vec.size()))
	{
		for (size_t i = 0; i < vec.size(); i++) {
			m->set(i, vec[i]);
		}
		set_leading();
	}

	template <typename T>
	T * ModularRow<T>::getArray() {
		return m->getArray();
	}
	template <typename T>
	const T* ModularRow<T>::getConstArray() const{
		return m->getConstArray();
	}
	template <typename T>
	TArray<T>* ModularRow<T>::getTArray() {
		return m;
	}
	template <typename T>
	TArray<T>* ModularRow<T>::getConstTArray() const{
		return new TArray<T>(*m);
	}
	template <typename T>
	void ModularRow<T>::setMat(TArray<T>* mat) {
		m = mat;
	}

	template <typename T>
	size_t ModularRow<T>::getSize() const {
		return m->getLength();
	}
	template <typename T>
	size_t ModularRow<T>::getLeadingIndex() const
	{
		return leading_index;
	}
	template <typename T>
	void ModularRow<T>::setLeadingIndex(size_t leadingIndex)
	{
		leading_index = leadingIndex;
	}
	template <typename T>
	size_t ModularRow<T>::getLeadingRank() const
	{
		return leading_rank;
	}
	template <typename T>
	void ModularRow<T>::setLeadingRank(size_t leadingRank)
	{
		leading_rank = leadingRank;
	}
	template <typename T>
	ModularRow<T> & ModularRow<T>::operator=(const ModularRow &a) {
		if (this != &a) {
			m = a.m;
			leading_index = a.leading_index;
			leading_rank = a.leading_rank;
		}
		return *this;
	}
	template <typename T>
	bool ModularRow<T>::operator== (const ModularRow & mat) const {
		if (m->getLength() != mat.m->getLength()) return false;
		if (leading_index != mat.leading_index) return false;
		if (leading_rank != mat.leading_rank) return false;
		size_t i;
		for (i = leading_index; i < m->getLength(); i++)
		{
			if (m->get(i) != mat.m->get(i)) return false;
		}
		return true;
	}
	template <typename T>
	ModularRow<T>::ModularRow(const TArray<T>* & _m, size_t n) {
		m = _m;

		const T * mat = m->getConstArray();
		size_t i = 0;
		while (i < n && mat[i] == 0) {
			i++;
		}
		if (i == n) {
			leading_index = n;
			leading_rank = T::highest_power();
		}
		else {
			leading_index = i;
			leading_rank = BitVector::compute_rank(mat[i]);
		}
	}
	template <typename T>
	ModularRow<T>::ModularRow(const TArray<T>* & _m, size_t _li, size_t _lr) {
		m = new TArray<T>(_m->getConstArray(), _m->getLength());
		leading_index = _li;
		leading_rank = _lr;
	}
	template <typename T>
	ModularRow<T>::~ModularRow() {
		m = NULL;
	}

	template <typename T>
	const T * ModularRow<T>::getConstMatrix() const {
		return m->getConstArray();
	}

	template <typename T>
	ModularRow<T>* ModularRow<T>::copy(const ModularRow<T>* a)
	{
		TArray<T>* tmp = a->m;
		size_t n = tmp->getLength();

		TArray<T>* ans = new TArray<T>(n);
		for (size_t i = 0; i<n; i++)
			ans->set(i, tmp->get(i));

		ModularRow* ret = new ModularRow(ans, n);
		return ret;
	}

	template <typename T>
	std::vector<ModularRow<T> > ModularRow<T>::copy(const std::vector<ModularRow<T> > & v)
	{
		typename std::vector<ModularRow> ans;
		typename std::vector<ModularRow>::const_iterator it;
		for (it = v.begin(); it != v.end(); it++) {
			ModularRow t = copy(*it);
			ans.push_back(t);
		}
		return ans;
	}

	template <typename T>
	std::vector<ModularRow<T>* > ModularRow<T>::copy(const std::vector<ModularRow<T>* > & v)
	{
		typename std::vector<ModularRow* > ans;
		typename std::vector<ModularRow* >::const_iterator it;

		for (it = v.begin(); it != v.end(); it++) {
			ModularRow* t = copy(*it);
			ans.push_back(t);
		}
		return ans;
	}

	template <typename T>
	inline T ModularRow<T>::get(size_t index) const {
		return m->get(index);
	}

	template <typename T>
	inline void ModularRow<T>::set(size_t index, T val) {
		m->set(index, val);
	}

	template <typename T>
	inline void ModularRow<T>::setVal(size_t begin, size_t length, const T val) {
		m->setRange(begin, length, val);
		if (leading_index >= begin) {
			if (val == 0) {
				leading_index += length;
			}
			else {
				leading_index = begin;
				leading_rank = BitVector::compute_rank((*m)[begin]);
			}
		}
	}

	template <typename T>
	inline void ModularRow<T>::setRange(size_t begin, size_t length,
		const ModularRow & val_array, size_t val_begin) {
		m->setRange(begin, length, val_array.m, val_begin);
		if (leading_index >= begin) { set_leading(begin); }
	}

	template <typename T>
	inline void ModularRow<T>::insertVal(const size_t begin, const size_t length, const T val) {
		m->insertVals(begin, length, val);

		// If we might have changed the leading index, fix it up.
		if (leading_index >= begin) {
			if (val == 0) {
				set_leading(begin + length);
			}
			else {
				leading_index = begin;
				leading_rank = BitVector::compute_rank((*m)[begin]);
			}
		}
	}

	template <typename T>
	inline void ModularRow<T>::insertRange(const size_t begin, const size_t length,
		const ModularRow& valArray, size_t val_begin) {
		m->insertRange(begin, length, *valArray.m, val_begin);

		// If we might have changed the leading index, fix it up.
		if (leading_index >= begin) { set_leading(begin); }

	}

	template <typename T>
	inline void ModularRow<T>::swapCols(const size_t a, const size_t b) {
		UWAssert::UWAssert::shouldNeverHappen(a > getSize() || b > getSize());
		if (a == b) return;

		swap((*m)[a], (*m)[b]);

		// Fix leading index and rank.
		size_t low = a<b ? a : b;
		size_t high = a<b ? b : a;
		if (leading_index == high) {
			leading_index = low;
		}
		else if (leading_index < low || leading_index > high) {
			return;
		}
		else { // low <= leading_index < high.
			set_leading(low);
		}
	}

	template <typename T>
	inline void ModularRow<T>::permuteCols(const std::vector<size_t> & from,
		const std::vector<size_t> & to) {
		UWAssert::UWAssert::shouldNeverHappen(from.size() != to.size());
		T *tmp = new T[m->getLength()];
		for (size_t i = 0; i < from.size(); ++i) {
			tmp[to[i]] = m->get(from[i]);
		}
		for (size_t i = 0; i < to.size(); ++i) {
			m->set(to[i], tmp[to[i]]);
		}
		set_leading();
		delete[] tmp;
	}

	template <typename T>
	std::ostream & ModularRow<T>::print(std::ostream & out) const {
		return m->print(out);
	}

	// Out = in1*mult1 + in2*mult2
	template <typename T>
	void vectorLinearCombination(ModularRow<T>& out,
		ModularRow<T>& in1, T mult1,
		ModularRow<T>& in2, T mult2) {
		UWAssert::UWAssert::shouldNeverHappen(out.getSize() != in1.getSize());
		UWAssert::UWAssert::shouldNeverHappen(in1.getSize() != in2.getSize());

		size_t i;
		T * out_vec = out.getArray();
		T * in1_vec = in1.getArray();
		T * in2_vec = in2.getArray();
		bool set = false;
		for (i = 0; i < out.getSize(); i++)
		{
			out_vec[i] = in1_vec[i] * mult1 + in2_vec[i] * mult2;
			if (out_vec[i] != 0 && !set)
			{
				set = true;
				out.setLeadingIndex(i);
				out.setLeadingRank(BitVector::compute_rank(out_vec[i]));
			}
		}
		if (!set)
		{
			out.setLeadingIndex(out.getSize());
			out.setLeadingRank(T::highest_power());
		}
	}
}
 
