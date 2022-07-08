//#include "uw_assert.h"
#include <algorithm> 
#include <iomanip>
#include "../bit_vector/bit_vector_1.h"
#include "../bit_vector/bit_vector_ops.h"
#include "sparse_tarray.h"

namespace HowellMatrix
{
	//-------------------------------------------------
	// SparseTArray: Sparse aray representation.
	//-------------------------------------------------

	template <typename T>
	SparseTArray<T>::SparseTArray(size_t l)
		: len_(l)
	{ }

	template <typename T>
	SparseTArray<T>::SparseTArray(const SparseTArray &ia)
		: len_(ia.len_)
	{
		arr_ = ia.arr_;
	}

	template <typename T>
	SparseTArray<T>::SparseTArray(const SparseTArray* ia)
		: len_(ia->len_)
	{
		arr_ = ia->arr_;
	}

	template <typename T>
	SparseTArray<T>::SparseTArray(const TArray<T> & arr)
		: len_(arr.getLength())
	{
		for (size_t i = 0; i < len_; i++)
		{
			if (arr[i] != T(0))
			{
				arr_.push_back(std::pair<size_t, T>(i, arr[i]));
			}
		}
	}

	template <typename T>
	SparseTArray<T>::SparseTArray(const T* ia, size_t length)
		: len_(length)
	{
		for (size_t i = 0; i < length; i++)
		{
			if (ia[i] != T(0))
			{
				arr_.push_back(std::pair<size_t, T>(i, ia[i]));
			}
		}
	}

	template <typename T>
	SparseTArray<T>::SparseTArray(const std::vector<T> & a)
		: len_(a.size())
	{
		for (size_t i = 0; i < len_; i++)
		{
			if (a[i] != T(0))
			{
				arr_.push_back(std::pair<size_t, T>(i, a[i]));
			}
		}
	}

	template <typename T>
	SparseTArray<T> SparseTArray<T>::operator=(const SparseTArray &ia) {
		if (this != &ia) {
			len_ = ia.len_;
			arr_ = ia.arr_;
		}
		return *this;
	}

	template <typename T>
	SparseTArray<T>::~SparseTArray() {
	}

	template <typename T>
	bool SparseTArray<T>::operator==(const SparseTArray &that) const {
		if (len_ != that.len_)
		{
			return false;
		}

		return arr_ == that.arr_;
	}

	template <typename T>
	bool SparseTArray<T>::operator!=(const SparseTArray &that) const {
		return !(*this == that);
	}

	template <typename T>
	inline void SparseTArray<T>::set(size_t index, T val) {
#ifdef UW_ASSERTS_ENABLED
		if (index >= len_)
			UWAssert::UWAssert::shouldNeverHappen();
#endif
		typename std::vector<std::pair<size_t, T> >::iterator it =
			std::lower_bound(arr_.begin(), arr_.end(), std::pair<size_t, T>(index, T(0)), firstKeyPairVectorComparator<T>());
		if (it != arr_.end() && it->first == index)
		{
			if (val != 0)
				it->second = val;
			else
				arr_.erase(it);
		}
		else
		{
			if (val != 0)
			{
				arr_.push_back(std::pair<size_t, T>(index, val));
				std::sort(arr_.begin(), arr_.end(), firstKeyPairVectorComparator<T>());
			}
		}
	}

	template <typename T>
	inline size_t SparseTArray<T>::Size() const {
		return arr_.size();
	}

	template <typename T>
	inline T SparseTArray<T>::Get(size_t index) const {
#ifdef UW_ASSERTS_ENABLED
		if (index >= len_)
			UWAssert::UWAssert::shouldNeverHappen();
#endif

#if 0
		// Using lower bound for smaller std::vector leads to lot of inefficiency
		typename std::vector<std::pair<size_t, T> >::const_iterator it =
			std::lower_bound(arr_.begin(), arr_.end(), std::make_pair<size_t, T>(index, T(0)), firstKeyPairVectorComparator<T>());

		if (it != arr_.end() && it->first == index)
		{
			return it->second;
		}
		else
		{
			return T(0);
		}
#else
		for (typename std::vector<std::pair<size_t, T> >::const_iterator it = arr_.begin(); it != arr_.end(); it++)
		{
			if (it->first == index) {
				return it->second;
			}
			else {
				if (it->first > index)
					return T(0);
			}
		}

		return T(0);
#endif
	}

	template <typename T>
	inline size_t SparseTArray<T>::GetLength() const {
		return len_;
	}

	template <typename T>
	inline size_t SparseTArray<T>::GetNumNonZeroEntries() const {
		return arr_.size();
	}

	template <typename T>
	std::pair<T, T> SparseTArray<T>::GetAbsoluteSumOfCoefficients() const {
		T carry = T(0);
		T abs_sum = T(0);
		for (typename std::vector<std::pair<size_t, T> >::const_iterator it = arr_.begin(); it != arr_.end(); it++)
		{
			T max_signed_int = (T(1) << (T::highest_power() - 1)) - 1;
			T abs_val_it = (it->second <= max_signed_int) ? it->second : (T(0) - (it->second));

			// Check for overflow
			if (abs_sum + abs_val_it > max_signed_int)
			{
				carry = carry + 1;
				abs_sum = abs_sum - max_signed_int;
			}

			abs_sum = abs_sum + abs_val_it;
		}
		UWAssert::UWAssert::shouldNeverHappen(abs_sum < T(0));
		return std::pair<T, T>(carry, abs_sum);
	}

	template <typename T>
	const std::vector<std::pair<size_t, T> > SparseTArray<T>::GetArray() const {
		return arr_;
	}

	template <typename T>
	inline size_t SparseTArray<T>::GetLeadingIndex() const {
		if (arr_.empty())
			return len_;

		typename std::vector<std::pair<size_t, T> >::const_iterator it = arr_.begin();
		return it->first;
	}

	template <typename T>
	inline size_t SparseTArray<T>::GetLeadingRank() const {
		// std::cout << "arr size: " << arr_.size() << std::endl;
		if (arr_.empty())
		{
			return BitVector::highest_power<T>();
		}

		return BitVector::compute_rank(arr_.begin()->second);
	}

	template <typename T>
	inline T SparseTArray<T>::operator [](size_t j) const {
#ifdef UW_ASSERTS_ENABLED
		if (j >= len_)
			UWAssert::UWAssert::shouldNeverHappen();
#endif
		return Get(j);
	}

	template <typename T>
	TArray<T>* SparseTArray<T>::GetTArray() const
	{
		TArray<T>* ret = new TArray<T>(len_);
		typename std::vector<std::pair<size_t, T> >::const_iterator it = arr_.begin();
		for (size_t i = 0; i < len_; i++)
		{
			if (it != arr_.end() && it->first == i)
			{
				ret->set(i, it->second);
				it++;
			}
			else
			{
				ret->set(i, 0);
			}
		}
		return ret;
	}

	// insertVals: insert `length` copies of `val` starting at index `begin`.
	// This inserts values, where setRange overwrites them.
	// So, this increases `len`.
	template <typename T>
	inline void SparseTArray<T>::insertVals(size_t beg_index, size_t length, T val) {
		if (length == 0)
			return;

		typename std::vector<std::pair<size_t, T> >::iterator it = (beg_index != 0) ? std::upper_bound(arr_.begin(), arr_.end(), std::make_pair<size_t, T>(beg_index - 1, 0), firstKeyPairVectorComparator<T>()) : arr_.begin();
		for (; it != arr_.end(); it++)
		{
			it->first = it->first + length;
		}

		len_ = len_ + length;

		//Add the new values
		if (val != 0)
		{
			for (size_t i = 0; i < length; i++)
			{
				arr_.push_back(std::pair<size_t, T>(beg_index + i, val));
			}
			std::sort(arr_.begin(), arr_.end(), firstKeyPairVectorComparator<T>());
		}
	}

	// insertVals: insert `length` copies of `val` starting at index `begin`.
	// This inserts values, where setRange overwrites them.
	// So, this increases `len`.
	template <typename T>
	inline void SparseTArray<T>::insertVals(size_t beg_index, size_t length, SparseTArray<T>& vals, size_t vals_beg_index) {
		if (length == 0)
			return;

		typename std::vector<std::pair<size_t, T> >::iterator it = (beg_index != 0) ? std::upper_bound(arr_.begin(), arr_.end(), std::pair<size_t, T>(beg_index - 1, 0), firstKeyPairVectorComparator<T>()) : arr_.begin();
		for (; it != arr_.end(); it++)
		{
			it->first = it->first + length;
		}

		len_ = len_ + length;

		//Add the new values
		typename std::vector<std::pair<size_t, T> >::const_iterator vals_it;
		for (vals_it = vals.arr_.begin(); vals_it != vals.arr_.end(); vals_it++)
		{
			if (vals_it->first >= vals_beg_index && vals_it->first < vals_beg_index + length)
				arr_.push_back(std::pair<size_t, T>(beg_index + (vals_it->first - vals_beg_index), vals_it->second));
		}
		std::sort(arr_.begin(), arr_.end(), firstKeyPairVectorComparator<T>());
	}

	//TODO: Inefficient because it actually creates a new map.
	template <typename T>
	void SparseTArray<T>::project(const std::vector<size_t> & projected_indices)
	{
		if (projected_indices.size() == 0)
		{
			len_ = 0;
			arr_.clear();
			return;
		}
		else if (projected_indices.size() == len_)
		{
			return;
		}

		std::vector<std::pair<size_t, T> > new_arr;

		typename std::vector<std::pair<size_t, T> >::iterator it = arr_.begin();

		size_t crt_index = 0;
		for (size_t i = 0; i < len_; i++)
		{
			if (std::find(projected_indices.begin(), projected_indices.end(), i)
				!= projected_indices.end())
			{
				if (it->first == i)
				{
					new_arr.push_back(std::pair<size_t, T>(crt_index, it->second));
				}
				crt_index++;
			}
			if (it->first == i)
				it++;
		}

		len_ = projected_indices.size();
#ifdef UW_ASSERTS_ENABLED
		UWAssert::UWAssert::shouldNeverHappen(len_ != crt_index);
#endif
		arr_ = new_arr;
	}

	template <typename T>
	void SparseTArray<T>::projectOut(const std::vector<size_t> & project_out_indices)
	{
		if (project_out_indices.size() == 0)
			return;

		std::vector<size_t> projected_indices;

		for (size_t i = 0; i < len_; i++)
		{
			if (std::find(project_out_indices.begin(), project_out_indices.end(), i) == project_out_indices.end())
				projected_indices.push_back(i);
		}

		project(projected_indices);
	}

	template <typename T>
	void SparseTArray<T>::projectRight(size_t beg_index)
	{
		if (beg_index == 0)
			return;

		typename std::vector<std::pair<size_t, T> >::iterator it = (beg_index != 0) ? std::upper_bound(arr_.begin(), arr_.end(), std::make_pair<size_t, T>(beg_index - 1, 0), firstKeyPairVectorComparator<T>()) : arr_.begin();

		for (; it != arr_.end(); it++)
		{
			it->first = it->first - beg_index;
		}

		len_ = len_ - beg_index;
	}

	template <typename T>
	void SparseTArray<T>::permuteCols(const std::vector<size_t> & from,
		const std::vector<size_t> & to)
	{
		UWAssert::UWAssert::shouldNeverHappen(from.size() != to.size());

		if (from.size() == 0)
			return;

		typename std::vector<std::pair<size_t, T> >::iterator it;

		for (it = arr_.begin(); it != arr_.end(); it++)
		{
			typename std::vector<size_t>::const_iterator from_loc =
				std::find(from.begin(), from.end(), it->first);
			if (from_loc != from.end())
			{
				it->first = to[from_loc - from.begin()];
			}
		}
		std::sort(arr_.begin(), arr_.end(), firstKeyPairVectorComparator<T>());
	}

	template <typename T>
	void SparseTArray<T>::vectorMul(SparseTArray<T>*& out, T mult)
	{
		if (mult == 1)
			return;

		if (mult == 0)
			out->arr_.clear();

		typename std::vector<std::pair<size_t, T> >::iterator it;

		for (it = out->arr_.begin(); it != out->arr_.end();)
		{
			T val = it->second * mult;
			it->second = val;

			if (val == 0)
				it = out->arr_.erase(it);
			else
				++it;
		}
	}


	template <typename T>
	void SparseTArray<T>::vectorLinearCombination(SparseTArray<T>*& out, const SparseTArray<T>& in1, T mult1, const SparseTArray<T>& in2, T mult2)
	{
#ifdef UW_ASSERTS_ENABLED
		UWAssert::UWAssert::shouldNeverHappen(in1.GetLength() != in2.GetLength());
		UWAssert::UWAssert::shouldNeverHappen(out.GetLength() != in1.GetLength());
#endif
		typename std::vector<std::pair<size_t, T> >::const_iterator it1 = in1.arr_.begin();
		typename std::vector<std::pair<size_t, T> >::const_iterator it2 = in2.arr_.begin();

		std::vector<std::pair<size_t, T> > new_arr;
		while (it1 != in1.arr_.end() || it2 != in2.arr_.end())
		{
			if (it1 != in1.arr_.end() && it2 != in2.arr_.end() && it1->first == it2->first)
			{
				T expr = mult1 * it1->second + mult2 * it2->second;
				if (expr != 0)
				{
					new_arr.push_back(std::pair<size_t, T>(it1->first, expr));
				}
				it1++;
				it2++;
			}
			else if ((it2 == in2.arr_.end()) || ((it1 != in1.arr_.end()) && (it1->first < it2->first)))
			{
				T expr = mult1 * it1->second;
				if (expr != 0)
				{
					new_arr.push_back(std::pair<size_t, T>(it1->first, expr));
				}
				it1++;
			}
			else
			{
				T expr = mult2 * it2->second;
				if (expr != 0)
				{
					new_arr.push_back(std::pair<size_t, T>(it2->first, expr));
				}
				it2++;
			}

		}
		out->arr_ = new_arr;
	}

	template <typename T>
	std::ostream & SparseTArray<T>::Print(std::ostream & out) const {
		//TODO: Can improve the performance if iterators are used directly.
		out << "  ";
		for (size_t i = 0; i < GetLength(); i++) {
			out << Get(i) << " ";
		}
		return out;
	}

}

