#ifndef MATRIX_SPARSE_TARRAY
#define MATRIX_SPARSE_TARRAY

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "ModularRow.h"


	//-------------------------------------------------
	// SparseTArray: Sparse representation of an Array over T.
	//-------------------------------------------------

namespace HowellMatrix
{
	template <typename T>
	class SparseTArray {

	public:
		SparseTArray(size_t l);
		SparseTArray(const SparseTArray &ia);  // Copy constructor.
		SparseTArray(const SparseTArray* ia);  // Copy constructor.
		SparseTArray(const TArray<T> & arr);
		SparseTArray(const std::vector<T> & arr);
		SparseTArray(const T* ia, size_t length);
		SparseTArray operator=(const SparseTArray &ia);
		~SparseTArray();

		void DeallocateMemory(){  }

		bool operator==(const SparseTArray &t) const;
		bool operator!=(const SparseTArray &t) const;

		size_t GetLength() const;
		size_t GetNumNonZeroEntries() const;
		std::pair<T, T> GetAbsoluteSumOfCoefficients() const;

		const std::vector<std::pair<size_t, T> > GetArray() const;

		void set(size_t index, T val);
		T Get(size_t index) const;
		size_t Size() const;
		T operator [](size_t j) const;

		size_t GetLeadingIndex() const;
		size_t GetLeadingRank() const;

		TArray<T>* GetTArray() const;
		const std::vector<std::pair<size_t, T> > GetSparseVector() const
		{
			return arr_;
		}

		// insertVals: Insert `length` copies of `val`, starting at index `begin`.
		// This inserts values, and so increases the size of the array.
		void insertVals(size_t begin, size_t length = 1,
			T val = 0);
		void insertVals(size_t begin, size_t length, SparseTArray<T>& vals, size_t vals_beg_index = 0);

		void project(const std::vector<size_t> & project_indices);
		void projectOut(const std::vector<size_t> & proj_out_indices);
		//Project on the variables on the right of i including i
		void projectRight(size_t beg_index);

		void permuteCols(const std::vector<size_t> & from,
			const std::vector<size_t> & to);

		// TODO: Add static members for arithmetic operations like addition and vectorCombination.
		static void vectorMul(SparseTArray<T>*& out, T mul);
		static void vectorLinearCombination(SparseTArray<T>*& out, const SparseTArray<T>& in1, T mult1, const SparseTArray<T>& in2, T mult2);

		std::ostream & Print(std::ostream & out) const;

		// RefCounter count;

	private:
		size_t len_;
		std::vector<std::pair<size_t, T> > arr_;
	};

	template <typename T>
	struct firstKeyPairVectorComparator
	{
		bool operator()(const std::pair<size_t, T>& a, const std::pair<size_t, T>& b) const
		{
			return a.first < b.first;
		}
	};
}


#include "sparse_tarray_impl.h"
#endif

