#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
// #include "../../../ref_ptr.h"
#include "sparse_tarray.h"

	namespace HowellMatrix
	{
		//Abstract class inteface for Matrix, used by HowellMatrix.
		template <typename T>
		class Matrix {

		public:
			virtual ~Matrix() {}
			virtual bool operator==(const Matrix &t) const = 0;
			virtual bool operator!=(const Matrix &t) const
			{
				return !(*this == t);
			}

			//Getter functions
			virtual size_t NumRows() const = 0;
			virtual size_t NumCols() const = 0;
			virtual T Get(size_t r, size_t c) const = 0;

			//Performance warning: This function call might be costly
			virtual const TArray<T>* GetRow(size_t r) const = 0;

			// virtual const std::unique_ptr<SparseTArray<T> > GetConstSparseRow(size_t r) const = 0;

			virtual SparseTArray<T>* GetSparseRow(size_t r) = 0;

			//Setter Function
			virtual void set(size_t r, size_t c, T val) = 0;

			//Modifier function
			virtual void clear() = 0;

			//Row additions: Unless specified, rows are always added at the end.
			virtual void addRow(const TArray<T>& row) = 0;
			virtual void addRow(const SparseTArray<T>& row) = 0;

			//Copy row at index i and add it to the matrix.
			virtual void addRow(size_t index) = 0;

			virtual void removeRows(const std::vector<size_t>& proj) = 0;
			virtual void removeCols(const std::vector<size_t>& proj) = 0;

			//More row manipulations functions
			virtual void vectorLinearCombination(size_t out_index, size_t in_index_1, T mult_1, size_t in_index_2, T mult_2) = 0;

			//Reference counting
			// RefCounter count;
		};
	}

#endif
