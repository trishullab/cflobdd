#ifndef SPARSE_ROW_MATRIX_H
#define SPARSE_ROW_MATRIX_H

#include <stdio.h>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include "../../../ref_ptr.h"
#include "matrix.h"
#include "sparse_tarray.h"

namespace CFL_OBDD
{
	namespace HowellMatrix
	{
		//Matrix with sparse representation for rows.
		template <typename T>
		class SparseRowMatrix : public Matrix<T> {

		public:
			template <typename BV_Type> friend class HowellMatrix;

			SparseRowMatrix(size_t col_size);
			SparseRowMatrix(const SparseRowMatrix& m);
			~SparseRowMatrix();

			void DeallocateMemory(){}
			SparseRowMatrix<T> operator=(const Matrix<T> &ia);
			bool operator==(const Matrix<T> &t) const;

			//Getter functions
			size_t NumRows() const;
			size_t NumCols() const;
			T Get(size_t r, size_t c) const;
			size_t GetLeadingIndex(size_t r) const;
			size_t GetLeadingRank(size_t r) const;

			//Performance warning: This function call might be costly
			const ref_ptr<TArray<T> > GetRow(size_t r) const;

			const ref_ptr<SparseTArray<T> > GetConstSparseRow(size_t r) const;
			ref_ptr<SparseTArray<T> > GetSparseRow(size_t r);
			const ref_ptr<SparseTArray<T> > operator [](size_t j) const;

			//Setter Function
			void set(size_t r, size_t c, T val);

			//Modifier function
			void clear();

			//Row additions: Unless specified, rows are always added at the end.
			void addRow(const T* arr);
			void addRow(const std::vector<T> & vec);
			void addRow(const TArray<T>& row);
			void addRow(const SparseTArray<T>& row);
			void addRow(ref_ptr<SparseTArray<T> >& row);
			//Copy row at index i and add it to the matrix.
			void addRow(size_t index);

			//Row removal
			void removeRows(const std::vector<size_t>& proj);

			//Column operations
			void addCols(size_t begin, size_t length);
			void removeCols(const std::vector<size_t>& proj);
			void permuteCols(const std::vector<size_t> & from,
				const std::vector<size_t> & to);
			void insertVals(size_t begin, size_t len,
				const ref_ptr<SparseRowMatrix<T> > & vals);

			//Put the columns in from at the end.
			void moveEnd(const std::vector<size_t> & m);

			//More row manipulations functions
			void vectorMul(size_t out_index, T mult);
			void vectorLinearCombination(size_t out_index, size_t in_index_1, T mult_1, size_t in_index_2, T mult_2);

			//Print functions
			std::ostream & Print(std::ostream & out) const;

			static SparseRowMatrix<T> parse(std::istream & buf);

		private:
			typedef std::vector<ref_ptr<SparseTArray<T> > > RowVector;
			RowVector mat_;
			size_t num_cols_;
		};
	}
}
#include "sparse_row_matrix_impl.h"
#endif
