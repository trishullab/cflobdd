
#include "ModularSquareMatrix.h"
#include "../bit_vector/bit_vector_ops.h"
#include "../assert/uw_assert.h"

#ifndef __MODULAR_SQUARE_MATRIX_IMPL_H
#define __MODULAR_SQUARE_MATRIX_IMPL_H


namespace HowellMatrix
{
	// ========================================
	// Row and column operations on matrices.
	// ========================================

	// col_swap: in-place column swap.
	template <typename T>
	void ModularSquareMatrix<T>::col_swap(unsigned pos, unsigned col1, unsigned col2) {
		if (col1 == col2) return;
		for (unsigned k = pos; k < dim; k++) {
			T tmp = at(k, col1);
			at(k, col1) = at(k, col2);
			at(k, col2) = tmp;
		}
	}

	// col_mult: in-place column multiplication. col = d * col
	template <typename T>
	void ModularSquareMatrix<T>::col_mult(unsigned pos, unsigned col, T d) {
		if (d == 1) return;
		for (unsigned k = pos; k < dim; k++) {
			at(k, col) *= d;
		}
	}

	// col_sub: in-place column operation: col1 = col1 - x*col2
	template <typename T>
	void ModularSquareMatrix<T>::col_sub(unsigned pos, unsigned col1, unsigned col2, T mult) {
		if (mult == 0) return;
		for (unsigned k = pos; k < dim; k++) {
			at(k, col1) -= mult * at(k, col2);
		}
	}

	// row_swap: in-place row swap.
	template <typename T>
	void ModularSquareMatrix<T>::row_swap(unsigned pos, unsigned row1, unsigned row2) {
		if (row1 == row2) return;
		for (unsigned k = pos; k < dim; k++) {
			T tmp = at(row1, k);
			at(row1, k) = at(row2, k);
			at(row2, k) = tmp;
		}
	}

	// row_mult: in-place row multiplication. row = d * row
	template <typename T>
	void ModularSquareMatrix<T>::row_mult(unsigned pos, unsigned row, T d) {
		if (d == 1) return;
		for (unsigned k = pos; k < dim; k++) {
			at(row, k) *= d;
		}
	}

	// row_sub: in-place row operation: row1 = row1 - x*row2
	template <typename T>
	void ModularSquareMatrix<T>::row_sub(unsigned pos, unsigned row1, unsigned row2, T mult) {
		if (mult == 0) return;
		for (unsigned k = pos; k < dim; k++) {
			at(row1, k) -= mult * at(row2, k);
		}
	}

#ifndef NORMAL_PIVOT_FINDING
	// is_smallest_row: part of find_pivot
	//
	// Determine whether (i,j) is a minimal-rank element in A(i,pos:N-1).
	// If yes, return true;
	// otherwise, set j to the column position of a minimal-rank element in row i, and return false.
	//
	// precondition: pos <= i,j < N
	template <typename T>
	bool is_smallest_row(ModularSquareMatrix<T>*& A,
		size_t pos,
		size_t i,
		size_t &j) {
		size_t minpos = j;
		size_t min = BitVector::compute_rank(A->at(i, j));

		for (size_t k = pos; k < A->getDim(); k++) {
			size_t r = BitVector::compute_rank(A->at(i, k));
			if (r < min) {
				min = r;
				minpos = k;
			}
		}
		if (minpos == j) return true;
		j = minpos;
		return false;
	}

	// is_smallest_col: part of find_pivot.
	//
	// Determine whether (i,j) is a minimal-rank element in A(pos:N-1,j).
	// If yes, return true;
	// otherwise, set i to the row position of a minimal-rank element in column j, and return false.
	//
	// precondition: pos <= i,j < N
	template <typename T>
	bool is_smallest_col(ModularSquareMatrix<T>*& A,
		size_t pos, size_t &i, size_t j) {
		size_t min = BitVector::compute_rank(A->at(i, j));
		size_t minpos = i;
		for (size_t k = pos; k < A->getDim(); k++) {
			size_t r = BitVector::compute_rank(A->at(k, j));
			if (r < min) {
				min = r;
				minpos = k;
			}
		}
		if (minpos == i) return true;
		i = minpos;
		return false;
	}

	//
	// find_pivot
	//
	// Identify the position of an element in A[pos:N-1, pos:N-1] that is of minimal rank
	// in both its row and column (to use, e.g., as a pivot in modular_svd)
	//
	template <typename T>
	void find_pivot(ModularSquareMatrix<T>*&A,
		size_t pos,
		size_t &i, //output
		size_t &j  //output
		) {
		UWAssert::UWAssert::shouldNeverHappen(pos >= A->getDim());

		// Repeatedly, set (i,j) to index the first minimal-rank element in current row,
		//       then, set (i,j) to index the first minimal-rank element in current col.
		//
		// Stop at a fixed point; that will be a min-rank element in its row and column.
		i = pos; j = pos;

		is_smallest_row(A, pos, i, j);
		while (true) {
			if (is_smallest_col(A, pos, i, j)) break;
			if (is_smallest_row(A, pos, i, j)) break;
		}
	}
#else
	//
	// find_pivot
	//
	// Identify the position of an element in A[pos:N-1, pos:N-1] that is of minimal rank
	// over all of the elements in A[pos:N-1, pos:N-1] (to use, e.g., as a pivot in modular_svd).
	//
	// Note: this procedure is over-kill because modular_svd only needs the position of an
	// element in A[pos:N-1, pos:N-1] that is of minimal rank in both its own row and column,
	// not one that is globally minimal in A[pos:N-1, pos:N-1].
	//
	template <typename T>
	void find_pivot(ModularSquareMatrix<T>*A,
		size_t pos,
		size_t &pivot_posi, //output
		size_t &pivot_posj //output 
		) {

		assert(pos < A->getDim());

		// Set (pivot_posi, pivot_posj) to the coordinates of a minimal-rank
		// element in A[pos..(N-1), pos..(N-1)]

		size_t min = highest_power<T>(); // min will be the min rank
		pivot_posi = pos;
		pivot_posj = pos;
		for (size_t i = pos; i<A->getDim(); i++) {
			for (size_t j = pos; j<N; j++) {
				size_t r = BitVector::compute_rank(A->at(i, j));
				if (r < min) { // new minimum-rank element found
					min = r;
					pivot_posi = i;
					pivot_posj = j;
				}
			}
		}
		return;
	}
#endif

	//ModularSquareMatrix implementation.

	//Constructors
	template <typename T>
	ModularSquareMatrix<T>::ModularSquareMatrix(unsigned d)
		: dim(d), m(new TArray<T>(d*d)) {}

	template <typename T>
	ModularSquareMatrix<T>::ModularSquareMatrix(const ModularSquareMatrix &a)
		: dim(a.dim) {
			auto t = a.m;
			SetM(t);
		}

	template <typename T>
	void ModularSquareMatrix<T>::SetM(TArray<T>*& am){
		m = am;
	}

	template <typename T>
	ModularSquareMatrix<T>::ModularSquareMatrix(ModularSquareMatrix* a)
		: dim(a->dim) {
			SetM(a->m);
		}

	// DEPRECATED. commented out to find uses.
	// template <typename T>
	// ModularSquareMatrix<T>::ModularSquareMatrix(const T* _m, unsigned d) 
	//     : m(new TArray<T>(_m, d*d)), dim(d) { }

	template <typename T>
	ModularSquareMatrix<T>::ModularSquareMatrix(TArray<T>* & _m, unsigned d)
		: dim(d) {
		UWAssert::UWAssert::shouldNeverHappen(_m->getLength() != d*d);
		TArray<T>* t = new TArray<T>(_m->getArray(), _m->getLength());
		m = t;
	}


	//Destructor
	template <typename T>
	ModularSquareMatrix<T>::~ModularSquareMatrix() { }

	//Operators
	template <typename T>
	ModularSquareMatrix<T> ModularSquareMatrix<T>::operator =(const ModularSquareMatrix<T> &a) {
		if (this != &a) {
			dim = a.dim;
			m = new TArray<T>(a.getConstMatrix(), dim*dim);
		}
		return *this;
	}

	template <typename T>
	bool ModularSquareMatrix<T>::operator == (const ModularSquareMatrix & mat) const {
		if (dim != mat.getDim()) return false;
		for (unsigned i = 0; i < dim*dim; i++) {
			if (m[i] != mat.m[i]) return false;
		}
		return true;
	}

	// Getters
	template <typename T>
	T * ModularSquareMatrix<T>::getMatrix() {
		return m->getArray();
	}

	template <typename T>
	const T * ModularSquareMatrix<T>::getConstMatrix() const {
		return m->getConstArray();
	}

	template <typename T>
	TArray<T>*& ModularSquareMatrix<T>::getArray() {
		return m;
	}

	template <typename T>
	const TArray<T>* ModularSquareMatrix<T>::getConstArray() const {
		return m;
	}

	template <typename T>
	unsigned ModularSquareMatrix<T>::getDim() const {
		return dim;
	}

	template <typename T>
	T & ModularSquareMatrix<T>::at(unsigned row, unsigned col) const {
		return (*m)[row*dim + col];
	}

	template <typename T>
	void ModularSquareMatrix<T>::set(unsigned r, unsigned c, T val) {
		UWAssert::UWAssert::shouldNeverHappen(r >= dim);
		UWAssert::UWAssert::shouldNeverHappen(c >= dim);
		m->set(r*dim + c, val);
	}

	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::copy(ModularSquareMatrix<T>* & a)
	{
		return new ModularSquareMatrix<T>(a->getArray(), a->getDim());
	}

	//Basic matrix operations.

	//Return a transposed matrix.
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::transpose(ModularSquareMatrix<T>* & mat)
	{
		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(mat->getDim());

		for (unsigned i = 0; i < mat->getDim(); i++) {
			for (unsigned j = 0; j < mat->getDim(); j++) {
				result->set(i, j, mat->at(j, i));
			}
		}
		return result;
	}

	//-------------------
	// Add two matrices
	//-------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::add(
		const ModularSquareMatrix<T>* & mat1,
		const ModularSquareMatrix<T>* & mat2) {
		UWAssert::UWAssert::shouldNeverHappen(mat1->getDim() != mat2->getDim());

		unsigned N = mat1->getDim();
		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);

		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				result->set(i, j, mat1->at(i, j) + mat2->at(i, j));
			}
		}
		return result;
	}

	//-------------------
	// Add two matrices
	//-------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::add_with_check(
		const ModularSquareMatrix<T>* & mat1,
		const ModularSquareMatrix<T>* & mat2)
	{
		UWAssert::UWAssert::shouldNeverHappen(mat1->getDim() != mat2->getDim());
		unsigned N = mat1->getDim();

		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				if (i == 0 && j == 0 && mat1->at(0, 0) == 1 && mat2->at(0, 0) == 1)
					result->set(0, 0, 1);
				else
					result->set(i, j, mat1->at(i, j) + mat2->at(i, j));
			}
		}

		return result;
	}

	//-------------------
	// Combine two matrices: 0 ~ n from mat1 and n+1 ~ N-1 from mat2
	// e.g., when combining regs and mem for ConcLARA or LARA, n is the number of registers
	//-------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::combine_with_check(
		const ModularSquareMatrix<T>* & mat1,
		const ModularSquareMatrix<T>* & mat2,
		unsigned n) {
		UWAssert::UWAssert::shouldNeverHappen(mat1->getDim() != mat2->getDim());
		unsigned N = mat1->getDim();
		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);

		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				if (i == 0 && j == 0 && mat1->at(0, 0) == 1 && mat2->at(0, 0) == 1)
					result->set(0, 0, 1);
				else {
					if (j <= n)
						result->set(i, j, mat1->at(i, j));
					else
						result->set(i, j, mat2->at(i, j));
				}
			}
		}
		return result;
	}

	//------------------------------------
	// Find the difference of two matrices
	//------------------------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::subtract(
		const ModularSquareMatrix<T>* & mat1,
		const ModularSquareMatrix<T>* & mat2)
	{
		UWAssert::UWAssert::shouldNeverHappen(mat1->getDim() != mat2->getDim());
		unsigned N = mat1->getDim();
		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);

		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				result->set(i, j, mat1->at(i, j) - mat2->at(i, j));
			}
		}
		return result;
	}

	//------------------------
	// Multiply two matrices
	//------------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::multiply(
		ModularSquareMatrix<T>* & mat1,
		ModularSquareMatrix<T>* & mat2)
	{
		UWAssert::UWAssert::shouldNeverHappen(mat1->getDim() != mat2->getDim());
		unsigned N = mat1->getDim();

		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);
		for (unsigned i = 0; i < N; ++i) {
			for (unsigned j = 0; j < N; ++j) {
				T sum = 0;
				for (unsigned k = 0; k < N; ++k) {
					sum += mat1->at(i, k) * mat2->at(k, j);
				}
				result->set(i, j, sum);
			}
		}
		return result;
	}

	// Added by Junghee (06/07/11)
	// Rewritten by Matt Elder (18 Sept 2012)
	// Returns col1 * col2^T
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::createNbyNMatrix(
		unsigned N,
		const std::vector<T> & col1,
		const std::vector<T> & col2)
	{
		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);
		for (unsigned i = 0; i<N; i++) {
			for (unsigned j = 0; j<N; j++) {
				result->set(i, j, col1.at(i) * col2.at(j));
			}
		}
		return result;
	}

	//------------------------------------
	// Return the identity matrix
	//------------------------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::mkId(unsigned N) {
		if (N == 0) return NULL;
		ModularSquareMatrix<T>* result(new ModularSquareMatrix<T>(N));
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				result->set(i, j, (i == j) ? 1 : 0);
			}
		}

		return result;
	}


	//------------------------------------
	// Return a matrix with 1's in diagonal from ith column to jth column
	// (0,0) is always 1
	//------------------------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::mkPartialId(unsigned N, unsigned start, unsigned stop) {
		UWAssert::UWAssert::shouldNeverHappen(!(start <= stop && stop < N));

		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				result->set(i, j, (i == j
					&& (i == 0
					|| (start <= i && i <= stop))) ? 1 : 0);

			}
		}

		return result;
	}

	//------------------------------------
	// Return Id, except for a 0 in one on-diagonal: e_i,i
	//------------------------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::mkSinglePointProjectionMatrix(
		unsigned N, unsigned omit) {

		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				result->set(i, j, (i == j && i != omit) ? 1 : 0);
			}
		}

		return result;
	}

	//------------------------------------
	// a 1 in some e_i,j; 0 everywhere else
	//------------------------------------
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::mkSinglePointMatrix(
		unsigned N, unsigned row, unsigned col) {

		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(N);
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < N; j++) {
				result->set(i, j, 0);
			}
		}
		result->set(row, col, 1);
		return result;
	}

	/* ----------------------------------------------------
	* inflateTo2K converts (n+1) x (n+1) matrix
	*    to (2k+1) x (2k+1) matrix in the following way
	*
	* -------------           ------------------
	* |j|   b^T   |           |j|  b^T  |  0   |
	* -------------           ------------------
	* | |         |     ==>   | |       |      |
	* |0|   A^T   |           |0|   0   |  0   |
	* | |         |           | |       |      |
	* | |         |           ------------------
	* -------------           | |       |      |
	*                         |0|  A^T  |  jxI |
	*                         | |       |      |
	*                         ------------------
	* junghee@cs.wisc.edu
	* edited to patch mem leaks by elder@cs.
	* ---------------------------------------------------
	*/
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::inflateTo2K(const ModularSquareMatrix<T>* & mat)
	{
		unsigned N = mat->getDim();
		unsigned K = N - 1;
		unsigned inflatedN = K * 2 + 1; // N + K

		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(inflatedN);
		// copy from (0,0) ~ (0,N-1)
		for (unsigned i = 0; i < N; i++) {
			result->set(0, i, mat->at(0, i));
		}
		// set 0 from (0,N) ~ (0,N+K-1)
		for (unsigned i = N; i < inflatedN; i++) {
			result->set(0, i, 0);
		}
		// set 0 from (1,0) ~ (N-1,N+K-1)
		for (unsigned i = 1; i < N; i++) {
			for (unsigned j = 0; j < inflatedN; j++) {
				result->set(i, j, 0);
			}
		}
		// set appropriate values from (N,0) ~ (N+K-1,N+K-1)
		for (unsigned i = N; i < inflatedN; i++) {
			for (unsigned j = 0; j < inflatedN; j++) {
				if (j == 0) {
					result->set(i, j, 0);
				}
				else if (j >= 1 && j < N) {
					result->set(i, j, mat->at(i - K, j));
				}
				else {// j >= N
					result->set(i, j, (i == j) ? mat->at(0, 0) : 0);
				}
			}
		}

		return result;
	}

	// Matrix operations

	// modular_svd
	//
	// Perform a modular-svd decomposition:
	// Find Linv, D, and Rinv such that A = LDR.
	//
	// The procedure destructively modifies A to create D, and produce Rinv and Linv.
	template <typename T>
	void ModularSquareMatrix<T>::modular_svd(ModularSquareMatrix<T>* & A,
		ModularSquareMatrix<T>* & Rinv,
		ModularSquareMatrix<T>* & Linv) {

		UWAssert::UWAssert::shouldNeverHappen(A->getDim() != Rinv->getDim() || Rinv->getDim() != Linv->getDim());

		size_t r;
		T d, dprime, x;

		size_t N = A->getDim(); // set class data member


		// T* A_mat = A.getMatrix();
		// T* Rinv_mat = Rinv.getMatrix();
		// T* Linv_mat = Linv.getMatrix();

		//ModularSquareMatrix<T>* original_A = ModularSquareMatrix<T>::copy(A);

		Rinv = mkId(N);
		Linv = mkId(N);

		// The factoring loop: similar to Gaussian elimination, but zeroes out both
		// (i) the column elements A[k][pos] below the pivot element A[pos][pos], and
		// (ii) the row elements of A[pos][k] to the right of the pivot element
		for (size_t pos = 0; pos < N - 1; pos++) {
			// Find the position of an element in A[pos:N-1, pos:N-1] that is
			// of minimal rank in both its row and column. Perform a row swap
			// and a column swap to bring the element to A[pos,pos].
			size_t i, j;

			find_pivot(A, pos, i, j);
			r = BitVector::compute_rank(A->at(i, j));

			// Bring pivot to diagonal
			A->col_swap(pos, pos, j);
			A->row_swap(pos, pos, i);

			// Do column operations on Rinv, row operations on Linv
			Rinv->col_swap(0, pos, j);
			Linv->row_swap(0, pos, i);

			if (r == BitVector::highest_power<T>()) { // A[pos,pos] is already zero, and hence so are
				continue;                 // all elements in both the pos^th row and pos^th column.
			}

			// Set d to the invertible "core" of the pivot element A[pos,pos].
			// That is, set d so that A[pos,pos] = d * 2^r.
			d = BitVector::div_pow2(A->at(pos, pos), r);

			// Perform row operations on A to zero out the elements of the pos^th column
			// below A[pos,pos]. Collect the row operations in Linv.
			for (unsigned k = pos + 1; k < N; k++) {
				unsigned rprime = BitVector::compute_rank(A->at(k, pos));
				if (rprime == BitVector::highest_power<T>()) continue; // already zero

				UWAssert::UWAssert::shouldNeverHappen(!(r <= rprime));

				T dprime = BitVector::div_pow2(A->at(k, pos), rprime);
				T x = BitVector::mult_pow2(dprime, rprime - r);

				// Perform a row operation on A to zero out A[k][pos]
				A->row_mult(pos, k, d);
				A->row_sub(pos, k, pos, x);

				// Perform the corresponding row operation on Linv
				Linv->row_mult(0, k, d);
				Linv->row_sub(0, k, pos, x);

				UWAssert::UWAssert::shouldNeverHappen(A->at(k, pos) != 0);
			}

			// Perform column operations on A and Rinv. The column operations on A are
			// simplified because the pos^th column [0, ..., 0, 2^r * d, 0, 0 ..., 0]^t
			// has just one non-zero element, namely, at the pivot position A[pos,pos].
			// Thus, we just zero out the elements of the pos^th row of A -- i.e.,
			// A[pos][k], for pos+1 <= k < N.
			// However, we still need to perform full column operations on Rinv.
			for (unsigned k = pos + 1; k < N; k++) {
				// Check whether we need to change Rinv
				unsigned rprime = BitVector::compute_rank(A->at(pos, k));
				if (rprime == BitVector::highest_power<T>()) continue; // already zero; no need to change Rinv

				// Perform a column operation on A to zero out A[pos][k]
				UWAssert::UWAssert::shouldNeverHappen(!(r <= rprime));

				T dprime = BitVector::div_pow2(A->at(pos, k), rprime);
				T x = BitVector::mult_pow2(dprime, rprime - r);

				A->col_mult(pos, k, d);
				A->set(pos, k, 0);   // Optimized col_sub(A, pos, k, pos, x);

				// Perform the corresponding column operation on Rinv
				Rinv->col_mult(0, k, d);
				Rinv->col_sub(0, k, pos, x);
			}
		} // factoring loop
	}

	// replace_column
	//
	//    Replace the column of the given index with the given col in the matrix
	//
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::replace_column(
		const ModularSquareMatrix<T>* & mat, unsigned index, const std::vector<T> & col
		) {
		UWAssert::UWAssert::shouldNeverHappen(index >= (int)col.size());

		unsigned N = mat->getDim();

		TArray<T>* m = new TArray<T>(mat->getConstMatrix(), N*N);
		ModularSquareMatrix<T>* result = new ModularSquareMatrix<T>(m, N);

		for (unsigned i = 0; i < N; i++) {
			result->set(i, index, col[i]);
		}
		return result;
	}

	// get_column_of_index
	// 
	//   Obtain the column of the given index in the matrix
	//
	template <typename T>
	std::vector<T> ModularSquareMatrix<T>::get_column_of_index(const ModularSquareMatrix<T>* & mat, int index)
	{
		std::vector<T> ans;

		for (unsigned i = 0; i < mat->getDim(); i++) {
			ans.push_back(mat->at(i, index));
		}

		return ans;
	}

	/*
	* dualize
	*
	* Do the dualization(perp operator) of the given matrix.
	* dualize(M) = (L^{-1})^t T (R^{-1})^t
	* Where L D R = M is the diagonal decomposition of M,
	* and T_{i,i} = 2^{w-rank(D_{i,i})}
	*/
	template <typename T>
	ModularSquareMatrix<T>* ModularSquareMatrix<T>::dualize(const ModularSquareMatrix<T> & M)
	{
		ModularSquareMatrix<T>* Linv = new ModularSquareMatrix<T>(M.getDim());
		ModularSquareMatrix<T>* Rinv = new ModularSquareMatrix<T>(M.getDim());
		ModularSquareMatrix<T>* D = new ModularSquareMatrix<T>(M);

		//Obtain D, Rinv, Linv by using modular_svd method.
		modular_svd(D, Rinv, Linv);

		ModularSquareMatrix<T>* Rinv_transpose = transpose(Rinv);
		ModularSquareMatrix<T>* Linv_transpose = transpose(Linv);

		//Obtain T from D
		for (size_t i = 0; i < D->getDim(); i++) {
			D->set(i, i, BitVector::mult_pow2<T>(1, BitVector::highest_power<T>() - BitVector::compute_rank(D->at(i, i))));
		}
		auto tmp = multiply(D, Rinv_transpose);
		return multiply(Linv_transpose, tmp);
	}


	template <typename T>
	std::ostream & ModularSquareMatrix<T>::print(std::ostream & out) const
	{
		out << "[" << std::endl;
		for (unsigned i = 0; i < getDim(); i++) {
			out << "  ";
			for (unsigned j = 0; j < getDim(); j++) {
				out << this->at(i, j) << ", ";
			}
			out << std::endl;
		}
		out << "]";
		return out;
	}

	//
	// validate
	//
	// Validate that the first entry of the first matrix is one.
	template <typename T>
	void ModularSquareMatrix<T>::validate() const {
		UWAssert::UWAssert::shouldNeverHappen(this->getConstArray()[0] != 1);
	}
}

#endif