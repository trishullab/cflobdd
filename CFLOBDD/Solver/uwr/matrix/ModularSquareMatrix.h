#ifndef MODULAR_SQUARE_MATRIX_H
#define MODULAR_SQUARE_MATRIX_H

#include <stdio.h>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include "../../../ref_ptr.h"
#include "../assert/uw_assert.h"
#include "ModularRow.h"
#include "HowellMatrix.h"

#define GEN_MODULE_DBGS(stmts) do{ if(false/*debug_module_space*/) { stmts; } }while(0)

namespace CFL_OBDD
{
	namespace HowellMatrix
	{
		/*
		*  ModularSquareMatrix
		*
		*  This is class used to represent Square Matrices over modular integers.
		*  It provides normal matrix operations like transpose and also a few complex operations
		*  like modular svd.
		*
		*/
		template <typename _T>
		class ModularSquareMatrix {

		private:
			//The dimension of the square matrix
			unsigned dim;
			//The matrix
			ref_ptr<TArray<_T> > m;

		public:
			RefCounter count;

			//Constructors
			ModularSquareMatrix(unsigned dimension);
			ModularSquareMatrix(const ModularSquareMatrix &a);
			// DEPRECATED ModularSquareMatrix(const _T* _m, unsigned dimension);
			ModularSquareMatrix(const ref_ptr<TArray<_T> > & _m, unsigned dimension);

			//Destructor
			~ModularSquareMatrix();
			void DeallocateMemory(){  }
			//Operators
			ModularSquareMatrix operator =(const ModularSquareMatrix &a);
			bool operator == (const ModularSquareMatrix & mat) const;

			// Getters



			// TODO: DEPRECATE
			// Deprecate these. Use inlined accessors instead; far cleaner interface.
			// I'm uneasy about these. getDim() is just fine, but we should eliminate the rest.
			// -- elder
			_T * getMatrix();
			const _T * getConstMatrix() const;
			ref_ptr<TArray<_T> > getArray();
			const ref_ptr<TArray<_T> > getConstArray() const;

			unsigned getDim() const;

			//Setter
			void set(unsigned r, unsigned col, _T val);
			// TODO: DEPRECATE
			// Why is this different from the copy constructor? Be suspicious everywhere you see this.
			// -- elder
			static ref_ptr<ModularSquareMatrix> copy(const ref_ptr<ModularSquareMatrix> & a);

			//Basic matrix operations.
			static ref_ptr<ModularSquareMatrix> transpose(const ref_ptr<ModularSquareMatrix> & mat);
			static ref_ptr<ModularSquareMatrix> add(const ref_ptr<ModularSquareMatrix<_T> > & mat1, const ref_ptr<ModularSquareMatrix<_T> > & mat2);
			static ref_ptr<ModularSquareMatrix> add_with_check(const ref_ptr<ModularSquareMatrix<_T> > & mat1, const ref_ptr<ModularSquareMatrix<_T> > & mat2);
			static ref_ptr<ModularSquareMatrix> combine_with_check(const ref_ptr<ModularSquareMatrix<_T> > & mat1, const ref_ptr<ModularSquareMatrix<_T> > & mat2, unsigned n);
			static ref_ptr<ModularSquareMatrix> subtract(const ref_ptr<ModularSquareMatrix<_T> > & mat1, const ref_ptr<ModularSquareMatrix<_T> > & mat2);
			static ref_ptr<ModularSquareMatrix> multiply(const ref_ptr<ModularSquareMatrix<_T> > & mat1, const ref_ptr<ModularSquareMatrix<_T> > & mat2);
			static ref_ptr<ModularSquareMatrix<_T> > createNbyNMatrix(unsigned N,
				const std::vector<_T> & col1, const std::vector<_T> & col2);
			static ref_ptr<ModularSquareMatrix<_T> > mkId(unsigned N);
			static ref_ptr<ModularSquareMatrix<_T> > mkPartialId(unsigned N, unsigned i, unsigned j);
			static ref_ptr<ModularSquareMatrix<_T> > mkSinglePointProjectionMatrix(unsigned N, unsigned i);
			static ref_ptr<ModularSquareMatrix<_T> > mkSinglePointMatrix(unsigned N, unsigned i, unsigned j);
			static ref_ptr<ModularSquareMatrix<_T> > inflateTo2K(const ref_ptr<ModularSquareMatrix<_T> > & mat);
			static void modular_svd(ref_ptr<ModularSquareMatrix<_T> > & A,    // input, modified to output.
				ref_ptr<ModularSquareMatrix<_T> > & Rinv, // output
				ref_ptr<ModularSquareMatrix<_T> > & Linv  // output
				);

			static ref_ptr<ModularSquareMatrix<_T> > replace_column(const ref_ptr<ModularSquareMatrix<_T> > & mat, unsigned index, const std::vector<_T> & col);
			static std::vector<_T> get_column_of_index(const ref_ptr<ModularSquareMatrix<_T> > & mat, int index);

			static ref_ptr<ModularSquareMatrix<_T> > dualize(const ModularSquareMatrix<_T> & M); //Dualization operator, for interconversion between AG and KS domain.
			std::ostream & print(std::ostream & out) const;
			void validate() const;

			_T & at(unsigned row, unsigned col) const;

		private:
			// All column operations (col_*) assume that rows above 'pos' are all zeroes.
			// All row operations (row_*) assume that columns left of 'pos' are all zeroes.
			// In all cases, it's always safe (but maybe inefficient) to use pos == 0).

			// col_swap: Swap columns 'col1' and 'col2'.
			void col_swap(unsigned pos, unsigned col1, unsigned col2);
			// col_mult: Multiply column 'col' by 'mult'.
			void col_mult(unsigned pos, unsigned col, _T mult);
			// col_sub: Column subtaction. Perform col1 -= mult * col2.
			void col_sub(unsigned pos, unsigned col1, unsigned col2, _T mult);

			// row_swap: Swap rowumns 'row1' and 'row2'.
			void row_swap(unsigned pos, unsigned row1, unsigned row2);
			// row_mult: Multiply rowumn 'row' by 'mult'.
			void row_mult(unsigned pos, unsigned row, _T mult);
			// row_sub: Rowumn subtaction. Perform row1 -= mult * row2.
			void row_sub(unsigned pos, unsigned row1, unsigned row2, _T mult);

		};
	}
}
#include "ModularSquareMatrix_impl.h"
#endif //MODULAR_SQUARE_MATRIX_H
