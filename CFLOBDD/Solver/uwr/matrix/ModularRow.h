#ifndef MODULAR_ROW_H
#define MODULAR_ROW_H

#include <stdio.h>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include "../../../ref_ptr.h"

namespace CFL_OBDD
{
	namespace HowellMatrix
	{
		using std::vector;

		//-------------------------------------------------
		// TArray: Just a wrapper around T* that does array
		// bound checking.
		//
		// Usage rule: This class has a private method for 
		// exposing the internal array. It is called from
		// ModuleSpace methods for performance reasons.
		// Care should be taken so that such arrays are never
		// returned outside the scope of ModuleSpace methods. 
		//-------------------------------------------------
		template <typename T>
		class TArray {

		public:
			TArray(size_t l);
			TArray(const TArray &ia);  // Copy constructor.
			TArray(const T* ia, size_t length); // Copies the contents of memory at ia.
			TArray operator=(const TArray &ia);
			~TArray();

			void DeallocateMemory(){  }

			bool operator==(const TArray &t) const;

			size_t getLength() const { return len; }
			const T * getConstArray() const { return m; }

			void set(size_t index, T val);
			T get(size_t index) const;
			size_t getLeadingIndex() const;
			T & operator [](size_t j) const;

			// setRange: Set m[begin] through m[begin + length - 1] to value val.
			// TODO: This should be named setVal. But this was here before I was, 
			// and I'll probably need to change the name at call sites, too.
			void setRange(size_t begin, size_t length, T val);

			// setRange: Set m[begin] through m[begin + length - 1] to values
			// valArray[val_begin] to valArray[val_begin + length - 1].
			void setRange(size_t begin, size_t length,
				const T* valArray, size_t val_begin);

			// insertVals: Insert `length` copies of `val`, starting at index `begin`.
			// This inserts values, and so increases the size of the array.
			void insertVals(const size_t begin, const size_t length = 1,
				const T val = 0);

			// insertRange: Insert a copy of T[val_begin:val_begin+length]
			// at index `begin` of `m`.
			// This inserts values, and so increases the size of the array.
			void insertRange(const size_t begin, const size_t length,
				const T* valArray, size_t val_begin = 0);
			void insertRange(const size_t begin, const size_t length,
				const TArray & valArray, size_t val_begin = 0);

			ref_ptr<TArray<T> > projectOut(const std::vector<size_t> & projectOutIndices) const;

			std::ostream & print(std::ostream & out) const;

			RefCounter count;
		private:
			template <typename _T> friend class ModularRow;
			template <typename _T> friend class ModularSquareMatrix;
			template <typename _T> friend class Mos;
			template <typename _T> friend class HowellMatrix;
			friend class ref_ptr<TArray<T> >;

			size_t len;
			T * m;

			// For use by Mos and ModularRow
			T * getArray() { return m; }
		};

		//Helper Functions
		template <typename T>
		void multiplyMatrices(ref_ptr<TArray<T> > multResult, const ref_ptr<TArray<T> > op1, const ref_ptr<TArray<T> > op2, unsigned N);


		/*
		*  ModularRow
		*
		*  It is used to represent a row in HowellMatrix.
		*  It is a wrapper over TArray but also maintains leading index
		*  and leading rank.
		*/

		template <typename T>
		class ModularRow  {
		public:
			//Constructors
			ModularRow(size_t l);
			ModularRow(const ModularRow &a);
			ModularRow(const T* a, size_t len);
			ModularRow(const vector<T> &);
			ModularRow(const ref_ptr<TArray<T> > & _m, size_t n);
			ModularRow(const ref_ptr<TArray<T> > & _m, size_t _li, size_t _lr);
			virtual ModularRow<T> & operator=(const ModularRow &a);
			virtual ~ModularRow();

			void DeallocateMemory(){ ~ModularRow(); }

			T * getArray();
			const T* getConstArray() const;
			const T * getConstMatrix() const;
			ref_ptr<TArray<T> > getTArray();
			ref_ptr<TArray<T> > getConstTArray() const;

			void setMat(ref_ptr<TArray<T> > mat);

			size_t getSize() const;
			size_t getLeadingIndex() const;
			void setLeadingIndex(size_t leadingIndex);
			size_t getLeadingRank() const;
			void setLeadingRank(size_t leadingRank);

			virtual bool operator== (const ModularRow & mat) const;


			static ref_ptr<ModularRow> copy(const ref_ptr<ModularRow> a);
			static std::vector<ModularRow<T> > copy(const std::vector<ModularRow<T> > & v);
			static std::vector<ref_ptr<ModularRow<T> > > copy(const std::vector<ref_ptr<ModularRow<T> > > & v);

			T get(size_t index) const;

			// set and insert functions
			// ========================
			// The set* functions leave the size of the array unchanged and
			// overwrite array values. In contrast, the insert* functions increase
			// the size of the array but do not overwrite any values.

			void set(size_t index, T val);

			// setVal: Set every value in m[begin:begin+length] to `val`.
			void setVal(size_t begin, size_t length = 1, const T val = 0);
			// setRange: For i with 0 <= i < length, m[i] = val_array[val_begin + i].
			void setRange(size_t begin, size_t length,
				const ModularRow & val_array, size_t val_begin = 0);

			// insertVal: insert `length` copies of `val` at m[begin].
			void insertVal(const size_t begin, const size_t length = 1, const T val = 0);
			// insertRange: insert the array val_array[val_begin:val_begin+length]
			// into this array at m[begin].
			void insertRange(const size_t begin, const size_t length,
				const ModularRow & val_array, size_t val_begin = 0);

			void swapCols(const size_t a, const size_t b);
			void permuteCols(const std::vector<size_t> & from,
				const std::vector<size_t> & to);

			std::ostream & print(std::ostream & out) const;
		private:
			template <typename _T> friend class ModularRow;
			template <typename _T> friend class HowellMatrix;
			template <typename _T> friend class Mos;
			template <typename _T> friend class Ks;

			// The matrix (row major order)
			ref_ptr<TArray<T> > m;

			// Index of the leading non-zero entry
			size_t leading_index;

			// Rank of the leading entry
			size_t leading_rank;

			friend class ref_ptr<ModularRow<T> >;
			RefCounter count;

			// I is hereby deprecate copyIntArray! -elder
			// T * copyIntArray( const T * p, size_t N );

			void init(const T* _m, size_t n);

			// set_leading: fixup leading rank and index, after breaking them.
			// i, if given, is an index no less than the new leading index.
			void set_leading(size_t i = 0);

		};

		template <typename T>
		std::ostream & operator << (std::ostream & out, const ModularRow<T> & row) {
			return row.print(out);
		}

		// Out = in1*mult1 + in2*mult2
		template <typename T>
		void vectorLinearCombination(ModularRow<T>& out,
			ModularRow<T>& in1, T mult1,
			ModularRow<T>& in2, T mult2);
	}
	
} 
#include "ModularRow_impl.h"
#endif //MODULAR_ROW

