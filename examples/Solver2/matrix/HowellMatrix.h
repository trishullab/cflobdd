#ifndef HOWELL_MATRIX_H
#define HOWELL_MATRIX_H

#include <algorithm>
#include <sstream>
#include "ModularRow.h"
#include "ModularSquareMatrix.h"
#include "sparse_row_matrix.h"
#include "../bit_vector/value.h"


namespace HowellMatrix
{
	// Generic HowellMatrix defines the interface that HowellMatrix must follow
	// This is primarily here to support matrix member is KsAv which has to be a
	// GenericHowellMatrix as it's size can be 64, 32, 16 or 8 bits.
	class GenericHowellMatrix {
	public:
		// Reference-count field needed to support use of std::unique_ptr<.>
		// RefCounter count;

		// constructor
		GenericHowellMatrix()         {}
		virtual ~GenericHowellMatrix() {}
		virtual void DeallocateMemory() const { }
		virtual GenericHowellMatrix* Copy() const = 0;
		virtual bool isEqual(const GenericHowellMatrix & op2) const = 0;
		virtual bool isEqual(const GenericHowellMatrix* & op2) const = 0;

		// Check if r1 is a logical consequence of r2
		virtual bool IsLogicalConsequence(size_t r1, size_t r2) const = 0;

		virtual size_t GetVectorSize() const = 0;
		virtual size_t GetLeadingIndex(size_t r) const = 0;
		virtual size_t GetLeadingRank(size_t r) const = 0;

		virtual std::vector<std::pair<size_t, BitVector::Value> > getRowAsValueCoeffs(size_t) const = 0;

		//Clears all the constraints. Used to get the top KS element.
		virtual void clear() = 0;

		// Output
		virtual std::ostream & print(std::ostream & out) const = 0;
		virtual std::string str() const = 0;

		virtual void reduceToHowellForm() = 0;

		virtual void setHowellized() = 0;

		virtual void negate() = 0;

		// Insert 'num' 0-columns at index 'index'. (index=0 => add to the left.)
		virtual void insertCols(const size_t index, const size_t num) = 0;

		// Permute columns. For each i, send column from[i] to to[i].
		virtual void permuteCols(const std::vector<size_t> & from,
			const std::vector<size_t> & to,
			bool perform_howellize = true) = 0;

		// Concatenate matrices side-by-side. Produce block matrix [this that].
		virtual void concatRight(const GenericHowellMatrix & that) = 0;
		void concatRight(const GenericHowellMatrix* that) { concatRight(*that); }

		// Concatenate matrices vertically. Produce block matrix:
		// [ this ]
		// [ that ]
		virtual void concatDown(const GenericHowellMatrix & that) = 0;
		void concatDown(const GenericHowellMatrix* that) { concatDown(*that); }

		// Project away (some) left-hand columns.
		// Really, this is shorthand for project([index,index+1,...,num_vars])
		virtual GenericHowellMatrix* generic_project_left(size_t index) const = 0;
		virtual GenericHowellMatrix* generic_project(const std::vector<size_t> & projectionVariables) const = 0;

		virtual GenericHowellMatrix* generic_havoc(std::vector<size_t> havocSet) const = 0;
		virtual size_t num_cols() const = 0;
		virtual size_t num_rows() const = 0;

		virtual BitVector::Value getEntryAsValue(size_t rowNum, size_t columnNum) const = 0;
		virtual void setEntryAsValue(size_t rowNum, size_t columnNum, BitVector::Value val) = 0;

	};

	/*
	* HowellMatrix
	*
	* This class implements modular-arithmetic matrices that are
	* maintained in Howell form. The class maintains as INVARIANT that
	* 'mat' is in HowellForm.
	*
	* See Definition 1 in
	* Elder, M., Lim, J., Sharma, T., Andersen, T., and Reps, T.,
	* Abstract domains of affine relations. To appear in Proc. Static
	* Analysis Symposium (SAS), 2011.
	*
	* A matrix M is in Howell form iff
	* 1. M is in row-echelon form,
	* 2. the leading value of every row is a power of two,
	* 3. each leading value is the largest value in its column, and
	* 4. for every row r of M, for any p in Z, if i is the leading index of (2^p)r,
	*    then (2^p)r in rowspace([M]_i), where [M]_i denotes the matrix that
	*    consists of all rows of M whose leading index is i or greater.
	*/
	template <typename T>
	class HowellMatrix : public GenericHowellMatrix{
	public:
		// Constructors and destructors
		HowellMatrix(size_t vectorSize, bool keep_zero_vector = false);// Creates an empty matrix
		HowellMatrix(const HowellMatrix & rhs);//Copy constructor
		HowellMatrix(HowellMatrix* rhs);//Copy constructor
		~HowellMatrix();      // Destructor
		void DeallocateMemory() const {}
		HowellMatrix(ModularSquareMatrix<T> ms_matrix);

		// Properties of the affine space
		bool isEqual(const HowellMatrix & op2) const;
		bool isEqual(const HowellMatrix* & op2) const;
		bool isEqual(const GenericHowellMatrix & op2) const;
		bool isEqual(const GenericHowellMatrix* & op2) const;
		bool operator==(const HowellMatrix & that) const;
		HowellMatrix& operator=(const HowellMatrix & that)
		{
			if (this == &that)
				return *this;
			mat_ = new SparseRowMatrix<T>(*(that.mat_));
			howellized_ = that.howellized_;
			return *this;
		}
		GenericHowellMatrix* Copy() const
		{
			return new HowellMatrix<T>(*this);
		}

		// Check if r1 is a logical consequence of r2
		bool IsLogicalConsequence(size_t r1, size_t r2) const;

		size_t GetVectorSize() const;
		size_t GetLeadingIndex(size_t r) const;
		size_t GetLeadingRank(size_t r) const;

		//Clears all the constraints. Used to get the top KS element.
		void clear();

		// Output
		std::ostream & print(std::ostream & out) const;
		std::string str() const;

		// Input
		static HowellMatrix of_str(const std::string & s);
		static HowellMatrix parse(std::istream & in);


		// Insert the set of vectors vectorSet.
		void insertSetOfRows(const std::vector<ModularSquareMatrix<T>* > & vectorSet, bool howellize = true);
		void insertSetOfRows(const std::vector<ModularRow<T>* > & vectorSet, bool howellize = true);
		void insertSetOfRows(const std::vector<ModularRow<T> > & vectorSet, bool howellize = true);
		void insertSetOfRows(const std::vector<TArray<T>* > & m, bool howellize = true);
		void insertSetOfRows(const std::vector<SparseTArray<T> > & m, bool howellize = true);
		void insertSetOfRows(const std::vector<TArray<T> > & m, bool howellize = true);
		void insertSetOfRows(const std::vector<const T *> & f, bool howellize = true);
		void insertSetOfRows(const std::vector<std::vector<T> > & f, bool howellize = true);
		void insertSetOfRows(std::vector<SparseTArray<T>* > & f, bool howellize = true);
		void insertSetOfRows(const std::vector<SparseTArray<T>* > & f, bool howellize = true);

		// Insert the vectorSize vector pointed. Allocates memory for vector
		void insertRow(const std::vector<T> & m, bool perform_howellize = true);
		void insertRow(const ModularRow<T> & m, bool perform_howellize = true);
		void insertRow(const TArray<T>* & m, bool perform_howellize = true);
		void insertRow(const TArray<T>& m, bool perform_howellize = true);
		void insertRow(const T * f, bool perform_howellize = true);
		void insertRow(const SparseTArray<T>& f, bool perform_howellize = true);
		void insertRow(const SparseTArray<T>* f, bool perform_howellize = true);

		// Insert 'num' 0-columns at index 'index'. (index=0 => add to the left.)
		void insertCols(const size_t index, const size_t num = 1);

		// Permute columns. For each i, send column from[i] to to[i].
		void permuteCols(const std::vector<size_t> & from,
			const std::vector<size_t> & to,
			bool perform_howellize = true);

		// Project away (some) left-hand columns.
		// Really, this is shorthand for project([index,index+1,...,num_vars])
		HowellMatrix* project_left(size_t index) const;
		GenericHowellMatrix* generic_project_left(size_t index) const;

		// Concatenate matrices side-by-side. Produce block matrix [this that].
		void concatRight(const GenericHowellMatrix & that);
		void concatRight(const HowellMatrix & that);
		void concatRight(const HowellMatrix* that);
		// Concatenate matrices vertically. Produce block matrix:
		// [ this ]
		// [ that ]
		void concatDown(const GenericHowellMatrix & that);
		void concatDown(const HowellMatrix & that);
		void concatDown(const HowellMatrix* that);

		void negate();

		void reduceToHowellForm();

		//The caller should only call this when it's sure that the matrix is already in howellForm
		void setHowellized() { howellized_ = true; }
		HowellMatrix* havoc(size_t i) const;
		HowellMatrix* havoc(std::vector<size_t> havocSet) const;
		GenericHowellMatrix* generic_havoc(std::vector<size_t> havocSet) const;
		HowellMatrix* project(const std::vector<size_t> & projectionVariables) const;
		GenericHowellMatrix* generic_project(const std::vector<size_t> & projectionVariables) const;
		size_t size() const; //Deprecated: Use num_cols instead.
		size_t num_cols() const;
		size_t num_rows() const;
		const ModularRow<T> * get(size_t) const;
		SparseTArray<T> * getSparseRow(size_t) const;

		std::vector<std::pair<size_t, BitVector::Value> > getRowAsValueCoeffs(size_t) const;
		const std::vector<ModularRow<T>* > getRows(size_t, size_t) const;

		T getEntry(size_t rowNum, size_t columnNum) const;
		BitVector::Value getEntryAsValue(size_t rowNum, size_t columnNum) const;

		// It's the responsibility of the caller to ensure that the
		// matrix is later reverted to HowellForm. The only purpose of
		// this setEntry interface is for the purposes of optimization.
		void setEntry(size_t rowNum, size_t columnNum, T val);
		void setEntryAsValue(size_t rowNum, size_t columnNum, BitVector::Value val);
		ModularSquareMatrix<T> getSquareMatrix() const;


	private:
		std::vector<std::pair<size_t, BitVector::Value> > getRowAsValueCoeffs(SparseTArray<T>* & sp_arr) const;

		//Data members
		SparseRowMatrix<T>* mat_;
		mutable bool howellized_;
	}; // class HowellMatrix

	template <typename T>
	std::ostream & operator << (std::ostream & out, const HowellMatrix<T> & hmat) {
		return hmat.print(out);
	}

	template <typename T>
	std::istream & operator>> (std::istream &, HowellMatrix<T> &);

}

#include "HowellMatrix_impl.h"
#endif  // HOWELL_MATRIX_H

