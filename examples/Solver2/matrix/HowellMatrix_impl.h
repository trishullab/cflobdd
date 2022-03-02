#include <stdio.h>
#include "HowellMatrix.h"
#include "../assert/uw_assert.h"
#include "../parsing/parsing.h"
#include "../bit_vector/bit_vector_ops.h"
#include <algorithm>    

#ifndef __HOWELL_MATRIX_IMPL_H
#define __HOWELL_MATRIX_IMPL_H



namespace HowellMatrix
{
	using std::vector;
	using std::string;

	//Implementation for HowellMatrix constructors and methods

	template <typename T>
	HowellMatrix<T>::HowellMatrix(size_t vecSize, bool)
		: mat_(new SparseRowMatrix<T>(vecSize)),
		howellized_(true)
	{ }
	//---------------------------------------------------
	// Copy constructor
	//---------------------------------------------------
	template <typename T>
	HowellMatrix<T>::HowellMatrix(const HowellMatrix<T> & rhs)
		: mat_(new SparseRowMatrix<T>(*(rhs.mat_))),
		howellized_(rhs.howellized_)
	{
	}

	template <typename T>
	HowellMatrix<T>::HowellMatrix(HowellMatrix<T>* rhs)
		: mat_(new SparseRowMatrix<T>(*(rhs->mat_))),
		howellized_(rhs->howellized_)
	{
	}

	// Construct this from a ModularSquareMatrix.
	template <typename T>
	HowellMatrix<T>::HowellMatrix(ModularSquareMatrix<T> ms_mat)
		: mat_(new SparseRowMatrix<T>(ms_mat.getDim())),
		howellized_(false)
	{
		const T* arr = ms_mat.getConstMatrix();
		for (size_t i = 0; i < ms_mat.getDim(); i++) {
			mat_->addRow(arr + i*GetVectorSize());
			//mat_->insertRow(arr + i*GetVectorSize(), ms_mat.getDim());
		}

		reduceToHowellForm();
	}

	//---------------------------------------------------
	// Destructor 
	//---------------------------------------------------
	template <typename T>
	HowellMatrix<T>::~HowellMatrix() {
		// Nothing needed
	}

	// clear constraints
	template <typename T>
	void HowellMatrix<T>::clear() {
		mat_->clear();
		howellized_ = true;
	}
	/*
	template <typename T>
	void HowellMatrix<T>::DeallocateMemory() const{

	}
	*/
	template <typename T>
	inline bool HowellMatrix<T>::isEqual(const GenericHowellMatrix & that) const {
		const HowellMatrix<T>* that_cast = dynamic_cast<const HowellMatrix<T>* >(&that);
		return isEqual(*that_cast);
	}

	template <typename T>
	inline bool HowellMatrix<T>::isEqual(const GenericHowellMatrix* & that) const {
		const HowellMatrix<T>* that_cast = dynamic_cast<const HowellMatrix<T>* >(that);
		return isEqual(*that_cast);
	}

	template <typename T>
	inline bool HowellMatrix<T>::isEqual(const HowellMatrix<T> & that) const {
		if (howellized_)
		{
			if (that.howellized_)
			{
				return ((this->mat_) == (that.mat_));
			}
			else
			{
				HowellMatrix<T> that_c(that);
				that_c.reduceToHowellForm();
				return ((this->mat_) == (that_c.mat_));
			}
		}
		else
		{
			HowellMatrix<T> this_c(*this);
			this_c.reduceToHowellForm();
			if (that.howellized_)
			{
				return ((this_c.mat_) == (that.mat_));
			}
			else
			{
				HowellMatrix<T> that_c(that);
				that_c.reduceToHowellForm();
				return ((this_c.mat_) == (that_c.mat_));
			}
		}

		return false;
	}

	template <typename T>
	inline bool HowellMatrix<T>::isEqual(const HowellMatrix<T>* & that) const {
		return isEqual(*that);
	}

	template <typename T>
	inline bool HowellMatrix<T>::operator==(const HowellMatrix<T> & that) const {
		return isEqual(that);
	}


	// Insert a row
	template <typename T>
	void HowellMatrix<T>::insertRow(const ModularRow<T> & f, bool perform_howellize) {
		UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != f.getSize());
		mat_->addRow(*(f.getConstTArray()));
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertRow(const TArray<T>& f, bool perform_howellize) {
		UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != f.getLength());
		mat_->addRow(f);
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertRow(const TArray<T>* & f, bool) {
		insertRow(*f);
	}
	template <typename T>
	void HowellMatrix<T>::insertRow(const SparseTArray<T> & f, bool perform_howellize) {
		UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != f.GetLength());
		mat_->addRow(f);
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertRow(const SparseTArray<T>* f, bool perform_howellize) {
		insertRow(*f, perform_howellize);
	}
	template <typename T>
	void HowellMatrix<T>::insertRow(const std::vector<T> & vec, bool perform_howellize) {
		UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vec.size());
		mat_->addRow(vec);
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertRow(const T * f, bool perform_howellize) {
		mat_->addRow(f);
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}

	template <typename T>
	void HowellMatrix<T>::insertCols(const size_t index, const size_t num) {
		mat_->addCols(index, num);
	}

	template <typename T>
	void HowellMatrix<T>::permuteCols(const std::vector<size_t> & from,
		const std::vector<size_t> & to,
		bool perform_howellize) {
		mat_->permuteCols(from, to);
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}

	template <typename T>
	GenericHowellMatrix* HowellMatrix<T>::generic_project_left(size_t index) const {
		return project_left(index);
	}

	// Project away (some) left-hand columns.
	// Really, this is shorthand for project([index,index+1,...,num_vars])
	template <typename T>
	HowellMatrix<T>* HowellMatrix<T>::project_left(size_t index) const {
		UWAssert::UWAssert::shouldNeverHappen(!howellized_);
		std::vector<SparseTArray<T>* > new_mat;
		for (size_t i = 0; i < size(); i++)
		{
			if (mat_->GetLeadingIndex(i) >= index)
			{
				SparseTArray<T>* projected_row = new SparseTArray<T>(*(mat_->GetSparseRow(i)));
				projected_row->projectRight(index);
				new_mat.push_back(projected_row);
			}
		}

		//Projected rows are already in howellForm
		HowellMatrix<T>* ret = new HowellMatrix<T>(GetVectorSize() - index);
		ret->insertSetOfRows(new_mat, false);
		ret->setHowellized();
		return ret;
	}

	// Concatenate matrices side-by-side. Produce block matrix [this that].
	template <typename T>
	void HowellMatrix<T>::concatRight(const GenericHowellMatrix & that) {
		const HowellMatrix<T>* that_cast = dynamic_cast<const HowellMatrix<T>* >(&that);
		concatRight(*that_cast);
	}

	template <typename T>
	inline void HowellMatrix<T>::concatRight(const HowellMatrix<T>* that) {
		concatRight(*that);
	}

	template <typename T>
	void HowellMatrix<T>::concatRight(const HowellMatrix<T> & that) {
		HowellMatrix<T> other(that);

		bool old_howellized = howellized_;
		SparseTArray<T>* zero_row = new SparseTArray<T>(mat_->NumCols());

		// Pad matrices with zeros to make them the same size, as needed.
		for (size_t i = size(); i < other.size(); ++i) {
			insertRow(zero_row, false);
		}

		SparseTArray<T>* that_zero_row = new SparseTArray<T>(that.mat_->NumCols());
		for (size_t i = other.size(); i < size(); ++i) {
			other.insertRow(that_zero_row, false);
		}

		// Concatenate rows.
		mat_->insertVals(mat_->NumCols(), other.mat_->NumCols(), other.mat_);

		// Concat right doesn't change the status of Howellization
		howellized_ = old_howellized;
	}

	// Concatenate matrices vertically. Produce block matrix:
	// [ this ]
	// [ that ]
	template <typename T>
	void HowellMatrix<T>::concatDown(const GenericHowellMatrix & that) {
		const HowellMatrix<T>* that_cast = dynamic_cast<const HowellMatrix<T>* >(&that);
		concatDown(*that_cast);
	}

	template <typename T>
	inline void HowellMatrix<T>::concatDown(const HowellMatrix<T>* that) {
		concatDown(*that);
	}
	template <typename T>
	inline void HowellMatrix<T>::concatDown(const HowellMatrix<T> & that) {
		const std::vector<SparseTArray<T>* > rows = that.mat_->mat_;
		insertSetOfRows(rows);
	}

	template <typename T>
	void HowellMatrix<T>::negate() {
		for (size_t i = 0; i < size(); i++) {
			mat_->vectorMul(i, -1);
		}
		howellized_ = false;
	}

	// Insert set of rows
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<ModularSquareMatrix<T>* > & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i]->getDim() * vectorSet[i]->getDim());
			SparseTArray<T>* r = new SparseTArray<T>(vectorSet[i]->getConstMatrix(), mat_->NumCols());
			mat_->addRow(r);
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<ModularRow<T>* > & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i]->getSize());
			mat_->addRow(*(vectorSet[i]->getConstTArray()));
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<ModularRow<T> > & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i].getSize());
			mat_->addRow(*(vectorSet[i].getConstTArray()));
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<TArray<T> > & vectorSet, bool perform_howellize) {
		for (unsigned int i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i].getLength());
			mat_->addRow(vectorSet[i]);
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<SparseTArray<T> > & vectorSet, bool perform_howellize) {
		for (unsigned int i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i].GetLength());
			mat_->addRow(vectorSet[i]);
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<TArray<T>* > & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i]->getLength());
			mat_->addRow(*(vectorSet[i]));
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(std::vector<SparseTArray<T>* > & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i]->GetLength());
			mat_->addRow(vectorSet[i]);
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<SparseTArray<T>* > & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			UWAssert::UWAssert::shouldNeverHappen(mat_->NumCols() != vectorSet[i]->GetLength());
			mat_->addRow(*(vectorSet[i]));
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<std::vector<T> > & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			mat_->addRow(vectorSet[i]);
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}
	template <typename T>
	void HowellMatrix<T>::insertSetOfRows(const std::vector<const T *> & vectorSet, bool perform_howellize) {
		for (size_t i = 0; i < vectorSet.size(); i++) {
			mat_->addRow(vectorSet[i]);
		}
		howellized_ = false;
		if (perform_howellize)
			reduceToHowellForm();
	}

	/*
	* reduceToHowellForm
	*
	* Convert a set of ModularRow matrices into a set in Howell form.
	* Each matrix in the set is considered to be a row of a larger matrix
	* (referred to as M in the definition below).
	*
	* A matrix M is in Howell form iff
	* 1. M is in row-echelon form,
	* 2. the leading value of every row is a power of two,
	* 3. each leading value is the largest value in its column, and
	* 4. for every row r of M, for any p in Z, if i is the leading index of (2^p)r,
	*    then (2^p)r in rowspace([M]_i), where [M]_i denotes the matrix that
	*    consists of all rows of M whose leading index is i or greater.
	*
	* See Definition 1 in
	* Elder, M., Lim, J., Sharma, T., Andersen, T., and Reps, T.,
	* Abstract domains of affine relations. To appear in Proc. Static
	* Analysis Symposium (SAS), 2011.
	*/
	template <typename T>
	void HowellMatrix<T>::reduceToHowellForm() {
		if (howellized_)
			return;

#ifdef GATHER_STATISTICS
		std::ofstream sparsity_out("howell_matrix_sparsity.stat", std::ios_base::app);
		sparsity_out << "Row";
		for (size_t i = 0; i < num_rows(); i++)
		{
			size_t num_non_zero = 0;
			for (size_t j = 0; j < GetVectorSize(); j++)
			{
				if (getEntry(i, j) != 0)
				{
					num_non_zero++;
				}
			}
			sparsity_out << " " << num_non_zero;
		}
		sparsity_out << "\nColumn";
		for (size_t i = 0; i < GetVectorSize(); i++)
		{
			size_t num_non_zero = 0;
			for (size_t j = 0; j < num_rows(); j++)
			{
				if (getEntry(j, i) != 0)
				{
					num_non_zero++;
				}
			}
			sparsity_out << " " << num_non_zero;
		}
		sparsity_out << "\n";
		sparsity_out.close();
#endif

		size_t i, j;

		// li is an (inverse) map that maps from leading-entry position to the row
		// that has that position as its leading entry. The map may be sparse because for
		// a given leading entry position there may not be such a row. However, for a given leading entry
		// position, there still might be more than one rows. So, the map is a vector of vectors.
		// i.e., if li[i] = j, then pBasis[j].leading_entry = i (unless j == -1)
		std::vector<std::vector<size_t> > li(GetVectorSize());

		//Create the reverse map
		for (i = 0; i < mat_->NumRows(); i++)
		{
			if (mat_->GetSparseRow(i)->GetLeadingIndex() != GetVectorSize())
				li[mat_->GetSparseRow(i)->GetLeadingIndex()].push_back(i);
		}

		//Note: We'll replace mat with howellFormKS vector at the end of the function
		SparseRowMatrix<T>* howell_form_ks = new SparseRowMatrix<T>(GetVectorSize());

		//Loop invariants:
		//1)Every row in howellFormKS is in howell form(note it's not reduced at the moment).
		//2)li[k], where 0<=k<=i, either has only one entry which is already stored in howellFormKS or has zero entries.

		size_t rowToAddToHowellForm, minRank, minRowNumber, r, rprime;
		T d, dprime;
		for (i = 0; i < GetVectorSize(); i++)
		{
			// Identify the pivot element, the pivot row, and decompose pivot as d*2^r
			if (li[i].size() == 0)
				//There is no constraint with i as the leading entry; move on to the next leading entry.
				continue;
			else if (li[i].size() == 1)
			{
				//No echelon operations need to be done as this is the only row with leading entry i
				rowToAddToHowellForm = li[i][0];
				// std::cout << "i: " << i << " hey " << mat_->GetLeadingRankWithArray(rowToAddToHowellForm, mat_->GetSparseRow(rowToAddToHowellForm)) << std::endl;
				// r = mat_->GetSparseRow(rowToAddToHowellForm)->GetLeadingRank();
				// d = BitVector::div_pow2(mat_->GetSparseRow(rowToAddToHowellForm)->Get(i), r);
				r = mat_->GetLeadingRank(rowToAddToHowellForm);//r = minRank
				d = BitVector::div_pow2(mat_->Get(rowToAddToHowellForm, i), r);
			}
			else
			{
				//Find the index of the row with the minimum rank in minRowNumber. This row will be added to howell_form_ks in this iteration i.
				//For all the other rows present in li[i] perform a row reduction wrt row minRowNumber.
				minRank = mat_->GetLeadingRank(li[i][0]);
				minRowNumber = li[i][0];
				for (j = 1; j < li[i].size(); j++)
				{
					if (mat_->GetLeadingRank(li[i][j]) < minRank)
					{
						minRank = mat_->GetLeadingRank(li[i][j]);
						minRowNumber = li[i][j];
					}
				}
				rowToAddToHowellForm = minRowNumber;
				r = mat_->GetLeadingRank(minRowNumber);
				d = BitVector::div_pow2(mat_->Get(minRowNumber, i), r);

				size_t numRowsWithLeadingIndexi = li[i].size();
				//Perform row reduction operation on the other rows through row minRowNumber
				for (j = 0; j < numRowsWithLeadingIndexi; j++)
				{
					//Make sure that you don't zero out the minRank row by subtracting it by itself.
					if (li[i][j] == minRowNumber)
						continue;
					//i is the leading index for rows in li[i].
					//Decompose leading entry of li[i][j] as dprime*2^rprime
					rprime = mat_->GetLeadingRank(li[i][j]);
					dprime = BitVector::div_pow2(mat_->Get(li[i][j], i), rprime);

					//Subtract an appropriate multiple of row minRowNumber from li[i][j] to zero out the leading entry
					T t = BitVector::mult_pow2(dprime, (rprime - r));
					size_t rowToBeSubtracted = li[i][j];


					//Do operation mat[rowToBeSubtracted] = d*mat[rowToBeSubtracted] - t*mat[minRowNumber]
					mat_->vectorLinearCombination(rowToBeSubtracted, rowToBeSubtracted, (T)d, minRowNumber, (T)-t);

					//Update li to reflect the change in the leading index of rowToBeSubtracted
					if (mat_->GetLeadingIndex(rowToBeSubtracted) != GetVectorSize())
						li[mat_->GetLeadingIndex(rowToBeSubtracted)].push_back(rowToBeSubtracted);
				}

				//Update li[i] to reflect the removal of all the previous rows with leading index i.
				int index;
				for (index = numRowsWithLeadingIndexi - 1; index >= 0; index--)
				{
					if (li[i][index] != minRowNumber)
					{
						li[i].erase(li[i].begin() + index);
					}
				}
			}

			// std::cout << "d: " << d << " i: " << i << " li_size: " << li[i].size()  << std::endl;

			// Get the logical consequent of row xmat and add it to the li map and mat
			// The logical consequent of a row is obtained by making its current leading entry zero by 
			// multiplying it by an appropriate power of 2.
			// Logical consequent is all zero when the leading rank is 0. So, don't consider that case
			if (r != 0)
			{	
				SparseTArray<T>* log_conseq = mat_->GetSparseRow(rowToAddToHowellForm);
				SparseTArray<T>::vectorMul(log_conseq, BitVector::mult_pow2<T>(1, BitVector::highest_power<T>() - r));
				mat_->addRow(log_conseq);

				if (log_conseq->GetLeadingIndex() != GetVectorSize())
					li[log_conseq->GetLeadingIndex()].push_back(mat_->NumRows() - 1);
			}

			SparseTArray<T>* rowToAdd = mat_->GetSparseRow(rowToAddToHowellForm);

			//Now make the leading entry of rowToAdd a power of 2.
			SparseTArray<T>::vectorMul(rowToAdd, BitVector::multiplicative_inverse<T>(d));

			//TODO: Avoid copying the rowToAdd in the constructor
			//Insert rowToAdd in howellFormKS
			howell_form_ks->addRow(rowToAdd);

			//Reduce the rows in howell_form_ks
			T s, d;
			const T kLeadingEntry = rowToAdd->Get(i);
			for (j = 0; j < howell_form_ks->NumRows() - 1; j++)
			{
				const T kEntryji = howell_form_ks->Get(j, i);
				// All the entries above the leading index must be less than the leading index's 
				// value(which is a power of two)..
				if (kEntryji >= kLeadingEntry)
				{
					//Change kEntryji to kEntryji mod 2^r and adjust other elements in the matrix above it..
					r = rowToAdd->GetLeadingRank();
					s = kEntryji >> r; //(s*2^r) + (kEntryji mod(2^r)) = kEntryji
					// This will make the ith entry of howellFormKS[j] less than 2^r.
					// Note: This wouldn't affect the leading entry and the matrix 
					// would still be in echelon form.
					howell_form_ks->vectorLinearCombination(j, j, T(1), howell_form_ks->NumRows() - 1, T(-s));
				}
			}
		}

		mat_ = howell_form_ks;
		// std::cout << "last: " << howell_form_ks->GetLeadingRank(0) << std::endl;

		howellized_ = true;
	}

	template <typename T>
	GenericHowellMatrix* HowellMatrix<T>::generic_project(const std::vector<size_t> & projectionVariables) const {
		return project(projectionVariables);
	}

	/* project(std::vector<int> projectionVariables
	* Projects the current ks-constraints-subspace into the variables.
	* The subspace returned has k dimensions where k is the number of variables projected upon.
	* Paramter: projectionVariables is the array of dimensions on which the projection has to be done.
	* Prerequistes: All the numbers in projectionVariables must be different and must be a valid dimension(>=0 and < TwoNPlus1)
	*/
	template <typename T>
	HowellMatrix<T>* HowellMatrix<T>::project(const std::vector<size_t> & projectionVariables) const
	{
		//Reorder columns to put projectionVariables at the end
		HowellMatrix<T> reorderedHowellMatrix(*this);
		reorderedHowellMatrix.mat_->moveEnd(projectionVariables);

		//Howellize the matrix now that projectionVariables are moved to the end
		reorderedHowellMatrix.howellized_ = false;
		reorderedHowellMatrix.reduceToHowellForm();

		/* Add the relevant rows to the mat of ret*/ // TODO: exploit howell form.
		/* If a row has the leading index as a variables being projected out(ie not in projectionVariables list), then
		* Case 1: If it's rank is 1, we simply don't add the row to projected_ks_matrix
		* Case 2: Else if the rank is k, we add it after removing the first vectorSize - projectionVariables
		*/
		std::vector<SparseTArray<T> > new_mat;
		for (size_t i = 0; i < reorderedHowellMatrix.size(); i++)
		{
			if (reorderedHowellMatrix.mat_->GetLeadingIndex(i) >= GetVectorSize() - projectionVariables.size())
			{
				SparseTArray<T>* projected_row = new SparseTArray<T>(*(reorderedHowellMatrix.mat_->GetSparseRow(i)));
				projected_row->projectRight(GetVectorSize() - projectionVariables.size());
				new_mat.push_back(projected_row);
			}
		}

		//Projected rows are already in howellForm
		HowellMatrix<T>* ret = 
			new HowellMatrix<T>(projectionVariables.size());
		ret->insertSetOfRows(new_mat, false);
		ret->setHowellized();
		return ret;
	}

	template <typename T>
	HowellMatrix<T>* HowellMatrix<T>::havoc(size_t i) const
	{
		UWAssert::UWAssert::shouldNeverHappen(!(i < GetVectorSize()));

		//Project out variable i
		std::vector<size_t> projectionVar;
		for (size_t k = 0; k < GetVectorSize(); k++)
		{
			if (k != i)
				projectionVar.push_back(k);
		}
		HowellMatrix<T>* iProjectedOutMatrix = project(projectionVar);

		HowellMatrix<T>* ret = iProjectedOutMatrix;

		//Add a zero column to iProjectedOutMatrix at position i.
		ret->mat_->addCols(i, 1);

		return ret;
	}

	template <typename T>
	GenericHowellMatrix* HowellMatrix<T>::generic_havoc(std::vector<size_t> havocSet) const {
		return havoc(havocSet);
	}

	template <typename T>
	HowellMatrix<T>* HowellMatrix<T>::havoc(std::vector<size_t> havocSet) const
	{
		if (havocSet.size() == 0u)
			return new HowellMatrix(*this);

		//Project out variables in havoc set
		std::vector<size_t> projectionVar;
#ifdef UW_ASSERTS_ENABLED
		for (size_t k = 0; k < havocSet.size(); k++)
		{
			UWAssert::UWAssert::shouldNeverHappen(!(havocSet[k] < GetVectorSize()));
		}
#endif

		//from, to determine the permutation to be done after project
		size_t num_havoc_indices = 0;
		for (size_t k = 0; k < GetVectorSize(); k++)
		{
			if (std::find(havocSet.begin(), havocSet.end(), k) == havocSet.end())
			{
				projectionVar.push_back(k);
			}
			else
			{
				num_havoc_indices++;
			}
		}

		const size_t kNumHavocVars = num_havoc_indices;
		const size_t kNumNonHavocVars = GetVectorSize() - kNumHavocVars;

#ifdef UW_ASSERTS_ENABLED
		UWAssert::UWAssert::shouldNeverHappen(num_havoc_indices != kNumHavocVars);
#endif
		HowellMatrix<T>* ret = project(projectionVar);

#ifdef UW_ASSERTS_ENABLED
		UWAssert::UWAssert::shouldNeverHappen(ret->GetVectorSize() != kNumNonHavocVars);
#endif

		// std::unique_ptr<HowellMatrix<T> > ret = proj_mat;

		//Add columns for elements that are projected away.
		ret->insertCols(ret->GetVectorSize(), kNumHavocVars);
#ifdef UW_ASSERTS_ENABLED
		UWAssert::UWAssert::shouldNeverHappen(ret->GetVectorSize() != GetVectorSize());
#endif
		//Ret looks like [(coefficent unhavoc vars)  0 0...kNumHavocVars]

		//Permute cols back to the correct position
		std::vector<size_t> from, to;
		size_t crt_havoc_index = 0;
		for (size_t k = 0; k < GetVectorSize(); k++)
		{
			to.push_back(k);
			if (std::find(havocSet.begin(), havocSet.end(), k) == havocSet.end())
			{
				from.push_back(k - crt_havoc_index);
			}
			else
			{
				from.push_back(kNumNonHavocVars + crt_havoc_index);
				crt_havoc_index++;
			}
		}

		ret->permuteCols(from, to);

		return ret;
	}

	template <typename T>
	size_t HowellMatrix<T>::size() const {
		return mat_->NumRows();
	}

	template <typename T>
	size_t HowellMatrix<T>::num_rows() const {
		return mat_->NumRows();
	}

	template <typename T>
	size_t HowellMatrix<T>::num_cols() const {
		return mat_->NumCols();
	}

	template <typename T>
	size_t HowellMatrix<T>::GetVectorSize() const {
		return mat_->NumCols();
	}

	template <typename T>
	size_t HowellMatrix<T>::GetLeadingIndex(size_t r) const {
		return mat_->GetLeadingIndex(r);
	}

	template <typename T>
	size_t HowellMatrix<T>::GetLeadingRank(size_t r) const {
		return mat_->GetLeadingRank(r);
	}

	template <typename T>
	const ModularRow<T>* HowellMatrix<T>::get(size_t i) const {
#ifdef UW_ASSERTS_ENABLED
		UWAssert::UWAssert::shouldNeverHappen(!(i < mat_->NumRows()));
#endif
		const ModularRow<T>* ret =
			new ModularRow<T>(mat_->GetRow(i), mat_->NumCols());
		return ret;
	}

	template <typename T>
	SparseTArray<T>* HowellMatrix<T>::getSparseRow(size_t i) const {
#ifdef UW_ASSERTS_ENABLED
		UWAssert::UWAssert::shouldNeverHappen(!(i < mat_->NumRows()));
#endif
		return mat_->GetSparseRow(i);
	}

	template <typename T>
	bool LessThanAbsSum(const std::pair<T, T> & a1, const std::pair<T, T> & a2)
	{
		bool cmp;
		if (a1.first == a2.first)
		{
			cmp = a1.second < a2.second;
		}
		else
		{
			cmp = a1.first < a2.first;
		}
		return cmp;
	}

	template <typename T>
	std::vector<std::pair<size_t, BitVector::Value> > HowellMatrix<T>::getRowAsValueCoeffs(size_t r) const {
		SparseTArray<T>* row = getSparseRow(r);
		return getRowAsValueCoeffs(row);
	}

	template <typename T>
	std::vector<std::pair<size_t, BitVector::Value> > HowellMatrix<T>::getRowAsValueCoeffs(SparseTArray<T>* & sp_arr) const {
		std::vector<std::pair<size_t, BitVector::Value> > ret;
		const std::vector<std::pair<size_t, T> > arr = sp_arr->GetSparseVector();
		for (typename std::vector<std::pair<size_t, T> >::const_iterator it = arr.begin(); it != arr.end(); it++) {
			ret.push_back(std::make_pair(it->first, BitVector::Value(it->second)));
		}
		return ret;
	}

	// Check if r1 is a logical consequence of r2
	template <typename T>
	bool HowellMatrix<T>::IsLogicalConsequence(size_t r1, size_t r2) const {
		SparseTArray<T>* r1_sp = getSparseRow(r1);
		SparseTArray<T>* r2_sp = getSparseRow(r2);

		HowellMatrix<T> mat1(GetVectorSize());
		mat1.insertRow(r1_sp, false/*howellize*/);
		mat1.insertRow(r2_sp, true);

		HowellMatrix<T> mat2(GetVectorSize());
		mat2.insertRow(r2_sp, true);

		return (mat1 == mat2);
	}

	template <typename T>
	const std::vector<ModularRow<T>* > HowellMatrix<T>::getRows(size_t i, size_t j) const {
		std::vector<ModularRow<T>* > ret;
		for (size_t k = i; k < j; k++)
		{
			ret.push_back(new ModularRow<T>(mat_->GetRow(k), mat_->NumCols()));
		}
		return ret;
	}


	//TODO: Optimize!
	template <typename T>
	ModularSquareMatrix<T> HowellMatrix<T>::getSquareMatrix() const
	{
		//UWAssert::UWAssert::shouldNeverHappen(!howellized_);

		//Allocate space for the square matrix
		TArray<T>* sq_mat = new TArray<T>(GetVectorSize() * GetVectorSize());
		size_t pos = 0;
		for (size_t i = 0; i < mat_->NumRows(); i++)
		{
			for (size_t j = 0; j < GetVectorSize(); j++)
			{
				sq_mat->set(pos, getEntry(i, j));
				pos++;
			}
		}

		//Pad zero rows
		for (size_t i = mat_->NumRows(); i < GetVectorSize(); i++)
		{
			for (size_t j = 0; j < GetVectorSize(); j++)
			{
				sq_mat->set(pos, 0);
				pos++;
			}
		}
		ModularSquareMatrix<T> ret(sq_mat, GetVectorSize());
		return ret;
	}

	/*
	*  setEntry
	*  replaces the value at 'rowNum' row and 'columnNum' column of HowellMatrix to 'value'.
	*  This procedure does not maintain Howell form.  It is the responsibility of the client caller
	*  to eventually reestablish Howell form.
	*/
	template <typename T>
	void HowellMatrix<T>::setEntry(size_t rowNum, size_t columnNum, T value)
	{
		mat_->set(rowNum, columnNum, value);
		howellized_ = false;
	}

	template <typename T>
	void HowellMatrix<T>::setEntryAsValue(size_t , size_t , BitVector::Value )
	{
		UWAssert::UWAssert::shouldNeverHappen();
	}

	/*template <>
	inline void HowellMatrix<BV64>::setEntryAsValue(size_t rowNum, size_t columnNum, Value value)
	{
		setEntry(rowNum, columnNum, value.bv64_);
	}

	template <>
	inline void HowellMatrix<BV32>::setEntryAsValue(size_t rowNum, size_t columnNum, Value value)
	{
		setEntry(rowNum, columnNum, value.bv32_);
	}

	template <>
	inline void HowellMatrix<BV16>::setEntryAsValue(size_t rowNum, size_t columnNum, Value value)
	{
		setEntry(rowNum, columnNum, value.bv16_);
	}

	template <>
	inline void HowellMatrix<BV8>::setEntryAsValue(size_t rowNum, size_t columnNum, Value value)
	{
		setEntry(rowNum, columnNum, value.bv8_);
	}
	*/
	template <>
	inline void HowellMatrix<BitVector::BV1>::setEntryAsValue(size_t rowNum, size_t columnNum, BitVector::Value value)
	{
		setEntry(rowNum, columnNum, value.bv1_);
	}

	/*
	*  getEntry
	*  gets the value at 'rowNum' row and 'columnNum' column of HowellMatrix.
	*/
	template <typename T>
	T HowellMatrix<T>::getEntry(size_t rowNum, size_t columnNum) const
	{
		return mat_->Get(rowNum, columnNum);
	}

	template <typename T>
	BitVector::Value HowellMatrix<T>::getEntryAsValue(size_t rowNum, size_t columnNum) const
	{
		return BitVector::Value(getEntry(rowNum, columnNum));
	}

	///////////////////////////////////////////
	//Print methods
	///////////////////////////////////////////
	template <typename T>
	std::ostream & HowellMatrix<T>::print(std::ostream & out) const {
		return mat_->Print(out);
	}

	template <typename T>
	std::string HowellMatrix<T>::str() const {
		std::stringstream buf;
		print(buf);
		return buf.str();
	}


	template <typename T>
	HowellMatrix<T> HowellMatrix<T>::parse(std::istream & buf) {
		//TODO: Avoid extra copying here.
		SparseRowMatrix<T>* mat = new SparseRowMatrix<T>(SparseRowMatrix<T>::parse(buf));
		HowellMatrix<T> output(mat->NumCols());
		output.mat_ = mat;

		//TODO: find out if howellized needs to be called or alternatively change the print function to print isHowellized information
		output.howellized_ = false;
		output.reduceToHowellForm();

		return output;
	}

	template <typename T>
	HowellMatrix<T> HowellMatrix<T>::of_str(const std::string & in) {
		std::istringstream in_stream(in);
		return parse(in_stream);
	}

	template <typename T>
	std::istream & operator>> (std::istream & in, HowellMatrix<T> & h_mat) {
		h_mat = HowellMatrix<T>::parse(in);
		return in;
	}

}


#endif