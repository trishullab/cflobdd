#include "../parsing/parsing.h"
#include "sparse_row_matrix.h"
#ifndef SPARSE_ROW_MATRIX_IMPL_H
#define SPARSE_ROW_MATRIX_IMPL_H
namespace CFL_OBDD
{
	namespace HowellMatrix
	{
		using std::vector;
		using std::string;

		//Constructors
		//---------------------------------------------------
		template <typename T>
		SparseRowMatrix<T>::SparseRowMatrix(size_t num_cols)
			: num_cols_(num_cols)
		{ }

		//---------------------------------------------------
		// Copy constructor
		//---------------------------------------------------
		template <typename T>
		SparseRowMatrix<T>::SparseRowMatrix(const SparseRowMatrix<T>& that)
			: num_cols_(that.num_cols_)
		{
			typename RowVector::const_iterator beg, end;
			beg = that.mat_.begin();
			end = that.mat_.end();

			while (beg != end) {
				ref_ptr<SparseTArray<T> > newCopy = new SparseTArray<T>(**beg);
				mat_.push_back(newCopy);
				beg++;
			}
		}

		//---------------------------------------------------
		// Destructor 
		//---------------------------------------------------
		template <typename T>
		SparseRowMatrix<T>::~SparseRowMatrix() {
			// Nothing needed
		}


		template <typename T>
		SparseRowMatrix<T> SparseRowMatrix<T>::operator=(const Matrix<T>& that_gen)
		{
			if (this == &that_gen)
				return *this;

			num_cols_ = that_gen.NumCols();
			mat_.clear();

			const SparseRowMatrix* that = dynamic_cast<const SparseRowMatrix*>
				(&that_gen);

			typename RowVector::const_iterator beg, end;
			beg = that->mat_.begin();
			end = that->mat_.end();

			while (beg != end) {
				ref_ptr<SparseTArray<T> > newCopy = SparseTArray<T>::copy(*beg);
				mat_.push_back(newCopy);
				beg++;
			}
			return *this;
		}

		template <typename T>
		inline bool SparseRowMatrix<T>::operator==(const Matrix<T> & that_gen) const {
			const SparseRowMatrix* that = dynamic_cast<const SparseRowMatrix*>
				(&that_gen);
			if (num_cols_ != that->num_cols_)
				return false;

			if (mat_.size() != that->mat_.size())
				return false;

			typename RowVector::const_iterator it_this;
			typename RowVector::const_iterator it_that;
			for (it_this = mat_.begin(), it_that = that->mat_.begin();
				it_this != mat_.end() && it_that != that->mat_.end();
				it_this++, it_that++)
			{
				ref_ptr<SparseTArray<T> > this_row = *it_this;
				ref_ptr<SparseTArray<T> > that_row = *it_that;
				if (*this_row != *that_row)
					return false;
			}

			return true;
		}

		//Getter functions
		template <typename T>
		size_t SparseRowMatrix<T>::NumRows() const {
			return mat_.size();
		}


		template <typename T>
		size_t SparseRowMatrix<T>::NumCols() const {
			return num_cols_;
		}

		template <typename T>
		inline T SparseRowMatrix<T>::Get(size_t r, size_t c) const {
			return mat_[r]->Get(c);
		}

		template <typename T>
		inline size_t SparseRowMatrix<T>::GetLeadingIndex(size_t r) const {
			return mat_[r]->GetLeadingIndex();
		}

		template <typename T>
		inline size_t SparseRowMatrix<T>::GetLeadingRank(size_t r) const {
			return mat_[r]->GetLeadingRank();
		}

		template <typename T>
		inline const ref_ptr<TArray<T> > SparseRowMatrix<T>::GetRow(size_t r) const {
			return mat_[r]->GetTArray();
		}

		template <typename T>
		inline const ref_ptr<SparseTArray<T> > SparseRowMatrix<T>::GetConstSparseRow(size_t r) const {
			return mat_[r];
		}

		template <typename T>
		inline ref_ptr<SparseTArray<T> > SparseRowMatrix<T>::GetSparseRow(size_t r) {
			return mat_[r];
		}

		template <typename T>
		inline const ref_ptr<SparseTArray<T> > SparseRowMatrix<T>::operator[](size_t r) const {
			return mat_[r];
		}

		//Setter Function    
		template <typename T>
		inline void SparseRowMatrix<T>::set(size_t r, size_t c, T val) {
			mat_[r]->set(c, val);
		}

		//Modifier function
		template <typename T>
		void SparseRowMatrix<T>::clear() {
			mat_.clear();
		}

		//Row addition functions
		template <typename T>
		void SparseRowMatrix<T>::addRow(const T* row) {
			ref_ptr<SparseTArray<T> > r = new SparseTArray<T>(row, NumCols());
			mat_.push_back(r);
		}

		template <typename T>
		void SparseRowMatrix<T>::addRow(const std::vector<T> & row) {
			ref_ptr<SparseTArray<T> > r = new SparseTArray<T>(row);
			mat_.push_back(r);
		}

		template <typename T>
		void SparseRowMatrix<T>::addRow(const TArray<T>& row) {
			ref_ptr<SparseTArray<T> > r = new SparseTArray<T>(row);
			mat_.push_back(r);
		}

		template <typename T>
		void SparseRowMatrix<T>::addRow(const SparseTArray<T>& row) {
			ref_ptr<SparseTArray<T> > r = new SparseTArray<T>(row);
			mat_.push_back(r);
		}

		template <typename T>
		void SparseRowMatrix<T>::addRow(ref_ptr<SparseTArray<T> > & row) {
			mat_.push_back(row);
		}

		template <typename T>
		void SparseRowMatrix<T>::addRow(size_t index) {
			ref_ptr<SparseTArray<T> > r = new SparseTArray<T>(*(mat_[index]));
			mat_.push_back(r);
		}

		template <typename T>
		void SparseRowMatrix<T>::removeRows(const std::vector<size_t>& proj) {
			std::vector<size_t> proj_sorted(proj);
			std::sort(proj_sorted.begin(), proj_sorted.end());
			typename std::vector<size_t>::reverse_iterator it = proj_sorted.rbegin();
			for (; it != proj_sorted.rend(); it++)
			{
				mat_.erase(mat_.begin() + *it);
			}
		}

		//Column operations
		template <typename T>
		void SparseRowMatrix<T>::addCols(size_t begin, size_t length) {
			num_cols_ = num_cols_ + length;
			typename RowVector::iterator it = mat_.begin();
			for (; it != mat_.end(); it++)
			{
				(*it)->insertVals(begin, length);
			}
		}

		template <typename T>
		void SparseRowMatrix<T>::removeCols(const std::vector<size_t>& proj) {
			num_cols_ = num_cols_ - proj.size();
			typename RowVector::iterator it = mat_.begin();
			for (; it != mat_.end(); it++)
			{
				(*it)->projectOut(proj);
			}

		}

		template <typename T>
		void SparseRowMatrix<T>::permuteCols(const std::vector<size_t>& from,
			const std::vector<size_t>& to) {
			typename RowVector::iterator it = mat_.begin();
			for (; it != mat_.end(); it++)
			{
				(*it)->permuteCols(from, to);
			}

		}

		template <typename T>
		void SparseRowMatrix<T>::insertVals(size_t begin, size_t len,
			const ref_ptr<SparseRowMatrix<T> > & vals)
		{
			typename RowVector::iterator it = mat_.begin();
			for (; it != mat_.end(); it++)
			{
				(*it)->insertVals(begin, len, *(vals->GetConstSparseRow(it - mat_.begin())));
			}
			num_cols_ += len;
		}

		template <typename T>
		void SparseRowMatrix<T>::moveEnd(const std::vector<size_t>& m) {
			std::vector<size_t> from, to;

			size_t t_count = 0;
			for (size_t f_count = 0; f_count < NumCols(); f_count++)
			{
				from.push_back(f_count);

				typename std::vector<size_t>::const_iterator count_loc =
					std::find(m.begin(), m.end(), f_count);

				if (count_loc == m.end())
				{
					to.push_back(t_count);
					t_count++;
				}
				else
				{
					//Add index to to depending on the location of index f_count in m.
					to.push_back(NumCols() - m.size() + count_loc - m.begin());
				}
			}

			UWAssert::UWAssert::shouldNeverHappen(t_count + m.size() != NumCols());

			permuteCols(from, to);
		}

		///////////////////////////////////////////
		//Print methods
		///////////////////////////////////////////
		template <typename T>
		std::ostream & SparseRowMatrix<T>::Print(std::ostream & out) const {
			std::ios::fmtflags oldsettings = out.flags();
			out.setf(std::ios::fixed, std::ios::floatfield);
			out.setf(std::ios::dec, std::ios::basefield);

			out << "[";
			out << " cols " << NumCols()
				<< " rows " << NumRows() << "\n";

			for (size_t i = 0; i < NumRows(); i++)
			{
				mat_[i]->Print(out);
				out << "\n";
			}

			out << "]";
			out << "\n";

			out.flags(oldsettings);

			return out;
		}

		template <typename T>
		SparseRowMatrix<T> SparseRowMatrix<T>::parse(std::istream & buf) {
			std::ios::fmtflags oldsettings = buf.flags();
			buf.setf(std::ios::fixed, std::ios::floatfield);
			buf.setf(std::ios::dec, std::ios::basefield);

			std::vector<ref_ptr<SparseTArray<T> > > rows;

			size_t num_cols, num_rows;

			Parsing::expect(buf, "[");
			Parsing::expect(buf, "cols"); buf >> num_cols;
			Parsing::expect(buf, "rows"); buf >> num_rows;

			for (size_t row = 0; row < num_rows && buf.good(); row++) {
				ref_ptr<SparseTArray<T> > this_row = new SparseTArray<T>(num_cols);
				size_t col;
				for (col = 0; col < num_cols && buf.good(); col++) {
					T value;
					buf >> value;
					this_row->set(col, value);
				}
				UWAssert::UWAssert::shouldNeverHappen(col != num_cols);
				rows.push_back(this_row);
			}
			Parsing::expect(buf, "]");
			SparseRowMatrix<T> output(num_cols);

			output.mat_ = rows;

			buf.flags(oldsettings);

			return output;
		}

		//More row manipulations functions
		template <typename T>
		void SparseRowMatrix<T>::vectorMul(size_t out_index, T mult)
		{
			SparseTArray<T>::vectorMul(*mat_[out_index], mult);
		}


		template <typename T>
		void SparseRowMatrix<T>::vectorLinearCombination(size_t out_index, size_t in_index_1, T mult_1, size_t in_index_2, T mult_2)
		{
			auto t = mat_[in_index_2];
			SparseTArray<T>::vectorLinearCombination(*mat_[out_index], *mat_[in_index_1], mult_1, *mat_[in_index_2], mult_2);
		}
	}

}

#endif
