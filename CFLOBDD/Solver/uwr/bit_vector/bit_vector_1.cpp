#include "bit_vector_1.h"
#include "bit_vector_ops.h"
#include "../matrix/HowellMatrix.h"
// Required because of dependency on highest_power
//#include "tsl/uwr/bit_vector/bit_vector_template.hpp"

namespace CFL_OBDD
{
	namespace BitVector
	{
		// --------------------------------------------------
		// Constructors
		// --------------------------------------------------
		BitVector1::BitVector1()
			: data(0)
		{
		}

		BitVector1::BitVector1(bool b)
			: data(b ? 1 : 0)
		{
		}

		BitVector1::BitVector1(unsigned char b)
			: data(b ? 1 : 0)
		{
		}

		BitVector1::BitVector1(int i)
			: data(i == 0 ? 0 : 1)
		{
		}

		BitVector1::BitVector1(unsigned i)
			: data(i == 0 ? 0 : 1)
		{
		}

		BitVector1::BitVector1(long long unsigned i)
			: data(i == 0 ? 0 : 1)
		{}
		// --------------------------------------------------
		// Destructor
		// --------------------------------------------------
		BitVector1::~BitVector1()
		{}

		// --------------------------------------------------
		// Arithmetic operators
		// --------------------------------------------------

		// Unary operators ----------------------------------
		const BitVector1& BitVector1::operator+() const {
			return *this;
		}

		const BitVector1 BitVector1::operator-() const {
			return *this;
			//return BitVector1(1 - data);
		}

		// Binary operations -------------------------------------------
		const BitVector1 BitVector1::operator+(const BitVector1& a) const {
			return BitVector1(data ^ a.data);
		}

		const BitVector1 BitVector1::operator-(const BitVector1& a) const {
			return BitVector1(data ^ a.data);
		}

		const BitVector1 BitVector1::operator*(const BitVector1& a) const {
			return BitVector1(data * a.data);
		}

		// --------------------------------------------------
		// Bitwise operators --------------------------------
		// --------------------------------------------------

		// Unary operators -----------------------------------
		const BitVector1 BitVector1::operator~() const {
			return BitVector1(1 - data);
		}

		// Binary operations --------------------------------
		const BitVector1 BitVector1::operator^(const BitVector1& a) const {
			return BitVector1(data ^ a.data);
		}

		const BitVector1 BitVector1::operator&(const BitVector1& a) const {
			return BitVector1(data & a.data);
		}

		const BitVector1 BitVector1::operator|(const BitVector1& a) const {
			return BitVector1(data | a.data);
		}

		// Shift operators ------------------------------------------
		const BitVector1 BitVector1::operator<<(const unsigned int a) const {
			return BitVector1((a > 0) ? 0 : data);
		}

		const BitVector1 BitVector1::operator>>(const unsigned int a) const {
			return BitVector1((a > 0) ? 0 : data);
		}

		// --------------------------------------------------
		// Relational operators -----------------------------
		// --------------------------------------------------

		// Unary operations ---------------------------------
		bool BitVector1::operator!() const {
			return data != 0;
		}

		// Binary operations --------------------------------
		bool BitVector1::operator==(const BitVector1& a) const {
			return data == a.data;
		}

		bool BitVector1::operator!=(const BitVector1& a) const {
			return data != a.data;
		}

		bool BitVector1::operator<(const BitVector1& a) const {
			return data < a.data;
		}

		bool BitVector1::operator>(const BitVector1& a) const {
			return data > a.data;
		}

		bool BitVector1::operator<=(const BitVector1& a) const {
			return data <= a.data;
		}

		bool BitVector1::operator>=(const BitVector1& a) const {
			return data >= a.data;
		}

		bool BitVector1::operator&&(const BitVector1& a) const {
			return data && a.data;
		}

		bool BitVector1::operator||(const BitVector1& a) const {
			return data || a.data;
		}

		// -----------------------------------------------------
		// Side-effecting operations (non-const member functions)
		// -----------------------------------------------------

		// Prefix increment
		const BitVector1& BitVector1::operator++() {
			data = 1 - data;
			return *this;
		}

		// Postfix increment
		const BitVector1 BitVector1::operator++(int) {
			BitVector1 before(data);
			data = 1 - data;
			return before;
		}

		// Prefix decrement
		const BitVector1& BitVector1::operator--() {
			data = 1 - data;
			return *this;
		}

		// Postfix decrement
		const BitVector1 BitVector1::operator--(int) {
			BitVector1 before(data);
			data = 1 - data;
			return before;
		}

		// -------------------------------------------------
		// Assignment operators
		// -------------------------------------------------

		BitVector1& BitVector1::operator=(const BitVector1& a) {
			if (this != &a) {
				data = a.data;
			}
			return *this;
		}

		// Arithmetic assignment operators ---------------------
		// Note: No checks for "this != &a" because there is no
		// potential for clobbering members in this case

		BitVector1& BitVector1::operator+=(const BitVector1& a) {
			data ^= a.data;
			return *this;
		}

		BitVector1& BitVector1::operator-=(const BitVector1& a) {
			data ^= a.data;
			return *this;
		}

		BitVector1& BitVector1::operator*=(const BitVector1& a) {
			data *= a.data;
			return *this;
		}

		// Bitwise assignment operators ------------------------
		// Note: No checks for "this != &a" because there is no
		// potential for clobbering members in this case

		BitVector1& BitVector1::operator^=(const BitVector1& a) {
			data ^= a.data;
			return *this;
		}

		BitVector1& BitVector1::operator&=(const BitVector1& a) {
			data &= a.data;
			return *this;
		}

		BitVector1& BitVector1::operator|=(const BitVector1& a) {
			data |= a.data;
			return *this;
		}

		// Shift-assignment operators
		BitVector1& BitVector1::operator>>=(const unsigned int a) {
			if (a > 0) data = 0;
			return *this;
		}

		BitVector1& BitVector1::operator<<=(const unsigned int a) {
			if (a > 0) data = 0;
			return *this;
		}

		// ---------------------------------------------------
		// Miscellaneous
		// ---------------------------------------------------

		// Address-of operator -------------------------------
		BitVector1* BitVector1::operator&() {
			return this;
		}

		// Output: Write the contents to an ostream -----------------
		std::ostream& BitVector1::print(std::ostream& os) const {
			os << (int)data;
			return os;
		}

		std::ostream& operator<<(std::ostream& o, const BitVector1 & bv) {
			return bv.print(o);
		};

		std::istream& operator>> (std::istream& in, BitVector1& val)
		{
			bool v;
			in >> v;
			val = v;
			return in;
		}

	}
}