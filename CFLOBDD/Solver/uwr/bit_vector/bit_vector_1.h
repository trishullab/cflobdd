#ifndef BIT_VECTOR_1_H
#define BIT_VECTOR_1_H

#include <iostream>

namespace CFL_OBDD
{
	namespace BitVector
	{
		class BitVector1 {
		public:
			// For compatibility with bit_vector_template's ::T typedef.
			typedef unsigned char T;

			// Constructors -------------------------------------
			BitVector1();
			BitVector1(bool b);
			BitVector1(unsigned char b);
			BitVector1(int i);
			BitVector1(unsigned i);
			BitVector1(long long unsigned i);

			// Destructor --------------------------------------
			~BitVector1();


			// --------------------------------------------------
			// Arithmetic operators
			// --------------------------------------------------

			// Unary operators ----------------------------------
			const BitVector1& operator+() const;
			const BitVector1 operator-() const;

			// Binary operations -------------------------------------------
			const BitVector1 operator+(const BitVector1& a) const;
			const BitVector1 operator-(const BitVector1& a) const;
			const BitVector1 operator*(const BitVector1& a) const;

			// --------------------------------------------------
			// Bitwise operators --------------------------------
			// --------------------------------------------------

			// Unary operators -----------------------------------
			const BitVector1 operator~() const;

			// Binary operations --------------------------------
			const BitVector1 operator^(const BitVector1& a) const;
			const BitVector1 operator&(const BitVector1& a) const;
			const BitVector1 operator|(const BitVector1& a) const;

			// Shift operators ------------------------------------------
			const BitVector1 operator<<(const unsigned int a) const;
			const BitVector1 operator>>(const unsigned int a) const;

			// --------------------------------------------------
			// Relational operators -----------------------------
			// --------------------------------------------------

			// Unary operations ---------------------------------
			bool operator!() const;

			// Binary operations --------------------------------
			bool operator==(const BitVector1& a) const;
			bool operator!=(const BitVector1& a) const;
			bool operator<(const BitVector1& a) const;
			bool operator>(const BitVector1& a) const;
			bool operator<=(const BitVector1& a) const;
			bool operator>=(const BitVector1& a) const;
			bool operator&&(const BitVector1& a) const;
			bool operator||(const BitVector1& a) const;

			// -----------------------------------------------------
			// Side-effecting operations (non-const member functions)
			// -----------------------------------------------------

			const BitVector1& operator++();   // Prefix increment
			const BitVector1 operator++(int); // Postfix increment
			const BitVector1& operator--();   // Prefix decrement
			const BitVector1 operator--(int);// Postfix decrement

			// -------------------------------------------------
			// Assignment operators
			// -------------------------------------------------

			BitVector1& operator=(const BitVector1& a);

			// Arithmetic assignment operators ---------------------
			BitVector1& operator+=(const BitVector1& a);
			BitVector1& operator-=(const BitVector1& a);
			BitVector1& operator*=(const BitVector1& a);

			// Bitwise assignment operators ------------------------
			BitVector1& operator^=(const BitVector1& a);
			BitVector1& operator&=(const BitVector1& a);
			BitVector1& operator|=(const BitVector1& a);

			// Shift-assignment operators
			BitVector1& operator>>=(const unsigned int a);
			BitVector1& operator<<=(const unsigned int a);

			// ---------------------------------------------------
			// Miscellaneous
			// ---------------------------------------------------

			// Address-of operator -------------------------------
			BitVector1* operator&();

			// TODO:
			BitVector1 operator%(const BitVector1& a) const { return BitVector1(data % a.data); }
			BitVector1 operator/(const BitVector1& a) const { return BitVector1(data / a.data); }
			bool isZero() const { return data == 0; }
			bool isOne() const { return data == 1; }
			static BitVector1 zero() { return BitVector1(0); }
			static BitVector1 one() { return BitVector1(1); }

			size_t pop_count() const { return size_t(data); }

			unsigned char get_value() const { return data; }

			// Output: Write the contents to an ostream -----------------
			std::ostream& print(std::ostream& os) const;

		private:
			unsigned char data;
		};

		std::ostream& operator<<(std::ostream& o, const BitVector1 & bv);
		std::istream & operator>> (std::istream& in, BitVector1& val);

		typedef BitVector1 BV1;
	}
}
#endif

