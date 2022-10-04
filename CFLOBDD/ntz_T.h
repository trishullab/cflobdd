#ifndef _NTZ_T
#define _NTZ_T

#include <cassert>
#include <cstdlib>
#include <iostream>

#define POWER_OF_TWO(x) ((x) == ((x) & (-((int)(x)))))

//
// ntz
//
// Number of trailing zeros in a number
//
// The algorithm is from Warren, "Hacker's Delight", p.78
//
template <typename BV>
BV ntz(const BV & _x)
{
	if (!_x)
		return 8 * sizeof(BV);

	BV x = _x;
	x = ~x & (x - 1);
	BV n = 0;
	while (x != 0) {
		n = n + 1;
		x = x >> 1;
	}
	return n;
}

#endif
