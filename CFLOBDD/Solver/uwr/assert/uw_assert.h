#ifndef ASSERT_H
#define ASSERT_H

#include <cstdlib> // for abort
#include <string>
#include <vector>
#include <cassert>

/*#ifdef _MSC_VER
#   define DECLSPEC_NORETURN __declspec(noreturn)
#else
#   define DECLSPEC_NORETURN
#endif
#ifdef __GNUC__
#   define ATTRIBUTE_NORETURN __attribute__ ((noreturn))
#else
#   define ATTRIBUTE_NORETURN
#endif
*/
namespace CFL_OBDD
{
	namespace UWAssert
	{
		class UWAssert
		{
		public:
			static void shouldNeverHappen()
			{
				abort();
			}

			static void shouldNeverHappen(bool cond)
			{
				if (cond) {
					abort();
				}
				assert(!cond);
			}
		};
	}
}

#endif
