#include "parsing.h"
#include "../assert/uw_assert.h"
#include <iostream>
#include <algorithm>

namespace Parsing
{
	using std::string;

	// expect(in,s):
	// Expect the string 's' to be the next whitespace-delimited token of 'in'.
	// Consume that token if it is 's'; cough and die if it is not 's'.

	void expect(std::istream &in, string s) {
		UWAssert::UWAssert::shouldNeverHappen(in.bad());

		string spot;
		in >> spot;

		if (spot != s) {
			std::cerr << "Expected \"" << s
				<< "\", grabbed \"" << spot << "\""
				<< std::endl;
			getline(in, spot);
			std::cerr << "Following line is:\n"
				<< spot << std::endl;

			UWAssert::UWAssert::shouldNeverHappen();
		}
	}

	std::string expect_one(std::istream &in, std::vector< std::string > choices) {
		UWAssert::UWAssert::shouldNeverHappen(in.bad());

		string spot;
		in >> spot;

		if (find(choices.begin(), choices.end(), spot) != choices.end())
			return spot;

		std::cerr << "Expected one of { ";
		std::vector<std::string>::iterator begin, end, i;
		begin = choices.begin(); end = choices.end();
		for (i = begin; i != end; i++) {
			std::cerr << *i << " ";
		}
		std::cerr << "}, grabbed " << spot << std::endl;
		getline(in, spot);
		std::cerr << "Following line is:\n"
			<< spot << std::endl;

		UWAssert::UWAssert::shouldNeverHappen();
		return "";
	}

}
