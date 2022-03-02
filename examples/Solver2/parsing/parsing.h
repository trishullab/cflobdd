// parsing.h:
// Utility functions for dealing with input streams.
// Author: Matt Elder <elder@cs.wisc.edu>
#ifndef PARSING_H
#define PARSING_H

#include <string>
#include <istream>
#include <vector>


namespace Parsing
{
	// expect(in, s):
	// Expect the string 's' to be the next
	// whitespace-delimited token of 'in'.
	// Consume that token if it is 's';
	// cough and die if it is not 's'.
	void expect(std::istream &in, std::string s);

	std::string expect_one(std::istream &in, std::vector< std::string > choices);


}
#endif
