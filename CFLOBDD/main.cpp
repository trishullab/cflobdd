//
//    Copyright (c) 1999 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "tests_cfl.h"
#include <iomanip>
#define _CRTDBG_MAP_ALLOC
#include <cstdlib>

static long seed_value = 27;

int main(int argc, char * argv[])
{
	// Parse flags (--verbose / -v) and collect positional arguments
	std::vector<std::string> posArgs;
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "--verbose" || arg == "-v") {
			CFL_OBDD::CFLTests::verbose = true;
		} else {
			posArgs.push_back(arg);
		}
	}

	// Supply a default argument for when invoking from Windows (e.g., for debugging)
	if (posArgs.empty()) {
		CFL_OBDD::CFLTests::runTests("And");
	}
	else if (posArgs.size() == 1) {
		CFL_OBDD::CFLTests::runTests(posArgs[0].c_str());
	}
	else if (posArgs.size() == 2) {
		CFL_OBDD::CFLTests::runTests(posArgs[0].c_str(), atoi(posArgs[1].c_str()));
	}
	else if (posArgs.size() == 3) {
		CFL_OBDD::CFLTests::runTests(posArgs[0].c_str(), atoi(posArgs[1].c_str()), atoi(posArgs[2].c_str()));
	}
	else if (posArgs.size() >= 4) {
		CFL_OBDD::CFLTests::runTests(posArgs[0].c_str(), atoi(posArgs[1].c_str()), atoi(posArgs[2].c_str()), atoi(posArgs[3].c_str()));
	}
}
