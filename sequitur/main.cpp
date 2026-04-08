#include <fstream>
#include <iostream>
#include <string>
#include "sequitur.hpp"

//
// This file just executes some random actions using my sequitur implementation, to help give you a feel for how it works.
//

using namespace std;
//my namespace gives access to sequitur and related:
using namespace jw;


int main(int argc, char* argv[])
    {
    //
    // Look at the first string passed in at the console, and treat it as a file path, opening said file:
    //
    if(argc != 2) cerr << "Need one argument (filename)" << endl;
    string filename = argv[1];

    ifstream input(filename, std::ios::binary);

    if(!input.is_open())
        {
        cerr << "File \"" << filename << "\" not found." << endl;
        return 1;
        }

    //
    // Make a new Sequitur for taking in the char type:
    //
    Sequitur<char> s;
    char temp_char;
    unsigned count = 0;
    while(input.get(temp_char))
        {
        //
        //add chars to Sequitur, using the familiar push_back syntax:
        //
        s.push_back(temp_char);
        //
        //record the count and print an output every 100,000 chars:
        //
        count++;
        if(count % 100000 == 0) cout << count << endl;
        }
    //clear the input flags and seek back to the beginning, so we can read it in again:
    input.clear();
    input.seekg(0);

    //
    //we can't copy sequitur, but we can move it using the C++11 std::move function.
    //This will move 's' to 's2', invalidating s:
    //
    Sequitur<char> s2(std::move(s));

    //
    //create a sequitur<char>::const_iter pointed to beginning char in s2:
    //
    auto seq_iter = s2.begin();

    //read through the file again, and compare each char with that which is stored in the sequitur container s2:
    count = 0;
    while(input.get(temp_char))
        {
        if(temp_char != *seq_iter)
            {
            cerr << "not equal at " << count << ", file: " << temp_char << " seq: " << *seq_iter << endl;
            }
        //
        // we can post or pre-increment the iterator to look at the next char:
        //
        ++seq_iter;
        ++count;
        }

    //### lets see some stats: ###

    //
    // get the rule table (of type std::unordered_map<Symbol*>):
    //
    auto rule_table = s2.getRules();
    auto it = rule_table.begin();
    unsigned symbol_total = 0;
    unsigned rule_count = 0;
    while(it != rule_table.end())
        {
        ++rule_count;
        //
        // consult baselist.hpp for list traversal commands given symbols.
        // consult symbols.hpp for the different symbol types and functions available on each
        //
        // baselist is just a linked list subclass, which the Symbol class inherits from, so
        // symbols can be linked together. Thus, the functionality in baselist.hpp will allow you
        // to traverse the symbols found.
        //
        // Symbol::next() and Symbol::prev() allow iteration through them.
        // Symbol::end()/begin() returns a pair<[final_symbol], distance>
        //
        // Here, we count how many hops it takes to get to the end. This will get us one less than
        // the total number of symbols. Since two of the symbols are RuleHead and RuleTail, we 
        // minus another one to get the number of symbols used in the rule.
        //
        symbol_total += it->second->end().second - 1;
        ++it;
        }

    //
    // a few print commands are included to quikly visualise the rules created:
    //
    s2.printRules();

    cout << "total symbols inserted: " << count << endl;
    cout << "symbols used in sequitur: " << symbol_total << endl;
    cout << "rules created: " << rule_count << endl;

    return 0;
}

