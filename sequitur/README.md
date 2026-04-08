On getting back into C++ a few years ago, I decided to have a go at reimplementing the Sequitur algorithm, which I had originally implemented in Javascript for use in my [compression application][compressor]. Despite the inherent advantages of Javascript as a dynamic language, it was never an ideal choice for an algorithm that requires lower level optimisation to perform at its best.

My C++ implementation therefore has a few advantages:

- It's significantly faster (somewhere around 10x).
- It's templated, and so can be used on any arbitrary type.
- Due to it's similarity to STL containers, it is also easy to use.

I have put some effort into optimising it speed-wise, although I'm sure it could be taken further. To this end, I have implemented my own doubly linked list (in the form of an invisible inheritable base class which provides the relevant functionality), and a custom memory pool for mitigating the slowdown caused by the frequent creation and deletion of `Symbols` as the algorithm reorganises things.

The template nature of this implementation means that it can be applied to any primitive type and any class that implements an equality (==) operator for equality testing. Running it with the `char` type, it can process binary data in the region of 800KB-1MB/sec, and works roughly in linear time (in other words, time taken is proportional to the number of inputs it receives). Memory consumption, while still better than that in my Javascript implementation by a fair degree, is still rather high, as between 32 and 48 bytes are required per symbol (for the `char` version, and on a 64 bit machine). That said, there is not a lot of room for improvement here if we are to preserve linear time.

Anyway, it is available for download here and can be freely used for any purpose you like.

# Basic Usage

The download includes a `main.cpp`, which creates a simple program that can be passed a filename to work on, and will process it and then print out some details. You should be able to compile using `clang` (tested on clang-602.0.49) with:

```
clang -lstdc++ -std=c++11 -Wc++11-extensions main.cpp sequitur/*.cpp
```

or `gcc` 4.8.* (untested on newer versions) using `g++` with:

```
g++ -std=c++11 main.cpp sequitur/*.cpp
```

Alternately, I have been kindly provided with a `CMakeLists.txt` file for `cmake`. To make use of this, do something like:

```
mkdir build
cd build
cmake ..
make -j4
```

Which handels generating a makefile, and builds into a seperate `build` directory utilising 4 cores.

To put to use in a custom project, just add `sequitur.hpp` and all of the files in `sequitur` to your project. Additionally, add `main.cpp` to compile a simple example. I have omitted unit tests in the download for simplicity.

```
//include the namespace:
using namespace jw;

int main()
	{
	//create a new sequitur instance:
	Sequitur<char> s;

	//add elements to it:
	s.push_back('l');
	s.push_back('d');

	//print how many elements have been added (2):
	cout << s.size() << endl;

	//iterate over elements, printing them:
	Sequitur<char>::const_iterator it = s.begin();
	while(it != s.end())
		{
		cout << *it << endl;
		++it;
		}
	
	//reverse iterators are also implemented:
	Sequitur<char>::const_reverse_iterator rit;
	while(rit != s.rend())
		{
		cout << *rit << endl;
		++rit;
		}
	
	return 0;	
	}
```

Internally, `s` will be creating rules to represent any repetition, whereby rule 0 is the original string. The main print functions are as follows, and allow for visualisation of what is happening internally:

```
//print all rules and digram index:
s.printAll();

//print rules:
s.printRules();

//print digram index:
s.printDigramIndex();

```

# Advanced Usage

For more complex tasks, you'll need to manually work with the rule index, which can be obtained as follows:

```
//get rules (of type std::unordered_map<unsigned, Symbol*>):
auto rules = s.getRules();
```

`Symbol` is an abstract type from which various symbols (`RuleSymbol`, `RuleHead`, `RuleTail`, and `ValueSymbol<type>`) are derived. `ValueSymbol` contains the value of whatever type we are working with (`char` in the above). The other symbols are used for internal rule manipulation and creation; notably, `RuleSymbol` is a symbol which points to a `RuleHead` in the sequence, and `RuleHead` is positioned at the start of a rule, which itself is a list of Symbols.

`Symbol` itself inherits from my `BaseList` implementation, which endows any inheriting class doubly linked list properties, so that they can be chained together/navigated over etc. It does this transparently by being passed the inheriting type as a template parameter, and outputting a pointer to this type in any functions that require a `BaseList*` type back.

Thus, to navigate through symbols, you'll need to inspect the BaseList implementation for details. Iterating over list elements can be done as follows, using the `next()` function to obtain a pointer to the next rule (or a nullptr if one does not exist):

```
auto rules = s.getRules();
for(auto rule : rules)
  {
  //rule is of type <unsigned,Symbol*>
  
  //print rule ID:
  cout << "rule ID: " << rule.first << endl;
  
  //get first rule symbol:
  auto current = rule.second;
  while(current)
      {
      //print the type of the symbol:
      cout << typeid(*current).name() << endl;
      
      //get the next symbol:
      current = current->next();
      }
  }
```

From `RuleSymbol`s and `RuleHead`s, you can get the number of times the rule occurs, and the rule ID. From `ValueSymbol`s, you can get the value stored at that location, as illustrated (expanding on the previous example):

```
auto rules = s.getRules();
for(auto rule : rules)
  {
  //rule is of type <unsigned,Symbol*>
  
  //print rule ID:
  cout << "rule ID: " << rule.first << endl;
  
  //get first rule symbol:
  auto current = rule.second;
  while(current)
      {
      if(typeid(*current) == s.RuleHeadType)
          {
          auto head = static_cast<RuleHead*>(current);
          cout << "Rule Head" <<endl;
          cout << "- count: " << head->getCount() << endl;
          cout << "- ID: " << head->getID() << endl;
          }
      else if(typeid(*current) == s.RuleSymbolType)
          {
          auto rule = static_cast<RuleSymbol*>(current);
          cout << "Pointer to rule: " << rule->getID() << endl;
          }
      else if(typeid(*current) == s.ValueType)
          {
          auto value = static_cast<ValueSymbol<char>*>(current);
          cout << "Value: " << value->getValue() << endl;
          }
      else
          {
          cout << "Rule Tail" << endl;
          }
      }
  }
```

As each symbol has a different interface (they are used for different purposes), we must establish which symbol we are looking at. We can then perform a cast to the correct type, and interact with it.

A compression algorithm that wished to use sequitur would need to work through the rules in this manner in order that it can form some form of binary representation of them.

Any other functionality is left for the intrepid user to explore. Let me know however if you are interested in implementing this code in your project, and I'll help you out wherever I can!


[compressor]: http://unbui.lt/projects/compressor
