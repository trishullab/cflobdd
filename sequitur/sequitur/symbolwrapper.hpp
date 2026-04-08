#ifndef SYMBOLWRAPPER_H
#define SYMBOLWRAPPER_H

#include <cassert>

#include "symbols.hpp"

//
// Symbol pointers and unique_ptrs can be implicitly cast to SymbolWrappers.
//
// SymbolWrapper enables:
// - hashing based on underlying Symbol
// - equality operator== on underlying Symbol
//
// Normally, we cannot use == operator on pointers, so we need to wrap them
// to enable this and maintain polymorphism.
//
// a little like std::reference_wrapper except works direct with what I need it for.

namespace jw
    {
    class Symbol;

    class SymbolWrapper
        {
        public:

        //takes ownership of raw pointer:
        SymbolWrapper(Symbol * a): val(a)
            {
            assert(val != nullptr && "###SymbolWrapper: should never be passed a nullptr!###");
            }

        //takes ownership of unique_ptr:
        template<typename T>
        SymbolWrapper(std::unique_ptr<T> && a): val(std::move(a))
            { }

        Symbol & get() const
            {
            return *val;
            }

        bool operator==(const SymbolWrapper & other) const
            {
            return val->isEqual(*(other.val));
            }

        private:
        //so we dont have to write move ctor or dtor:
        std::unique_ptr<Symbol> val;
        };


    }//end jw namespace


//### Now, let's override the custom hash function for SymbolWrapper:
namespace std
    {
    using namespace jw;

    //overload of hash<t> for Symbol to get precalculated Hash:
    template <>
    struct hash<SymbolWrapper>
        {
        size_t operator()(const SymbolWrapper & x) const
            {
            return x.get().getHash();
            }
        };


    };

#endif // SYMBOLWRAPPER_H
