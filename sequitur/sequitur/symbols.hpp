#ifndef SYMBOLS_HPP
#define SYMBOLS_HPP

#include <functional>
#include <typeinfo>
#include <stdexcept>
#include <cassert>
#include <memory>

#include "id.hpp"
#include "baselist.hpp"
#include "objectpool.hpp"

namespace jw
    {

    //### Base symbol class (which inherits list functionality): ###
    class Symbol: public BaseList<Symbol>
        {
        public:

        virtual std::unique_ptr<Symbol> clone() const = 0;

        unsigned int getHash() const
            {
            return hash_value;
            }

        //polymorphic equality check:
        virtual bool isEqual(const Symbol & other) const =0;

        Symbol(): hash_value(0)
            {}
        Symbol(unsigned int hash_val): hash_value(hash_val)
            {}

        //virtual destructor so can polymorphically delete:
        virtual ~Symbol()
            {}

        private:
        unsigned int hash_value;
        };


    //forward declaration of RuleHead for use in RuleSymbol:
    class RuleHead;

    //symbol to denote a rule, with a pointer to it:
    class RuleSymbol: public Symbol, UseObjectPool(RuleSymbol)
        {
        public:
        //RuleHead is a friend, so it can call constructor:
        friend class RuleHead;

        unsigned int getID() const;
        unsigned int getCount() const;

        bool isEqual(const Symbol &other) const
            {
            if(typeid(other) != typeid(RuleSymbol)) return false;
            auto & o = static_cast<const RuleSymbol&>(other);
            return o.getID() == this->getID();
            }

        RuleHead * getRule() const
            {
            return rule_ptr;
            }

        //clone the RuleSymbol:
        std::unique_ptr<Symbol> clone() const;

        protected:
        //can't construct element directly:
        //must be made from RuleHead::makeRuleSymbol()
        RuleSymbol(RuleHead * rule);
        RuleSymbol(const std::unique_ptr<RuleHead> rule): RuleSymbol(rule.get()) {}

        //don't allow assignment (shouldnt need too):
        RuleSymbol& operator=(const RuleSymbol&) = delete;
        RuleSymbol(const RuleSymbol&) = delete;

        private:
        //points to rulehead (doesn't need managing):
        RuleHead * rule_ptr;
        };


    template<typename Type>
    class ValueSymbol: public Symbol, UseObjectPool(ValueSymbol<Type>)
        {
        public:
        ValueSymbol(Type val):
            Symbol(std::hash<Type>()(val)), value(val)
            {}

        Type & getValue()
            {
            return value;
            }

        const Type & getValue() const
            {
            return value;
            }

        std::unique_ptr<Symbol> clone() const
            {
            return std::unique_ptr<Symbol>(new ValueSymbol(value));
            }

        bool isEqual(const Symbol &other) const
            {
            if(typeid(other) != typeid(ValueSymbol)) return false;
            auto & o = static_cast<const ValueSymbol&>(other);
            return o.value == value;
            }

        private:
        Type value;
        };

    //##### RULE TAIL #####
    class RuleTail: public Symbol, UseObjectPool(RuleTail)
        {
        bool isEqual(const Symbol &other) const
            {
            return (typeid(other) == typeid(RuleTail));
            }

        std::unique_ptr<Symbol> clone() const
            {
            return std::unique_ptr<Symbol>(new RuleTail());
            }
        };

    //##### RULE HEAD #####
    class RuleHead: public Symbol, UseObjectPool(RuleHead)
        {
        public:
        friend class RuleSymbol;

        RuleHead(unsigned int id, RuleTail * tail_in):
            count(0), rule_id(id), tail(tail_in)
            {}
        std::unique_ptr<RuleSymbol> makeRuleSymbol()
            {
            return std::unique_ptr<RuleSymbol>(new RuleSymbol(this));
            }
        unsigned int getID() const
            {
            return rule_id;
            }
        unsigned int getCount() const
            {
            return count;
            }
        RuleTail * getTail() const
            {
            return tail;
            }

        unsigned int increment()
            {
            return ++count;
            }
        unsigned int decrement()
            {
            if(!count) throw std::range_error("count not allowed to drop below 0.");
            return --count;
            }

        std::unique_ptr<Symbol> clone() const
            {
            return std::unique_ptr<Symbol>(new RuleHead(rule_id, tail));
            }

        bool isEqual(const Symbol &other) const
            {
            if(typeid(other) != typeid(RuleHead)) return false;
            auto & o = static_cast<const RuleHead&>(other);
            return rule_id == o.rule_id;
            }

        private:
        unsigned int count;
        unsigned int rule_id;
        RuleTail * tail;
        };


    }//end sequitur namespace


#endif // SYMBOLS_HPP
