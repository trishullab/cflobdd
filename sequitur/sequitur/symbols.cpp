#include "symbols.hpp"

namespace jw
    {

    RuleSymbol::RuleSymbol(RuleHead *rule):
        Symbol(std::hash<unsigned>()(rule->getID())),
        rule_ptr(rule)
        { }

    std::unique_ptr<Symbol> RuleSymbol::clone() const
        {
        return rule_ptr->makeRuleSymbol();
        }

    //get count of RuleHead through RuleSymbol:
    unsigned int RuleSymbol::getCount() const
        {
        return rule_ptr->getCount();
        }

    //get ID of RuleHead through RuleSymbol:
    unsigned int RuleSymbol::getID() const
        {
        return rule_ptr->getID();
        }


    }//end jw namespace
