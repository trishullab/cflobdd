#ifndef SEQUITUR_H
#define SEQUITUR_H

#include <unordered_map>
#include <cassert>
#include <list>
#include <tuple>
#include <stack>
#include <memory>
#include "sequitur/symbols.hpp"
#include "sequitur/symbolwrapper.hpp"
#include "sequitur/hashing.hpp"
#include "sequitur/id.hpp"

namespace jw
    {
    //global aliases to help me along:
    using uint = unsigned int;

    //declare iter class before sequitur:
    template<typename Type> class SequiturIter;

    template<typename Type>
    class Sequitur
        {
        //### iterator classes ###
        template<typename ChildIter>
        class SequiturIter
            {
            public:
            using value_type = Type;
            using const_value_type = const Type;

            SequiturIter(const Sequitur<Type> * in_parent): parent(in_parent) {}

            const_value_type & operator* ();
            const_value_type * operator-> ();

            bool operator==(const ChildIter & other) const;
            bool operator!=(const ChildIter & other) const;

            ChildIter& operator++();
            ChildIter operator++(int);
            ChildIter& operator--();
            ChildIter operator--(int);

            protected:
            virtual void forward()=0;
            virtual void backward()=0;

            const Symbol * resolveForward(const Symbol * in);
            const Symbol * resolveBackward(const Symbol * in);

            std::stack<const Symbol*> pointer_stack;
            const Symbol* current_item;
            const Sequitur<Type> * parent;
            };

        struct ForwardIter: public SequiturIter<ForwardIter>
            {
            ForwardIter(const Sequitur<Type> * in_parent, Symbol* c): SequiturIter<ForwardIter>(in_parent)
                {
                this->current_item = this->resolveForward(c);
                }
            void forward();
            void backward();
            };

        struct ReverseIter: public SequiturIter<ReverseIter>
            {
            ReverseIter(const Sequitur<Type> * in_parent, Symbol* c): SequiturIter<ReverseIter>(in_parent)
                {
                this->current_item = this->resolveBackward(c);
                }
            void forward();
            void backward();
            };
        //### end iterator classes ###

        public:

        //some types for comparison:
        const std::type_info & RuleHeadType = typeid(RuleHead);
        const std::type_info & RuleTailType = typeid(RuleTail);
        const std::type_info & RuleSymbolType = typeid(RuleSymbol);
        const std::type_info & ValueType = typeid(Value);

        //let's simplify some names:
        using DigramIndex = std::unordered_map<std::pair<SymbolWrapper,SymbolWrapper>,Symbol*>;
        using RuleIndex = std::unordered_map<uint, Symbol*>;
        using Value = ValueSymbol<Type>;

        using const_iterator = ForwardIter;
        using const_reverse_iterator = ReverseIter;

        //may be useful outside the class:
        using value_type = Type;
        using const_value_type = const Type;

        //add new symbol to list.
        //linkMade if symbol not first one added.
        void push_back(Type);

        //get const iterators:
        const_iterator begin() const {return const_iterator(this, rule_index.at(0));}
        const_iterator end() const {return const_iterator(this, sequence_end);}
        const_reverse_iterator rbegin() const {return const_reverse_iterator(this, sequence_end);}
        const_reverse_iterator rend() const {return const_reverse_iterator(this, rule_index.at(0));}

        //return const references to rules for deep inspection:
        const RuleIndex & getRules() const { return rule_index; }

        //print things:
        void printList(const Symbol *, unsigned int number) const;
        void printAll() const;
        void printSequence() const;
        void printRules() const;
        void printDigramIndex() const;

        unsigned size() const { return length; }

        //constructor:
        Sequitur();
        //move constructor:
        Sequitur(Sequitur<Type> &&)=default;
        //destructor to clean up:
        ~Sequitur();

        private:

        //waht to do when a link is made between two symbols:
        void linkMade(Symbol * first);

        //return Iter pointing to digram location, OR end of sequence if none:
        Symbol * findAndAddDigram(Symbol *first);

        //make a digram pair for use in the digram index:
        std::pair<SymbolWrapper,SymbolWrapper> makeDigramPair(Symbol * first);

        //remove a digram from the digram index:
        void removeDigramFromIndex(Symbol * first);

        //return ptr to the RuleHead if digram is complete rule:
        RuleHead * getCompleteRule(Symbol * first);

        //swaps two occurrences with a new rule:
        std::pair<Symbol*,Symbol*> swapForNewRule(Symbol * first, Symbol * other);

        //swaps a digram with an existing rule:
        Symbol *swapForExistingRule(Symbol * first, RuleHead *rule_head);

        //decrement Item if it's a rule:
        bool decrementIfRule(Symbol * item);

        //increment if it's a rule:
        bool incrementIfRule(Symbol * item);

        //check links formed around rules if they are valid etc:
        void checkNewLinks(Symbol *rule1, Symbol *rule2);
        void checkNewLinks(Symbol *rule1);

        //check for rules only used once in rule_container and expand if found:
        void expandRuleIfNecessary(Symbol * potential_rule);

        //swapForRule(symbol it 1, symbol it 2)
        //Each time a digram is replaced by a symbol representing a rule
        //    For each symbol in the digram
        //        If the symbol represents a rule
        //            If the rule represented only occurs once now
        //                Replace the symbol with contents of the rule
        //                remove the rule

        ID id_generator;
        Symbol * sequence_end;
        unsigned int length = 0;
        DigramIndex digram_index;
        RuleIndex rule_index;
        };

    //CONSTRUCTOR
    template<typename Type>
    inline Sequitur<Type>::Sequitur()
        {
        RuleTail * start_tail = new RuleTail();
        RuleHead * start_head = new RuleHead(id_generator.get(), start_tail);
        start_head->insertAfter(start_tail);
        sequence_end = start_tail;
        rule_index[start_head->getID()] = start_head;
        }

    //DESTRUCTOR
    template<typename Type>
    inline Sequitur<Type>::~Sequitur()
        {
        //delete symbols:
        auto it = rule_index.begin();
        while(it != rule_index.end())
            {
            it->second->forUntil([](Symbol * item)
                {
                delete item;
                return true;
                });
            it = rule_index.erase(it);
            }

        //clear digrams:
        digram_index.clear();

        //clear the object pools, freeing up all available memory:
        //(won't free up used memory in the case of SinglePool)
        ObjectPool<Value>::clear();
        ObjectPool<RuleHead>::clear();
        ObjectPool<RuleTail>::clear();
        ObjectPool<RuleSymbol>::clear();
        }

    template<typename Type>
    void Sequitur<Type>::push_back(Type s)
        {
        //add new symbol:
        Symbol * val = sequence_end->insertBefore(new Value(s));
        if(++length > 1)
            {
            auto one_from_end = val->prev();
            linkMade(one_from_end);
            }
        }


    //linkMade(symbol it 1, symbol it 2)
    //Each time a new link is made between two symbols
    //    Check the digram index for another instance of this digram
    //    If another instance is found and does not overlap with this one
    //        If the other instance is a complete rule
    //            swapForRule(1, 2)
    //        Otherwise
    //            Create a new rule
    //            Replace both digrams with a new symbol representing this rule
    //    Else if this digram does not exist
    //        Add this digram to the index
    //    Else
    //        Do nothing (because other digram is overlapping).
    template<typename Type>
    void Sequitur<Type>::linkMade(Symbol * first)
        {
        assert(first != nullptr && "###linkMade: No nullptr expected here###");
        assert(first->isNext() && "###linkMade: digram has only one symbol###");

        //return iter to match if there's anything to do with it:
        // - match does not overlap
        // - match exists (digram added if not)
        Symbol * match_location = findAndAddDigram(first);

        //if there was a match, and is no overlap:
        if(match_location != nullptr)
            {
            //get ruleHead it's pointing too, or nullptr if not rule:
            RuleHead * rule_head = getCompleteRule(match_location);

            //if other occurrence is not a rule:
            if(!rule_head)
                {
                //match_location already in digram index, so swap that for new rule:
                auto locations = swapForNewRule(first, match_location);
                checkNewLinks(locations.first, locations.second);
                }
            //if it is a rule...
            else
                {
                //replace digram with rule symbol:
                Symbol * location = swapForExistingRule(first, rule_head);
                checkNewLinks(location);
                }
            }
        }

    template<typename Type>
    Symbol * Sequitur<Type>::findAndAddDigram(Symbol * first)
        {
        assert(first->isNext() && "###Digram is invalid!###");

        //place this pair into digram_index if it doesnt exist:
        auto out_pair = digram_index.emplace(makeDigramPair(first),first);

        //get bool indicating whether insertion took place, and iter to location:
        bool inserted = out_pair.second;
        Symbol * other_first = out_pair.first->second;

        //if already inserted, return end:
        if(inserted) return nullptr;

        //check for overlap:
        if(other_first->next() == first || other_first == first->next())
            {
            return nullptr;
            }
        else return other_first;
        }

    template<typename Type>
    std::pair<SymbolWrapper,SymbolWrapper> Sequitur<Type>::makeDigramPair(Symbol *first)
        {
        //while we can, we should not ever be making digram pairs out of ruleheads or ruletails:
        assert(first->isNext());
        assert(typeid(*first) != typeid(RuleHead*));
        assert(typeid(first->next()) != typeid(RuleTail*));

        //make a copy whatever symbols are in our digram:
        std::unique_ptr<Symbol> first_copy = first->clone();
        std::unique_ptr<Symbol> second_copy = first->next()->clone();

        //and add these to a pair, moving them to pass ownerhsip into container:
        return std::make_pair(std::move(first_copy), std::move(second_copy));
        }

    template<typename Type>
    void Sequitur<Type>::removeDigramFromIndex(Symbol *first)
        {
        if(typeid(*first) == RuleHeadType) return;
        if(typeid(*(first->next())) == RuleTailType) return;

        //digram must be pointed at Item to be removed:
        auto iter = digram_index.find(makeDigramPair(first));
        if(iter != digram_index.end() && iter->second == first)
            {
            digram_index.erase(iter);
            }
        }


    template<typename Type>
    RuleHead * Sequitur<Type>::getCompleteRule(Symbol * first)
        {
        assert(first->isNext() && "should be at least one symbol following this");

        //surrounding symbols:
        auto tail = first->next(2);
        auto head = first->prev();

        if(typeid(*head) == RuleHeadType && typeid(*tail) == RuleTailType)
            return static_cast<RuleHead*>(head);
        else return nullptr;
        }

    template<typename Type>
    std::pair<Symbol *, Symbol *> Sequitur<Type>::swapForNewRule(Symbol *match1, Symbol *match2)
        {
        assert(match1->next() && "first should be part of a digram");
        assert(match2->next() && "other should be part of digram");
        assert(match1->prev() != match2 && "should be no overlap");
        assert(match2->next() != match1 && "should be no overlap");
        assert(makeDigramPair(match1) == makeDigramPair(match2) && "digrams should be equal");

        Symbol * match1_second = match1->next();

        //make a new rule representing digram:
        RuleTail * rule_tail = new RuleTail();
        RuleHead * rule_head = new RuleHead(id_generator.get(), rule_tail);

        Symbol * rule_item1 = rule_head->insertAfter(match1->clone().release());
        Symbol * rule_item2 = rule_item1->insertAfter(match1_second->clone().release());
        rule_item2->insertAfter(rule_tail);

        //point digram_index to rule now:
        digram_index[makeDigramPair(match1)] = rule_item1;

        //point rule index to rule too:
        rule_index[rule_head->getID()] = rule_head;

        //increment count of any rules in digram, as we've added a copy:
        incrementIfRule(match1);
        incrementIfRule(match1_second);

        //swap digrams for this new rule (expands rules out if needbe):
        Symbol * loc1 = swapForExistingRule(match1, rule_head);
        Symbol * loc2 = swapForExistingRule(match2, rule_head);

        //return pair of locations:
        return std::make_pair(loc1, loc2);
        }

    template<typename Type>
    Symbol * Sequitur<Type>::swapForExistingRule(Symbol *first, RuleHead *rule_head)
        {
        assert(first->isPrev() && "should ALWAYS be one symbol before.");
        assert(first->isNext() && "incomplete digram.");
        assert(first->next()->isNext() && "should always be a tail after this digram.");
        assert(rule_index.find(rule_head->getID()) != rule_index.end() && "rule should exist in index");

        //replaces digram at first, with rule at rule_start (contains RuleHead)
        //- remove digrams around first (not first itself)
        //- decrement the count of any rules in digram to be removed
        //- increment count of new rule
        //- unlink and then delete first digram
        //- place new RuleSymbol where first was.

        Symbol * second = first->next();
        Symbol * before_digram = first->prev();

        //remove digrams around match if they exist and point to same location:
        removeDigramFromIndex(second);
        removeDigramFromIndex(before_digram);

        //remove digram itself, decrementing count of any rule symbols:
        first->unlink(2);

        decrementIfRule(first);
        decrementIfRule(second);

        //now, we can delete the original digram elements entirely:
        delete first;
        delete second;

        //insert rule in it's place, incrementing its count:
        RuleSymbol * new_rule = rule_head->makeRuleSymbol().release();
        new_rule->getRule()->increment();

        //expand any rules contained within this rule now if needbe:
        Symbol * rule_item1 = rule_head->next();
        Symbol * rule_item2 = rule_item1->next();
        expandRuleIfNecessary(rule_item1);
        expandRuleIfNecessary(rule_item2);

        //return position of rule in sequence:
        return before_digram->insertAfter(new_rule);
        }

    //decrement Item if it's a rule:
    template<typename Type>
    bool Sequitur<Type>::decrementIfRule(Symbol *item)
        {
        if(typeid(*item) == RuleSymbolType)
            {
            static_cast<RuleSymbol*>(item)->getRule()->decrement();
            return true;
            }
        else return false;
        }

    //increment if it's a rule:
    template<typename Type>
    bool Sequitur<Type>::incrementIfRule(Symbol *item)
        {
        if(typeid(*item) == RuleSymbolType)
            {
            static_cast<RuleSymbol*>(item)->getRule()->increment();
            return true;
            }
        else return false;
        }

    template<typename Type>
    void Sequitur<Type>::checkNewLinks(Symbol * rule1, Symbol * rule2)
        {
        assert(typeid(rule1) != typeid(RuleTail*) && "rule1 should never point to a RuleTail");
        assert(typeid(rule2) != typeid(RuleTail*) && "rule2 should never point to a RuleTail");

        //get any UNIQUE iterators to starts of digrams given that we inserted rule1 and rule2
        //and check them:

        //1. check digram at position rule1 if there is something valid following it.
        auto rule1_next = rule1->next();
        if(typeid(*rule1_next) != RuleTailType && typeid(*rule1) != RuleHeadType)
            linkMade(rule1);

        //2. check digram at position rule2 if there is something valid following it:
        auto rule2_next = rule2->next();
        if(typeid(*rule2_next) != RuleTailType && typeid(*rule2) != RuleHeadType)
            linkMade(rule2);

        //3. check digram at position --rule2 if it exists and is not equal to rule1:
        auto rule2_prev = rule2->prev();
        if(rule2_prev != rule1 && typeid(*rule2_prev) != RuleHeadType)
            linkMade(rule2_prev);

        //4. check digram at position --rule1 if it exists and is not equal to rule2:
        auto rule1_prev = rule1->prev();
        if(rule1_prev != rule2 && typeid(*rule1_prev) != RuleHeadType)
            linkMade(rule1_prev);
        }

    template<typename Type>
    void Sequitur<Type>::checkNewLinks(Symbol *rule1)
        {
        assert(typeid(*rule1) != typeid(RuleTail) && "rule should never point to a RuleTail");

        //1. check digram at position rule if something valid following it:
        auto rule1_next = rule1->next();
        if(typeid(*rule1_next) != RuleTailType && typeid(*rule1) != RuleHeadType)
            linkMade(rule1);

        //2. check digram at position --rule if it exists and is valid:
        auto rule1_prev = rule1->prev();
        if(typeid(*rule1_prev) != RuleHeadType)
            linkMade(rule1_prev);

        }

    template<typename Type>
    void Sequitur<Type>::expandRuleIfNecessary(Symbol *potential_rule)
        {
        assert(typeid(*potential_rule) != typeid(RuleHead));
        assert(typeid(*potential_rule) != typeid(RuleTail));

        //if it's not a rule, ignore it:
        if(typeid(*potential_rule) != RuleSymbolType) return;

        RuleSymbol * rule_symbol = static_cast<RuleSymbol*>(potential_rule);

        //a couple of checks to make sure nothing is broken:
        assert(rule_symbol->getCount() && "count should never be 0");
        assert(rule_index.find(rule_symbol->getID()) != rule_index.end() && "rule no exist!");

        //if rule symbol count is not 1 (or 0) leave it be:
        if(rule_symbol->getCount() != 1) return;

        //if we've got this far, expand the rule:
        //(by making rules hold ptr's to the Item* containing the rule head,
        // and rule heads holding ptr's to the rule tail, we could improve this)

        RuleHead * rule_head_item = rule_symbol->getRule();
        Symbol * rule_first_item = rule_head_item->next();

        Symbol * rule_tail_item = rule_head_item->getTail();
        Symbol * rule_last_item = rule_tail_item->prev();

        //check that indeed the rule_tail is a RuleTail (it always should be):
        assert(typeid(*rule_head_item) == typeid(RuleHead) && "should always be a head at start.");
        assert(typeid(*rule_tail_item) == typeid(RuleTail) && "should always be a tail at end.");

        //delete rule from rule_index and free up ID for future use:
        unsigned int rule_id = rule_head_item->getID();
        rule_index.erase(rule_id);
        id_generator.free(rule_id);

        //get items surrounding the rule symbol we wanna swap out:
        Symbol * before_potential_rule = potential_rule->prev();
        Symbol * after_potential_rule = potential_rule->next();

        //delete the digrams pointing to these items (relies on rule head):
        removeDigramFromIndex(before_potential_rule);
        removeDigramFromIndex(potential_rule);

        //separate rule from its head and tail:
        rule_head_item->splitAfter();
        rule_tail_item->splitBefore();

        //now we no longer need them, we can delete them:
        delete rule_head_item;
        delete rule_tail_item;

        //unlink and delete the rule symbol:
        potential_rule->splitBefore();
        potential_rule->splitAfter();
        delete potential_rule;

        //join up the pieces:
        before_potential_rule->joinAfter(rule_first_item);
        after_potential_rule->joinBefore(rule_last_item);

        //now, check new digrams made if they don't contain rule heads or tails:
        if(typeid(*before_potential_rule) != RuleHeadType) linkMade(before_potential_rule);
        if(typeid(*(rule_last_item->next())) != RuleTailType) linkMade(rule_last_item);
        }

    template<typename Type>
    void Sequitur<Type>::printList(const Symbol * list, unsigned int number) const
        {
        list->forUntil([&number,this](const Symbol * item)
            {
            if(typeid(*item) == ValueType)
                {
                std::cout << static_cast<const Value*>(item)->getValue() << " ";
                }
            else if(typeid(*item) == RuleSymbolType)
                {
                std::cout << "[" << static_cast<const RuleSymbol*>(item)->getID() << "] ";
                }
            else if(typeid(*item) == RuleHeadType)
                {
                std::cout << "[" << static_cast<const RuleHead*>(item)->getID()
                          << "(" << static_cast<const RuleHead*>(item)->getCount() << ")]: < ";
                }
            else if(typeid(*item) == RuleTailType)
                {
                std::cout << ">";
                }

            if(!--number) return false;
            else return true;
            });

        }

    template<typename Type>
    void Sequitur<Type>::printSequence() const
        {
        printList(rule_index.at(0), 0);
        std::cout << std::endl;
        }

    template<typename Type>
    void Sequitur<Type>::printRules() const
        {
        for(const auto & rule_pair : rule_index )
            {
            std::cout << rule_pair.first << ": ";
            printList(rule_pair.second, 0);
            std::cout << std::endl;
            }
        }

    template<typename Type>
    void Sequitur<Type>::printAll() const
        {
        //print out rules:
        for(const auto & rule_pair : rule_index )
            {
            std::cout << rule_pair.first << ": ";
            printList(rule_pair.second, 0);
            std::cout << std::endl;
            }

        printDigramIndex();
        std::cout << std::endl;
        }

    template<typename Type>
    void Sequitur<Type>::printDigramIndex() const
        {
        for(auto & i : digram_index)
            {
            printList(i.second, 2);
            std::cout << ", ";
            }
        std::cout << std::endl;
        }


    //###############################
    //### Sequitur Iterator Class ###
    //###############################

    template<typename Type>
    template<typename ChildIter>
    typename Sequitur<Type>::const_value_type &
    Sequitur<Type>::SequiturIter<ChildIter>::operator*()
        {
        return static_cast<const Value*>(current_item)->getValue();
        }

    template<typename Type>
    template<typename ChildIter>
    typename Sequitur<Type>::const_value_type *
    Sequitur<Type>::SequiturIter<ChildIter>::operator->()
        {
        return &(static_cast<const Value*>(current_item)->getValue());
        }

    template<typename Type>
    template<typename ChildIter>
    bool Sequitur<Type>::SequiturIter<ChildIter>::operator==(const ChildIter & other) const
        {
        return current_item == other.current_item &&
                pointer_stack == other.pointer_stack;
        }

    template<typename Type>
    template<typename ChildIter>
    bool Sequitur<Type>::SequiturIter<ChildIter>::operator!=(const ChildIter & other) const
        {
        return !(*this == other);
        }

    template<typename Type>
    template<typename ChildIter>
    const Symbol * Sequitur<Type>::SequiturIter<ChildIter>::resolveForward(const Symbol * in)
        {
        const std::type_info & type = typeid(*in);
        const Symbol * output;

        if(type == parent->ValueType)
            {
            output = in;
            }
        else if(type == parent->RuleSymbolType)
            {
            //go down one level:
            pointer_stack.push(in);
            output = resolveForward(static_cast<const RuleSymbol*>(in)->getRule()->next());
            }
        else if(type == parent->RuleHeadType)
            {
            output = resolveForward(in->next());
            }
        else //in == parent->RuleTailType
            {
            //###only hit when going forwards###
            if(pointer_stack.empty()) return in;
            //otherwise, go up one level:
            const Symbol * back = pointer_stack.top();
            pointer_stack.pop();
            output = resolveForward(back->next());
            }

        return output;
        }

    template<typename Type>
    template<typename ChildIter>
    const Symbol * Sequitur<Type>::SequiturIter<ChildIter>::resolveBackward(const Symbol * in)
        {
        const std::type_info & type = typeid(*in);
        const Symbol * output;

        if(type == parent->ValueType)
            {
            output = in;
            }
        else if(type == parent->RuleSymbolType)
            {
            //go down one level:
            pointer_stack.push(in);
            output = resolveBackward(static_cast<const RuleSymbol*>(in)->getRule()->getTail()->prev());
            }
        else if(type == parent->RuleTailType)
            {
            output = resolveBackward(in->prev());
            }
        else //in == parent->RuleHeadType
            {
            //if at the end:
            if(pointer_stack.empty()) return in;
            //otherwise, go up one level:
            const Symbol * back = pointer_stack.top();
            pointer_stack.pop();
            output = resolveBackward(back->prev());
            }

        return output;
        }

    template<typename Type>
    template<typename ChildIter>
    ChildIter& Sequitur<Type>::SequiturIter<ChildIter>::operator++()
        {
        this->forward();
        return *static_cast<ChildIter*>(this);
        }

    template<typename Type>
    template<typename ChildIter>
    ChildIter Sequitur<Type>::SequiturIter<ChildIter>::operator++(int)
        {
        ChildIter tmp(*static_cast<ChildIter*>(this));
        this->forward();
        return tmp;
        }

    template<typename Type>
    template<typename ChildIter>
    ChildIter & Sequitur<Type>::SequiturIter<ChildIter>::operator--()
        {
        this->backward();
        return *static_cast<ChildIter*>(this);
        }

    template<typename Type>
    template<typename ChildIter>
    ChildIter Sequitur<Type>::SequiturIter<ChildIter>::operator--(int)
        {
        ChildIter tmp(*static_cast<ChildIter*>(this));
        this->backward();
        return tmp;
        }

    // ##############################################################
    // ### function definitions for forward and reverse iterator: ###
    // ##############################################################

    template<typename Type>
    void Sequitur<Type>::ForwardIter::forward()
        {
        this->current_item = this->resolveForward(this->current_item->next());
        }

    template<typename Type>
    void Sequitur<Type>::ForwardIter::backward()
        {
        this->current_item = this->resolveBackward(this->current_item->prev());
        }

    template<typename Type>
    void Sequitur<Type>::ReverseIter::forward()
        {
        this->current_item = this->resolveBackward(this->current_item->prev());
        }

    template<typename Type>
    void Sequitur<Type>::ReverseIter::backward()
        {
        this->current_item = this->resolveForward(this->current_item->next());
        }


    }//end of jw namespace;
#endif // SEQUITUR_H
