#ifndef BaseList_HPP
#define BaseList_HPP

// Inherit from this, passing the inherited type to the base, to endow your class with linked list properties.

#include <functional>
#include <stdexcept>
#include <utility>

namespace jw
    {

    template<typename Child>
    class BaseList
        {
        public:

        virtual ~BaseList() {};

        //run a function (returning bool) for each item until function returns false.
        // - returns pointer to false returning element or nullptr if finished.
        template<typename Function> Child * forUntil(const Function & f);
        template<typename Function> Child * reverseForUntil(const Function & f);

        template<typename Function> const Child * forUntil(const Function & f) const;
        template<typename Function> const Child * reverseForUntil(const Function & f) const;

        //unlink items, and return last item in unlinked chain:
        Child * unlink(unsigned int number=1);
        Child * reverseUnlink(unsigned int number=1);

        //sever links, splitting the list in two:
        Child * splitBefore();
        Child * splitAfter();

        //join two lists together (must be working with ends of lists):
        void joinBefore(Child * other);
        void joinAfter(Child * other);

        //insert items after/before current item, and return last item inserted:
        Child * insertAfter(Child * list, unsigned int number=1);
        Child * insertBefore(Child * list, unsigned int number=1);

        //emplace items using forward to child constructor:
        template<typename... Args>
        Child * emplaceAfter(Args && ...args)
            { return this->insertAfter(new Child(std::forward<Args>(args)...)); }

        template<typename... Args>
        Child * emplaceBefore(Args && ...args)
            { return this->insertBefore(new Child(std::forward<Args>(args)...)); }

        //return pointer to next/prev element:
        const Child * next() const;
        const Child * prev() const;
        Child * next();
        Child * prev();

        //return pointer to element some number away:
        const Child * next(unsigned int number) const;
        const Child * prev(unsigned int number) const;
        Child * next(unsigned int number);
        Child * prev(unsigned int number);

        //is there a next/prev item from this?
        bool isNext() const { return next_ptr? true : false; }
        bool isPrev() const { return prev_ptr? true : false; }

        //return a pair containing pointer to begin/end item, and count of advances made:
        std::pair<const Child *, unsigned> end() const;
        std::pair<const Child *,unsigned> begin() const;

        std::pair<Child *, unsigned> end();
        std::pair<Child *,unsigned> begin();


        private:

        BaseList * next_ptr = nullptr;
        BaseList * prev_ptr = nullptr;
        };

    //definitions in here:
    #include "baselist.tpp"
    } //end jw namespace

#endif // BaseList_HPP
