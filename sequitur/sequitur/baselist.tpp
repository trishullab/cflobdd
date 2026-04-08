template<typename Child>
template<typename Function> Child * BaseList<Child>::forUntil(const Function & f)
    {
    BaseList * item = this;
    BaseList * next;
    while(item)
        {
        next = item->next_ptr;
        if(!f(static_cast<Child*>(item))) break;
        item = next;
        }
    return static_cast<Child*>(item);
    }

template<typename Child>
template<typename Function> Child * BaseList<Child>::reverseForUntil(const Function & f)
    {
    BaseList * item = this;
    BaseList * prev;
    while(item)
        {
        prev = item->prev_ptr;
        if(!f(static_cast<Child*>(item))) break;
        item = prev;
        }
    return static_cast<Child*>(item);
    }

template<typename Child>
template<typename Function> const Child * BaseList<Child>::forUntil(const Function & f) const
    {
    const BaseList * item = this;
    const BaseList * next;
    while(item)
        {
        next = item->next_ptr;
        if(!f(static_cast<const Child*>(item))) break;
        item = next;
        }
    return static_cast<const Child*>(item);
    }

template<typename Child>
template<typename Function> const Child * BaseList<Child>::reverseForUntil(const Function & f) const
    {
    const BaseList * item = this;
    const BaseList * prev;
    while(item)
        {
        prev = item->prev_ptr;
        if(f(static_cast<const Child*>(item))) break;
        item = prev;
        }
    return static_cast<const Child*>(item);
    }

template<typename Child>
Child * BaseList<Child>::unlink(unsigned int number)
    {
    if(!number) throw std::range_error("Cannot unlick 0 items.");

    //find the required items:
    BaseList * this_prev = this->prev_ptr; //possibly nullptr
    BaseList * this_last = this->next(number-1);
    BaseList * this_after = this_last->next_ptr; //possibly nullptr

    //do the unlinking:
    this->prev_ptr = nullptr;
    this_last->next_ptr = nullptr;
    if(this_prev) this_prev->next_ptr = this_after;
    if(this_after) this_after->prev_ptr = this_prev;

    //return last item in unlinked chain:
    return static_cast<Child*>(this_last);
    }

template<typename Child>
Child * BaseList<Child>::reverseUnlink(unsigned int number)
    {
    if(!number) throw std::range_error("Cannot unlick 0 items.");

    //find the required items:
    BaseList * this_after = this->next_ptr; //possibly nullptr
    BaseList * this_first = this->prev(number-1);
    BaseList * this_before = this_first->prev_ptr; //possibly nullptr

    //do the unlinking:
    this_first->prev_ptr = nullptr;
    this->next_ptr = nullptr;
    if(this_before) this_before->next_ptr = this_after;
    if(this_after) this_after->prev_ptr = this_before;

    //return first item in unlinked chain:
    return static_cast<Child*>(this_first);
    }

template<typename Child>
Child * BaseList<Child>::splitBefore()
    {
    if(!this->prev_ptr) throw std::range_error("no items before to split from.");
    BaseList * prev = this->prev_ptr;
    prev->next_ptr = nullptr;
    this->prev_ptr = nullptr;
    return static_cast<Child*>(prev);
    }

template<typename Child>
Child * BaseList<Child>::splitAfter()
    {
    if(!this->next_ptr) throw std::range_error("no items after to split from.");
    BaseList * next = this->next_ptr;
    next->prev_ptr = nullptr;
    this->next_ptr = nullptr;
    return static_cast<Child*>(next);
    }

template<typename Child>
void BaseList<Child>::joinBefore(Child * other)
    {
    if(this->prev_ptr || other->next_ptr)
        throw std::runtime_error("joining must join ends of two lists.");
    this->prev_ptr = other;
    other->next_ptr = this;
    }

template<typename Child>
void BaseList<Child>::joinAfter(Child * other)
    {
    if(this->next_ptr || other->prev_ptr)
        throw std::runtime_error("joining must join ends of two lists.");
    this->next_ptr = other;
    other->prev_ptr = this;
    }

template<typename Child>
Child * BaseList<Child>::insertAfter(Child * list, unsigned int number)
    {
    if(!number) throw std::range_error("insertAfter: number items to add should be > 0");

    //get a couple of pointers:
    BaseList * this_after = this->next_ptr; //potentially nullptr
    BaseList * list_prev = list->prev_ptr; //potentially nullptr

    //link first item in new list with current item:
    this->next_ptr = list;
    list->prev_ptr = this;

    //find last and one past last item in new list:
    BaseList * list_last = list->next(number-1);
    BaseList * list_after = list_last->next_ptr; //potentially nullptr

    //link last item in new list up:
    list_last->next_ptr = this_after;
    if(this_after) this_after->prev_ptr = list_last;

    //link input list together again (we've just cut a chunk out):
    if(list_prev) list_prev->next_ptr = list_after;
    if(list_after) list_after->prev_ptr = list_prev;

    return static_cast<Child*>(list_last);
    }

template<typename Child>
Child * BaseList<Child>::insertBefore(Child *list, unsigned int number)
    {
    //get a couple of pointers:
    BaseList * this_prev = this->prev_ptr; //potentialls nullptr
    BaseList * list_prev = list->prev_ptr; //potentially nullptr

    //link first item in new list with this_prev:
    if(this_prev) this_prev->next_ptr = list;
    list->prev_ptr = this_prev;

    //find last item in new list:
    BaseList * list_last = list->next(number-1);
    BaseList * list_after = list_last->next_ptr; //potentially nullptr

    //link last item in new list with this:
    list_last->next_ptr = this;
    this->prev_ptr = list_last;

    //link input list together again:
    if(list_prev) list_prev->next_ptr = list_after;
    if(list_after) list_after->prev_ptr = list_prev;

    return static_cast<Child*>(list);
    }

template<typename Child>
const Child * BaseList<Child>::next() const
    {
    return static_cast<const Child*>(this->next_ptr);
    }

template<typename Child>
const Child * BaseList<Child>::prev() const
    {
    return static_cast<const Child*>(this->prev_ptr);
    }

template<typename Child>
Child * BaseList<Child>::next()
    {
    return static_cast<Child*>(this->next_ptr);
    }

template<typename Child>
Child * BaseList<Child>::prev()
    {
    return static_cast<Child*>(this->prev_ptr);
    }

template<typename Child>
const Child * BaseList<Child>::next(unsigned int number) const
    {
    if(!number) return static_cast<const Child*>(this);

    const BaseList * output = this;
    while(number--)
        {
        if(!output) throw std::range_error("next: trying to access item which does not exist.");
        output = output->next_ptr;
        }
    return static_cast<const Child*>(output);
    }

template<typename Child>
const Child * BaseList<Child>::prev(unsigned int number) const
    {
    if(!number) return static_cast<const Child*>(this);

    const BaseList * output = this;
    while(number--)
        {
        if(!output) throw std::range_error("prev: trying to access item which does not exist.");
        output = output->prev_ptr;
        }
    return static_cast<const Child*>(output);
    }

//non const versions of next and prev:
template<typename Child>
Child * BaseList<Child>::next(unsigned int number)
    {
    return (Child*)(const_cast<const BaseList*>(this)->next(number));
    }

template<typename Child>
Child * BaseList<Child>::prev(unsigned int number)
    {
    return (Child*)(const_cast<const BaseList*>(this)->prev(number));
    }


template<typename Child>
std::pair<const Child*,unsigned> BaseList<Child>::end() const
    {
    unsigned int count = 0;
    const BaseList * output = this;
    while(output->next_ptr)
        {
        ++count;
        output = output->next_ptr;
        }
    return std::make_pair(static_cast<const Child*>(output), count);
    }

template<typename Child>
std::pair<const Child*,unsigned> BaseList<Child>::begin() const
    {
    unsigned int count = 0;
    const BaseList * output = this;
    while(output->prev_ptr)
        {
        ++count;
        output = output->prev_ptr;
        }
    return std::make_pair(static_cast<const Child*>(output), count);
    }

template<typename Child>
std::pair<Child*,unsigned> BaseList<Child>::end()
    {
    auto output = const_cast<const BaseList*>(this)->end();
    return std::make_pair((Child*)(output.first), output.second);
    }

template<typename Child>
std::pair<Child*,unsigned> BaseList<Child>::begin()
    {
    auto output = const_cast<const BaseList*>(this)->begin();
    return std::make_pair((Child*)(output.first), output.second);
    }
