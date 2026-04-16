#ifndef JAMDAWG_ID_HPP_INCLUDED
#define JAMDAWG_ID_HPP_INCLUDED

#include <stack>
#include <iostream>

//generates a unique value when used:
class ID
    {
    private:

    unsigned int val = 0;
    std::stack<unsigned int> free_ids;

    public:

    unsigned int get()
        {
        if(!free_ids.empty())
            {
            unsigned int return_val = free_ids.top();
            free_ids.pop();
            return return_val;
            }

        return val++;
        }
    void free(unsigned int id)
        {
        free_ids.push(id);
        }
    };

#endif // JAMDAWG_ID_HPP_INCLUDED
