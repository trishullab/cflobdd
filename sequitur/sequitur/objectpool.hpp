#ifndef OBJECTPOOL_HPP
#define OBJECTPOOL_HPP

//simple templated object pool. deals with allocating and deleting objects of type Type.
// - can inherit from UseObjectPool(classname) to use object pool with a class
// - otherwise, ObjectPool<classname>.create and ObjectPool<classname>.remove
//   act as new and delete.

#include <stack>
#include <algorithm>
#include <utility>

namespace jw{

    template<typename Child>
    class ChunkPool
        {
        public:

        ChunkPool(unsigned int size): block_size(size)
            {
            addBlock();
            }
        ChunkPool(): ChunkPool(100) {}

        ChunkPool(ChunkPool && other)
            {
            this->blocks(std::move(other.blocks));
            this->free_spaces(std::move(other.free_spaces));
            other.current_position = 0;
            addBlock();
            }

        void * allocate()
            {
            Child * free_memory;
            if(!free_spaces.empty())
                {
                free_memory = free_spaces.top();
                free_spaces.pop();
                }
            else
                {
                if(current_position == block_size)
                    {
                    free_memory = addBlock();
                    }
                else
                    {
                    free_memory = blocks.top() + current_position;
                    }
                ++current_position;
                }
            return static_cast<void*>(free_memory);
            }

        void deallocate(void * in)
            {
            free_spaces.push(static_cast<Child*>(in));
            }

        template<typename... Args>
        Child * create(Args && ...args)
            {
            Child * free_memory = static_cast<Child*>(allocate());
            try { new (free_memory) Child(std::forward<Args>(args)...); }
            catch(...) { free_spaces.push(free_memory); throw; }
            return free_memory;
            }

        void remove(Child * o)
            {
            o->~Child();
            free_spaces.push(o);
            }

        void clear()
            {
            while(!blocks.empty())
                {
                ::operator delete(blocks.top());
                blocks.pop();
                }
            while(!free_spaces.empty()) free_spaces.pop();
            addBlock();
            current_position = 0;
            }

        ~ChunkPool()
            {
            while(!blocks.empty())
                {
                ::operator delete(blocks.top());
                blocks.pop();
                }
            //other stuff handled now dynamically allcoated memory dealt with.
            }

        private:

        ChunkPool(const ChunkPool & other)=delete;

        Child * addBlock()
            {
            //allocate new block of memory:
            Child * block = static_cast<Child*>(::operator new(sizeof(Child)*block_size));
            blocks.push(block);
            current_position = 0;
            return block;
            }

        std::stack<Child*> blocks;
        std::stack<Child*> free_spaces;
        unsigned current_position;
        const unsigned block_size;
        };




    template<typename Type>
    class SinglePool
        {
        public:

        SinglePool() {}
        SinglePool(unsigned) {}
        SinglePool(SinglePool && other)
            {
            this->free(std::move(other.free));
            }

        void * allocate()
            {
            void * place;
            if(!free.empty())
                {
                place = static_cast<void*>(free.top());
                free.pop();
                }
            else
                {
                place = operator new(sizeof(Type));
                }
            return place;
            }

        void deallocate(void * o)
            {
            free.push(static_cast<Type*>(o));
            }

        template<typename... Args>
        Type * create(Args && ...args)
            {
            Type * place = (Type*)(allocate());
            try{ new (place) Type(std::forward<Args>(args)...); }
            catch(...) { free.push(place); throw; }
            return place;
            }

        void remove(Type * o)
            {
            o->~Type();
            free.push(o);
            }

        void clear()
            {
            while(!free.empty())
                {
                ::operator delete(free.top());
                free.pop();
                }
            }

        ~SinglePool()
            {
            clear();
            }

        private:
        SinglePool(const SinglePool &)=delete;
        std::stack<Type*> free;
        };


    template<typename Child, template<typename=Child> class PoolClass = SinglePool >
    class ObjectPool
        {
        public:
        //provide access through an object class to clear the pool:
        static void clear()
            {
            pool.clear();
            }

        static void * operator new(size_t)
            {
            return pool.allocate();
            }
        static void operator delete(void* in)
            {
            pool.deallocate(in);
            }

        private:
        //protected ctors so only Child can make instances:
        ObjectPool() {}
        ObjectPool(const ObjectPool &) {}

        static PoolClass<Child> pool;

        //Child class can make instances:
        friend Child;
        };

    template<typename Child, template<typename> class PoolClass>
    PoolClass<Child> ObjectPool<Child,PoolClass>::pool(1000);

    //inherit from UseObjectPool(classname} to use the object pool:
    // - disallows any inheritance from this class (would cause issues)
    // - virtual keyword means it's constructed from the most derived
    //   class, rather than the next one down. Which is impossible.
    // - second parameter specifies the pool to use (hidden beyond here)
    #define UseObjectPool( x ) public virtual ObjectPool<x>

    }//end of jw namespace.
#endif // OBJECTPOOL_HPP
