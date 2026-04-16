#ifndef PAIR_HASH_HPP
#define PAIR_HASH_HPP

#include <tuple>
#include <cstdlib>

namespace std
    {
    namespace
        {
        // Code from boost:
        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
            {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }

        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
            {
            static void apply(size_t& seed, Tuple const& tuple)
                {
                HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
                hash_combine(seed, get<Index>(tuple));
                }
            };

        template <class Tuple>
        struct HashValueImpl<Tuple,0>
            {
            static void apply(size_t& seed, Tuple const& tuple)
                {
                hash_combine(seed, get<0>(tuple));
                }
            };
        }

    //hash class overload for a generic tuple:
    template <typename ... TT>
    struct hash<std::tuple<TT...>>
        {
        size_t operator()(std::tuple<TT...> const& tt) const
            {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
            return seed;
            }
        };

    //hash class overload for a pair:
    template<typename S, typename T>
    struct hash<pair<S, T>>
        {
        inline size_t operator()(const pair<S, T> & v) const
            {
            size_t seed = 0;
            hash_combine(seed, v.first);
            hash_combine(seed, v.second);
            return seed;
            }
        };
    }

#endif // PAIR_HASH_HPP
