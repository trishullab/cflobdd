
#ifndef WEIGHTED_VALUES_GUARD
#define WEIGHTED_VALUES_GUARD

#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace CFL_OBDD {
    typedef boost::multiprecision::cpp_dec_float_100 FLOAT_BOOST;
}


namespace CFL_OBDD {

    // *******************************************************************************
    // Identity Value
    // ******************************************************************************* 

    template <typename T, typename Op>
    T getIdentityValue();

    // *******************************************************************************
    // Annhilator Value
    // *******************************************************************************

    template <typename T, typename Op>
    T getAnnhilatorValue();

    // *******************************************************************************
    // Annhilator Value
    // *******************************************************************************
    
    template <typename T, typename Op>
    std::tuple<T,T,T> computeInverseValue(T lw, T rw);

    template <typename T, typename Op>
    T computeComposition(T lw, T rw);

}

#endif