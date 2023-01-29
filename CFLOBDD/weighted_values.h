
#ifndef WEIGHTED_VALUES_GUARD
#define WEIGHTED_VALUES_GUARD

#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include "fourier_semiring.h"

namespace CFL_OBDD {
    typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
    typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
    // typedef double BIG_FLOAT;
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

    template <typename T, typename Op>
    long double computeProbabilityFromAmplitude(T w);

}

#endif