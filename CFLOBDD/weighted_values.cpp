
#include "weighted_values.h"

namespace CFL_OBDD {

    template <>
    BIG_FLOAT getIdentityValue<BIG_FLOAT, std::multiplies<BIG_FLOAT>>()
    {
        return 1.0;
    }

    template <>
    BIG_FLOAT getAnnhilatorValue<BIG_FLOAT, std::multiplies<BIG_FLOAT>>()
    {
        return 0.0;
    }

    template <>
    std::tuple<BIG_FLOAT, BIG_FLOAT, BIG_FLOAT> computeInverseValue<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(BIG_FLOAT lw, BIG_FLOAT rw)
    {
        if (lw == 0.0)
            return std::make_tuple(rw, 0.0, 1.0);
        if (rw == 0.0)
            return std::make_tuple(lw, 1.0, 0.0);
        return std::make_tuple(lw, 1.0, rw/lw);
    }

    template <>
    BIG_FLOAT computeComposition<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(BIG_FLOAT c, BIG_FLOAT w)
    {
        return c * w;
    }

    template <>
    BIG_FLOAT computeProbabilityFromAmplitude<BIG_FLOAT,std::multiplies<BIG_FLOAT>>(BIG_FLOAT w)
    {
        return w * w;
    }

    template <>
    BIG_COMPLEX_FLOAT getIdentityValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>()
    {
        return 1.0;
    }

    template <>
    BIG_COMPLEX_FLOAT getAnnhilatorValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>()
    {
        return 0.0;
    }

    template <>
    std::tuple<BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT> computeInverseValue<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(BIG_COMPLEX_FLOAT lw, BIG_COMPLEX_FLOAT rw)
    {
        if (lw == 0.0)
            return std::make_tuple(rw, 0.0, 1.0);
        if (rw == 0.0)
            return std::make_tuple(lw, 1.0, 0.0);
        return std::make_tuple(lw, 1.0, rw/lw);
    }

    template <>
    BIG_COMPLEX_FLOAT computeComposition<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(BIG_COMPLEX_FLOAT c, BIG_COMPLEX_FLOAT w)
    {
        return c * w;
    }

    template <>
    BIG_COMPLEX_FLOAT computeProbabilityFromAmplitude<BIG_COMPLEX_FLOAT,std::multiplies<BIG_COMPLEX_FLOAT>>(BIG_COMPLEX_FLOAT w)
    {
        BIG_COMPLEX_FLOAT ans(w.real() * w.real() + w.imag() * w.imag());
        return ans;
    }

}