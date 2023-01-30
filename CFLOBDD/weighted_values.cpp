
#include "weighted_values.h"

namespace CFL_OBDD {

    template <>
    long double getIdentityValue<long double, std::multiplies<long double>>()
    {
        return 1.0;
    }

    template <>
    long double getAnnhilatorValue<long double, std::multiplies<long double>>()
    {
        return 0.0;
    }

    template <>
    std::tuple<long double, long double, long double> computeInverseValue<long double, std::multiplies<long double>>(long double lw, long double rw)
    {
        if (lw == 0.0)
            return std::make_tuple(rw, 0.0, (rw == 0.0) ? 0.0 : 1.0);
        if (rw == 0.0)
            return std::make_tuple(lw, 1.0, 0.0);
        return std::make_tuple(lw, 1.0, rw/lw);
    }

    template <>
    long double computeComposition<long double, std::multiplies<long double>>(long double c, long double w)
    {
        return c * w;
    }

    template <>
    long double computeProbabilityFromAmplitude<long double,std::multiplies<long double>>(long double w)
    {
        long double ans = (w * w);
        return ans;
    }

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
            return std::make_tuple(rw, 0.0, (rw == 0.0) ? 0.0 : 1.0);
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
    long double computeProbabilityFromAmplitude<BIG_FLOAT,std::multiplies<BIG_FLOAT>>(BIG_FLOAT w)
    {
        long double ans = (w * w).convert_to<long double>();
        return ans;
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
            return std::make_tuple(rw, 0.0, (rw == 0.0) ? 0.0 : 1.0);
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
    long double computeProbabilityFromAmplitude<BIG_COMPLEX_FLOAT,std::multiplies<BIG_COMPLEX_FLOAT>>(BIG_COMPLEX_FLOAT w)
    {
        auto ans = (w.real() * w.real() + w.imag() * w.imag()).convert_to<long double>();
        return ans;
    }

    template <>
    fourierSemiring getIdentityValue<fourierSemiring, std::multiplies<fourierSemiring>>()
    {
        return fourierSemiring(1, 1);
    }

    template <>
    fourierSemiring getAnnhilatorValue<fourierSemiring, std::multiplies<fourierSemiring>>()
    {
        return fourierSemiring(0, 1);
    }

    template <>
    std::tuple<fourierSemiring, fourierSemiring, fourierSemiring> computeInverseValue<fourierSemiring, std::multiplies<fourierSemiring>>(fourierSemiring lw, fourierSemiring rw)
    {
        if (lw == fourierSemiring(0, 1))
            return std::make_tuple(rw, fourierSemiring(0, 1), (rw == fourierSemiring(0, 1)) ? fourierSemiring(0, 1) : fourierSemiring(1, 1)); 
        if (rw == fourierSemiring(0, 1))
            return std::make_tuple(lw, fourierSemiring(1, 1), fourierSemiring(0, 1));
        return std::make_tuple(lw, fourierSemiring(1, 1), rw/lw);
    }

    template <>
    fourierSemiring computeComposition<fourierSemiring, std::multiplies<fourierSemiring>>(fourierSemiring c, fourierSemiring w)
    {
        return c * w;
    }

    template <>
    long double computeProbabilityFromAmplitude<fourierSemiring,std::multiplies<fourierSemiring>>(fourierSemiring w)
    {
        abort();
    }

}