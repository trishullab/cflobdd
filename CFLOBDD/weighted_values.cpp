
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


}