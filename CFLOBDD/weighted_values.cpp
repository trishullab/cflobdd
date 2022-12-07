
#include "weighted_values.h"

namespace CFL_OBDD {

    template <>
    FLOAT_BOOST getIdentityValue<FLOAT_BOOST, std::multiplies<FLOAT_BOOST>>()
    {
        return 1.0;
    }

    template <>
    FLOAT_BOOST getAnnhilatorValue<FLOAT_BOOST, std::multiplies<FLOAT_BOOST>>()
    {
        return 0.0;
    }

    template <>
    std::tuple<FLOAT_BOOST, FLOAT_BOOST, FLOAT_BOOST> computeInverseValue<FLOAT_BOOST, std::multiplies<FLOAT_BOOST>>(FLOAT_BOOST lw, FLOAT_BOOST rw)
    {
        if (lw == 0.0)
            return std::make_tuple(rw, 0.0, 1.0);
        return std::make_tuple(lw, 1.0, rw/lw);
    }

    template <>
    FLOAT_BOOST computeComposition<FLOAT_BOOST, std::multiplies<FLOAT_BOOST>>(FLOAT_BOOST c, FLOAT_BOOST w)
    {
        return c * w;
    }


}