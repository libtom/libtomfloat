#include "tomfloat.h"

int mpf_ldexp(mp_float * a, long exp, mp_float * b)
{
    if (mpf_isinf(a) || mpf_isnan(a) || mpf_iszero(a)) {
	if ((err = mpf_copy(a, b)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }
    if ((err = mp_copy(&a->mantissa, &b->mantissa)) != MP_OKAY) {
	return err;
    }
    b->exp = exp - a->radix;
    b->radix = a->radix;
    return MP_OKAY;
}

