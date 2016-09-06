#include "tomfloat.h"

int mpf_frexp(mp_float * a, mp_float * b, long *exp)
{
    int err;
    if (mpf_isinf(a) || mpf_isnan(a) || mpf_iszero(a)) {
	if ((err = mpf_copy(a, b)) != MP_OKAY) {
	    return err;
	}
	*exp = 0;
	return MP_OKAY;
    }
    *exp = a->exp + a->radix;
    b->exp = -a->radix;
    b->radix = a->radix;

    if ((err = mp_copy(&a->mantissa, &b->mantissa)) != MP_OKAY) {
	return err;
    }
    return MP_OKAY;
}

