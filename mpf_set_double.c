#include "tomfloat.h"
#include <math.h>

#if __STDC_VERSION__ >= 199901L
#   include <fenv.h>
#endif
#include <float.h>
// set the double to the value of the mp_float if possible
// this method relies on the IEEE-754 complience of the libc!
int mpf_set_double(mp_float * a, double *d)
{
    mp_float cp;
    int err;
    //mp_digit high,low;
    int digits_needed, i, limit;
    double multiplier;


    if ((err = mpf_init_copy(a, &cp)) != MP_OKAY) {
	return err;
    }
    // cheap check for boundaries
    if ((err = mpf_normalize_to(&cp, DBL_MANT_DIG)) != MP_OKAY) {
	mpf_clear(&cp);
	return err;
    }
#if __STDC_VERSION__ >= 199901L
    if (a->exp < -1022) {
	*d = -INFINITY;
	goto END;
    } else if (a->exp > 1023) {
	*d = INFINITY;
	goto END;
    }
#else
    if (a->exp < -1022) {
	*d = -HUGE_VAL;
	goto END;
    } else if (a->exp > 1023) {
	*d = HUGE_VAL;
	goto END;
    }
#endif
    else if (mpf_isnan(a)) {
	*d = sqrt(-1);
	goto END;
    } else if (mpf_iszero(a)) {
	*d = (a->mantissa.sign == MP_NEG) ? -0.0 : 0.0;
	goto END;
    }
    // TODO: feed the double directly (but beware of endianess!)

    // Shamelessly ported from the libbtommath function mp_get_double

    digits_needed = ((DBL_MANT_DIG + MP_DIGIT_BIT) / MP_DIGIT_BIT) + 1;

    if ((err = mpf_copy(a, &cp)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_normalize_to(&cp, DBL_MANT_DIG + MP_DIGIT_BIT)) != MP_OKAY) {
	mpf_clear(&cp);
	return err;
    }
    multiplier = (double) (MP_MASK + 1);
    *d = 0.0;
    /* Could be assumed, couldn't it? */
    mp_clamp(&cp.mantissa);
    i = cp.mantissa.used;
    limit = (i <= digits_needed) ? 0 : i - digits_needed;

    while (i-- > limit) {
	*d += cp.mantissa.dp[i];
	*d *= multiplier;
    }

    if (cp.mantissa.sign == MP_NEG) {
	*d *= -1.0;
    }

    *d *= pow(2.0, i * DIGIT_BIT);
    // The content of the mantissa is in *d now, apply exponent
    *d *= pow(2.0, (double) (cp.exp));

    /* Handle overflow */
#if __STDC_VERSION__ >= 199901L
    if (*d == INFINITY) {
	feraiseexcept(FE_OVERFLOW);
#else
    if (*d == HUGE_VAL) {
#endif
	return MP_RANGE;
    }
  END:
    mpf_clear(&cp);
    return MP_OKAY;
}

