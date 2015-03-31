#include <tomfloat.h>

/* IEEE-754 sec. 3.4: quiet NaN */
int mpf_const_nan(mp_float * a)
{
    int err;
    // IEEE-754 says in sec. 6.2.1: set most significant bit to one
    // We do something similar but more brutal and set the most
    // siginificant _limb_ to one; easy to find, easy to check
    if (a->mantissa.used == 0) {
	if ((err = mpf_const_d(a, 1)) != MP_OKAY) {
	    return err;
	}
    }
    a->mantissa.dp[a->mantissa.used] = (mp_digit)(1);
    // set exponent to all-ones
    a->exp = LONG_MAX;
    // keep sign as it is and do not normalize
    return MP_OKAY;
}

