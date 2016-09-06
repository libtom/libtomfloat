#include <tomfloat.h>
/* IEEE-754 sec. 3.4: +/- infinity */
int mpf_const_inf(mp_float * a, int sign)
{
    int err;
    // set the most siginificant _limb_ to zero; easy to find,
    // easy to check
    if (a->mantissa.used == 0) {
	if ((err = mpf_const_d(a, 1)) != MP_OKAY) {
	    return err;
	}
    }
    a->mantissa.dp[a->mantissa.used - 1] = (mp_digit) (0);
    // set exponent to LONG_MAX
    a->exp = LONG_MAX;
    // keep sign as it is and do not normalize
    a->mantissa.sign = sign;
    return MP_OKAY;
}

