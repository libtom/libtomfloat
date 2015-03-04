#include <tomfloat.h>

/* IEEE-754 sec. 3.4: quiet NaN */
int mpf_const_nan(mp_float * a)
{
    long r, ilog2;
    int err;
    // IEEE-754 says in sec. 6.2.1: set most significant bit to one
    // We do something similar but more brutal and set the most
    // siginificant _limb_ to all-ones; easy to find, easy to check
    ilog2 = 0;
    if (a->mantissa.used == 0) {
	if ((err = mpf_const_d(a, 1)) != MP_OKAY) {
	    return err;
	}
    }
    r = a->mantissa.dp[a->mantissa.used - 1];
    while (r >>= 1) {
	ilog2++;
    }
    a->mantissa.dp[a->mantissa.used] = (1 << (ilog2 + 1)) - 1;
    // set exponent to all-ones
    a->exp = 1 << ((sizeof(long) * CHAR_BIT) - 1);
    // keep sign as it is and do not normalize
    return MP_OKAY;
}

