#include <tomfloat.h>
/* IEEE-754 sec. 3.4: quiet NaN */
int mpf_const_nan(mp_float * a)
{
    // IEEE-754 says in sec. 6.2.1: set most significant bit to one
    // We do something similar but more brutal and set the most
    // siginificant _limb_ to all-ones; easy to find, easy to check
    a->mantissa.dp[a->mantissa.used - 1] = (mp_digit) (~0U);
    // set exponent to all-ones
    a->exp = ~0L;
    // keep sign as it is and do not normalize
    return MP_OKAY;
}

