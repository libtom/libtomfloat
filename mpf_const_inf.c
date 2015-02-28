#include <tomfloat.h>
/* IEEE-754 sec. 3.4: +/- infinity */
int mpf_const_inf(mp_float * a, int sign)
{
    // set the most siginificant _limb_ to zero; easy to find,
    // easy to check
    a->mantissa.dp[a->mantissa.used - 1] = (mp_digit) (0U);
    // set exponent to all-ones
    a->exp = ~0L;
    a->mantissa.sign = sign;
    return MP_OKAY;
}
