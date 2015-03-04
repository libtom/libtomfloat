#include <tomfloat.h>
/* IEEE-754 sec. 3.4: +/- infinity */
int mpf_const_inf(mp_float * a, int sign)
{
    // set the most siginificant _limb_ to zero; easy to find,
    // easy to check
    if (a->mantissa.used == 0) {
	mpf_const_d(a, 1);
    }
    a->mantissa.dp[a->mantissa.used] = (mp_digit) (0);
    // set exponent to all-ones
    a->exp = 1 << ((sizeof(long) * CHAR_BIT) - 1);
    // keep sign as it is and do not normalize
    return MP_OKAY;
}

