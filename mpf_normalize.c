/* LibTomFloat, multiple-precision floating-point library
 *
 * LibTomFloat is a library that provides multiple-precision
 * floating-point artihmetic as well as trigonometric functionality.
 *
 * This library requires the public domain LibTomMath to be installed.
 * 
 * This library is free for all purposes without any express
 * gurantee it works
 *
 * Tom St Denis, tomstdenis@iahu.ca, http://float.libtomcrypt.org
 */
#include <tomfloat.h>

int mpf_normalize(mp_float * a)
{
    long cb, diff;
    int err;
    mp_digit c;

    /* sanity */
    if (a->radix < 2) {
	return MP_VAL;
    }
  loop:
    cb = mp_count_bits(&(a->mantissa));
    if (cb > a->radix) {
	diff = cb - a->radix;
	a->exp += diff;

	/* round it, add 1 after shift if diff-1'th bit is 1 */
	c = a->mantissa.dp[diff / DIGIT_BIT] & (1U << (diff % DIGIT_BIT));
	if ((err =
	     mp_div_2d(&(a->mantissa), diff, &(a->mantissa),
		       NULL)) != MP_OKAY) {
	    return err;
	}

	if (c != 0) {
	    if((err = mp_add_d(&(a->mantissa), 1, &(a->mantissa))) != MP_OKAY){
	        return err;
	    }
	    // in case of a carry: shift one right; rinse and repeat
	    if (mp_count_bits(&(a->mantissa)) > cb) {
		if ((err =
		     mp_div_2d(&(a->mantissa), 1, &(a->mantissa),
			       NULL)) != MP_OKAY) {
		    return err;
		}
                a->exp += 1;
                goto loop;
	    } else {
	        return MP_OKAY;
	    }
	} else {
	    return MP_OKAY;
	}
    } else if (cb < a->radix) {
	if (mp_iszero(&(a->mantissa)) == MP_YES) {
	    return mpf_const_0(a);
	} else {
	    diff = a->radix - cb;
	    a->exp -= diff;
	    return mp_mul_2d(&(a->mantissa), diff, &(a->mantissa));
	}
    }
    return MP_OKAY;
}

