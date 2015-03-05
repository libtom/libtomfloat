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

/*
    Argument reduction for atan(x) with |x| > 1

      atan(x) = Pi/2 - atan(1/x)       \; x > 0
      atan(x) = -1 * (Pi/2 + atan(1/x) \; x < 0

    Argument reduction for atan(x) with |x| < 1

      atan(x) = atan( 1 ) + atan( (t-1)/(1+t) ) \; good for .5<x<1

    Sign:

      atan(x) = -atan(-x)

    Fixed values:
      atan(-1) = -Pi/4
      atan(1) = Pi/4
      atan(0) = 0
      atan(Inf) = P/2;
*/

int mpf_atan(mp_float * a, mp_float * b)
{
    mp_float one, x, pi, t1, t2, EPS;
    long eps, oldeps;
    int err, sign;

    if (mpf_isnan(a) || mpf_iszero(a)) {
	return mpf_copy(a, b);
    }

    err = MP_OKAY;
    oldeps = a->radix;
    // TODO: series evaluation has guard digits already
    eps = oldeps + 3;

    if ((err = mpf_init_multi(oldeps, &one, &x, NULL)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_init_multi(eps, &pi, &t1, &t2, &EPS, NULL)) != MP_OKAY) {
	mpf_clear_multi(&one, &x, NULL);
	return err;
    }

    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_pi(&pi)) != MP_OKAY) {
	goto _ERR;
    }

    sign = a->mantissa.sign;

    if ((err = mpf_abs(a, &x)) != MP_OKAY) {
	goto _ERR;
    }
    // we have to do this check in the precision of the input
    // TODO: use EPS instead
    if (mpf_cmp(a, &one) == MP_EQ) {
	pi.exp -= 2;
	if (sign == MP_NEG) {
	    pi.mantissa.sign = MP_NEG;
	}
	if ((err = mpf_normalize_to(&pi, oldeps)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_copy(&pi, b)) != MP_OKAY) {
	    goto _ERR;
	}
	goto _ERR;
    }

    if ((err = mpf_normalize_to(&x, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&one, eps)) != MP_OKAY) {
	goto _ERR;
    }
    // a new one; would get hit by rounding errors otherwise
    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }
    // The smaller x is, the better the series works
    // but foremost because of the instability near one.
    if (mpf_cmp(&x, &one) == MP_LT) {
	// atan(x) = atan( 1 ) + atan( (x-1)/(1+x) ) \; .5<x<1

	// Actually, the exact limit is sqrt(2)-1, the positive solution of
	// x + (x - 1)/(1 + x)
	// but .5 is easier to handle
	one.exp -= 1;
	if (mpf_cmp(&x, &one) != MP_LT) {
	    // atan(1) = pi/4
	    pi.exp -= 2;
	    one.exp += 1;
	    // x = (x - 1)/(1 + x)
	    if ((err = mpf_sub(&x, &one, &t1)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_add(&one, &x, &t2)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_div(&t1, &t2, &x)) != MP_OKAY) {
		goto _ERR;
	    }
	    // x will be negative, make it positive again
	    x.mantissa.sign = MP_ZPOS;
	    // evaluate series for x and subtract from pi/4
	    if ((err = mpf_kernel_atan(&x, &x, 0)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_sub(&pi, &x, &x)) != MP_OKAY) {
		goto _ERR;
	    }
	    if (sign == MP_NEG) {
		x.mantissa.sign = MP_NEG;
	    }
	    if ((err = mpf_normalize_to(&x, oldeps)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_copy(&x, b)) != MP_OKAY) {
		goto _ERR;
	    }
	    goto _ERR;
	} else {
	    if ((err = mpf_kernel_atan(&x, &x, 0)) != MP_OKAY) {
		goto _ERR;
	    }
	    if (sign == MP_NEG) {
		x.mantissa.sign = MP_NEG;
	    }
	    if ((err = mpf_normalize_to(&x, oldeps)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_copy(&x, b)) != MP_OKAY) {
		goto _ERR;
	    }
	    goto _ERR;
	}
    } else {
	// atan(x) = Pi/2 - atan(1/x)       \; x > 0
	// atan(x) = -1 * (Pi/2 + atan(1/x) \; x < 0
	// the latter is not necessary because of
	// atan(x) = -atan(-x)

	// The inverse can get very small if x is very large, obviously.
	// Cutoff depends on actual precision and absolute size of x, hence the
	// check against EPS.

	if ((err = mpf_inv(&x, &x)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	    goto _ERR;
	}
	if (mpf_cmp(&x, &EPS) == MP_GT) {
	    pi.exp -= 1;
	    if ((err = mpf_kernel_atan(&x, &x, 0)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_sub(&pi, &x, &x)) != MP_OKAY) {
		goto _ERR;
	    }
	    if (sign == MP_NEG) {
		x.mantissa.sign = MP_NEG;
	    }
	    if ((err = mpf_normalize_to(&x, oldeps)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_copy(&x, b)) != MP_OKAY) {
		goto _ERR;
	    }
	    goto _ERR;
	} else {
	    // \lim x\to\infty tan^{-1} x = \frac{\pi}{2}
	    pi.exp -= 1;
	    if (sign == MP_NEG) {
		pi.mantissa.sign = MP_NEG;
	    }
	    if ((err = mpf_normalize_to(&pi, oldeps)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_copy(&pi, b)) != MP_OKAY) {
		goto _ERR;
	    }
	    goto _ERR;
	}
    }


  _ERR:
    mpf_clear_multi(&one, &x, &pi, &t1, &t2, &EPS, NULL);
    return err;
}

