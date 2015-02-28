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

#include <tommath.h>

/* compute 1/x */
int mpf_inv(mp_float * a, mp_float * b)
{
    int err, sign;
    double d;
    long expnt, oldeps, nloops, maxrounds;
    mp_float frac, one, x0, xn, A, hn;

    /* get sign of input */
    sign = a->mantissa.sign;

    /* force to positive */
    a->mantissa.sign = MP_ZPOS;

    oldeps = a->radix;

    if ((err = mpf_init(&frac, a->radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_frexp(a, &frac, &expnt)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_set_double(&frac, &d)) != MP_OKAY) {
	return err;
    }
    d = 1.0 / d;
    if ((err = mpf_get_double(d, &frac)) != MP_OKAY) {
	return err;
    }
    // TODO: checks & balances
    expnt = -expnt + 1;
    if ((err = mpf_ldexp(&frac, expnt, &frac)) != MP_OKAY) {
	return err;
    }
    // TODO calculate guard digits more exactly and do it loop-wise
    if ((err = mpf_normalize_to(&frac, frac.radix + MP_DIGIT_BIT)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_init_copy(&frac, &xn)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_init(&x0, frac.radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_init(&one, frac.radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_const_d(&one, 1L)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_init(&hn, frac.radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_init_copy(a, &A)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_normalize_to(&A, frac.radix)) != MP_OKAY) {
	return err;
    }
    maxrounds = A.radix;
    nloops = 0L;
    do {
	if ((err = mpf_copy(&xn, &x0)) != MP_OKAY) {
	    return err;
	}
	// hn = 1 - (A * xn);
	if ((err = mpf_mul(&A, &xn, &hn)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_sub(&one, &hn, &hn)) != MP_OKAY) {
	    return err;
	}
	// TODO: check against EPS instead
	if (mpf_iszero(&hn)) {
	    break;
	}
	// xn = xn + (xn * hn);
	if ((err = mpf_mul(&xn, &hn, &hn)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_add(&xn, &hn, &xn)) != MP_OKAY) {
	    return err;
	}
	nloops++;
	if (nloops >= maxrounds) {
	    // it might be a bug elsewhere, please report
	    fprintf(stderr, "mpf_inv did not converge in %ld rounds", nloops);
	    return MP_RANGE;
	}
    } while (mpf_cmp(&x0, &xn) != MP_EQ);
    if ((err = mpf_normalize_to(&xn, oldeps)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_copy(&xn, b)) != MP_OKAY) {
	return err;
    }

    /* now restore the sign */
    b->mantissa.sign = sign;

    return err;
}

