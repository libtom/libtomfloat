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
// http://www.jjj.de/fxt/#fxtbook
static int machin_ln_10(mp_float * a)
{
    int err;
    long oldeps, eps;
    mp_int c1, c2, c3, a1, a2, a3, sum, zero, t1, t2, P, Q, R, EPS;

    oldeps = a->radix;
    eps = (oldeps + MP_DIGIT_BIT);

    // TODO: error calculation!

    err = MP_OKAY;
    if ((err =
	 mp_init_multi(&c1, &c2, &c3, &a1, &a2, &a3, &sum, &zero, &t1, &t2,
		       &EPS, &P, &Q, &R, NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mp_set_int(&EPS, (unsigned long) (eps))) != MP_OKAY) {
	goto _ERR;
    }
    // using acoth instead of atanh to make things simpler 
    // atanh(x) = acoth(1/x) for x > 1

    // coefficients due to P. Sebah "Machin like formulae for logarithm" (1997)
    // http://numbers.computation.free.fr/Constants/PiProgram/userconstants.html
    mp_set(&c1, 46);
    mp_set(&c2, 34);
    mp_set(&c3, 20);

    mp_set(&a1, 31);
    mp_set(&a2, 49);
    mp_set(&a3, 161); // http://oeis.org/A175607

    mp_zero(&sum);
    mp_zero(&zero);

    if ((err =
	 mp_acoth_binary_splitting(&a1, &zero, &EPS, &P, &Q, &R)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_add(&P, &Q, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul_2d(&t1, eps, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul(&Q, &a1, &t2)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_div(&t1, &t2, &t1, NULL)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul(&t1, &c1, &sum)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err =
	 mp_acoth_binary_splitting(&a2, &zero, &EPS, &P, &Q, &R)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_add(&P, &Q, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul_2d(&t1, eps, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul(&Q, &a2, &t2)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_div(&t1, &t2, &t1, NULL)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul(&t1, &c2, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_add(&sum, &t1, &sum)) != MP_OKAY) {
	goto _ERR;
    }


    if ((err =
	 mp_acoth_binary_splitting(&a3, &zero, &EPS, &P, &Q, &R)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_add(&P, &Q, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul_2d(&t1, eps, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul(&Q, &a3, &t2)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_div(&t1, &t2, &t1, NULL)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul(&t1, &c3, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_add(&sum, &t1, &sum)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_from_mp_int(&sum, a)) != MP_OKAY) {
	goto _ERR;
    }
    a->exp -= eps;

_ERR:
    mp_clear_multi(&c1, &c2, &c3, &a1, &a2, &a3, &sum, &zero, &t1, &t2, &EPS,
		   &P, &Q, &R, NULL);
    return err;
}

static mp_float mpf_le10;

static long mpf_le10_precision;

int mpf_const_le10(mp_float * a)
{
    int err;
    long eps;

    err = MP_OKAY;

    if (mpf_le10_precision > 0 && a == NULL) {
	mpf_clear(&mpf_le10);
	mpf_le10_precision = 0;
	return err;
    }
    if (mpf_le10_precision >= a->radix) {
	eps = a->radix;
	if ((err = mpf_copy(&mpf_le10, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, eps);
    } else {
	if (mpf_le10_precision == 0) {
	    if ((err = mpf_init(&mpf_le10, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = machin_ln_10(&mpf_le10)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_copy(&mpf_le10, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}

