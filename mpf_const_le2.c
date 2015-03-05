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

static mp_float mpf_le2;

static long mpf_le2_precision;

int mpf_const_le2(mp_float * a)
{
    int err;
    long eps;

    err = MP_OKAY;

    if (mpf_le2_precision > 0 && a == NULL) {
	mpf_clear(&mpf_le2);
	mpf_le2_precision = 0;
	return err;
    }
    if (mpf_le2_precision >= a->radix) {
        eps = a->radix;
	if ((err = mpf_copy(&mpf_le2, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, eps);
    } else {
	if (mpf_le2_precision == 0) {
	    if ((err = mpf_init(&mpf_le2, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = mpf_const_ln_d(&mpf_le2, 2)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_copy(&mpf_le2, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}

		/* log_e 2       */

