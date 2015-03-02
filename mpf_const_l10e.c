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

static mp_float mpf_l10e;

static long mpf_l10e_precision;

int  mpf_const_l10e(mp_float *a)
{
    int err;

    err = MP_OKAY;

    if (mpf_l10e_precision > 0 && a == NULL) {
	mpf_clear(&mpf_l10e);
	mpf_l10e_precision = 0;
	return err;
    }
    if (mpf_l10e_precision >= a->radix) {
	if ((err = mpf_copy(&mpf_l10e, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, a->radix);
    } else {
	if (mpf_l10e_precision == 0) {
	    if ((err = mpf_init(&mpf_l10e, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = mpf_const_ln_d(&mpf_l10e, 10)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_inv(&mpf_l10e, &mpf_l10e)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_copy(&mpf_l10e, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}
                /* log_10 e == 1/ln(10)      */

