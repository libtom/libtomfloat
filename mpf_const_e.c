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

static mp_float mpf_e;

static long mpf_e_precision;

int mpf_const_e(mp_float * a)
{
    int err;

    err = MP_OKAY;

    if (mpf_e_precision > 0 && a == NULL) {
	mpf_clear(&mpf_e);
	mpf_e_precision = 0;
	return err;
    }
    if (mpf_e_precision >= a->radix) {
	if ((err = mpf_copy(&mpf_e, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, a->radix);
    } else {
	if (mpf_e_precision == 0) {
	    if ((err = mpf_init(&mpf_e, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = mpf_const_d(&mpf_e, 1)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_exp(&mpf_e, &mpf_e)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_copy(&mpf_e, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}

