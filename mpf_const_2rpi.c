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

static mp_float mpf_2rpi;

static long mpf_2rpi_precision;

int mpf_const_2rpi(mp_float * a)
{

    int err;
    long eps;

    err = MP_OKAY;

    if (mpf_2rpi_precision > 0 && a == NULL) {
	mpf_clear(&mpf_2rpi);
	mpf_2rpi_precision = 0;
	return err;
    }
    if (mpf_2rpi_precision >= a->radix) {
        eps = a->radix;
	if ((err = mpf_copy(&mpf_2rpi, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, eps);
    } else {
	if (mpf_2rpi_precision == 0) {
	    if ((err = mpf_init(&mpf_2rpi, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = mpf_const_pi(&mpf_2rpi)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_invsqrt(&mpf_2rpi, &mpf_2rpi)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_mul_d(&mpf_2rpi, 2, &mpf_2rpi)) != MP_OKAY) {
	    return err;
	}
	mpf_2rpi.exp += 1;
	if ((err = mpf_copy(&mpf_2rpi, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}

		/* 2/sqrt(Pi)     */

