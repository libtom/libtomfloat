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

int mpf_sqrt(mp_float * a, mp_float * b)
{
    int err;
    long oldeps, eps;
    mp_float A, B;

    oldeps = a->radix;
    eps = oldeps + 10;

    err = MP_OKAY;

    if ((err = mpf_init_multi(eps, &A, &B, NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&A, eps)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_invsqrt(&A, &B)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_mul(&A, &B, &B)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_normalize_to(&B, oldeps)) != MP_OKAY) {
	goto _ERR;
    }

    mpf_exch(&B, b);
  _ERR:
    mpf_clear_multi(&A, &B, NULL);
    return err;
}
