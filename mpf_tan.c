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

/* tan = sin/cos, prolly fix this up later */
/* CZ: Still tan = sin/cos, Tom,  but only very high precisions shall fear it */
int mpf_tan(mp_float * a, mp_float * b)
{
    mp_float x;
    int err, k, sign;
    long size;

    if (mpf_iszero(a)) {
	return mpf_const_0(b);
    }

    size = a->radix + a->exp;
    if (size > a->radix) {
	// raise invalid (TLOSS in SysV)?
	// feraiseexcept(FE_INVALID)
	// mpf_errno = E_DOM
	// or inexact?
	if ((err = mpf_const_nan(b)) != MP_OKAY) {
	    return err;
	}
	return MP_VAL;
    }
    sign = a->mantissa.sign;
    if ((err = mpf_init(&x, a->radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_trig_arg_reduct(a, &x, &k)) != MP_OKAY) {
	goto _ERR;
    }
    k = k % 4;

    if ((err = mpf_sincos(&x, &x, 1, 1, 0)) != MP_OKAY) {
	goto _ERR;
    }

    switch (k) {
    case 0:
	break;
    case 1:
	if ((err = mpf_inv(&x, &x)) != MP_OKAY) {
	    goto _ERR;
	}
	sign = (sign == MP_NEG) ? MP_ZPOS : MP_NEG;
	break;
    case 2:
	break;
    default:
	if ((err = mpf_inv(&x, &x)) != MP_OKAY) {
	    goto _ERR;
	}
	sign = (sign == MP_NEG) ? MP_ZPOS : MP_NEG;
	break;
    }
    x.mantissa.sign = sign;
    if ((err = mpf_copy(&x, b)) != MP_OKAY) {
	goto _ERR;
    }
  _ERR:
    mpf_clear(&x);
    return err;
}

