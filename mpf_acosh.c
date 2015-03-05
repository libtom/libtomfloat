#include <tomfloat.h>

int mpf_acosh(mp_float * a, mp_float * b)
{
    int err, sign;
    mp_float one, t, t2;

    err = MP_OKAY;
    sign = MP_ZPOS;
    if ((err = mpf_init_multi(a->radix, &one, &t, NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_sqr(a, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sub(&t, &one, &t)) != MP_OKAY) {
	goto _ERR;
    }
    // simplest way to check for domain errors, although not the cheapest
    if (t.mantissa.sign == MP_NEG) {
	// E_DOM
	if ((err = mpf_const_nan(b)) != MP_OKAY) {
	    goto _ERR;
	}
	err = MP_VAL;
	goto _ERR;
    }
    // acosh(1) = 0
    if (mpf_iszero(&t)) {
	if ((err = mpf_const_d(b, 0)) != MP_OKAY) {
	    goto _ERR;
	}
	goto _ERR;
    }

    if ((err = mpf_sqrt(&t, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(a, &t, b)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_ln(b, b)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mpf_clear_multi(&one, &t, NULL);
    return err;
}

