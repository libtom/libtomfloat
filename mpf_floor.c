#include <tomfloat.h>
/* floor function, rounds away from zero to +/- infinity */
int mpf_floor(mp_float * a, mp_float * b)
{
    mp_int q, r;
    mp_float ret;
    int err;

    err = MP_OKAY;
    // integer
    if (a->exp > 0) {
	return mpf_copy(a, b);
    }
    // fraction, return zero
    if (a->exp <= -a->radix) {
	return mpf_const_0(b);
    }

    if ((err = mp_init_multi(&q, &r, NULL)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_init(&ret, a->radix)) != MP_OKAY) {
	mp_clear_multi(&q, &r, NULL);
	return err;
    }
    //mixed/integer
    if ((err = mp_div_2d(&a->mantissa, abs(a->exp), &q, &r)) != MP_OKAY) {
	goto _ERR;
    }
    // no remainder: it is an integer
    if (mp_iszero(&r)) {
	if ((err = mpf_copy(a, b)) != MP_OKAY) {
	    goto _ERR;
	}
	goto _ERR;
    }
    // should not happen
    if (mp_iszero(&q)) {
	err = mpf_const_0(b);
	goto _ERR;
    }
    if ((err = mp_copy(&q, &ret.mantissa)) != MP_OKAY) {
	goto _ERR;
    }
    ret.exp = 0;
    if ((err = mpf_normalize(&ret)) != MP_OKAY) {
	goto _ERR;
    }
    // round to -infinity, e.g.: floor(-1.9) = -2
    if (a->mantissa.sign == MP_NEG) {
	if ((err = mpf_sub_d(&ret, 1, &ret)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_copy(&ret, b)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mp_clear_multi(&q, &r,NULL);
    mpf_clear(&ret);
    return err;
}

