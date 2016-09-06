#include <tomfloat.h>

int mpf_asinh(mp_float * a, mp_float * b)
{
    int err;
    mp_float one, t;

    err = MP_OKAY;
    if ((err = mpf_init_multi(a->radix, &one, &t, NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_sqr(a, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(&t, &one, &t)) != MP_OKAY) {
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

