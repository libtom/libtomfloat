#include <tomfloat.h>

/* 
                    (         x          )
    asin(x) = 2 atan( ------------------ )
                    (  1 + sqrt(1 -x^2)  )

*/
int mpf_asin(mp_float * a, mp_float * b)
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
    if ((err = mpf_sub(&one, &t, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sqrt(&t, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(&one, &t, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_div(a, &t, &t)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_atan(&t, b)) != MP_OKAY) {
	goto _ERR;
    }

    b->exp += 1;
  _ERR:
    mpf_clear_multi(&one, &t, NULL);
    return err;
}

