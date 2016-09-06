#include <tomfloat.h>

int mpf_atanh(mp_float * a, mp_float * b)
{
    mp_float one, x, t;
    int err;
    /*
     * Placeholder (for now) but log() is quite fast
     * 
     * 1      (  1 + x  )
     * atanh(x) = -  log( ------- )
     * 2     (  1 - x  )
     * 
     * The series at mpf_kernel_atan() is able to evaluate atanh, too.
     * TODO: evaluate which is faster and when.
     */
    err = MP_OKAY;

    if (mpf_iszero(a)) {
	return mpf_copy(a, b);
    }

    if ((err = mpf_init_multi(a->radix, &one, &x, &t, NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_copy(a, &x)) != MP_OKAY) {
	goto _ERR;
    }
    x.mantissa.sign = MP_ZPOS;

    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }
    // TODO: do it within the limits of EPS
    if (mpf_cmp(a, &one) == MP_EQ) {
	// pole
	if ((err = mpf_const_inf(&x, a->mantissa.sign)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_copy(&x, b)) != MP_OKAY) {
	    goto _ERR;
	}
	goto _ERR;
    }

    if (mpf_isinf(a) || mpf_cmp(a, &one) == MP_GT) {
	// E_DOM
	if ((err = mpf_const_nan(&x)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_copy(&x, b)) != MP_OKAY) {
	    goto _ERR;
	}
	goto _ERR;
    }

    x.mantissa.sign = a->mantissa.sign;

    if ((err = mpf_add(&one, &x, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sub(&one, &x, &x)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_div(&t, &x, &x)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_ln(&x, &x)) != MP_OKAY) {
	goto _ERR;
    }
    x.exp -= 1;
    if ((err = mpf_copy(&x, b)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mpf_clear_multi(&one, &x, &t, NULL);
    return err;
}

