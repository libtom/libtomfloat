#include <tomfloat.h>

/* Argument reduction for sine, cosine and tangent to x <= pi/4 */

int mpf_trig_arg_reduct(mp_float * a, mp_float * b, int *k)
{
    int err, sign;
    mp_float pi, pihalf, piquart, r, x, three, one, K;
    long size, oldeps, eps, oldprec, newprec;

    if (mpf_iszero(a)) {
	*k = 0;
	return mpf_copy(a, b);
    }
    oldeps = a->radix;
    eps = oldeps + 10;
    err = MP_OKAY;
    if ((err =
	 mpf_init_multi(eps, &pi, &pihalf, &piquart, &r, &x, &three, &one, &K,
			NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_abs(a, &x)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&x, eps)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_const_d(&three, 3)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_pi(&pi)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_pi(&pihalf)) != MP_OKAY) {
	goto _ERR;
    }
    pihalf.exp -= 1;
    if ((err = mpf_const_pi(&piquart)) != MP_OKAY) {
	goto _ERR;
    }
    piquart.exp -= 1;

    // nothing to do if it is already small enough
    if (mpf_cmp(&x, &piquart) != MP_GT) {
	if ((err = mpf_copy(a, b)) != MP_OKAY) {
	    goto _ERR;
	}
	*k = 0;
	goto _ERR;
    }
    // it starts to get tricky for x < 3pi/4 especially around pi/2
    // but not for the reduction part
    if ((err = mpf_mul(&three, &piquart, &three)) != MP_OKAY) {
	goto _ERR;
    }

    if (mpf_cmp(&x, &three) == MP_LT) {
	if ((err = mpf_sub(&x, &pihalf, &x)) != MP_OKAY) {
	    goto _ERR;
	}
	x.mantissa.sign = a->mantissa.sign;
	if ((err = mpf_copy(&x, b)) != MP_OKAY) {
	    goto _ERR;
	}
	*k = 1;
	goto _ERR;
    }

    sign = a->mantissa.sign;
    // size of integer part in bits
    size = a->radix + a->exp;

    if (size < 0) {
	size = 0;
    }
    // Reduction must be done in precision
    //     work_precision_(base 2) +  log_2(x)
    // if we have an integer part
    // But the main problem with such large numbers lies in the loss of accuracy
    // in the original number. That needs to get caught in the calling function.
    oldprec = eps;
    if (size > 0) {
	newprec = oldprec + size + 3;
    } else {
	newprec = oldprec + 3;
    }
    // we need as many digits of pi as there are digits in the integer part of the
    // number, so for e.g.: 1e308 we need 308 decimal digits of pi.
    // Computing that much may take a while.

    // Compute remainder of x/Pi/2
    if ((err = mpf_normalize_to(&pi, newprec)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_pi(&pi)) != MP_OKAY) {
	goto _ERR;
    }
    // k = round(x * 2/Pi)
    if ((err = mpf_normalize_to(&x, newprec)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_normalize_to(&K, newprec)) != MP_OKAY) {
	goto _ERR;
    }
    // multiply by two
    x.exp += 1;
    if ((err = mpf_div(&x, &pi, &K)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_round(&K, &K)) != MP_OKAY) {
	goto _ERR;
    }
    // undo multiplying by two
    x.exp -= 1;

    // r = x - k * Pi/2
    if ((err = mpf_normalize_to(&r, newprec)) != MP_OKAY) {
	goto _ERR;
    }
    pi.exp -= 1;
    if ((err = mpf_mul(&K, &pi, &pi)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sub(&x, &pi, &r)) != MP_OKAY) {
	goto _ERR;
    }
    // TODO: check for size of "size" and don't delete pi if "size" is moderate
    // mpf_const_pi(NULL);

    if ((err =
	 mp_div_2d(&K.mantissa, abs(K.exp), &K.mantissa, NULL)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&r, oldprec)) != MP_OKAY) {
	goto _ERR;
    }
    // we need the last two bits only
    *k = K.mantissa.dp[0];
    r.mantissa.sign = a->mantissa.sign;
    mpf_copy(&r, b);

    err = MP_OKAY;

  _ERR:
    mpf_clear_multi(&pi, &pihalf, &piquart, &r, &x, &three, &one, &K, NULL);
    return err;
}

