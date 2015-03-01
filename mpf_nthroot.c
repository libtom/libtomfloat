#include <tomfloat.h>
static int mpf_print(mp_float * a)
{
    // one for EOS and one for sign
    mp_fput(&(a->mantissa), 10, stdout);
    printf(" * 2^%ld * 1.0\n", a->exp);

    return MP_OKAY;
}

int mpf_nthroot(mp_float * a, long n, mp_float * b)
{
    int err, sign;
    mp_float one, t2, t3, invb, x0;
    mp_float ret;
    mp_float frac;
    long oldeps, nm1, loops, expnt;
    double d;

    // TODO: IEEE-754 has a quite complicated system here, adapt
    if (mpf_iszero(a)) {
	if (n > 0) {
	    mpf_const_0(b);
	    return MP_OKAY;
	} else {
	    mpf_const_nan(b);
	    return MP_OKAY;
	}
    }
    if (mpf_isnan(a)) {
	mpf_const_nan(b);
	return MP_OKAY;
    }
    if (mpf_isinf(a)) {
	mpf_const_inf(b, MP_ZPOS);
	return MP_OKAY;
    }
    if (n < 0) {
	// 1/(x^(1/-n))
	n = -n;
	mpf_init_copy(&ret, a);
	mpf_inv(&ret, &ret);
	mpf_nthroot(&ret, n, b);
	mpf_clear(&ret);
	return MP_OKAY;
    }
    if (n == 1) {
	// a^(1/1)
	mpf_copy(a, b);
	return MP_OKAY;
    }
    /* input must be positive if b is even */
    if ((n & 1) == 0 && a->mantissa.sign == MP_NEG) {
	return MP_VAL;
    }
    if (n == 2) {
	return mpf_sqrt(a, b);
    }
    mpf_init(&one, a->radix);
    mpf_const_d(&one, 1);

    if (mpf_cmp(a, &one) == MP_EQ) {
	if (n == 0) {
	    mpf_const_nan(b);
	    mpf_clear(&one);
	    return MP_VAL;
	} else {
	    mpf_const_0(b);
	    mpf_clear(&one);
	    return MP_OKAY;
	}
    }

    mpf_clear(&one);

    oldeps = a->radix;
    mpf_init_copy(a, &ret);
    mpf_normalize_to(&ret, oldeps + 10);

    sign = a->mantissa.sign;
    ret.mantissa.sign = MP_ZPOS;

    /*
     * Compute initial value (assuming n is positive)
     *
     * x^n = exp(log(x) * n)
     * x^n = 2^(log_2(x) * n)
     *
     */
    mpf_init(&t2, ret.radix);
    mpf_init(&frac, ret.radix);
    mpf_frexp(&ret, &frac, &expnt);
    expnt = expnt / n;
    mpf_set_double(&frac, &d);
    d = pow(d, (double) (n));
    mpf_get_double(d, &frac);
    mpf_ldexp(&frac, expnt, &t2);
    // we should have a good initial value in t2 now.
    // It can be done better if we check if the input has a magnitude that
    // allows to do it all in double precision which can safe one round.
    mpf_init(&invb, ret.radix);
    mpf_const_d(&invb, n);
    mpf_inv(&invb, &invb);
    nm1 = n - 1;
    // safe guard for the worst case
    loops = ret.radix + 1;
    mpf_init_multi(ret.radix, &t3, &x0, NULL);
    do {
	// x0 = t2
	mpf_copy(&t2, &x0);
	// t3 = t2^nm1
	mpf_pow_d(&t2, nm1, &t3);
	// t3 = a / t3
	mpf_div(&ret, &t3, &t3);
	// t3 = t3 - t2
	mpf_sub(&t3, &t2, &t3);
	// t3 = t3 * invb
	mpf_mul(&t3, &invb, &t3);
	// t2 = t2 + t3
	mpf_add(&t2, &t3, &t2);
	if (mpf_iszero(&t2)) {
	    break;
	}
	if (loops-- == 0) {
	    fprintf(stderr, "mp_float nthroot did not converge in %ld rounds\n",
		    (ret.radix + 1));
	    err = MP_VAL;
	    goto _ERR;
	}
    } while (mpf_cmp(&t2, &x0) != MP_EQ);
//#ifdef DEBUG
    fprintf(stderr, "nth-root loops =  %ld\n", (ret.radix + 1) - loops);
//#endif
    /* set the sign of the result */
    t2.mantissa.sign = sign;
    mpf_normalize_to(&t2, a->radix);
    mpf_copy(&t2, b);
    err = MP_OKAY;

  _ERR:
    mpf_clear_multi(&t2, &t3, &x0, &invb, &ret, NULL);
    return err;
}

