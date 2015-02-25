#include <tomfloat.h>
/*
   mp_float to the power of a long. Will return the reciprocal if the exponent
   is negative.
*/
int mpf_pow_d(mp_float * a, long e, mp_float * c)
{
    mp_float ret, b;
    long ex, prec;
    int err = MP_OKAY;

    // x^1 = x
    if (e == 1) {
	// save one copying
	if (mpf_iszero(a)) {
	    return MP_OKAY;
	}
	if ((err = mpf_copy(a, c)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }
    // x^-1 = 1/x
    if (e == -1) {
	// TODO: raise DivisionByZero?
	if (mpf_iszero(a)) {
	    return MP_VAL;
	}
	if ((err = mpf_inv(a, c)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }
    // x^0 = 1, -x^0 = -1
    if (e == 0) {
	if ((err = mpf_const_d(c, 1)) != MP_OKAY) {
	    return err;
	}
	c->mantissa.sign = a->mantissa.sign;
	return MP_OKAY;
    }
    // some guard digits (one limb at least)
    // TODO: adjust according to size of exponent 
    prec = c->radix + MP_DIGIT_BIT;
    if ((err = mpf_init(&ret, prec)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_init_copy(a, &b)) != MP_OKAY) {
	mpf_clear(&ret);
	return err;
    }
    if ((err = mpf_normalize_to(&b, prec)) != MP_OKAY) {
	mpf_clear_multi(&ret,&b,NULL);
	return err;
    }
    if ((err = mpf_const_d(&ret, 1)) != MP_OKAY) {
	mpf_clear_multi(&ret, &b, NULL);
	return err;
    }
    // work on abs. values, do the reciprocal later
    if (e < 0) {
	ex = -e;
    } else {
	ex = e;
    }
    while (ex != 0) {
	if ((ex & 1) == 1) {
	    if ((err = mpf_mul(&ret, &b, &ret)) != MP_OKAY) {
		mpf_clear_multi(&ret, &b, NULL);
		return err;
	    }
	}
	if ((err = mpf_sqr(&b, &b)) != MP_OKAY) {
	    mpf_clear_multi(&ret, &b, NULL);
	    return err;
	};
	ex >>= 1;
    }
    if (e < 0) {
	if ((err = mpf_inv(&ret, &ret)) != MP_OKAY) {
	    mpf_clear_multi(&ret, &b, NULL);
	    return err;
	}
    }
    if ((err = mpf_normalize_to(&ret, c->radix)) != MP_OKAY) {
	mpf_clear_multi(&ret,&b,NULL);
	return err;
    }
    if ((err = mpf_copy(&ret, c)) != MP_OKAY) {
	mpf_clear_multi(&ret, &b, NULL);
	return err;
    }
    mpf_clear_multi(&ret, &b, NULL);
    return err;
}

