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
/* log base e (series for small precisions, AGM for the rest) */
int mpf_ln(mp_float * a, mp_float * b)
{
    mp_float ret, one, argred, frac, x, xc, exponent, t, diff, x0;
    long eps, oldeps, expnt, ar, n;
    int err;

    err = MP_OKAY;

    if (mpf_isnan(a)) {
	if ((err = mpf_copy(a, b)) != MP_OKAY) {
	    return err;
	}
	return MP_VAL;
    }

    if (mpf_iszero(a)) {
	if ((err = mpf_const_inf(b, MP_NEG)) != MP_OKAY) {
	    return err;
	}
	// TODO: raise FE_DIVBYZERO
	return MP_RANGE;
    }

    if (mpf_isinf(a)) {
	if ((err = mpf_const_inf(b, MP_ZPOS)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }

    if (a->mantissa.sign == MP_NEG) {
	if ((err = mpf_const_nan(b)) != MP_OKAY) {
	    return err;
	}
	b->mantissa.sign = MP_NEG;
	// TODO: raise FE_INVALID
	return MP_VAL;
    }

    if ((err = mpf_init(&one, a->radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	mpf_clear(&one);
	return err;
    }
    if (mpf_cmp(a, &one) == MP_EQ) {
	if ((err = mpf_const_0(b)) != MP_OKAY) {
	    return err;
	}
	mpf_clear(&one);
	return MP_OKAY;
    }
    mpf_clear(&one);
    // see mpf_global_variables.c (original was set at 4096)
    ar = MPF_LOG_AGM_REDUX_CUTOFF;

    oldeps = a->radix;
    // one limb extra precision
    eps = oldeps + MP_DIGIT_BIT;
    // Log by the AGM is faster already at some very low cutoff. Somebody
    // (Borwein?) said it would be faster at about 16 decimal digits precision.
    // The original value of the cutoff is set at 2,000 bits radix and
    // 100 bits absolute size (~10^30) but your version might have changed that
    // see mpf_global_variables.c for the exact values and adjust according to
    // your needs. Please contact the author if you did so succesfully.
    if (a->radix >= MPF_LOG_AGM_1_CUTOFF ||
        a->exp + a->radix >=  MPF_LOG_AGM_2_CUTOFF) {
        return mpf_ln_agm(a, b);
    }
    // TODO: make series an extra function (useful for e.g.: benchmarking)
    if ((err = mpf_init_copy(a, &x)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_normalize_to(&x, eps)) != MP_OKAY) {
	mpf_clear(&x);
	return err;
    }

    // argument reduction
    if ((err =
	 mpf_init_multi(eps, &argred, &frac, &xc, &ret, &exponent, &one, &t,
			&diff, &x0, NULL)) != MP_OKAY) {
	mpf_clear(&x);
	return err;
    }
    // seems small but more will need more guard bits, too
    if ((err = mpf_const_d(&argred, ar)) != MP_OKAY) {
	goto _ERR;
    }
    // work on mantissa which has the guaranteed magnitude .5<=x<1

    if ((err = mpf_frexp(&x, &frac, &expnt)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_nthroot(&frac, ar, &x)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sub(&x, &one, &x)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&x, &xc)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&x, &ret)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&exponent, expnt)) != MP_OKAY) {
	goto _ERR;
    }
    n = 2;
    do {
	// x0 = ret
	if ((err = mpf_copy(&ret, &x0)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = mpf_float(n)
	if ((err = mpf_const_d(&t, n)) != MP_OKAY) {
	    goto _ERR;
	}
	// n = 1/n
	if ((err = mpf_inv(&t, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	// x^(n-1) -> x^(n)
	if ((err = mpf_mul(&x, &xc, &x)) != MP_OKAY) {
	    goto _ERR;
	}
	// 1/n*x^(n)
	if ((err = mpf_mul(&t, &x, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	// series + 1/n*x^(n)
	if ((n & 1) == 1) {
	    if ((err = mpf_add(&ret, &t, &ret)) != MP_OKAY) {
		goto _ERR;
	    }
	} else {
	    if ((err = mpf_sub(&ret, &t, &ret)) != MP_OKAY) {
		goto _ERR;
	    }
	}
	if ((err = mpf_sub(&x0, &ret, &diff)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_abs(&diff, &diff)) != MP_OKAY) {
	    goto _ERR;
	}
	// TODO: check against EPS here, would make no sense otherwise
	if (mpf_iszero(&diff)) {
	    break;
	}
	n++;
    } while (mpf_cmp(&x0, &ret) != MP_EQ);

    if ((err = mpf_const_le2(&one)) != MP_OKAY) {
	goto _ERR;
    }
    // ret = ret * argred + (exponent * log(2))

    if ((err = mpf_mul(&exponent, &one, &exponent)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_mul(&ret, &argred, &ret)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(&ret, &exponent, &ret)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_normalize_to(&ret, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&ret, b)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mpf_clear_multi(&argred, &frac, &one, &xc, &ret, &exponent, &t, &diff, &one,
		    NULL);
    return err;
}

int mpf_ln_agm(mp_float * a, mp_float * b)
{
    int err, fract;
    long eps, oldeps, p, ilog2, sublog2, m;
    mp_float x, pi, minv, one, fourx;

    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;
    if ((err =
	 mpf_init_multi(eps, &one, &fourx, &x, &minv, &pi, NULL)) != MP_OKAY) {
	return err;
    }
    err = MP_OKAY;
    if ((err = mpf_copy(a, &x)) != MP_OKAY) {
	goto _ERR;
    }
    // if x < 1 compute log(1/x) = -log(x)
    fract = 0;
    if (mpf_isfraction(a)) {
	if ((err = mpf_inv(&x, &x)) != MP_OKAY) {
	    goto _ERR;
	}
	fract = 1;
    }
    // it said: x^m > 2^(precision/2) but it
    // is not that simple
    p = x.radix >> 1;
    ilog2 = (x.radix - abs(x.exp)) - 1;
    sublog2 = 0;
    if (ilog2 == 0) {
	sublog2 = 1;
	x.exp += 1;
	// uhm...well...
	ilog2 = (x.radix - abs(x.exp)) - 1;
    }
    m = p / ilog2;

    // a lower limit is simpler
    if (m < 10) {
	m = 10;
    }

    if ((err = mpf_const_d(&minv, m)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_inv(&minv, &minv)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_pow_d(&x, m, &x)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&fourx, 4)) != MP_OKAY) {
	goto _ERR;
    }
    mpf_div(&fourx, &x, &fourx);

    if ((err = mpf_agm(&one, &fourx, &x)) != MP_OKAY) {
	goto _ERR;
    }

    // (pi/2)/(AGM(1,4/x))
    if ((err = mpf_const_pi2(&pi)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_div(&pi, &x, &x)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_mul(&x, &minv, &x)) != MP_OKAY) {
	goto _ERR;
    }
    if (sublog2 == 1) {
	if ((err = mpf_const_ln_d(&one, 2)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_sub(&x, &one, &x)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_normalize_to(&x, a->radix)) != MP_OKAY) {
	goto _ERR;
    }
    if (fract == 1) {
	x.mantissa.sign = MP_NEG;
    }
    if ((err = mpf_copy(&x, b)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mpf_clear_multi(&x, &pi, &minv, &one, &fourx, NULL);
    return err;
}

