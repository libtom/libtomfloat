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


/* compute b = e^a using e^x == \sum_{n=0}^{\infty} {1 \over n!}x^n */
int mpf_exp(mp_float * a, mp_float * b)
{
    long n, oldeps, eps, loops, decexpo;
    mp_float to, t, tx, ret, x0, one, two, nt, diff;
    int err, sign, m, i;

    sign = MP_ZPOS;
    err = MP_OKAY;

    // TODO: more checks & balances

    if (mpf_iszero(a)) {
	return mpf_const_d(a, 1);
    }

    oldeps = a->radix;

    decexpo = a->exp + a->radix;
    if (decexpo > 0) {
	decexpo = mpf_getdecimalexponent(decexpo) - 1;
	// TODO: evaluate the cutoffs and guard digits more exactly
	// (Mmh... I'm curious how the compiler optimizes it)
        eps = oldeps + 20;
	m = 1 << 4;
        if (decexpo > 1) {
	    eps = oldeps + 30;
	    m = 1 << 5;
	} else	if (decexpo > 2) {
	    eps = oldeps + 30;
	    m = 1 << 5;
	} else if (decexpo > 3) {
	    eps = oldeps + 60;
	    m = 1 << 6;
	} else if (decexpo > 4) {
	    eps = oldeps + 120;
	    m = 1 << 7;
	} else if (decexpo > 5) {
	    eps = oldeps + 240;
	    m = 1 << 8;
	} else if (decexpo > 6) {
	    eps = oldeps + 480;
	    m = 1 << 9;
	} else if (decexpo > 7) {
	    // exp(1e7) ~ 6.592 e4342944
	    // raise exponent overflow
	    err = MP_RANGE;
	    goto _ERR;
	}
    } else {
	m = 1 << 4;
	eps = oldeps + 20;
    }

    if ((err =
	 mpf_init_multi(eps, &to, &t, &tx, &ret, &x0, &one, &two, &nt, &diff,
			NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_copy(a, &nt)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&nt, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if (a->mantissa.sign == MP_NEG) {
	sign = MP_NEG;
	nt.mantissa.sign = MP_ZPOS;
    }

    if ((err = mpf_const_d(&ret, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&to, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&tx, 1)) != MP_OKAY) {
	goto _ERR;
    }
    n = 1;

    // NOTE: calculate a bit more precisely
    loops = eps + 10;
    // argument reduction by 1/2^m

    nt.exp -= m;
    do {
	// x0 = ret
	if ((err = mpf_copy(&ret, &x0)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = n++
	if ((err = mpf_const_d(&t, n++)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = 1/t
	if ((err = mpf_inv(&t, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	// to = t * to
	if ((err = mpf_mul(&t, &to, &to)) != MP_OKAY) {
	    goto _ERR;
	}
	// tx = tx * nt
	if ((err = mpf_mul(&tx, &nt, &tx)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = to * tx
	if ((err = mpf_mul(&to, &tx, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	// ret = t + ret
	if ((err = mpf_add(&t, &ret, &ret)) != MP_OKAY) {
	    goto _ERR;
	}
	if (loops-- == 0) {
	    fprintf(stderr, "exp did not converge in %ld rounds\n", eps + 10);
	    goto _ERR;
	}
    } while (mpf_cmp(&x0, &ret) != MP_EQ);
#ifdef DEBUG
    fprintf(stderr, "loops = %ld\n", (eps + 10) - loops);
#endif
    // we used the standard series to compute exp(z/2^m) + 1
    // ret = ret - 1;
    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sub(&ret, &one, &ret)) != MP_OKAY) {
	goto _ERR;
    }
    // reverse argument reduction
    for (i = 0; i < m; i++) {
	// to = ret^2
	if ((err = mpf_sqr(&ret, &to)) != MP_OKAY) {
	    goto _ERR;
	}
	// ret = ret * 2
	ret.exp += 1;
	// ret = ret + to
	if ((err = mpf_add(&ret, &to, &ret)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    // we have exp(z) - 1 now, add one unit
    if ((err = mpf_add(&ret, &one, &ret)) != MP_OKAY) {
	goto _ERR;
    }
    // exp(-z) = 1/exp(z)
    if (sign == MP_NEG) {
	if ((err = mpf_inv(&ret, &ret)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_normalize_to(&ret, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&ret, b)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mpf_clear_multi(&to, &t, &tx, &ret, &x0, &one, &two, &nt, &diff, NULL);
    return err;
}

