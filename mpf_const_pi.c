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

static mp_float mpf_pi;
// This is automatically initialized to 0 according to C99
// TODO: check the other standard versions
static long mpf_pi_precision;

/* Pi by the AGM (Brent-Salamin) */
int mpf_const_pi(mp_float * a)
{
    int err;
    mp_float aa, b, d, d2, s, t, p, two, twoinv;
    long eps, oldeps, extra, k, loops, r;

    err = MP_OKAY;

    // Sometimes memory is short, even today
    if (mpf_pi_precision > 0 && a == NULL) {
	mpf_clear(&mpf_pi);
	mpf_pi_precision = 0;
	return err;
    }

    oldeps = a->radix;

    // two percent plus  angst-allowance
    // TODO: compute correct value
    extra = (oldeps / 100) * 2 + MP_DIGIT_BIT;
    if (mpf_pi_precision >= oldeps) {
	if ((err = mpf_copy(&mpf_pi, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, oldeps);
    } else {
	eps = oldeps + extra;
	if (mpf_pi_precision <= 0) {
	    if ((err = mpf_init(&mpf_pi, eps)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = mpf_normalize_to(&mpf_pi, eps)) != MP_OKAY) {
	    return err;
	}
	if ((err =
	     mpf_init_multi(eps, &aa, &b, &d, &d2, &s, &t, &p, &two, &twoinv,
			    NULL)) != MP_OKAY) {
	    return err;
	}
	// this algorithm is quadratic, produces twice the amount of digits
	// every round, so the maximum number of rounds is log_2(radix)
	loops = 0;
	// if eps is negative we have a really large problem but somewhere else
	r = eps;
	while (r >>= 1) {
	    loops++;
	}
	// ceil(log_2(radix)) + angst-allowance
	loops += 3;
	/*
	 * Initialize:
	 * a = 1, b = 1/sqrt(2), t = 1/2, k = 1
	 */
	// aa = 1
	if ((err = mpf_const_d(&aa, 1)) != MP_OKAY) {
	    goto _ERR;
	}
	// two = 2
	if ((err = mpf_const_d(&two, 2)) != MP_OKAY) {
	    goto _ERR;
	}
	// twoinv = 1/two
	if ((err = mpf_inv(&two, &twoinv)) != MP_OKAY) {
	    goto _ERR;
	}
	// b = sqrt(two)
	if ((err = mpf_sqrt(&two, &b)) != MP_OKAY) {
	    goto _ERR;
	}
	// b = 1/b
	if ((err = mpf_inv(&b, &b)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = twoinv
	if ((err = mpf_copy(&twoinv, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	k = 1;

	do {
	    /* s = 1/2 * (a + b) */
	    // s = aa + b
	    if ((err = mpf_add(&aa, &b, &s)) != MP_OKAY) {
		goto _ERR;
	    }
	    // s *= twoinv
	    if ((err = mpf_mul(&s, &twoinv, &s)) != MP_OKAY) {
		goto _ERR;
	    }
	    /* d = aa - s */
	    if ((err = mpf_sub(&aa, &s, &d)) != MP_OKAY) {
		goto _ERR;
	    }
	    // d = d^2
	    if ((err = mpf_sqr(&d, &d)) != MP_OKAY) {
		goto _ERR;
	    }
	    // aa = s
	    if ((err = mpf_copy(&s, &aa)) != MP_OKAY) {
		goto _ERR;
	    }
	    // s = s^2
	    if ((err = mpf_sqr(&s, &s)) != MP_OKAY) {
		goto _ERR;
	    }
	    /* t = t- 2^k * d */
	    // d2 = d * 2^k
	    if ((err = mpf_copy(&d, &d2)) != MP_OKAY) {
		goto _ERR;
	    }
	    d2.exp += k;
	    // t = t - d2
	    if ((err = mpf_sub(&t, &d2, &t)) != MP_OKAY) {
		goto _ERR;
	    }
	    /* b = sqrt(s - d) */
	    // b = s - d
	    if ((err = mpf_sub(&s, &d, &b)) != MP_OKAY) {
		goto _ERR;
	    }
	    // b = sqrt(b);
	    if ((err = mpf_sqrt(&b, &b)) != MP_OKAY) {
		goto _ERR;
	    }
	    k++;
	    if (k == loops) {
		err = MP_VAL;
		fprintf(stderr, "Pi by AGM did not converge in %ld rounds\n",
			k);
printf("aa ");mp_put(&(aa.mantissa),16);puts("");
printf("bb ");mp_put(&(b.mantissa),16);puts("");
		goto _ERR;
	    }
	} while (mpf_cmp(&aa, &b) != MP_EQ);
	/*
	 *       (a + b)^2
	 * pi = -----------
	 *          2t
	 */
	// t = t * 2
	t.exp += 1;
	// t = 1/t
	if ((err = mpf_inv(&t, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	// p = aa + b but aa = b, so p = aa * 2
	aa.exp += 1;
	// p = p^2
	if ((err = mpf_sqr(&aa, &p)) != MP_OKAY) {
	    goto _ERR;
	}
	// p = p * t
	if ((err = mpf_mul(&p, &t, &mpf_pi)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mpf_copy(&mpf_pi, a)) != MP_OKAY) {
	    goto _ERR;
	}
	mpf_pi_precision = eps;
	return mpf_normalize_to(a, oldeps);
    }
  _ERR:
    // Make it clear that it was a failure and that the value of mpf_pi cannot
    // be trusted anymore by setting mpf_pi_precision to a value that
    // cannot happen naturally.
    mpf_pi_precision = -1;
    mpf_clear_multi(&aa, &b, &d, &d2, &s, &t, &p, &two, &twoinv, NULL);
    return err;
}

