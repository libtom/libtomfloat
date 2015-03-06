#include <tomfloat.h>

int mpf_kernel_atan(mp_float * a, mp_float * b, int hyper)
{
    mp_float sum1, sum2, EPS, t2, t4, tt, tmp;
    long n, eps, oldeps;
    int err;

    // checks & balances
    // check for -1 < this < 1 must get done by caller
    // atan(h)(0) = 0
    if (mpf_iszero(a)) {
	mpf_copy(a, b);
	return MP_OKAY;
    }

    err = MP_OKAY;
    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;

    if ((err =
	 mpf_init_multi(eps, &sum1, &sum2, &EPS, &t2, &t4, &tt, &tmp,
			NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_copy(a, &sum1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&sum1, eps)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_sqr(&sum1, &t2)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sqr(&t2, &t4)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_div_d(&sum1, 3, &sum2)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_mul(&sum1, &t4, &tmp)) != MP_OKAY) {
	goto _ERR;
    }
    // we skipped the first terms, hence...
    n = 5;
    //...a while-loop instead of my usual do-while-loops because it could
    // already stop here if all the stars are aligned right.
    while (mpf_cmp(&tmp, &EPS) != MP_LT) {
	// Doing it this way saves some branching for the hyperbolic case and
	// because of the rule "legibility wins if it is cheaply to have" it is
	// done here.
	// It can be done for the cosine, too.
	if ((err = mpf_div_d(&tmp, n, &tt)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_add(&sum1, &tt, &sum1)) != MP_OKAY) {
	    goto _ERR;
	}
	n += 2;
	if ((err = mpf_div_d(&tmp, n, &tt)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_add(&sum2, &tt, &sum2)) != MP_OKAY) {
	    goto _ERR;
	}
	n += 2;
	if ((err = mpf_mul(&tmp, &t4, &tmp)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_mul(&sum2, &t2, &sum2)) != MP_OKAY) {
	goto _ERR;
    }
    if (hyper == 1) {
	if ((err = mpf_add(&sum1, &sum2, &sum1)) != MP_OKAY) {
	    goto _ERR;
	}
    } else {
	if ((err = mpf_sub(&sum1, &sum2, &sum1)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_normalize_to(&sum1, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&sum1, b)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mpf_clear_multi(&sum1, &sum2, &EPS, &t2, &t4, &tt, &tmp, NULL);
    return err;
}

