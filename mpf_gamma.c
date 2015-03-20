#include <tomfloat.h>

#include <math.h>

static mp_float *spougecache;
static size_t spougecache_len;
static long spougecache_eps;

// flag for handling precision for the cache
static int reflection;

static int fill_spougecache(size_t A, int accuracy, long eps)
{
    int err, n;
    mp_float factrl, e, t1, t2, t3, pi;
    size_t start;

    if ((err =
	 mpf_init_multi(eps, &factrl, &e, &t1, &t2, &t3, NULL)) != MP_OKAY) {
	return err;
    }

    err = MP_OKAY;

    if (spougecache_len < A || spougecache_eps < accuracy) {
	//puts("filling spougecache");
	if ((err = mpf_const_d(&factrl, 1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_const_e(&e)) != MP_OKAY) {
	    goto _ERR;
	}

	if (spougecache_len != 0) {
	    spougecache = realloc(spougecache, (A + 1) * sizeof(mp_float));
	    if (spougecache == NULL) {
		return MP_MEM;
	    }
	    start = spougecache_len;
	} else {
	    spougecache = malloc((A + 1) * sizeof(mp_float));
	    if (spougecache == NULL) {
		return MP_MEM;
	    }
	    start = 1;
	    if ((err = mpf_init(&pi, eps)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_const_pi(&pi)) != MP_OKAY) {
		goto _ERR;
	    }
	    pi.exp += 1;
	    if ((err = mpf_init(&(spougecache[0]), eps)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_sqrt(&pi, &(spougecache[0]))) != MP_OKAY) {
		goto _ERR;
	    }
	    mpf_clear(&pi);
	}
	for (n = start; n < (int) A; n++) {
	    // to avoid the more expensive exp(log(a-n)*(n-0.5))
	    // TODO: check if exp() is fast enough now
	    if ((err = mpf_set_int(&t1, (int) (A - n))) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_pow_d(&t1, (int) (n - 1), &t2)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_sqrt(&t1, &t1)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_mul(&t2, &t1, &t2)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_pow_d(&e, (int) (A - n), &t3)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_mul(&t2, &t3, &t2)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_init(&(spougecache[n]), eps)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_div(&t2, &factrl, &(spougecache[n]))) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_set_int(&t1, -n)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_mul(&factrl, &t1, &factrl)) != MP_OKAY) {
		goto _ERR;
	    }
	}
	spougecache_eps = accuracy;
	spougecache_len = A;
    }
_ERR:
    mpf_clear_multi(&factrl, &e, &t1, &t2, &t3, NULL);
    return err;
}
// only (exp(log(gamma(a)))) for now
int mpf_lngamma(mp_float * a, mp_float * b)
{
    int oldeps, eps, err, n, accuracy;
    size_t A;
    mp_float z, ONE, factrl, e, pi, t1, t2, t3, sum;

    err = MP_OKAY;

    if (mpf_iszero(a)) {
	// raise singularity error
	return MP_VAL;
    }
    // Check is expensive and inexact.
    // if(mpf_isint(a){...}

    // very near zero
    if (a->exp + a->radix < -(a->radix)) {
	if ((err = mpf_ln(a, b)) != MP_OKAY) {
	    return err;
	}
	return err;
    }

    oldeps = a->radix;

    if (reflection == 1) {
	// angst-allowance seemed not necessary
	// TODO: check if that is really true
	eps = oldeps;
    } else {
	accuracy = mpf_getdecimalexponent(oldeps + MP_DIGIT_BIT);
	A = (size_t) ceil(1.252850440912568095810522965 * accuracy);
	// We need to compute the coefficients as exact as possible, so
	// increase the working precision acccording to the largest coefficient
	// TODO: probably too much
	eps = oldeps + ((A * 144) / 100);
	if ((err = fill_spougecache(A, accuracy, eps)) != MP_OKAY) {
	    return err;
	}
    }


    if ((err =
	 mpf_init_multi(eps, &z, &ONE, &factrl, &e, &t1, &t2, &t3, &sum,
			NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_copy(a, &z)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&z, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_set_int(&ONE, 1)) != MP_OKAY) {
	goto _ERR;
    }
    //   printf("splen = %u, A = %u, seps = %ld, acc = %d\n",
    //    spougecache_len , A , spougecache_eps , accuracy);

    if (a->mantissa.sign == MP_NEG) {
	//   eps = oldeps + 10;
	if ((err = mpf_copy(a, &z)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_normalize_to(&z, eps)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_init(&pi, eps)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_const_pi(&pi)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_normalize_to(&ONE, eps)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_sub(&ONE, &z, &t1)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mpf_mul(&pi, &z, &t2)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_sin(&t2, &t2)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_div(&pi, &t2, &t2)) != MP_OKAY) {
	    goto _ERR;
	}
	t2.mantissa.sign = MP_ZPOS;
	if ((err = mpf_ln(&t2, &t2)) != MP_OKAY) {
	    goto _ERR;
	}
	reflection = 1;
	if ((err = mpf_lngamma(&t1, &t1)) != MP_OKAY) {
	    goto _ERR;
	}
	reflection = 0;
	if ((err = mpf_sub(&t2, &t1, &t1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_normalize_to(&t1, oldeps)) != MP_OKAY) {
	    goto _ERR;
	}
	mpf_exch(&t1, b);
	goto _ERR;
    }

    if ((err = mpf_copy(&(spougecache[0]), &sum)) != MP_OKAY) {
	goto _ERR;
    }
    for (n = 1; n < (int) spougecache_len; n++) {
	if ((err = mpf_const_d(&t1, n)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_add(&z, &t1, &t1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_div(&(spougecache[n]), &t1, &t1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_add(&sum, &t1, &sum)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_const_d(&t1, (int) spougecache_len)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(&z, &t1, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    // 1/2
    ONE.exp -= 1;

    //ret = log(sum) + (-(t1)) + log(t1) * (z + 1/2);
    //    = log(sum) + t2 + log(t1) * (z + 1/2);
    //    = log(sum) + t2 + log(t1) * t3
    if ((err = mpf_neg(&t1, &t2)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(&z, &ONE, &t3)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_ln(&sum, &sum)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_ln(&t1, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_mul(&t1, &t3, &t3)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(&sum, &t2, &sum)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_add(&sum, &t3, &sum)) != MP_OKAY) {
	goto _ERR;
    }
    // ret = ret - log(z)
    if ((err = mpf_ln(&z, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_sub(&sum, &t1, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&t1, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    mpf_exch(&t1, b);
_ERR:
    mpf_clear_multi(&z, &ONE, &factrl, &e, &t1, &t2, &t3, &sum, NULL);
    return err;
}

void mpf_clear_spougecache()
{
    size_t i;
    for (i = 0; i < spougecache_len; i++) {
	mpf_clear(&(spougecache[i]));
    }
    free(spougecache);
}

void mpf_dump_spougecache()
{
    size_t i;
    for (i = 0; i < spougecache_len; i++) {
	printf("%u\n", i);
	mpf_dump(&(spougecache[i]));
    }
    printf("eps = %ld\n", spougecache_eps);
}
// TODO: adjust sign
int mpf_gamma(mp_float * a, mp_float * b)
{
    int err;
    long oldeps, eps;
    mp_float t;

    err = MP_OKAY;

    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;
    if ((err = mpf_init(&t, oldeps)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_copy(a, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&t, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_lngamma(&t, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_exp(&t, &t)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&t, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    mpf_exch(&t, b);
_ERR:
    mpf_clear(&t);
    return err;
}

