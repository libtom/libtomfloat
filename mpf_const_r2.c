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
static mp_float mpf_r2;

static long mpf_r2_precision;

int mpf_const_r2(mp_float * a)
{
    int err;
    long eps;

    err = MP_OKAY;

    if (mpf_r2_precision > 0 && a == NULL) {
	mpf_clear(&mpf_r2);
	mpf_r2_precision = 0;
	return err;
    }
    if (mpf_r2_precision >= a->radix) {
	eps = a->radix;
	if ((err = mpf_copy(&mpf_r2, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, eps);
    } else {
	if (mpf_r2_precision == 0) {
	    if ((err = mpf_init(&mpf_r2, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = mpf_const_d(&mpf_r2, 2)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_sqrt(&mpf_r2, &mpf_r2)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_copy(&mpf_r2, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}

// TODO: is slower but still lacks some optimizsations
// (275807/195025) * (1-1/76069501250)^(-1/2); nope, because 76069501250>2^31
// (47321/33461) * (1-1/2239277042)^(-1/2); 2239277042/2 = 1119638521 < 2^31
//                                          can be done directly with c-99
int mpf_const_r2_new(mp_float * a)
{
    int err;
    long eps;
    mp_float t1, t2, one;

    err = MP_OKAY;

    if (mpf_r2_precision > 0 && a == NULL) {
	mpf_clear(&mpf_r2);
	mpf_r2_precision = 0;
	return err;
    }
    if (mpf_r2_precision >= a->radix) {
	eps = a->radix;
	if ((err = mpf_copy(&mpf_r2, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, eps);
    } else {
	if ((err = mpf_init_multi(a->radix, &t1, &t2, &one, NULL)) != MP_OKAY) {
	    return err;
	}
	if (mpf_r2_precision == 0) {
	    if ((err = mpf_init(&mpf_r2, a->radix)) != MP_OKAY) {
		goto _ERR;
	    }
	}
	if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_const_d(&t1, 47321)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_const_d(&t2, 33461)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mpf_div(&t1, &t2, &t1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_const_d(&t2, 1119638521)) != MP_OKAY) {
	    goto _ERR;
	}
	t2.exp += 1;

	if ((err = mpf_inv(&t2, &t2)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mpf_sub(&one, &t2, &t2)) != MP_OKAY) {
	    goto _ERR;
	}


	if ((err = mpf_sqrt(&t2, &mpf_r2)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mpf_inv(&mpf_r2, &mpf_r2)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_mul(&t1, &mpf_r2, &mpf_r2)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mpf_copy(&mpf_r2, a)) != MP_OKAY) {
	    goto _ERR;
	}
      _ERR:
	mpf_clear_multi(&t1, &t2, &one, NULL);

    }
    return err;
}

