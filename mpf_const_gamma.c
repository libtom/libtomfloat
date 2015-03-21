#include <tomfloat.h>
#include <math.h>

// Xavier Gourdon & Pascal Sebah, The Euler constant: gamma
// http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf
// Jonathan Borwein & David Bailey, Mathematics by Experiment, A. K. Peters, 2003
// Algorithm 6 (3.8.65) (on page 173 in the authors edition), nearly verbatim
static int brent_macmillan_gamma(mp_float * a)
{
    int err;
    long eps, oldeps;
    unsigned long tmp;
    mp_int A, B, U, V, nsquare, k, t1;
    mp_float E, u, v;
    mp_digit dec;

    err = MP_OKAY;

    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;

    // the limit to end the loop, see paper of Borwein et al.
    dec = (mp_digit) ((mpf_getdecimalexponent(eps) * 208) / 100) + 3;
    if ((err =
	 mp_init_multi(&A, &B, &U, &V, &nsquare, &k, &t1, NULL)) != MP_OKAY) {
	return err;
    }
    // compute the n from the error O(e^(-4n))
    tmp = (unsigned long) (floor(log((eps / 4.0) * log(2)) / log(2))) + 1;
    // tmp = 31 for eps = 10^10 (It is quite a brave assumption that somebody
    //                           might use this lib for 10^10 bit long numbers)
    mp_set(&nsquare, 1);
    if ((err = mp_mul_2d(&nsquare, tmp, &nsquare)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_sqr(&nsquare, &nsquare)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_init(&E, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_le2(&E)) != MP_OKAY) {
	mpf_clear(&E);
	goto _ERR;
    }
    if ((err = mp_copy(&E.mantissa, &A)) != MP_OKAY) {
	mpf_clear(&E);
	goto _ERR;
    }
    mpf_clear(&E);

    if ((err = mp_mul_d(&A, tmp, &A)) != MP_OKAY) {
	goto _ERR;
    }
    A.sign = MP_NEG;
    if ((err = mp_copy(&A, &U)) != MP_OKAY) {
	goto _ERR;
    }
    mp_set(&B, 1);
    if ((err = mp_mul_2d(&B, eps, &B)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_copy(&B, &V)) != MP_OKAY) {
	goto _ERR;
    }
    mp_set(&k, 1);

    for (;;) {
	if ((err = mp_sqr(&k, &t1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_mul(&B, &nsquare, &B)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_div(&B, &t1, &B, NULL)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mp_mul(&A, &nsquare, &A)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_div(&A, &k, &A, NULL)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_add(&A, &B, &A)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_div(&A, &k, &A, NULL)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mp_add(&U, &A, &U)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_add(&V, &B, &V)) != MP_OKAY) {
	    goto _ERR;
	}
	if (k.dp[0] == dec) {
	    break;
	}
	if ((err = mp_add_d(&k, 1, &k)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_init_multi(eps, &u, &v, NULL)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_from_mp_int(&U, &u)) != MP_OKAY) {
	goto _SERR;
    }
    if ((err = mpf_from_mp_int(&V, &v)) != MP_OKAY) {
	goto _SERR;
    }
    if ((err = mpf_div(&u, &v, a)) != MP_OKAY) {
	goto _SERR;
    }
    if ((err = mpf_normalize_to(a, oldeps)) != MP_OKAY) {
	goto _SERR;
    }
  _SERR:
    mpf_clear_multi(&u, &v, NULL);
  _ERR:
    mp_clear_multi(&A, &B, &U, &V, &nsquare, &k, &t1, NULL);
    return err;
}


static mp_float mpf_egamma;

static long mpf_gamma_precision;

int mpf_const_gamma(mp_float * a)
{
    int err;
    long eps;

    err = MP_OKAY;

    if (mpf_gamma_precision > 0 && a == NULL) {
	mpf_clear(&mpf_egamma);
	mpf_gamma_precision = 0;
	return err;
    }
    if (mpf_gamma_precision >= a->radix) {
	eps = a->radix;
	if ((err = mpf_copy(&mpf_egamma, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, eps);
    } else {
	if (mpf_gamma_precision == 0) {
	    if ((err = mpf_init(&mpf_egamma, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = brent_macmillan_gamma(&mpf_egamma)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_copy(&mpf_egamma, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}

