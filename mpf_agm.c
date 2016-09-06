#include <tomfloat.h>
// TODO: At least hint at the compiler to inline these two functions.
//       Or do it by hand?
// arithmetic mean am = (a + b)/2
static int mpf_alpha(mp_float * a, mp_float * b, mp_float * c)
{
    int err;
    if ((err = mpf_add(a, b, c)) != MP_OKAY) {
	return err;
    }
    c->exp -= 1;
    return MP_OKAY;
}

// geometric mean gm = sqrt(a * b)
static int mpf_beta(mp_float * a, mp_float * b, mp_float * c)
{
    int err;
    if ((err = mpf_mul(a, b, c)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_sqrt(c, c)) != MP_OKAY) {
	return err;
    }
    return MP_OKAY;
}

int mpf_agm(mp_float * a, mp_float * b, mp_float * c)
{
    mp_float A, aa, B, bb, diff, EPS;
    int err;
    long oldeps, eps;

    err = MP_OKAY;
    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;
    if ((err = mpf_init_multi(eps, &A, &aa, &B, &bb, &diff, &EPS, NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_alpha(a, b, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_beta(a, b, &B)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&EPS, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    do {
	if ((err = mpf_copy(&A, &aa)) != MP_OKAY) {err = -10;
	    goto _ERR;
	}
	if ((err = mpf_copy(&B, &bb)) != MP_OKAY) {err = -11;
	    goto _ERR;
	}
	if ((err = mpf_alpha(&aa, &bb, &A)) != MP_OKAY) {err = -12;
	    goto _ERR;
	}
	if ((err = mpf_beta(&aa, &bb, &B)) != MP_OKAY) {err = -13;
	    goto _ERR;
	}
	if ((err = mpf_sub(&A, &B, &diff)) != MP_OKAY) {err = -14;
	    goto _ERR;
	}
	if ((err = mpf_abs(&diff, &diff)) != MP_OKAY) {err = -15;
	    goto _ERR;
	}
    } while (mpf_cmp(&diff,&EPS) == MP_GT);

    if ((err = mpf_normalize_to(&A, oldeps)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_copy(&A, c)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mpf_clear_multi(&A, &aa, &B, &bb, &diff, &EPS, NULL);
    return err;
}

