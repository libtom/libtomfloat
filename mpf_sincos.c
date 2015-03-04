#include <tomfloat.h>
/*
   More of an internal function but if you have use for it:

  int cosine = compute cosine if set to one
  int tan    = compute tangent if cosine _and_ tan are set to one
  int hyper  = compute hyperbolic functions instead

  This functions assumes that the argument is already reduced!

*/
static int mpf_print(mp_float * a)
{
    mp_fput(&(a->mantissa),10,stdout);
    printf(" * 2^%ld * 1.0\n", a->exp);
    return MP_OKAY;
}
int mpf_sincos(mp_float * a, mp_float * b, int cosine, int tan, int hyper)
{
    int err;
    mp_float sin, cos, tmp, atmp, A, EPS;
    long n, eps, oldeps;

    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;

    if ((err =
	 mpf_init_multi(eps, &sin, &cos, &tmp, &atmp, &A, &EPS,
			NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&A, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&A, &sin)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&cos, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&sin, &tmp)) != MP_OKAY) {
	goto _ERR;
    }

    n = 2;
    do {
	// sin/cos have alternating sums, the hyperbolic functions don't,
	// vice versa with the tangent
	if (hyper == 1) {
	    if ((err = mpf_mul(&tmp, &A, &tmp)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_div_d(&tmp, n, &tmp)) != MP_OKAY) {
		goto _ERR;
	    }
	} else {
	    if ((err = mpf_mul(&tmp, &A, &tmp)) != MP_OKAY) {
		goto _ERR;
	    }
	    if ((err = mpf_div_d(&tmp, -n, &tmp)) != MP_OKAY) {
		goto _ERR;
	    }
	}
	if (cosine == 1 || tan == 1) {
	    if ((err = mpf_add(&cos, &tmp, &cos)) != MP_OKAY) {
		goto _ERR;
	    }
	}
	n++;
	if ((err = mpf_mul(&tmp, &A, &tmp)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_div_d(&tmp, n, &tmp)) != MP_OKAY) {
	    goto _ERR;
	}

	if (cosine != 1 || tan == 1) {
	    if ((err = mpf_add(&sin, &tmp, &sin)) != MP_OKAY) {
		goto _ERR;
	    }
	}
	n++;
	if ((err = mpf_abs(&tmp, &atmp)) != MP_OKAY) {
	    goto _ERR;
	}
    } while (mpf_cmp(&atmp, &EPS) == MP_GT);

    if (tan == 1) {
	if ((err = mpf_div(&sin, &cos, b)) != MP_OKAY) {
	    goto _ERR;
	}
    } else if (cosine == 1 && tan != 1) {
	if ((err = mpf_copy(&cos, b)) != MP_OKAY) {
	    goto _ERR;
	}

    } else {
	if ((err = mpf_copy(&sin, b)) != MP_OKAY) {
	    goto _ERR;
	}

    }
    if ((err = mpf_normalize_to(b, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
  _ERR:
    mpf_clear_multi(&sin, &cos, &tmp, &atmp, &A, &EPS, NULL);
    return err;
}

