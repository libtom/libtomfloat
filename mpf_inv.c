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

#include <tommath.h>

/* compute 1/x */
int mpf_inv_newton(mp_float * a, mp_float * b)
{
    int err, sign, sign2;
    double d;
    long expnt, oldeps, eps, nloops, maxrounds;
    mp_float frac, one, x0, xn, A, hn, EPS;

    /* get sign of input */
    sign = a->mantissa.sign;

    /* force to positive */
    a->mantissa.sign = MP_ZPOS;

    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;

    err = MP_OKAY;

    if(mpf_iszero(a)){
      // raise DivisionByZero
      return MP_VAL;
    }


    if ((err = mpf_init_multi(eps, &frac, &one, &x0, &xn, &A, &hn,NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_frexp(a, &frac, &expnt)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_set_double(&frac, &d)) != MP_OKAY) {
	goto _ERR;
    }
    d = 1.0 / d;
    if ((err = mpf_get_double(d, &frac)) != MP_OKAY) {
	goto _ERR;
    }
    expnt = -expnt + 1;
    if ((err = mpf_ldexp(&frac, expnt, &frac)) != MP_OKAY) {
	goto _ERR;
    }
    // TODO calculate guard digits more exactly and do it loop-wise
    if ((err = mpf_normalize_to(&frac, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&frac, &xn)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&one, 1L)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&A, eps)) != MP_OKAY) {
	goto _ERR;
    }
    maxrounds = A.radix;
    nloops = 0L;
    if ((err = mpf_init(&EPS, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	goto _ERR;
    }
    do {
	if ((err = mpf_copy(&xn, &x0)) != MP_OKAY) {
	    goto _ERR;
	}
	// hn = 1 - (A * xn);
	if ((err = mpf_mul(&A, &xn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_sub(&one, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	sign2 = hn.mantissa.sign;
	hn.mantissa.sign = MP_ZPOS;
	// It makes more sense to compare after that limit is reached
	if (hn.exp <= EPS.exp) {
	    if (mpf_cmp(&hn, &EPS) != MP_GT) {
		break;
	    }
	}
	hn.mantissa.sign = sign2;
	// xn = xn + (xn * hn);
	if ((err = mpf_mul(&xn, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_add(&xn, &hn, &xn)) != MP_OKAY) {
	    goto _ERR;
	}
	nloops++;
	if (nloops >= maxrounds) {
	    // it might be a bug elsewhere, please report
	    fprintf(stderr, "mpf_inv did not converge in %ld rounds", nloops);
	    return MP_RANGE;
	}
      // Comparing is not necessary here but might save an iteration
      } while (mpf_cmp(&x0, &xn) != MP_EQ);
//      } while (1);
    if ((err = mpf_normalize_to(&xn, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    mpf_exch(&xn, b);
    /* now restore the signs */
    a->mantissa.sign = sign;
    b->mantissa.sign = sign;

_ERR:

    mpf_clear_multi(&frac, &one, &x0, &xn, &A, &hn, &EPS ,NULL);
    return err;
}

// not faster?
int mpf_inv(mp_float * a, mp_float * b)
{
    int err, sign, sign2;
    double d;
    long expnt, oldeps, eps, nloops, maxrounds, starteps, maxeps;
    mp_float frac, one, x0, xn, A, hn,hn2, EPS;

    /* get sign of input */
    sign = a->mantissa.sign;

    /* force to positive */
    a->mantissa.sign = MP_ZPOS;

    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;

    starteps = 2 * MP_DIGIT_BIT;
    maxeps = oldeps + MP_DIGIT_BIT;

    err = MP_OKAY;

    if(mpf_iszero(a)){
      // raise DivisionByZero
      return MP_VAL;
    }


    if ((err = mpf_init_multi(starteps, &frac, &one, &x0, &xn, &A, &hn,&hn2,NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_frexp(a, &frac, &expnt)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_set_double(&frac, &d)) != MP_OKAY) {
	goto _ERR;
    }
    d = 1.0 / d;
    if ((err = mpf_get_double(d, &frac)) != MP_OKAY) {
	goto _ERR;
    }
    expnt = -expnt + 1;
    if ((err = mpf_ldexp(&frac, expnt, &frac)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_copy(&frac, &xn)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&one, 1L)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&A, starteps)) != MP_OKAY) {
	goto _ERR;
    }
    maxrounds = A.radix;
    nloops = 0L;
    if ((err = mpf_init(&EPS, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	goto _ERR;
    }
    do {
	starteps = starteps * 2;
	if (starteps > maxeps) {
	    // do one round with full precision
	    // or...well...die from exhaustion
	    starteps = maxeps;
	}
	if ((err =
	     mpf_normalize_to_multi(starteps, &one, &x0, &xn, &A, &hn,&hn2,
				    &EPS, NULL)) != MP_OKAY) {
	    return err;
	}

	if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	    goto _ERR;
	}
        if ((err = mpf_normalize_to(&A, starteps)) != MP_OKAY) {
	    goto _ERR;
        }
	if ((err = mpf_copy(&xn, &x0)) != MP_OKAY) {
	    goto _ERR;
	}

	// hn = 1 - (A * xn);
	if ((err = mpf_mul(&A, &xn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_sub(&one, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	sign2 = hn.mantissa.sign;
	hn.mantissa.sign = MP_ZPOS;
	// It makes more sense to compare after that limit is reached
	if (hn.exp <= EPS.exp) {
	    if (mpf_cmp(&hn, &EPS) != MP_GT) {
		break;
	    }
	}
	hn.mantissa.sign = sign2;

	// xn = xn + xn( hn + hn^2)

	if ((err = mpf_sqr(&hn, &hn2)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_add(&hn, &hn2, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
        if ((err = mpf_mul(&xn, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_add(&xn, &hn, &xn)) != MP_OKAY) {
	    goto _ERR;
	}

	nloops++;
	if (nloops >= maxrounds) {
	    // it might be a bug elsewhere, please report
	    fprintf(stderr, "mpf_inv did not converge in %ld rounds", nloops);
	    return MP_RANGE;
	}
      // Comparing is not necessary here but might save an iteration
      } while (mpf_cmp(&x0, &xn) != MP_EQ);
//      } while (1);
    if ((err = mpf_normalize_to(&xn, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    mpf_exch(&xn, b);
    /* now restore the signs */
    a->mantissa.sign = sign;
    b->mantissa.sign = sign;

_ERR:

    mpf_clear_multi(&frac, &one, &x0, &xn, &A, &hn, &EPS ,NULL);
    return err;
}
