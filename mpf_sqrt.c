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
// as elegant as this is, but it is to slow for larger precisions and the AGM
int mpf_sqrt_old(mp_float * a, mp_float * b)
{
    int err;
    long oldeps, eps;
    mp_float A, B;

    oldeps = a->radix;
    eps = 2 * oldeps + MP_DIGIT_BIT;

    err = MP_OKAY;

    if ((err = mpf_init_multi(eps, &A, &B, NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&A, eps)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_invsqrt(&A, &B)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_mul(&A, &B, &B)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_normalize_to(&B, oldeps)) != MP_OKAY) {
	goto _ERR;
    }

    mpf_exch(&B, b);
  _ERR:
    mpf_clear_multi(&A, &B, NULL);
    return err;
}


int mpf_sqrt_newton(mp_float * a, mp_float * b)
{
    int err, sign;
    long oldeps, eps, nloops, maxrounds, expnt, rest;
    double d;
    mp_float one, x0, xn, A, hn, EPS, frac;

    if (mpf_isnan(a) || mpf_isinf(a)) {
	return mpf_copy(a, b);
    }

    if (a->mantissa.sign == MP_NEG) {
	return MP_VAL;
    }
    if (mpf_iszero(a)) {
	return mpf_const_0(b);
    }
    oldeps = a->radix;
    eps = oldeps + 4 * MP_DIGIT_BIT ;
    err = MP_OKAY;

    if ((err = mpf_init_multi(eps, &one, &x0, &xn, &A, &hn, &frac, NULL)) != MP_OKAY) {
	return err;
    }


    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&A, eps)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_copy(&A, &xn)) != MP_OKAY) {
	goto _ERR;
    }

    // compute seed
    if ((err = mpf_frexp(&xn, &frac, &expnt)) != MP_OKAY) {
        goto _ERR;
    }
    if ((err = mpf_set_double(&frac, &d)) != MP_OKAY) {
        goto _ERR;
    }
    // (f) ^(1/n) 
    d = sqrt(d);
    // integer part floor(e/2)
    expnt = expnt / 2;
    // add the fractional part {e/2}
    rest = ((double)(expnt)) / (2.0) - expnt;
    rest = pow(2,rest);
    d *= rest;
    if ((err = mpf_get_double(d, &frac)) != MP_OKAY) {
        goto _ERR;
    }
    // f^(1/2) * 2^( e/2 )
    if ((err = mpf_ldexp(&frac, expnt, &xn)) != MP_OKAY) {
        goto _ERR;
    }
    if (d >= 1.0){
        if( (err = mpf_normalize(&xn) ) != MP_OKAY){
            goto _ERR;
        }
    }

    if ((err = mpf_inv(&xn, &xn)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_const_d(&one, 1L)) != MP_OKAY) {
	goto _ERR;
    }

    maxrounds = (long)(log(eps)/log(2)) + 5;
    nloops = 0L;
    if ((err = mpf_init(&EPS, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	goto _ERR;
    }
    // TODO: work with increasing precision, starting at radix = 50 (or so)
    do {
	if ((err = mpf_copy(&xn, &x0)) != MP_OKAY) {
	    goto _ERR;
	}
	// hn = 1 - (A * xn^2);
	if ((err = mpf_sqr(&xn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_mul(&A, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_sub(&one, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	sign = hn.mantissa.sign;
	hn.mantissa.sign = MP_ZPOS;
	// It makes more sense to compare after that limit is reached
	if (hn.exp <= EPS.exp) {
	    if (mpf_cmp(&hn, &EPS) != MP_GT) {
		break;
	    }
	}
	hn.mantissa.sign = sign;
	// x(n+1) = xn + xn/2 * hn
	xn.exp -= 1;
	if ((err = mpf_mul(&xn, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	xn.exp += 1;
	if ((err = mpf_add(&xn, &hn, &xn)) != MP_OKAY) {
	    goto _ERR;
	}
	nloops++;
	if (nloops >= maxrounds) {
	    // it might be a bug elsewhere, please report
	    fprintf(stderr, "mpf_sqrt did not converge in %ld rounds\n",
		    nloops);
	    return MP_RANGE;
	}

    } while (mpf_cmp(&x0, &xn) != MP_EQ);
    if ((err = mpf_mul(&A, &xn, &xn)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&xn, a->radix)) != MP_OKAY) {
	goto _ERR;
    }
    mpf_exch(&xn, b);
  _ERR:
    mpf_clear_multi(&one, &x0, &xn, &A, &hn, &EPS, NULL);
    return err;
}

//TODO: working with increasing precision does not gain very much but something
//      and this something gets significant with high precision.
//      On the other side: the overhead of doing it is not insignificant and
//      influences the operation time on low precisions to the worse.
//      Example for sqrt(2) 1,000 times on my good old 1GHz Duron and 112 bits
//      precision.
//      old:  ~200ms
//      new:  ~220ms
//      It already wins at 333 bits precision (ca 100 dec. digits)
//      old:  ~320ms
//      new:  ~270ms
//      Significant at 3333 bits precision (ca 1000 dec. digits)
//      old:  ~7.0s
//      new:  ~3.5s
//      Alternatives: check size and branch or raise starting precision
//      NOTE: halfs time of computing log by the AGM at 100 dec. digits already
int mpf_sqrt(mp_float * a, mp_float * b)
{
    int err, sign;
    long oldeps, starteps, maxeps, nloops, maxrounds;
    mp_float one, x0, xn, A, hn, EPS;

    if (mpf_isnan(a) || mpf_isinf(a)) {
	return mpf_copy(a, b);
    }

    if (a->mantissa.sign == MP_NEG) {
	return MP_VAL;
    }
    if (mpf_iszero(a)) {
	return mpf_const_0(b);
    }
    oldeps = a->radix;
    starteps = 64;//2 * MP_DIGIT_BIT;
    maxeps = oldeps + MP_DIGIT_BIT;

    err = MP_OKAY;

    if ((err =
	 mpf_init_multi(starteps, &one, &x0, &xn, &A, &hn, NULL)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&A, starteps)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_copy(&A, &xn)) != MP_OKAY) {
	goto _ERR;
    }

    xn.exp -= (xn.exp + xn.radix) / 2  - 1;

    if ((err = mpf_inv(&xn, &xn)) != MP_OKAY) {	goto _ERR;    }

    if ((err = mpf_const_d(&one, 1L)) != MP_OKAY) {
	goto _ERR;
    }

    maxrounds = (long)(log(maxeps)/log(2)) + 5;
    nloops = 0L;
    // hammer a nail into the floor, it shall be the calibration point
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
	     mpf_normalize_to_multi(starteps, &one, &x0, &xn, &A, &hn,
				    &EPS, NULL)) != MP_OKAY) {
	    return err;
	}
	// TODO: check if "one" always holds one (it must, but one never knows)
	if ((err = mpf_copy(&xn, &x0)) != MP_OKAY) {
	    goto _ERR;
	}
	// hn = 1 - (A * xn^2);
	if ((err = mpf_sqr(&xn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	    goto _ERR;
	}
        if ((err = mpf_normalize_to(&A, starteps)) != MP_OKAY) {
	    goto _ERR;
        }
	if ((err = mpf_mul(&A, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mpf_sub(&one, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	sign = hn.mantissa.sign;
	hn.mantissa.sign = MP_ZPOS;
//fprintf(stderr,"SQRT 5 %ld %ld\n",hn.exp, EPS.exp);
        if (mpf_iszero(&hn))
            break;
	// It makes more sense to compare after that limit is reached
	if (hn.exp <= EPS.exp) {
	    if (mpf_cmp(&hn, &EPS) != MP_GT) {
		break;
	    }
	}
	hn.mantissa.sign = sign;
	// x(n+1) = xn + xn/2 * hn
	xn.exp -= 1;
	if ((err = mpf_mul(&xn, &hn, &hn)) != MP_OKAY) {
	    goto _ERR;
	}
	xn.exp += 1;
	if ((err = mpf_add(&xn, &hn, &xn)) != MP_OKAY) {
	    goto _ERR;
	}
	nloops++;
	if (nloops >= maxrounds) {
	    // it might be a bug elsewhere, please report
	    fprintf(stderr, "mpf_sqrt did not converge in %ld rounds\n",
		    nloops);
	    return MP_RANGE;
	}

    } while (mpf_cmp(&x0, &xn) != MP_EQ);
    // TODO: renormalize here if maxeps is very large
    //if ((err = mpf_normalize_to(&xn, oldeps + MP_DIGIT_BIT)) != MP_OKAY) {	goto _ERR;    }
    //if ((err = mpf_normalize_to(&A, oldeps + MP_DIGIT_BIT)) != MP_OKAY) {	goto _ERR;    }

    if ((err = mpf_mul(&A, &xn, &xn)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&xn, oldeps)) != MP_OKAY) {	goto _ERR;    }

    mpf_exch(&xn, b);
  _ERR:
    mpf_clear_multi(&one, &x0, &xn, &A, &hn, &EPS, NULL);
    return err;
}

