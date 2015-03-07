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
#include <math.h>
int  mpf_invsqrt(mp_float *a, mp_float *b)
{
    int err, sign;
    long oldeps, eps, nloops, maxrounds;
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
    eps = oldeps + MP_DIGIT_BIT;
    err = MP_OKAY;

    if ((err =
	 mpf_init_multi(eps, &one, &x0, &xn, &A, &hn,
			NULL)) != MP_OKAY) {
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
    xn.exp /= 2; 
    if ((err = mpf_inv(&xn, &xn)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_const_d(&one, 1L)) != MP_OKAY) {
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
    // TODO: work with increasing precision, starting at radix = 50, the
    //       accuracy of the initial value and double each round
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
	    fprintf(stderr, "mpf_invsqrt did not converge in %ld rounds\n",
		    nloops);
	    return MP_RANGE;
	}

    } while (mpf_cmp(&x0, &xn) != MP_EQ);

    if ((err = mpf_normalize_to(&xn, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    mpf_exch(&xn, b);
  _ERR:
    mpf_clear_multi(&one, &x0, &xn, &A, &hn, &EPS, NULL);
    return err;
}

