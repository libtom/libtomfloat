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

static int mpf_exp_newton(mp_float * a, mp_float * b);

/* compute b = e^a using e^x == \sum_{n=0}^{\infty} {1 \over n!}x^n */
int mpf_exp(mp_float * a, mp_float * b)
{
    long n, oldeps, eps, loops, decexpo;
    mp_float to, t, tx, ret, x0, one, two, nt, diff;
    int err, sign, m, i;

    sign = MP_ZPOS;
    err = MP_OKAY;

    // TODO: more checks & balances

    if (mpf_iszero(a)) {
	return mpf_const_d(a, 1);
    }



    oldeps = a->radix;

    if (oldeps > 2 * MP_DIGIT_BIT) {
puts("not more than one time?");
        return mpf_exp_newton(a,b);
    }

    decexpo = a->exp + a->radix;
    if (decexpo > 0) {
	decexpo = mpf_getdecimalexponent(decexpo) - 1;
	// TODO: evaluate the cutoffs and guard digits more exactly
	switch (decexpo) {
        case 0:
	case 1:
	    eps = oldeps + 20;
	    m = 1 << 4;
	    break;
	case 2:
	case 3:
	    eps = oldeps + 30;
	    m = 1 << 5;
	    break;
	case 4:
	    eps = oldeps + 60;
	    m = 1 << 6;
	    break;
	case 5:
	    eps = oldeps + 120;
	    m = 1 << 7;
	    break;
	case 6:
	    eps = oldeps + 240;
	    m = 1 << 8;
	    break;
	case 7:
	    eps = oldeps + 480;
	    m = 1 << 9;
	    break;
	case 8:
	    eps = oldeps + 960;
	    m = 1 << 10;
	    break;
/*
        This gives ~3.33561e434294481 for 1e9 instead of ~8.003e434294481
        There seems to be a still undetected flaw (over/underflow?), maybe
        even somewhere else.

	case 9:
	    eps = oldeps + 2000;
	    m = 1 << 11;
	    break;
*/
	default:
	    err = MP_RANGE;
	    break;
	}
    } else {
	m = 1 << 4;
	eps = oldeps + 20;
    }
    if ((err =
	 mpf_init_multi(eps, &to, &t, &tx, &ret, &x0, &one, &two, &nt,
			&diff, NULL)) != MP_OKAY) { 
	return err;
    }
    if ((err = mpf_copy(a, &nt)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_normalize_to(&nt, eps)) != MP_OKAY) {
	goto _ERR;
    }
    if (a->mantissa.sign == MP_NEG) {
	sign = MP_NEG;
	nt.mantissa.sign = MP_ZPOS;
    }

    if ((err = mpf_const_0(&ret)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&to, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&tx, 1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
	goto _ERR;
    }
    n = 1;

    // NOTE: calculate a bit more precisely
    loops = eps + 10;
    // argument reduction by 1/2^m

    nt.exp -= m;
    do {
	// x0 = ret
	if ((err = mpf_copy(&ret, &x0)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = n++
	if ((err = mpf_const_d(&t, n++)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = 1/t
	if ((err = mpf_inv(&t, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	// to = t * to
	if ((err = mpf_mul(&t, &to, &to)) != MP_OKAY) {
	    goto _ERR;
	}
	// tx = tx * nt
	if ((err = mpf_mul(&tx, &nt, &tx)) != MP_OKAY) {
	    goto _ERR;
	}
	// t = to * tx
	if ((err = mpf_mul(&to, &tx, &t)) != MP_OKAY) {
	    goto _ERR;
	}
	// ret = t + ret
	if ((err = mpf_add(&t, &ret, &ret)) != MP_OKAY) {
	    goto _ERR;
	}
	if (loops-- == 0) {
	    fprintf(stderr, "exp did not converge in %ld rounds\n", eps + 10);
	    goto _ERR;
	}
    } while (mpf_cmp(&x0, &ret) != MP_EQ);

#ifdef DEBUG
    fprintf(stderr, "loops = %ld\n", (eps + 10) - loops);
#endif
    // reverse argument reduction
    for (i = 0; i < m; i++) {
	// to = ret^2
	if ((err = mpf_sqr(&ret, &to)) != MP_OKAY) {
	    goto _ERR;
	}
	// ret = ret * 2
	ret.exp += 1;
	// ret = ret + to
	if ((err = mpf_add(&ret, &to, &ret)) != MP_OKAY) {
	    goto _ERR;
	}
    }

    // we have exp(z) - 1 now, add one unit
    if ((err = mpf_add(&ret, &one, &ret)) != MP_OKAY) {
	goto _ERR;
    }
    // exp(-z) = 1/exp(z)
    if (sign == MP_NEG) {
	if ((err = mpf_inv(&ret, &ret)) != MP_OKAY) {
	    goto _ERR;
	}
    }
    if ((err = mpf_normalize_to(&ret, oldeps)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&ret, b)) != MP_OKAY) {
	goto _ERR;
    }
  _ERR:
    mpf_clear_multi(&to, &t, &tx, &ret, &x0, &one, &two, &nt, &diff, NULL);
    return err;
}

// $x_{n+1} = x_n ( 1 + a - \log( x_n)  )$


static int mpf_exp_newton(mp_float * a, mp_float * b){
    long n, oldeps, eps, nloops, maxrounds, starteps, maxeps;
    mp_float t, x0, xn, one, A, EPS;
    int err, sign, m, i;

    sign = MP_ZPOS;
    err = MP_OKAY;

    oldeps = a->radix;
    eps = oldeps + 4 * MP_DIGIT_BIT;

    if ((err = mpf_init(&A, oldeps)) != MP_OKAY) {
        return err;
    }
    if ((err = mpf_copy(a, &A)) != MP_OKAY) {
	mpf_clear(&A);
        return err;
    }
    // get the seed with the series
    if ((err = mpf_normalize_to(&A, 2 * MP_DIGIT_BIT)) != MP_OKAY) {
	mpf_clear(&A);
        return err;
    }

    mpf_exp(&A,&A);

    maxrounds = A.radix;
    nloops = 0L;

    if ((err = mpf_init(&EPS, oldeps)) != MP_OKAY) {
	mpf_clear(&A);
        return err;
    }
    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	mpf_clear(&A);
	mpf_clear(&EPS);
        return err;
    }

    oldeps = a->radix;
    starteps = 2 * MP_DIGIT_BIT;
    maxeps = oldeps + MP_DIGIT_BIT;

    if ((err =
	 mpf_init_multi(starteps, &EPS, &t, &x0, &xn,&one, NULL)) != MP_OKAY) { 
	return err;
    }
    if ((err = mpf_normalize_to(&A, starteps)) != MP_OKAY) {
        goto _ERR;
    }
    if ((err = mpf_copy(&A, &xn)) != MP_OKAY) {
        goto _ERR;
    }
    if ((err = mpf_const_d(&one, 1L)) != MP_OKAY) {
        goto _ERR;
    }

    // $x_{n+1} = x_n ( 1 + a - \log( x_n)  )$
    do {
        if ((err = mpf_copy(&xn, &x0)) != MP_OKAY) {
            goto _ERR;
        }
        starteps = starteps * 2;
        if (starteps > maxeps) {
            // do one round with full precision
            starteps = maxeps;
        }
        if ((err =
             mpf_normalize_to_multi(starteps, &one, &x0, &xn, &A, &t,
                                    &EPS,NULL)) != MP_OKAY) {
            goto _ERR;
        }

        // t = log(x_n)
        if ((err = mpf_ln(&xn, &t)) != MP_OKAY) {
            goto _ERR;
        }

        if ((err = mpf_copy(a, &A)) != MP_OKAY) {
            goto _ERR;
        }
        if ((err = mpf_normalize_to(&A, starteps)) != MP_OKAY) {
            goto _ERR;
        }
        // t = a - t
        if ((err = mpf_sub(&A, &t, &t)) != MP_OKAY) {
            goto _ERR;
        }
        // t = 1 + t
        if ((err = mpf_add(&t, &one, &t)) != MP_OKAY) {
            goto _ERR;
        }
        // xn = xn * t
        if ((err = mpf_mul(&xn, &t, &xn)) != MP_OKAY) {
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

    if ((err = mpf_normalize_to(&xn, oldeps)) != MP_OKAY) {
       goto _ERR;    
    }

    mpf_exch(&xn, b);

  _ERR:
    mpf_clear_multi(&t, &x0, &xn, &one, NULL);
    return err;
}
