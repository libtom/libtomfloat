#include <tomfloat.h>
int mpf_nthroot(mp_float * a, long n, mp_float * b)
{
    int err, sign;
    mp_float one, t2, t3, invb, x0;
    mp_float ret;
    mp_float frac;
    long oldeps, nm1, loops, expnt;
    double d,rest;

    // TODO: IEEE-754 has a quite complicated system here, adapt
    if (mpf_iszero(a)) {
	if (n > 0) {
	    if ((err = mpf_const_0(b)) != MP_OKAY) {
		return err;
	    }
	    return MP_OKAY;
	} else {
	    if ((err = mpf_const_nan(b)) != MP_OKAY) {
		return err;
	    }
	    return MP_OKAY;
	}
    }
    if (mpf_isnan(a)) {
	if ((err = mpf_const_nan(b)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }
    if (mpf_isinf(a)) {
	if ((err = mpf_const_inf(b, MP_ZPOS)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }
    if (n < 0) {
	// 1/(x^(1/-n))
	n = -n;
	if ((err = mpf_init_copy(&ret, a)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_inv(&ret, &ret)) != MP_OKAY) {
            mpf_clear(&ret);
	    return err;
	}
	if ((err = mpf_nthroot(&ret, n, b)) != MP_OKAY) {
            mpf_clear(&ret);
	    return err;
	}
	mpf_clear(&ret);
	return MP_OKAY;
    }
    if (n == 1) {
	// a^(1/1)
	if ((err = mpf_copy(a, b)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }
    /* input must be positive if b is even */
    if ((n & 1) == 0 && a->mantissa.sign == MP_NEG) {
	return MP_VAL;
    }
    if (n == 2) {
	return mpf_sqrt(a, b);
    }
    if ((err = mpf_init(&one, a->radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_const_d(&one, 1)) != MP_OKAY) {
        mpf_clear(&one);
	return err;
    }

    if (mpf_cmp(a, &one) == MP_EQ) {
	if (n == 0) {
	    if ((err = mpf_const_nan(b)) != MP_OKAY) {
                mpf_clear(&one);
		return err;
	    }
	    mpf_clear(&one);
	    return MP_VAL;
	} else {
	    if ((err = mpf_const_0(b)) != MP_OKAY) {
                mpf_clear(&one);
		return err;
	    }
	    mpf_clear(&one);
	    return MP_OKAY;
	}
    }

    mpf_clear(&one);

    oldeps = a->radix;
    if ((err = mpf_init_copy(a, &ret)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_normalize_to(&ret, oldeps + MP_DIGIT_BIT)) != MP_OKAY) {
        mpf_clear(&ret);
	return err;
    }

    sign = a->mantissa.sign;
    ret.mantissa.sign = MP_ZPOS;

    /*
     * Compute initial value (assuming n is positive)
     * 
     * x^n = exp(log(x) * n)
     * x^n = 2^(log_2(x) * n)
     * 
     * (f * 2^e)^n = 2^(log_2(f * 2^e) * n)
     *             = 2^( (log_2(f) + e)   * n)
     *             = 2^( log_2(f) * n + e * n )
     * 
     * and with n the fraction 1/n
     * 
     * (f * 2^e)^(1/n)    = 2^( log_2(f) / n + e / n )
     * 
     * handling both factors individually:
     * 
     * (f * 2^e)^(1/n)  = f^(1/n) * (2^e)^(1/n)
     * 
     * f^(1/n) = 2^(log_2(f)/n)
     * 
     * (2^e)^(1/n) = 2^( log_2(2^e)/n )
     *             = 2^(e / n)
     * We need an integer for the exponent so we use the integer
     * part for the big exponent and add the small fractional part f
     * the big fractional part F by doing F*2^f
     */
    if (mpf_isdouble(a)) {
	if ((err = mpf_init(&t2, ret.radix)) != MP_OKAY) {
            mpf_clear(&ret);
	    return err;
	}
	if ((err = mpf_init(&frac, ret.radix)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,NULL);
	    return err;
	}
	if ((err = mpf_set_double(&ret, &d)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,&frac,NULL);
	    return err;
	}
        // TODO: this can underflow for very large n and very small d.
        //       Unlikely, but not impossible even for n = 2^32 and d = 1e-300
        //       More likely if the long type is larger than 32 bit
	d = pow(d, 1.0 / n);
	if ((err = mpf_get_double(d, &t2)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,&frac,NULL);
	    return err;
	}

    } else {
	if ((err = mpf_init(&t2, ret.radix)) != MP_OKAY) {
            mpf_clear(&ret);
	    return err;
	}
	if ((err = mpf_init(&frac, ret.radix)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,NULL);
	    return err;
	}
	if ((err = mpf_frexp(&ret, &frac, &expnt)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,&frac,NULL);
	    return err;
	}
	
	if ((err = mpf_set_double(&frac, &d)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,&frac,NULL);
	    return err;
	}
	// (f) ^(1/n) 
	d = pow(d, 1.0 / n);
        // integer part floor(e/n)
	expnt = expnt / n;
        // add the fractional part {e/n}
        rest = ((double)(expnt)) / ((double)(n)) ;
        rest = pow(2,rest);
        d *= rest;
	if ((err = mpf_get_double(d, &frac)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,&frac,NULL);
	    return err;
	}
	// f^(1/n) * 2^( log_2(2^e)/n )
	// f^(1/n) * 2^( e/n )
	if ((err = mpf_ldexp(&frac, expnt, &t2)) != MP_OKAY) {
            mpf_clear_multi(&ret,&t2,&frac,NULL);
	    return err;
	}
        if (d >= 1){
            if( (err = mpf_normalize(&t2) ) != MP_OKAY){
                mpf_clear_multi(&ret,&t2,&frac,NULL);
	        return err;
            }
        }
    }
    // we should have a good enough initial value in t2 now.
    if ((err = mpf_init(&invb, ret.radix)) != MP_OKAY) {
        mpf_clear_multi(&ret,&t2,&frac,NULL);
	return err;
    }
    if ((err = mpf_const_d(&invb, n)) != MP_OKAY) {
        mpf_clear_multi(&ret,&t2,&frac,&invb,NULL);
	return err;
    }
    if ((err = mpf_inv(&invb, &invb)) != MP_OKAY) {
        mpf_clear_multi(&ret,&t2,&frac,&invb,NULL);
	return err;
    }
    nm1 = n - 1;
    // safe guard for the worst case
    loops = ret.radix;
    if ((err = mpf_init_multi(ret.radix, &t3, &x0, NULL)) != MP_OKAY) {
        mpf_clear_multi(&ret,&t2,&frac,&invb,NULL);
	return err;
    }
    do {
	// x0 = t2
	if ((err = mpf_copy(&t2, &x0)) != MP_OKAY) {
	    goto _ERR;
	}
	// t3 = t2^nm1
	if ((err = mpf_pow_d(&t2, nm1, &t3)) != MP_OKAY) {
	    goto _ERR;
	}
	// t3 = a / t3
	if ((err = mpf_div(&ret, &t3, &t3)) != MP_OKAY) {
	    goto _ERR;
	}
	// t3 = t3 - t2
	if ((err = mpf_sub(&t3, &t2, &t3)) != MP_OKAY) {
	    goto _ERR;
	}
	// t3 = t3 * invb
	if ((err = mpf_mul(&t3, &invb, &t3)) != MP_OKAY) {
	    goto _ERR;
	}
	// t2 = t2 + t3
	if ((err = mpf_add(&t2, &t3, &t2)) != MP_OKAY) {
	    goto _ERR;
	}
	if (mpf_iszero(&t2)) {
	    break;
	}
	if (loops-- == 0) {
	    fprintf(stderr, "mp_float nthroot did not converge in some %ld rounds\n",
		    (ret.radix));
	    err = MP_VAL;
	    goto _ERR;
	}
    } while (mpf_cmp(&t2, &x0) != MP_EQ);
#ifdef DEBUG
    fprintf(stderr, "nth-root loops =  %ld\n", (ret.radix + 1) - loops);
#endif
    /* set the sign of the result */
    t2.mantissa.sign = sign;
    if ((err = mpf_normalize_to(&t2, a->radix)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mpf_copy(&t2, b)) != MP_OKAY) {
	goto _ERR;
    }
    err = MP_OKAY;

  _ERR:        
    mpf_clear_multi(&t2, &t3, &x0, &invb, &ret, &frac,NULL);
    return err;
}

