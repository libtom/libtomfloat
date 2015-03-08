#include <tomfloat.h>

#include <math.h>

static double lambertw_series(double x)
{   // probably too many but memory consumption is neglible here
    double b[20] = {
	-1.0,
	1.0,
	-0.3333333333333333333333333333333333333333,
	0.1527777777777777777777777777777777777778,
	-0.07962962962962962962962962962962962962963,
	0.04450231481481481481481481481481481481481,
	-0.02598471487360376249265138154027042915932,
	0.01563563253233392122281011169900058788948,
	-0.009616892024299431706839114246521653929061,
	0.006014543252956117860953251899753957367126,
	-0.003811298034891999226704302150118336675264,
	0.002440877991143982665896858528643753021570,
	-0.001576930344686784253923409539931411597316,
	0.001026263320507607154437548153390686105647,
	-0.0006720616311561362040020200434190752159122,
	0.0004424730618146209099302076085847372647923,
	-0.0002926772247296274448499290408144925243472,
	0.0001943872760545393178222591050496671973510,
	-0.0001295742668527488188824085734788824435592,
	0.00008665035805208127166045159046239029319060
    };
    double y, ret;
    int blen, i;
    y = sqrt(2*(1+ exp(1)*x));
    blen = 20;
    ret = b[blen - 1];
    for (i = blen - 2; i >= 0; i--) {
	ret = ret * y + b[i];
    }
    return ret;
}

static double lambertw_asymp(double x, int br)
{
    double ret, a, b;
    if (br == 0) {
        a = log(x);
        b = log(a);
    } else {
        a = log(-x);
        b = log(-a);
    }
    ret = a - b
          + b / a
          + (b * (-2.0 + b)) / (2.0*a*a)
          + (b * (6.0 - 9.0*b + 2.0*b*b)) / (6.0*a*a*a)
          + (b * (-12.0 + 36.0*b - 22.0*b*b + 3.0*b*b*b)) / (12.0*a*a*a*a)
          + (b * (60.0 - 300.0*b + 350.0*b*b - 125.0*b*b*b + 12.0*b*b*b*b))
                                                            / (60.0*a*a*a*a*a);
    return ret;
}
static double lambertw_Q01(double x)
{
    double ret, num, den;
    double a[] = {
	5.931375839364438,
	11.39220550532913,
	7.33888339911111,
	0.653449016991959
    };
    double b[] = {
	6.931373689597704,
	16.82349461388016,
	16.43072324143226,
	5.115235195211697
    };
    num =
	1.0 + a[0] * x + a[1] * x * x + a[2] * x * x * x + a[3] * x * x * x * x;
    den =
	1.0 + b[0] * x + b[1] * x * x + b[2] * x * x * x + b[3] * x * x * x * x;
    ret = x * num / den;
    return ret;
}

static double lambertw_Q02(double x)
{
    double ret, num, den;
    double a[] = {
	2.445053070726557,
	1.343664225958226,
	0.148440055397592,
	0.0008047501729130
    };
    double b[] = {
	3.444708986486002,
	3.292489857371952,
	0.916460018803122,
	0.0530686404483322
    };
    num =
	1.0 + a[0] * x + a[1] * x * x + a[2] * x * x * x + a[3] * x * x * x * x;
    den =
	1.0 + b[0] * x + b[1] * x * x + b[2] * x * x * x + b[3] * x * x * x * x;
    ret = x * num / den;
    return ret;
}
static double lambertw_Qm1(double x)
{
    double ret, num, den;
    double a[] = {
	-7.81417672390744,
	253.88810188892484,
	657.9493176902304
    };
    double b[] = {
	-60.43958713690808,
	99.9856708310761,
	682.6073999909428,
	962.1784396969866,
	1477.9341280760887
    };
    num = a[0] + a[1] * x + a[2] * x * x;
    den = 1 + b[0] * x + b[1] * x * x + b[2] * x * x * x + b[3] * x * x * x * x
	+ b[4] * x * x * x * x * x;
    ret = num / den;
    return ret;
}

static double lambertw_R(double x, int n)
{
    double ret;
    if (n == 0) {
	return log(-x) - log(-(log(-x)));
    } else {
	n--;
	ret = log(-x) - log(-lambertw_R(x, n));
	return ret;
    }
}


/**
   The Lambert-W function (aka. ProductLog) for real results.<br>

   Darko Veberic "Having Fun with Lambert W(x) Function", CoRR, 2010
   {@link http://arxiv.org/abs/1003.1628} with Halley's iteration instead
   of Fritsch's (for now).<br>
   This function is slow but usable.

   The branch must be either 0 or -1
*/
int mpf_lambertw(mp_float * a, mp_float * b, int branch)
{
    int err, k, i, sign, cmp;
    long oldeps, eps;
    double numx;
    mp_float w, wn, z, ret, branchpoint, negone, t1, t2, t3, t4, ONE, TWO, EPS;
    mp_float E, log2half, pihalf;

    k = branch;
    if (!mpf_isdouble(a)) {
	fprintf(stderr, "lambertw restricted to double magnitude for now\n");
	return MP_VAL;
    }
    oldeps = a->radix;
    eps = oldeps + MP_DIGIT_BIT;

    err = MP_OKAY;

    if ((err =
	 mpf_init_multi(eps, &w, &wn, &z, &ret, &branchpoint, &negone, &t1, &t2,
			&t3, &t4, &ONE, &TWO, &EPS, &E, &log2half, &pihalf,
			NULL)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_const_d(&ONE, 1)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_const_d(&TWO, 2)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_const_eps(&EPS)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_normalize_to(&EPS, oldeps)) != MP_OKAY) {
	return err;
    }

    if ((err = mpf_copy(a, &z)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_normalize_to(&z, eps)) != MP_OKAY) {
	return err;
    }

    if (branch == 0) {
	if ((err = mpf_const_e(&E)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_const_le2(&log2half)) != MP_OKAY) {
	    return err;
	}
	log2half.exp -= 1;
	log2half.mantissa.sign = MP_NEG;
	if ((err = mpf_const_pi(&pihalf)) != MP_OKAY) {
	    return err;
	}
	pihalf.exp -= 1;
	pihalf.mantissa.sign = MP_NEG;
	sign = z.mantissa.sign;
	z.mantissa.sign = MP_ZPOS;
	if (mpf_cmp(&z, &EPS) != MP_GT) {
	    if ((err = mpf_const_0(b)) != MP_OKAY) {
		return err;
	    }
	    goto _ERR;
	}
	if ((err = mpf_sub(&z, &E, &E)) != MP_OKAY) {
	    return err;
	}
	E.mantissa.sign = MP_ZPOS;
	if (mpf_cmp(&E, &EPS) != MP_GT) {
	    if ((err = mpf_copy(&ONE, b)) != MP_OKAY) {
		return err;
	    }
	    goto _ERR;
	}
	if ((err = mpf_sub(&z, &log2half, &log2half)) != MP_OKAY) {
	    return err;
	}
	log2half.mantissa.sign = MP_ZPOS;
	if (mpf_cmp(&log2half, &EPS) != MP_GT) {
	    if ((err = mpf_const_le2(b)) != MP_OKAY) {
		return err;
	    }
	    goto _ERR;
	}
	if ((err = mpf_sub(&z, &pihalf, &pihalf)) != MP_OKAY) {
	    return err;
	}
	pihalf.mantissa.sign = MP_ZPOS;
	if (mpf_cmp(&pihalf, &EPS) != MP_GT) {
	    // actually pi()/2 * I;
	    if ((err = mpf_const_nan(b)) != MP_OKAY) {
		return err;
	    }
	    err = MP_VAL;
	    goto _ERR;
	}
	z.mantissa.sign = sign;
    }

    if ((err = mpf_set_double(&z, &numx)) != MP_OKAY) {
	return err;
    }
    // check if z is >= the branchpoint at -1/e
    // next digit in -1/e would be a 7 (seven)
    if (numx < -0.3678794411) {
	if ((err = mpf_const_e(&E)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_inv(&E, &E)) != MP_OKAY) {
	    return err;
	}
	E.mantissa.sign = MP_NEG;
	cmp = mpf_cmp(&z, &E);
	if (cmp == MP_EQ) {
	    if ((err = mpf_const_d(b, 1)) != MP_OKAY) {
		return err;
	    }
	    b->mantissa.sign = MP_NEG;
	    goto _ERR;
	}
	if (cmp != MP_GT) {
	    err = MP_VAL;
	    goto _ERR;
	}
    }
    if (k == 0) {
	if (numx < -0.32358170806015724) {
	    numx = lambertw_series(numx);
	} else if (numx < 0.14546954290661823) {
	    numx = lambertw_Q01(numx);
	} else if (numx < 8.706658967856612) {
	    numx = lambertw_Q02(numx);
	} else {
	    numx = lambertw_asymp(numx, 0);
	}
    } else if (k == -1) {
	if (numx < -0.30298541769) {
	    numx = lambertw_series(numx);
	} else if (numx < -0.051012917658221676) {
	    numx = lambertw_Qm1(numx);
	} else {
	    numx = lambertw_R(numx, 10);
	}
    } else {
	// only branches with at least some real results supported
	err = MP_VAL;
	goto _ERR;
    }

    if ((err = mpf_get_double(numx, &w)) != MP_OKAY) {
	return err;
    }
    // Halley's iteration looses much accuracy near the singularities
    // at -1/e and 0
    if ((err = mpf_copy(&ONE, &negone)) != MP_OKAY) {
	return err;
    }
    negone.mantissa.sign = MP_NEG;
    for (i = 0; i < eps; i++) {
	// t1 = exp(w)
	if ((err = mpf_exp(&w, &t1)) != MP_OKAY) {
	    return err;
	}
	// t2 = w * t1 - z
	if ((err = mpf_mul(&w, &t1, &t3)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_sub(&t3, &z, &t2)) != MP_OKAY) {
	    return err;
	}
	if (mpf_cmp(&w, &negone) != MP_EQ) {
	    // t4 = w + 1
	    if ((err = mpf_add(&w, &ONE, &t4)) != MP_OKAY) {
		return err;
	    }
	    // t1 = t1 * (w + 1)
	    if ((err = mpf_mul(&t1, &t4, &t1)) != MP_OKAY) {
		return err;
	    }
	    // t3 = (w + 2.0)*t2/2
	    if ((err = mpf_add(&w, &TWO, &t3)) != MP_OKAY) {
		return err;
	    }
	    if ((err = mpf_mul(&t3, &t2, &t3)) != MP_OKAY) {
		return err;
	    }
	    t3.exp -= 1;
	    // t1*(w + 1) - 0.5*(w + 2.0)*t2/(w + 1)
	    // t1 = t1 - t3
	    if ((err = mpf_sub(&t1, &t3, &t1)) != MP_OKAY) {
		return err;
	    }
	    // t1 = t2/t1
	    if ((err = mpf_div(&t2, &t1, &t1)) != MP_OKAY) {
		return err;
	    }
	    sign = t1.mantissa.sign;
	    t1.mantissa.sign = MP_ZPOS;
	    if (mpf_cmp(&t1, &EPS) != MP_GT) {
		break;
	    }
	    t1.mantissa.sign = sign;
	}
	if ((err = mpf_sub(&w, &t1, &w)) != MP_OKAY) {
	    return err;
	}

    }

    if (i == eps) {
	fprintf(stderr, "Max. rounds in lambertw reached for z = ");
	mpf_dump(&z);
	fprintf(stderr, "Result after %d iterations = ", i);
	mpf_dump(&w);
    }
    if ((err = mpf_normalize_to(&w, oldeps)) != MP_OKAY) {
	return err;
    }
    mpf_exch(&w, b);

  _ERR:
    mpf_clear_multi(&w, &wn, &z, &ret, &branchpoint, &negone, &t1, &t2, &t3,
		    &t4, &ONE, &TWO, &EPS, &E, &log2half, &pihalf, NULL);
    return err;
}



