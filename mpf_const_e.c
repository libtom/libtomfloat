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



#include <tomfloat.h>
//http://numbers.computation.free.fr/Constants/Algorithms/splitting.html
static int mp_e_binary_splitting(mp_int * a, mp_int * b, mp_int * P, mp_int * Q)
{
    int err;
    mp_int p1, q1, p2, q2, t1, one;
    if ((err = mp_init_multi(&p1, &q1, &p2, &q2, &t1, &one, NULL)) != MP_OKAY) {
	return err;
    }

    err = MP_OKAY;
    mp_set(&one, 1);
    if ((err = mp_sub(b, a, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if (mp_cmp(&t1, &one) == MP_EQ) {
	if ((err = mp_set_int(P, 1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_copy(b, Q)) != MP_OKAY) {
	    goto _ERR;
	}
	goto _ERR;
    }

    if ((err = mp_add(a, b, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_div_2d(&t1, 1, &t1, NULL)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mp_e_binary_splitting(a, &t1, &p1, &q1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_e_binary_splitting(&t1, b, &p2, &q2)) != MP_OKAY) {
	goto _ERR;
    }
    //P = q2*p1 + p2
    if ((err = mp_mul(&q2, &p1, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_add(&t1, &p2, P)) != MP_OKAY) {
	goto _ERR;
    }
    //Q =  q1*q2
    if ((err = mp_mul(&q1, &q2, Q)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mp_clear_multi(&p1, &q1, &p2, &q2, &t1, &one, NULL);
    return err;
}

static int mpf_exp1(mp_float * a)
{
    int err;
    long eps;
    mp_int p, q, zero, EPS;

    eps = a->radix + MP_DIGIT_BIT;
    err = MP_OKAY;

    if ((err = mp_init_multi(&p, &q, &zero, &EPS, NULL)) != MP_OKAY) {
	return err;
    }
    if ((err = mp_set_int(&EPS, eps)) != MP_OKAY) {
	goto _ERR;
    }
    mp_zero(&zero);

    if ((err = mp_e_binary_splitting(&zero, &EPS, &p, &q)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mp_add(&p, &q, &p)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul_2d(&p, eps, &p)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_div(&p, &q, &p, NULL)) != MP_OKAY) {
	goto _ERR;
    }

    if ((err = mpf_from_mp_int(&p, a)) != MP_OKAY) {
	goto _ERR;
    }
    a->exp -= eps;

  _ERR:
    mp_clear_multi(&p, &q, &zero, &EPS, NULL);
    return err;
}

static mp_float mpf_e;

static long mpf_e_precision;

int mpf_const_e(mp_float * a)
{
    int err;
    long eps;

    err = MP_OKAY;

    if (mpf_e_precision > 0 && a == NULL) {
	mpf_clear(&mpf_e);
	mpf_e_precision = 0;
	return err;
    }
    if (mpf_e_precision >= a->radix) {
	eps = a->radix;
	if ((err = mpf_copy(&mpf_e, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, eps);
    } else {
	if (mpf_e_precision == 0) {
	    if ((err = mpf_init(&mpf_e, a->radix)) != MP_OKAY) {
		return err;
	    }
	}
	if ((err = mpf_exp1(&mpf_e)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_copy(&mpf_e, a)) != MP_OKAY) {
	    return err;
	}
    }
    return MP_OKAY;
}

