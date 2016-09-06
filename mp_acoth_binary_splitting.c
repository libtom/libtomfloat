#include <tomfloat.h>
//http://numbers.computation.free.fr/Constants/Algorithms/splitting.html
int mp_acoth_binary_splitting(mp_int * q, mp_int * a, mp_int * b, mp_int * P,
			      mp_int * Q, mp_int * R)
{
    int err;
    mp_int p1, q1, r1, p2, q2, r2, t1, t2, one;
    if ((err =
	 mp_init_multi(&p1, &q1, &r1, &p2, &q2, &r2, &t1, &t2, &one,
		       NULL)) != MP_OKAY) {
	return err;
    }

    err = MP_OKAY;
    mp_set(&one, 1);
    if ((err = mp_sub(b, a, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if (mp_cmp(&t1, &one) == MP_EQ) {
	if ((err = mp_mul_2d(a, 1, &t1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_add_d(&t1, 3, &t1)) != MP_OKAY) {
	    goto _ERR;
	}

	if ((err = mp_set_int(P, 1)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_sqr(q, &t2)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_mul(&t1, &t2, Q)) != MP_OKAY) {
	    goto _ERR;
	}
	if ((err = mp_copy(&t1, R)) != MP_OKAY) {
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

    if ((err = mp_acoth_binary_splitting(q, a, &t1, &p1, &q1, &r1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_acoth_binary_splitting(q, &t1, b, &p2, &q2, &r2)) != MP_OKAY) {
	goto _ERR;
    }
    //P = q2*p1 + r1*p2
    if ((err = mp_mul(&q2, &p1, &t1)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_mul(&r1, &p2, &t2)) != MP_OKAY) {
	goto _ERR;
    }
    if ((err = mp_add(&t1, &t2, P)) != MP_OKAY) {
	goto _ERR;
    }
    //Q =  q1*q2
    if ((err = mp_mul(&q1, &q2, Q)) != MP_OKAY) {
	goto _ERR;
    }
    //R = r1*r2
    if ((err = mp_mul(&r1, &r2, R)) != MP_OKAY) {
	goto _ERR;
    }

  _ERR:
    mp_clear_multi(&p1, &q1, &r1, &p2, &q2, &r2, &t1, &t2, &one, NULL);
    return err;
}

