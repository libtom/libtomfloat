#include <tomfloat.h>

int mpf_sinh(mp_float * a, mp_float * b)
{

    mp_float x, t1, t2;
    int err, sign;
    long size;

    if (mpf_iszero(a)) {
	return mpf_const_0(b);
    }
    err = MP_OKAY;
    // check size of input (sinh grows rapidly!)
    // global limit is at about +/- 1e11 but YMMV and if it varies: please 
    // contact the author of this file with the machine data
    size = a->radix + a->exp;
    if (size > 0) {
	// upper limit of exp(x) which is about 2.1144e43,429,448,190
	if (size / 3.3 > 1e11) {
	    // raise invalid (TLOSS in SysV)?
	    // feraiseexcept(FE_INVALID)
	    // mpf_errno = E_DOM
	    // or inexact?
	    if ((err = mpf_const_nan(b)) != MP_OKAY) {
		return err;
	    }
	    return MP_VAL;
	}
	// use (exp(x) - exp(-x))/2 otherwise
	if ((err = mpf_init_multi(a->radix, &t1, &t2, NULL)) != MP_OKAY) {
	    return err;
	}
	if ((err = mpf_exp(a, &t1)) != MP_OKAY) {
	    return err;
	}
        // The limit is much smaller, actually, but to be on the safer side...
	if (size <= a->radix) {
	    a->mantissa.sign = (a->mantissa.sign == MP_NEG) ? MP_ZPOS : MP_NEG;
	    if ((err = mpf_exp(a, &t2)) != MP_OKAY) {
		a->mantissa.sign =
		    (a->mantissa.sign == MP_NEG) ? MP_ZPOS : MP_NEG;
		goto _SERR;
	    }
	    a->mantissa.sign = (a->mantissa.sign == MP_NEG) ? MP_ZPOS : MP_NEG;
	    if ((err = mpf_sub(&t1, &t2, &t1)) != MP_OKAY) {
		goto _SERR;
	    }
	}
	t1.exp -= 1;
	t1.mantissa.sign = a->mantissa.sign;
	if ((err = mpf_copy(&t1, b)) != MP_OKAY) {
	    goto _SERR;
	}

      _SERR:
	mpf_clear_multi(&t1, &t2, NULL);
	return err;
    }
    sign = a->mantissa.sign;
    if ((err = mpf_init(&x, a->radix)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_abs(a, &x)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_sincos(&x, &x, 0, 0, 1)) != MP_OKAY) {
	goto _ERR;
    }

    x.mantissa.sign = sign;
    if ((err = mpf_copy(&x, b)) != MP_OKAY) {
	goto _ERR;
    }
  _ERR:
    mpf_clear(&x);
    return err;
}

