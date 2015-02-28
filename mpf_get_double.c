#include "tomfloat.h"
#include <math.h>
#include <float.h>
// set the mp_float to the value of the double
int mpf_get_double(double d, mp_float * a)
{

    int err, sign, e;
    mp_int tmp;

    sign = (d < 0) ? MP_NEG : MP_ZPOS;
    d = fabs(d);
#if __STDC_VERSION__ >= 199901L
    if (isinf(d) || isnan(d)) {
	if (isinf(d)) {
	    if ((err = mpf_const_inf(a, sign)) != MP_OKAY) {
		return err;
	    }
	    return MP_OKAY;
	}
	if (isnan(d)) {
	    if ((err = mpf_const_nan(a)) != MP_OKAY) {
		return err;
	    }
	    a->mantissa.sign = sign;
	    return MP_OKAY;
	}
    }
#endif
    if (d == 0) {
	a->mantissa.sign = sign;
	return MP_OKAY;
    }

    if ((err = mp_init(&tmp)) != MP_OKAY) {
	return err;
    }
    // part the fraction from the exponent
    d = frexp(d,&e);
    // multiply by 2^m whith m the number of bits in the mantissa including the
    // mutual bit
    d *= pow(2.0, DBL_MANT_DIG);
    d = round(d);
    // convert the double to a mp_int (rounding mode 0 = to nearest)-unexpensive
    if ((err = mp_set_double(&tmp, d, 0) ) != MP_OKAY) {
	return err;
    }
    // convert the mp_int to a mp_float-cheap
    if ((err =  mpf_from_mp_int(&tmp, a)) != MP_OKAY) {	return err;    }
    //  correct exponent-costs next to nothing
    a->exp -= DBL_MANT_DIG - e;
    mp_clear(&tmp);
    // don't forget the sign which seems to happen quite often
    a->mantissa.sign = sign;
    //TODO: the way we do it is, although very simpel, very sensitive to
    //      rounding errors. We could do it directly from the innards of the
    //      double (as has been done in the JavaScript version of this lib)
    //      but endinaness will need a large amount of error-prone branches.
    //      Especially for the case of some modern processors which can switch
    //      endianess willi-nilly. That needs something in software to check
    //      it at run time.
    return MP_OKAY;
}

