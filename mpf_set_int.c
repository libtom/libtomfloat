#include <tomfloat.h>

int mpf_set_int(mp_float * a, int d)
{
    int err;
    err = MP_OKAY;
    if ((err = mp_set_int(&(a->mantissa), abs(d))) != MP_OKAY) {
	return err;
    }
    a->mantissa.sign = (d < 0) ? MP_NEG : MP_ZPOS;
    if ((err = mpf_from_mp_int(&(a->mantissa), a)) != MP_OKAY) {
	return err;
    }
    return err;
}

