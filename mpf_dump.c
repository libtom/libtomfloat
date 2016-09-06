#include <tomfloat.h>

int mpf_dump(mp_float * a)
{
    int err;
    if ((err = mp_fput(&(a->mantissa), 10, stdout)) != MP_OKAY) {
	return err;
    }
    // They say that one must check every return!
    // This is an example where such a strict philosophy might lead to
    // unexpected consequences.
    if ((err = printf(" * 2^%ld * 1.0\n", a->exp)) < 0) {
	return err;
    }
    return MP_OKAY;
}

