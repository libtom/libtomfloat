#include <tomfloat.h>
/* Converts a mp_int to a mp_float */
// TODO: checks and balances
int mpf_from_mp_int(mp_int * a, mp_float * b)
{
    int err;
    if ((err = mp_copy(a, &b->mantissa)) != MP_OKAY) {
	return err;
    }
    b->exp = 0;
    mpf_normalize(b);
    return MP_OKAY;
}
