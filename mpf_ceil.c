#include <tomfloat.h>
/* ceil function: floor(a) + 1 */
int mpf_ceil(mp_float * a, mp_float * b)
{
    int err;
    if ((err = mpf_floor(a, b)) != MP_OKAY) {
	return err;
    }
    if ((err = mpf_add_d(b, 1, b)) != MP_OKAY) {
	return err;
    }
    return MP_OKAY;
}

