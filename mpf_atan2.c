#include <tomfloat.h>

/*
    atan2(a,b) = atan(a/b);

*/
int mpf_atan2(mp_float * a, mp_float * b, mp_float * c)
{
    int err;
    // TODO: implement IEEE-754. But be warned: it's a looong list!
    if ((err = mpf_div(a, b, c)) != MP_OKAY) {
	return err;
    }
    return mpf_atan(c, c);
}

