#include <tomfloat.h>

static mp_float mpf_eps;

static long mpf_eps_precision;

int mpf_const_eps(mp_float * a)
{
    int err;
    long eps;

    err = MP_OKAY;

    if (mpf_eps_precision > 0 && a == NULL) {
        mpf_clear(&mpf_eps);
        mpf_eps_precision = 0;
        return err;
    }
    if (mpf_eps_precision >= a->radix) {
        eps = a->radix;
        if ((err = mpf_copy(&mpf_eps, a)) != MP_OKAY) {
            return err;
        }
        return mpf_normalize_to(a, eps);
    } else {
        if (mpf_eps_precision == 0) {
            if (mpf_eps.mantissa.dp != NULL) {
        		mpf_clear(&mpf_eps);
        	}
            if ((err = mpf_init(&mpf_eps, a->radix)) != MP_OKAY) {
                return err;
            }
        }
        if ((err = mpf_const_d(&mpf_eps, 1)) != MP_OKAY) {
            return err;
        }
        mpf_eps.exp = -(2 * mpf_eps.radix);
        if ((err = mpf_copy(&mpf_eps, a)) != MP_OKAY) {
            return err;
        }
    }
    return MP_OKAY;
}

