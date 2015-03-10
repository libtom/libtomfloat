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
#include <stdarg.h>

int  mpf_normalize_to_multi(long radix, mp_float *a, ...)
{
    mp_err res = MP_OKAY;      /* Assume ok until proven otherwise */
    int n = 0;                 /* Number of ok inits */
    mp_float* cur_arg = a;
    va_list args;

    va_start(args, a);        /* init args to next argument from caller */
    while (cur_arg != NULL) {
        if (mpf_normalize_to(cur_arg, radix) != MP_OKAY) {
            res = MP_MEM;
            break;
        }
        n++;
        cur_arg = va_arg(args, mp_float*);
    }
    va_end(args);
    return res;                /* Assumed ok, if error flagged above. */
}
