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

int  mpf_init_multi(long radix, mp_float *a, ...)
{
    mp_err res = MP_OKAY;      /* Assume ok until proven otherwise */
    int n = 0;                 /* Number of ok inits */
    mp_float* cur_arg = a;
    va_list args;

    va_start(args, a);        /* init args to next argument from caller */
    while (cur_arg != NULL) {
        if (mpf_init(cur_arg, radix) != MP_OKAY) {
            /* Oops - error! Back-track and mp_clear what we already
               succeeded in init-ing, then return error.
            */
            va_list clean_args;
            
            /* end the current list */
            va_end(args);
            
            /* now start cleaning up */            
            cur_arg = a;
            va_start(clean_args, a);
            while (n--) {
                mpf_clear(cur_arg);
                cur_arg = va_arg(clean_args, mp_float*);
            }
            va_end(clean_args);
            res = MP_MEM;
            break;
        }
        n++;
        cur_arg = va_arg(args, mp_float*);
    }
    va_end(args);
    return res;                /* Assumed ok, if error flagged above. */
}
