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

void mpf_clear_multi(mp_float *a, ...)
{
    mp_float* next_mp = a;
    va_list args;
    va_start(args, a);
    while (next_mp != NULL) {
        mpf_clear(next_mp);
        next_mp = va_arg(args, mp_float*);
    }
    va_end(args);
}
