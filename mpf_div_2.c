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

int  mpf_div_2(mp_float *a, mp_float *b)
{
   int err;
   
   /* |b| = |a|/2 */
   if ((err = mp_copy(&(a->mantissa), &(b->mantissa))) != MP_OKAY) {
      return err;
   }
   b->exp = a->exp - 1;
   return mpf_normalize(b);
}
