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

int  mpf_mul(mp_float *a, mp_float *b, mp_float *c)
{
   int err;
 
   if ((err = mp_mul(&(a->mantissa), &(b->mantissa), &(c->mantissa))) != MP_OKAY) {
      return err;
   }
   c->exp = a->exp + b->exp;
   return mpf_normalize(c);
}
