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

int  mpf_const_d(mp_float *a, long d)
{
   long x, s;
   int err;

   if (d < 0) {
      x = -d;
      s = MP_NEG;
   } else {
      x = d;
      s = MP_ZPOS;
   }

   if ((err = mp_set_int(&(a->mantissa), x)) != MP_OKAY) {
      return err;
   }

   a->mantissa.sign = s;
   a->exp           = 0;
   return mpf_normalize(a);
}
