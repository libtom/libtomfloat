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

int  mpf_normalize(mp_float *a)
{
   long cb, diff;

   /* sanity */
   if (a->radix < 2) {
      return MP_VAL;
   }

   cb = mp_count_bits(&(a->mantissa));
   if (cb > a->radix) {
      diff    = cb - a->radix;
      a->exp += diff;
      return mp_div_2d(&(a->mantissa), diff, &(a->mantissa), NULL);
   } else if (cb < a->radix) {
      if (mp_iszero(&(a->mantissa)) == MP_YES) {
         return mpf_const_0(a);
      } else {
         diff    = a->radix - cb;
         a->exp -= diff;
         return mp_mul_2d(&(a->mantissa), diff, &(a->mantissa));
      }
   }
   return MP_OKAY;
}

