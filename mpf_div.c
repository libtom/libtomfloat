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

int  mpf_div(mp_float *a, mp_float *b, mp_float *c)
{
   mp_float tmp;
   int err;

   /* ensure b is not zero */
   if (mp_iszero(&(b->mantissa)) == MP_YES) {
      return MP_VAL;
   }

   /* find 1/b */
   if ((err = mpf_init(&tmp, c->radix)) != MP_OKAY) {
      return err;
   }
   if ((err = mpf_inv(b, &tmp)) != MP_OKAY)                           { goto __ERR; }
   
   /* now multiply */
   err = mpf_mul(&tmp, a, c);

__ERR: mpf_clear(&tmp);
   return err;
}

