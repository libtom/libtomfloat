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

/* compute 1/x by (1/sqrt(x))^2 */
int  mpf_inv(mp_float *a, mp_float *b)
{
   int err, sign;

   /* get sign of input */
   sign = a->mantissa.sign;

   /* force to positive */
   a->mantissa.sign = MP_ZPOS;

   /* compute 1/sqrt(a) */
   if ((err = mpf_invsqrt(a, b)) != MP_OKAY) {
      return err;
   }

   /* square 1/sqrt(a) to get 1/a */
   err = mpf_sqr(b, b);

   /* now restore the sign */
   b->mantissa.sign = a->mantissa.sign = sign;

   return err;
}

