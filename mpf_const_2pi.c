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

int  mpf_const_2pi(mp_float *a)
{
   int err;
   if ((err = mpf_const_pi(a)) != MP_OKAY) {
      return err;
   }
   if ((err = mpf_inv(a, a)) != MP_OKAY) {
      return err;
   }
   return mpf_mul_2(a, a);
}
                /* 2/Pi     */

