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

/* we have e^x, so why not write a^b as e^(lna * b) ;-) w00t w00t */
int  mpf_pow(mp_float *a, mp_float *b, mp_float *c)
{
   mp_float  exponent;
   int       err;

   if ((err = mpf_init(&exponent, c->radix)) != MP_OKAY) {
      return err;
   }

   /* get ln of a */
   if ((err = mpf_ln(a, &exponent)) != MP_OKAY)                      { goto __ERR; }

   /* multiply it by b */
   if ((err = mpf_mul(&exponent, b, &exponent)) != MP_OKAY)          { goto __ERR; }

   /* now evaluate it */
   err = mpf_exp(&exponent, c);

__ERR: mpf_clear(&exponent);
   return err;
}
