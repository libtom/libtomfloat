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

/* tan = sin/cos, prolly fix this up later */
int  mpf_tan(mp_float *a, mp_float *b)
{
   int       err;
   mp_float  tmp;

   /* get cos of a */
   if ((err = mpf_init(&tmp, b->radix)) != MP_OKAY) {
      return err;
   }
   if ((err = mpf_cos(a, &tmp)) != MP_OKAY)                     { goto __ERR; }

   /* now make it upside down ;-) (this should catch domain errors too) */
   if ((err = mpf_inv(&tmp, &tmp)) != MP_OKAY)                  { goto __ERR; }

   /* put sin in b */
   if ((err = mpf_sin(a, b)) != MP_OKAY)                        { goto __ERR; }

   /* multiply the two and we done */
	   err = mpf_mul(b, &tmp, b);

__ERR: mpf_clear(&tmp);
   return err;
}

