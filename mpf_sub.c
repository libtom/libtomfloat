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

int  mpf_sub(mp_float *a, mp_float *b, mp_float *c)
{
   int err;
   mp_float tmp;
   long diff;

   if (a->exp < b->exp) {
      /* tmp == a normalize to b's exp */
      if ((err = mpf_init_copy(a, &tmp)) != MP_OKAY) {
         return err;
      }

      /* now make tmp.exp == b.exp by dividing tmp by 2^(b.exp - tmp.exp) */
      diff = b->exp - tmp.exp;
      tmp.exp = b->exp;
      if ((err = mp_div_2d(&(tmp.mantissa), diff, (&tmp.mantissa), NULL)) != MP_OKAY)  { goto __TMP; }
      if ((err = mp_sub(&(tmp.mantissa), &(b->mantissa), &(c->mantissa))) != MP_OKAY)  { goto __TMP; }
      c->exp = b->exp;
   } else {
      /* tmp == b normalize to a's radix */
      if ((err = mpf_init_copy(b, &tmp)) != MP_OKAY) {
         return err;
      }

      /* now make tmp.exp == a.exp by dividing tmp by 2^(a.exp - tmp.exp) */
      diff = a->exp - tmp.exp;
      tmp.exp = a->exp;
      if ((err = mp_div_2d(&(tmp.mantissa), diff, (&tmp.mantissa), NULL)) != MP_OKAY)   { goto __TMP; }
      if ((err = mp_sub(&(a->mantissa), &(tmp.mantissa), &(c->mantissa))) != MP_OKAY)  { goto __TMP; }
      c->exp = a->exp;
   }

   err    = mpf_normalize(c);

__TMP:  mpf_clear(&tmp);
   return err;
}
