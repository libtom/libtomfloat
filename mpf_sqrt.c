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

/* compute using sqrt(x) = y_{n+1} = 0.5 * (y_{n} + x/y_{n}) */
int  mpf_sqrt(mp_float *a, mp_float *b)
{
   mp_float oldval, tmp, res;
   int      err, itts;

   /* make sure it's positive */
   if (a->mantissa.sign == MP_NEG) {
      return MP_VAL;
   }

   /* let's roll */
   if ((err = mpf_init_multi(b->radix, &oldval, &tmp, &res, NULL)) != MP_OKAY) {
      return err;
   }

   /* get copy of a and make reasonable guestimte towards the sqrt */
   if ((err = mpf_copy(a, &res)) != MP_OKAY)                                     { goto __ERR; }
   mp_rshd(&(res.mantissa), res.mantissa.used >> 1);
   res.exp   = res.exp / 2;
   if ((err = mpf_normalize_to(&res, b->radix)) != MP_OKAY)                      { goto __ERR; }

   /* number of iterations */
   itts = mpf_iterations(b);

   while (itts--) {
        if ((err = mpf_copy(&res, &oldval)) != MP_OKAY)                          { goto __ERR; }

        /* compute x/res */
        if ((err = mpf_div(a, &res, &tmp)) != MP_OKAY)                           { goto __ERR; }
        /* res + x/res */
        if ((err = mpf_add(&res, &tmp, &res)) != MP_OKAY)                        { goto __ERR; }
        /* 0.5 * (res + x/res) */
        if ((err = mpf_div_2(&res, &res)) != MP_OKAY)                            { goto __ERR; }

        if (mpf_cmp(&res, &oldval) == MP_EQ) {
           break;
        }
   }

   mpf_exch(&res, b);
__ERR: mpf_clear_multi(&oldval, &tmp, &res, NULL);
   return err;
}

