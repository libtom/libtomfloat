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


/* compute b = e^a using e^x == \sum_{n=0}^{\infty} {1 \over n!}x^n */
int  mpf_exp(mp_float *a, mp_float *b)
{
   mp_float  tmpx, tmpovern, tmp, res; 
   int       err, itts;
   long      n;

   /* initialize temps */
   if ((err = mpf_init_multi(b->radix, &tmpx, &tmpovern, &tmp, &res, NULL)) != MP_OKAY) {
      return err;
   }

   /* initlialize temps */
   /* all three start at one */
   if ((err = mpf_const_d(&res, 1)) != MP_OKAY)                                         { goto __ERR; }
   if ((err = mpf_const_d(&tmpovern, 1)) != MP_OKAY)                                    { goto __ERR; }
   if ((err = mpf_const_d(&tmpx, 1)) != MP_OKAY)                                        { goto __ERR; }
   n = 1;

   /* get number of iterations */
   itts = mpf_iterations(b);

   while (itts-- > 0) {
       /* compute 1/n! as 1/(n-1)! * 1/n */
// hack: this won't be portable for n>127
       if ((err = mpf_const_d(&tmp, n++)) != MP_OKAY)                                   { goto __ERR; }
       if ((err = mpf_inv(&tmp, &tmp)) != MP_OKAY)                                      { goto __ERR; }
       if ((err = mpf_mul(&tmp, &tmpovern, &tmpovern)) != MP_OKAY)                      { goto __ERR; }
 
       /* compute x^n as x^(n-1) * x */
       if ((err = mpf_mul(&tmpx, a, &tmpx)) != MP_OKAY)                                 { goto __ERR; }

       /* multiply and sum them */
       if ((err = mpf_mul(&tmpovern, &tmpx, &tmp)) != MP_OKAY)                          { goto __ERR; }
       if ((err = mpf_add(&tmp, &res, &res)) != MP_OKAY)                                { goto __ERR; }
   }

   mpf_exch(&res, b);
__ERR: mpf_clear_multi(&tmpx, &tmpovern, &tmp, &res, NULL);
   return err;
}

