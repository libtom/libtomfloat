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

/* using sin x == \sum_{n=0}^{\infty} ((-1)^n/(2n+1)!) * x^(2n+1) */
int  mpf_sin(mp_float *a, mp_float *b)
{
   mp_float tmpovern, tmp, tmpx, res, sqr;
   int      oddeven, err, itts;
   long     n;

   /* initialize temps */
   if ((err = mpf_init_multi(b->radix, &tmpx, &tmpovern, &tmp, &res, &sqr, NULL)) != MP_OKAY) {
      return err;
   }

   /* initlialize temps */
   /* this is the 1/n! which starts at 1 */
   if ((err = mpf_const_d(&tmpovern, 1)) != MP_OKAY)                                    { goto __ERR; }

   /* the square of the input, used to save multiplications later */
   if ((err = mpf_sqr(a, &sqr)) != MP_OKAY)                                             { goto __ERR; }

   /* tmpx starts at the input, so we copy and normalize */
   if ((err = mpf_copy(a, &tmpx)) != MP_OKAY)                                           { goto __ERR; }
   tmpx.radix = b->radix;
   if ((err = mpf_normalize(&tmpx)) != MP_OKAY)                                         { goto __ERR; }

   /* the result starts off at a as we skip a term in the series */
   if ((err = mpf_copy(&tmpx, &res)) != MP_OKAY)                                        { goto __ERR; }

   /* this is the denom counter.  Goes up by two per pass */
   n       = 1;

   /* we alternate between adding and subtracting */
   oddeven = 1;

   /* get number of iterations */
   itts = mpf_iterations(b);

   while (itts-- > 0) {
       /* compute 1/(2n)! from 1/(2(n-1))! by multiplying by (1/n)(1/(n+1)) */
       if ((err = mpf_const_d(&tmp, ++n)) != MP_OKAY)                                   { goto __ERR; }
       if ((err = mpf_inv(&tmp, &tmp)) != MP_OKAY)                                      { goto __ERR; }
       if ((err = mpf_mul(&tmpovern, &tmp, &tmpovern)) != MP_OKAY)                      { goto __ERR; }
       /* we do this twice */
       if ((err = mpf_const_d(&tmp, ++n)) != MP_OKAY)                                   { goto __ERR; }
       if ((err = mpf_inv(&tmp, &tmp)) != MP_OKAY)                                      { goto __ERR; }
       if ((err = mpf_mul(&tmpovern, &tmp, &tmpovern)) != MP_OKAY)                      { goto __ERR; }

       /* now multiply sqr into tmpx */
       if ((err = mpf_mul(&tmpx, &sqr, &tmpx)) != MP_OKAY)                              { goto __ERR; }

       /* now multiply the two */
       if ((err = mpf_mul(&tmpx, &tmpovern, &tmp)) != MP_OKAY)                          { goto __ERR; }

       /* now depending on if this is even or odd we add/sub */
       oddeven ^= 1;
       if (oddeven  == 1) {
          if ((err = mpf_add(&res, &tmp, &res)) != MP_OKAY)                             { goto __ERR; }
       } else {
          if ((err = mpf_sub(&res, &tmp, &res)) != MP_OKAY)                             { goto __ERR; }
       }
   }
   mpf_exch(&res, b);
__ERR: mpf_clear_multi(&tmpx, &tmpovern, &tmp, &res, &sqr, NULL);
   return err;
}
