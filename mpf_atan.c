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

/* y = y - (tan(y) - x)/(tan(y+0.1)-tan(x)) */

int  mpf_atan(mp_float *a, mp_float *b)
{
   mp_float oldval, tmp, tmpx, res, sqr;
   int      oddeven, ires, err, itts;
   long     n;

   /* ensure -1 <= a <= 1 */
   if ((err = mpf_cmp_d(a, -1, &ires)) != MP_OKAY) {
      return err; 
   }
   if (ires == MP_LT) {
      return MP_VAL;
   }

   if ((err = mpf_cmp_d(a, 1, &ires)) != MP_OKAY) {
      return err; 
   }
   if (ires == MP_GT) {
      return MP_VAL;
   }

   /* easy out if a == 0 */
   if (mpf_iszero(a) == MP_YES) {
      return mpf_const_d(b, 1);
   }
   
   /* now a != 0 */

   /* initialize temps */
   if ((err = mpf_init_multi(b->radix, &oldval, &tmpx, &tmp, &res, &sqr, NULL)) != MP_OKAY) {
      return err;
   }

   /* initlialize temps */
   /* res = 0 */
   /* tmpx = 1/a */
   if ((err = mpf_inv(a, &tmpx)) != MP_OKAY)                                            { goto __ERR; }

   /* sqr = a^2 */
   if ((err = mpf_sqr(a, &sqr)) != MP_OKAY)                                             { goto __ERR; }

   /* this is the denom counter.  Goes up by two per pass */
   n       = 1;

   /* we alternate between adding and subtracting */
   oddeven = 0;

   /* get number of iterations */
   itts = mpf_iterations(b);

   while (itts-- > 0) {
       if ((err = mpf_copy(&res, &oldval)) != MP_OKAY)                                  { goto __ERR; }

       /* compute 1/(2n-1) */
       if ((err = mpf_const_d(&tmp, (2*n++ - 1))) != MP_OKAY)                           { goto __ERR; }
       if ((err = mpf_inv(&tmp, &tmp)) != MP_OKAY)                                      { goto __ERR; }

       /* now multiply a into tmpx twice */
       if ((err = mpf_mul(&tmpx, &sqr, &tmpx)) != MP_OKAY)                              { goto __ERR; }

       /* now multiply the two */
       if ((err = mpf_mul(&tmpx, &tmp, &tmp)) != MP_OKAY)                               { goto __ERR; }

       /* now depending on if this is even or odd we add/sub */
       oddeven ^= 1;
       if (oddeven  == 1) {
          if ((err = mpf_add(&res, &tmp, &res)) != MP_OKAY)                             { goto __ERR; }
       } else {
          if ((err = mpf_sub(&res, &tmp, &res)) != MP_OKAY)                             { goto __ERR; }
       }

       if (mpf_cmp(&oldval, &res) == MP_EQ) {
          break;
       }
   }
   mpf_exch(&res, b);
__ERR: mpf_clear_multi(&oldval, &tmpx, &tmp, &res, &sqr, NULL);
   return err;

}
