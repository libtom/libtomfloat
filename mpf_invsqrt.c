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

/* using newtons method we have 1/sqrt(x) = Y_{n+1} = y_n * ((3 - xy^2_n)/2) */
int  mpf_invsqrt(mp_float *a, mp_float *b)
{
   mp_float oldval, tmp1, tmp2, const_3;
   int err, itts;
   long expo, nbits;

   /* ensure a is not zero or negative */
   if ((mp_iszero(&(a->mantissa)) == MP_YES) || (a->mantissa.sign == MP_NEG)) {
      return MP_VAL;
   }

   /* get number of iterations */
   itts = mpf_iterations(b);
 
   /* init temps */
   if ((err = mpf_init_multi(b->radix, &oldval, &tmp1, &tmp2, &const_3, NULL)) != MP_OKAY) {
      return err;
   }

   /* const_3 */
   if ((err = mpf_const_d(&const_3, 3)) != MP_OKAY)                               { goto __ERR; }

   /* tmp1 == reasonable guess at sqrt */
   if ((err = mpf_copy(a, &tmp1)) != MP_OKAY)                                     { goto __ERR; }
   expo = tmp1.exp;
   if (expo < -tmp1.radix) {
      nbits = -tmp1.radix - expo;
      tmp1.exp = (expo + ( nbits + (nbits/2) ) );
   }
   else if (expo > -tmp1.radix) {
      nbits = expo - -tmp1.radix ;
      tmp1.exp = (expo - ( nbits + (nbits/2) ) );
   }
   else {
      /* do nothing, near enough? */
      tmp1.exp = expo;
   }

   while (itts-- > 0) {
       /* grap copy of tmp1 for early out */
       if ((err = mpf_copy(&tmp1, &oldval)) != MP_OKAY)                           { goto __ERR; }

       /* first tmp2 = y^2 == tmp1^2 */
       if ((err = mpf_sqr(&tmp1, &tmp2)) != MP_OKAY)                              { goto __ERR; }
       /* multiply by x, tmp1 * a */
       if ((err = mpf_mul(&tmp2, a, &tmp2)) != MP_OKAY)                           { goto __ERR; }
       /* 3 - xy^2_n == 3 - tmp1 */
       if ((err = mpf_sub(&const_3, &tmp2, &tmp2)) != MP_OKAY)                    { goto __ERR; }
       /* halve it */
       if ((err = mpf_div_2(&tmp2, &tmp2)) != MP_OKAY)                            { goto __ERR; }
       /* multiply by y_n and feedback */
       if ((err = mpf_mul(&tmp1, &tmp2, &tmp1)) != MP_OKAY)                       { goto __ERR; }

       /* early out if stable */
       if (mpf_cmp(&oldval, &tmp1) == MP_EQ) {
          break;
       }
   }

   mpf_exch(&tmp1, b);
__ERR: mpf_clear_multi(&oldval, &tmp1, &tmp2, &const_3, NULL);
   return err;
}

