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

/* 

Using the newton approximation y1 = y - (e^y - x)/e^y

Which converges quickly and we can reuse e^y so we only calc it once per loop

*/
int  mpf_ln(mp_float *a, mp_float *b)
{
    mp_float  oldval, tmpey, tmpy, val;
    int       itts, err;
    long      k;

    /* ensure positive */
    if (a->mantissa.sign == MP_NEG) {
       return MP_VAL;
    }

    /* easy out for 0 */
    if (mpf_iszero(a) == MP_YES) {
       return mpf_const_d(b, 1);
    }

    /* initialize temps */
    if ((err = mpf_init_multi(b->radix, &oldval, &tmpey, &tmpy, &val, NULL)) != MP_OKAY) {
       return err;
    }

    /* initial guess */
    if ((err = mpf_const_e(&val)) != MP_OKAY)                                           { goto __ERR; }
    if ((err = mpf_sqr(&val, &val)) != MP_OKAY)                                         { goto __ERR; }
    if ((err = mpf_inv(&val, &tmpey)) != MP_OKAY)                                       { goto __ERR; }
    if ((err = mpf_copy(a, &tmpy)) != MP_OKAY)                                          { goto __ERR; }
    if ((err = mpf_normalize_to(&tmpy, b->radix)) != MP_OKAY)                           { goto __ERR; }

    /* divide out e's */
    k = 0;
    while (mpf_cmp(&tmpy, &val) == MP_GT) {
        ++k;
        if ((err = mpf_mul(&tmpy, &tmpey, &tmpy)) != MP_OKAY)                           { goto __ERR; }
    }
    if ((err = mpf_const_d(&tmpy, k*2)) != MP_OKAY)                                     { goto __ERR; }

    /* number of iterations */
    itts = mpf_iterations(b);

    while (itts--) {
        if ((err = mpf_copy(&tmpy, &oldval)) != MP_OKAY)                                { goto __ERR; }

        /* get e^y and save it */
        if ((err = mpf_exp(&tmpy, &tmpey)) != MP_OKAY)                                  { goto __ERR; }
    
        /* now compute e^y - x */
        if ((err = mpf_sub(&tmpey, a, &val)) != MP_OKAY)                                { goto __ERR; }

        /* now compute (e^y - x) / e^y */
        if ((err = mpf_div(&val, &tmpey, &val)) != MP_OKAY)                             { goto __ERR; }

        /* y = y - (e^y - x)/e^y */
        if ((err = mpf_sub(&tmpy, &val, &tmpy)) != MP_OKAY)                             { goto __ERR; }

        if (mpf_cmp(&tmpy, &oldval) == MP_EQ) {
           break;
        }
    }
    mpf_exch(&tmpy, b);
__ERR: mpf_clear_multi(&oldval, &tmpey, &tmpy, &val, NULL);
   return err;
}
