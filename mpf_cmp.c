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

int  mpf_cmp(mp_float *a,   mp_float *b)
{
   int za, zb, sa, sb;

   /* if one is zero than we early out */
   za = mp_iszero(&(a->mantissa));
   sa = a->mantissa.sign;
   zb = mp_iszero(&(b->mantissa));
   sb = b->mantissa.sign;

   if (za == MP_YES && zb == MP_NO) {
      /* result depends on b */
      if (sb == MP_NEG) {
         return MP_GT;
      } else {
         return MP_LT;
      }
   } else if (za == MP_NO && zb == MP_YES) {
      /* result depends on a */
      if (sa == MP_NEG) {
         return MP_LT;
      } else {
         return MP_GT;
      }
   }

   /* compare the signs */
   if (sa == MP_NEG && sb == MP_ZPOS) {
      return MP_LT;
   } else if (sa == MP_ZPOS && sb == MP_NEG) {
      return MP_GT;
   }

   /* they're both non-zero, the same sign and normalized, compare the exponents */
   if (a->exp > b->exp) {
      return (sa == MP_NEG) ? MP_LT : MP_GT;
   } else if (a->exp < b->exp) {
      return (sa == MP_NEG) ? MP_GT : MP_LT;
   }

   /* same exponent and sign, compare mantissa */
   return mp_cmp(&(a->mantissa), &(b->mantissa));
}

