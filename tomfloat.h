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
#ifndef TF_H_
#define TF_H_

#include <tommath.h>

/* this is mp_float type */
typedef struct {
     mp_int mantissa;
     long   radix,       /* how many bits for mantissa */
            exp;         /* current exponent, e.g. mantissa * 2^exp == number  */
} mp_float;

/* initializers */
int  mpf_init(mp_float *a, long radix);
void mpf_clear(mp_float *a);

int  mpf_init_multi(long radix, mp_float *a, ...);
void mpf_clear_multi(mp_float *a, ...);

int  mpf_init_copy(mp_float *a, mp_float *b);

int  mpf_copy(mp_float *src, mp_float *dest);
void mpf_exch(mp_float *a, mp_float *b);

/* maintainers */
int  mpf_normalize(mp_float *a);
int  mpf_normalize_to(mp_float *a, long radix);
int  mpf_iterations(mp_float *a);

/* constants */
int  mpf_const_0(mp_float *a);                  /* valid zero */
int  mpf_const_d(mp_float *a, long d);          /* valid d */
int  mpf_const_ln_d(mp_float *a, long b);       /* a = ln(b)     */
int  mpf_const_sqrt_d(mp_float *a, long b);     /* a = sqrt(b);  */

/* math constants as they appear in math.h */
int  mpf_const_e(mp_float *a);                  /* e             */
int  mpf_const_l2e(mp_float *a);                /* log_2 e       */
int  mpf_const_l10e(mp_float *a);               /* log_10 e      */
int  mpf_const_le2(mp_float *a);                /* log_e 2       */
int  mpf_const_pi(mp_float *a);                 /* Pi            */
int  mpf_const_pi2(mp_float *a);                /* Pi/2          */
int  mpf_const_pi4(mp_float *a);                /* Pi/4          */
int  mpf_const_1pi(mp_float *a);                /* 1/Pi          */
int  mpf_const_2pi(mp_float *a);                /* 2/Pi          */
int  mpf_const_2rpi(mp_float *a);               /* 2/sqrt(Pi)    */
int  mpf_const_r2(mp_float *a);                 /* sqrt(2)       */
int  mpf_const_1r2(mp_float *a);                /* 1/sqrt(2)     */

/* sign operators */
int  mpf_abs(mp_float *a, mp_float *b);         /* absolute */
int  mpf_neg(mp_float *a, mp_float *b);         /* negation */

/* basic math */
int  mpf_mul_2(mp_float *a, mp_float *b);              /* b = 2a       */
int  mpf_div_2(mp_float *a, mp_float *b);              /* b = a/2      */
int  mpf_add(mp_float *a, mp_float *b, mp_float *c);   /* c = a + b    */
int  mpf_sub(mp_float *a, mp_float *b, mp_float *c);   /* c = a - b    */
int  mpf_mul(mp_float *a, mp_float *b, mp_float *c);   /* c = a * b    */
int  mpf_div(mp_float *a, mp_float *b, mp_float *c);   /* c = a / b    */
int  mpf_sqr(mp_float *a, mp_float *b);                /* b = a^2      */

int  mpf_add_d(mp_float *a, long b, mp_float *c);      /* c = a + b    */
int  mpf_sub_d(mp_float *a, long b, mp_float *c);      /* c = a - b    */
int  mpf_mul_d(mp_float *a, long b, mp_float *c);      /* c = a * b    */
int  mpf_div_d(mp_float *a, long b, mp_float *c);      /* c = a / b    */

/* compares */
int  mpf_cmp(mp_float *a,   mp_float *b);
int  mpf_cmp_d(mp_float *a, long b, int *res);
#define mpf_iszero(a) mp_iszero(&((a)->mantissa))

/* Algebra */
int  mpf_exp(mp_float *a, mp_float *b);                /* b = e^a       */
int  mpf_pow(mp_float *a, mp_float *b, mp_float *c);   /* c = a^b       */
int  mpf_ln(mp_float *a, mp_float *b);                 /* b = ln a      */
int  mpf_invsqrt(mp_float *a, mp_float *b);            /* b = 1/sqrt(a) */
int  mpf_inv(mp_float *a, mp_float *b);                /* b = 1/a       */
int  mpf_sqrt(mp_float *a, mp_float *b);               /* b = sqrt(a)   */

/* Trig */
int  mpf_cos(mp_float *a, mp_float *b);                /* b = cos(a)    */
int  mpf_sin(mp_float *a, mp_float *b);                /* b = sin(a)    */
int  mpf_tan(mp_float *a, mp_float *b);                /* b = tan(a)    */
int  mpf_acos(mp_float *a, mp_float *b);               /* b = acos(a)   */
int  mpf_asin(mp_float *a, mp_float *b);               /* b = asin(a)   */
int  mpf_atan(mp_float *a, mp_float *b);               /* b = atan(a)   */

/* ASCII <=> mp_float conversions */
char *mpf_to_string(mp_float *a, mp_digit radix);
int mpf_from_string(mp_float *a, const char *str, mp_digit radix);

#endif
