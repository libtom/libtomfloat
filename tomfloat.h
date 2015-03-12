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
// The math-functions gets used a lot, especially as initial values for the
// Newton based algorithms
#include <math.h>

/* this is mp_float type */
typedef struct {
     mp_int mantissa;
     long   radix,       /* how many bits for mantissa */
            exp;         /* current exponent, e.g. mantissa * 2^exp == number  */
} mp_float;



/* Global radix value */
extern long mpf_global_radix;
/* Global error holder?  */
extern int mpf_errno;

/* Some cutoffs */
extern int MPF_LOG_AGM_1_CUTOFF;
extern int MPF_LOG_AGM_2_CUTOFF;
extern long MPF_LOG_AGM_REDUX_CUTOFF;
/* Handling of the precision set above */
long mpf_getprecision();
long mpf_getdecimalprecision();
/* A threadsafe variant useing a mutex is available if MPF_USE_THREADS
   is defined */
void mpf_setprecision(long r);

/* A helper for the radix conversion */
long mpf_getdecimalexponent(long exp);

/* initializers */
int  mpf_init(mp_float *a, long radix);
void mpf_clear(mp_float *a);

int  mpf_init_multi(long radix, mp_float *a, ...);
void mpf_clear_multi(mp_float *a, ...);

int  mpf_init_copy(mp_float *a, mp_float *b);

int  mpf_copy(mp_float *src, mp_float *dest);
void mpf_exch(mp_float *a, mp_float *b);

/* maintainers/misc */
int  mpf_normalize(mp_float *a);
int  mpf_normalize_to(mp_float *a, long radix);
int  mpf_normalize_to_multi(long radix, mp_float *a, ...);
int  mpf_iterations(mp_float *a);
int  mpf_from_mp_int(mp_int * a, mp_float * b);

/* constants */
int  mpf_const_inf(mp_float *a, int sign);      /* set to +/- infinity */
int  mpf_const_nan(mp_float *a);                /* set to NaN */
int  mpf_const_0(mp_float *a);                  /* valid zero */
int  mpf_const_d(mp_float *a, long d);          /* valid d */
int  mpf_const_ln_d(mp_float *a, long b);       /* a = ln(b)     */
int  mpf_const_sqrt_d(mp_float *a, long b);     /* a = sqrt(b);  */

int  mpf_const_gamma(mp_float * a);             /* Euler-Gamma ~0.577215 */

/* math constants as they appear in math.h */
int  mpf_const_e(mp_float *a);                  /* e             */
int  mpf_const_l2e(mp_float *a);                /* log_2 e       */
int  mpf_const_l10e(mp_float *a);               /* log_10 e      */
int  mpf_const_le2(mp_float *a);                /* log_e 2       */
int  mpf_const_le10(mp_float *a);               /* log_e 10      */
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

#define mpf_isnan(a) \
( \
a->mantissa.dp[a->mantissa.used - 1] ==  (mp_digit) (1) \
&& \
a->exp == 1 << ((sizeof(long) * CHAR_BIT) - 1)\
)

#define mpf_isinf(a) \
( \
a->mantissa.dp[a->mantissa.used - 1] ==  (mp_digit) (0) \
&& \
a->exp == 1 << ((sizeof(long) * CHAR_BIT) - 1)\
)

#define mpf_isdouble(a) (-(1021 + a->radix) <= a->exp && a->exp <= (1024 - a->radix))

#define mpf_isfraction(a) (mpf_iszero(a) || (a->exp <= -a->radix))

/* Algebra */
int  mpf_exp(mp_float *a, mp_float *b);                /* b = e^a       */
int  mpf_pow(mp_float *a, mp_float *b, mp_float *c);   /* c = a^b       */
int  mpf_pow_d(mp_float *a, long e, mp_float *c);      /* c = a^b       */
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

/* Bases 2, 10, and 16 from IEEE-754 formated strings to mp_float */
int mpf_set_str(const char *str, mp_float * c);
/* Conversion of mp_float to string, base 10 only (for now) */
/* This function allocates memory for the string by itself, freeing it is
   left to the caller */
int mpf_get_str(mp_float * a, char **str, int base);

int mpf_get_double(double d, mp_float * a);
int mpf_set_double(mp_float * a, double *d);

int mpf_frexp(mp_float * a, mp_float * b, long *exp);
int mpf_ldexp(mp_float * a, long exp, mp_float * b);

long mpf_digits(mp_float *a);

int mpf_nthroot(mp_float * a, long n, mp_float * b);

int mpf_agm(mp_float * a, mp_float * b, mp_float * c);


int mpf_floor(mp_float *a, mp_float *b);
int mpf_ceil(mp_float *a, mp_float *b);
int mpf_round(mp_float *a, mp_float *b);

int mpf_trig_arg_reduct(mp_float *a, mp_float *b, int *k);

int mpf_sincos(mp_float *a, mp_float *b, int cosine, int tan, int hyper);
int mpf_const_eps(mp_float *a);
int mpf_dump(mp_float * a);

int mpf_sinh(mp_float * a, mp_float * b);
int mpf_cosh(mp_float * a, mp_float * b);
int mpf_tanh(mp_float * a, mp_float * b);

int mpf_kernel_atan(mp_float * a, mp_float * b, int hyper);
int mpf_atan(mp_float * a, mp_float * b);
int mpf_atan2(mp_float * a, mp_float * b, mp_float *c);
int mpf_atanh(mp_float * a, mp_float * b);
int mpf_asinh(mp_float * a, mp_float * b);
int mpf_acosh(mp_float * a, mp_float * b);

// integer trig. functions (for Machin's fomrulas)
#define mp_isone(a) ( (a)->used == 1 && (a)->dp[0] == 1 && (a)->sign == MP_ZPOS )
int mp_acoth_binary_splitting(mp_int * q, mp_int * a, mp_int * b, mp_int * P, mp_int * Q, mp_int * R);
int mp_acot_binary_splitting(mp_int * q, mp_int * a, mp_int * b, mp_int * P, mp_int * Q, mp_int * R);

// Lambert-W aka ProductLog
int mpf_lambertw(mp_float *a, mp_float *b, int branch);



/* temporary functions. if this is non-empty you have just entered
   development area and are on your own risk from now on */
int  mpf_inv_old(mp_float *a, mp_float *b);
int mpf_ln_agm(mp_float * a, mp_float * b);
int mpf_exp_new(mp_float * a, mp_float * b);
int  mpf_cos_old(mp_float *a, mp_float *b);


#endif
