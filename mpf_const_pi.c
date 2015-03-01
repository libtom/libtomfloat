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

static mp_float mpf_pi;
// This is automatically initialized to 0 according to C99
// TODO: check the other standard versions
static long mpf_pi_precision;

/* Pi by the AGM (Brent-Salamin) */
int mpf_const_pi(mp_float * a)
{
    int err;
    mp_float aa, b, d, d2, s, t, p, two, twoinv;
    long eps, oldeps, extra, k, loops,r;

    oldeps = a->radix;
    err = MP_OKAY;

    // Sometimes memory is short, even today
    if(mpf_pi_precision > 0 && a == NULL){
      mpf_clear(&mpf_pi);
      mpf_pi_precision = 0;
      return err;
    }

    // five percent plus 3 bit angst-allowance
    // TODO: compute correct value
    extra = (oldeps / 100) * 5 + 3;
    if (mpf_pi_precision >= oldeps) {
	if ((err = mpf_copy(&mpf_pi, a)) != MP_OKAY) {
	    return err;
	}
	return mpf_normalize_to(a, oldeps);
    } else {
	eps = oldeps + extra;
	if (mpf_pi_precision == 0) {
	    mpf_init(&mpf_pi, eps);
	}

	mpf_init_multi(eps, &aa, &b, &d, &d2, &s, &t, &p, &two, &twoinv, NULL);

        // this algorithm is quadratic, produces twice the amount of digits
        // every round, so the maximum number of rounds is log_2(radix)
        loops = 0;
        // if eps is negative we have a really large problem but somewhere else
        r = eps;
        while (r >>= 1){
           loops++;
        }
        // ceil(log_2(radix)) + angst-allowance
        loops += 3;

        /*
            Initialize:
            a = 1, b = 1/sqrt(2), t = 1/2, k = 1
         */
	// aa = 1
	mpf_const_d(&aa, 1);
	// two = 2
	mpf_const_d(&two, 2);
	// twoinv = 1/two
	mpf_inv(&two, &twoinv);
	// b = sqrt(two)
	mpf_sqrt(&two, &b);
	// b = 1/b
	mpf_inv(&b, &b);
	// t = twoinv
	mpf_copy(&twoinv, &t);
	k = 1;

	do {
            /* s = 1/2 * (a + b) */
	    // s = aa + b
	    mpf_add(&aa, &b, &s);
	    // s *= twoinv
	    mpf_mul(&s, &twoinv, &s);
	    /* d = aa - s */
	    mpf_sub(&aa, &s, &d);
	    // d = d^2
	    mpf_sqr(&d, &d);
	    // aa = s
	    mpf_copy(&s, &aa);
	    // s = s^2
	    mpf_sqr(&s, &s);
            /* t = t- 2^k * d */
	    // d2 = d * 2^k
	    mpf_copy(&d, &d2);
	    d2.exp += k;
	    // t = t - d2
	    mpf_sub(&t, &d2, &t);
            /* b = sqrt(s - d) */
	    // b = s - d
	    mpf_sub(&s, &d, &b);
	    // b = sqrt(b);
	    mpf_sqrt(&b, &b);
	    k++;
            if(k == loops){
              err = MP_VAL;
#ifdef DEBUG
              fprintf(stderr,"Pi by AGM did not converge in %ld rounds\n",k);
#endif
              goto _ERR;
            }
	} while (mpf_cmp(&aa, &b) != MP_EQ);
        /*
                  (a + b)^2
            pi = -----------
                    2t
        */
	// t = t * 2
	t.exp += 1;
	// t = 1/t
	mpf_inv(&t, &t);
	// p = aa + b but aa = b, so p = aa * 2
	aa.exp += 1;
	// p = p^2
	mpf_sqr(&aa, &p);
	// p = p * t
	mpf_mul(&p, &t, &mpf_pi);

	mpf_copy(&mpf_pi, a);
        mpf_pi_precision = oldeps;
	return mpf_normalize_to(a, oldeps);
    }
  _ERR:
    // Make it clear that it was a failure and that the value of mpf_pi cannot
    // be trusted anymore by setting mpf_pi_precision to a value that
    // cannot happen naturally.
    mpf_pi_precision = -1;
    mpf_clear_multi(&aa, &b, &d, &d2, &s, &t, &p, &two, &twoinv, NULL);
    return err;
}

