#include <tomfloat.h>

static int mpf_print(mp_float * a)
{
    int size, err;
    // one for EOS and one for sign
    mp_fput(&(a->mantissa),10,stdout);
    printf(" * 2^%ld * 1.0\n", a->exp);

    return MP_OKAY;
}

/* round to next integer (towards +infinity)*/
int mpf_round(mp_float *a, mp_float *b){
    int err, cmp;
    mp_float two, ret;

    if (mpf_isnan(a) || mpf_isinf(a)|| mpf_iszero(a)) {
        return mpf_copy(a,b);
    }
    // integer
    if (a->exp > 0) {
        return mpf_copy(a,b);
    }

    mpf_init_multi(a->radix,&two,&ret,NULL);
    mpf_const_d(&two,2);
    // using inv() would cause rounding errors
    two.exp -= 2;

    // fraction
    if (a->exp <= -a->radix) {
        mpf_abs(a,&ret);
        cmp = mpf_cmp(&ret, &two);
        // >= .5 round up
        if (cmp != MP_LT) {
            if(cmp == MP_EQ){
               // -.5 rounds up to zero
               if(a->mantissa.sign == MP_NEG){
                  mpf_const_0(b);
                  goto _ERR;
               }
            }
            if (a->mantissa.sign == MP_NEG) {
                mpf_const_d(b,-1);
                goto _ERR;
            }
            mpf_const_d(b,1);
            goto _ERR;
        } else {
             mpf_const_0(b);
             goto _ERR;
        }
    }
    // mixed/integer
    // just add .5 for now, other rounding methods might come later
    mpf_add(a,&two,b);
    mpf_floor(b,b);
    err = MP_OKAY;

_ERR:
 mpf_clear_multi(&two,&ret,NULL);
   return err;
}

