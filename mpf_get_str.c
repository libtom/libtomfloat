#include <tomfloat.h>

static long get_mp_digits(mp_int * a, int base)
{
    double log210;
    long bits = mp_count_bits(a);
    switch (base) {
    case 2:
        log210 = 1.0;
        break;
    case 8:
        log210 = 3.0;
        break;
    case 10:
        log210 = 3.321928094887362347870319429489390175865;
        break;
    case 16:
        log210 = 4.0;
        break;
    default:
        log210 = log(bits) / log(2.0);
        break;
    }
    return (long) (round(bits / log210));
}


/* Convert a mp_float into a string. Scientific notation and base 10 only
  (for now) 
  The string gets allocated by the function, freeing is left to the caller.
  This is not the fastest way to do it
*/
int mpf_get_str(mp_float * a, char **str, int base)
{
  mp_int ret, ten;
  mp_float ret2, ten2, anew, log10, loga;
  mp_digit grs;
  int err, digits, tmpdigits;
  int decprec, decexpo;
  int sign, signbit;
  int rlen;
  char *tmp, *s;
  size_t len, offset;
  long eps, move;
  double d;
  // rounding bit
  mp_digit c;

  if (base != 2 && base != 10 && base != 16) {
    fprintf(stderr, "Base must be one of 2, 10, or 16\n");
    return MP_VAL;
  }
  if (base != 10) {
    fprintf(stderr, "Base is (temporarily) restricted to 10\n");
    return MP_VAL;
  }
  // number of decimal digits to print
  decprec = mpf_getdecimalprecision(a->radix);
  // some guard digits
  eps = (a->radix) + MP_DIGIT_BIT;

  if (mpf_iszero(a)) {
    goto print_zero;
  }


  if ((err = mpf_init_copy(a, &anew)) != MP_OKAY) {
    mp_clear(&ret);
    return err;
  }
  if ((err = mpf_normalize_to(&anew, eps)) != MP_OKAY) {
    mp_clear(&ret);
    mpf_clear(&anew);
    return err;
  }
  if ((err = mpf_init_multi(anew.radix, &log10, &loga, NULL)) != MP_OKAY) {
    mp_clear(&ret);
    mpf_clear(&anew);
    return err;
  }


  if ((&a->mantissa)->sign == MP_NEG) {
    sign = MP_NEG;
    anew.mantissa.sign = MP_ZPOS;
  } else {
    sign = MP_ZPOS;
  }


  // it must be an integer
  if (a->exp >= 0) {
    if ((err = mp_init_copy(&ret, &a->mantissa)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    ret.sign = MP_ZPOS;

    // TODO: implement the other two bases
    if ((err = mp_mul_2d(&ret, (int) a->exp, &ret)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }

    if( (err = mpf_const_le10(&log10) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mpf_ln(&anew, &loga) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mpf_div(&loga, &log10, &loga) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mpf_floor(&loga, &loga) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mpf_set_double(&loga, &d) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    // exponent is always in base ten
    digits = (int) (d);
    decexpo = digits;

    if ((err = mpf_init(&ret2, eps)) != MP_OKAY) {
      mp_clear(&ret);
      return err;
    }

    if ((err = mpf_from_mp_int(&ret, &ret2)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear(&ret2);
      return err;
    }

    if ((err = mpf_init(&ten2, eps)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear(&ret2);
      return err;
    }
    if ((err = mpf_const_d(&ten2, base)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear_multi(&ret2, &ten2, NULL);
      return err;
    }
    if ((err = mpf_pow_d(&ten2, digits - decprec, &ten2)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear_multi(&ret2, &ten2, NULL);
      return err;
    }
    if ((err = mpf_div(&ret2, &ten2, &ret2)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear_multi(&ret2, &ten2, NULL);
      return err;
    }
    if ((err = mpf_normalize_to(&ret2, a->radix)) != MP_OKAY) {
      mp_clear(&ret);
      mpf_clear(&anew);
      return err;
    }

    if (ret2.exp < 0) {
      c = ret2.mantissa.dp[abs(ret2.exp) /
			   MP_DIGIT_BIT] & (1U << (abs(ret2.exp) % DIGIT_BIT));
      if ((err =
	   mp_div_2d(&ret2.mantissa, abs(ret2.exp), &ret2.mantissa,
		     NULL)) != MP_OKAY) {
	mp_clear(&ret);
	mpf_clear_multi(&ret2, &ten2, NULL);
	return err;
      }
      // TODO round according to fenv
      if (c != 0) {
        rlen = mp_count_bits(&ret2.mantissa);
	if ((err = mp_add_d(&ret2.mantissa, 1, &ret2.mantissa)) != MP_OKAY) {
	  mp_clear(&ret);
	  mpf_clear_multi(&ret2, &ten2, NULL);
	  return err;
	}
        // rounding may have added one bit
        if(rlen < mp_count_bits(&ret)){
          if( (err = mp_div_2(&ret,&ret) ) != MP_OKAY){
            mp_clear(&ret);
	    mpf_clear_multi(&ret2, &ten2, NULL);
            return err;
          }
        }
      }
    } else {
      if ((err =
	   mp_mul_2d(&ret2.mantissa, abs(ret2.exp),
		     &ret2.mantissa)) != MP_OKAY) {
	mp_clear(&ret);
	mpf_clear_multi(&ret2, &ten2, NULL);
	return err;
      }
    }
    mp_exch(&ret2.mantissa, &ret);
    mpf_clear_multi(&ret2, &ten2, NULL);
    if (mp_iszero(&ret)) {
      mp_clear(&ret);
      goto print_zero;
    }

    // We have to allocate memory for the resulting string
    // but do not know how much. Finding it out involves a
    // lot of guesswork, hence the large amount of angst-allowance
    tmpdigits = get_mp_digits(&ret, base) + 5;
    // digits in ret + sign + expo.-mark + expo + angst-allowance = 50
    *str = malloc((tmpdigits + 50) * sizeof(char));
    if (*str == NULL) {
      fprintf(stderr, "malloc failed to allocate %lu bytes\n",
	      (tmpdigits + 50) * sizeof(char));
      mp_clear_multi(&ret, &ten, NULL);
      return err;
    }
    if (sign == MP_NEG) {
      ret.sign = MP_NEG;
    }
    if ((err = mp_toradix(&ret, *str, base)) != MP_OKAY) {
      mp_clear_multi(&ret, &ten, NULL);
      free(str);
      return err;
    }
    goto print_number;

  } else {
    // it is not necessarily an integer
    if ((err = mp_init_copy(&ret, &anew.mantissa)) != MP_OKAY) {
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    ret.sign = MP_ZPOS;
    decexpo = mpf_getdecimalexponent(a->exp);
    tmpdigits = get_mp_digits(&ret, base) + 5;

    // TODO: implement the other two bases
    if( (err = mpf_const_le10(&log10) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mpf_ln(&anew, &loga) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mpf_div(&loga, &log10, &loga) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    signbit = mpf_signbit(&loga);
    if( (err = mpf_set_double(&loga, &d) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    d = fabs(d);
    
    digits = (int) (d);
    decexpo = digits;

    if( (err = mp_init_set_int(&ten,10L) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mp_expt_d(&ten,digits + decprec,&ten) ) != MP_OKAY){
      mp_clear_multi(&ret,&ten,NULL);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    if( (err = mp_mul(&ret,&ten,&ret) ) != MP_OKAY){
      mp_clear_multi(&ret,&ten,NULL);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }

    mp_clear(&ten);

    // TODO: round according to fenv!
    move = -((&anew)->exp);
    // We need at least four bits for rounding. The simplest thing is
    // to shift "ret" for "move - 4" first to get the bits.
    if( (err = mp_div_2d(&ret,move - 4,&ret,NULL) ) != MP_OKAY){
      mp_clear(&ret);
      mpf_clear_multi(&anew,&loga,&log10,NULL);
      return err;
    }
    // get LSL (last significant limb) and the last 4 bits of it
    grs = ret.dp[0] & 0xF;
    // round nearest to even
    /*
                     decimal
        x0xx down       -       
        1100 up        12
        0100 down       4
        x101 up         5
        x110 up         6
        x111 up         7
    */
    //    0100  ||      x0xx
    if(grs == 4 || ((grs>>2)&1) == 0){
      // just truncate
      if( (err = mp_div_2d(&ret,4,&ret,NULL) ) != MP_OKAY){
        mp_clear(&ret);
        mpf_clear_multi(&anew,&loga,&log10,NULL);
        return err;
      }
    } else {
      // adding might carry, save length
      rlen = mp_count_bits(&ret);
      if( (err = mp_add_d(&ret, 4, &ret) ) != MP_OKAY){
        mp_clear(&ret);
        mpf_clear_multi(&anew,&loga,&log10,NULL);
        return err;
      }
      if(rlen < mp_count_bits(&ret)){
        if( (err = mp_div_2d(&ret,3,&ret,NULL) ) != MP_OKAY){
          mp_clear(&ret);
          mpf_clear_multi(&anew,&loga,&log10,NULL);
          return err;
        }
      } else {
        if( (err = mp_div_2d(&ret,4,&ret,NULL) ) != MP_OKAY){
          mp_clear(&ret);
          mpf_clear_multi(&anew,&loga,&log10,NULL);
          return err;
        }
      }
    }

    // digits in ret + sign + expo.-mark + expo + angst-allowance = 50
    *str = malloc((tmpdigits + 50) * sizeof(char));
    if (*str == NULL) {
      fprintf(stderr, "malloc failed to allocate %lu bytes\n",
	      (tmpdigits + 50) * sizeof(char));
      mp_clear_multi(&ret, &ten, NULL);
      return err;
    }
    if (sign == MP_NEG) {
      ret.sign = MP_NEG;
    }
    if ((err = mp_toradix(&ret, *str, base)) != MP_OKAY) {
      mp_clear_multi(&ret, &ten, NULL);
      free(str);
      return err;
    }

    // round away from zero
    if (signbit == MP_NEG) {
      digits = -ceil(d);
    } else {
      digits = floor(d);
    }
    decexpo = digits;
    if (base == 10) {
      if (sign == MP_NEG) {
	*((*str) + decprec + 1) = '\0';
      } else {
	*((*str) + decprec) = '\0';
      }
    }
    goto print_number;
  }

print_zero:
  if (base == 10) {
    *str = malloc(6 * sizeof(char));
    if (NULL == *str) {
      fprintf(stderr, "malloc failed to allocate %lu bytes\n",
	      6 * sizeof(char));
    }
    *str = strcpy(*str, "0.0E0");
  } else {
    *str = malloc(8 * sizeof(char));
    if (NULL == *str) {
      fprintf(stderr, "malloc failed to allocate %lu bytes\n",
	      8 * sizeof(char));
    }
    if (base == 16) {
      *str = strcpy(*str, "0x0.0P0");
    } else {
      *str = strcpy(*str, "0b0.0P0");
    }
  }
  return MP_OKAY;

print_number:
  // assuming 8 bit char < 4 decimal digits (plus sign, EOS and a.-a.)
  tmp = malloc((sizeof(int) * 4 + 3) * sizeof(char));
  if (tmp == NULL) {
    fprintf(stderr, "malloc failed to allocate %lu bytes\n",
	    (sizeof(long) * 4 + 3) * sizeof(char));
    mp_clear_multi(&ret, &ten, NULL);
    free(str);
    return err;
  }
  // TODO: check return of snprintf(3)
  snprintf(tmp, (sizeof(int) * 4 + 3) * sizeof(char), "%d", decexpo);
  /*
   * The digits are allocated, the exponent is allocated, leaves the
   * sign, the decimal mark (we do scientific notation only), the
   * exponent mark and the EOS; makes 4. So lets take 5 because I
   * like the rythm. And don't forget to add 3 more for a prefix.
   */
  if (strlen(*str) + strlen(tmp) + 8 >= (tmpdigits + 50) * sizeof(char)) {
    // TODO: use realloc() ?
    // Bail out. For now.
    fprintf(stderr, "Not enough memory allocated for output. Please report\n");
    mp_clear_multi(&ret, &ten, NULL);
    free(str);
    free(tmp);
    return err;
  }
  // a bit of pointer juggling, sorry.
  // string length plus EOS
  len = strlen(*str) + 1;
  // adjust for a sign if necessary
  if (ret.sign == MP_NEG) {
    offset = 2;
  } else {
    offset = 1;
  }
  // move the whole thing by one byte
  s = memmove((*str) + 1, *str, len);
  // set first digit
  *((*str) + 1) = *(s + 1);
  // set decimal point
  *((*str) + offset) = '.';

  // put all together: prefix, exponent mark, and exponent
  if (base == 2 || base == 16) {
    len = strlen(*str) + 1;
    // move the whole thing by two bytes
    s = memmove((*str) + 2, *str, len);
    if (ret.sign == MP_NEG) {
      // set sign
      *((*str)) = *s;
      // set prefix
      *((*str) + 1) = '0';
      *((*str) + 2) = (base == 2) ? 'b' : 'x';
    } else {
      *(*str) = '0';
      *((*str) + 1) = (base == 2) ? 'b' : 'x';
    }
    *str = strcat(*str, "p");
    *str = strcat(*str, tmp);
  } else {
    *str = strcat(*str, "e");
    *str = strcat(*str, tmp);
  }
  mp_clear(&ret);
  free(tmp);
  return MP_OKAY;
}

