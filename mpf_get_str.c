#include <tomfloat.h>


/* Convert a mp_float into a string. Scientific notation and base 10 only
  (for now) 
  The string gets allocated by the function, freeing is left to the caller.
  This is not the fastest way to do it
*/
int mpf_get_str(mp_float * a, char **str, int base)
{
    mp_int ret, ten;
    mp_float ret2, ten2, anew, log10,loga;
    int err, digits, tmpdigits;
    int decprec, decexpo;
    int sign, signbit;
    char *tmp, *s;
    size_t len, offset;
    long eps;
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
    decprec = mpf_getdecimalprecision(a->radix);
    // some guard digits
    eps = (a->radix) + MP_DIGIT_BIT;

    if (mpf_iszero(a)) {
	goto print_zero;
    }

    if ((err = mp_init_copy(&ret, &a->mantissa)) != MP_OKAY) {
	return err;
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
    if ((err = mpf_init_multi(anew.radix,&log10,&loga, NULL)) != MP_OKAY) {
	mp_clear(&ret);
	mpf_clear(&anew);
	return err;
    }


    if (ret.sign == MP_NEG) {
	sign = MP_NEG;
	ret.sign = MP_ZPOS;
	anew.mantissa.sign = MP_ZPOS;
    } else {
	sign = MP_ZPOS;
    }


    // it must be an integer
    if (a->exp >= 0) {
	// TODO: implement the other two bases
	if ((err = mp_mul_2d(&ret, (int) a->exp, &ret)) != MP_OKAY) {
	    mp_clear(&ret);
	    return err;
	}

        mp_radix_size(&ret,10,&digits);
	decexpo = digits - 2;
// TODO: log is fast enough now, use it
/*
        mpf_const_le10(&log10);
        mpf_ln(&anew,&loga);
        mpf_div(&loga,&log10,&loga);
        mpf_floor(&loga,&loga);
        mpf_set_double(&loga, &d);
	// exponent is always in base ten
        digits = (int)(d);
	decexpo = digits;
*/


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
				 MP_DIGIT_BIT] & (1U << (abs(ret2.exp) %
							 DIGIT_BIT));
	    if ((err =
		 mp_div_2d(&ret2.mantissa, abs(ret2.exp), &ret2.mantissa,
			   NULL)) != MP_OKAY) {
		mp_clear(&ret);
		mpf_clear_multi(&ret2, &ten2, NULL);
		return err;
	    }
            // TODO: rounding may have added a bit
	    if (c != 0) {
//printf("bits = %d\n",mp_count_bits(&ret));
		if ((err = mp_add_d(&ret2.mantissa, 1, &ret2.mantissa)) != MP_OKAY) {
		    mp_clear(&ret);
		    mpf_clear_multi(&ret2, &ten2, NULL);
		    return err;
		}
//printf("bits = %d\n",mp_count_bits(&ret));
	    }
	} else {
	    if ((err =
		 mp_mul_2d(&ret2.mantissa, abs(ret2.exp), &ret2.mantissa)) != MP_OKAY) {
		mp_clear(&ret);
		mpf_clear_multi(&ret2, &ten2, NULL);
		return err;
	    }
	}
        mp_exch(&ret2.mantissa,&ret);
	mpf_clear_multi(&ret2, &ten2, NULL);
	if (mp_iszero(&ret)) {
	    mp_clear(&ret);
	    goto print_zero;
	}
        // TODO: the code below is used twice with very small chances
	// We have to allocate memory for the resulting string
	// but do not know how much. Finding it out involves a
	// lot of guesswork, hence the large amount of angst-allowance
	tmpdigits = mp_digits(&ret, base) + 5;
	// digits in ret + sign + expo.-mark + expo + angst-allowance = 50
	*str = malloc((tmpdigits + 50) * sizeof(char));
	if (*str == NULL) {
	    fprintf(stderr, "malloc failed to allocate %d bytes\n",
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

	// assuming 8 bit char < 4 decimal digits (plus sign, EOS and a.-a.)
	tmp = malloc((sizeof(int) * 4 + 3) * sizeof(char));
	if (tmp == NULL) {
	    fprintf(stderr, "malloc failed to allocate %d bytes\n",
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
	    fprintf(stderr,
		    "Not enough memory allocated for output. Please report\n");
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

    } else {
	// it is not necessarily an integer
	decexpo = mpf_getdecimalexponent(a->exp);
	if ((err = mpf_init(&ret2, eps)) != MP_OKAY) {
	    mp_clear(&ret);
	    return err;
	}
	if ((err = mpf_init(&ten2, eps)) != MP_OKAY) {
	    mp_clear(&ret);
	    mpf_clear(&ret2);
	    return err;
	}
	if ((err = mpf_const_d(&ten2, 10)) != MP_OKAY) {
	    mp_clear(&ret);
	    mpf_clear_multi(&ret2, &ten2, NULL);
	    return err;
	}
	if ((err = mpf_pow_d(&ten2, decexpo, &ten2)) != MP_OKAY) {
	    mp_clear(&ret);
	    mpf_clear_multi(&ret2, &ten2, NULL);
	    return err;
	}
	if ((err = mpf_mul(&anew, &ten2, &ret2)) != MP_OKAY) {
	    mp_clear(&ret);
	    mpf_clear_multi(&ret2, &ten2, NULL);
	    return err;
	}
	if ((err = mpf_normalize_to(&ret2, a->radix)) != MP_OKAY) {
	    mp_clear(&ret);
	    mpf_clear(&anew);
	    return err;
	}
	// round
	if (ret2.exp < 0) {
	    c = ret2.mantissa.dp[(abs(ret2.exp)) /
				 MP_DIGIT_BIT] & (1U << ((abs(ret2.exp)) %
							 DIGIT_BIT));
	    if ((err =
		 mp_div_2d(&ret2.mantissa, abs(ret2.exp), &ret,
			   NULL)) != MP_OKAY) {
		mp_clear(&ret);
		mpf_clear_multi(&ret2, &ten2, NULL);
		return err;
	    }
	    if (c != 0) {
		if ((err = mp_add_d(&ret, 1, &ret)) != MP_OKAY) {
		    mp_clear(&ret);
		    mpf_clear_multi(&ret2, &ten2, NULL);
		    return err;
		}
	    }
	} else {
	    if ((err =
		 mp_mul_2d(&ret2.mantissa, abs(ret2.exp), &ret)) != MP_OKAY) {
		mp_clear(&ret);
		mpf_clear_multi(&ret2, &ten2, NULL);
		return err;
	    }
	    if (mp_count_bits(&ret) > mp_count_bits(&ret2.mantissa)) {
		ret2.exp--;
	    }
	}
	mpf_clear_multi(&ret2, &ten2, NULL);

	if (mp_iszero(&ret)) {
	    mp_clear(&ret);
	    goto print_zero;
	}

	tmpdigits = mp_digits(&ret, base) + 3;

	// digits in ret + sign + expo.-mark + expo + angst-allowance = 50
	*str = malloc((tmpdigits + 50) * sizeof(char));
	if (*str == NULL) {
	    fprintf(stderr, "malloc failed to allocate %d bytes\n",
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
        // TODO: implement the other two bases
        mpf_const_le10(&log10);
        mpf_ln(&anew,&loga);
        mpf_div(&loga,&log10,&loga);
        signbit = mpf_signbit(&loga);
        mpf_set_double(&loga, &d);
        d = fabs(d);
	// round away from zero
        if(signbit == MP_NEG){
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
	tmp = malloc((sizeof(int) * 4 + 3) * sizeof(char));
	if (tmp == NULL) {
	    fprintf(stderr, "malloc failed to allocate %d bytes\n",
		    (sizeof(long) * 4 + 3) * sizeof(char));
	    mp_clear_multi(&ret, &ten, NULL);
	    free(str);
	    return err;
	}
	snprintf(tmp, (sizeof(int) * 4 + 3) * sizeof(char), "%d", decexpo);
	if (strlen(*str) + strlen(tmp) + 5 >= (tmpdigits + 50) * sizeof(char)) {
	    fprintf(stderr,
		    "Not enough memory allocated for output. Please report\n");
	    mp_clear_multi(&ret, &ten, NULL);
	    free(str);
	    free(tmp);
	    return err;
	}
	len = strlen(*str) + 1;
	if (ret.sign == MP_NEG) {
	    offset = 2;
	} else {
	    offset = 1;
	}
	s = memmove((*str) + 1, *str, len);
	*((*str) + 1) = *(s + 1);
	*((*str) + offset) = '.';

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
	return MP_OKAY;
    }
  print_zero:
    if (base == 10) {
	*str = malloc(6 * sizeof(char));
	if (NULL == *str) {
	    fprintf(stderr, "malloc failed to allocate %d bytes\n",
		    6 * sizeof(char));
	}
	*str = strcpy(*str, "0.0E0");
    } else {
	*str = malloc(8 * sizeof(char));
	if (NULL == *str) {
	    fprintf(stderr, "malloc failed to allocate %d bytes\n",
		    8 * sizeof(char));
	}
	if (base == 16) {
	    *str = strcpy(*str, "0x0.0P0");
	} else {
	    *str = strcpy(*str, "0b0.0P0");
	}
    }
    return MP_OKAY;
}



