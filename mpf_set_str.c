#include <tomfloat.h>

// unistd.h maight not exist in non posix systems
#if !defined(_WIN32)
#   include <unistd.h>
#endif
// for str(n)casecmp if on a unix-type OS
#include <strings.h>
// for, who might have guessed, errno
#include <errno.h>

/* Convert an IEEE-754 formated string (base 2, 20, and 16) into a mp_float */

/* str(n)casecmp is not in STD-C but in POSIX.1-2001 (and 4.4BSD)*/
/* Used for detecting the strings "Infinity",  "Inf", and "NaN" */
#if (!(defined _POSIX_VERSION ) || (_POSIX_VERSION < 200112L))
static int strncasecmp(const char *s1, const char *s2, int n)
{
    if (n == 0) {
	return 0;
    }
    while (tolower(*s1) == tolower(*s2) && n-- != 0) {
	if (n == 0 || *s1 == '\0' || *s2 == '\0') {
	    break;
	}
	s1++;
	s2++;
    }
    return tolower(*(unsigned char *) s1) - tolower(*(unsigned char *) s2);
}

static int strcasecmp(const char *s1, const char *s2)
{
    while (tolower(*s1) == tolower(*s2)) {
	if (*s1 == '\0' || *s2 == '\0') {
	    break;
	}
	s1++;
	s2++;
    }
    return tolower(*(unsigned char *) s1) - tolower(*(unsigned char *) s2);
}
#endif
/* 
 Case ignorant version of strchr(3)
 Returns a pointer to the findings and sets len to the position of it
 Used to detect an exponent mark and to put a pointer at that position
 to be able to use strtol(3) directly;
*/
static char *strcasechar(const char *s, int c, size_t * len)
{
    size_t l = 0;
    while (tolower(*s) != c) {
	if (*s == '\0') {
	    return NULL;
	}
	s++;
	l++;
    }
    *len = l;
    return (char *) s;
}

int mpf_set_str(const char *str, mp_float * c)
{
    mp_float tmp, ret;
    int sign = MP_ZPOS;
    int base = 10;
    char *exponent, *endptr;
    size_t *slen = 0, k = 0;
    long expo = 0, decimalpoint = -1, fdigs = 0, eps;
    int err;


    if (NULL == str) {
	fprintf(stderr, "Input is NULL in parsenumber\n");
	if ((err = mpf_const_nan(c)) != MP_OKAY) {
	    return err;
	}
	return MP_VAL;
    }
    if (strlen(str) == 0) {
	fprintf(stderr, "Input has length zero in parsenumber\n");
	if ((err = mpf_const_nan(c)) != MP_OKAY) {
	    return err;
	}
	return MP_VAL;
    }
    if (*str == '\0') {
	fprintf(stderr, "Input has no characters in parsenumber\n");
	if ((err = mpf_const_nan(c)) != MP_OKAY) {
	    return err;
	}
	return MP_VAL;
    }

    if (*str == '-' || *str == '+') {
	if (*str == '-') {
	    sign = MP_NEG;
	}
	str++;
    }
    if (!strncasecmp(str, "inf", 3)) {
	if ((err = mpf_const_inf(c, c->mantissa.sign)) != MP_OKAY) {
	    return err;
	}
	return MP_OKAY;
    }

    if (!strncasecmp(str, "nan", 3)) {
	if ((err = mpf_const_nan(c)) != MP_OKAY) {
	    return err;
	}
	// return an error here?
	return MP_OKAY;
    }

    if (*str == '0') {
	if (*(str + 1) == 'x' || *(str + 1) == 'b') {
	    if (*(str + 1) == 'x') {
		base = 16;
	    }
	    if (*(str + 1) == 'b') {
		base = 2;
	    }
	    str += 2;
	}
    }
    // strip leading zeros
    while (*str == '0') {
	str++;
    }
    if (*str == '\0') {
	// Zero is the special number 0e1
	c->exp = 1;
	return MP_OKAY;
    }

    /* 
     * Grab the exponent first.
     * Allowing bases 2, 10 and 16 only, the corresponding exponent marks
     * are the case ignorant letters "E" for base 10 and and "P" for the
     * two other bases.
     */
    // simply count length instead of illegible pointer juggling
    slen = malloc(sizeof(size_t));
    if (NULL == slen) {
	fprintf(stderr, "malloc failed to allocate a single size_t\n");
    }
    *slen = strlen(str);
    exponent = NULL;
    if (base == 10) {
	exponent = strcasechar(str, 'e', slen);
    }
    if (base == 2 || base == 16) {
	exponent = strcasechar(str, 'p', slen);
    }
    if (exponent != NULL) {
	// skip the exponent mark
	exponent++;
	// The IEEE-754 exponent is always encoded in base ten
	expo = strtol(exponent, &endptr, 10);
	if ((errno == ERANGE && (expo == LONG_MAX || expo == LONG_MIN))
	    || (errno != 0 && expo == 0)) {
	    fprintf(stderr, "An error occured in parsing the exponent\n");
	    free(slen);
	    if ((err = mpf_const_nan(c)) != MP_OKAY) {
		return err;
	    }
	    return MP_VAL;
	}
	if (endptr == exponent) {
	    fprintf(stderr, "No digits in the exponent\n");
	    free(slen);
	    if ((err = mpf_const_nan(c)) != MP_OKAY) {
		return err;
	    }
	    return MP_VAL;
	}
    }
    // Add some guard bits
    // adding MP_DIGIT_BIT adds on limb at least
    eps = c->radix + MP_DIGIT_BIT;
    if ((err = mpf_init(&ret, eps)) != MP_OKAY) {
	goto free_slen;
    }


    while (*str != '\0' && k++ < *slen) {
	switch (tolower(*str)) {
	case '0':
	case '1':
	case '2':
	    if ((err = mpf_mul_d(&ret, base, &ret)) != MP_OKAY) {
		goto clear_ret;
	    }
	    if ((err = mpf_add_d(&ret, ((int) (*str)) - 48, &ret)) != MP_OKAY) {
		goto clear_ret;
	    }
	    if (decimalpoint > 0) {
		fdigs++;
	    }
	    break;
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
	    if (base == 2) {
		fprintf(stderr, "wrong digit for base %i found\n", base);
		mpf_clear(&ret);
		free(slen);
		if ((err = mpf_const_nan(c)) != MP_OKAY) {
		    return err;
		}
		return MP_VAL;
	    }
	    if ((err = mpf_mul_d(&ret, base, &ret)) != MP_OKAY) {
		goto clear_ret;
	    }
	    if ((err = mpf_add_d(&ret, ((int) (*str)) - 48, &ret)) != MP_OKAY) {
		goto clear_ret;
	    }
	    if (decimalpoint > 0) {
		fdigs++;
	    }
	    break;
	case 'a':
	case 'b':
	case 'c':
	case 'd':
	case 'e':
	case 'f':
	    if (base == 2 || base == 10) {
		fprintf(stderr, "wrong digit for base %i found\n", base);
		free(slen);
		mpf_clear(&ret);
		if ((err = mpf_const_nan(c)) != MP_OKAY) {
		    return err;
		}
		return MP_VAL;
	    }
	    if ((err = mpf_mul_d(&ret, base, &ret)) != MP_OKAY) {
		goto clear_ret;
		return err;
	    }
	    if ((err =
		 mpf_add_d(&ret, ((int) tolower(*str)) - 87,
			   &ret)) != MP_OKAY) {
		goto clear_ret;
	    }
	    if (decimalpoint > 0) {
		fdigs++;
	    }
	    break;
	case '.':
	    // TODO: decimalpoint is only a flag now, change code accordingly
	    decimalpoint = (long) k;
	    break;
	default:
	    fprintf(stderr, "Unknown character %c found\n", *str);
	    free(slen);
	    mpf_clear(&ret);
	    if ((err = mpf_const_nan(c)) != MP_OKAY) {
		return err;
	    }
	    return MP_VAL;
	}
	str++;
    }
    // take every chance to get out early
    if (mpf_iszero(&ret)) {
	mpf_clear(&ret);
	free(slen);
	return MP_OKAY;
    }
    if ((err = mpf_init(&tmp, eps)) != MP_OKAY) {
	goto clear_ret;
    }
    // adjust fractional part
    if (fdigs != 0) {
	if ((err = mpf_const_d(&tmp, base)) != MP_OKAY) {
	    goto clear_tmp;
	}
	if ((err = mpf_pow_d(&tmp, fdigs, &tmp)) != MP_OKAY) {
	    goto clear_tmp;
	}
	if ((err = mpf_div(&ret, &tmp, &ret)) != MP_OKAY) {
	    goto clear_tmp;
	}
    }
    // weave exponent in
    if (exponent != NULL) {
	// the only exponents allowed here are all of base 10
	if ((err = mpf_const_d(&tmp, 10)) != MP_OKAY) {
	    goto clear_tmp;
	}
	// pow_d knows how to handle negative exponents
	if ((err = mpf_pow_d(&tmp, expo, &tmp)) != MP_OKAY) {
	    goto clear_tmp;
	}
	if ((err = mpf_mul(&ret, &tmp, &ret)) != MP_OKAY) {
	    goto clear_tmp;
	}
    }
    mpf_normalize_to(&ret, c->radix);

    if (sign == MP_NEG) {
	ret.mantissa.sign = MP_NEG;
    }

    if ((err = mpf_copy(&ret, c)) != MP_OKAY) {
	goto clear_tmp;
    }

    return MP_OKAY;

  clear_tmp:
    mpf_clear(&tmp);
  clear_ret:
    mpf_clear(&ret);
  free_slen:
    free(slen);
    return err;
}

