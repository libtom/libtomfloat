#include <tomfloat.h>
/* approximate(!) the number of base 10 digits of the mp_float.
   either the number of decimal digits in the integer part if it exists
   with a positive value or the number of decimal digits in the fractional part
   with a negative value. An approximation of the decimal exponent used with the
   scientific format.
 */
long mpf_digits(mp_float * a)
{
    long ret;
    double log102 = 3.321928094887362347870319429489390175864831393;

    ret = ceil((a->exp + a->radix) / log102);
    if (ret < 0) {
	ret++;
    }
    return ret;
}

