#include <tomfloat.h>
/* Get the argument converted to decimal */
long mpf_getdecimalexponent(long exp)
{
    double log210 = 3.321928094887362347870319429489390175865;
    long bits = abs(exp);
    return (lround(bits / log210) + 1);
}
