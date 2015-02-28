#include <tomfloat.h>
/* Get the global precision in decimal */
long mpf_getdecimalprecision()
{
    double log210 = 3.321928094887362347870319429489390175865;
    long bits = mpf_global_radix;
    return (floor(bits / log210) + 1);
}