#include <tomfloat.h>
/* Get the global precision in decimal */
long mpf_getdecimalprecision(long radix)
{
    double log210 = 3.321928094887362347870319429489390175865;
    long bits;
    if(radix <= 0){
        bits = mpf_global_radix;
    } else {
        bits = radix;
    }
    return (floor(bits / log210) + 1);
}