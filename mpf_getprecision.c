#include <tomfloat.h>
/* get the global precision in bits  */
long mpf_getprecision()
{
    return mpf_global_radix;
}