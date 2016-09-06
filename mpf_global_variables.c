#include <tomfloat.h>
/* Global radix value */
long mpf_global_radix = 4 * MP_DIGIT_BIT;
/* Extra bits for rounding */
long mpf_guard_bits   = 3;


/* global error */
int mpf_errno = 0;

/* cut-off values in log (in bits, if not indicated otherwise)*/
int MPF_LOG_AGM_1_CUTOFF = 100;
// the n in nth-root for argument reduction in the log-series algo.
long MPF_LOG_AGM_REDUX_CUTOFF = 4096; 

