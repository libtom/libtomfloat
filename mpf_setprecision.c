#include <tomfloat.h>
/* Set the global precision */ 
#ifdef MPF_USE_THREADS
#   include <pthread.h>
// don't forget to do a pthread_mutex_init(&radix_mutex, NULL);
pthread_mutex_t radix_mutex;
void mpf_setprecision(long r)
{
    pthread_mutex_lock(&radix_mutex);
    mpf_global_radix = r;
    pthread_mutex_unlock(&radix_mutex);
}
#else
void mpf_setprecision(long r)
{
    mpf_global_radix = r;
}
#endif