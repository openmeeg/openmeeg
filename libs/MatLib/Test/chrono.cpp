/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#if WIN32
#include <time.h>
#else
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#endif

static double ref;
void timer_start()
{
#if WIN32
    ref=double(clock())/CLOCKS_PER_SEC;
#else
    struct tms b;
    times(&b);
    ref=double(b.tms_utime)/sysconf(_SC_CLK_TCK);
#endif
}

double timer_lap()
{
#if WIN32
    double lap=double(clock())/CLOCKS_PER_SEC-ref;
#else
    struct tms b;
    times(&b);
    double lap=double(b.tms_utime)/sysconf(_SC_CLK_TCK)-ref;
#endif
    return lap;
}
