#include <stdlib.h>
#include <stdint.h>

#include "lfa.h"

#ifdef _WIN32

#define WINVER 0x0500
#include <windows.h>

int get_num_cpus (void) {
  SYSTEM_INFO si;

  GetSystemInfo(&si);

  return(si.dwNumberOfProcessors);
}

uint64_t get_total_mem (void) {
  MEMORYSTATUSEX stat;

  stat.dwLength = sizeof(stat);

  if(GlobalMemoryStatusEx(&stat))
    return((uint64_t) stat.ullTotalPhys);
  else
    return(0);
}

#else  /* assume POSIX */

#include <sys/types.h>
#include <unistd.h>

#if defined(_SC_NPROCESSORS_ONLN)

int get_num_cpus (void) {
  return(sysconf(_SC_NPROCESSORS_ONLN));
}

#elif defined(_SC_NPROC_ONLN)

int get_num_cpus (void) {
  return(sysconf(_SC_NPROC_ONLN));
}

#elif defined(__FreeBSD__) || defined(__NetBSD__) || \
      defined(__OpenBSD__) || defined(__APPLE__)

#include <sys/param.h>
#include <sys/sysctl.h>

int get_num_cpus (void) {
  int mib[4] = {
    CTL_HW,
#if defined(HW_AVAILCPU)
    HW_AVAILCPU
#else
    HW_NCPU
#endif
  };
  int rv;
  size_t len = sizeof(rv);

  if(sysctl(mib, 2, &rv, &len, NULL, 0))
    return(-1);

  return(rv);
}

#else  /* don't know */

int get_num_cpus (void) {
  return(-1);
}

#endif

#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)

uint64_t get_total_mem (void) {
  long npages, pagesize;

  npages = sysconf(_SC_PHYS_PAGES);
  pagesize = sysconf(_SC_PAGE_SIZE);

  if(npages == -1 || pagesize == -1)
    return(0);
  else
    return(((uint64_t) npages) * ((uint64_t) pagesize));
}

#elif defined(__FreeBSD__) || defined(__NetBSD__) || \
      defined(__OpenBSD__) || defined(__APPLE__)

#include <sys/param.h>
#include <sys/sysctl.h>

uint64_t get_total_mem (void) {
  int mib[4] = { CTL_HW, HW_MEMSIZE };
  uint64_t rv;
  size_t len = sizeof(rv);

  if(sysctl(mib, 2, &rv, &len, NULL, 0))
    return(-1);

  return(rv);
}

#else  /* don't know */

uint64_t get_total_mem (void) {
  return(0);
}

#endif

#endif  /* _WIN32 */

