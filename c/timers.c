/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifdef _WIN32
#include <windows.h>
double cpuClock() {
    	LARGE_INTEGER fi,ti;
    	double fd,td;
    	QueryPerformanceFrequency( &fi);
    	fd=(double) fi.QuadPart;
    	QueryPerformanceCounter(&ti);
	td=(double) ti.QuadPart;
	return td/fd;
}
void getWallTime(double* wcTime) {
    LARGE_INTEGER ti,fi;
    double fd,td;
    FILETIME ft;
    QueryPerformanceFrequency( &fi);
    fd=(double) fi.QuadPart;
    GetSystemTimeAsFileTime(&ft);
    ti.LowPart=ft.dwLowDateTime;
    ti.HighPart=ft.dwHighDateTime;
    td=(double) ti.QuadPart;
    return td/fd;
}

#else
#include <sys/time.h>
#include <time.h>
/*High-resolution per-process timer from the CPU 
 * -lrt library should be linked*/
double cpuClock() {
    struct timespec t0;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
    return (double) t0.tv_sec+1.0e-9*(double) t0.tv_nsec;
}
void getWallTime(double* wcTime) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    *wcTime = (double)(tp.tv_sec + tp.tv_usec*1.0e-6);
}
#endif

void getWallTime_(double* wcTime) {
    getWallTime(wcTime);
}


