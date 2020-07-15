#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#if DEBUG_LEVEL > 0
    // http://support.raisonance.com/content/how-remove-file-path
    #define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
    #define SSA_LOG(setlevel, ...) \
        do { ssa_printf(setlevel, DEBUG_LEVEL, __DATE__, __TIME__, __FILENAME__, __LINE__, ## __VA_ARGS__); } while(0)
#else
    #define SSA_LOG(setlevel, fmt, ...)
#endif

const char *ssa_level(int setlevel);
extern void ssa_printf(int setlevel, int debuglevel, const char *date, const char* time, const char *filename, int line, const char *fmt, ...);
