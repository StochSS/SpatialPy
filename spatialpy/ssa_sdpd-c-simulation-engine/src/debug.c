#include "debug.h"

const char *ssa_level(int setlevel){
    switch (setlevel){
        case 1:
            return "INFO";
        case 2:
            return "DEBUG";
        default:
            return "UNKNOWN LEVEL";
    }
}

extern void ssa_printf(int setlevel, int debuglevel, const char *date, const char *time, const char *filename, int line, const char *fmt, ...){
    if (debuglevel >= setlevel){
        va_list vaargs;
        fprintf(stdout, "%s %s - ssa-sdpd.%s - %s - %i - ", date, time, filename, ssa_level(setlevel), line);
        va_start(vaargs, fmt);
        vfprintf(stdout, fmt, vaargs);
        va_end(vaargs);
        fflush(stdout);
    }
}
