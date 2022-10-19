#ifdef CBLAS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "common.h"

void CNAME(blasint p, char *rout, char *form, ...)
{
   va_list args;

   va_start(args, form);

   if (p)
      FPRINTF(stderr, "Parameter %d to routine %s was incorrect\n", p, rout);
   VFPRINTF(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

