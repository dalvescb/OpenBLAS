#include <stdio.h>
#include <stdarg.h>
#include <sys/time.h>
#include <time.h>
#include "cblas.h"
#include "common.h"

// User print functions

  int my_printf (const char *fmt, ...) {
    va_list args;
    va_start (args, fmt);
    printf ("my_printf: ");
    vprintf (fmt, args);
    va_end (args);
    return 0;
  }

  int my_fprintf (FILE *fp, const char *fmt, ...) {
    va_list args;
    va_start (args, fmt);
    fprintf (fp, "my_fprintf: ");
    vfprintf (fp, fmt, args);
    va_end (args);
    return 0;
  }

  int my_vfprintf (FILE *fp, const char *fmt, char **args) {
    fprintf (fp, "my_vfprintf: ");
    vfprintf (fp, fmt, args);
    return 0;
  }

  void my_perror (const char *msg) {
    extern int errno;
    printf ("my_perror: %s errno=%d\n", msg, errno);
  }

  // Wrapper to call vfprintf through function pointer but with ... args.
  void wrap_vfprintf (FILE *fp, const char *fmt, ...) {
    va_list args;
    va_start (args, fmt);
    (*__openblas_vfprintf_ptr) (fp, fmt, args);
    va_end (args);
  }

main(){
  FILE *fp;
  fp = fopen ("foo","w");

  // Call function pointers - should get default print funcs
  (*__openblas_printf_ptr)("printf: %08X %f\n", 0xDEADBEEF, 3.14);
  (*__openblas_fprintf_ptr)(fp, "fprintf: %08X %f\n", 0xDEADBEEF, 3.14);
  wrap_vfprintf (fp, "vfprintf: %08X %f\n", 0xDEADBEEF, 3.14);
  (*__openblas_perror_ptr)("perror: an error message");

  // Set to use user print functions
  __openblas_set_printf (&my_printf);
  __openblas_set_fprintf (&my_fprintf);
  __openblas_set_vfprintf (&my_vfprintf);
  __openblas_set_perror (&my_perror);

  // Call function pointers - should get user print funcs
  (*__openblas_printf_ptr)("printf: %08X %f\n", 0xDEADBEEF, 3.14);
  (*__openblas_fprintf_ptr)(fp, "fprintf: %08X %f\n", 0xDEADBEEF, 3.14);
  wrap_vfprintf (fp, "vfprintf: %08X %f\n", 0xDEADBEEF, 3.14);
  (*__openblas_perror_ptr)("perror: an error message");

  fclose(fp);
}
