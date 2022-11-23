////////////////////////////////////////////////////////////////////////////
// Test program to exercise print function called through function pointers
// by calling OpenBLAS functions with invalid arguments.  Since some of these
// calls stop the program, use a script to set the following macros for each
// running of this program, and run once per test.
//
// DEFAULT:
//   defined means to use default print functions,
//   undefined means to use user print functions.
// TEST:
//   value from 1 to 5 to select which test to run.
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdarg.h>
#include <sys/time.h>
#include <time.h>
#include "cblas.h"
#include "common.h"

extern int errno;

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
    printf ("my_perror: %s errno=%d\n", msg, errno);
  }

main(){
  int i;

/*
  printf ("sizeof(blasint)=%d\n", sizeof(blasint));
  printf ("sizeof(long)=%d\n", sizeof(long));
  printf ("sizeof(long long)=%d\n", sizeof(long long));
  printf ("sizeof(size_t)=%d\n", sizeof(size_t));
*/

  #ifdef DEFAULT
    // Use default print functions
    printf ("\nTEST=%1d DEFAULT\n", TEST);
  #else
    // Set to use user print functions
    printf ("\nTEST=%1d USER\n", TEST);
    __openblas_set_printf (&my_printf);
    __openblas_set_fprintf (&my_fprintf);
    __openblas_set_vfprintf (&my_vfprintf);
    __openblas_set_perror (&my_perror);
  #endif

  ////////////////////////////////////////////////////////////////////////////
  // cblas_xerbla()
  ////////////////////////////////////////////////////////////////////////////

  #if TEST==1
    printf ("Calling cblas_xerbla()\n");
    cblas_xerbla (3, "foo", "parm=%f\n", 1.23);  // calling this ends program!
  #endif


  ////////////////////////////////////////////////////////////////////////////
  // cblas_sgemv() to test xerbla()
  // void cblas_sgemv(OPENBLAS_CONST enum CBLAS_ORDER order,
  //             OPENBLAS_CONST enum CBLAS_TRANSPOSE trans,
  //             OPENBLAS_CONST blasint m, OPENBLAS_CONST blasint n,
  //             OPENBLAS_CONST float alpha, OPENBLAS_CONST float *a,
  //             OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *x,
  //             OPENBLAS_CONST blasint incx, OPENBLAS_CONST float beta,
  //             float *y, OPENBLAS_CONST blasint incy);
  ////////////////////////////////////////////////////////////////////////////
  #if TEST==2
  {
    CBLAS_ORDER order = CblasRowMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    blasint m=-1, n=1;  // invalid m
    float alpha=1.0;
    float a[1][1]={{1.0}};
    blasint lda=m;
    float x[1]={1.0};
    blasint incx=1;
    float beta=1.0;
    float y[1];
    blasint incy=1;

    printf ("Calling cblas_sgemv() with m=-1 to test xerbla()\n");
    cblas_sgemv (order, trans, m, n, alpha, &(a[0][0]), lda, x, incx, beta, y, incy);
  }

  // Calling Fortran interface sgemv_()
  {
    float a[1][1]={{1.0}}, x[1]={1.0}, y[1];
    float alpha=1.0;
    float beta=1.0;
    char trans='N';
    blasint m=-1, n=1;  // invalid m
    blasint lda=m;
    blasint inc_x=1,inc_y=1;

    printf ("Calling sgemv_() with m=-1 to test xerbla()\n");
    sgemv_ (&trans, &m, &n, &alpha, &(a[0][0]), &lda, x, &inc_x, &beta, y, &inc_y);
  }
  #endif

  ////////////////////////////////////////////////////////////////////////////
  // cblas_simatcopy() to test PRINTF("Memory alloc failed\n");
  // void cblas_simatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
  //             OPENBLAS_CONST enum CBLAS_TRANSPOSE ctrans,
  //             OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
  //             OPENBLAS_CONST float calpha, float *a,
  //             OPENBLAS_CONST blasint clda, OPENBLAS_CONST blasint cldb);
  ////////////////////////////////////////////////////////////////////////////
  #if TEST==3
  {
    CBLAS_ORDER order = CblasRowMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    // Need rows != cols to enter imatcopy.c code with malloc()
    blasint rows=0x0FFFFFFF, cols=0x0FFFFFFE;  // rows*cols > 2^61
    float alpha=1.0;
    float a[1][1]={{1.0}};
    blasint lda=rows, ldb=cols;

    printf ("Calling cblas_simatcopy() with huge matrix to test Memory alloc fail\n");
    // This error cause program to stop!
    cblas_simatcopy (order, trans, rows, cols, alpha, &(a[0][0]), lda, ldb);
  }
  #endif

  ////////////////////////////////////////////////////////////////////////////
  // cblas_zimatcopy() to test PRINTF("Memory alloc failed\n");
  // void cblas_zimatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER,
  //             OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS,
  //             OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols,
  //             OPENBLAS_CONST double *calpha, double *a,
  //             OPENBLAS_CONST blasint clda, OPENBLAS_CONST blasint cldb);
  ////////////////////////////////////////////////////////////////////////////
  #if TEST==4
  {
    CBLAS_ORDER order = CblasRowMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    // Need rows != cols to enter imatcopy.c code with malloc()
    blasint rows=0x0FFFFFFF, cols=0x0FFFFFFE;  // rows*cols > 2^61
    double alpha[]={1.0,0.0};
    double a[1][1]={{1.0}};
    blasint lda=rows, ldb=cols;

    printf ("Calling cblas_zimatcopy() with huge matrix to test Memory alloc fail\n");
    // This error cause program to stop!
    cblas_zimatcopy (order, trans, rows, cols, alpha, &(a[0][0]), lda, ldb);
  }
  #endif


  ////////////////////////////////////////////////////////////////////////////
  //  openblas_warning()
  ////////////////////////////////////////////////////////////////////////////
  #if TEST==5
  {
    extern void openblas_warning (int, const char *);
    int verbose=0;
    printf ("Calling openblas_warning()\n");
    openblas_warning (verbose, "error msg from openblas_warning()\n");
  }
  #endif

  ////////////////////////////////////////////////////////////////////////////
  //  driver/others/memory.c has several error calls.  Not including any tests
  //  here because most (if not all) are under #ifdef DEBUG or an arch other
  //  than z/OS.
  ////////////////////////////////////////////////////////////////////////////
}
