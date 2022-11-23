#define AUTHORIZED_ZARCH_C

#ifdef __OPENBLAS_AUTHORIZED_ZOS
  #include "common.h"  // For typedefs

  // Definitions of pointers to functions not provided in
  // z/OS Authorized LE environment.  They are initialized statically
  // at program loading to the LE defaults.
  
  __openblas_printf_ptr_t __openblas_printf_ptr = printf;
  __openblas_fprintf_ptr_t __openblas_fprintf_ptr = fprintf;
  __openblas_vfprintf_ptr_t __openblas_vfprintf_ptr =vfprintf;
  __openblas_perror_ptr_t __openblas_perror_ptr = perror;

  // User-callable functions to override defaults with their own functions that
  // work in Authorized LE environment.

  void __openblas_set_printf (__openblas_printf_ptr_t ptr) {
    __openblas_printf_ptr = ptr;
  }
  void __openblas_set_fprintf (__openblas_fprintf_ptr_t ptr) {
    __openblas_fprintf_ptr = ptr;
  }
  void __openblas_set_vfprintf (__openblas_vfprintf_ptr_t ptr) {
    __openblas_vfprintf_ptr = ptr;
  }
  void __openblas_set_perror (__openblas_perror_ptr_t ptr) {
    __openblas_perror_ptr = ptr;
  }
#endif
