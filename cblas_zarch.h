#ifndef CBLAS_ZARCH
#define CBLAS_ZARCH

#ifdef ALLOW_AUTHORIZED_LE
  // User-callable functions to override defaults with their own functions that
  // work in Authorized LE environment.

  extern void __openblas_set_printf (__openblas_printf_ptr_t);
  extern void __openblas_set_fprintf (__openblas_fprintf_ptr_t);
  extern void __openblas_set_vfprintf (__openblas_vfprintf_ptr_t);
  extern void __openblas_set_perror (__openblas_perror_ptr_t);
#endif

#endif
