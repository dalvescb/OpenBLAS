/***************************************************************************
Copyright (c) 2013-2019, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

static void zdot_kernel_8(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *d) {
__asm__ (
" vzero 24\n"
" vzero 25\n"
" vzero 26\n"
" vzero 27\n"
" vzero 28\n"
" vzero 29\n"
" vzero 30\n"
" vzero 31\n"
" srlg %[n],%[n],3\n"
" xgr 1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" pfd 1,1024(1,%[y])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  0,0(1,%[y])\n"
" vl  1,16(1,%[y])\n"
" vl  2,32(1,%[y])\n"
" vl  3,48(1,%[y])\n"
" vpdi 20,16,16,4\n"
" vpdi 21,17,17,4\n"
" vpdi 22,18,18,4\n"
" vpdi 23,19,19,4\n"
" vfmadb    24,16,0,24\n"
" vfmadb    25,20,0,25\n"
" vfmadb    26,17,1,26\n"
" vfmadb    27,21,1,27\n"
" vfmadb    28,18,2,28\n"
" vfmadb    29,22,2,29\n"
" vfmadb    30,19,3,30\n"
" vfmadb    31,23,3,31\n"
" vl  16,64(1,%[x])\n"
" vl  17,80(1,%[x])\n"
" vl  18,96(1,%[x])\n"
" vl  19,112(1,%[x])\n"
" vl  0,64(1,%[y])\n"
" vl  1,80(1,%[y])\n"
" vl  2,96(1,%[y])\n"
" vl  3,112(1,%[y])\n"
" vpdi 20,16,16,4\n"
" vpdi 21,17,17,4\n"
" vpdi 22,18,18,4\n"
" vpdi 23,19,19,4\n"
" vfmadb    24,16,0,24\n"
" vfmadb    25,20,0,25\n"
" vfmadb    26,17,1,26\n"
" vfmadb    27,21,1,27\n"
" vfmadb    28,18,2,28\n"
" vfmadb    29,22,2,29\n"
" vfmadb    30,19,3,30\n"
" vfmadb    31,23,3,31\n"
" agfi   1,128\n"
" brctg  %[n],label0\n"
" vfadb  24,24,26\n"
" vfadb  24,24,28\n"
" vfadb  24,24,30\n"
" vfadb  25,25,27\n"
" vfadb  25,25,29\n"
" vfadb  25,25,31\n"
" vsteg  24,0(%[d]),0\n"
" vsteg  24,8(%[d]),1\n"
" vsteg  25,16(%[d]),1\n"
" vsteg  25,24(%[d]),0\n"
    : "=m"(*(FLOAT (*)[4]) d),[n] "+&r"(n)
    : [d] "a"(d), "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x),
       "m"(*(const FLOAT (*)[n * 2]) y),[y] "a"(y)
    : "cc", "r1", "v0", "v1", "v2", "v3", "v16", "v17", "v18", "v19", "v20",
       "v21", "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

OPENBLAS_COMPLEX_FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y,
                             BLASLONG inc_y) {
  BLASLONG i;
  BLASLONG ix, iy;
  OPENBLAS_COMPLEX_FLOAT result;
  FLOAT dot[4] __attribute__ ((aligned(16))) = {
  0.0, 0.0, 0.0, 0.0};

  if (n <= 0) {
    CREAL(result) = 0.0;
    CIMAG(result) = 0.0;
    return (result);

  }

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -8;

    if (n1)
      zdot_kernel_8(n1, x, y, dot);

    i = n1;
    BLASLONG j = i * 2;

    while (i < n) {

      dot[0] += x[j] * y[j];
      dot[1] += x[j + 1] * y[j + 1];
      dot[2] += x[j] * y[j + 1];
      dot[3] += x[j + 1] * y[j];

      j += 2;
      i++;

    }

  } else {
    i = 0;
    ix = 0;
    iy = 0;
    inc_x <<= 1;
    inc_y <<= 1;
    while (i < n) {

      dot[0] += x[ix] * y[iy];
      dot[1] += x[ix + 1] * y[iy + 1];
      dot[2] += x[ix] * y[iy + 1];
      dot[3] += x[ix + 1] * y[iy];

      ix += inc_x;
      iy += inc_y;
      i++;

    }
  }

#if !defined(CONJ)
  CREAL(result) = dot[0] - dot[1];
  CIMAG(result) = dot[2] + dot[3];
#else
  CREAL(result) = dot[0] + dot[1];
  CIMAG(result) = dot[2] - dot[3];

#endif

  return (result);

}
