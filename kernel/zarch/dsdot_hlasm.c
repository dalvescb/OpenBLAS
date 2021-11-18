/***************************************************************************
Copyright (c) 2013-2019,The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms,with or without
modification,are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice,this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice,this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES,INCLUDING,BUT NOT LIMITED TO,THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT,INDIRECT,INCIDENTAL,SPECIAL,EXEMPLARY,OR CONSEQUENTIAL
DAMAGES (INCLUDING,BUT NOT LIMITED TO,PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE,DATA,OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY,WHETHER IN CONTRACT,STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE,EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

static double dsdot_kernel_16(BLASLONG n, FLOAT *x, FLOAT *y) {
  double dot;
__asm__ (
" vzero 0\n"
" vzero 1\n"
" vzero 2\n"
" vzero 3\n"
" vzero 4\n"
" vzero 5\n"
" vzero 6\n"
" vzero 7\n"
" srlg  %[n],%[n],4\n"
" xgr   1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" pfd 1,1024(1,%[y])\n"
" vlef  16,0(1,%[x]),0\n"
" vlef  16,4(1,%[x]),2\n"
" vlef  17,8(1,%[x]),0\n"
" vlef  17,12(1,%[x]),2\n"
" vlef  18,16(1,%[x]),0\n"
" vlef  18,20(1,%[x]),2\n"
" vlef  19,24(1,%[x]),0\n"
" vlef  19,28(1,%[x]),2\n"
" vlef  20,32(1,%[x]),0\n"
" vlef  20,36(1,%[x]),2\n"
" vlef  21,40(1,%[x]),0\n"
" vlef  21,44(1,%[x]),2\n"
" vlef  22,48(1,%[x]),0\n"
" vlef  22,52(1,%[x]),2\n"
" vlef  23,56(1,%[x]),0\n"
" vlef  23,60(1,%[x]),2\n"
" vflls   16,16\n"
" vflls   17,17\n"
" vflls   18,18\n"
" vflls   19,19\n"
" vflls   20,20\n"
" vflls   21,21\n"
" vflls   22,22\n"
" vflls   23,23\n"
" vlef    24,0(1,%[y]),0\n"
" vlef    24,4(1,%[y]),2\n"
" vflls   24,24\n"
" vfmadb  0,16,24,0\n"
" vlef    25,8(1,%[y]),0\n"
" vlef    25,12(1,%[y]),2\n"
" vflls   25,25\n"
" vfmadb  1,17,25,1\n"
" vlef    26,16(1,%[y]),0\n"
" vlef    26,20(1,%[y]),2\n"
" vflls   26,26\n"
" vfmadb  2,18,26,2\n"
" vlef    27,24(1,%[y]),0\n"
" vlef    27,28(1,%[y]),2\n"
" vflls   27,27\n"
" vfmadb  3,19,27,3\n"
" vlef    28,32(1,%[y]),0\n"
" vlef    28,36(1,%[y]),2\n"
" vflls   28,28\n"
" vfmadb  4,20,28,4\n"
" vlef    29,40(1,%[y]),0\n"
" vlef    29,44(1,%[y]),2\n"
" vflls   29,29\n"
" vfmadb  5,21,29,5\n"
" vlef    30,48(1,%[y]),0\n"
" vlef    30,52(1,%[y]),2\n"
" vflls   30,30\n"
" vfmadb  6,22,30,6\n"
" vlef    31,56(1,%[y]),0\n"
" vlef    31,60(1,%[y]),2\n"
" vflls   31,31\n"
" vfmadb  7,23,31,7\n"
" agfi   1,64\n"
" brctg  %[n],label0\n"
" vfadb  0,0,1\n"
" vfadb  0,0,2\n"
" vfadb  0,0,3\n"
" vfadb  0,0,4\n"
" vfadb  0,0,5\n"
" vfadb  0,0,6\n"
" vfadb  0,0,7\n"
" vrepg  1,0,1\n"
" adbr   0,1\n"
" ldr    %[dot],0\n"
    : [dot] "=f"(dot),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n]) x),[x] "a"(x),
       "m"(*(const FLOAT (*)[n]) y),[y] "a"(y)
    : "cc", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v16",
       "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26",
       "v27", "v28", "v29", "v30", "v31");
  return dot;
}

double CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;

  double dot = 0.0;

  if (n <= 0)
    return (dot);

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -16;

    if (n1)
      dot = dsdot_kernel_16(n1, x, y);

    i = n1;
    while (i < n) {

      dot += (double) y[i] * (double) x[i];
      i++;

    }
    return (dot);

  }

  BLASLONG n1 = n & -2;

  while (i < n1) {

    dot += (double) y[iy] * (double) x[ix];
    dot += (double) y[iy + inc_y] * (double) x[ix + inc_x];
    ix += inc_x * 2;
    iy += inc_y * 2;
    i += 2;

  }

  while (i < n) {

    dot += (double) y[iy] * (double) x[ix];
    ix += inc_x;
    iy += inc_y;
    i++;

  }
  return (dot);

}
