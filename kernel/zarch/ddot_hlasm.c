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

static FLOAT ddot_kernel_16(BLASLONG n, FLOAT *x, FLOAT *y) {
  FLOAT dot;
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
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,64(1,%[x])\n"
" vl  21,80(1,%[x])\n"
" vl  22,96(1,%[x])\n"
" vl  23,112(1,%[x])\n"
" vl  24,0(1,%[y])\n"
" vl  25,16(1,%[y])\n"
" vl  26,32(1,%[y])\n"
" vl  27,48(1,%[y])\n"
" vl  28,64(1,%[y])\n"
" vl  29,80(1,%[y])\n"
" vl  30,96(1,%[y])\n"
" vl  31,112(1,%[y])\n"
" vfmadb   0,16,24,0\n"
" vfmadb   1,17,25,1\n"
" vfmadb   2,18,26,2\n"
" vfmadb   3,19,27,3\n"
" vfmadb   4,20,28,4\n"
" vfmadb   5,21,29,5\n"
" vfmadb   6,22,30,6\n"
" vfmadb   7,23,31,7\n"
" agfi   1,128\n"
" brctg  %[n],label0\n"
" vfadb   0,0,1\n"
" vfadb   0,0,2\n"
" vfadb   0,0,3\n"
" vfadb   0,0,4\n"
" vfadb   0,0,5\n"
" vfadb   0,0,6\n"
" vfadb   0,0,7\n"
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

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;

  FLOAT dot = 0.0;

  if (n <= 0)
    return (dot);

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -16;

    if (n1)
      dot = ddot_kernel_16(n1, x, y);

    i = n1;
    while (i < n) {

      dot += y[i] * x[i];
      i++;

    }
    return (dot);

  }

  FLOAT temp1 = 0.0;
  FLOAT temp2 = 0.0;

  BLASLONG n1 = n & -4;

  while (i < n1) {

    FLOAT m1 = y[iy] * x[ix];
    FLOAT m2 = y[iy + inc_y] * x[ix + inc_x];

    FLOAT m3 = y[iy + 2 * inc_y] * x[ix + 2 * inc_x];
    FLOAT m4 = y[iy + 3 * inc_y] * x[ix + 3 * inc_x];

    ix += inc_x * 4;
    iy += inc_y * 4;

    temp1 += m1 + m3;
    temp2 += m2 + m4;

    i += 4;

  }

  while (i < n) {

    temp1 += y[iy] * x[ix];
    ix += inc_x;
    iy += inc_y;
    i++;

  }
  dot = temp1 + temp2;
  return (dot);

}
