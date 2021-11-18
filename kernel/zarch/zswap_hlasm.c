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

static void zswap_kernel_16(BLASLONG n, FLOAT *x, FLOAT *y) {
__asm__ (
" srlg %[n],%[n],4\n"
" xgr  1,1\n"
"label0 pfd 2,1024(1,%[x])\n"
" pfd 2,1024(1,%[y])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,64(1,%[x])\n"
" vl  21,80(1,%[x])\n"
" vl  22,96(1,%[x])\n"
" vl  23,112(1,%[x])\n"
" vl  24,128(1,%[x])\n"
" vl  25,144(1,%[x])\n"
" vl  26,160(1,%[x])\n"
" vl  27,176(1,%[x])\n"
" vl  28,192(1,%[x])\n"
" vl  29,208(1,%[x])\n"
" vl  30,224(1,%[x])\n"
" vl  31,240(1,%[x])\n"
" vl  0,0(1,%[y])\n"
" vl  1,16(1,%[y])\n"
" vl  2,32(1,%[y])\n"
" vl  3,48(1,%[y])\n"
" vl  4,64(1,%[y])\n"
" vl  5,80(1,%[y])\n"
" vl  6,96(1,%[y])\n"
" vl  7,112(1,%[y])\n"
" vst  0,0(1,%[x])\n"
" vst  1,16(1,%[x])\n"
" vst  2,32(1,%[x])\n"
" vst  3,48(1,%[x])\n"
" vst  4,64(1,%[x])\n"
" vst  5,80(1,%[x])\n"
" vst  6,96(1,%[x])\n"
" vst  7,112(1,%[x])\n"
" vl  0,128(1,%[y])\n"
" vl  1,144(1,%[y])\n"
" vl  2,160(1,%[y])\n"
" vl  3,176(1,%[y])\n"
" vl  4,192(1,%[y])\n"
" vl  5,208(1,%[y])\n"
" vl  6,224(1,%[y])\n"
" vl  7,240(1,%[y])\n"
" vst  0,128(1,%[x])\n"
" vst  1,144(1,%[x])\n"
" vst  2,160(1,%[x])\n"
" vst  3,176(1,%[x])\n"
" vst  4,192(1,%[x])\n"
" vst  5,208(1,%[x])\n"
" vst  6,224(1,%[x])\n"
" vst  7,240(1,%[x])\n"
" vst  16,0(1,%[y])\n"
" vst  17,16(1,%[y])\n"
" vst  18,32(1,%[y])\n"
" vst  19,48(1,%[y])\n"
" vst  20,64(1,%[y])\n"
" vst  21,80(1,%[y])\n"
" vst  22,96(1,%[y])\n"
" vst  23,112(1,%[y])\n"
" vst  24,128(1,%[y])\n"
" vst  25,144(1,%[y])\n"
" vst  26,160(1,%[y])\n"
" vst  27,176(1,%[y])\n"
" vst  28,192(1,%[y])\n"
" vst  29,208(1,%[y])\n"
" vst  30,224(1,%[y])\n"
" vst  31,240(1,%[y])\n"
" agfi   1,256\n"
" brctg  %[n],label0\n"
    : "+m"(*(FLOAT (*)[n * 2]) x),
       "+m"(*(FLOAT (*)[n * 2]) y),[n] "+&r"(n)
    : [x] "a"(x),[y] "a"(y)
    : "cc", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v16",
       "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26",
       "v27", "v28", "v29", "v30", "v31");
}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT dummy3,
          FLOAT dummy4, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y,
          FLOAT *dummy, BLASLONG dummy2) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;
  FLOAT temp[2];
  BLASLONG inc_x2, inc_y2;

  if (n <= 0)
    return (0);

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -16;
    if (n1 > 0) {
      zswap_kernel_16(n1, x, y);
      i = n1;
      ix = 2 * n1;
      iy = 2 * n1;
    }

    while (i < n) {

      temp[0] = x[ix];
      temp[1] = x[ix + 1];
      x[ix] = y[iy];
      x[ix + 1] = y[iy + 1];
      y[iy] = temp[0];
      y[iy + 1] = temp[1];

      ix += 2;
      iy += 2;
      i++;

    }

  } else {

    inc_x2 = 2 * inc_x;
    inc_y2 = 2 * inc_y;

    while (i < n) {

      temp[0] = x[ix];
      temp[1] = x[ix + 1];
      x[ix] = y[iy];
      x[ix + 1] = y[iy + 1];
      y[iy] = temp[0];
      y[iy + 1] = temp[1];

      ix += inc_x2;
      iy += inc_y2;
      i++;

    }

  }
  return (0);

}
