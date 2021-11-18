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

static void srot_kernel_64(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *c, FLOAT *s) {
__asm__ (
" vlrepf 0,%[c]\n"
" vlrepf 1,%[s]\n"
" srlg   %[n],%[n],6\n"
" xgr    1,1\n"
"label0 pfd 2,1024(1,%[x])\n"
" pfd 2,1024(1,%[y])\n"
" vl  24,0(1,%[x])\n"
" vl  25,16(1,%[x])\n"
" vl  26,32(1,%[x])\n"
" vl  27,48(1,%[x])\n"
" vl  16,0(1,%[y])\n"
" vl  17,16(1,%[y])\n"
" vl  18,32(1,%[y])\n"
" vl  19,48(1,%[y])\n"
" vfmsb 28,24,0\n"
" vfmsb 29,25,0\n"
" vfmsb 20,24,1\n"
" vfmsb 21,25,1\n"
" vfmsb 30,26,0\n"
" vfmsb 22,26,1\n"
" vfmsb 31,27,0\n"
" vfmsb 23,27,1\n"
" vfmasb 28,16,1,28\n"
" vfmssb 20,16,0,20\n"
" vfmasb 29,17,1,29\n"
" vfmssb 21,17,0,21\n"
" vfmasb 30,18,1,30\n"
" vfmssb 22,18,0,22\n"
" vfmasb 31,19,1,31\n"
" vfmssb 23,19,0,23\n"
" vst  28,0(1,%[x])\n"
" vst  29,16(1,%[x])\n"
" vst  30,32(1,%[x])\n"
" vst  31,48(1,%[x])\n"
" vst  20,0(1,%[y])\n"
" vst  21,16(1,%[y])\n"
" vst  22,32(1,%[y])\n"
" vst  23,48(1,%[y])\n"
" vl  24,64(1,%[x])\n"
" vl  25,80(1,%[x])\n"
" vl  26,96(1,%[x])\n"
" vl  27,112(1,%[x])\n"
" vl  16,64(1,%[y])\n"
" vl  17,80(1,%[y])\n"
" vl  18,96(1,%[y])\n"
" vl  19,112(1,%[y])\n"
" vfmsb 28,24,0\n"
" vfmsb 29,25,0\n"
" vfmsb 20,24,1\n"
" vfmsb 21,25,1\n"
" vfmsb 30,26,0\n"
" vfmsb 22,26,1\n"
" vfmsb 31,27,0\n"
" vfmsb 23,27,1\n"
" vfmasb 28,16,1,28\n"
" vfmssb 20,16,0,20\n"
" vfmasb 29,17,1,29\n"
" vfmssb 21,17,0,21\n"
" vfmasb 30,18,1,30\n"
" vfmssb 22,18,0,22\n"
" vfmasb 31,19,1,31\n"
" vfmssb 23,19,0,23\n"
" vst  28,64(1,%[x])\n"
" vst  29,80(1,%[x])\n"
" vst  30,96(1,%[x])\n"
" vst  31,112(1,%[x])\n"
" vst  20,64(1,%[y])\n"
" vst  21,80(1,%[y])\n"
" vst  22,96(1,%[y])\n"
" vst  23,112(1,%[y])\n"
" vl  24,128(1,%[x])\n"
" vl  25,144(1,%[x])\n"
" vl  26,160(1,%[x])\n"
" vl  27,176(1,%[x])\n"
" vl  16,128(1,%[y])\n"
" vl  17,144(1,%[y])\n"
" vl  18,160(1,%[y])\n"
" vl  19,176(1,%[y])\n"
" vfmsb 28,24,0\n"
" vfmsb 29,25,0\n"
" vfmsb 20,24,1\n"
" vfmsb 21,25,1\n"
" vfmsb 30,26,0\n"
" vfmsb 22,26,1\n"
" vfmsb 31,27,0\n"
" vfmsb 23,27,1\n"
" vfmasb 28,16,1,28\n"
" vfmssb 20,16,0,20\n"
" vfmasb 29,17,1,29\n"
" vfmssb 21,17,0,21\n"
" vfmasb 30,18,1,30\n"
" vfmssb 22,18,0,22\n"
" vfmasb 31,19,1,31\n"
" vfmssb 23,19,0,23\n"
" vst  28,128(1,%[x])\n"
" vst  29,144(1,%[x])\n"
" vst  30,160(1,%[x])\n"
" vst  31,176(1,%[x])\n"
" vst  20,128(1,%[y])\n"
" vst  21,144(1,%[y])\n"
" vst  22,160(1,%[y])\n"
" vst  23,176(1,%[y])\n"
" vl  24,192(1,%[x])\n"
" vl  25,208(1,%[x])\n"
" vl  26,224(1,%[x])\n"
" vl  27,240(1,%[x])\n"
" vl  16,192(1,%[y])\n"
" vl  17,208(1,%[y])\n"
" vl  18,224(1,%[y])\n"
" vl  19,240(1,%[y])\n"
" vfmsb 28,24,0\n"
" vfmsb 29,25,0\n"
" vfmsb 20,24,1\n"
" vfmsb 21,25,1\n"
" vfmsb 30,26,0\n"
" vfmsb 22,26,1\n"
" vfmsb 31,27,0\n"
" vfmsb 23,27,1\n"
" vfmasb 28,16,1,28\n"
" vfmssb 20,16,0,20\n"
" vfmasb 29,17,1,29\n"
" vfmssb 21,17,0,21\n"
" vfmasb 30,18,1,30\n"
" vfmssb 22,18,0,22\n"
" vfmasb 31,19,1,31\n"
" vfmssb 23,19,0,23\n"
" vst  28,192(1,%[x])\n"
" vst  29,208(1,%[x])\n"
" vst  30,224(1,%[x])\n"
" vst  31,240(1,%[x])\n"
" vst  20,192(1,%[y])\n"
" vst  21,208(1,%[y])\n"
" vst  22,224(1,%[y])\n"
" vst  23,240(1,%[y])\n"
" agfi  1,256\n"
" brctg %[n],label0\n"
    : "+m"(*(FLOAT (*)[n]) x), "+m"(*(FLOAT (*)[n]) y),
       [n] "+&r"(n)
    : [x] "a"(x),[y] "a"(y),[c] "Q"(*c),[s] "Q"(*s)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

int CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y,
          FLOAT c, FLOAT s) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;

  FLOAT temp;

  if (n <= 0)
    return (0);

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -64;
    if (n1 > 0) {
      FLOAT cosa, sina;
      cosa = c;
      sina = s;
      srot_kernel_64(n1, x, y, &cosa, &sina);
      i = n1;
    }

    while (i < n) {
      temp = c * x[i] + s * y[i];
      y[i] = c * y[i] - s * x[i];
      x[i] = temp;

      i++;

    }

  } else {

    while (i < n) {
      temp = c * x[ix] + s * y[iy];
      y[iy] = c * y[iy] - s * x[ix];
      x[ix] = temp;

      ix += inc_x;
      iy += inc_y;
      i++;

    }

  }
  return (0);

}
