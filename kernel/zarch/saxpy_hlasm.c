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

static void saxpy_kernel_64(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *alpha) {
__asm__ (
" vlrepf 0,%[alpha]\n"
" srlg  %[n],%[n],6\n"
" xgr   1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" pfd 2,1024(1,%[y])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,0(1,%[y])\n"
" vl  21,16(1,%[y])\n"
" vl  22,32(1,%[y])\n"
" vl  23,48(1,%[y])\n"
" vl  24,64(1,%[x])\n"
" vl  25,80(1,%[x])\n"
" vl  26,96(1,%[x])\n"
" vl  27,112(1,%[x])\n"
" vl  28,64(1,%[y])\n"
" vl  29,80(1,%[y])\n"
" vl  30,96(1,%[y])\n"
" vl  31,112(1,%[y])\n"
" vfmasb   16,0,16,20\n"
" vfmasb   17,0,17,21\n"
" vfmasb   18,0,18,22\n"
" vfmasb   19,0,19,23\n"
" vfmasb   24,0,24,28\n"
" vfmasb   25,0,25,29\n"
" vfmasb   26,0,26,30\n"
" vfmasb   27,0,27,31\n"
" vst  16,0(1,%[y])\n"
" vst  17,16(1,%[y])\n"
" vst  18,32(1,%[y])\n"
" vst  19,48(1,%[y])\n"
" vst  24,64(1,%[y])\n"
" vst  25,80(1,%[y])\n"
" vst  26,96(1,%[y])\n"
" vst  27,112(1,%[y])\n"
" vl  16,128(1,%[x])\n"
" vl  17,144(1,%[x])\n"
" vl  18,160(1,%[x])\n"
" vl  19,176(1,%[x])\n"
" vl  20,128(1,%[y])\n"
" vl  21,144(1,%[y])\n"
" vl  22,160(1,%[y])\n"
" vl  23,176(1,%[y])\n"
" vl  24,192(1,%[x])\n"
" vl  25,208(1,%[x])\n"
" vl  26,224(1,%[x])\n"
" vl  27,240(1,%[x])\n"
" vl  28,192(1,%[y])\n"
" vl  29,208(1,%[y])\n"
" vl  30,224(1,%[y])\n"
" vl  31,240(1,%[y])\n"
" vfmasb   16,0,16,20\n"
" vfmasb   17,0,17,21\n"
" vfmasb   18,0,18,22\n"
" vfmasb   19,0,19,23\n"
" vfmasb   24,0,24,28\n"
" vfmasb   25,0,25,29\n"
" vfmasb   26,0,26,30\n"
" vfmasb   27,0,27,31\n"
" vst  16,128(1,%[y])\n"
" vst  17,144(1,%[y])\n"
" vst  18,160(1,%[y])\n"
" vst  19,176(1,%[y])\n"
" vst  24,192(1,%[y])\n"
" vst  25,208(1,%[y])\n"
" vst  26,224(1,%[y])\n"
" vst  27,240(1,%[y])\n"
" agfi  1,256\n"
" brctg %[n],label0\n"
    : "+m"(*(FLOAT (*)[n]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n]) x),[x] "a"(x),
       [alpha] "Q"(*alpha)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");
}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x,
          BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy,
          BLASLONG dummy2) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;

  if (n <= 0)
    return 0;

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -64;

    if (n1)
      saxpy_kernel_64(n1, x, y, &da);

    i = n1;
    while (i < n) {

      y[i] += da * x[i];
      i++;

    }
    return 0;

  }

  BLASLONG n1 = n & -4;

  while (i < n1) {

    FLOAT m1 = da * x[ix];
    FLOAT m2 = da * x[ix + inc_x];
    FLOAT m3 = da * x[ix + 2 * inc_x];
    FLOAT m4 = da * x[ix + 3 * inc_x];

    y[iy] += m1;
    y[iy + inc_y] += m2;
    y[iy + 2 * inc_y] += m3;
    y[iy + 3 * inc_y] += m4;

    ix += inc_x * 4;
    iy += inc_y * 4;
    i += 4;

  }

  while (i < n) {

    y[iy] += da * x[ix];
    ix += inc_x;
    iy += inc_y;
    i++;

  }
  return 0;

}
