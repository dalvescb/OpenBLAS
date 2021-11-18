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

static void sscal_kernel_32(BLASLONG n, FLOAT da, FLOAT *x) {
__asm__ (
" vlrepf 0,%[da]\n"
" srlg  %[n],%[n],5\n"
" xgr   1,1\n"
"label0 pfd 2,1024(1,%[x])\n"
" vl    24,0(1,%[x])\n"
" vfmsb 24,24,0\n"
" vst   24,0(1,%[x])\n"
" vl    25,16(1,%[x])\n"
" vfmsb 25,25,0\n"
" vst   25,16(1,%[x])\n"
" vl    26,32(1,%[x])\n"
" vfmsb 26,26,0\n"
" vst   26,32(1,%[x])\n"
" vl    27,48(1,%[x])\n"
" vfmsb 27,27,0\n"
" vst   27,48(1,%[x])\n"
" vl    28,64(1,%[x])\n"
" vfmsb 28,28,0\n"
" vst   28,64(1,%[x])\n"
" vl    29,80(1,%[x])\n"
" vfmsb 29,29,0\n"
" vst   29,80(1,%[x])\n"
" vl    30,96(1,%[x])\n"
" vfmsb 30,30,0\n"
" vst   30,96(1,%[x])\n"
" vl    31,112(1,%[x])\n"
" vfmsb 31,31,0\n"
" vst   31,112(1,%[x])\n"
" agfi   1,128\n"
" brctg  %[n],label0\n"
    : "+m"(*(FLOAT (*)[n]) x),[n] "+&r"(n)
    : [x] "a"(x),[da] "Q"(da)
    : "cc", "r1", "v0", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

static void sscal_kernel_32_zero(BLASLONG n, FLOAT *x) {
__asm__ (
" vzero 0\n"
" srlg %[n],%[n],5\n"
" xgr   1,1\n"
"label1 pfd 2,1024(1,%[x])\n"
" vst  0,0(1,%[x])\n"
" vst  0,16(1,%[x])\n"
" vst  0,32(1,%[x])\n"
" vst  0,48(1,%[x])\n"
" vst  0,64(1,%[x])\n"
" vst  0,80(1,%[x])\n"
" vst  0,96(1,%[x])\n"
" vst  0,112(1,%[x])\n"
" agfi  1,128\n"
" brctg %[n],label1\n"
    : "=m"(*(FLOAT (*)[n]) x),[n] "+&r"(n)
    : [x] "a"(x)
    : "cc", "r1", "v0");
}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x,
          BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy,
          BLASLONG dummy2) {
  BLASLONG i = 0, j = 0;
  if (n <= 0 || inc_x <= 0)
    return (0);

  if (inc_x == 1) {

    if (da == 0.0) {

      BLASLONG n1 = n & -32;
      if (n1 > 0) {

        sscal_kernel_32_zero(n1, x);
        j = n1;
      }

      while (j < n) {

        x[j] = 0.0;
        j++;
      }

    } else {

      BLASLONG n1 = n & -32;
      if (n1 > 0) {
        sscal_kernel_32(n1, da, x);
        j = n1;
      }
      while (j < n) {

        x[j] = da * x[j];
        j++;
      }
    }

  } else {

    if (da == 0.0) {

      BLASLONG n1 = n & -2;

      while (j < n1) {

        x[i] = 0.0;
        x[i + inc_x] = 0.0;

        i += inc_x * 2;
        j += 2;

      }
      while (j < n) {

        x[i] = 0.0;
        i += inc_x;
        j++;
      }

    } else {
      BLASLONG n1 = n & -2;

      while (j < n1) {

        x[i] = da * x[i];
        x[i + inc_x] = da * x[i + inc_x];

        i += inc_x * 2;
        j += 2;

      }

      while (j < n) {

        x[i] = da * x[i];
        i += inc_x;
        j++;
      }
    }

  }
  return 0;

}
