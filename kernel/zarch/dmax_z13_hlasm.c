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

static FLOAT dmax_kernel_32(BLASLONG n, FLOAT *x) {
  FLOAT max;

__asm__ (
" vl    0,0(%[x])\n"
" srlg  %[n],%[n],5\n"
" xgr 1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,64(1,%[x])\n"
" vl  21,80(1,%[x])\n"
" vl  22,96(1,%[x])\n"
" vl  23,112(1,%[x])\n"
" vfchdb  24,16,17\n"
" vfchdb  25,18,19\n"
" vfchdb  26,20,21\n"
" vfchdb  27,22,23\n"
" vsel    24,16,17,24\n"
" vsel    25,18,19,25\n"
" vsel    26,20,21,26\n"
" vsel    27,22,23,27\n"
" vfchdb  28,24,25\n"
" vfchdb  29,26,27\n"
" vsel    28,24,25,28\n"
" vsel    29,26,27,29\n"
" vfchdb  30,28,29\n"
" vsel    30,28,29,30\n"
" vfchdb  31,30,0\n"
" vsel    0,30,0,31\n"
" vl  16,128(1,%[x])\n"
" vl  17,144(1,%[x])\n"
" vl  18,160(1,%[x])\n"
" vl  19,176(1,%[x])\n"
" vl  20,192(1,%[x])\n"
" vl  21,208(1,%[x])\n"
" vl  22,224(1,%[x])\n"
" vl  23,240(1,%[x])\n"
" vfchdb  24,16,17\n"
" vfchdb  25,18,19\n"
" vfchdb  26,20,21\n"
" vfchdb  27,22,23\n"
" vsel    24,16,17,24\n"
" vsel    25,18,19,25\n"
" vsel    26,20,21,26\n"
" vsel    27,22,23,27\n"
" vfchdb  28,24,25\n"
" vfchdb  29,26,27\n"
" vsel    28,24,25,28\n"
" vsel    29,26,27,29\n"
" vfchdb  30,28,29\n"
" vsel    30,28,29,30\n"
" vfchdb  31,30,0\n"
" vsel    0,30,0,31\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" vrepg  16,0,1\n"
" wfchdb 17,0,16\n"
" vsel   0,0,16,17\n"
" ldr    %[max],0\n"
    : [max] "=f"(max),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");

  return max;
}

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG j = 0;
  FLOAT maxf = 0.0;

  if (n <= 0 || inc_x <= 0)
    return (maxf);

  if (inc_x == 1) {

    BLASLONG n1 = n & -32;
    if (n1 > 0) {

      maxf = dmax_kernel_32(n1, x);

      i = n1;
    } else {
      maxf = x[0];
      i++;
    }

    while (i < n) {
      if (x[i] > maxf) {
        maxf = x[i];
      }
      i++;
    }
    return (maxf);

  } else {

    maxf = x[0];

    BLASLONG n1 = n & -4;
    while (j < n1) {

      if (x[i] > maxf) {
        maxf = x[i];
      }
      if (x[i + inc_x] > maxf) {
        maxf = x[i + inc_x];
      }
      if (x[i + 2 * inc_x] > maxf) {
        maxf = x[i + 2 * inc_x];
      }
      if (x[i + 3 * inc_x] > maxf) {
        maxf = x[i + 3 * inc_x];
      }

      i += inc_x * 4;

      j += 4;

    }

    while (j < n) {
      if (x[i] > maxf) {
        maxf = x[i];
      }
      i += inc_x;
      j++;
    }
    return (maxf);
  }
}
