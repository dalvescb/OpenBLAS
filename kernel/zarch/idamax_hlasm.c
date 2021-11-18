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
#include <math.h>

#define ABS fabs

static BLASLONG idamax_kernel_32(BLASLONG n, FLOAT *x, FLOAT *amax) {
  BLASLONG iamax;
__asm__ (
" vl     0,0(%[x])\n"
" vflpdb 0,0\n"
" vleig  1,0,0\n"
" vleig  1,1,1\n"
" vrepig 2,16\n"
" vzero  3\n"
" vleig  24,0,0\n"
" vleig  24,1,1\n"
" vleig  25,2,0\n"
" vleig  25,3,1\n"
" vleig  26,4,0\n"
" vleig  26,5,1\n"
" vleig  27,6,0\n"
" vleig  27,7,1\n"
" vleig  28,8,0\n"
" vleig  28,9,1\n"
" vleig  29,10,0\n"
" vleig  29,11,1\n"
" vleig  30,12,0\n"
" vleig  30,13,1\n"
" vleig  31,14,0\n"
" vleig  31,15,1\n"
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
" vflpdb  16,16\n"
" vflpdb  17,17\n"
" vflpdb  18,18\n"
" vflpdb  19,19\n"
" vflpdb  20,20\n"
" vflpdb  21,21\n"
" vflpdb  22,22\n"
" vflpdb  23,23\n"
" vfchedb  4,16,17\n"
" vfchedb  5,18,19\n"
" vfchedb  6,20,21\n"
" vfchedb  7,22,23\n"
" vsel    16,16,17,4\n"
" vsel    4,24,25,4\n"
" vsel    17,18,19,5\n"
" vsel    5,26,27,5\n"
" vsel    18,20,21,6\n"
" vsel    6,28,29,6\n"
" vsel    19,22,23,7\n"
" vsel    7,30,31,7\n"
" vfchedb  20,16,17\n"
" vfchedb  21,18,19\n"
" vsel    16,16,17,20\n"
" vsel    4,4,5,20\n"
" vsel    17,18,19,21\n"
" vsel    5,6,7,21\n"
" vfchedb  18,16,17\n"
" vsel    16,16,17,18\n"
" vsel    4,4,5,18\n"
" vag     4,4,3\n"
" vfchedb  5,0,16\n"
" vsel    0,0,16,5\n"
" vsel    1,1,4,5\n"
" vag     3,3,2\n"
" vl  16,128(1,%[x])\n"
" vl  17,144(1,%[x])\n"
" vl  18,160(1,%[x])\n"
" vl  19,176(1,%[x])\n"
" vl  20,192(1,%[x])\n"
" vl  21,208(1,%[x])\n"
" vl  22,224(1,%[x])\n"
" vl  23,240(1,%[x])\n"
" vflpdb  16,16\n"
" vflpdb  17,17\n"
" vflpdb  18,18\n"
" vflpdb  19,19\n"
" vflpdb  20,20\n"
" vflpdb  21,21\n"
" vflpdb  22,22\n"
" vflpdb  23,23\n"
" vfchedb  4,16,17\n"
" vfchedb  5,18,19\n"
" vfchedb  6,20,21\n"
" vfchedb  7,22,23\n"
" vsel    16,16,17,4\n"
" vsel    4,24,25,4\n"
" vsel    17,18,19,5\n"
" vsel    5,26,27,5\n"
" vsel    18,20,21,6\n"
" vsel    6,28,29,6\n"
" vsel    19,22,23,7\n"
" vsel    7,30,31,7\n"
" vfchedb  20,16,17\n"
" vfchedb  21,18,19\n"
" vsel    16,16,17,20\n"
" vsel    4,4,5,20\n"
" vsel    17,18,19,21\n"
" vsel    5,6,7,21\n"
" vfchedb  18,16,17\n"
" vsel    16,16,17,18\n"
" vsel    4,4,5,18\n"
" vag     4,4,3\n"
" vfchedb  5,0,16\n"
" vsel    0,0,16,5\n"
" vsel    1,1,4,5\n"
" vag     3,3,2\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" vrepg  2,0,1\n"
" vrepg  3,1,1\n"
" wfcdb  2,0\n"
" jne label1\n"
" vsteg  0,%[amax],0\n"
" vmnlg  0,1,3\n"
" vlgvg  %[iamax],0,0\n"
" j label2\n"
"label1 wfchdb 4,2,0\n"
" vsel   1,3,1,4\n"
" vsel   0,2,0,4\n"
" std    0,%[amax]\n"
" vlgvg  %[iamax],1,0\n"
"label2 nop 0\n"
    : [iamax] "=r"(iamax),[amax] "=Q"(*amax),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v16",
       "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26",
       "v27", "v28", "v29", "v30", "v31");
  return iamax;
}

BLASLONG CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG j = 0;
  FLOAT maxf = 0.0;
  BLASLONG max = 0;

  if (n <= 0 || inc_x <= 0)
    return (max);

  if (inc_x == 1) {

    BLASLONG n1 = n & -32;
    if (n1 > 0) {

      max = idamax_kernel_32(n1, x, &maxf);

      i = n1;
    } else {
      maxf = ABS(x[0]);
      i++;
    }

    while (i < n) {
      if (ABS(x[i]) > maxf) {
        max = i;
        maxf = ABS(x[i]);
      }
      i++;
    }
    return (max + 1);

  } else {

    max = 0;
    maxf = ABS(x[0]);

    BLASLONG n1 = n & -4;
    while (j < n1) {

      if (ABS(x[i]) > maxf) {
        max = j;
        maxf = ABS(x[i]);
      }
      if (ABS(x[i + inc_x]) > maxf) {
        max = j + 1;
        maxf = ABS(x[i + inc_x]);
      }
      if (ABS(x[i + 2 * inc_x]) > maxf) {
        max = j + 2;
        maxf = ABS(x[i + 2 * inc_x]);
      }
      if (ABS(x[i + 3 * inc_x]) > maxf) {
        max = j + 3;
        maxf = ABS(x[i + 3 * inc_x]);
      }

      i += inc_x * 4;

      j += 4;

    }

    while (j < n) {
      if (ABS(x[i]) > maxf) {
        max = j;
        maxf = ABS(x[i]);
      }
      i += inc_x;
      j++;
    }
    return (max + 1);
  }
}
