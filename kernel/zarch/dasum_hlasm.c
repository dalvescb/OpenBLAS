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

static FLOAT dasum_kernel_32(BLASLONG n, FLOAT *x) {
  FLOAT asum;
__asm__ (
" vzero   24\n"
" vzero   25\n"
" vzero   26\n"
" vzero   27\n"
" vzero   28\n"
" vzero   29\n"
" vzero   30\n"
" vzero   31\n"
" srlg  %[n],%[n],5\n"
" xgr 1,1\n"
"label0 pfd  1,1024(1,%[x])\n"
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
" vfadb   24,24,16\n"
" vfadb   25,25,17\n"
" vfadb   26,26,18\n"
" vfadb   27,27,19\n"
" vfadb   28,28,20\n"
" vfadb   29,29,21\n"
" vfadb   30,30,22\n"
" vfadb   31,31,23\n"
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
" vfadb   24,24,16\n"
" vfadb   25,25,17\n"
" vfadb   26,26,18\n"
" vfadb   27,27,19\n"
" vfadb   28,28,20\n"
" vfadb   29,29,21\n"
" vfadb   30,30,22\n"
" vfadb   31,31,23\n"
" agfi  1,256\n"
" brctg %[n],label0\n"
" vfadb   24,24,25\n"
" vfadb   24,24,26\n"
" vfadb   24,24,27\n"
" vfadb   24,24,28\n"
" vfadb   24,24,29\n"
" vfadb   24,24,30\n"
" vfadb   24,24,31\n"
" vrepg   25,24,1\n"
" vfadb   24,24,25\n"
" vsteg   24,%[asum],0\n"
    : [asum] "=Q"(asum),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n]) x),[x] "a"(x)
    : "cc", "r1", "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23",
       "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");
  return asum;
}

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG j = 0;
  FLOAT sumf = 0.0;
  BLASLONG n1;

  if (n <= 0 || inc_x <= 0)
    return sumf;

  if (inc_x == 1) {

    n1 = n & -32;

    if (n1 > 0) {

      sumf = dasum_kernel_32(n1, x);
      i = n1;
    }

    while (i < n) {
      sumf += ABS(x[i]);
      i++;
    }

  } else {
    BLASLONG n1 = n & -4;
    register FLOAT sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;
    while (j < n1) {

      sum1 += ABS(x[i]);
      sum2 += ABS(x[i + inc_x]);
      sum1 += ABS(x[i + 2 * inc_x]);
      sum2 += ABS(x[i + 3 * inc_x]);

      i += inc_x * 4;
      j += 4;

    }
    sumf = sum1 + sum2;
    while (j < n) {

      sumf += ABS(x[i]);
      i += inc_x;
      j++;
    }

  }
  return sumf;
}
