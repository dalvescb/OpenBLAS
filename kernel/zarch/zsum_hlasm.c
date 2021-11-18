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


static FLOAT zsum_kernel_16(BLASLONG n, FLOAT *x) {
  FLOAT sum;
__asm__ (
" vzero   24\n"
" vzero   25\n"
" vzero   26\n"
" vzero   27\n"
" vzero   28\n"
" vzero   29\n"
" vzero   30\n"
" vzero   31\n"
" srlg  %[n],%[n],4\n"
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
" vsteg   24,%[sum],0\n"
    : [sum] "=Q"(sum),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x)
    : "cc", "r1", "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23",
       "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");
  return sum;
}

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG ip = 0;
  FLOAT sumf = 0.0;
  BLASLONG n1;
  BLASLONG inc_x2;

  if (n <= 0 || inc_x <= 0)
    return (sumf);

  if (inc_x == 1) {

    n1 = n & -16;
    if (n1 > 0) {

      sumf = zsum_kernel_16(n1, x);
      i = n1;
      ip = 2 * n1;
    }

    while (i < n) {
      sumf += x[ip] + x[ip + 1];
      i++;
      ip += 2;
    }

  } else {
    inc_x2 = 2 * inc_x;

    while (i < n) {
      sumf += x[ip] + x[ip + 1];
      ip += inc_x2;
      i++;
    }

  }
  return (sumf);
}
