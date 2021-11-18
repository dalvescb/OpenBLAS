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

static FLOAT smax_kernel_64(BLASLONG n, FLOAT *x) {
  FLOAT max;
__asm__ (
" vl     0,0(%[x])\n"
" srlg   %[n],%[n],6\n"
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
" vl  24,128(1,%[x])\n"
" vl  25,144(1,%[x])\n"
" vl  26,160(1,%[x])\n"
" vl  27,176(1,%[x])\n"
" vl  28,192(1,%[x])\n"
" vl  29,208(1,%[x])\n"
" vl  30,224(1,%[x])\n"
" vl  31,240(1,%[x])\n"
" vfmaxsb  16,16,24,0\n"
" vfmaxsb  17,17,25,0\n"
" vfmaxsb  18,18,26,0\n"
" vfmaxsb  19,19,27,0\n"
" vfmaxsb  20,20,28,0\n"
" vfmaxsb  21,21,29,0\n"
" vfmaxsb  22,22,30,0\n"
" vfmaxsb  23,23,31,0\n"
" vfmaxsb  16,16,20,0\n"
" vfmaxsb  17,17,21,0\n"
" vfmaxsb  18,18,22,0\n"
" vfmaxsb  19,19,23,0\n"
" vfmaxsb  16,16,18,0\n"
" vfmaxsb  17,17,19,0\n"
" vfmaxsb  16,16,17,0\n"
" vfmaxsb  0,0,16,0\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" veslg   16,0,32\n"
" vfmaxsb 0,0,16,0\n"
" vrepf   16,0,2\n"
" wfmaxsb 0,0,16,0\n"
" ler    %[max],0\n"
    : [max] "=f"(max),[n] "+&r"(n)
    : "m"(*(const  FLOAT(*)[n]) x),[x] "a"(x)
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

    BLASLONG n1 = n & -64;
    if (n1 > 0) {

      maxf = smax_kernel_64(n1, x);

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
