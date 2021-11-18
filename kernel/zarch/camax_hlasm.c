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

#define CABS1(x,i) (fabsf(x[i]) + fabsf(x[i + 1]))

static FLOAT camax_kernel_32(BLASLONG n, FLOAT *x) {
  FLOAT amax;
__asm__ (
" vlef   0,0(%[x]),0\n"
" vlef   16,4(%[x]),0\n"
" vlef   0,8(%[x]),1\n"
" vlef   16,12(%[x]),1\n"
" vlef   0,16(%[x]),2\n"
" vlef   16,20(%[x]),2\n"
" vlef   0,24(%[x]),3\n"
" vlef   16,28(%[x]),3\n"
" vflpsb 0,0\n"
" vflpsb 16,16\n"
" vfasb  0,0,16\n"
" vleib  1,0,0\n"
" vleib  1,1,1\n"
" vleib  1,2,2\n"
" vleib  1,3,3\n"
" vleib  1,8,4\n"
" vleib  1,9,5\n"
" vleib  1,10,6\n"
" vleib  1,11,7\n"
" vleib  1,16,8\n"
" vleib  1,17,9\n"
" vleib  1,18,10\n"
" vleib  1,19,11\n"
" vleib  1,24,12\n"
" vleib  1,25,13\n"
" vleib  1,26,14\n"
" vleib  1,27,15\n"
" srlg  %[n],%[n],5\n"
" xgr 1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" vl    16,0(1,%[x])\n"
" vl    2,16(1,%[x])\n"
" vpkg  17,16,2\n"
" vperm 16,16,2,1\n"
" vl    18,32(1,%[x])\n"
" vl    2,48(1,%[x])\n"
" vpkg  19,18,2\n"
" vperm 18,18,2,1\n"
" vl    20,64(1,%[x])\n"
" vl    2,80(1,%[x])\n"
" vpkg  21,20,2\n"
" vperm 20,20,2,1\n"
" vl    22,96(1,%[x])\n"
" vl    2,112(1,%[x])\n"
" vpkg  23,22,2\n"
" vperm 22,22,2,1\n"
" vl    24,128(1,%[x])\n"
" vl    2,144(1,%[x])\n"
" vpkg  25,24,2\n"
" vperm 24,24,2,1\n"
" vl    26,160(1,%[x])\n"
" vl    2,176(1,%[x])\n"
" vpkg  27,26,2\n"
" vperm 26,26,2,1\n"
" vl    28,192(1,%[x])\n"
" vl    2,208(1,%[x])\n"
" vpkg  29,28,2\n"
" vperm 28,28,2,1\n"
" vl    30,224(1,%[x])\n"
" vl    2,240(1,%[x])\n"
" vpkg  31,30,2\n"
" vperm 30,30,2,1\n"
" vflpsb  16,16\n"
" vflpsb  17,17\n"
" vflpsb  18,18\n"
" vflpsb  19,19\n"
" vflpsb  20,20\n"
" vflpsb  21,21\n"
" vflpsb  22,22\n"
" vflpsb  23,23\n"
" vflpsb  24,24\n"
" vflpsb  25,25\n"
" vflpsb  26,26\n"
" vflpsb  27,27\n"
" vflpsb  28,28\n"
" vflpsb  29,29\n"
" vflpsb  30,30\n"
" vflpsb  31,31\n"
" vfasb 16,16,17\n"
" vfasb 18,18,19\n"
" vfasb 20,20,21\n"
" vfasb 22,22,23\n"
" vfasb 24,24,25\n"
" vfasb 26,26,27\n"
" vfasb 28,28,29\n"
" vfasb 30,30,31\n"
" vfmaxsb  16,16,24,0\n"
" vfmaxsb  18,18,26,0\n"
" vfmaxsb  20,20,28,0\n"
" vfmaxsb  22,22,30,0\n"
" vfmaxsb  16,16,20,0\n"
" vfmaxsb  18,18,22,0\n"
" vfmaxsb  16,16,18,0\n"
" vfmaxsb  0,0,16,0\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" veslg   16,0,32\n"
" vfmaxsb 0,0,16,0\n"
" vrepf  16,0,2\n"
" wfmaxsb 0,0,16,0\n"
" ler    %[amax],0\n"
    : [amax] "=f"(amax),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v2", "v16", "v17", "v18", "v19", "v20",
       "v21", "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
  return amax;
}

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG ix = 0;
  FLOAT maxf = 0.0;
  BLASLONG inc_x2;

  if (n <= 0 || inc_x <= 0)
    return (maxf);

  if (inc_x == 1) {

    BLASLONG n1 = n & -32;
    if (n1 > 0) {

      maxf = camax_kernel_32(n1, x);
      ix = n1 * 2;
      i = n1;
    } else {
      maxf = CABS1(x, 0);
      ix += 2;
      i++;
    }

    while (i < n) {
      if (CABS1(x, ix) > maxf) {
        maxf = CABS1(x, ix);
      }
      ix += 2;
      i++;
    }
    return (maxf);

  } else {

    maxf = CABS1(x, 0);
    inc_x2 = 2 * inc_x;

    BLASLONG n1 = n & -4;
    while (i < n1) {

      if (CABS1(x, ix) > maxf) {
        maxf = CABS1(x, ix);
      }
      if (CABS1(x, ix + inc_x2) > maxf) {
        maxf = CABS1(x, ix + inc_x2);
      }
      if (CABS1(x, ix + inc_x2 * 2) > maxf) {
        maxf = CABS1(x, ix + inc_x2 * 2);
      }
      if (CABS1(x, ix + inc_x2 * 3) > maxf) {
        maxf = CABS1(x, ix + inc_x2 * 3);
      }

      ix += inc_x2 * 4;

      i += 4;

    }

    while (i < n) {
      if (CABS1(x, ix) > maxf) {
        maxf = CABS1(x, ix);
      }
      ix += inc_x2;
      i++;
    }
    return (maxf);
  }
}
