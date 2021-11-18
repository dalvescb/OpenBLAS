/***************************************************************************
Copyright (c) 2019, The OpenBLAS Project
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

static BLASLONG icamin_kernel_32(BLASLONG n, FLOAT *x, FLOAT *amin) {
  BLASLONG iamin;
__asm__ (
" vlef   0,0(%[x]),0\n"
" vlef   1,4(%[x]),0\n"
" vlef   0,8(%[x]),1\n"
" vlef   1,12(%[x]),1\n"
" vlef   0,16(%[x]),2\n"
" vlef   1,20(%[x]),2\n"
" vlef   0,24(%[x]),3\n"
" vlef   1,28(%[x]),3\n"
" vflpsb 0,0\n"
" vflpsb 1,1\n"
" vfasb  0,0,1\n"
" vleig  1,0,0\n"
" vleig  1,2,1\n"
" vleig  2,1,0\n"
" vleig  2,3,1\n"
" vrepig 3,16\n"
" vzero  4\n"
" vleib  9,0,0\n"
" vleib  9,1,1\n"
" vleib  9,2,2\n"
" vleib  9,3,3\n"
" vleib  9,8,4\n"
" vleib  9,9,5\n"
" vleib  9,10,6\n"
" vleib  9,11,7\n"
" vleib  9,16,8\n"
" vleib  9,17,9\n"
" vleib  9,18,10\n"
" vleib  9,19,11\n"
" vleib  9,24,12\n"
" vleib  9,25,13\n"
" vleib  9,26,14\n"
" vleib  9,27,15\n"
" vleif  24,0,0\n"
" vleif  24,1,1\n"
" vleif  24,2,2\n"
" vleif  24,3,3\n"
" vleif  25,4,0\n"
" vleif  25,5,1\n"
" vleif  25,6,2\n"
" vleif  25,7,3\n"
" vleif  26,8,0\n"
" vleif  26,9,1\n"
" vleif  26,10,2\n"
" vleif  26,11,3\n"
" vleif  27,12,0\n"
" vleif  27,13,1\n"
" vleif  27,14,2\n"
" vleif  27,15,3\n"
" srlg  %[n],%[n],5\n"
" xgr 1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" vl    16,0(1,%[x])\n"
" vl    28,16(1,%[x])\n"
" vpkg  17,16,28\n"
" vperm 16,16,28,9\n"
" vl    18,32(1,%[x])\n"
" vl    29,48(1,%[x])\n"
" vpkg  19,18,29\n"
" vperm 18,18,29,9\n"
" vl    20,64(1,%[x])\n"
" vl    30,80(1,%[x])\n"
" vpkg  21,20,30\n"
" vperm 20,20,30,9\n"
" vl    22,96(1,%[x])\n"
" vl    31,112(1,%[x])\n"
" vpkg  23,22,31\n"
" vperm 22,22,31,9\n"
" vflpsb  16,16\n"
" vflpsb  17,17\n"
" vflpsb  18,18\n"
" vflpsb  19,19\n"
" vflpsb  20,20\n"
" vflpsb  21,21\n"
" vflpsb  22,22\n"
" vflpsb  23,23\n"
" vfasb 16,16,17\n"
" vfasb 17,18,19\n"
" vfasb 18,20,21\n"
" vfasb 19,22,23\n"
" vfchesb  5,17,16\n"
" vfchesb  6,19,18\n"
" vsel    16,16,17,5\n"
" vsel    5,24,25,5\n"
" vsel    17,18,19,6\n"
" vsel    6,26,27,6\n"
" vfchesb  18,17,16\n"
" vsel    16,16,17,18\n"
" vsel    5,5,6,18\n"
" vsegf   6,5\n"
" vesrlg  5,5,32\n"
" vag     5,5,4\n"
" vag     6,6,4\n"
" vfchesb  7,16,0\n"
" vsel    0,0,16,7\n"
" vsegf   8,7\n"
" vesrlg  7,7,32\n"
" vsegf   7,7\n"
" vsel    1,1,5,7\n"
" vsel    2,2,6,8\n"
" vag     4,4,3\n"
" vl    16,128(1,%[x])\n"
" vl    28,144(1,%[x])\n"
" vpkg  17,16,28\n"
" vperm 16,16,28,9\n"
" vl    18,160(1,%[x])\n"
" vl    29,176(1,%[x])\n"
" vpkg  19,18,29\n"
" vperm 18,18,29,9\n"
" vl    20,192(1,%[x])\n"
" vl    30,208(1,%[x])\n"
" vpkg  21,20,30\n"
" vperm 20,20,30,9\n"
" vl    22,224(1,%[x])\n"
" vl    31,240(1,%[x])\n"
" vpkg  23,22,31\n"
" vperm 22,22,31,9\n"
" vflpsb  16,16\n"
" vflpsb  17,17\n"
" vflpsb  18,18\n"
" vflpsb  19,19\n"
" vflpsb  20,20\n"
" vflpsb  21,21\n"
" vflpsb  22,22\n"
" vflpsb  23,23\n"
" vfasb 16,16,17\n"
" vfasb 17,18,19\n"
" vfasb 18,20,21\n"
" vfasb 19,22,23\n"
" vfchesb  5,17,16\n"
" vfchesb  6,19,18\n"
" vsel    16,16,17,5\n"
" vsel    5,24,25,5\n"
" vsel    17,18,19,6\n"
" vsel    6,26,27,6\n"
" vfchesb  18,17,16\n"
" vsel    16,16,17,18\n"
" vsel    5,5,6,18\n"
" vsegf   6,5\n"
" vesrlg  5,5,32\n"
" vag     5,5,4\n"
" vag     6,6,4\n"
" vfchesb  7,16,0\n"
" vsel    0,0,16,7\n"
" vsegf   8,7\n"
" vesrlg  7,7,32\n"
" vsegf   7,7\n"
" vsel    1,1,5,7\n"
" vsel    2,2,6,8\n"
" vag     4,4,3\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" veslg   3,0,32\n"
" vfchsb  4,3,0\n"
" vchlg   5,2,1\n"
" vfcesb  6,0,3\n"
" vn      5,5,6\n"
" vo      4,4,5\n"
" vsel    0,0,3,4\n"
" vesrlg  4,4,32\n"
" vsegf   4,4\n"
" vsel    1,1,2,4\n"
" vrepf  2,0,2\n"
" vrepg  3,1,1\n"
" wfcsb  2,0\n"
" jne label1\n"
" vstef  0,%[amin],0\n"
" vmnlg  0,1,3\n"
" vlgvg  %[iamin],0,0\n"
" j label2\n"
"label1 wfchsb 4,0,2\n"
" vesrlg 4,4,32\n"
" vsegf  4,4\n"
" vsel   1,3,1,4\n"
" vsel   0,2,0,4\n"
" ste    0,%[amin]\n"
" vlgvg  %[iamin],1,0\n"
"label2 nop 0\n"
    : [iamin] "=r"(iamin),[amin] "=Q"(*amin),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8",
       "v9", "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24",
       "v25", "v26", "v27", "v28", "v29", "v30", "v31");
  return iamin;
}

BLASLONG CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG ix = 0;
  FLOAT minf = 0;
  BLASLONG min = 0;
  BLASLONG inc_x2;

  if (n <= 0 || inc_x <= 0)
    return (min);

  if (inc_x == 1) {

    BLASLONG n1 = n & -32;
    if (n1 > 0) {

      min = icamin_kernel_32(n1, x, &minf);
      ix = n1 * 2;
      i = n1;
    } else {
      minf = CABS1(x, 0);
      ix += 2;
      i++;
    }

    while (i < n) {
      if (CABS1(x, ix) < minf) {
        min = i;
        minf = CABS1(x, ix);
      }
      ix += 2;
      i++;
    }
    return (min + 1);

  } else {

    min = 0;
    minf = CABS1(x, 0);
    inc_x2 = 2 * inc_x;

    BLASLONG n1 = n & -4;
    while (i < n1) {

      if (CABS1(x, ix) < minf) {
        min = i;
        minf = CABS1(x, ix);
      }
      if (CABS1(x, ix + inc_x2) < minf) {
        min = i + 1;
        minf = CABS1(x, ix + inc_x2);
      }
      if (CABS1(x, ix + 2 * inc_x2) < minf) {
        min = i + 2;
        minf = CABS1(x, ix + 2 * inc_x2);
      }
      if (CABS1(x, ix + 3 * inc_x2) < minf) {
        min = i + 3;
        minf = CABS1(x, ix + 3 * inc_x2);
      }

      ix += inc_x2 * 4;

      i += 4;

    }

    while (i < n) {
      if (CABS1(x, ix) < minf) {
        min = i;
        minf = CABS1(x, ix);
      }
      ix += inc_x2;
      i++;
    }
    return (min + 1);
  }
}
