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

#define CABS1(x,i) (fabs(x[i]) + fabs(x[i + 1]))

static BLASLONG izamin_kernel_16(BLASLONG n, FLOAT *x, FLOAT *amin) {
  BLASLONG iamin;
__asm__ (
" vleg   0,0(%[x]),0\n"
" vleg   1,8(%[x]),0\n"
" vleg   0,16(%[x]),1\n"
" vleg   1,24(%[x]),1\n"
" vflpdb 0,0\n"
" vflpdb 1,1\n"
" vfadb  0,0,1\n"
" vleig  1,0,0\n"
" vleig  1,1,1\n"
" vrepig 2,8\n"
" vzero  3\n"
" vleig  24,0,0\n"
" vleig  24,1,1\n"
" vleig  25,2,0\n"
" vleig  25,3,1\n"
" vleig  26,4,0\n"
" vleig  26,5,1\n"
" vleig  27,6,0\n"
" vleig  27,7,1\n"
" srlg  %[n],%[n],4\n"
" xgr 1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" vleg  16,0(1,%[x]),0\n"
" vleg  17,8(1,%[x]),0\n"
" vleg  16,16(1,%[x]),1\n"
" vleg  17,24(1,%[x]),1\n"
" vleg  18,32(1,%[x]),0\n"
" vleg  19,40(1,%[x]),0\n"
" vleg  18,48(1,%[x]),1\n"
" vleg  19,56(1,%[x]),1\n"
" vleg  20,64(1,%[x]),0\n"
" vleg  21,72(1,%[x]),0\n"
" vleg  20,80(1,%[x]),1\n"
" vleg  21,88(1,%[x]),1\n"
" vleg  22,96(1,%[x]),0\n"
" vleg  23,104(1,%[x]),0\n"
" vleg  22,112(1,%[x]),1\n"
" vleg  23,120(1,%[x]),1\n"
" vflpdb  16,16\n"
" vflpdb  17,17\n"
" vflpdb  18,18\n"
" vflpdb  19,19\n"
" vflpdb  20,20\n"
" vflpdb  21,21\n"
" vflpdb  22,22\n"
" vflpdb  23,23\n"
" vfadb 16,16,17\n"
" vfadb 17,18,19\n"
" vfadb 18,20,21\n"
" vfadb 19,22,23\n"
" vfchedb  4,17,16\n"
" vfchedb  5,19,18\n"
" vsel    16,16,17,4\n"
" vsel    4,24,25,4\n"
" vsel    17,18,19,5\n"
" vsel    5,26,27,5\n"
" vfchedb  18,17,16\n"
" vsel    16,16,17,18\n"
" vsel    4,4,5,18\n"
" vag     4,4,3\n"
" vfchedb  5,16,0\n"
" vsel    0,0,16,5\n"
" vsel    1,1,4,5\n"
" vag     3,3,2\n"
" vleg  16,128(1,%[x]),0\n"
" vleg  17,136(1,%[x]),0\n"
" vleg  16,144(1,%[x]),1\n"
" vleg  17,152(1,%[x]),1\n"
" vleg  18,160(1,%[x]),0\n"
" vleg  19,168(1,%[x]),0\n"
" vleg  18,176(1,%[x]),1\n"
" vleg  19,184(1,%[x]),1\n"
" vleg  20,192(1,%[x]),0\n"
" vleg  21,200(1,%[x]),0\n"
" vleg  20,208(1,%[x]),1\n"
" vleg  21,216(1,%[x]),1\n"
" vleg  22,224(1,%[x]),0\n"
" vleg  23,232(1,%[x]),0\n"
" vleg  22,240(1,%[x]),1\n"
" vleg  23,248(1,%[x]),1\n"
" vflpdb  16,16\n"
" vflpdb  17,17\n"
" vflpdb  18,18\n"
" vflpdb  19,19\n"
" vflpdb  20,20\n"
" vflpdb  21,21\n"
" vflpdb  22,22\n"
" vflpdb  23,23\n"
" vfadb 16,16,17\n"
" vfadb 17,18,19\n"
" vfadb 18,20,21\n"
" vfadb 19,22,23\n"
" vfchedb  4,17,16\n"
" vfchedb  5,19,18\n"
" vsel    16,16,17,4\n"
" vsel    4,24,25,4\n"
" vsel    17,18,19,5\n"
" vsel    5,26,27,5\n"
" vfchedb  18,17,16\n"
" vsel    16,16,17,18\n"
" vsel    4,4,5,18\n"
" vag     4,4,3\n"
" vfchedb  5,16,0\n"
" vsel    0,0,16,5\n"
" vsel    1,1,4,5\n"
" vag     3,3,2\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" vrepg  2,0,1\n"
" vrepg  3,1,1\n"
" wfcdb  2,0\n"
" jne label1\n"
" vsteg  0,%[amin],0\n"
" vmnlg  0,1,3\n"
" vlgvg  %[iamin],0,0\n"
" j label2\n"
"label1 wfchdb 4,0,2\n"
" vsel   1,3,1,4\n"
" vsel   0,2,0,4\n"
" std    0,%[amin]\n"
" vlgvg  %[iamin],1,0\n"
"label2 nop 0\n"
    : [iamin] "=r"(iamin),[amin] "=Q"(*amin),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v16", "v17", "v18",
       "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26", "v27");
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

    BLASLONG n1 = n & -16;
    if (n1 > 0) {

      min = izamin_kernel_16(n1, x, &minf);
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
