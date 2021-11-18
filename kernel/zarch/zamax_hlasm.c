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

#define CABS1(x,i) (fabs(x[i]) + fabs(x[i + 1]))

static FLOAT zamax_kernel_16(BLASLONG n, FLOAT *x) {
  FLOAT amax;
__asm__ (
" vleg   0,0(%[x]),0\n"
" vleg   16,8(%[x]),0\n"
" vleg   0,16(%[x]),1\n"
" vleg   16,24(%[x]),1\n"
" vflpdb 0,0\n"
" vflpdb 16,16\n"
" vfadb  0,0,16\n"
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
" vleg  24,128(1,%[x]),0\n"
" vleg  25,136(1,%[x]),0\n"
" vleg  24,144(1,%[x]),1\n"
" vleg  25,152(1,%[x]),1\n"
" vleg  26,160(1,%[x]),0\n"
" vleg  27,168(1,%[x]),0\n"
" vleg  26,176(1,%[x]),1\n"
" vleg  27,184(1,%[x]),1\n"
" vleg  28,192(1,%[x]),0\n"
" vleg  29,200(1,%[x]),0\n"
" vleg  28,208(1,%[x]),1\n"
" vleg  29,216(1,%[x]),1\n"
" vleg  30,224(1,%[x]),0\n"
" vleg  31,232(1,%[x]),0\n"
" vleg  30,240(1,%[x]),1\n"
" vleg  31,248(1,%[x]),1\n"
" vflpdb  16,16\n"
" vflpdb  17,17\n"
" vflpdb  18,18\n"
" vflpdb  19,19\n"
" vflpdb  20,20\n"
" vflpdb  21,21\n"
" vflpdb  22,22\n"
" vflpdb  23,23\n"
" vflpdb  24,24\n"
" vflpdb  25,25\n"
" vflpdb  26,26\n"
" vflpdb  27,27\n"
" vflpdb  28,28\n"
" vflpdb  29,29\n"
" vflpdb  30,30\n"
" vflpdb  31,31\n"
" vfadb 16,16,17\n"
" vfadb 18,18,19\n"
" vfadb 20,20,21\n"
" vfadb 22,22,23\n"
" vfadb 24,24,25\n"
" vfadb 26,26,27\n"
" vfadb 28,28,29\n"
" vfadb 30,30,31\n"
" vfmaxdb  16,16,24,0\n"
" vfmaxdb  18,18,26,0\n"
" vfmaxdb  20,20,28,0\n"
" vfmaxdb  22,22,30,0\n"
" vfmaxdb  16,16,20,0\n"
" vfmaxdb  18,18,22,0\n"
" vfmaxdb  16,16,18,0\n"
" vfmaxdb  0,0,16,0\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" vrepg  16,0,1\n"
" wfmaxdb 0,0,16,0\n"
" ldr    %[amax],0\n"
    : [amax] "=f"(amax),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");
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

    BLASLONG n1 = n & -16;
    if (n1 > 0) {

      maxf = zamax_kernel_16(n1, x);
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
