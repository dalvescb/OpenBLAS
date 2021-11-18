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

static FLOAT damin_kernel_32(BLASLONG n, FLOAT *x) {
  FLOAT amin;
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
" vl  24,128(1,%[x])\n"
" vl  25,144(1,%[x])\n"
" vl  26,160(1,%[x])\n"
" vl  27,176(1,%[x])\n"
" vl  28,192(1,%[x])\n"
" vl  29,208(1,%[x])\n"
" vl  30,224(1,%[x])\n"
" vl  31,240(1,%[x])\n"
" vfmindb  16,16,24,8\n"
" vfmindb  17,17,25,8\n"
" vfmindb  18,18,26,8\n"
" vfmindb  19,19,27,8\n"
" vfmindb  20,20,28,8\n"
" vfmindb  21,21,29,8\n"
" vfmindb  22,22,30,8\n"
" vfmindb  23,23,31,8\n"
" vfmindb  16,16,20,8\n"
" vfmindb  17,17,21,8\n"
" vfmindb  18,18,22,8\n"
" vfmindb  19,19,23,8\n"
" vfmindb  16,16,18,8\n"
" vfmindb  17,17,19,8\n"
" vfmindb  16,16,17,8\n"
" vfmindb  0,0,16,8\n"
" agfi    1,256\n"
" brctg   %[n],label0\n"
" vrepg   16,0,1\n"
" wfmindb 0,0,16,8\n"
" lpdr    %[amin],0\n"
    : [amin] "=f"(amin),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");
  return amin;
}

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG j = 0;
  FLOAT minf = 0.0;

  if (n <= 0 || inc_x <= 0)
    return (minf);

  if (inc_x == 1) {

    BLASLONG n1 = n & -32;
    if (n1 > 0) {

      minf = damin_kernel_32(n1, x);

      i = n1;
    } else {
      minf = ABS(x[0]);
      i++;
    }

    while (i < n) {
      if (ABS(x[i]) < minf) {
        minf = ABS(x[i]);
      }
      i++;
    }
    return (minf);

  } else {

    minf = ABS(x[0]);

    BLASLONG n1 = n & -4;
    while (j < n1) {

      if (ABS(x[i]) < minf) {
        minf = ABS(x[i]);
      }
      if (ABS(x[i + inc_x]) < minf) {
        minf = ABS(x[i + inc_x]);
      }
      if (ABS(x[i + 2 * inc_x]) < minf) {
        minf = ABS(x[i + 2 * inc_x]);
      }
      if (ABS(x[i + 3 * inc_x]) < minf) {
        minf = ABS(x[i + 3 * inc_x]);
      }

      i += inc_x * 4;

      j += 4;

    }

    while (j < n) {
      if (ABS(x[i]) < minf) {
        minf = ABS(x[i]);
      }
      i += inc_x;
      j++;
    }
    return (minf);
  }
}
