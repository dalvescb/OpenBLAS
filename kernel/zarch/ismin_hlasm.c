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

static BLASLONG ismin_kernel_64(BLASLONG n, FLOAT *x, FLOAT *min) {
  BLASLONG imin;
__asm__ (
" vl     0,0(%[x])\n"
" vleig  1,0,0\n"
" vleig  1,2,1\n"
" vleig  2,1,0\n"
" vleig  2,3,1\n"
" vrepig 3,32\n"
" vzero  4\n"
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
" vleif  28,16,0\n"
" vleif  28,17,1\n"
" vleif  28,18,2\n"
" vleif  28,19,3\n"
" vleif  29,20,0\n"
" vleif  29,21,1\n"
" vleif  29,22,2\n"
" vleif  29,23,3\n"
" vleif  30,24,0\n"
" vleif  30,25,1\n"
" vleif  30,26,2\n"
" vleif  30,27,3\n"
" vleif  31,28,0\n"
" vleif  31,29,1\n"
" vleif  31,30,2\n"
" vleif  31,31,3\n"
" srlg  %[n],%[n],6\n"
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
" vfchesb  5,17,16\n"
" vfchesb  6,19,18\n"
" vfchesb  7,21,20\n"
" vfchesb  8,23,22\n"
" vsel    16,16,17,5\n"
" vsel    5,24,25,5\n"
" vsel    17,18,19,6\n"
" vsel    6,26,27,6\n"
" vsel    18,20,21,7\n"
" vsel    7,28,29,7\n"
" vsel    19,22,23,8\n"
" vsel    8,30,31,8\n"
" vfchesb  20,17,16\n"
" vfchesb  21,19,18\n"
" vsel    16,16,17,20\n"
" vsel    5,5,6,20\n"
" vsel    17,18,19,21\n"
" vsel    6,7,8,21\n"
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
" vl  16,128(1,%[x])\n"
" vl  17,144(1,%[x])\n"
" vl  18,160(1,%[x])\n"
" vl  19,176(1,%[x])\n"
" vl  20,192(1,%[x])\n"
" vl  21,208(1,%[x])\n"
" vl  22,224(1,%[x])\n"
" vl  23,240(1,%[x])\n"
" vfchesb  5,17,16\n"
" vfchesb  6,19,18\n"
" vfchesb  7,21,20\n"
" vfchesb  8,23,22\n"
" vsel    16,16,17,5\n"
" vsel    5,24,25,5\n"
" vsel    17,18,19,6\n"
" vsel    6,26,27,6\n"
" vsel    18,20,21,7\n"
" vsel    7,28,29,7\n"
" vsel    19,22,23,8\n"
" vsel    8,30,31,8\n"
" vfchesb  20,17,16\n"
" vfchesb  21,19,18\n"
" vsel    16,16,17,20\n"
" vsel    5,5,6,20\n"
" vsel    17,18,19,21\n"
" vsel    6,7,8,21\n"
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
" vstef  0,%[min],0\n"
" vmnlg  0,1,3\n"
" vlgvg  %[imin],0,0\n"
" j label2\n"
"label1 wfchsb 4,0,2\n"
" vesrlg 4,4,32\n"
" vsegf  4,4\n"
" vsel   1,3,1,4\n"
" vsel   0,2,0,4\n"
" ste    0,%[min]\n"
" vlgvg  %[imin],1,0\n"
"label2 nop 0\n"
    : [imin] "=r"(imin),[min] "=Q"(*min),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v2", "v4", "v5", "v6", "v7", "v8", "v16",
       "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26",
       "v27", "v28", "v29", "v30", "v31");
  return imin;
}

BLASLONG CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
  BLASLONG i = 0;
  BLASLONG j = 0;
  FLOAT minf = 0.0;
  BLASLONG min = 0;

  if (n <= 0 || inc_x <= 0)
    return (min);

  if (inc_x == 1) {

    BLASLONG n1 = n & -64;
    if (n1 > 0) {

      min = ismin_kernel_64(n1, x, &minf);

      i = n1;
    } else {
      minf = x[0];
      i++;
    }

    while (i < n) {
      if (x[i] < minf) {
        min = i;
        minf = x[i];
      }
      i++;
    }
    return (min + 1);

  } else {

    min = 0;
    minf = x[0];

    BLASLONG n1 = n & -4;
    while (j < n1) {

      if (x[i] < minf) {
        min = j;
        minf = x[i];
      }
      if (x[i + inc_x] < minf) {
        min = j + 1;
        minf = x[i + inc_x];
      }
      if (x[i + 2 * inc_x] < minf) {
        min = j + 2;
        minf = x[i + 2 * inc_x];
      }
      if (x[i + 3 * inc_x] < minf) {
        min = j + 3;
        minf = x[i + 3 * inc_x];
      }

      i += inc_x * 4;

      j += 4;

    }

    while (j < n) {
      if (x[i] < minf) {
        min = j;
        minf = x[i];
      }
      i += inc_x;
      j++;
    }
    return (min + 1);
  }
}
