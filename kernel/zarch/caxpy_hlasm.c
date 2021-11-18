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

static void caxpy_kernel_16(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *alpha) {
__asm__ (
#if !defined(CONJ)
" vlrepf 0,0(%[alpha])\n"
" vlef   1,4(%[alpha]),0\n"
" vlef   1,4(%[alpha]),2\n"
" vflcsb 1,1\n"
" vlef   1,4(%[alpha]),1\n"
" vlef   1,4(%[alpha]),3\n"
#else
" vlef   0,0(%[alpha]),1\n"
" vlef   0,0(%[alpha]),3\n"
" vflcsb 0,0\n"
" vlef   0,0(%[alpha]),0\n"
" vlef   0,0(%[alpha]),2\n"
" vlrepf 1,4(%[alpha])\n"
#endif
" srlg %[n],%[n],4\n"
" xgr  1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" pfd 2,1024(1,%[y])\n"
" vl   8,0(1,%[x])\n"
" vl   9,16(1,%[x])\n"
" vl   10,32(1,%[x])\n"
" vl   11,48(1,%[x])\n"
" vl   12,0(1,%[y])\n"
" vl   13,16(1,%[y])\n"
" vl   14,32(1,%[y])\n"
" vl   15,48(1,%[y])\n"
" vl   16,64(1,%[x])\n"
" vl   17,80(1,%[x])\n"
" vl   18,96(1,%[x])\n"
" vl   19,112(1,%[x])\n"
" vl   20,64(1,%[y])\n"
" vl   21,80(1,%[y])\n"
" vl   22,96(1,%[y])\n"
" vl   23,112(1,%[y])\n"
" verllg   24,8,32\n"
" verllg   25,9,32\n"
" verllg   26,10,32\n"
" verllg   27,11,32\n"
" verllg   28,16,32\n"
" verllg   29,17,32\n"
" verllg   30,18,32\n"
" verllg   31,19,32\n"
" vfmasb 8,8,0,12\n"
" vfmasb 9,9,0,13\n"
" vfmasb 10,10,0,14\n"
" vfmasb 11,11,0,15\n"
" vfmasb 16,16,0,20\n"
" vfmasb 17,17,0,21\n"
" vfmasb 18,18,0,22\n"
" vfmasb 19,19,0,23\n"
" vfmasb 8,24,1,8\n"
" vfmasb 9,25,1,9\n"
" vfmasb 10,26,1,10\n"
" vfmasb 11,27,1,11\n"
" vfmasb 16,28,1,16\n"
" vfmasb 17,29,1,17\n"
" vfmasb 18,30,1,18\n"
" vfmasb 19,31,1,19\n"
" vst 8,0(1,%[y])\n"
" vst 9,16(1,%[y])\n"
" vst 10,32(1,%[y])\n"
" vst 11,48(1,%[y])\n"
" vst 16,64(1,%[y])\n"
" vst 17,80(1,%[y])\n"
" vst 18,96(1,%[y])\n"
" vst 19,112(1,%[y])\n"
" agfi  1,128\n"
" brctg %[n],label0\n"
    : "+m"(*(FLOAT (*)[n * 2]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x),
       "m"(*(const FLOAT (*)[2]) alpha),[alpha] "a"(alpha)
    : "cc", "r1", "v0", "v1", "v8", "v9", "v10", "v11", "v12", "v13",
       "v14", "v15", "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23",
       "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");
}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da_r, FLOAT da_i,
          FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy,
          BLASLONG dummy2) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;
  FLOAT da[2] __attribute__ ((aligned(16)));

  if (n <= 0)
    return (0);

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -16;

    if (n1) {
      da[0] = da_r;
      da[1] = da_i;
      caxpy_kernel_16(n1, x, y, da);
      ix = 2 * n1;
    }
    i = n1;
    while (i < n) {
#if !defined(CONJ)
      y[ix] += (da_r * x[ix] - da_i * x[ix + 1]);
      y[ix + 1] += (da_r * x[ix + 1] + da_i * x[ix]);
#else
      y[ix] += (da_r * x[ix] + da_i * x[ix + 1]);
      y[ix + 1] -= (da_r * x[ix + 1] - da_i * x[ix]);
#endif
      i++;
      ix += 2;

    }
    return (0);

  }

  inc_x *= 2;
  inc_y *= 2;

  while (i < n) {

#if !defined(CONJ)
    y[iy] += (da_r * x[ix] - da_i * x[ix + 1]);
    y[iy + 1] += (da_r * x[ix + 1] + da_i * x[ix]);
#else
    y[iy] += (da_r * x[ix] + da_i * x[ix + 1]);
    y[iy + 1] -= (da_r * x[ix + 1] - da_i * x[ix]);
#endif
    ix += inc_x;
    iy += inc_y;
    i++;

  }
  return (0);

}
