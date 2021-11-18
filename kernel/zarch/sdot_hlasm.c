/***************************************************************************
Copyright (c) 2013-2019,The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms,with or without
modification,are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice,this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice,this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES,INCLUDING,BUT NOT LIMITED TO,THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT,INDIRECT,INCIDENTAL,SPECIAL,EXEMPLARY,OR CONSEQUENTIAL
DAMAGES (INCLUDING,BUT NOT LIMITED TO,PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE,DATA,OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY,WHETHER IN CONTRACT,STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE,EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

static FLOAT sdot_kernel_32(BLASLONG n, FLOAT *x, FLOAT *y) {
  FLOAT dot;
//#if 0
__asm__ (
" vzero 0\n"
" vzero 1\n"
" vzero 2\n"
" vzero 3\n"
" vzero 4\n"
" vzero 5\n"
" vzero 6\n"
" vzero 7\n"
" srlg  %[n],%[n],5\n"
" xgr   1,1\n"
"label0 pfd 1,1024(1,%[x])\n"
" pfd 1,1024(1,%[y])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,64(1,%[x])\n"
" vl  21,80(1,%[x])\n"
" vl  22,96(1,%[x])\n"
" vl  23,112(1,%[x])\n"
" vl  24,0(1,%[y])\n"
" vl  25,16(1,%[y])\n"
" vl  26,32(1,%[y])\n"
" vl  27,48(1,%[y])\n"
" vl  28,64(1,%[y])\n"
" vl  29,80(1,%[y])\n"
" vl  30,96(1,%[y])\n"
" vl  31,112(1,%[y])\n"
" vfmasb   0,16,24,0\n"
" vfmasb   1,17,25,1\n"
" vfmasb   2,18,26,2\n"
" vfmasb   3,19,27,3\n"
" vfmasb   4,20,28,4\n"
" vfmasb   5,21,29,5\n"
" vfmasb   6,22,30,6\n"
" vfmasb   7,23,31,7\n"
" agfi   1,128\n"
" brctg  %[n],label0\n"
" vfasb   0,0,1\n"
" vfasb   0,0,2\n"
" vfasb   0,0,3\n"
" vfasb   0,0,4\n"
" vfasb   0,0,5\n"
" vfasb   0,0,6\n"
" vfasb   0,0,7\n"
" vrepf  1,0,1\n"
" vrepf  2,0,2\n"
" vrepf  3,0,3\n"
" aebr   0,1\n"
" aebr   0,2\n"
" aebr   0,3\n"
" ler    %[dot],0\n"
    : [dot] "=f"(dot),[n] "+&r"(n)
    : "m"(*(const FLOAT (*)[n]) x),[x] "a"(x),
       "m"(*(const FLOAT (*)[n]) y),[y] "a"(y)
    : "cc", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v16",
       "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26",
       "v27", "v28", "v29", "v30", "v31");
//#endif
  return dot;
}

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;

  FLOAT dot = 0.0;

  if (n <= 0)
    return (dot);

  if ((inc_x == 1) && (inc_y == 1)) {

    BLASLONG n1 = n & -32;

    if (n1)
      dot = sdot_kernel_32(n1, x, y);

    i = n1;
    while (i < n) {

      dot += y[i] * x[i];
      i++;

    }
    return (dot);

  }

  BLASLONG n1 = n & -2;

  while (i < n1) {

    dot += y[iy] * x[ix] + y[iy + inc_y] * x[ix + inc_x];
    ix += inc_x * 2;
    iy += inc_y * 2;
    i += 2;

  }

  while (i < n) {

    dot += y[iy] * x[ix];
    ix += inc_x;
    iy += inc_y;
    i++;

  }
  return (dot);

}
