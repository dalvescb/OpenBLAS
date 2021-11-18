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

#define NBMAX 2048

static void sgemv_kernel_4x4(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
  register FLOAT *ap2 = ap[2];
  register FLOAT *ap3 = ap[3];
__asm__ (
" vzero 0\n"
" vzero 1\n"
" vzero 2\n"
" vzero 3\n"
" vzero 4\n"
" vzero 5\n"
" vzero 6\n"
" vzero 7\n"
" xgr   1,1\n"
" lghi    0,-32\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelone%=\n"
" srlg  0,0,5\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n"
" pfd 1,1024(1,%[ap1])\n"
" pfd 1,1024(1,%[ap2])\n"
" pfd 1,1024(1,%[ap3])\n"
" pfd 1,1024(1,%[x])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,64(1,%[x])\n"
" vl  21,80(1,%[x])\n"
" vl  22,96(1,%[x])\n"
" vl  23,112(1,%[x])\n"
" vl  24,0(1,%[ap0])\n"
" vfmasb   0,16,24,0\n"
" vl  25,0(1,%[ap1])\n"
" vfmasb   1,16,25,1\n"
" vl  26,0(1,%[ap2])\n"
" vfmasb   2,16,26,2\n"
" vl  27,0(1,%[ap3])\n"
" vfmasb   3,16,27,3\n"
" vl  28,16(1,%[ap0])\n"
" vfmasb   4,17,28,4\n"
" vl  29,16(1,%[ap1])\n"
" vfmasb   5,17,29,5\n"
" vl  30,16(1,%[ap2])\n"
" vfmasb   6,17,30,6\n"
" vl  31,16(1,%[ap3])\n"
" vfmasb   7,17,31,7\n"
" vl  24,32(1,%[ap0])\n"
" vfmasb   0,18,24,0\n"
" vl  25,32(1,%[ap1])\n"
" vfmasb   1,18,25,1\n"
" vl  26,32(1,%[ap2])\n"
" vfmasb   2,18,26,2\n"
" vl  27,32(1,%[ap3])\n"
" vfmasb   3,18,27,3\n"
" vl  28,48(1,%[ap0])\n"
" vfmasb   4,19,28,4\n"
" vl  29,48(1,%[ap1])\n"
" vfmasb   5,19,29,5\n"
" vl  30,48(1,%[ap2])\n"
" vfmasb   6,19,30,6\n"
" vl  31,48(1,%[ap3])\n"
" vfmasb   7,19,31,7\n"
" vl  24,64(1,%[ap0])\n"
" vfmasb   0,20,24,0\n"
" vl  25,64(1,%[ap1])\n"
" vfmasb   1,20,25,1\n"
" vl  26,64(1,%[ap2])\n"
" vfmasb   2,20,26,2\n"
" vl  27,64(1,%[ap3])\n"
" vfmasb   3,20,27,3\n"
" vl  28,80(1,%[ap0])\n"
" vfmasb   4,21,28,4\n"
" vl  29,80(1,%[ap1])\n"
" vfmasb   5,21,29,5\n"
" vl  30,80(1,%[ap2])\n"
" vfmasb   6,21,30,6\n"
" vl  31,80(1,%[ap3])\n"
" vfmasb   7,21,31,7\n"
" vl  24,96(1,%[ap0])\n"
" vfmasb   0,22,24,0\n"
" vl  25,96(1,%[ap1])\n"
" vfmasb   1,22,25,1\n"
" vl  26,96(1,%[ap2])\n"
" vfmasb   2,22,26,2\n"
" vl  27,96(1,%[ap3])\n"
" vfmasb   3,22,27,3\n"
" vl  28,112(1,%[ap0])\n"
" vfmasb   4,23,28,4\n"
" vl  29,112(1,%[ap1])\n"
" vfmasb   5,23,29,5\n"
" vl  30,112(1,%[ap2])\n"
" vfmasb   6,23,30,6\n"
" vl  31,112(1,%[ap3])\n"
" vfmasb   7,23,31,7\n"
" agfi   1,128\n"
" brctg  0,labelzero%=\n"
"labelone%= lghi    0,28\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelthree%=\n"
" srlg  0,0,2\n"
"labeltwo%= vl  16,0(1,%[x])\n"
" vl  24,0(1,%[ap0])\n"
" vfmasb   0,16,24,0\n"
" vl  25,0(1,%[ap1])\n"
" vfmasb   1,16,25,1\n"
" vl  26,0(1,%[ap2])\n"
" vfmasb   2,16,26,2\n"
" vl  27,0(1,%[ap3])\n"
" vfmasb   3,16,27,3\n"
" agfi   1,16\n"
" brctg  0,labeltwo%=\n"
"labelthree%= vfasb  0,0,4\n"
" vfasb  1,1,5\n"
" vfasb  2,2,6\n"
" vfasb  3,3,7\n"
" veslg  4,0,32\n"
" vfasb  0,0,4\n"
" vrepg  4,0,1\n"
" aebr   0,4\n"
" ste    0,0(%[y])\n"
" veslg  4,1,32\n"
" vfasb  1,1,4\n"
" vrepg  4,1,1\n"
" aebr   1,4\n"
" ste    1,4(%[y])\n"
" veslg  4,2,32\n"
" vfasb  2,2,4\n"
" vrepg  4,2,1\n"
" aebr   2,4\n"
" ste    2,8(%[y])\n"
" veslg  4,3,32\n"
" vfasb  3,3,4\n"
" vrepg  4,3,1\n"
" aebr   3,4\n"
" ste    3,12(%[y])\n"
    : "=m"(*(FLOAT (*)[4]) y)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n]) ap0),[ap0] "a"(ap0),
       "m"(*(const FLOAT (*)[n]) ap1),[ap1] "a"(ap1),
       "m"(*(const FLOAT (*)[n]) ap2),[ap2] "a"(ap2),
       "m"(*(const FLOAT (*)[n]) ap3),[ap3] "a"(ap3),
       "m"(*(const FLOAT (*)[n]) x),[x] "a"(x),[n] "r"(n)
    : "cc", "r0", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7",
       "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25",
       "v26", "v27", "v28", "v29", "v30", "v31");
}

static void sgemv_kernel_4x2(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
__asm__ (
" vzero 0\n"
" vzero 1\n"
" vzero 2\n"
" vzero 3\n"
" vzero 4\n"
" vzero 5\n"
" vzero 6\n"
" vzero 7\n"
" xgr   1,1\n"
" lghi    0,-32\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelone%=\n"
" srlg  0,0,5\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n"
" pfd 1,1024(1,%[ap1])\n"
" pfd 1,1024(1,%[x])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,64(1,%[x])\n"
" vl  21,80(1,%[x])\n"
" vl  22,96(1,%[x])\n"
" vl  23,112(1,%[x])\n"
" vl  24,0(1,%[ap0])\n"
" vfmasb   0,16,24,0\n"
" vl  25,0(1,%[ap1])\n"
" vfmasb   1,16,25,1\n"
" vl  26,16(1,%[ap0])\n"
" vfmasb   2,17,26,2\n"
" vl  27,16(1,%[ap1])\n"
" vfmasb   3,17,27,3\n"
" vl  28,32(1,%[ap0])\n"
" vfmasb   4,18,28,4\n"
" vl  29,32(1,%[ap1])\n"
" vfmasb   5,18,29,5\n"
" vl  30,48(1,%[ap0])\n"
" vfmasb   6,19,30,6\n"
" vl  31,48(1,%[ap1])\n"
" vfmasb   7,19,31,7\n"
" vl  24,64(1,%[ap0])\n"
" vfmasb   0,20,24,0\n"
" vl  25,64(1,%[ap1])\n"
" vfmasb   1,20,25,1\n"
" vl  26,80(1,%[ap0])\n"
" vfmasb   2,21,26,2\n"
" vl  27,80(1,%[ap1])\n"
" vfmasb   3,21,27,3\n"
" vl  28,96(1,%[ap0])\n"
" vfmasb   4,22,28,4\n"
" vl  29,96(1,%[ap1])\n"
" vfmasb   5,22,29,5\n"
" vl  30,112(1,%[ap0])\n"
" vfmasb   6,23,30,6\n"
" vl  31,112(1,%[ap1])\n"
" vfmasb   7,23,31,7\n"
" agfi   1,128\n"
" brctg  0,labelzero%=\n"
"labelone%= lghi    0,28\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelthree%=\n"
" srlg  0,0,2\n"
"labeltwo%= vl  16,0(1,%[x])\n"
" vl  24,0(1,%[ap0])\n"
" vfmasb   0,16,24,0\n"
" vl  25,0(1,%[ap1])\n"
" vfmasb   1,16,25,1\n"
" agfi   1,16\n"
" brctg  0,labeltwo%=\n"
"labelthree%= vfasb  0,0,2\n"
" vfasb  0,0,4\n"
" vfasb  0,0,6\n"
" vfasb  1,1,3\n"
" vfasb  1,1,5\n"
" vfasb  1,1,7\n"
" veslg  2,0,32\n"
" vfasb  0,0,2\n"
" vrepg  2,0,1\n"
" aebr   0,2\n"
" ste    0,0(%[y])\n"
" veslg  2,1,32\n"
" vfasb  1,1,2\n"
" vrepg  2,1,1\n"
" aebr   1,2\n"
" ste    1,4(%[y])\n"
    : "=m"(*(FLOAT (*)[2]) y)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n]) ap0),[ap0] "a"(ap0),
       "m"(*(const FLOAT (*)[n]) ap1),[ap1] "a"(ap1),
       "m"(*(const FLOAT (*)[n]) x),[x] "a"(x),[n] "r"(n)
    : "cc", "r0", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7",
       "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25",
       "v26", "v27", "v28", "v29", "v30", "v31");
}

static void sgemv_kernel_4x1(BLASLONG n, FLOAT *a0, FLOAT *x, FLOAT *y) {
__asm__ (
" vzero 0\n"
" vzero 1\n"
" vzero 2\n"
" vzero 3\n"
" vzero 4\n"
" vzero 5\n"
" vzero 6\n"
" vzero 7\n"
" xgr   1,1\n"
" lghi    0,-32\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelone%=\n"
" srlg  0,0,5\n"
"labelzero%= pfd 1,1024(1,%[a0])\n"
" pfd 1,1024(1,%[x])\n"
" vl  16,0(1,%[x])\n"
" vl  17,16(1,%[x])\n"
" vl  18,32(1,%[x])\n"
" vl  19,48(1,%[x])\n"
" vl  20,64(1,%[x])\n"
" vl  21,80(1,%[x])\n"
" vl  22,96(1,%[x])\n"
" vl  23,112(1,%[x])\n"
" vl  24,0(1,%[a0])\n"
" vfmasb   0,16,24,0\n"
" vl  25,16(1,%[a0])\n"
" vfmasb   1,17,25,1\n"
" vl  26,32(1,%[a0])\n"
" vfmasb   2,18,26,2\n"
" vl  27,48(1,%[a0])\n"
" vfmasb   3,19,27,3\n"
" vl  28,64(1,%[a0])\n"
" vfmasb   4,20,28,4\n"
" vl  29,80(1,%[a0])\n"
" vfmasb   5,21,29,5\n"
" vl  30,96(1,%[a0])\n"
" vfmasb   6,22,30,6\n"
" vl  31,112(1,%[a0])\n"
" vfmasb   7,23,31,7\n"
" agfi   1,128\n"
" brctg  0,labelzero%=\n"
"labelone%= lghi    0,28\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelthree%=\n"
" srlg  0,0,2\n"
"labeltwo%= vl  16,0(1,%[x])\n"
" vl  24,0(1,%[a0])\n"
" vfmasb   0,16,24,0\n"
" agfi   1,16\n"
" brctg  0,labeltwo%=\n"
"labelthree%= vfasb  0,0,1\n"
" vfasb  0,0,2\n"
" vfasb  0,0,3\n"
" vfasb  0,0,4\n"
" vfasb  0,0,5\n"
" vfasb  0,0,6\n"
" vfasb  0,0,7\n"
" veslg  1,0,32\n"
" vfasb  0,0,1\n"
" vrepg  1,0,1\n"
" aebr   0,1\n"
" ste    0,0(%[y])\n"
    : "=m"(*(FLOAT (*)[1]) y)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n]) a0),[a0] "a"(a0),
       "m"(*(const FLOAT (*)[n]) x),[x] "a"(x),[n] "r"(n)
    : "cc", "r0", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7",
       "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25",
       "v26", "v27", "v28", "v29", "v30", "v31");
}

static void copy_x(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src) {
  BLASLONG i;
  for (i = 0; i < n; i++) {
    dest[i] = *src;
    src += inc_src;
  }
}

static void add_y_kernel_4(BLASLONG n, FLOAT da, FLOAT *src, FLOAT *dest) {
__asm__ (
" vlrepf 0,%[da]\n"
" xgr   1,1\n"
" lghi    0,-32\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelone%=\n"
" srlg  0,0,5\n"
"labelzero%= pfd 1,1024(1,%[src])\n"
" pfd 2,1024(1,%[dest])\n"
" vl  16,0(1,%[src])\n"
" vl  17,16(1,%[src])\n"
" vl  18,32(1,%[src])\n"
" vl  19,48(1,%[src])\n"
" vl  20,64(1,%[src])\n"
" vl  21,80(1,%[src])\n"
" vl  22,96(1,%[src])\n"
" vl  23,112(1,%[src])\n"
" vl  24,0(1,%[dest])\n"
" vfmasb   24,16,0,24\n"
" vst  24,0(1,%[dest])\n"
" vl  25,16(1,%[dest])\n"
" vfmasb   25,17,0,25\n"
" vst  25,16(1,%[dest])\n"
" vl  26,32(1,%[dest])\n"
" vfmasb   26,18,0,26\n"
" vst  26,32(1,%[dest])\n"
" vl  27,48(1,%[dest])\n"
" vfmasb   27,19,0,27\n"
" vst  27,48(1,%[dest])\n"
" vl  28,64(1,%[dest])\n"
" vfmasb   28,20,0,28\n"
" vst  28,64(1,%[dest])\n"
" vl  29,80(1,%[dest])\n"
" vfmasb   29,21,0,29\n"
" vst  29,80(1,%[dest])\n"
" vl  30,96(1,%[dest])\n"
" vfmasb   30,22,0,30\n"
" vst  30,96(1,%[dest])\n"
" vl  31,112(1,%[dest])\n"
" vfmasb   31,23,0,31\n"
" vst  31,112(1,%[dest])\n"
" agfi   1,128\n"
" brctg  0,labelzero%=\n"
"labelone%= lghi    0,28\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelthree%=\n"
" srlg  0,0,2\n"
"labeltwo%= vl  16,0(1,%[src])\n"
" vl  24,0(1,%[dest])\n"
" vfmasb   24,16,0,24\n"
" vst  24,0(1,%[dest])\n"
" agfi   1,16\n"
" brctg  0,labeltwo%=\n"
"labelthree%= nop 0\n"
    : "+m"(*(FLOAT (*)[n]) dest)
    : [dest] "a"(dest),[da] "Q"(da), "m"(*(const FLOAT (*)[n]) src),
       [src] "a"(src),[n] "r"(n)
    : "cc", "r0", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}
static void add_y(BLASLONG n, FLOAT da, FLOAT *src, FLOAT *dest,
                  BLASLONG inc_dest) {
  if (inc_dest == 1)
    add_y_kernel_4(n, da, src, dest);
  else {
    BLASLONG i;
    for (i = 0; i < n; i++) {
      *dest += src[i] * da;
      dest += inc_dest;
    }
  }
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a,
          BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y,
          FLOAT *buffer) {
  BLASLONG register i;
  BLASLONG register j;
  FLOAT *a_ptr;
  FLOAT *x_ptr;
  FLOAT *y_ptr;
  BLASLONG n0;
  BLASLONG n1;
  BLASLONG m1;
  BLASLONG m2;
  BLASLONG m3;
  BLASLONG n2;
  FLOAT ybuffer[2] __attribute__ ((aligned(16)));
  FLOAT *xbuffer;
  FLOAT *ytemp;

  if (m < 1)
    return (0);
  if (n < 1)
    return (0);

  xbuffer = buffer;
  ytemp = buffer + (m < NBMAX ? m : NBMAX);

  n0 = n / NBMAX;
  n1 = (n % NBMAX) >> 2;
  n2 = n & 3;

  m3 = m & 3;
  m1 = m & -4;
  m2 = (m & (NBMAX - 1)) - m3;

  BLASLONG NB = NBMAX;

  while (NB == NBMAX) {
    m1 -= NB;
    if (m1 < 0) {
      if (m2 == 0)
        break;
      NB = m2;
    }

    y_ptr = y;
    a_ptr = a;
    x_ptr = x;

    if (inc_x == 1)
      xbuffer = x_ptr;
    else
      copy_x(NB, x_ptr, xbuffer, inc_x);

    FLOAT *ap[4];
    FLOAT *yp;
    BLASLONG register lda4 = 4 * lda;
    ap[0] = a_ptr;
    ap[1] = a_ptr + lda;
    ap[2] = ap[1] + lda;
    ap[3] = ap[2] + lda;

    if (n0 > 0) {
      BLASLONG nb1 = NBMAX / 4;
      for (j = 0; j < n0; j++) {

        yp = ytemp;
        for (i = 0; i < nb1; i++) {
          sgemv_kernel_4x4(NB, ap, xbuffer, yp);
          ap[0] += lda4;
          ap[1] += lda4;
          ap[2] += lda4;
          ap[3] += lda4;
          yp += 4;
        }
        add_y(nb1 * 4, alpha, ytemp, y_ptr, inc_y);
        y_ptr += nb1 * inc_y * 4;
        a_ptr += nb1 * lda4;

      }

    }

    yp = ytemp;

    for (i = 0; i < n1; i++) {
      sgemv_kernel_4x4(NB, ap, xbuffer, yp);
      ap[0] += lda4;
      ap[1] += lda4;
      ap[2] += lda4;
      ap[3] += lda4;
      yp += 4;
    }
    if (n1 > 0) {
      add_y(n1 * 4, alpha, ytemp, y_ptr, inc_y);
      y_ptr += n1 * inc_y * 4;
      a_ptr += n1 * lda4;
    }

    if (n2 & 2) {

      sgemv_kernel_4x2(NB, ap, xbuffer, ybuffer);
      a_ptr += lda * 2;
      *y_ptr += ybuffer[0] * alpha;
      y_ptr += inc_y;
      *y_ptr += ybuffer[1] * alpha;
      y_ptr += inc_y;

    }

    if (n2 & 1) {

      sgemv_kernel_4x1(NB, a_ptr, xbuffer, ybuffer);
      // a_ptr  += lda;
      *y_ptr += ybuffer[0] * alpha;
      // y_ptr  += inc_y;

    }
    a += NB;
    x += NB * inc_x;
  }

  if (m3 == 0)
    return (0);

  x_ptr = x;
  a_ptr = a;
  if (m3 == 3) {
    FLOAT xtemp0 = *x_ptr * alpha;
    x_ptr += inc_x;
    FLOAT xtemp1 = *x_ptr * alpha;
    x_ptr += inc_x;
    FLOAT xtemp2 = *x_ptr * alpha;

    FLOAT *aj = a_ptr;
    y_ptr = y;

    if (lda == 3 && inc_y == 1) {

      for (j = 0; j < (n & -4); j += 4) {

        y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1 + aj[2] * xtemp2;
        y_ptr[j + 1] += aj[3] * xtemp0 + aj[4] * xtemp1 + aj[5] * xtemp2;
        y_ptr[j + 2] += aj[6] * xtemp0 + aj[7] * xtemp1 + aj[8] * xtemp2;
        y_ptr[j + 3] += aj[9] * xtemp0 + aj[10] * xtemp1 + aj[11] * xtemp2;
        aj += 12;
      }

      for (; j < n; j++) {
        y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1 + aj[2] * xtemp2;
        aj += 3;
      }

    } else {

      if (inc_y == 1) {

        BLASLONG register lda2 = lda << 1;
        BLASLONG register lda4 = lda << 2;
        BLASLONG register lda3 = lda2 + lda;

        for (j = 0; j < (n & -4); j += 4) {

          y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1 + *(aj + 2) * xtemp2;
          y_ptr[j + 1] +=
            *(aj + lda) * xtemp0 + *(aj + lda + 1) * xtemp1 + *(aj + lda +
                                                                2) * xtemp2;
          y_ptr[j + 2] +=
            *(aj + lda2) * xtemp0 + *(aj + lda2 + 1) * xtemp1 + *(aj + lda2 +
                                                                  2) * xtemp2;
          y_ptr[j + 3] +=
            *(aj + lda3) * xtemp0 + *(aj + lda3 + 1) * xtemp1 + *(aj + lda3 +
                                                                  2) * xtemp2;
          aj += lda4;
        }

        for (; j < n; j++) {

          y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1 + *(aj + 2) * xtemp2;
          aj += lda;
        }

      } else {

        for (j = 0; j < n; j++) {
          *y_ptr += *aj * xtemp0 + *(aj + 1) * xtemp1 + *(aj + 2) * xtemp2;
          y_ptr += inc_y;
          aj += lda;
        }

      }

    }
    return (0);
  }

  if (m3 == 2) {
    FLOAT xtemp0 = *x_ptr * alpha;
    x_ptr += inc_x;
    FLOAT xtemp1 = *x_ptr * alpha;

    FLOAT *aj = a_ptr;
    y_ptr = y;

    if (lda == 2 && inc_y == 1) {

      for (j = 0; j < (n & -4); j += 4) {
        y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1;
        y_ptr[j + 1] += aj[2] * xtemp0 + aj[3] * xtemp1;
        y_ptr[j + 2] += aj[4] * xtemp0 + aj[5] * xtemp1;
        y_ptr[j + 3] += aj[6] * xtemp0 + aj[7] * xtemp1;
        aj += 8;

      }

      for (; j < n; j++) {
        y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1;
        aj += 2;
      }

    } else {
      if (inc_y == 1) {

        BLASLONG register lda2 = lda << 1;
        BLASLONG register lda4 = lda << 2;
        BLASLONG register lda3 = lda2 + lda;

        for (j = 0; j < (n & -4); j += 4) {

          y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1;
          y_ptr[j + 1] += *(aj + lda) * xtemp0 + *(aj + lda + 1) * xtemp1;
          y_ptr[j + 2] += *(aj + lda2) * xtemp0 + *(aj + lda2 + 1) * xtemp1;
          y_ptr[j + 3] += *(aj + lda3) * xtemp0 + *(aj + lda3 + 1) * xtemp1;
          aj += lda4;
        }

        for (; j < n; j++) {

          y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1;
          aj += lda;
        }

      } else {
        for (j = 0; j < n; j++) {
          *y_ptr += *aj * xtemp0 + *(aj + 1) * xtemp1;
          y_ptr += inc_y;
          aj += lda;
        }
      }

    }
    return (0);

  }

  FLOAT xtemp = *x_ptr * alpha;
  FLOAT *aj = a_ptr;
  y_ptr = y;
  if (lda == 1 && inc_y == 1) {
    for (j = 0; j < (n & -4); j += 4) {
      y_ptr[j] += aj[j] * xtemp;
      y_ptr[j + 1] += aj[j + 1] * xtemp;
      y_ptr[j + 2] += aj[j + 2] * xtemp;
      y_ptr[j + 3] += aj[j + 3] * xtemp;
    }
    for (; j < n; j++) {
      y_ptr[j] += aj[j] * xtemp;
    }

  } else {
    if (inc_y == 1) {

      BLASLONG register lda2 = lda << 1;
      BLASLONG register lda4 = lda << 2;
      BLASLONG register lda3 = lda2 + lda;
      for (j = 0; j < (n & -4); j += 4) {
        y_ptr[j] += *aj * xtemp;
        y_ptr[j + 1] += *(aj + lda) * xtemp;
        y_ptr[j + 2] += *(aj + lda2) * xtemp;
        y_ptr[j + 3] += *(aj + lda3) * xtemp;
        aj += lda4;
      }

      for (; j < n; j++) {
        y_ptr[j] += *aj * xtemp;
        aj += lda;
      }

    } else {
      for (j = 0; j < n; j++) {
        *y_ptr += *aj * xtemp;
        y_ptr += inc_y;
        aj += lda;
      }

    }
  }

  return (0);
}
