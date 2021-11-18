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

static void dgemv_kernel_4x4(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y,
                             FLOAT *alpha) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
  register FLOAT *ap2 = ap[2];
  register FLOAT *ap3 = ap[3];
__asm__ (
" vlrepg 0,0(%[x])\n"
" vlrepg 1,8(%[x])\n"
" vlrepg 2,16(%[x])\n"
" vlrepg 3,24(%[x])\n"
" vlrepg 4,%[alpha]\n"
" vfmdb  0,0,4\n"
" vfmdb  1,1,4\n"
" vfmdb  2,2,4\n"
" vfmdb  3,3,4\n"
" xgr   1,1\n"
" lghi    0,-16\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelone%=\n" // Jump to label1
" srlg  0,0,4\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n" // label0
" pfd 1,1024(1,%[ap1])\n"
" pfd 1,1024(1,%[ap2])\n"
" pfd 1,1024(1,%[ap3])\n"
" pfd 2,1024(1,%[y])\n"
" vl  16,0(1,%[ap0])\n"
" vl  17,0(1,%[ap1])\n"
" vl  18,0(1,%[ap2])\n"
" vl  19,0(1,%[ap3])\n"
" vl  20,16(1,%[ap0])\n"
" vl  21,16(1,%[ap1])\n"
" vl  22,16(1,%[ap2])\n"
" vl  23,16(1,%[ap3])\n"
" vl  24,32(1,%[ap0])\n"
" vl  25,32(1,%[ap1])\n"
" vl  26,32(1,%[ap2])\n"
" vl  27,32(1,%[ap3])\n"
" vl  28,48(1,%[ap0])\n"
" vl  29,48(1,%[ap1])\n"
" vl  30,48(1,%[ap2])\n"
" vl  31,48(1,%[ap3])\n"
" vl  4,0(1,%[y])\n"
" vl  5,16(1,%[y])\n"
" vl  6,32(1,%[y])\n"
" vl  7,48(1,%[y])\n"
" vfmadb   4,16,0,4\n"
" vfmadb   5,20,0,5\n"
" vfmadb   6,24,0,6\n"
" vfmadb   7,28,0,7\n"
" vfmadb   4,17,1,4\n"
" vfmadb   5,21,1,5\n"
" vfmadb   6,25,1,6\n"
" vfmadb   7,29,1,7\n"
" vfmadb   4,18,2,4\n"
" vfmadb   5,22,2,5\n"
" vfmadb   6,26,2,6\n"
" vfmadb   7,30,2,7\n"
" vfmadb   4,19,3,4\n"
" vfmadb   5,23,3,5\n"
" vfmadb   6,27,3,6\n"
" vfmadb   7,31,3,7\n"
" vst 4,0(1,%[y])\n"
" vst 5,16(1,%[y])\n"
" vst 6,32(1,%[y])\n"
" vst 7,48(1,%[y])\n"
" vl  16,64(1,%[ap0])\n"
" vl  17,64(1,%[ap1])\n"
" vl  18,64(1,%[ap2])\n"
" vl  19,64(1,%[ap3])\n"
" vl  20,80(1,%[ap0])\n"
" vl  21,80(1,%[ap1])\n"
" vl  22,80(1,%[ap2])\n"
" vl  23,80(1,%[ap3])\n"
" vl  24,96(1,%[ap0])\n"
" vl  25,96(1,%[ap1])\n"
" vl  26,96(1,%[ap2])\n"
" vl  27,96(1,%[ap3])\n"
" vl  28,112(1,%[ap0])\n"
" vl  29,112(1,%[ap1])\n"
" vl  30,112(1,%[ap2])\n"
" vl  31,112(1,%[ap3])\n"
" vl  4,64(1,%[y])\n"
" vl  5,80(1,%[y])\n"
" vl  6,96(1,%[y])\n"
" vl  7,112(1,%[y])\n"
" vfmadb   4,16,0,4\n"
" vfmadb   5,20,0,5\n"
" vfmadb   6,24,0,6\n"
" vfmadb   7,28,0,7\n"
" vfmadb   4,17,1,4\n"
" vfmadb   5,21,1,5\n"
" vfmadb   6,25,1,6\n"
" vfmadb   7,29,1,7\n"
" vfmadb   4,18,2,4\n"
" vfmadb   5,22,2,5\n"
" vfmadb   6,26,2,6\n"
" vfmadb   7,30,2,7\n"
" vfmadb   4,19,3,4\n"
" vfmadb   5,23,3,5\n"
" vfmadb   6,27,3,6\n"
" vfmadb   7,31,3,7\n"
" vst 4,64(1,%[y])\n"
" vst 5,80(1,%[y])\n"
" vst 6,96(1,%[y])\n"
" vst 7,112(1,%[y])\n"
" agfi   1,128\n"
" brctg  0,labelzero%=\n" // Jump to label0
"labelone%= lghi    0,12\n" // label1
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelthree%=\n" // Jump to label3
" srlg  0,0,2\n"
"labeltwo%= vl  16,0(1,%[ap0])\n" // label2
" vl  17,0(1,%[ap1])\n"
" vl  18,0(1,%[ap2])\n"
" vl  19,0(1,%[ap3])\n"
" vl  20,16(1,%[ap0])\n"
" vl  21,16(1,%[ap1])\n"
" vl  22,16(1,%[ap2])\n"
" vl  23,16(1,%[ap3])\n"
" vl  4,0(1,%[y])\n"
" vl  5,16(1,%[y])\n"
" vfmadb   4,16,0,4\n"
" vfmadb   5,20,0,5\n"
" vfmadb   4,17,1,4\n"
" vfmadb   5,21,1,5\n"
" vfmadb   4,18,2,4\n"
" vfmadb   5,22,2,5\n"
" vfmadb   4,19,3,4\n"
" vfmadb   5,23,3,5\n"
" vst 4,0(1,%[y])\n"
" vst 5,16(1,%[y])\n"
" agfi   1,32\n"
" brctg  0,labeltwo%=\n" // Jump to label2
"labelthree%= nop 0\n" // label3
    : "+m"(*(FLOAT (*)[n]) y)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n]) ap0),[ap0] "a"(ap0),
       "m"(*(const FLOAT (*)[n]) ap1),[ap1] "a"(ap1),
       "m"(*(const FLOAT (*)[n]) ap2),[ap2] "a"(ap2),
       "m"(*(const FLOAT (*)[n]) ap3),[ap3] "a"(ap3),
       "m"(*(const FLOAT (*)[4]) x),[x] "a"(x),[alpha] "Q"(*alpha),
       [n] "r"(n)
    : "cc", "r0", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7",
       "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23", "v24", "v25",
       "v26", "v27", "v28", "v29", "v30", "v31");
}

static void dgemv_kernel_4x2(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y,
                             FLOAT *alpha) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
__asm__ (
" vlrepg 0,0(%[x])\n"
" vlrepg 1,8(%[x])\n"
" vlrepg 2,%[alpha]\n"
" vfmdb  0,0,2\n"
" vfmdb  1,1,2\n"
" xgr   1,1\n"
" lghi    0,-16\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelone%=\n"	// Jump	to label1
" srlg  0,0,4\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n" // label0
" pfd 1,1024(1,%[ap1])\n"
" pfd 2,1024(1,%[y])\n"
" vl  16,0(1,%[ap0])\n"
" vl  17,0(1,%[ap1])\n"
" vl  18,16(1,%[ap0])\n"
" vl  19,16(1,%[ap1])\n"
" vl  20,32(1,%[ap0])\n"
" vl  21,32(1,%[ap1])\n"
" vl  22,48(1,%[ap0])\n"
" vl  23,48(1,%[ap1])\n"
" vl  24,64(1,%[ap0])\n"
" vl  25,64(1,%[ap1])\n"
" vl  26,80(1,%[ap0])\n"
" vl  27,80(1,%[ap1])\n"
" vl  28,96(1,%[ap0])\n"
" vl  29,96(1,%[ap1])\n"
" vl  30,112(1,%[ap0])\n"
" vl  31,112(1,%[ap1])\n"
" vl  2,0(1,%[y])\n"
" vl  3,16(1,%[y])\n"
" vl  4,32(1,%[y])\n"
" vl  5,48(1,%[y])\n"
" vl  6,64(1,%[y])\n"
" vl  7,80(1,%[y])\n"
" vl  8,96(1,%[y])\n"
" vl  9,112(1,%[y])\n"
" vfmadb   2,16,0,2\n"
" vfmadb   3,18,0,3\n"
" vfmadb   4,20,0,4\n"
" vfmadb   5,22,0,5\n"
" vfmadb   6,24,0,6\n"
" vfmadb   7,26,0,7\n"
" vfmadb   8,28,0,8\n"
" vfmadb   9,30,0,9\n"
" vfmadb   2,17,1,2\n"
" vfmadb   3,19,1,3\n"
" vfmadb   4,21,1,4\n"
" vfmadb   5,23,1,5\n"
" vfmadb   6,25,1,6\n"
" vfmadb   7,27,1,7\n"
" vfmadb   8,29,1,8\n"
" vfmadb   9,31,1,9\n"
" vst 2,0(1,%[y])\n"
" vst 3,16(1,%[y])\n"
" vst 4,32(1,%[y])\n"
" vst 5,48(1,%[y])\n"
" vst 6,64(1,%[y])\n"
" vst 7,80(1,%[y])\n"
" vst 8,96(1,%[y])\n"
" vst 9,112(1,%[y])\n"
" agfi   1,128\n"
" brctg  0,labelzero%=\n" // Jump to label0
"labelone%= lghi    0,12\n" // label1
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelthree%=\n"	// Jump	to label3
" srlg  0,0,2\n"
"labeltwo%= vl  16,0(1,%[ap0])\n" // label2
" vl  17,0(1,%[ap1])\n"
" vl  18,16(1,%[ap0])\n"
" vl  19,16(1,%[ap1])\n"
" vl  2,0(1,%[y])\n"
" vl  3,16(1,%[y])\n"
" vfmadb   2,16,0,2\n"
" vfmadb   3,18,0,3\n"
" vfmadb   2,17,1,2\n"
" vfmadb   3,19,1,3\n"
" vst 2,0(1,%[y])\n"
" vst 3,16(1,%[y])\n"
" agfi   1,32\n"
" brctg  0,labeltwo%=\n" // Jump to label2
"labelthree%= nop 0\n" // label3
    : "+m"(*(FLOAT (*)[n]) y)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n]) ap0),[ap0] "a"(ap0),
       "m"(*(const FLOAT (*)[n]) ap1),[ap1] "a"(ap1),
       "m"(*(const FLOAT (*)[2]) x),[x] "a"(x),[alpha] "Q"(*alpha),
       [n] "r"(n)
    : "cc", "r0", "r1", "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7",
       "v8", "v9", "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23",
       "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31");
}

static void dgemv_kernel_4x1(BLASLONG n, FLOAT *a0, FLOAT *x, FLOAT *y,
                             FLOAT *alpha) {
__asm__ (
" vlrepg 0,0(%[x])\n"
" vlrepg 16,%[alpha]\n"
" vfmdb  0,0,16\n"
" xgr   1,1\n"
" lghi    0,-16\n"
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelone%=\n"	// Jump	to label1
" srlg  0,0,4\n"
"labelzero%= pfd 1,1024(1,%[a0])\n" // label0
" pfd 2,1024(1,%[y])\n"
" vl  16,0(1,%[a0])\n"
" vl  17,16(1,%[a0])\n"
" vl  18,32(1,%[a0])\n"
" vl  19,48(1,%[a0])\n"
" vl  20,64(1,%[a0])\n"
" vl  21,80(1,%[a0])\n"
" vl  22,96(1,%[a0])\n"
" vl  23,112(1,%[a0])\n"
" vl  24,0(1,%[y])\n"
" vl  25,16(1,%[y])\n"
" vl  26,32(1,%[y])\n"
" vl  27,48(1,%[y])\n"
" vl  28,64(1,%[y])\n"
" vl  29,80(1,%[y])\n"
" vl  30,96(1,%[y])\n"
" vl  31,112(1,%[y])\n"
" vfmadb   24,16,0,24\n"
" vfmadb   25,17,0,25\n"
" vfmadb   26,18,0,26\n"
" vfmadb   27,19,0,27\n"
" vfmadb   28,20,0,28\n"
" vfmadb   29,21,0,29\n"
" vfmadb   30,22,0,30\n"
" vfmadb   31,23,0,31\n"
" vst 24,0(1,%[y])\n"
" vst 25,16(1,%[y])\n"
" vst 26,32(1,%[y])\n"
" vst 27,48(1,%[y])\n"
" vst 28,64(1,%[y])\n"
" vst 29,80(1,%[y])\n"
" vst 30,96(1,%[y])\n"
" vst 31,112(1,%[y])\n"
" agfi   1,128\n"
" brctg  0,labelzero%=\n" // Jump to label0
"labelone%= lghi    0,12\n" // label1
" ngr     0,%[n]\n"
" ltgr    0,0\n"
" jz      labelthree%=\n"	// Jump	to label3
" srlg  0,0,2\n"
"labeltwo%= vl  16,0(1,%[a0])\n" // label2
" vl  17,16(1,%[a0])\n"
" vl  18,0(1,%[y])\n"
" vl  19,16(1,%[y])\n"
" vfmadb   18,16,0,18\n"
" vfmadb   19,17,0,19\n"
" vst 18,0(1,%[y])\n"
" vst 19,16(1,%[y])\n"
" agfi   1,32\n"
" brctg  0,labeltwo%=\n" // Jump to label2
"labelthree%= nop 0\n" // label3
    : "+m"(*(FLOAT (*)[n]) y)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n]) a0),[a0] "a"(a0),
       "m"(*(const FLOAT (*)[1]) x),[x] "a"(x),[alpha] "Q"(*alpha),
       [n] "r"(n)
    : "cc", "r0", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

static void add_y(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_dest) {
  BLASLONG i;
  for (i = 0; i < n; i++) {
    *dest += src[i];
    dest += inc_dest;
  }
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a,
          BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y,
          FLOAT *buffer) {
  BLASLONG i;
  FLOAT *a_ptr;
  FLOAT *x_ptr;
  FLOAT *y_ptr;
  FLOAT *ap[4];
  BLASLONG n1;
  BLASLONG m1;
  BLASLONG m2;
  BLASLONG m3;
  BLASLONG n2;
  BLASLONG lda4 = lda << 2;
  FLOAT xbuffer[8], *ybuffer;

  if (m < 1)
    return (0);
  if (n < 1)
    return (0);

  ybuffer = buffer;

  n1 = n >> 2;
  n2 = n & 3;

  m3 = m & 3;
  m1 = m & -4;
  m2 = (m & (NBMAX - 1)) - m3;

  y_ptr = y;

  BLASLONG NB = NBMAX;

  while (NB == NBMAX) {

    m1 -= NB;
    if (m1 < 0) {
      if (m2 == 0)
        break;
      NB = m2;
    }

    a_ptr = a;
    x_ptr = x;

    ap[0] = a_ptr;
    ap[1] = a_ptr + lda;
    ap[2] = ap[1] + lda;
    ap[3] = ap[2] + lda;

    if (inc_y != 1)
      memset(ybuffer, 0, NB * 8);
    else
      ybuffer = y_ptr;

    if (inc_x == 1) {

      for (i = 0; i < n1; i++) {
        dgemv_kernel_4x4(NB, ap, x_ptr, ybuffer, &alpha);
        ap[0] += lda4;
        ap[1] += lda4;
        ap[2] += lda4;
        ap[3] += lda4;
        a_ptr += lda4;
        x_ptr += 4;
      }

      if (n2 & 2) {
        dgemv_kernel_4x2(NB, ap, x_ptr, ybuffer, &alpha);
        a_ptr += lda * 2;
        x_ptr += 2;
      }

      if (n2 & 1) {
        dgemv_kernel_4x1(NB, a_ptr, x_ptr, ybuffer, &alpha);
        /* a_ptr += lda;
           x_ptr += 1; */

      }

    } else {

      for (i = 0; i < n1; i++) {
        xbuffer[0] = x_ptr[0];
        x_ptr += inc_x;
        xbuffer[1] = x_ptr[0];
        x_ptr += inc_x;
        xbuffer[2] = x_ptr[0];
        x_ptr += inc_x;
        xbuffer[3] = x_ptr[0];
        x_ptr += inc_x;
        dgemv_kernel_4x4(NB, ap, xbuffer, ybuffer, &alpha);
        ap[0] += lda4;
        ap[1] += lda4;
        ap[2] += lda4;
        ap[3] += lda4;
        a_ptr += lda4;
      }

      for (i = 0; i < n2; i++) {
        xbuffer[0] = x_ptr[0];
        x_ptr += inc_x;
        dgemv_kernel_4x1(NB, a_ptr, xbuffer, ybuffer, &alpha);
        a_ptr += lda;

      }

    }

    a += NB;
    if (inc_y != 1) {
      add_y(NB, ybuffer, y_ptr, inc_y);
      y_ptr += NB * inc_y;
    } else
      y_ptr += NB;

  }

  if (m3 == 0)
    return (0);

  if (m3 == 3) {
    a_ptr = a;
    x_ptr = x;
    FLOAT temp0 = 0.0;
    FLOAT temp1 = 0.0;
    FLOAT temp2 = 0.0;
    if (lda == 3 && inc_x == 1) {

      for (i = 0; i < (n & -4); i += 4) {

        temp0 += a_ptr[0] * x_ptr[0] + a_ptr[3] * x_ptr[1];
        temp1 += a_ptr[1] * x_ptr[0] + a_ptr[4] * x_ptr[1];
        temp2 += a_ptr[2] * x_ptr[0] + a_ptr[5] * x_ptr[1];

        temp0 += a_ptr[6] * x_ptr[2] + a_ptr[9] * x_ptr[3];
        temp1 += a_ptr[7] * x_ptr[2] + a_ptr[10] * x_ptr[3];
        temp2 += a_ptr[8] * x_ptr[2] + a_ptr[11] * x_ptr[3];

        a_ptr += 12;
        x_ptr += 4;
      }

      for (; i < n; i++) {
        temp0 += a_ptr[0] * x_ptr[0];
        temp1 += a_ptr[1] * x_ptr[0];
        temp2 += a_ptr[2] * x_ptr[0];
        a_ptr += 3;
        x_ptr++;
      }

    } else {

      for (i = 0; i < n; i++) {
        temp0 += a_ptr[0] * x_ptr[0];
        temp1 += a_ptr[1] * x_ptr[0];
        temp2 += a_ptr[2] * x_ptr[0];
        a_ptr += lda;
        x_ptr += inc_x;

      }

    }
    y_ptr[0] += alpha * temp0;
    y_ptr += inc_y;
    y_ptr[0] += alpha * temp1;
    y_ptr += inc_y;
    y_ptr[0] += alpha * temp2;
    return (0);
  }

  if (m3 == 2) {
    a_ptr = a;
    x_ptr = x;
    FLOAT temp0 = 0.0;
    FLOAT temp1 = 0.0;
    if (lda == 2 && inc_x == 1) {

      for (i = 0; i < (n & -4); i += 4) {
        temp0 += a_ptr[0] * x_ptr[0] + a_ptr[2] * x_ptr[1];
        temp1 += a_ptr[1] * x_ptr[0] + a_ptr[3] * x_ptr[1];
        temp0 += a_ptr[4] * x_ptr[2] + a_ptr[6] * x_ptr[3];
        temp1 += a_ptr[5] * x_ptr[2] + a_ptr[7] * x_ptr[3];
        a_ptr += 8;
        x_ptr += 4;

      }

      for (; i < n; i++) {
        temp0 += a_ptr[0] * x_ptr[0];
        temp1 += a_ptr[1] * x_ptr[0];
        a_ptr += 2;
        x_ptr++;
      }

    } else {

      for (i = 0; i < n; i++) {
        temp0 += a_ptr[0] * x_ptr[0];
        temp1 += a_ptr[1] * x_ptr[0];
        a_ptr += lda;
        x_ptr += inc_x;

      }

    }
    y_ptr[0] += alpha * temp0;
    y_ptr += inc_y;
    y_ptr[0] += alpha * temp1;
    return (0);
  }

  if (m3 == 1) {
    a_ptr = a;
    x_ptr = x;
    FLOAT temp = 0.0;
    if (lda == 1 && inc_x == 1) {

      for (i = 0; i < (n & -4); i += 4) {
        temp +=
          a_ptr[i] * x_ptr[i] + a_ptr[i + 1] * x_ptr[i + 1] + a_ptr[i +
                                                                    2] *
          x_ptr[i + 2] + a_ptr[i + 3] * x_ptr[i + 3];

      }

      for (; i < n; i++) {
        temp += a_ptr[i] * x_ptr[i];
      }

    } else {

      for (i = 0; i < n; i++) {
        temp += a_ptr[0] * x_ptr[0];
        a_ptr += lda;
        x_ptr += inc_x;
      }

    }
    y_ptr[0] += alpha * temp;
    return (0);
  }

  return (0);
}
