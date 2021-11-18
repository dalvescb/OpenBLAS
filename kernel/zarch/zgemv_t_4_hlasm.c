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

#define NBMAX 1024

static void zgemv_kernel_4x4(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y,
                             FLOAT *alpha) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
  register FLOAT *ap2 = ap[2];
  register FLOAT *ap3 = ap[3];
__asm__ (
" vzero  16\n"
" vzero  17\n"
" vzero  18\n"
" vzero  19\n"
" vzero  20\n"
" vzero  21\n"
" vzero  22\n"
" vzero  23\n"
" xgr   1,1\n"
" srlg  %[n],%[n],1\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n"
" pfd 1,1024(1,%[ap1])\n"
" pfd 1,1024(1,%[ap2])\n"
" pfd 1,1024(1,%[ap3])\n"
" pfd 1,1024(1,%[x])\n"
" vl     0,0(1,%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   1,8(1,%[x]),0\n"
" wflcdb 1,1\n"
" vleg   1,0(1,%[x]),1\n"
#else
" vleg   1,0(1,%[x]),1\n"
" vflcdb 1,1\n"
" vleg   1,8(1,%[x]),0\n"
#endif
" vlrepg 24,0(1,%[ap0])\n"
" vlrepg 25,8(1,%[ap0])\n"
" vlrepg 26,0(1,%[ap1])\n"
" vlrepg 27,8(1,%[ap1])\n"
" vlrepg 28,0(1,%[ap2])\n"
" vlrepg 29,8(1,%[ap2])\n"
" vlrepg 30,0(1,%[ap3])\n"
" vlrepg 31,8(1,%[ap3])\n"
" vfmadb   16,24,0,16\n"
" vfmadb   20,25,1,20\n"
" vfmadb   17,26,0,17\n"
" vfmadb   21,27,1,21\n"
" vfmadb   18,28,0,18\n"
" vfmadb   22,29,1,22\n"
" vfmadb   19,30,0,19\n"
" vfmadb   23,31,1,23\n"
" vl     0,16(1,%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   1,24(1,%[x]),0\n"
" wflcdb 1,1\n"
" vleg   1,16(1,%[x]),1\n"
#else
" vleg   1,16(1,%[x]),1\n"
" vflcdb 1,1\n"
" vleg   1,24(1,%[x]),0\n"
#endif
" vlrepg 24,16(1,%[ap0])\n"
" vlrepg 25,24(1,%[ap0])\n"
" vlrepg 26,16(1,%[ap1])\n"
" vlrepg 27,24(1,%[ap1])\n"
" vlrepg 28,16(1,%[ap2])\n"
" vlrepg 29,24(1,%[ap2])\n"
" vlrepg 30,16(1,%[ap3])\n"
" vlrepg 31,24(1,%[ap3])\n"
" vfmadb   16,24,0,16\n"
" vfmadb   20,25,1,20\n"
" vfmadb   17,26,0,17\n"
" vfmadb   21,27,1,21\n"
" vfmadb   18,28,0,18\n"
" vfmadb   22,29,1,22\n"
" vfmadb   19,30,0,19\n"
" vfmadb   23,31,1,23\n"
" agfi   1,32\n"
" brctg  %[n],labelzero%=\n"
" vfadb  16,16,20\n"
" vfadb  17,17,21\n"
" vfadb  18,18,22\n"
" vfadb  19,19,23\n"
" vpdi 20,16,16,4\n"
" vpdi 21,17,17,4\n"
" vpdi 22,18,18,4\n"
" vpdi 23,19,19,4\n"
#if !defined(XCONJ)
" vlrepg 24,0(%[alpha])\n"
" vleg   25,8(%[alpha]),0\n"
" wflcdb 25,25\n"
" vleg   25,8(%[alpha]),1\n"
#else
" vleg   24,0(%[alpha]),1\n"
" vflcdb 24,24\n"
" vleg   24,0(%[alpha]),0\n"
" vlrepg 25,8(%[alpha])\n"
#endif
" vl  26,0(%[y])\n"
" vl  27,16(%[y])\n"
" vl  28,32(%[y])\n"
" vl  29,48(%[y])\n"
" vfmadb   26,16,24,26\n"
" vfmadb   26,20,25,26\n"
" vfmadb   27,17,24,27\n"
" vfmadb   27,21,25,27\n"
" vfmadb   28,18,24,28\n"
" vfmadb   28,22,25,28\n"
" vfmadb   29,19,24,29\n"
" vfmadb   29,23,25,29\n"
" vst  26,0(%[y])\n"
" vst  27,16(%[y])\n"
" vst  28,32(%[y])\n"
" vst  29,48(%[y])\n"
    : "+m"(*(FLOAT (*)[8]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n * 2]) ap0),[ap0] "a"(ap0),
       "m"(*(const FLOAT (*)[n * 2]) ap1),[ap1] "a"(ap1),
       "m"(*(const FLOAT (*)[n * 2]) ap2),[ap2] "a"(ap2),
       "m"(*(const FLOAT (*)[n * 2]) ap3),[ap3] "a"(ap3),
       "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x),
       "m"(*(const FLOAT (*)[2]) alpha),[alpha] "a"(alpha)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

static void zgemv_kernel_4x2(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y,
                             FLOAT *alpha) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
__asm__ (
" vzero  16\n"
" vzero  17\n"
" vzero  18\n"
" vzero  19\n"
" xgr   1,1\n"
" srlg  %[n],%[n],1\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n"
" pfd 1,1024(1,%[ap1])\n"
" pfd 1,1024(1,%[x])\n"
" vl     0,0(1,%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   1,8(1,%[x]),0\n"
" wflcdb 1,1\n"
" vleg   1,0(1,%[x]),1\n"
#else
" vleg   1,0(1,%[x]),1\n"
" vflcdb 1,1\n"
" vleg   1,8(1,%[x]),0\n"
#endif
" vlrepg 20,0(1,%[ap0])\n"
" vlrepg 21,8(1,%[ap0])\n"
" vlrepg 22,0(1,%[ap1])\n"
" vlrepg 23,8(1,%[ap1])\n"
" vfmadb   16,20,0,16\n"
" vfmadb   18,21,1,18\n"
" vfmadb   17,22,0,17\n"
" vfmadb   19,23,1,19\n"
" vl     0,16(1,%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   1,24(1,%[x]),0\n"
" wflcdb 1,1\n"
" vleg   1,16(1,%[x]),1\n"
#else
" vleg   1,16(1,%[x]),1\n"
" vflcdb 1,1\n"
" vleg   1,24(1,%[x]),0\n"
#endif
" vlrepg 20,16(1,%[ap0])\n"
" vlrepg 21,24(1,%[ap0])\n"
" vlrepg 22,16(1,%[ap1])\n"
" vlrepg 23,24(1,%[ap1])\n"
" vfmadb   16,20,0,16\n"
" vfmadb   18,21,1,18\n"
" vfmadb   17,22,0,17\n"
" vfmadb   19,23,1,19\n"
" agfi   1,32\n"
" brctg  %[n],labelzero%=\n"
" vfadb  16,16,18\n"
" vfadb  17,17,19\n"
" vpdi 18,16,16,4\n"
" vpdi 19,17,17,4\n"
#if !defined(XCONJ)
" vlrepg 20,0(%[alpha])\n"
" vleg   21,8(%[alpha]),0\n"
" wflcdb 21,21\n"
" vleg   21,8(%[alpha]),1\n"
#else
" vleg   20,0(%[alpha]),1\n"
" vflcdb 20,20\n"
" vleg   20,0(%[alpha]),0\n"
" vlrepg 21,8(%[alpha])\n"
#endif
" vl  22,0(%[y])\n"
" vl  23,16(%[y])\n"
" vfmadb   22,16,20,22\n"
" vfmadb   22,18,21,22\n"
" vfmadb   23,17,20,23\n"
" vfmadb   23,19,21,23\n"
" vst  22,0(%[y])\n"
" vst  23,16(%[y])\n"
    : "+m"(*(FLOAT (*)[4]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n * 2]) ap0),[ap0] "a"(ap0),
       "m"(*(const FLOAT (*)[n * 2]) ap1),[ap1] "a"(ap1),
       "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x),
       "m"(*(const FLOAT (*)[2]) alpha),[alpha] "a"(alpha)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23");
}

static void zgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y,
                             FLOAT *alpha) {
__asm__ (
" vzero  16\n"
" vzero  17\n"
" xgr   1,1\n"
" srlg  %[n],%[n],1\n"
"labelzero%= pfd 1,1024(1,%[ap])\n"
" pfd 1,1024(1,%[x])\n"
" vl     0,0(1,%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   1,8(1,%[x]),0\n"
" wflcdb 1,1\n"
" vleg   1,0(1,%[x]),1\n"
#else
" vleg   1,0(1,%[x]),1\n"
" vflcdb 1,1\n"
" vleg   1,8(1,%[x]),0\n"
#endif
" vlrepg 18,0(1,%[ap])\n"
" vlrepg 19,8(1,%[ap])\n"
" vfmadb   16,18,0,16\n"
" vfmadb   17,19,1,17\n"
" vl     0,16(1,%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   1,24(1,%[x]),0\n"
" wflcdb 1,1\n"
" vleg   1,16(1,%[x]),1\n"
#else
" vleg   1,16(1,%[x]),1\n"
" vflcdb 1,1\n"
" vleg   1,24(1,%[x]),0\n"
#endif
" vlrepg 18,16(1,%[ap])\n"
" vlrepg 19,24(1,%[ap])\n"
" vfmadb   16,18,0,16\n"
" vfmadb   17,19,1,17\n"
" agfi   1,32\n"
" brctg  %[n],labelzero%=\n"
" vfadb  16,16,17\n"
" vpdi 17,16,16,4\n"
#if !defined(XCONJ)
" vlrepg 18,0(%[alpha])\n"
" vleg   19,8(%[alpha]),0\n"
" wflcdb 19,19\n"
" vleg   19,8(%[alpha]),1\n"
#else
" vleg   18,0(%[alpha]),1\n"
" vflcdb 18,18\n"
" vleg   18,0(%[alpha]),0\n"
" vlrepg 19,8(%[alpha])\n"
#endif
" vl  0,0(%[y])\n"
" vfmadb   0,16,18,0\n"
" vfmadb   0,17,19,0\n"
" vst  0,0(%[y])\n"
    : "+m"(*(FLOAT (*)[2]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n * 2]) ap),[ap] "a"(ap),
       "m"(*(const FLOAT (*)[n * 2]) x),[x] "a"(x),
       "m"(*(const FLOAT (*)[2]) alpha),[alpha] "a"(alpha)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19");
}

static void copy_x(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src) {
  BLASLONG i;
  for (i = 0; i < n; i++) {
    *dest = *src;
    *(dest + 1) = *(src + 1);
    dest += 2;
    src += inc_src;
  }
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i,
          FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y,
          BLASLONG inc_y, FLOAT *buffer) {
  BLASLONG i;
  BLASLONG j;
  FLOAT *a_ptr;
  FLOAT *x_ptr;
  FLOAT *y_ptr;
  FLOAT *ap[8];
  BLASLONG n1;
  BLASLONG m1;
  BLASLONG m2;
  BLASLONG m3;
  BLASLONG n2;
  BLASLONG lda4;
  FLOAT ybuffer[8], *xbuffer;
  FLOAT alpha[2];

  if (m < 1)
    return (0);
  if (n < 1)
    return (0);

  inc_x <<= 1;
  inc_y <<= 1;
  lda <<= 1;
  lda4 = lda << 2;

  xbuffer = buffer;

  n1 = n >> 2;
  n2 = n & 3;

  m3 = m & 3;
  m1 = m - m3;
  m2 = (m & (NBMAX - 1)) - m3;

  alpha[0] = alpha_r;
  alpha[1] = alpha_i;

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
    ap[0] = a_ptr;
    ap[1] = a_ptr + lda;
    ap[2] = ap[1] + lda;
    ap[3] = ap[2] + lda;
    if (inc_x != 2)
      copy_x(NB, x_ptr, xbuffer, inc_x);
    else
      xbuffer = x_ptr;

    if (inc_y == 2) {

      for (i = 0; i < n1; i++) {
        zgemv_kernel_4x4(NB, ap, xbuffer, y_ptr, alpha);
        ap[0] += lda4;
        ap[1] += lda4;
        ap[2] += lda4;
        ap[3] += lda4;
        a_ptr += lda4;
        y_ptr += 8;

      }

      if (n2 & 2) {
        zgemv_kernel_4x2(NB, ap, xbuffer, y_ptr, alpha);
        a_ptr += lda * 2;
        y_ptr += 4;

      }

      if (n2 & 1) {
        zgemv_kernel_4x1(NB, a_ptr, xbuffer, y_ptr, alpha);
        /* a_ptr += lda;
           y_ptr += 2; */

      }

    } else {

      for (i = 0; i < n1; i++) {
        memset(ybuffer, 0, sizeof(ybuffer));
        zgemv_kernel_4x4(NB, ap, xbuffer, ybuffer, alpha);
        ap[0] += lda4;
        ap[1] += lda4;
        ap[2] += lda4;
        ap[3] += lda4;
        a_ptr += lda4;

        y_ptr[0] += ybuffer[0];
        y_ptr[1] += ybuffer[1];
        y_ptr += inc_y;
        y_ptr[0] += ybuffer[2];
        y_ptr[1] += ybuffer[3];
        y_ptr += inc_y;
        y_ptr[0] += ybuffer[4];
        y_ptr[1] += ybuffer[5];
        y_ptr += inc_y;
        y_ptr[0] += ybuffer[6];
        y_ptr[1] += ybuffer[7];
        y_ptr += inc_y;

      }

      for (i = 0; i < n2; i++) {
        memset(ybuffer, 0, sizeof(ybuffer));
        zgemv_kernel_4x1(NB, a_ptr, xbuffer, ybuffer, alpha);
        a_ptr += lda;
        y_ptr[0] += ybuffer[0];
        y_ptr[1] += ybuffer[1];
        y_ptr += inc_y;

      }

    }
    a += 2 * NB;
    x += NB * inc_x;
  }

  if (m3 == 0)
    return (0);

  x_ptr = x;
  j = 0;
  a_ptr = a;
  y_ptr = y;

  if (m3 == 3) {

    FLOAT temp_r;
    FLOAT temp_i;
    FLOAT x0 = x_ptr[0];
    FLOAT x1 = x_ptr[1];
    x_ptr += inc_x;
    FLOAT x2 = x_ptr[0];
    FLOAT x3 = x_ptr[1];
    x_ptr += inc_x;
    FLOAT x4 = x_ptr[0];
    FLOAT x5 = x_ptr[1];
    while (j < n) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
      temp_r = a_ptr[0] * x0 - a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 + a_ptr[1] * x0;
      temp_r += a_ptr[2] * x2 - a_ptr[3] * x3;
      temp_i += a_ptr[2] * x3 + a_ptr[3] * x2;
      temp_r += a_ptr[4] * x4 - a_ptr[5] * x5;
      temp_i += a_ptr[4] * x5 + a_ptr[5] * x4;
#else

      temp_r = a_ptr[0] * x0 + a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 - a_ptr[1] * x0;
      temp_r += a_ptr[2] * x2 + a_ptr[3] * x3;
      temp_i += a_ptr[2] * x3 - a_ptr[3] * x2;
      temp_r += a_ptr[4] * x4 + a_ptr[5] * x5;
      temp_i += a_ptr[4] * x5 - a_ptr[5] * x4;
#endif

#if !defined(XCONJ)
      y_ptr[0] += alpha_r * temp_r - alpha_i * temp_i;
      y_ptr[1] += alpha_r * temp_i + alpha_i * temp_r;
#else
      y_ptr[0] += alpha_r * temp_r + alpha_i * temp_i;
      y_ptr[1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif

      a_ptr += lda;
      y_ptr += inc_y;
      j++;
    }
    return (0);
  }

  if (m3 == 2) {

    FLOAT temp_r;
    FLOAT temp_i;
    FLOAT temp_r1;
    FLOAT temp_i1;
    FLOAT x0 = x_ptr[0];
    FLOAT x1 = x_ptr[1];
    x_ptr += inc_x;
    FLOAT x2 = x_ptr[0];
    FLOAT x3 = x_ptr[1];
    FLOAT ar = alpha[0];
    FLOAT ai = alpha[1];

    while (j < (n & -2)) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
      temp_r = a_ptr[0] * x0 - a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 + a_ptr[1] * x0;
      temp_r += a_ptr[2] * x2 - a_ptr[3] * x3;
      temp_i += a_ptr[2] * x3 + a_ptr[3] * x2;
      a_ptr += lda;
      temp_r1 = a_ptr[0] * x0 - a_ptr[1] * x1;
      temp_i1 = a_ptr[0] * x1 + a_ptr[1] * x0;
      temp_r1 += a_ptr[2] * x2 - a_ptr[3] * x3;
      temp_i1 += a_ptr[2] * x3 + a_ptr[3] * x2;
#else

      temp_r = a_ptr[0] * x0 + a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 - a_ptr[1] * x0;
      temp_r += a_ptr[2] * x2 + a_ptr[3] * x3;
      temp_i += a_ptr[2] * x3 - a_ptr[3] * x2;
      a_ptr += lda;
      temp_r1 = a_ptr[0] * x0 + a_ptr[1] * x1;
      temp_i1 = a_ptr[0] * x1 - a_ptr[1] * x0;
      temp_r1 += a_ptr[2] * x2 + a_ptr[3] * x3;
      temp_i1 += a_ptr[2] * x3 - a_ptr[3] * x2;
#endif

#if !defined(XCONJ)
      y_ptr[0] += ar * temp_r - ai * temp_i;
      y_ptr[1] += ar * temp_i + ai * temp_r;
      y_ptr += inc_y;
      y_ptr[0] += ar * temp_r1 - ai * temp_i1;
      y_ptr[1] += ar * temp_i1 + ai * temp_r1;
#else
      y_ptr[0] += ar * temp_r + ai * temp_i;
      y_ptr[1] -= ar * temp_i - ai * temp_r;
      y_ptr += inc_y;
      y_ptr[0] += ar * temp_r1 + ai * temp_i1;
      y_ptr[1] -= ar * temp_i1 - ai * temp_r1;
#endif

      a_ptr += lda;
      y_ptr += inc_y;
      j += 2;
    }

    while (j < n) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
      temp_r = a_ptr[0] * x0 - a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 + a_ptr[1] * x0;
      temp_r += a_ptr[2] * x2 - a_ptr[3] * x3;
      temp_i += a_ptr[2] * x3 + a_ptr[3] * x2;
#else

      temp_r = a_ptr[0] * x0 + a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 - a_ptr[1] * x0;
      temp_r += a_ptr[2] * x2 + a_ptr[3] * x3;
      temp_i += a_ptr[2] * x3 - a_ptr[3] * x2;
#endif

#if !defined(XCONJ)
      y_ptr[0] += ar * temp_r - ai * temp_i;
      y_ptr[1] += ar * temp_i + ai * temp_r;
#else
      y_ptr[0] += ar * temp_r + ai * temp_i;
      y_ptr[1] -= ar * temp_i - ai * temp_r;
#endif

      a_ptr += lda;
      y_ptr += inc_y;
      j++;
    }

    return (0);
  }

  if (m3 == 1) {

    FLOAT temp_r;
    FLOAT temp_i;
    FLOAT temp_r1;
    FLOAT temp_i1;
    FLOAT x0 = x_ptr[0];
    FLOAT x1 = x_ptr[1];
    FLOAT ar = alpha[0];
    FLOAT ai = alpha[1];

    while (j < (n & -2)) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
      temp_r = a_ptr[0] * x0 - a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 + a_ptr[1] * x0;
      a_ptr += lda;
      temp_r1 = a_ptr[0] * x0 - a_ptr[1] * x1;
      temp_i1 = a_ptr[0] * x1 + a_ptr[1] * x0;
#else

      temp_r = a_ptr[0] * x0 + a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 - a_ptr[1] * x0;
      a_ptr += lda;
      temp_r1 = a_ptr[0] * x0 + a_ptr[1] * x1;
      temp_i1 = a_ptr[0] * x1 - a_ptr[1] * x0;
#endif

#if !defined(XCONJ)
      y_ptr[0] += ar * temp_r - ai * temp_i;
      y_ptr[1] += ar * temp_i + ai * temp_r;
      y_ptr += inc_y;
      y_ptr[0] += ar * temp_r1 - ai * temp_i1;
      y_ptr[1] += ar * temp_i1 + ai * temp_r1;
#else
      y_ptr[0] += ar * temp_r + ai * temp_i;
      y_ptr[1] -= ar * temp_i - ai * temp_r;
      y_ptr += inc_y;
      y_ptr[0] += ar * temp_r1 + ai * temp_i1;
      y_ptr[1] -= ar * temp_i1 - ai * temp_r1;
#endif

      a_ptr += lda;
      y_ptr += inc_y;
      j += 2;
    }

    while (j < n) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
      temp_r = a_ptr[0] * x0 - a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 + a_ptr[1] * x0;
#else

      temp_r = a_ptr[0] * x0 + a_ptr[1] * x1;
      temp_i = a_ptr[0] * x1 - a_ptr[1] * x0;
#endif

#if !defined(XCONJ)
      y_ptr[0] += ar * temp_r - ai * temp_i;
      y_ptr[1] += ar * temp_i + ai * temp_r;
#else
      y_ptr[0] += ar * temp_r + ai * temp_i;
      y_ptr[1] -= ar * temp_i - ai * temp_r;
#endif

      a_ptr += lda;
      y_ptr += inc_y;
      j++;
    }
    return (0);
  }

  return (0);
}
