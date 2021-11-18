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

static void zgemv_kernel_4x4(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
  register FLOAT *ap2 = ap[2];
  register FLOAT *ap3 = ap[3];
__asm__ (
" vl     16,0(%[x])\n"
" vl     17,16(%[x])\n"
" vl     18,32(%[x])\n"
" vl     19,48(%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   20,8(%[x]),0\n"
" wflcdb 20,20\n"
" vleg   20,0(%[x]),1\n"
" vleg   21,24(%[x]),0\n"
" wflcdb 21,21\n"
" vleg   21,16(%[x]),1\n"
" vleg   22,40(%[x]),0\n"
" wflcdb 22,22\n"
" vleg   22,32(%[x]),1\n"
" vleg   23,56(%[x]),0\n"
" wflcdb 23,23\n"
" vleg   23,48(%[x]),1\n"
#else
" vleg   20,0(%[x]),1\n"
" vflcdb 20,20\n"
" vleg   20,8(%[x]),0\n"
" vleg   21,16(%[x]),1\n"
" vflcdb 21,21\n"
" vleg   21,24(%[x]),0\n"
" vleg   22,32(%[x]),1\n"
" vflcdb 22,22\n"
" vleg   22,40(%[x]),0\n"
" vleg   23,48(%[x]),1\n"
" vflcdb 23,23\n"
" vleg   23,56(%[x]),0\n"
#endif
" xgr   1,1\n"
" srlg  %[n],%[n],1\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n"
" pfd 1,1024(1,%[ap1])\n"
" pfd 1,1024(1,%[ap2])\n"
" pfd 1,1024(1,%[ap3])\n"
" pfd 2,1024(1,%[y])\n"
" vl  0,0(1,%[y])\n"
" vl  1,16(1,%[y])\n"
" vlrepg 24,0(1,%[ap0])\n"
" vlrepg 25,8(1,%[ap0])\n"
" vlrepg 26,0(1,%[ap1])\n"
" vlrepg 27,8(1,%[ap1])\n"
" vlrepg 28,16(1,%[ap0])\n"
" vlrepg 29,24(1,%[ap0])\n"
" vlrepg 30,16(1,%[ap1])\n"
" vlrepg 31,24(1,%[ap1])\n"
" vfmadb   0,24,16,0\n"
" vfmadb   1,28,16,1\n"
" vfmadb   0,25,20,0\n"
" vfmadb   1,29,20,1\n"
" vfmadb   0,26,17,0\n"
" vfmadb   1,30,17,1\n"
" vfmadb   0,27,21,0\n"
" vfmadb   1,31,21,1\n"
" vlrepg 24,0(1,%[ap2])\n"
" vlrepg 25,8(1,%[ap2])\n"
" vlrepg 26,0(1,%[ap3])\n"
" vlrepg 27,8(1,%[ap3])\n"
" vlrepg 28,16(1,%[ap2])\n"
" vlrepg 29,24(1,%[ap2])\n"
" vlrepg 30,16(1,%[ap3])\n"
" vlrepg 31,24(1,%[ap3])\n"
" vfmadb   0,24,18,0\n"
" vfmadb   1,28,18,1\n"
" vfmadb   0,25,22,0\n"
" vfmadb   1,29,22,1\n"
" vfmadb   0,26,19,0\n"
" vfmadb   1,30,19,1\n"
" vfmadb   0,27,23,0\n"
" vfmadb   1,31,23,1\n"
" vst 0,0(1,%[y])\n"
" vst 1,16(1,%[y])\n"
" agfi   1,32\n"
" brctg  %[n],labelzero%=\n"
    : "+m"(*(FLOAT (*)[n * 2]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n * 2]) ap0),[ap0] "a"(ap0),
       "m"(*(const FLOAT (*)[n * 2]) ap1),[ap1] "a"(ap1),
       "m"(*(const FLOAT (*)[n * 2]) ap2),[ap2] "a"(ap2),
       "m"(*(const FLOAT (*)[n * 2]) ap3),[ap3] "a"(ap3),
       "m"(*(const FLOAT (*)[8]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

static void zgemv_kernel_4x2(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y) {
  register FLOAT *ap0 = ap[0];
  register FLOAT *ap1 = ap[1];
__asm__ (
" vl     16,0(%[x])\n"
" vl     17,16(%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   18,8(%[x]),0\n"
" wflcdb 18,18\n"
" vleg   18,0(%[x]),1\n"
" vleg   19,24(%[x]),0\n"
" wflcdb 19,19\n"
" vleg   19,16(%[x]),1\n"
#else
" vleg   18,0(%[x]),1\n"
" vflcdb 18,18\n"
" vleg   18,8(%[x]),0\n"
" vleg   19,16(%[x]),1\n"
" vflcdb 19,19\n"
" vleg   19,24(%[x]),0\n"
#endif
" xgr   1,1\n"
" srlg  %[n],%[n],1\n"
"labelzero%= pfd 1,1024(1,%[ap0])\n"
" pfd 1,1024(1,%[ap1])\n"
" pfd 2,1024(1,%[y])\n"
" vl  0,0(1,%[y])\n"
" vl  1,16(1,%[y])\n"
" vlrepg 20,0(1,%[ap0])\n"
" vlrepg 21,8(1,%[ap0])\n"
" vlrepg 22,0(1,%[ap1])\n"
" vlrepg 23,8(1,%[ap1])\n"
" vlrepg 24,16(1,%[ap0])\n"
" vlrepg 25,24(1,%[ap0])\n"
" vlrepg 26,16(1,%[ap1])\n"
" vlrepg 27,24(1,%[ap1])\n"
" vfmadb   0,20,16,0\n"
" vfmadb   1,24,16,1\n"
" vfmadb   0,21,18,0\n"
" vfmadb   1,25,18,1\n"
" vfmadb   0,22,17,0\n"
" vfmadb   1,26,17,1\n"
" vfmadb   0,23,19,0\n"
" vfmadb   1,27,19,1\n"
" vst 0,0(1,%[y])\n"
" vst 1,16(1,%[y])\n"
" agfi   1,32\n"
" brctg  %[n],labelzero%=\n"
    : "+m"(*(FLOAT (*)[n * 2]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n * 2]) ap0),[ap0] "a"(ap0),
         "m"(*(const FLOAT (*)[n * 2]) ap1),[ap1] "a"(ap1),
         "m"(*(const FLOAT (*)[4]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27");
}

static void zgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y) {
__asm__ (
" vl     16,0(%[x])\n"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
" vleg   17,8(%[x]),0\n"
" wflcdb 17,17\n"
" vleg   17,0(%[x]),1\n"
#else
" vleg   17,0(%[x]),1\n"
" vflcdb 17,17\n"
" vleg   17,8(%[x]),0\n"
#endif
" xgr   1,1\n"
" srlg  %[n],%[n],1\n"
"labelzero%= pfd 1,1024(1,%[ap])\n"
" pfd 2,1024(1,%[y])\n"
" vl  0,0(1,%[y])\n"
" vl  1,16(1,%[y])\n"
" vlrepg 18,0(1,%[ap])\n"
" vlrepg 19,8(1,%[ap])\n"
" vlrepg 20,16(1,%[ap])\n"
" vlrepg 21,24(1,%[ap])\n"
" vfmadb   0,18,16,0\n"
" vfmadb   1,20,16,1\n"
" vfmadb   0,19,17,0\n"
" vfmadb   1,21,17,1\n"
" vst 0,0(1,%[y])\n"
" vst 1,16(1,%[y])\n"
" agfi   1,32\n"
" brctg  %[n],labelzero%=\n"
    : "+m"(*(FLOAT (*)[n * 2]) y),[n] "+&r"(n)
    : [y] "a"(y), "m"(*(const FLOAT (*)[n * 2]) ap),[ap] "a"(ap),
       "m"(*(const FLOAT (*)[2]) x),[x] "a"(x)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21");
}

static void add_y_4(BLASLONG n, FLOAT *src, FLOAT *dest, FLOAT alpha_r,
                    FLOAT alpha_i) {
  __asm__(
#if !defined(XCONJ)
" vlrepg 0,%[alpha_r]\n"
" vleg   1,%[alpha_i],0\n"
" wflcdb 1,1\n"
" vleg   1,%[alpha_i],1\n"
#else
" vleg   0,%[alpha_r],1\n"
" vflcdb 0,0\n"
" vleg   0,%[alpha_r],0\n"
" vlrepg 1,%[alpha_i]\n"
#endif
" xgr   1,1\n"
" srlg  %[n],%[n],2\n"
"labelzero%= pfd 1,1024(1,%[src])\n"
" pfd 2,1024(1,%[dest])\n"
" vl   16,0(1,%[src])\n"
" vl   17,16(1,%[src])\n"
" vl   18,32(1,%[src])\n"
" vl   19,48(1,%[src])\n"
" vl   20,0(1,%[dest])\n"
" vl   21,16(1,%[dest])\n"
" vl   22,32(1,%[dest])\n"
" vl   23,48(1,%[dest])\n"
" vpdi 24,16,16,4\n"
" vpdi 25,17,17,4\n"
" vpdi 26,18,18,4\n"
" vpdi 27,19,19,4\n"
" vfmadb 28,16,0,20\n"
" vfmadb 29,17,0,21\n"
" vfmadb 30,18,0,22\n"
" vfmadb 31,19,0,23\n"
" vfmadb 28,24,1,28\n"
" vfmadb 29,25,1,29\n"
" vfmadb 30,26,1,30\n"
" vfmadb 31,27,1,31\n"
" vst 28,0(1,%[dest])\n"
" vst 29,16(1,%[dest])\n"
" vst 30,32(1,%[dest])\n"
" vst 31,48(1,%[dest])\n"
" agfi   1,64\n"
" brctg  %[n],labelzero%=\n"
    : "+m"(*(FLOAT (*)[n * 2]) dest),[n] "+&r"(n)
    : [dest] "a"(dest), "m"(*(const FLOAT (*)[n * 2]) src),
       [src] "a"(src),[alpha_r] "Q"(alpha_r),[alpha_i] "Q"(alpha_i)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

static void add_y(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_dest,
                  FLOAT alpha_r, FLOAT alpha_i) {
  BLASLONG i;

  if (inc_dest != 2) {

    FLOAT temp_r;
    FLOAT temp_i;
    for (i = 0; i < n; i++) {
#if !defined(XCONJ)
      temp_r = alpha_r * src[0] - alpha_i * src[1];
      temp_i = alpha_r * src[1] + alpha_i * src[0];
#else
      temp_r = alpha_r * src[0] + alpha_i * src[1];
      temp_i = -alpha_r * src[1] + alpha_i * src[0];
#endif

      *dest += temp_r;
      *(dest + 1) += temp_i;

      src += 2;
      dest += inc_dest;
    }
    return;
  }

  add_y_4(n, src, dest, alpha_r, alpha_i);
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i,
          FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y,
          BLASLONG inc_y, FLOAT *buffer) {
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
  BLASLONG lda4;
  FLOAT xbuffer[8], *ybuffer;

  if (m < 1)
    return (0);
  if (n < 1)
    return (0);

  ybuffer = buffer;

  inc_x *= 2;
  inc_y *= 2;
  lda *= 2;
  lda4 = 4 * lda;

  n1 = n / 4;
  n2 = n % 4;

  m3 = m % 4;
  m1 = m - (m % 4);
  m2 = (m % NBMAX) - (m % 4);

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
    ap[0] = a_ptr;
    ap[1] = a_ptr + lda;
    ap[2] = ap[1] + lda;
    ap[3] = ap[2] + lda;
    x_ptr = x;
    //zero_y(NB,ybuffer);
    memset(ybuffer, 0, NB * 16);

    if (inc_x == 2) {

      for (i = 0; i < n1; i++) {
        zgemv_kernel_4x4(NB, ap, x_ptr, ybuffer);
        ap[0] += lda4;
        ap[1] += lda4;
        ap[2] += lda4;
        ap[3] += lda4;
        a_ptr += lda4;
        x_ptr += 8;
      }

      if (n2 & 2) {
        zgemv_kernel_4x2(NB, ap, x_ptr, ybuffer);
        x_ptr += 4;
        a_ptr += 2 * lda;

      }

      if (n2 & 1) {
        zgemv_kernel_4x1(NB, a_ptr, x_ptr, ybuffer);
        /* x_ptr += 2;  
           a_ptr += lda; */

      }
    } else {

      for (i = 0; i < n1; i++) {

        xbuffer[0] = x_ptr[0];
        xbuffer[1] = x_ptr[1];
        x_ptr += inc_x;
        xbuffer[2] = x_ptr[0];
        xbuffer[3] = x_ptr[1];
        x_ptr += inc_x;
        xbuffer[4] = x_ptr[0];
        xbuffer[5] = x_ptr[1];
        x_ptr += inc_x;
        xbuffer[6] = x_ptr[0];
        xbuffer[7] = x_ptr[1];
        x_ptr += inc_x;

        zgemv_kernel_4x4(NB, ap, xbuffer, ybuffer);
        ap[0] += lda4;
        ap[1] += lda4;
        ap[2] += lda4;
        ap[3] += lda4;
        a_ptr += lda4;
      }

      for (i = 0; i < n2; i++) {
        xbuffer[0] = x_ptr[0];
        xbuffer[1] = x_ptr[1];
        x_ptr += inc_x;
        zgemv_kernel_4x1(NB, a_ptr, xbuffer, ybuffer);
        a_ptr += 1 * lda;

      }

    }

    add_y(NB, ybuffer, y_ptr, inc_y, alpha_r, alpha_i);
    a += 2 * NB;
    y_ptr += NB * inc_y;
  }

  if (m3 == 0)
    return (0);

  if (m3 == 1) {
    a_ptr = a;
    x_ptr = x;
    FLOAT temp_r = 0.0;
    FLOAT temp_i = 0.0;

    if (lda == 2 && inc_x == 2) {

      for (i = 0; i < (n & -2); i += 2) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        temp_r += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
        temp_r += a_ptr[2] * x_ptr[2] - a_ptr[3] * x_ptr[3];
        temp_i += a_ptr[2] * x_ptr[3] + a_ptr[3] * x_ptr[2];
#else
        temp_r += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
        temp_r += a_ptr[2] * x_ptr[2] + a_ptr[3] * x_ptr[3];
        temp_i += a_ptr[2] * x_ptr[3] - a_ptr[3] * x_ptr[2];
#endif

        a_ptr += 4;
        x_ptr += 4;
      }

      for (; i < n; i++) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        temp_r += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
#else
        temp_r += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
#endif

        a_ptr += 2;
        x_ptr += 2;
      }

    } else {

      for (i = 0; i < n; i++) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        temp_r += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
#else
        temp_r += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
#endif

        a_ptr += lda;
        x_ptr += inc_x;
      }

    }
#if !defined(XCONJ)
    y_ptr[0] += alpha_r * temp_r - alpha_i * temp_i;
    y_ptr[1] += alpha_r * temp_i + alpha_i * temp_r;
#else
    y_ptr[0] += alpha_r * temp_r + alpha_i * temp_i;
    y_ptr[1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif
    return (0);
  }

  if (m3 == 2) {
    a_ptr = a;
    x_ptr = x;
    FLOAT temp_r0 = 0.0;
    FLOAT temp_i0 = 0.0;
    FLOAT temp_r1 = 0.0;
    FLOAT temp_i1 = 0.0;

    if (lda == 4 && inc_x == 2) {

      for (i = 0; i < (n & -2); i += 2) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )

        temp_r0 += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] - a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] + a_ptr[3] * x_ptr[0];

        temp_r0 += a_ptr[4] * x_ptr[2] - a_ptr[5] * x_ptr[3];
        temp_i0 += a_ptr[4] * x_ptr[3] + a_ptr[5] * x_ptr[2];
        temp_r1 += a_ptr[6] * x_ptr[2] - a_ptr[7] * x_ptr[3];
        temp_i1 += a_ptr[6] * x_ptr[3] + a_ptr[7] * x_ptr[2];

#else
        temp_r0 += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] + a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] - a_ptr[3] * x_ptr[0];

        temp_r0 += a_ptr[4] * x_ptr[2] + a_ptr[5] * x_ptr[3];
        temp_i0 += a_ptr[4] * x_ptr[3] - a_ptr[5] * x_ptr[2];
        temp_r1 += a_ptr[6] * x_ptr[2] + a_ptr[7] * x_ptr[3];
        temp_i1 += a_ptr[6] * x_ptr[3] - a_ptr[7] * x_ptr[2];

#endif

        a_ptr += 8;
        x_ptr += 4;
      }

      for (; i < n; i++) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        temp_r0 += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] - a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] + a_ptr[3] * x_ptr[0];
#else
        temp_r0 += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] + a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] - a_ptr[3] * x_ptr[0];
#endif

        a_ptr += 4;
        x_ptr += 2;
      }

    } else {

      for (i = 0; i < n; i++) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        temp_r0 += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] - a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] + a_ptr[3] * x_ptr[0];
#else
        temp_r0 += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] + a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] - a_ptr[3] * x_ptr[0];
#endif

        a_ptr += lda;
        x_ptr += inc_x;
      }

    }
#if !defined(XCONJ)
    y_ptr[0] += alpha_r * temp_r0 - alpha_i * temp_i0;
    y_ptr[1] += alpha_r * temp_i0 + alpha_i * temp_r0;
    y_ptr += inc_y;
    y_ptr[0] += alpha_r * temp_r1 - alpha_i * temp_i1;
    y_ptr[1] += alpha_r * temp_i1 + alpha_i * temp_r1;
#else
    y_ptr[0] += alpha_r * temp_r0 + alpha_i * temp_i0;
    y_ptr[1] -= alpha_r * temp_i0 - alpha_i * temp_r0;
    y_ptr += inc_y;
    y_ptr[0] += alpha_r * temp_r1 + alpha_i * temp_i1;
    y_ptr[1] -= alpha_r * temp_i1 - alpha_i * temp_r1;
#endif
    return (0);
  }

  if (m3 == 3) {
    a_ptr = a;
    x_ptr = x;
    FLOAT temp_r0 = 0.0;
    FLOAT temp_i0 = 0.0;
    FLOAT temp_r1 = 0.0;
    FLOAT temp_i1 = 0.0;
    FLOAT temp_r2 = 0.0;
    FLOAT temp_i2 = 0.0;

    if (lda == 6 && inc_x == 2) {

      for (i = 0; i < n; i++) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        temp_r0 += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] - a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] + a_ptr[3] * x_ptr[0];
        temp_r2 += a_ptr[4] * x_ptr[0] - a_ptr[5] * x_ptr[1];
        temp_i2 += a_ptr[4] * x_ptr[1] + a_ptr[5] * x_ptr[0];
#else
        temp_r0 += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] + a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] - a_ptr[3] * x_ptr[0];
        temp_r2 += a_ptr[4] * x_ptr[0] + a_ptr[5] * x_ptr[1];
        temp_i2 += a_ptr[4] * x_ptr[1] - a_ptr[5] * x_ptr[0];
#endif

        a_ptr += 6;
        x_ptr += 2;
      }

    } else {

      for (i = 0; i < n; i++) {
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        temp_r0 += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] - a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] + a_ptr[3] * x_ptr[0];
        temp_r2 += a_ptr[4] * x_ptr[0] - a_ptr[5] * x_ptr[1];
        temp_i2 += a_ptr[4] * x_ptr[1] + a_ptr[5] * x_ptr[0];
#else
        temp_r0 += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
        temp_i0 += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
        temp_r1 += a_ptr[2] * x_ptr[0] + a_ptr[3] * x_ptr[1];
        temp_i1 += a_ptr[2] * x_ptr[1] - a_ptr[3] * x_ptr[0];
        temp_r2 += a_ptr[4] * x_ptr[0] + a_ptr[5] * x_ptr[1];
        temp_i2 += a_ptr[4] * x_ptr[1] - a_ptr[5] * x_ptr[0];
#endif

        a_ptr += lda;
        x_ptr += inc_x;
      }

    }
#if !defined(XCONJ)
    y_ptr[0] += alpha_r * temp_r0 - alpha_i * temp_i0;
    y_ptr[1] += alpha_r * temp_i0 + alpha_i * temp_r0;
    y_ptr += inc_y;
    y_ptr[0] += alpha_r * temp_r1 - alpha_i * temp_i1;
    y_ptr[1] += alpha_r * temp_i1 + alpha_i * temp_r1;
    y_ptr += inc_y;
    y_ptr[0] += alpha_r * temp_r2 - alpha_i * temp_i2;
    y_ptr[1] += alpha_r * temp_i2 + alpha_i * temp_r2;
#else
    y_ptr[0] += alpha_r * temp_r0 + alpha_i * temp_i0;
    y_ptr[1] -= alpha_r * temp_i0 - alpha_i * temp_r0;
    y_ptr += inc_y;
    y_ptr[0] += alpha_r * temp_r1 + alpha_i * temp_i1;
    y_ptr[1] -= alpha_r * temp_i1 - alpha_i * temp_r1;
    y_ptr += inc_y;
    y_ptr[0] += alpha_r * temp_r2 + alpha_i * temp_i2;
    y_ptr[1] -= alpha_r * temp_i2 - alpha_i * temp_r2;
#endif
    return (0);
  }

  return (0);
}
