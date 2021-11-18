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

/*
 * Avoid contraction of floating point operations, specifically fused
 * multiply-add, because they can cause unexpected results in complex
 * multiplication.
 */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif

#if defined(__clang__)
#pragma clang fp contract(off)
#endif

#include "common.h"
#include "vector-common.h"

static void cscal_kernel_16(BLASLONG n, FLOAT da_r, FLOAT da_i, FLOAT *x) {
    vector_float da_r_vec = vec_splats(da_r);
    vector_float da_i_vec = { -da_i, da_i, -da_i, da_i };

    vector_float *x_vec_ptr = (vector_float *)x;

#pragma GCC unroll 16
    for (size_t i = 0; i < n/2; i++) {
	vector_float x_vec = vec_load_hinted(x + i * VLEN_FLOATS);
	vector_float x_swapped = {x_vec[1], x_vec[0], x_vec[3], x_vec[2]};

	x_vec_ptr[i] =  x_vec * da_r_vec + x_swapped * da_i_vec;
    }
}

static void cscal_kernel_16_zero_r(BLASLONG n, FLOAT *alpha, FLOAT *x) {
__asm__ (
" vlef   0,4(%[alpha]),0\n"
" vlef   0,4(%[alpha]),2\n"
" vflcsb 0,0\n"
" vlef   0,4(%[alpha]),1\n"
" vlef   0,4(%[alpha]),3\n"
" srlg %[n],%[n],4\n"
" xgr   1,1\n"
"label0 pfd 2,1024(1,%[x])\n"
" vl   16,0(1,%[x])\n"
" vl   17,16(1,%[x])\n"
" vl   18,32(1,%[x])\n"
" vl   19,48(1,%[x])\n"
" vl   20,64(1,%[x])\n"
" vl   21,80(1,%[x])\n"
" vl   22,96(1,%[x])\n"
" vl   23,112(1,%[x])\n"
" verllg   16,16,32\n"
" verllg   17,17,32\n"
" verllg   18,18,32\n"
" verllg   19,19,32\n"
" verllg   20,20,32\n"
" verllg   21,21,32\n"
" verllg   22,22,32\n"
" verllg   23,23,32\n"
" vfmsb 16,16,0\n"
" vfmsb 17,17,0\n"
" vfmsb 18,18,0\n"
" vfmsb 19,19,0\n"
" vfmsb 20,20,0\n"
" vfmsb 21,21,0\n"
" vfmsb 22,22,0\n"
" vfmsb 23,23,0\n"
" vst 16,0(1,%[x])\n"
" vst 17,16(1,%[x])\n"
" vst 18,32(1,%[x])\n"
" vst 19,48(1,%[x])\n"
" vst 20,64(1,%[x])\n"
" vst 21,80(1,%[x])\n"
" vst 22,96(1,%[x])\n"
" vst 23,112(1,%[x])\n"
" agfi  1,128\n"
" brctg %[n],label0\n"
    : "+m"(*(FLOAT (*)[n * 2]) x),[n] "+&r"(n)
    : [x] "a"(x), "m"(*(const FLOAT (*)[2]) alpha),
       [alpha] "a"(alpha)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23");
}

static void cscal_kernel_16_zero_i(BLASLONG n, FLOAT *alpha, FLOAT *x) {
__asm__ (
" vlrepf 0,0(%[alpha])\n"
" srlg %[n],%[n],4\n"
" xgr   1,1\n"
"label1 pfd 2,1024(1,%[x])\n"
" vl   16,0(1,%[x])\n"
" vl   17,16(1,%[x])\n"
" vl   18,32(1,%[x])\n"
" vl   19,48(1,%[x])\n"
" vl   20,64(1,%[x])\n"
" vl   21,80(1,%[x])\n"
" vl   22,96(1,%[x])\n"
" vl   23,112(1,%[x])\n"
" vfmsb 16,16,0\n"
" vfmsb 17,17,0\n"
" vfmsb 18,18,0\n"
" vfmsb 19,19,0\n"
" vfmsb 20,20,0\n"
" vfmsb 21,21,0\n"
" vfmsb 22,22,0\n"
" vfmsb 23,23,0\n"
" vst 16,0(1,%[x])\n"
" vst 17,16(1,%[x])\n"
" vst 18,32(1,%[x])\n"
" vst 19,48(1,%[x])\n"
" vst 20,64(1,%[x])\n"
" vst 21,80(1,%[x])\n"
" vst 22,96(1,%[x])\n"
" vst 23,112(1,%[x])\n"
" agfi  1,128\n"
" brctg %[n],label1\n"
    : "+m"(*(FLOAT (*)[n * 2]) x),[n] "+&r"(n)
    : [x] "a"(x), "m"(*(const FLOAT (*)[2]) alpha),
       [alpha] "a"(alpha)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23");
}

static void cscal_kernel_16_zero(BLASLONG n, FLOAT *x) {
__asm__ (
" vzero 0\n"
" srlg %[n],%[n],4\n"
" xgr   1,1\n"
"label2 pfd 2,1024(1,%[x])\n"
" vst  0,0(1,%[x])\n"
" vst  0,16(1,%[x])\n"
" vst  0,32(1,%[x])\n"
" vst  0,48(1,%[x])\n"
" vst  0,64(1,%[x])\n"
" vst  0,80(1,%[x])\n"
" vst  0,96(1,%[x])\n"
" vst  0,112(1,%[x])\n"
" agfi  1,128\n"
" brctg %[n],label2\n"
    : "=m"(*(FLOAT (*)[n * 2]) x),[n] "+&r"(n)
    : [x] "a"(x)
    : "cc", "r1", "v0");
}

static void cscal_kernel_inc_8(BLASLONG n, FLOAT da_r, FLOAT da_i, FLOAT *x,
                               BLASLONG inc_x) {
  BLASLONG i;
  BLASLONG inc_x2 = 2 * inc_x;
  BLASLONG inc_x3 = inc_x2 + inc_x;
  FLOAT t0, t1, t2, t3;

  for (i = 0; i < n; i += 4) {
    t0 = da_r * x[0] - da_i * x[1];
    t1 = da_r * x[inc_x] - da_i * x[inc_x + 1];
    t2 = da_r * x[inc_x2] - da_i * x[inc_x2 + 1];
    t3 = da_r * x[inc_x3] - da_i * x[inc_x3 + 1];

    x[1] = da_i * x[0] + da_r * x[1];
    x[inc_x + 1] = da_i * x[inc_x] + da_r * x[inc_x + 1];
    x[inc_x2 + 1] = da_i * x[inc_x2] + da_r * x[inc_x2 + 1];
    x[inc_x3 + 1] = da_i * x[inc_x3] + da_r * x[inc_x3 + 1];

    x[0] = t0;
    x[inc_x] = t1;
    x[inc_x2] = t2;
    x[inc_x3] = t3;

    x += 4 * inc_x;
  }
}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da_r, FLOAT da_i,
          FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy,
          BLASLONG dummy2) {
  BLASLONG i = 0, j = 0;
  FLOAT temp0;
  FLOAT temp1;
  FLOAT alpha[2] __attribute__ ((aligned(16)));

  if (inc_x != 1) {
    inc_x <<= 1;

    if (da_r == 0.0) {

      BLASLONG n1 = n & -2;

      if (da_i == 0.0) {

        while (j < n1) {

          x[i] = 0.0;
          x[i + 1] = 0.0;
          x[i + inc_x] = 0.0;
          x[i + 1 + inc_x] = 0.0;
          i += 2 * inc_x;
          j += 2;

        }

        while (j < n) {

          x[i] = 0.0;
          x[i + 1] = 0.0;
          i += inc_x;
          j++;

        }

      } else {

        while (j < n1) {

          temp0 = -da_i * x[i + 1];
          x[i + 1] = da_i * x[i];
          x[i] = temp0;
          temp1 = -da_i * x[i + 1 + inc_x];
          x[i + 1 + inc_x] = da_i * x[i + inc_x];
          x[i + inc_x] = temp1;
          i += 2 * inc_x;
          j += 2;

        }

        while (j < n) {

          temp0 = -da_i * x[i + 1];
          x[i + 1] = da_i * x[i];
          x[i] = temp0;
          i += inc_x;
          j++;

        }

      }

    } else {

      if (da_i == 0.0) {
        BLASLONG n1 = n & -2;

        while (j < n1) {

          temp0 = da_r * x[i];
          x[i + 1] = da_r * x[i + 1];
          x[i] = temp0;
          temp1 = da_r * x[i + inc_x];
          x[i + 1 + inc_x] = da_r * x[i + 1 + inc_x];
          x[i + inc_x] = temp1;
          i += 2 * inc_x;
          j += 2;

        }

        while (j < n) {

          temp0 = da_r * x[i];
          x[i + 1] = da_r * x[i + 1];
          x[i] = temp0;
          i += inc_x;
          j++;

        }

      } else {

        BLASLONG n1 = n & -8;
        if (n1 > 0) {
          cscal_kernel_inc_8(n1, da_r, da_i, x, inc_x);
          j = n1;
          i = n1 * inc_x;
        }

        while (j < n) {

          temp0 = da_r * x[i] - da_i * x[i + 1];
          x[i + 1] = da_r * x[i + 1] + da_i * x[i];
          x[i] = temp0;
          i += inc_x;
          j++;

        }

      }

    }

    return (0);
  }

  BLASLONG n1 = n & -16;
  if (n1 > 0) {

    alpha[0] = da_r;
    alpha[1] = da_i;

    if (da_r == 0.0)
      if (da_i == 0)
        cscal_kernel_16_zero(n1, x);
      else
        cscal_kernel_16_zero_r(n1, alpha, x);
    else if (da_i == 0)
      cscal_kernel_16_zero_i(n1, alpha, x);
    else
      cscal_kernel_16(n1, da_r, da_i, x);

    i = n1 << 1;
    j = n1;
  }

  if (da_r == 0.0) {

    if (da_i == 0.0) {

      while (j < n) {

        x[i] = 0.0;
        x[i + 1] = 0.0;
        i += 2;
        j++;

      }

    } else {

      while (j < n) {

        temp0 = -da_i * x[i + 1];
        x[i + 1] = da_i * x[i];
        x[i] = temp0;
        i += 2;
        j++;

      }

    }

  } else {

    if (da_i == 0.0) {

      while (j < n) {

        temp0 = da_r * x[i];
        x[i + 1] = da_r * x[i + 1];
        x[i] = temp0;
        i += 2;
        j++;

      }

    } else {

      while (j < n) {

        temp0 = da_r * x[i] - da_i * x[i + 1];
        x[i + 1] = da_r * x[i + 1] + da_i * x[i];
        x[i] = temp0;
        i += 2;
        j++;

      }

    }

  }

  return (0);
}
