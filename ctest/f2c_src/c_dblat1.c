/* c_dblat1.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer icase, n, incx, incy, mode;
    logical pass;
} combla_;

#define combla_1 combla_

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static doublereal c_b34 = 1.;

/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static doublereal sfac = 9.765625e-4;

    /* Format strings */
    static char fmt_99999[] = "(\002 Real CBLAS Test Program Results\002,/1x)"
	    ;
    static char fmt_99998[] = "(\002                                    ----"
	    "- PASS -----\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int header_(void), check0_(doublereal *), check1_(
	    doublereal *), check2_(doublereal *), check3_(doublereal *);
    static integer ic;

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 6, 0, fmt_99999, 0 };
    static cilist io___4 = { 0, 6, 0, fmt_99998, 0 };


/*     Test program for the DOUBLE PRECISION Level 1 CBLAS. */
/*     Based upon the original CBLAS test routine together with: */
/*     F06EAF Example Program Text */
/*     .. Parameters .. */
/*     .. Scalars in Common .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Common blocks .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */
    s_wsfe(&io___2);
    e_wsfe();
    for (ic = 1; ic <= 11; ++ic) {
	combla_1.icase = ic;
	header_();

/*        .. Initialize  PASS,  INCX,  INCY, and MODE for a new case. .. */
/*        .. the value 9999 for INCX, INCY or MODE will appear in the .. */
/*        .. detailed  output, if any, for cases  that do not involve .. */
/*        .. these parameters .. */

	combla_1.pass = TRUE_;
	combla_1.incx = 9999;
	combla_1.incy = 9999;
	combla_1.mode = 9999;
	if (combla_1.icase == 3) {
	    check0_(&sfac);
	} else if (combla_1.icase == 7 || combla_1.icase == 8 || 
		combla_1.icase == 9 || combla_1.icase == 10) {
	    check1_(&sfac);
	} else if (combla_1.icase == 1 || combla_1.icase == 2 || 
		combla_1.icase == 5 || combla_1.icase == 6) {
	    check2_(&sfac);
	} else if (combla_1.icase == 4 || combla_1.icase == 11) {
	    check3_(&sfac);
	}
/*        -- Print */
	if (combla_1.pass) {
	    s_wsfe(&io___4);
	    e_wsfe();
	}
/* L20: */
    }
    s_stop("", (ftnlen)0);

    return 0;
} /* MAIN__ */

/* Subroutine */ int header_(void)
{
    /* Initialized data */

    static char l[15*11] = "CBLAS_DDOT     " "CBLAS_DAXPY    " "CBLAS_DROTG "
	    "   " "CBLAS_DROT     " "CBLAS_DCOPY    " "CBLAS_DSWAP    " "CBLA"
	    "S_DNRM2    " "CBLAS_DASUM    " "CBLAS_DSCAL    " "CBLAS_IDAMAX   "
	     "CBLAS_DROTM    ";

    /* Format strings */
    static char fmt_99999[] = "(/\002 Test of subprogram number\002,i3,9x,a1"
	    "5)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 6, 0, fmt_99999, 0 };


/*     .. Parameters .. */
/*     .. Scalars in Common .. */
/*     .. Local Arrays .. */
/*     .. Common blocks .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */
    s_wsfe(&io___6);
    do_fio(&c__1, (char *)&combla_1.icase, (ftnlen)sizeof(integer));
    do_fio(&c__1, l + (0 + (0 + (combla_1.icase - 1) * 15)), (ftnlen)15);
    e_wsfe();
    return 0;

} /* header_ */

/* Subroutine */ int check0_(doublereal *sfac)
{
    /* Initialized data */

    static doublereal ds1[8] = { .8,.6,.8,-.6,.8,0.,1.,0. };
    static doublereal datrue[8] = { .5,.5,.5,-.5,-.5,0.,1.,1. };
    static doublereal dbtrue[8] = { 0.,.6,0.,-.6,0.,0.,1.,0. };
    static doublereal da1[8] = { .3,.4,-.3,-.4,-.3,0.,0.,1. };
    static doublereal db1[8] = { .4,.3,.4,.3,-.4,0.,1.,0. };
    static doublereal dc1[8] = { .6,.8,-.6,.8,.6,1.,0.,1. };

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int stest1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal sa, sb, sc, ss;
    extern /* Subroutine */ int drotgtest_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___18 = { 0, 6, 0, 0, 0 };


/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Scalars in Common .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Common blocks .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */

/*     Compute true values which cannot be prestored */
/*     in decimal notation */

    dbtrue[0] = 1.6666666666666667;
    dbtrue[2] = -1.6666666666666667;
    dbtrue[4] = 1.6666666666666667;

    for (k = 1; k <= 8; ++k) {
/*        .. Set N=K for identification in output if any .. */
	combla_1.n = k;
	if (combla_1.icase == 3) {
/*           .. DROTGTEST .. */
	    if (k > 8) {
		goto L40;
	    }
	    sa = da1[k - 1];
	    sb = db1[k - 1];
	    drotgtest_(&sa, &sb, &sc, &ss);
	    stest1_(&sa, &datrue[k - 1], &datrue[k - 1], sfac);
	    stest1_(&sb, &dbtrue[k - 1], &dbtrue[k - 1], sfac);
	    stest1_(&sc, &dc1[k - 1], &dc1[k - 1], sfac);
	    stest1_(&ss, &ds1[k - 1], &ds1[k - 1], sfac);
	} else {
	    s_wsle(&io___18);
	    do_lio(&c__9, &c__1, " Shouldn't be here in CHECK0", (ftnlen)28);
	    e_wsle();
	    s_stop("", (ftnlen)0);
	}
/* L20: */
    }
L40:
    return 0;
} /* check0_ */

/* Subroutine */ int check1_(doublereal *sfac)
{
    /* Initialized data */

    static doublereal sa[10] = { .3,-1.,0.,1.,.3,.3,.3,.3,.3,.3 };
    static doublereal dv[80]	/* was [8][5][2] */ = { .1,2.,2.,2.,2.,2.,2.,
	    2.,.3,3.,3.,3.,3.,3.,3.,3.,.3,-.4,4.,4.,4.,4.,4.,4.,.2,-.6,.3,5.,
	    5.,5.,5.,5.,.1,-.3,.5,-.1,6.,6.,6.,6.,.1,8.,8.,8.,8.,8.,8.,8.,.3,
	    9.,9.,9.,9.,9.,9.,9.,.3,2.,-.4,2.,2.,2.,2.,2.,.2,3.,-.6,5.,.3,2.,
	    2.,2.,.1,4.,-.3,6.,-.5,7.,-.1,3. };
    static doublereal dtrue1[5] = { 0.,.3,.5,.7,.6 };
    static doublereal dtrue3[5] = { 0.,.3,.7,1.1,1. };
    static doublereal dtrue5[80]	/* was [8][5][2] */ = { .1,2.,2.,2.,
	    2.,2.,2.,2.,-.3,3.,3.,3.,3.,3.,3.,3.,0.,0.,4.,4.,4.,4.,4.,4.,.2,
	    -.6,.3,5.,5.,5.,5.,5.,.03,-.09,.15,-.03,6.,6.,6.,6.,.1,8.,8.,8.,
	    8.,8.,8.,8.,.09,9.,9.,9.,9.,9.,9.,9.,.09,2.,-.12,2.,2.,2.,2.,2.,
	    .06,3.,-.18,5.,.09,2.,2.,2.,.03,4.,-.09,6.,-.15,7.,-.03,3. };
    static integer itrue2[5] = { 0,1,2,2,3 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int itest1_(integer *, integer *), stest1_(
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal sx[8];
    static integer len;
    extern /* Subroutine */ int dscaltest_(integer *, doublereal *, 
	    doublereal *, integer *);
    static integer np1;
    extern doublereal dasumtest_(integer *, doublereal *, integer *), 
	    dnrm2test_(integer *, doublereal *, integer *);
    extern integer idamaxtest_(integer *, doublereal *, integer *);
    static doublereal stemp[1], strue[8];
    extern /* Subroutine */ int stest_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___31 = { 0, 6, 0, 0, 0 };


/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Scalars in Common .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Common blocks .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */
    for (combla_1.incx = 1; combla_1.incx <= 2; ++combla_1.incx) {
	for (np1 = 1; np1 <= 5; ++np1) {
	    combla_1.n = np1 - 1;
	    len = max(combla_1.n,1) << 1;
/*           .. Set vector arguments .. */
	    i__1 = len;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		sx[i__ - 1] = dv[i__ + (np1 + combla_1.incx * 5 << 3) - 49];
/* L20: */
	    }

	    if (combla_1.icase == 7) {
/*              .. DNRM2TEST .. */
		stemp[0] = dtrue1[np1 - 1];
		d__1 = dnrm2test_(&combla_1.n, sx, &combla_1.incx);
		stest1_(&d__1, stemp, stemp, sfac);
	    } else if (combla_1.icase == 8) {
/*              .. DASUMTEST .. */
		stemp[0] = dtrue3[np1 - 1];
		d__1 = dasumtest_(&combla_1.n, sx, &combla_1.incx);
		stest1_(&d__1, stemp, stemp, sfac);
	    } else if (combla_1.icase == 9) {
/*              .. DSCALTEST .. */
		dscaltest_(&combla_1.n, &sa[(combla_1.incx - 1) * 5 + np1 - 1]
			, sx, &combla_1.incx);
		i__1 = len;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    strue[i__ - 1] = dtrue5[i__ + (np1 + combla_1.incx * 5 << 
			    3) - 49];
/* L40: */
		}
		stest_(&len, sx, strue, strue, sfac);
	    } else if (combla_1.icase == 10) {
/*              .. IDAMAXTEST .. */
		i__1 = idamaxtest_(&combla_1.n, sx, &combla_1.incx);
		itest1_(&i__1, &itrue2[np1 - 1]);
	    } else {
		s_wsle(&io___31);
		do_lio(&c__9, &c__1, " Shouldn't be here in CHECK1", (ftnlen)
			28);
		e_wsle();
		s_stop("", (ftnlen)0);
	    }
/* L60: */
	}
/* L80: */
    }
    return 0;
} /* check1_ */

/* Subroutine */ int check2_(doublereal *sfac)
{
    /* Initialized data */

    static doublereal sa = .3;
    static integer incxs[4] = { 1,2,-2,-1 };
    static integer incys[4] = { 1,-2,1,-2 };
    static integer lens[8]	/* was [4][2] */ = { 1,1,2,4,1,1,3,7 };
    static integer ns[4] = { 0,1,2,4 };
    static doublereal dx1[7] = { .6,.1,-.5,.8,.9,-.3,-.4 };
    static doublereal dy1[7] = { .5,-.9,.3,.7,-.6,.2,.8 };
    static doublereal dt7[16]	/* was [4][4] */ = { 0.,.3,.21,.62,0.,.3,-.07,
	    .85,0.,.3,-.79,-.74,0.,.3,.33,1.27 };
    static doublereal dt8[112]	/* was [7][4][4] */ = { .5,0.,0.,0.,0.,0.,0.,
	    .68,0.,0.,0.,0.,0.,0.,.68,-.87,0.,0.,0.,0.,0.,.68,-.87,.15,.94,0.,
	    0.,0.,.5,0.,0.,0.,0.,0.,0.,.68,0.,0.,0.,0.,0.,0.,.35,-.9,.48,0.,
	    0.,0.,0.,.38,-.9,.57,.7,-.75,.2,.98,.5,0.,0.,0.,0.,0.,0.,.68,0.,
	    0.,0.,0.,0.,0.,.35,-.72,0.,0.,0.,0.,0.,.38,-.63,.15,.88,0.,0.,0.,
	    .5,0.,0.,0.,0.,0.,0.,.68,0.,0.,0.,0.,0.,0.,.68,-.9,.33,0.,0.,0.,
	    0.,.68,-.9,.33,.7,-.75,.2,1.04 };
    static doublereal dt10x[112]	/* was [7][4][4] */ = { .6,0.,0.,0.,
	    0.,0.,0.,.5,0.,0.,0.,0.,0.,0.,.5,-.9,0.,0.,0.,0.,0.,.5,-.9,.3,.7,
	    0.,0.,0.,.6,0.,0.,0.,0.,0.,0.,.5,0.,0.,0.,0.,0.,0.,.3,.1,.5,0.,0.,
	    0.,0.,.8,.1,-.6,.8,.3,-.3,.5,.6,0.,0.,0.,0.,0.,0.,.5,0.,0.,0.,0.,
	    0.,0.,-.9,.1,.5,0.,0.,0.,0.,.7,.1,.3,.8,-.9,-.3,.5,.6,0.,0.,0.,0.,
	    0.,0.,.5,0.,0.,0.,0.,0.,0.,.5,.3,0.,0.,0.,0.,0.,.5,.3,-.6,.8,0.,
	    0.,0. };
    static doublereal dt10y[112]	/* was [7][4][4] */ = { .5,0.,0.,0.,
	    0.,0.,0.,.6,0.,0.,0.,0.,0.,0.,.6,.1,0.,0.,0.,0.,0.,.6,.1,-.5,.8,
	    0.,0.,0.,.5,0.,0.,0.,0.,0.,0.,.6,0.,0.,0.,0.,0.,0.,-.5,-.9,.6,0.,
	    0.,0.,0.,-.4,-.9,.9,.7,-.5,.2,.6,.5,0.,0.,0.,0.,0.,0.,.6,0.,0.,0.,
	    0.,0.,0.,-.5,.6,0.,0.,0.,0.,0.,-.4,.9,-.5,.6,0.,0.,0.,.5,0.,0.,0.,
	    0.,0.,0.,.6,0.,0.,0.,0.,0.,0.,.6,-.9,.1,0.,0.,0.,0.,.6,-.9,.1,.7,
	    -.5,.2,.8 };
    static doublereal ssize1[4] = { 0.,.3,1.6,3.2 };
    static doublereal ssize2[28]	/* was [14][2] */ = { 0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,1.17,1.17,1.17,1.17,1.17,1.17,1.17,
	    1.17,1.17,1.17,1.17,1.17,1.17,1.17 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int stest1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ki, kn, mx, my;
    static doublereal sx[7], sy[7];
    extern doublereal ddottest_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal stx[7], sty[7];
    extern /* Subroutine */ int dcopytest_(integer *, doublereal *, integer *,
	     doublereal *, integer *), dswaptest_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), daxpytest_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer lenx, leny, ksize;
    extern /* Subroutine */ int stest_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 6, 0, 0, 0 };


/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Scalars in Common .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Common blocks .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */

    for (ki = 1; ki <= 4; ++ki) {
	combla_1.incx = incxs[ki - 1];
	combla_1.incy = incys[ki - 1];
	mx = abs(combla_1.incx);
	my = abs(combla_1.incy);

	for (kn = 1; kn <= 4; ++kn) {
	    combla_1.n = ns[kn - 1];
	    ksize = min(2,kn);
	    lenx = lens[kn + (mx << 2) - 5];
	    leny = lens[kn + (my << 2) - 5];
/*           .. Initialize all argument arrays .. */
	    for (i__ = 1; i__ <= 7; ++i__) {
		sx[i__ - 1] = dx1[i__ - 1];
		sy[i__ - 1] = dy1[i__ - 1];
/* L20: */
	    }

	    if (combla_1.icase == 1) {
/*              .. DDOTTEST .. */
		d__1 = ddottest_(&combla_1.n, sx, &combla_1.incx, sy, &
			combla_1.incy);
		stest1_(&d__1, &dt7[kn + (ki << 2) - 5], &ssize1[kn - 1], 
			sfac);
	    } else if (combla_1.icase == 2) {
/*              .. DAXPYTEST .. */
		daxpytest_(&combla_1.n, &sa, sx, &combla_1.incx, sy, &
			combla_1.incy);
		i__1 = leny;
		for (j = 1; j <= i__1; ++j) {
		    sty[j - 1] = dt8[j + (kn + (ki << 2)) * 7 - 36];
/* L40: */
		}
		stest_(&leny, sy, sty, &ssize2[ksize * 14 - 14], sfac);
	    } else if (combla_1.icase == 5) {
/*              .. DCOPYTEST .. */
		for (i__ = 1; i__ <= 7; ++i__) {
		    sty[i__ - 1] = dt10y[i__ + (kn + (ki << 2)) * 7 - 36];
/* L60: */
		}
		dcopytest_(&combla_1.n, sx, &combla_1.incx, sy, &
			combla_1.incy);
		stest_(&leny, sy, sty, ssize2, &c_b34);
	    } else if (combla_1.icase == 6) {
/*              .. DSWAPTEST .. */
		dswaptest_(&combla_1.n, sx, &combla_1.incx, sy, &
			combla_1.incy);
		for (i__ = 1; i__ <= 7; ++i__) {
		    stx[i__ - 1] = dt10x[i__ + (kn + (ki << 2)) * 7 - 36];
		    sty[i__ - 1] = dt10y[i__ + (kn + (ki << 2)) * 7 - 36];
/* L80: */
		}
		stest_(&lenx, sx, stx, ssize2, &c_b34);
		stest_(&leny, sy, sty, ssize2, &c_b34);
	    } else {
		s_wsle(&io___58);
		do_lio(&c__9, &c__1, " Shouldn't be here in CHECK2", (ftnlen)
			28);
		e_wsle();
		s_stop("", (ftnlen)0);
	    }
/* L100: */
	}
/* L120: */
    }
    return 0;
} /* check2_ */

/* Subroutine */ int check3_(doublereal *sfac)
{
    /* Initialized data */

    static integer incxs[7] = { 1,1,2,2,-2,-1,-2 };
    static integer incys[7] = { 1,2,2,-2,1,-2,-2 };
    static integer ns[5] = { 0,1,2,4,5 };
    static doublereal dx[10] = { .6,.1,-.5,.8,.9,-.3,-.4,.7,.5,.2 };
    static doublereal dy[10] = { .5,-.9,.3,.7,-.6,.2,.8,-.5,.1,-.3 };
    static doublereal sc = .8;
    static doublereal ss = .6;
    static integer len = 10;
    static doublereal param[20]	/* was [5][4] */ = { -2.,1.,0.,0.,1.,-1.,.2,
	    .3,.4,.5,0.,1.,.3,.4,1.,1.,.2,-1.,1.,.5 };
    static doublereal ssize2[20]	/* was [10][2] */ = { 0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,1.17,1.17,1.17,1.17,1.17,1.17,1.17,1.17,1.17,1.17 }
	    ;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static doublereal dparam[5];
    static integer i__, k, ki, kn;
    static doublereal sx[10], sy[10];
    extern /* Subroutine */ int drottest_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal stx[10], sty[10];
    extern /* Subroutine */ int drotmtest_(integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *), drot_(integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *);
    static integer ksize;
    extern /* Subroutine */ int drotm_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *), stest_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___79 = { 0, 6, 0, 0, 0 };


/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Scalars in Common .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Common blocks .. */
/*     .. Data statements .. */
/*     .. Executable Statements .. */

    for (ki = 1; ki <= 7; ++ki) {
	combla_1.incx = incxs[ki - 1];
	combla_1.incy = incys[ki - 1];

	for (kn = 1; kn <= 5; ++kn) {
	    combla_1.n = ns[kn - 1];
	    ksize = min(2,kn);

	    if (combla_1.icase == 4) {
/*              .. DROTTEST .. */
		for (i__ = 1; i__ <= 10; ++i__) {
		    sx[i__ - 1] = dx[i__ - 1];
		    sy[i__ - 1] = dy[i__ - 1];
		    stx[i__ - 1] = dx[i__ - 1];
		    sty[i__ - 1] = dy[i__ - 1];
/* L20: */
		}
		drottest_(&combla_1.n, sx, &combla_1.incx, sy, &combla_1.incy,
			 &sc, &ss);
		drot_(&combla_1.n, stx, &combla_1.incx, sty, &combla_1.incy, &
			sc, &ss);
		stest_(&len, sx, stx, &ssize2[ksize * 10 - 10], sfac);
		stest_(&len, sy, sty, &ssize2[ksize * 10 - 10], sfac);
	    } else if (combla_1.icase == 11) {
/*              .. DROTMTEST .. */
		for (i__ = 1; i__ <= 10; ++i__) {
		    sx[i__ - 1] = dx[i__ - 1];
		    sy[i__ - 1] = dy[i__ - 1];
		    stx[i__ - 1] = dx[i__ - 1];
		    sty[i__ - 1] = dy[i__ - 1];
/* L90: */
		}
		for (i__ = 1; i__ <= 4; ++i__) {
		    for (k = 1; k <= 5; ++k) {
			dparam[k - 1] = param[k + i__ * 5 - 6];
/* L80: */
		    }
		    drotmtest_(&combla_1.n, sx, &combla_1.incx, sy, &
			    combla_1.incy, dparam);
		    drotm_(&combla_1.n, stx, &combla_1.incx, sty, &
			    combla_1.incy, dparam);
		    stest_(&len, sx, stx, &ssize2[ksize * 10 - 10], sfac);
		    stest_(&len, sy, sty, &ssize2[ksize * 10 - 10], sfac);
/* L70: */
		}
	    } else {
		s_wsle(&io___79);
		do_lio(&c__9, &c__1, " Shouldn't be here in CHECK3", (ftnlen)
			28);
		e_wsle();
		s_stop("", (ftnlen)0);
	    }
/* L40: */
	}
/* L60: */
    }
    return 0;
} /* check3_ */

/* Subroutine */ int stest_(integer *len, doublereal *scomp, doublereal *
	strue, doublereal *ssize, doublereal *sfac)
{
    /* Format strings */
    static char fmt_99999[] = "(\002                                       F"
	    "AIL\002)";
    static char fmt_99998[] = "(/\002 CASE  N INCX INCY MODE  I             "
	    "               \002,\002 COMP(I)                             TRU"
	    "E(I)  DIFFERENCE\002,\002     SIZE(I)\002,/1x)";
    static char fmt_99997[] = "(1x,i4,i3,3i5,i3,2d36.8,2d12.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal sd;
    extern doublereal sdiff_(doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___82 = { 0, 6, 0, fmt_99999, 0 };
    static cilist io___83 = { 0, 6, 0, fmt_99998, 0 };
    static cilist io___84 = { 0, 6, 0, fmt_99997, 0 };


/*     ********************************* STEST ************************** */

/*     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO */
/*     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE */
/*     NEGLIGIBLE. */

/*     C. L. LAWSON, JPL, 1974 DEC 10 */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Scalars in Common .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Common blocks .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ssize;
    --strue;
    --scomp;

    /* Function Body */
    i__1 = *len;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sd = scomp[i__] - strue[i__];
	d__4 = (d__1 = ssize[i__], abs(d__1)) + (d__2 = *sfac * sd, abs(d__2))
		;
	d__5 = (d__3 = ssize[i__], abs(d__3));
	if (sdiff_(&d__4, &d__5) == 0.) {
	    goto L40;
	}

/*                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I). */

	if (! combla_1.pass) {
	    goto L20;
	}
/*                             PRINT FAIL MESSAGE AND HEADER. */
	combla_1.pass = FALSE_;
	s_wsfe(&io___82);
	e_wsfe();
	s_wsfe(&io___83);
	e_wsfe();
L20:
	s_wsfe(&io___84);
	do_fio(&c__1, (char *)&combla_1.icase, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&combla_1.n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&combla_1.incx, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&combla_1.incy, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&combla_1.mode, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&scomp[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&strue[i__], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ssize[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
L40:
	;
    }
    return 0;

} /* stest_ */

/* Subroutine */ int stest1_(doublereal *scomp1, doublereal *strue1, 
	doublereal *ssize, doublereal *sfac)
{
    static doublereal scomp[1], strue[1];
    extern /* Subroutine */ int stest_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*     ************************* STEST1 ***************************** */

/*     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN */
/*     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE */
/*     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT. */

/*     C.L. LAWSON, JPL, 1978 DEC 6 */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ssize;

    /* Function Body */
    scomp[0] = *scomp1;
    strue[0] = *strue1;
    stest_(&c__1, scomp, strue, &ssize[1], sfac);

    return 0;
} /* stest1_ */

doublereal sdiff_(doublereal *sa, doublereal *sb)
{
    /* System generated locals */
    doublereal ret_val;

/*     ********************************* SDIFF ************************** */
/*     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15 */

/*     .. Scalar Arguments .. */
/*     .. Executable Statements .. */
    ret_val = *sa - *sb;
    return ret_val;
} /* sdiff_ */

/* Subroutine */ int itest1_(integer *icomp, integer *itrue)
{
    /* Format strings */
    static char fmt_99999[] = "(\002                                       F"
	    "AIL\002)";
    static char fmt_99998[] = "(/\002 CASE  N INCX INCY MODE                "
	    "               \002,\002 COMP                                TRU"
	    "E     DIFFERENCE\002,/1x)";
    static char fmt_99997[] = "(1x,i4,i3,3i5,2i36,i12)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer id;

    /* Fortran I/O blocks */
    static cilist io___87 = { 0, 6, 0, fmt_99999, 0 };
    static cilist io___88 = { 0, 6, 0, fmt_99998, 0 };
    static cilist io___90 = { 0, 6, 0, fmt_99997, 0 };


/*     ********************************* ITEST1 ************************* */

/*     THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR */
/*     EQUALITY. */
/*     C. L. LAWSON, JPL, 1974 DEC 10 */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Scalars in Common .. */
/*     .. Local Scalars .. */
/*     .. Common blocks .. */
/*     .. Executable Statements .. */

    if (*icomp == *itrue) {
	goto L40;
    }

/*                            HERE ICOMP IS NOT EQUAL TO ITRUE. */

    if (! combla_1.pass) {
	goto L20;
    }
/*                             PRINT FAIL MESSAGE AND HEADER. */
    combla_1.pass = FALSE_;
    s_wsfe(&io___87);
    e_wsfe();
    s_wsfe(&io___88);
    e_wsfe();
L20:
    id = *icomp - *itrue;
    s_wsfe(&io___90);
    do_fio(&c__1, (char *)&combla_1.icase, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&combla_1.n, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&combla_1.incx, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&combla_1.incy, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&combla_1.mode, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*icomp), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*itrue), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
    e_wsfe();
L40:
    return 0;

} /* itest1_ */

/* Subroutine */ int drot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ix, iy;
    static doublereal dtemp;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     applies a plane rotation. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*     .. Local Scalars .. */
/*     .. */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[ix] + *s * dy[iy];
	dy[iy] = *c__ * dy[iy] - *s * dx[ix];
	dx[ix] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[i__] + *s * dy[i__];
	dy[i__] = *c__ * dy[i__] - *s * dx[i__];
	dx[i__] = dtemp;
/* L30: */
    }
    return 0;
} /* drot_ */

/* Subroutine */ int drotm_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer nsteps, i__;
    static doublereal w, z__;
    static integer kx, ky;
    static doublereal dflag, dh11, dh12, dh21, dh22;


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Data statements .. */
    /* Parameter adjustments */
    --dparam;
    --dy;
    --dx;

    /* Function Body */
/*     .. */

    dflag = dparam[1];
    if (*n <= 0 || dflag + two == zero) {
	return 0;
    }
    if (*incx == *incy && *incx > 0) {

	nsteps = *n * *incx;
	if (dflag < zero) {
	    dh11 = dparam[2];
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    dh22 = dparam[5];
	    i__1 = nsteps;
	    i__2 = *incx;
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		w = dx[i__];
		z__ = dy[i__];
		dx[i__] = w * dh11 + z__ * dh12;
		dy[i__] = w * dh21 + z__ * dh22;
	    }
	} else if (dflag == zero) {
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    i__2 = nsteps;
	    i__1 = *incx;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
		w = dx[i__];
		z__ = dy[i__];
		dx[i__] = w + z__ * dh12;
		dy[i__] = w * dh21 + z__;
	    }
	} else {
	    dh11 = dparam[2];
	    dh22 = dparam[5];
	    i__1 = nsteps;
	    i__2 = *incx;
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		w = dx[i__];
		z__ = dy[i__];
		dx[i__] = w * dh11 + z__;
		dy[i__] = -w + dh22 * z__;
	    }
	}
    } else {
	kx = 1;
	ky = 1;
	if (*incx < 0) {
	    kx = (1 - *n) * *incx + 1;
	}
	if (*incy < 0) {
	    ky = (1 - *n) * *incy + 1;
	}

	if (dflag < zero) {
	    dh11 = dparam[2];
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    dh22 = dparam[5];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w = dx[kx];
		z__ = dy[ky];
		dx[kx] = w * dh11 + z__ * dh12;
		dy[ky] = w * dh21 + z__ * dh22;
		kx += *incx;
		ky += *incy;
	    }
	} else if (dflag == zero) {
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w = dx[kx];
		z__ = dy[ky];
		dx[kx] = w + z__ * dh12;
		dy[ky] = w * dh21 + z__;
		kx += *incx;
		ky += *incy;
	    }
	} else {
	    dh11 = dparam[2];
	    dh22 = dparam[5];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w = dx[kx];
		z__ = dy[ky];
		dx[kx] = w * dh11 + z__;
		dy[ky] = -w + dh22 * z__;
		kx += *incx;
		ky += *incy;
	    }
	}
    }
    return 0;
} /* drotm_ */

/* Main program alias */ int dcblat1_ () { MAIN__ (); return 0; }
