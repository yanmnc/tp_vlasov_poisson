      SUBROUTINE DPTTRF( N, D, E, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTTRF computes the L*D*L' factorization of a real symmetric
!  positive definite tridiagonal matrix A.  The factorization may also
!  be regarded as having the form A = U'*D*U.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix
!          A.  On exit, the n diagonal elements of the diagonal matrix
!          D from the L*D*L' factorization of A.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix A.  On exit, the (n-1) subdiagonal elements of the
!          unit bidiagonal factor L from the L*D*L' factorization of A.
!          E can also be regarded as the superdiagonal of the unit
!          bidiagonal factor U from the U'*D*U factorization of A.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite; if k < N, the factorization could not
!               be completed, while if k = N, the factorization was
!               completed, but D(N) = 0.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I4
      DOUBLE PRECISION   EI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         WRITE( *, * )'DPTTRF', -INFO
         STOP
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
        RETURN
!
!     Compute the L*D*L' (or U'*D*U) factorization of A.
!
      I4 = MOD( N-1, 4 )
      DO 10 I = 1, I4
         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
   10 CONTINUE
!
      DO 20 I = I4 + 1, N - 4, 4
!
!        Drop out of the loop if d(i) <= 0: the matrix is not positive
!        definite.
!
         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
!
!        Solve for e(i) and d(i+1).
!
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
!
         IF( D( I+1 ).LE.ZERO ) THEN
            INFO = I + 1
            GO TO 30
         END IF
!
!        Solve for e(i+1) and d(i+2).
!
         EI = E( I+1 )
         E( I+1 ) = EI / D( I+1 )
         D( I+2 ) = D( I+2 ) - E( I+1 )*EI
!
         IF( D( I+2 ).LE.ZERO ) THEN
            INFO = I + 2
            GO TO 30
         END IF
!
!        Solve for e(i+2) and d(i+3).
!
         EI = E( I+2 )
         E( I+2 ) = EI / D( I+2 )
         D( I+3 ) = D( I+3 ) - E( I+2 )*EI
!
         IF( D( I+3 ).LE.ZERO ) THEN
            INFO = I + 3
            GO TO 30
         END IF
!
!        Solve for e(i+3) and d(i+4).
!
         EI = E( I+3 )
         E( I+3 ) = EI / D( I+3 )
         D( I+4 ) = D( I+4 ) - E( I+3 )*EI
   20 CONTINUE
!
!     Check d(n) for positive definiteness.
!
      IF( D( N ).LE.ZERO ) &
        INFO = N
!
   30 CONTINUE
      RETURN
!
!     End of DPTTRF
!
      END


      SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTTRS solves a tridiagonal system of the form
!     A * X = B
!  using the L*D*L' factorization of A computed by DPTTRF.  D is a
!  diagonal matrix specified in the vector D, L is a unit bidiagonal
!  matrix whose subdiagonal is specified in the vector E, and X and B
!  are N by NRHS matrices.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the diagonal matrix D from the
!          L*D*L' factorization of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) subdiagonal elements of the unit bidiagonal factor
!          L from the L*D*L' factorization of A.  E can also be regarded
!          as the superdiagonal of the unit bidiagonal factor U from the
!          factorization A = U'*D*U.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            J, JB, NB
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DPTTS2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         WRITE( *, * )'DPTTRS', -INFO
         STOP
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
        RETURN
!
!     Determine the number of right-hand sides to solve at a time.
!
      IF( NRHS.EQ.1 ) THEN
         NB = 1
      ELSE
         !NB = MAX( 1, ILAENV( 1, 'DPTTRS', ' ', N, NRHS, -1, -1 ) )
         NB = 1
      END IF
!
      IF( NB.GE.NRHS ) THEN
         CALL DPTTS2( N, NRHS, D, E, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL DPTTS2( N, JB, D, E, B( 1, J ), LDB )
   10    CONTINUE
      END IF
!
      RETURN
!
!     End of DPTTRS
!
      END


      SUBROUTINE DPTTS2( N, NRHS, D, E, B, LDB )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTTS2 solves a tridiagonal system of the form
!     A * X = B
!  using the L*D*L' factorization of A computed by DPTTRF.  D is a
!  diagonal matrix specified in the vector D, L is a unit bidiagonal
!  matrix whose subdiagonal is specified in the vector E, and X and B
!  are N by NRHS matrices.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the diagonal matrix D from the
!          L*D*L' factorization of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) subdiagonal elements of the unit bidiagonal factor
!          L from the L*D*L' factorization of A.  E can also be regarded
!          as the superdiagonal of the unit bidiagonal factor U from the
!          factorization A = U'*D*U.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 ) THEN
         IF( N.EQ.1 ) &
           CALL DSCAL( NRHS, 1.D0 / D( 1 ), B, LDB )
         RETURN
      END IF
!
!     Solve A * X = B using the factorization A = L*D*L',
!     overwriting each right hand side vector with its solution.
!
      DO 30 J = 1, NRHS
!
!           Solve L * x = b.
!
         DO 10 I = 2, N
            B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   10    CONTINUE
!
!           Solve D * L' * x = b.
!
         B( N, J ) = B( N, J ) / D( N )
         DO 20 I = N - 1, 1, -1
            B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   20    CONTINUE
   30 CONTINUE
!
      RETURN
!
!     End of DPTTS2
!
      END


      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end


!**************************************************************************
!      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!**************************************************************************
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     January 2007
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  ILAENV returns an INTEGER
!  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines (DEPRECATED)
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR method
!               for nonsymmetric eigenvalue problems (DEPRECATED)
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!          12 <= ISPEC <= 16:
!               xHSEQR or one of its subroutines,
!               see IPARMQ for detailed explanation
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
!     ..
!     .. Executable Statements ..
!
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
             130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
   10 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.   &
            ( IC.GE.145 .AND. IC.LE.153 ) .OR.    &
            ( IC.GE.162 .AND. IC.LE.169 ) ) THEN  
              SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )               
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
        RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
!
      GO TO ( 50, 60, 70 )ISPEC
!
   50 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
   60 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
   70 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
               'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
   80 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
   90 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  100 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  110 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  120 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  130 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
  140 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
!
  150 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
!
  160 CONTINUE
!
!     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
!
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
!
!     End of ILAENV
!
      END


!**************************************************************************
!      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!**************************************************************************
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                         NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 ) RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END


!**************************************************************************
!      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!**************************************************************************
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!     
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
!
!  Purpose
!  =======
!
!       This program sets problem and machine dependent parameters
!       useful for xHSEQR and its subroutines. It is called whenever 
!       ILAENV is called with 12 <= ISPEC <= 16
!
!  Arguments
!  =========
!
!       ISPEC  (input) integer scalar
!              ISPEC specifies which tunable parameter IPARMQ should
!              return.
!
!              ISPEC=12: (INMIN)  Matrices of order nmin or less
!                        are sent directly to xLAHQR, the implicit
!                        double shift QR algorithm.  NMIN must be
!                        at least 11.
!
!              ISPEC=13: (INWIN)  Size of the deflation window.
!                        This is best set greater than or equal to
!                        the number of simultaneous shifts NS.
!                        Larger matrices benefit from larger deflation
!                        windows.
!
!              ISPEC=14: (INIBL) Determines when to stop nibbling and
!                        invest in an (expensive) multi-shift QR sweep.
!                        If the aggressive early deflation subroutine
!                        finds LD converged eigenvalues from an order
!                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
!                        then the next QR sweep is skipped and early
!                        deflation is applied immediately to the
!                        remaining active diagonal block.  Setting
!                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!                        multi-shift QR sweep whenever early deflation
!                        finds a converged eigenvalue.  Setting
!                        IPARMQ(ISPEC=14) greater than or equal to 100
!                        prevents TTQRE from skipping a multi-shift
!                        QR sweep.
!
!              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!                        a multi-shift QR iteration.
!
!              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!                        following meanings.
!                        0:  During the multi-shift QR sweep,
!                            xLAQR5 does not accumulate reflections and
!                            does not use matrix-matrix multiply to
!                            update the far-from-diagonal matrix
!                            entries.
!                        1:  During the multi-shift QR sweep,
!                            xLAQR5 and/or xLAQRaccumulates reflections and uses
!                            matrix-matrix multiply to update the
!                            far-from-diagonal matrix entries.
!                        2:  During the multi-shift QR sweep.
!                            xLAQR5 accumulates reflections and takes
!                            advantage of 2-by-2 block structure during
!                            matrix-matrix multiplies.
!                        (If xTRMM is slower than xGEMM, then
!                        IPARMQ(ISPEC=16)=1 may be more efficient than
!                        IPARMQ(ISPEC=16)=2 despite the greater level of
!                        arithmetic work implied by the latter choice.)
!
!       NAME    (input) character string
!               Name of the calling subroutine
!
!       OPTS    (input) character string
!               This is a concatenation of the string arguments to
!               TTQRE.
!
!       N       (input) integer scalar
!               N is the order of the Hessenberg matrix H.
!
!       ILO     (input) INTEGER
!       IHI     (input) INTEGER
!               It is assumed that H is already upper triangular
!               in rows and columns 1:ILO-1 and IHI+1:N.
!
!       LWORK   (input) integer scalar
!               The amount of workspace available.
!
!  Further Details
!  ===============
!
!       Little is known about how best to choose these parameters.
!       It is possible to use different values of the parameters
!       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!
!       It is probably best to choose different parameters for
!       different matrices and different parameters at different
!       times during the iteration, but this has not been
!       implemented --- yet.
!
!
!       The best choices of most of the parameters depend
!       in an ill-understood way on the relative execution
!       rate of xLAQR3 and xLAQR5 and on the nature of each
!       particular eigenvalue problem.  Experiment may be the
!       only practical way to determine which choices are most
!       effective.
!
!       Following is a list of default values supplied by IPARMQ.
!       These defaults may be adjusted in order to attain better
!       performance in any particular computational environment.
!
!       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!                        Default: 75. (Must be at least 11.)
!
!       IPARMQ(ISPEC=13) Recommended deflation window size.
!                        This depends on ILO, IHI and NS, the
!                        number of simultaneous shifts returned
!                        by IPARMQ(ISPEC=15).  The default for
!                        (IHI-ILO+1).LE.500 is NS.  The default
!                        for (IHI-ILO+1).GT.500 is 3*NS/2.
!
!       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!
!       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!                        a multi-shift QR iteration.
!
!                        If IHI-ILO+1 is ...
!
!                        greater than      ...but less    ... the
!                        or equal to ...      than        default is
!
!                                0               30       NS =   2+
!                               30               60       NS =   4+
!                               60              150       NS =  10
!                              150              590       NS =  **
!                              590             3000       NS =  64
!                             3000             6000       NS = 128
!                             6000             infinity   NS = 256
!
!                    (+)  By default matrices of this order are
!                         passed to the implicit double shift routine
!                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!                         values of NS are used only in case of a rare
!                         xLAHQR failure.
!
!                    (**) The asterisks (**) indicate an ad-hoc
!                         function increasing from 10 to 64.
!
!       IPARMQ(ISPEC=16) Select structured matrix multiply.
!                        (See ISPEC=16 above for details.)
!                        Default: 3.
!
!     ================================================================
!     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, &
                         ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, &
                         NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
!     ..
!     .. Local Scalars ..
      INTEGER            NH, NS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
!     ..
!     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. &
          ( ISPEC.EQ.IACC22 ) ) THEN
!
!        ==== Set the number simultaneous shifts ====
!
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 ) &
            NS = 4
         IF( NH.GE.60 ) &
            NS = 10
         IF( NH.GE.150 ) &
            NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 ) &
            NS = 64
         IF( NH.GE.3000 ) &
            NS = 128
         IF( NH.GE.6000 ) &
            NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
!
      IF( ISPEC.EQ.INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         IPARMQ = NMIN
!
      ELSE IF( ISPEC.EQ.INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         IPARMQ = NIBBLE
!
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         IPARMQ = NS
!
      ELSE IF( ISPEC.EQ.INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
!
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
         IPARMQ = 0
         IF( NS.GE.KACMIN ) &
            IPARMQ = 1
         IF( NS.GE.K22MIN ) &
            IPARMQ = 2
!
      ELSE
!        ===== invalid value of ispec =====
         IPARMQ = -1
!
      END IF
!
!     ==== End of IPARMQ ====
!
      END
