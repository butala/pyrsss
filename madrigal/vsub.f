C     $Id: vsub.f 3304 2011-01-17 15:25:59Z brideout $
C
      SUBROUTINE VSUB(A,B,C)
C
C     jmh - 11/79  ans fortran 66
C
C     VSUB calculates the difference of two vectors A and B, C = A - B.
C
C     Input:
C        A - floating point vector of dimension 3
C        B - floating point vector of dimension 3
C
C     Output:
C        C - floating point vector of dimension 3
C
C     .. Array Arguments ..
      DOUBLE PRECISION A(3),B(3),C(3)
C     ..
      C(1) = A(1) - B(1)
      C(2) = A(2) - B(2)
      C(3) = A(3) - B(3)
      RETURN
C
      END
