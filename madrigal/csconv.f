C     $Id: csconv.f 3304 2011-01-17 15:25:59Z brideout $
C
      SUBROUTINE CSCONV(X,Y,Z,R,THETA,PHI,IMODE)
C
C     jmh - 11/79  ans fortran 66
C
C     Converts between cartesian coordinates x,y,z and spherical
C     coordinates r,theta,phi.  theta and phi are in degrees.
C
C       Input:
C         IMODE - 1 (x,y,z) -> (r,theta,phi)
C                 2 (r,theta,phi) -> (x,y,z)
C
C       Input, Output:
C              X, Y, Z - cartesian coordinates
C        R, THETA, PHI - spherical coordinates (degrees)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION PHI,R,THETA,X,Y,Z
      INTEGER IMODE
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CP,CT,DTR,P,RHO2,RLMIN,SP,ST,T
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DATAN2,DCOS,DMAX1,DSIN,DSQRT,SIGN
C     ..
C     .. Data statements ..
      DATA DTR/.0174532925199D0/
      DATA RLMIN/1.D-38/
C     ..
C
      IF (IMODE.EQ.2) THEN
C
         T = DTR*THETA
         P = DTR*PHI
         CT = DCOS(T)
         ST = DSIN(T)
         CP = DCOS(P)
         SP = DSIN(P)
         X = R*ST*CP
         Y = R*ST*SP
         Z = R*CT
      ELSE
C
         RHO2 = X*X + Y*Y
         R = DSQRT(RHO2+Z*Z)
         THETA = DATAN2(DSQRT(RHO2),SIGN(1.D0,Z)*DMAX1(ABS(Z),RLMIN))/
     *           DTR
         PHI = DATAN2(Y,SIGN(1.D0,X)*DMAX1(ABS(X),RLMIN))/DTR
      END IF
      RETURN
C
      END
