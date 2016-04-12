C     $Id: vctcnv.f 3304 2011-01-17 15:25:59Z brideout $
C
      SUBROUTINE VCTCNV(FX,FY,FZ,X,Y,Z,FR,FT,FP,R,THETA,PHI,IMODE)
C
C     jmh - 11/79  ans fortran 66
C
C     VCTCNV converts between the cartesian and spherical coordinate
C     representations of a vector field f. (FX,FY,FZ) are the
C     components of the field at (X,Y,Z). (FR,FT,FP) are the
C     components of the field at (R,THETA,PHI) in the directions of
C     increasing R, increasing THETA and increasing PHI. if IMODE=1,
C     (FX,FY,FZ,X,Y,Z) -> (FR,FT,FP,R,THETA,PHI). if IMODE=2,
C     (FR,FT,FP,R,THETA,PHI) -> (FX,FY,FZ,X,Y,Z). THETA and PHI are
C     in degrees.
C
C       Input:
C         IMODE - 1 cartesian to spherical
C                 2 spherical to cartesian
C
C     Input, output:
C      FX,FY,FZ - cartesian vector field components
C         X,Y,Z - cartesian vector field coordinates
C      FR,FT,FP - spherical vector field components
C   R,THETA,PHI - spherical vector field coordinates
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION FP,FR,FT,FX,FY,FZ,PHI,R,THETA,X,Y,Z
      INTEGER IMODE
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CP,CT,DTR,P,RHO2,RLMIN,SP,ST,T
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DATAN2,DCOS,DMAX1,DSIN,DSQRT,SIGN
C     ..
C     .. Data statements ..
      DATA DTR/.0174532925199D0/,RLMIN/1.D-38/
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
         FX = ST*CP*FR + CT*CP*FT - SP*FP
         FY = ST*SP*FR + CT*SP*FT + CP*FP
         FZ = CT*FR - ST*FT
      ELSE
C
         RHO2 = X*X + Y*Y
         R = DSQRT(RHO2+Z*Z)
         THETA = DATAN2(DSQRT(RHO2),SIGN(1.D0,Z)*DMAX1(ABS(Z),RLMIN))
         PHI = DATAN2(Y,SIGN(1.D0,X)*DMAX1(ABS(X),RLMIN))
         CT = DCOS(THETA)
         ST = DSIN(THETA)
         CP = DCOS(PHI)
         SP = DSIN(PHI)
         THETA = THETA/DTR
         PHI = PHI/DTR
         FR = ST*CP*FX + ST*SP*FY + CT*FZ
         FT = CT*CP*FX + CT*SP*FY - ST*FZ
         FP = -SP*FX + CP*FY
      END IF
      RETURN
C
      END
