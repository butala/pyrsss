C     $Id: convrt.f 3304 2011-01-17 15:25:59Z brideout $
C
      SUBROUTINE CONVRT(I,GDLAT,GDALT,GCLAT,RKM)
C
C     jmh - 11/79  ans fortran 66
C
C     Converts between geodetic and geocentric coordinates. the
C     reference geoid is that adopted by the iau in 1964. a=6378.16,
C     b=6356.7746, f=1/298.25. the equations for conversion from
C     geocentric to geodetic are from astron. j., vol 66, 1961, p. 15.
C
C       Input:
C             I - 1 geodetic to geocentric, 2 geocentric to geodetic
C       Input, Output:
C         GDLAT - geodetic latitude (degrees)
C         GDALT - altitude above geoid (km)
C         GCLAT - geocentric latitude (degrees)
C           RKM - geocentric radial distance (km)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION GCLAT,GDALT,GDLAT,RKM
      INTEGER I
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,A2,A4,A6,A8,AB2,C2CL,C4CL,CCL,CL2,COSBET,
     *                 COSLAT,DLTCL,DTR,EP2,GCL,GDL,RER,RGEOID,S2CL,
     *                 S4CL,S6CL,S8CL,SB2,SCL,SINBET,SINLAT,SL2,X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN2,DCOS,DMIN1,DSIN,DSQRT
C     ..
C     .. Data statements ..
C
      DATA A/6378.16D0/,AB2/1.0067397D0/,EP2/.0067397D0/
      DATA DTR/.0174532925199D0/
C     ..
C
      IF (I.EQ.2) THEN
C
C        .....geocentic to geodetic.....
         RER = RKM/A
         A2 = ((-1.4127348D-8/RER+.94339131D-8)/RER+.33523288D-2)/RER
         A4 = (((-1.2545063D-10/RER+.11760996D-9)/RER+.11238084D-4)/RER-
     *        .2814244D-5)/RER
         A6 = ((54.939685D-9/RER-28.301730D-9)/RER+3.5435979D-9)/RER
         A8 = (((320.D0/RER-252.D0)/RER+64.D0)/RER-5.D0)/RER*
     *        .98008304D-12
         GCL = DTR*GCLAT
         CCL = DCOS(GCL)
         SCL = DSIN(GCL)
         S2CL = 2.D0*SCL*CCL
         C2CL = 2.D0*CCL*CCL - 1.D0
         S4CL = 2.D0*S2CL*C2CL
         C4CL = 2.D0*C2CL*C2CL - 1.D0
         S8CL = 2.D0*S4CL*C4CL
         S6CL = S2CL*C4CL + C2CL*S4CL
         DLTCL = S2CL*A2 + S4CL*A4 + S6CL*A6 + S8CL*A8
         GDLAT = GCLAT + DLTCL/DTR
         GDALT = RKM - A/DSQRT(1.D0+EP2*SCL*SCL)
      ELSE
C
C        .....geodetic to geocentric.....
         GDL = DTR*GDLAT
         SINLAT = DSIN(GDL)
         COSLAT = DCOS(GDL)
         SL2 = SINLAT*SINLAT
         CL2 = AB2*COSLAT
         CL2 = CL2*CL2
         SINBET = SINLAT/DSQRT(SL2+CL2)
         SB2 = DMIN1(SINBET*SINBET,1.D0)
         COSBET = DSQRT(1.D0-SB2)
         RGEOID = A/DSQRT(1.D0+EP2*SB2)
         X = RGEOID*COSBET + GDALT*COSLAT
         Y = RGEOID*SINBET + GDALT*SINLAT
         RKM = DSQRT(X*X+Y*Y)
         GCLAT = DATAN2(Y,X)/DTR
      END IF
      RETURN
C
      END
