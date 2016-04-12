C     $Id: point.f 3304 2011-01-17 15:25:59Z brideout $
C
      SUBROUTINE POINT(SR,SLAT,SLON,AZ,EL,RANGE,PR,GLAT,GLON)
C
C     jmh - 1/80  ans fortran 66
C
C     POINT calculates the position of a point defined by the radar
C     line-of sight vector to that point.
C
C     Input:
C       SR    - distance of station from center of earth (km)
C       SLAT  - geocentric latitude of station (deg)
C       SLON  - longitude of station (deg)
C       AZ    - radar azimuth (deg)
C       EL    - radar elevation (deg)
C       RANGE - radar range (km)
C
C     Output:
C       PR    - distance from center of earth of observation point (km)
C       GLAT  - observation point geocentric latitude (deg)
C       GLON  - observation point longitude (deg)
C
C     ...calculate "line-of-sight" station centered cartesian coords...
C     .. Scalar Arguments ..
      DOUBLE PRECISION AZ,EL,GLAT,GLON,PR,RANGE,SLAT,SLON,SR
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RP,RR,RT,T
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION P(3),R(3),S(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL CSCONV,VADD,VCTCNV
C     ..
      CALL CSCONV(RT,RP,RR,RANGE,90.0D0-EL,180.0D0-AZ,2)
C
C     ...calculate "line-of-sight" earth centered cartesian coords
C          and "station" earth centered cartesian coords...
      CALL VCTCNV(R(1),R(2),R(3),S(1),S(2),S(3),RR,RT,RP,SR,90.0D0-SLAT,
     *            SLON,2)
C
C     ...calculate "observation-point" earth centered cartesian coords..
      CALL VADD(S,R,P)
C
C     ...calculate "observation-point" earth centered spherical coords..
      CALL CSCONV(P(1),P(2),P(3),PR,T,GLON,1)
      GLAT = 90.0D0 - T
      RETURN
C
      END
