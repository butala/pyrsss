C     $Id: look.f 3304 2011-01-17 15:25:59Z brideout $
C
      SUBROUTINE LOOK(SR,SLAT,SLON,PR,GLAT,GLON,AZ,EL,RANGE)
C
C     jmh - 1/80  ans fortran 66
C
C     LOOK calculates the azimuth, elevation and range from a radar
C     of a specified point.
C
C       Input:
C        SR    - distance of station from center of earth (km)
C        SLAT  - geocentric latitude of station (deg)
C        SLON  - longitude of station (deg)
C        PR    - distance from center of earth of observation point (km)
C        GLAT  - observation point geocentric latitude (deg)
C        GLON  - observation point longitude (deg)
C
C      Output:
C        AZ    - radar azimuth (deg)
C        EL    - radar elevation (deg)
C        RANGE - radar range (km)
C
C     ...calculate "observation-point" earth centered cartesian coords..
C     .. Scalar Arguments ..
      DOUBLE PRECISION AZ,EL,GLAT,GLON,PR,RANGE,SLAT,SLON,SR
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RP,RR,RT,SP1,SR1,ST2
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION P(3),R(3),S(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL CSCONV,VCTCNV,VSUB
C     ..
      CALL CSCONV(P(1),P(2),P(3),PR,90.0D0-GLAT,GLON,2)
C
C     ...calculate "station" earth centered cartesian coordinates...
      CALL CSCONV(S(1),S(2),S(3),SR,90.0D0-SLAT,SLON,2)
C
C     ...calculate "line-of-sight" earth centered cartesian coords...
      CALL VSUB(P,S,R)
C
C     ...calculate "line-of-sight" station centered cartesian coords...
      CALL VCTCNV(R(1),R(2),R(3),S(1),S(2),S(3),RR,RT,RP,SR1,ST2,SP1,1)
C
C     ...calculate "line-of-sight" station centered spherical coords...
      CALL CSCONV(RT,RP,RR,RANGE,EL,AZ,1)
      EL = 90.D0 - EL
      AZ = 180.D0 - AZ
      RETURN
C
      END
