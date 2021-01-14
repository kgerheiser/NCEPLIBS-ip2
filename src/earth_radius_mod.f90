module earth_radius_mod
  implicit none

  private
  public :: earth_radius

contains

  SUBROUTINE EARTH_RADIUS(IGDTMPL, IGDTLEN, RADIUS, ECCEN_SQUARED)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  EARTH_RADIUS   DETERMINE EARTH RADIUS AND SHAPE
    !   PRGMMR: GAYNO    ORG: W/NMC23     DATE: 2015-07-14
    !
    ! ABSTRACT: DETERMINE THE RADIUS AND SHAPE OF THE EARTH FROM
    !   THE GRIB 2 GRID DEFINITION TEMPLATE ARRAY - SECTION 3
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-07-14  GAYNO
    !
    ! USAGE:   CALL EARTH_RADIUS(IGDTMPL, IGDTLEN, &
    !                            RADIUS, ECCEN_SQUARED)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTMPL  - INTEGER (IGDTLEN) GRID DEFINITION TEMPLATE ARRAY.
    !                CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                FOR ALL MAP PROJECTIONS RECOGNIZED BY IPLIB,
    !                THE ENTRIES USE BY THIS ROUTINE ARE:
    !                 (1) - SHAPE OF EARTH, SECTION 3, OCTET 15
    !                 (2) - SCALE FACTOR OF SPHERICAL EARTH RADIUS,
    !                       OCTET 16
    !                 (3) - SCALED VALUE OF RADIUS OF SPHERICAL EARTH,
    !                       OCTETS 17-20
    !                 (4) - SCALE FACTOR OF MAJOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTET 21
    !                 (5) - SCALED VALUE OF MAJOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTETS 22-25
    !                 (6) - SCALE FACTOR OF MINOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTET 26
    !                 (7) - SCALED VALUE OF MINOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTETS 27-30
    !     IGDTLEN  - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD
    !                DATA STRUCTURE.
    !
    !   OUTPUT ARGUMENT LIST:
    !     RADIUS        - REAL EARTH RADIUS IN METERS.
    !                     FOR ELLIPITICAL EARTHS, THIS IS THE
    !                     SEMI MAJOR AXIS.  SEE "MAP PROJECTSIONS -
    !                     A WORKING MANUAL" BY SNYDER (1987)
    !                     FOR DETAILS.
    !     ECCEN_SQUARED - REAL EARTH ECCENTRICITY SQUARED
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    INTEGER,                INTENT(IN   ) :: IGDTLEN
    INTEGER,                INTENT(IN   ) :: IGDTMPL(IGDTLEN)
    !
    REAL,                   INTENT(  OUT) :: ECCEN_SQUARED
    REAL,                   INTENT(  OUT) :: RADIUS
    !
    REAL                                  :: FLAT
    REAL                                  :: MAJOR_AXIS, MINOR_AXIS
    !
    SELECT CASE (IGDTMPL(1))
    CASE (0)
       RADIUS        = 6367470.0
       ECCEN_SQUARED = 0.0
    CASE (1)  ! USER SPECIFIED SPHERICAL
       RADIUS        = FLOAT(IGDTMPL(3))/FLOAT(10**IGDTMPL(2))
       ECCEN_SQUARED = 0.0
    CASE (2)  ! IAU 1965
       RADIUS        = 6378160.0      ! SEMI MAJOR AXIS
       FLAT          = 1.0/297.0      ! FLATTENING
       ECCEN_SQUARED = (2.0*FLAT) - (FLAT**2)
    CASE (3)  ! USER SPECIFIED ELLIPTICAL (KM)
       MAJOR_AXIS    = FLOAT(IGDTMPL(5))/FLOAT(10**IGDTMPL(4))
       MAJOR_AXIS    = MAJOR_AXIS * 1000.0
       MINOR_AXIS    = FLOAT(IGDTMPL(7))/FLOAT(10**IGDTMPL(6))
       MINOR_AXIS    = MINOR_AXIS * 1000.0
       ECCEN_SQUARED = 1.0 - (MINOR_AXIS**2 / MAJOR_AXIS**2)
       RADIUS        = MAJOR_AXIS
    CASE (4)  ! IAG-GRS80 MODEL
       RADIUS        = 6378137.0      ! SEMI MAJOR AXIS
       FLAT          = 1.0/298.2572   ! FLATTENING
       ECCEN_SQUARED = (2.0*FLAT) - (FLAT**2)
    CASE (5)  ! WGS84 DATUM
       RADIUS        = 6378137.0      ! SEMI MAJOR AXIS
       ECCEN_SQUARED = 0.00669437999013
    CASE (6)
       RADIUS        = 6371229.0
       ECCEN_SQUARED = 0.0
    CASE (7)  ! USER SPECIFIED ELLIPTICAL (M)
       MAJOR_AXIS    = FLOAT(IGDTMPL(5))/FLOAT(10**IGDTMPL(4))
       MINOR_AXIS    = FLOAT(IGDTMPL(7))/FLOAT(10**IGDTMPL(6))
       ECCEN_SQUARED = 1.0 - (MINOR_AXIS**2 / MAJOR_AXIS**2)
       RADIUS        = MAJOR_AXIS
    CASE (8)
       RADIUS        = 6371200.0
       ECCEN_SQUARED = 0.0
    CASE DEFAULT
       RADIUS        = -9999.
       ECCEN_SQUARED = -9999.
    END SELECT
    !
    RETURN
    !
  END SUBROUTINE EARTH_RADIUS
end module earth_radius_mod
