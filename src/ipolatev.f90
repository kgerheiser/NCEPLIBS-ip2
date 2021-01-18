module ipolatev_mod
  use polatev0_mod
  use polatev1_mod
  use polatev2_mod
  use polatev3_mod
  use polatev4_mod
  use polatev6_mod

  use ip_grid_factory_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod

  implicit none

  private
  public :: ipolatev

  interface ipolatev
     module procedure ipolatev_grib1
     module procedure ipolatev_grib2
  end interface ipolatev

contains

  SUBROUTINE IPOLATEV_grib2(IP,IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       IGDTNUMO,IGDTMPLO,IGDTLENO, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  IPOLATEV   IREDELL'S POLATE FOR VECTOR FIELDS
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM INTERPOLATES VECTOR FIELDS
    !           FROM ANY GRID TO ANY GRID (JOE IRWIN'S DREAM).
    !           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
    !           THE FOLLOWING INTERPOLATION METHODS ARE POSSIBLE:
    !             (IP=0) BILINEAR
    !             (IP=1) BICUBIC
    !             (IP=2) NEIGHBOR
    !             (IP=3) BUDGET
    !             (IP=4) SPECTRAL
    !             (IP=6) NEIGHBOR-BUDGET
    !           SOME OF THESE METHODS HAVE INTERPOLATION OPTIONS AND/OR
    !           RESTRICTIONS ON THE INPUT OR OUTPUT GRIDS, BOTH OF WHICH
    !           ARE DOCUMENTED MORE FULLY IN THEIR RESPECTIVE SUBPROGRAMS.
    !
    !           THE INPUT AND OUTPUT GRIDS ARE DEFINED BY THEIR GRIB 2 GRID
    !           DEFINITION TEMPLATE AS DECODED BY THE NCEP G2 LIBRARY.  THE
    !           CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
    !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
    !           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
    !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
    !             (IGDTNUMI/O=01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
    !                             NON-"E" STAGGERED
    !             (IGDTNUMI/O=10) MERCATOR CYLINDRICAL
    !             (IGDTNUMI/O=20) POLAR STEREOGRAPHIC AZIMUTHAL
    !             (IGDTNUMI/O=30) LAMBERT CONFORMAL CONICAL
    !             (IGDTNUMI/O=40) GAUSSIAN CYLINDRICAL
    !
    !           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
    !           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
    !           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
    !           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
    !           AS DESIGNATED BY THEIR RESPECTIVE GRID DEFINITION SECTIONS.
    !
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
    !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
    !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !           ON THE OTHER HAND, THE DATA MAY BE INTERPOLATED TO A SET OF
    !           STATION POINTS IF IGDTNUMO<0 (IGDTNUMO-255 FOR THE BUDGET OPTION),
    !           IN WHICH CASE THE NUMBER OF POINTS AND THEIR LATITUDES AND
    !           LONGITUDES MUST BE INPUT ALONG WITH THEIR VECTOR ROTATION 
    !           PARAMETERS.
    !
    !           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
    !           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
    !           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
    !           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    ! 2003-06-23  IREDELL  STAGGERING FOR GRID TYPE 203
    ! 2015-01-27  GAYNO    REMOVE REFERENCES TO OBSOLETE NCEP GRIDS 201
    !                      AND 202.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL IPOLATEV(IP,IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,UI,VI, &
    !                         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IP       - INTEGER INTERPOLATION METHOD
    !                (IP=0 FOR BILINEAR;
    !                 IP=1 FOR BICUBIC;
    !                 IP=2 FOR NEIGHBOR;
    !                 IP=3 FOR BUDGET;
    !                 IP=4 FOR SPECTRAL;
    !                 IP=6 FOR NEIGHBOR-BUDGET)
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                (IP=0: MIN MASK %
    !                 IP=1: CONSTRAINT OPTION, MIN MASK %
    !                 IP=2: SEARCH RADIUS
    !                 IP=3: NUMBER IN RADIUS, RADIUS WEIGHTS, MIN MASK %
    !                 IP=4: SPECTRAL SHAPE, SPECTRAL TRUNCATION
    !                 IP=6: NUMBER IN RADIUS, RADIUS WEIGHTS, MIN MASK %
    !     IGDTNUMI - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE:
    !                  00 - EQUIDISTANT CYLINDRICAL
    !                  01 - ROTATED EQUIDISTANT CYLINDRICAL.  "E"
    !                       AND NON-"E" STAGGERED
    !                  10 - MERCATOR CYCLINDRICAL
    !                  20 - POLAR STEREOGRAPHIC AZIMUTHAL
    !                  30 - LAMBERT CONFORMAL CONICAL
    !                  40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
    !     IGDTMPLI - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY -
    !                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE
    !                (SECTION 3 INFO):
    !                ALL MAP PROJECTIONS:
    !                 (1):  SHAPE OF EARTH, OCTET 15
    !                 (2):  SCALE FACTOR OF SPHERICAL EARTH RADIUS,
    !                       OCTET 16
    !                 (3):  SCALED VALUE OF RADIUS OF SPHERICAL EARTH,
    !                       OCTETS 17-20
    !                 (4):  SCALE FACTOR OF MAJOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTET 21
    !                 (5):  SCALED VALUE OF MAJOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTETS 22-25
    !                 (6):  SCALE FACTOR OF MINOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTET 26
    !                 (7):  SCALED VALUE OF MINOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTETS 27-30
    !                EQUIDISTANT CYCLINDRICAL:
    !                 (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
    !                 (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
    !                 (10): BASIC ANGLE OF INITIAL PRODUCTION DOMAIN,
    !                       OCTETS 39-42.
    !                 (11): SUBDIVISIONS OF BASIC ANGLE, OCTETS 43-46
    !                 (12): LATITUDE OF FIRST GRID POINT, OCTETS 47-50
    !                 (13): LONGITUDE OF FIRST GRID POINT, OCTETS 51-54
    !                 (14): RESOLUTION AND COMPONENT FLAGS, OCTET 55
    !                 (15): LATITUDE OF LAST GRID POINT, OCTETS 56-59
    !                 (16): LONGITUDE OF LAST GRID POINT, OCTETS 60-63
    !                 (17): I-DIRECTION INCREMENT, OCTETS 64-67
    !                 (18): J-DIRECTION INCREMENT, OCTETS 68-71
    !                 (19): SCANNING MODE, OCTET 72
    !                MERCATOR CYCLINDRICAL:
    !                 (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
    !                 (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
    !                 (10): LATITUDE OF FIRST POINT, OCTETS 39-42
    !                 (11): LONGITUDE OF FIRST POINT, OCTETS 43-46
    !                 (12): RESOLUTION AND COMPONENT FLAGS, OCTET 47
    !                 (13): TANGENT LATITUDE, OCTETS 48-51
    !                 (14): LATITUDE OF LAST POINT, OCTETS 52-55
    !                 (15): LONGITUDE OF LAST POINT, OCTETS 56-59
    !                 (16): SCANNING MODE FLAGS, OCTET 60
    !                 (17): ORIENTATION OF GRID, OCTETS 61-64
    !                 (18): LONGITUDINAL GRID LENGTH, OCTETS 65-68
    !                 (19): LATITUDINAL GRID LENGTH, OCTETS 69-72
    !                LAMBERT CONFORMAL CONICAL:
    !                 (8):  NUMBER OF POINTS ALONG X-AXIS, OCTS 31-34
    !                 (9):  NUMBER OF POINTS ALONG Y-AXIS, OCTS 35-38
    !                 (10): LATITUDE OF FIRST POINT, OCTETS 39-42
    !                 (11): LONGITUDE OF FIRST POINT, OCTETS 43-46
    !                 (12): RESOLUTION OF COMPONENT FLAG, OCTET 47
    !                 (13): LATITUDE WHERE GRID LENGTHS SPECIFIED,
    !                       OCTETS 48-51
    !                 (14): LONGITUDE OF MERIDIAN THAT IS PARALLEL TO
    !                       Y-AXIS, OCTETS 52-55
    !                 (15): X-DIRECTION GRID LENGTH, OCTETS 56-59
    !                 (16): Y-DIRECTION GRID LENGTH, OCTETS 60-63
    !                 (17): PROJECTION CENTER FLAG, OCTET 64
    !                 (18): SCANNING MODE, OCTET 65
    !                 (19): FIRST TANGENT LATITUDE FROM POLE, OCTETS 66-69
    !                 (20): SECOND TANGENT LATITUDE FROM POLE, OCTETS 70-73
    !                 (21): LATITUDE OF SOUTH POLE OF PROJECTION,
    !                       OCTETS 74-77
    !                 (22): LONGITUDE OF SOUTH POLE OF PROJECTION,
    !                       OCTETS 78-81
    !                GAUSSIAN CYLINDRICAL:
    !                 (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
    !                 (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
    !                 (10): BASIC ANGLE OF INITIAL PRODUCTION DOMAIN,
    !                       OCTETS 39-42
    !                 (11): SUBDIVISIONS OF BASIC ANGLE, OCTETS 43-46
    !                 (12): LATITUDE OF FIRST GRID POINT, OCTETS 47-50
    !                 (13): LONGITUDE OF FIRST GRID POINT, OCTETS 51-54
    !                 (14): RESOLUTION AND COMPONENT FLAGS, OCTET 55
    !                 (15): LATITUDE OF LAST GRID POINT, OCTETS 56-59
    !                 (16): LONGITUDE OF LAST GRID POINT, OCTETS 60-63
    !                 (17): I-DIRECTION INCREMENT, OCTETS 64-67
    !                 (18): NUMBER OF PARALLELS BETWEEN POLE AND EQUATOR,
    !                       OCTETS 68-71
    !                 (19): SCANNING MODE, OCTET 72
    !                POLAR STEREOGRAPHIC AZIMUTHAL:
    !                 (8):  NUMBER OF POINTS ALONG X-AXIS, OCTETS 31-34
    !                 (9):  NUMBER OF POINTS ALONG Y-AXIS, OCTETS 35-38
    !                 (10): LATITUDE OF FIRST GRID POINT, OCTETS 39-42
    !                 (11): LONGITUDE OF FIRST GRID POINT, OCTETS 43-46
    !                 (12): RESOLUTION AND COMPONENT FLAGS, OCTET 47
    !                 (13): TRUE LATITUDE, OCTETS 48-51
    !                 (14): ORIENTATION LONGITUDE, OCTETS 52-55
    !                 (15): X-DIRECTION GRID LENGTH, OCTETS 56-59
    !                 (16): Y-DIRECTION GRID LENGTH, OCTETS 60-63
    !                 (17): PROJECTION CENTER FLAG, OCTET 64
    !                 (18): SCANNING MODE FLAGS, OCTET 65
    !                ROTATED EQUIDISTANT CYCLINDRICAL:
    !                 (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
    !                 (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
    !                 (10): BASIC ANGLE OF INITIAL PRODUCTION DOMAIN,
    !                       OCTETS 39-42
    !                 (11): SUBDIVISIONS OF BASIC ANGLE, OCTETS 43-46
    !                 (12): LATITUDE OF FIRST GRID POINT, OCTETS 47-50
    !                 (13): LONGITUDE OF FIRST GRID POINT, OCTETS 51-54
    !                 (14): RESOLUTION AND COMPONENT FLAGS, OCTET 55
    !                 (15): LATITUDE OF LAST GRID POINT, OCTETS 56-59
    !                 (16): LONGITUDE OF LAST GRID POINT, OCTETS 60-63
    !                 (17): I-DIRECTION INCREMENT, OCTETS 64-67
    !                 (18): J-DIRECTION INCREMENT, OCTETS 68-71
    !                 (19): SCANNING MODE, OCTET 72
    !                 (20): LATITUDE OF SOUTHERN POLE OF PROJECTION,
    !                       OCTETS 73-76
    !                 (21): LONGITUDE OF SOUTHERN POLE OF PROJECTION,
    !                       OCTETS 77-80
    !                 (22): ANGLE OF ROTATION OF PROJECTION, OCTS 81-84
    !     IGDTLENI - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IGDTNUMO - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE. SEE "IGDTNUMI"
    !                FOR SPECIFIC TEMPLATE DEFINITIONS. NOTE: IGDTNUMO<0
    !                MEANS INTERPOLATE TO RANDOM STATION POINTS.
    !     IGDTMPLO - INTEGER (IGDTLENO) GRID DEFINITION TEMPLATE ARRAY -
    !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                SEE "IGDTMPLI" FOR DEFINITION OF ARRAY ELEMENTS.
    !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
    !     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
    !     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
    !     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF RESPECTIVE IBI(K)=1)
    !     UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
    !     VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO<0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO<0)
    !     CROT     - REAL (MO) VECTOR ROTATION COSINES (IF IGDTNUMO<0)
    !     SROT     - REAL (MO) VECTOR ROTATION SINES (IF IGDTNUMO<0)
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !
    !   OUTPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO>=0)
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     CROT     - REAL (MO) VECTOR ROTATION COSINES (IF IGDTNUMO>=0)
    !     SROT     - REAL (MO) VECTOR ROTATION SINES (IF IGDTNUMO>=0)
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
    !     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
    !     UO       - REAL (MO,KM) OUTPUT U-COMPONENT FIELDS INTERPOLATED
    !     VO       - REAL (MO,KM) OUTPUT V-COMPONENT FIELDS INTERPOLATED
    !     IRET     - INTEGER RETURN CODE
    !                0    SUCCESSFUL INTERPOLATION
    !                1    UNRECOGNIZED INTERPOLATION METHOD
    !                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
    !                3    UNRECOGNIZED OUTPUT GRID
    !                1X   INVALID BICUBIC METHOD PARAMETERS
    !                3X   INVALID BUDGET METHOD PARAMETERS
    !                4X   INVALID SPECTRAL METHOD PARAMETERS
    !
    ! SUBPROGRAMS CALLED:
    !   POLATEV0     INTERPOLATE VECTOR FIELDS (BILINEAR)
    !   POLATEV1     INTERPOLATE VECTOR FIELDS (BICUBIC)
    !   POLATEV2     INTERPOLATE VECTOR FIELDS (NEIGHBOR)
    !   POLATEV3     INTERPOLATE VECTOR FIELDS (BUDGET)
    !   POLATEV4     INTERPOLATE VECTOR FIELDS (SPECTRAL)
    !   POLATEV6     INTERPOLATE VECTOR FIELDS (NEIGHBOR-BUDGET)
    !
    ! REMARKS: EXAMPLES DEMONSTRATING RELATIVE CPU COSTS.
    !   THIS EXAMPLE IS INTERPOLATING 12 LEVELS OF WINDS
    !   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
    !   TO THE 93 X 68 HAWAIIAN MERCATOR GRID (NCEP GRID 204).
    !   THE EXAMPLE TIMES ARE FOR THE C90.  AS A REFERENCE, THE CP TIME
    !   FOR UNPACKING THE GLOBAL 12 PAIRS OF WIND FIELDS IS 0.07 SECONDS.
    !
    !   BILINEAR    0                   0.05
    !   BICUBIC     1   0               0.16
    !   BICUBIC     1   1               0.17
    !   NEIGHBOR    2                   0.02
    !   BUDGET      3   -1,-1           0.94
    !   SPECTRAL    4   0,40            0.31
    !   SPECTRAL    4   1,40            0.33
    !   SPECTRAL    4   0,-1            0.59
    !   N-BUDGET    6   0,-1            0.31
    !
    !   THE SPECTRAL INTERPOLATION IS FAST FOR THE MERCATOR GRID.
    !   HOWEVER, FOR SOME GRIDS THE SPECTRAL INTERPOLATION IS SLOW.
    !   THE FOLLOWING EXAMPLE IS INTERPOLATING 12 LEVELS OF WINDS
    !   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
    !   TO THE 93 X 65 CONUS LAMBERT CONFORMAL GRID (NCEP GRID 211).
    !
    !   METHOD      IP  IPOPT          CP SECONDS
    !   --------    --  -------------  ----------
    !   BILINEAR    0                   0.05
    !   BICUBIC     1   0               0.15
    !   BICUBIC     1   1               0.16
    !   NEIGHBOR    2                   0.02
    !   BUDGET      3   -1,-1           0.92
    !   SPECTRAL    4   0,40            4.51
    !   SPECTRAL    4   1,40            5.77
    !   SPECTRAL    4   0,-1           12.60
    !   N-BUDGET    6   0,-1            0.33
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    INTEGER,               INTENT(IN   ) :: IP, IPOPT(20), IBI(KM)
    INTEGER,               INTENT(IN   ) :: KM, MI, MO
    INTEGER,               INTENT(IN   ) :: IGDTNUMI, IGDTLENI
    INTEGER,               INTENT(IN   ) :: IGDTMPLI(IGDTLENI)
    INTEGER,               INTENT(IN   ) :: IGDTNUMO, IGDTLENO
    INTEGER,               INTENT(IN   ) :: IGDTMPLO(IGDTLENO)
    INTEGER,               INTENT(  OUT) :: IBO(KM), IRET, NO
    !
    LOGICAL*1,             INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,             INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,                  INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,                  INTENT(INOUT) :: CROT(MO),SROT(MO)
    REAL,                  INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,                  INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    !
    INTEGER                              :: K, N

    type(grib2_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(igdtnumi, igdtleni, igdtmpli)
    desc_out = init_descriptor(igdtnumo, igdtleno, igdtmplo)

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  BILINEAR INTERPOLATION
    IF(IP.EQ.0) THEN
       ! CALL POLATEV0(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       !      IGDTNUMO,IGDTMPLO,IGDTLENO, &
       !      MI,MO,KM,IBI,LI,UI,VI,&
       !      NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)

       CALL interpolate_bilinear_vector(IPOPT,grid_in,grid_out, &
            MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  BICUBIC INTERPOLATION
    ELSEIF(IP.EQ.1) THEN
      CALL interpolate_bicubic_vector(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
           NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  NEIGHBOR INTERPOLATION
    ELSEIF(IP.EQ.2) THEN
       ! CALL POLATEV2(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       !      IGDTNUMO,IGDTMPLO,IGDTLENO, &
       !      MI,MO,KM,IBI,LI,UI,VI,&
       !      NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       CALL interpolate_neighbor_vector(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
           NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  BUDGET INTERPOLATION
    ELSEIF(IP.EQ.3) THEN
       CALL interpolate_budget_vector(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  SPECTRAL INTERPOLATION
    ELSEIF(IP.EQ.4) THEN
       CALL POLATEV4(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
            IGDTNUMO,IGDTMPLO,IGDTLENO, &
            MI,MO,KM,IBI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  NEIGHBOR-BUDGET INTERPOLATION
    ELSEIF(IP.EQ.6) THEN
       CALL POLATEV6(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
            IGDTNUMO,IGDTMPLO,IGDTLENO, &
            MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  UNRECOGNIZED INTERPOLATION METHOD
    ELSE
       IF(IGDTNUMO.GE.0) NO=0
       DO K=1,KM
          IBO(K)=1
          DO N=1,NO
             LO(N,K)=.FALSE.
             UO(N,K)=0.
             VO(N,K)=0.
          ENDDO
       ENDDO
       IRET=1
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE IPOLATEV_GRIB2


  !> @file
  !! IREDELL'S POLATE FOR VECTOR FIELDS
  !! @author IREDELL @date 96-04-10
  !
  !> THIS SUBPROGRAM INTERPOLATES VECTOR FIELDS
  !!           FROM ANY GRID TO ANY GRID (JOE IRWIN'S DREAM).
  !!           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
  !!           THE FOLLOWING INTERPOLATION METHODS ARE POSSIBLE:
  !!             (IP=0) BILINEAR
  !!             (IP=1) BICUBIC
  !!             (IP=2) NEIGHBOR
  !!             (IP=3) BUDGET
  !!             (IP=4) SPECTRAL
  !!             (IP=6) NEIGHBOR-BUDGET
  !!           SOME OF THESE METHODS HAVE INTERPOLATION OPTIONS AND/OR
  !!           RESTRICTIONS ON THE INPUT OR OUTPUT GRIDS, BOTH OF WHICH
  !!           ARE DOCUMENTED MORE FULLY IN THEIR RESPECTIVE SUBPROGRAMS.
  !!           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
  !!           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
  !!           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
  !!             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
  !!             (KGDS(1)=001) MERCATOR CYLINDRICAL
  !!             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
  !!             (KGDS(1)=004) GAUSSIAN CYLINDRICAL
  !!             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
  !!             (KGDS(1)=203) ROTATED EQUIDISTANT CYLINDRICAL - E-STAGGER
  !!             (KGDS(1)=205) ROTATED EQUIDISTANT CYLINDRICAL - B-STAGGER
  !!           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
  !!           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
  !!           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
  !!           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
  !!           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
  !!           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
  !!           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
  !!           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
  !!           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  !!           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
  !!           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
  !!           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT 
  !!           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  !!           NOTE: FOR THE BUDGET APPROACH, A SUBSECTION OF THE GRID MAY
  !!           BE OUTPUT BY SUBTRACTING KGDSO(1) FROM 255 AND PASSING
  !!           IN THE LATITUDES AND LONGITUDES OF THE POINTS.
  !!           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
  !!           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
  !!           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
  !!           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
  !!        
  !! PROGRAM HISTORY LOG:
  !!   96-04-10  IREDELL
  !! 2003-06-23  IREDELL  STAGGERING FOR GRID TYPE 203
  !! 2015-01-27  GAYNO    REMOVE REFERENCES TO OBSOLETE NCEP GRIDS 201
  !!                      AND 202.
  !!
  !! USAGE:    CALL IPOLATEV(IP,IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI,
  !!    &                    NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
  !!
  !!   INPUT ARGUMENT LIST:
  !!     IP       - INTEGER INTERPOLATION METHOD
  !!                (IP=0 FOR BILINEAR;
  !!                 IP=1 FOR BICUBIC;
  !!                 IP=2 FOR NEIGHBOR;
  !!                 IP=3 FOR BUDGET;
  !!                 IP=4 FOR SPECTRAL;
  !!                 IP=6 FOR NEIGHBOR-BUDGET)
  !!     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
  !!                (IP=0: (NO OPTIONS)
  !!                 IP=1: CONSTRAINT OPTION
  !!                 IP=2: (NO OPTIONS)
  !!                 IP=3: NUMBER IN RADIUS, RADIUS WEIGHTS ...
  !!                 IP=4: SPECTRAL SHAPE, SPECTRAL TRUNCATION
  !!                 IP=6: NUMBER IN RADIUS, RADIUS WEIGHTS ...)
  !!     KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
  !!                NOTE: IF KGDSI(1)=203, THEN THE 9TH BIT OF  KGDSI(11)
  !!                      IS TEMPORARILY SET TO 1 TO ALERT THE GDS WIZARD
  !!                      THAT THESE FIELDS ARE STAGGERED ETA WINDS.
  !!     KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
  !!                NOTE: IF KGDSO(1)=203, THEN THE 9TH BIT OF KGDSO(11)
  !!                      IS TEMPORARILY SET TO 1 TO ALERT THE GDS WIZARD
  !!                      THAT THESE FIELDS ARE STAGGERED ETA WINDS.
  !!     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
  !!                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
  !!     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
  !!                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
  !!     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
  !!     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
  !!     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF RESPECTIVE IBI(K)=1)
  !!     UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
  !!     VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
  !!     NO       - INTEGER NUMBER OF OUTPUT POINTS (IF KGDSO(1)<0)
  !!     RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
  !!     RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
  !!     CROT     - REAL (NO) VECTOR ROTATION COSINES (IF KGDSO(1)<0)
  !!     SROT     - REAL (NO) VECTOR ROTATION SINES (IF KGDSO(1)<0)
  !!                (UGRID=CROT*UEARTH-SROT*VEARTH;
  !!                 VGRID=SROT*UEARTH+CROT*VEARTH)
  !!
  !!   OUTPUT ARGUMENT LIST:
  !!     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)>=0)
  !!     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)>=0)
  !!     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)>=0)
  !!     CROT     - REAL (MO) VECTOR ROTATION COSINES (IF KGDSO(1)>=0)
  !!     SROT     - REAL (MO) VECTOR ROTATION SINES (IF KGDSO(1)>=0)
  !!                (UGRID=CROT*UEARTH-SROT*VEARTH;
  !!                 VGRID=SROT*UEARTH+CROT*VEARTH)
  !!     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
  !!     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
  !!     UO       - REAL (MO,KM) OUTPUT U-COMPONENT FIELDS INTERPOLATED
  !!     VO       - REAL (MO,KM) OUTPUT V-COMPONENT FIELDS INTERPOLATED
  !!     IRET     - INTEGER RETURN CODE
  !!                0    SUCCESSFUL INTERPOLATION
  !!                1    UNRECOGNIZED INTERPOLATION METHOD
  !!                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
  !!                3    UNRECOGNIZED OUTPUT GRID
  !!                1X   INVALID BICUBIC METHOD PARAMETERS
  !!                3X   INVALID BUDGET METHOD PARAMETERS
  !!                4X   INVALID SPECTRAL METHOD PARAMETERS
  !!
  !! SUBPROGRAMS CALLED:
  !!   POLATEV0     INTERPOLATE VECTOR FIELDS (BILINEAR)
  !!   POLATEV1     INTERPOLATE VECTOR FIELDS (BICUBIC)
  !!   POLATEV2     INTERPOLATE VECTOR FIELDS (NEIGHBOR)
  !!   POLATEV3     INTERPOLATE VECTOR FIELDS (BUDGET)
  !!   POLATEV4     INTERPOLATE VECTOR FIELDS (SPECTRAL)
  !!   POLATEV6     INTERPOLATE VECTOR FIELDS (NEIGHBOR-BUDGET)
  !!
  !! REMARKS: EXAMPLES DEMONSTRATING RELATIVE CPU COSTS.
  !!   THIS EXAMPLE IS INTERPOLATING 12 LEVELS OF WINDS
  !!   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
  !!   TO THE 93 X 68 HAWAIIAN MERCATOR GRID (NCEP GRID 204).
  !!   THE EXAMPLE TIMES ARE FOR THE C90.  AS A REFERENCE, THE CP TIME
  !!   FOR UNPACKING THE GLOBAL 12 PAIRS OF WIND FIELDS IS 0.07 SECONDS.
  !!
  !!   BILINEAR    0                   0.05
  !!   BICUBIC     1   0               0.16
  !!   BICUBIC     1   1               0.17
  !!   NEIGHBOR    2                   0.02
  !!   BUDGET      3   -1,-1           0.94
  !!   SPECTRAL    4   0,40            0.31
  !!   SPECTRAL    4   1,40            0.33
  !!   SPECTRAL    4   0,-1            0.59
  !!   N-BUDGET    6   0,-1            0.31
  !!
  !!   THE SPECTRAL INTERPOLATION IS FAST FOR THE MERCATOR GRID.
  !!   HOWEVER, FOR SOME GRIDS THE SPECTRAL INTERPOLATION IS SLOW.
  !!   THE FOLLOWING EXAMPLE IS INTERPOLATING 12 LEVELS OF WINDS
  !!   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
  !!   TO THE 93 X 65 CONUS LAMBERT CONFORMAL GRID (NCEP GRID 211).
  !!
  !!   METHOD      IP  IPOPT          CP SECONDS
  !!   --------    --  -------------  ----------
  !!   BILINEAR    0                   0.05
  !!   BICUBIC     1   0               0.15
  !!   BICUBIC     1   1               0.16
  !!   NEIGHBOR    2                   0.02
  !!   BUDGET      3   -1,-1           0.92
  !!   SPECTRAL    4   0,40            4.51
  !!   SPECTRAL    4   1,40            5.77
  !!   SPECTRAL    4   0,-1           12.60
  !!   N-BUDGET    6   0,-1            0.33
  !!
  SUBROUTINE IPOLATEV_grib1(IP,IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    IMPLICIT NONE
    !
    INTEGER,               INTENT(IN   ):: IP, IPOPT(20), IBI(KM)
    INTEGER,               INTENT(IN   ):: KM, MI, MO
    INTEGER,               INTENT(INOUT):: KGDSI(200), KGDSO(200)
    INTEGER,               INTENT(  OUT):: IBO(KM), IRET, NO
    !
    LOGICAL*1,             INTENT(IN   ):: LI(MI,KM)
    LOGICAL*1,             INTENT(  OUT):: LO(MO,KM)
    !
    REAL,                  INTENT(IN   ):: UI(MI,KM),VI(MI,KM)
    REAL,                  INTENT(INOUT):: CROT(MO),SROT(MO)
    REAL,                  INTENT(INOUT):: RLAT(MO),RLON(MO)
    REAL,                  INTENT(  OUT):: UO(MO,KM),VO(MO,KM)
    !
    INTEGER                             :: K, N, KGDSI11, KGDSO11

    type(grib1_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(KGDSI(1).EQ.203) THEN
       KGDSI11=KGDSI(11)
       KGDSI(11)=IOR(KGDSI(11),256)
    ENDIF
    IF(KGDSO(1).EQ.203) THEN
       KGDSO11=KGDSO(11)
       KGDSO(11)=IOR(KGDSO(11),256)
    ENDIF

    desc_in = init_descriptor(kgdsi)
    desc_out = init_descriptor(kgdso)

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  BILINEAR INTERPOLATION
    IF(IP.EQ.0) THEN
       ! CALL POLATEV0(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI,&
       !      NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       CALL interpolate_bilinear_vector(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  BICUBIC INTERPOLATION
    ELSEIF(IP.EQ.1) THEN
       CALL interpolate_bicubic_vector(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  NEIGHBOR INTERPOLATION
    ELSEIF(IP.EQ.2) THEN
      CALL interpolate_neighbor_vector(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
           NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  BUDGET INTERPOLATION
    ELSEIF(IP.EQ.3) THEN
       CALL interpolate_budget_vector(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  SPECTRAL INTERPOLATION
    ELSEIF(IP.EQ.4) THEN
       CALL POLATEV4(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  NEIGHBOR-BUDGET INTERPOLATION
    ELSEIF(IP.EQ.6) THEN
       CALL POLATEV6(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  UNRECOGNIZED INTERPOLATION METHOD
    ELSE
       IF(KGDSO(1).GE.0) NO=0
       DO K=1,KM
          IBO(K)=1
          DO N=1,NO
             LO(N,K)=.FALSE.
             UO(N,K)=0.
             VO(N,K)=0.
          ENDDO
       ENDDO
       IRET=1
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(KGDSI(1).EQ.203) THEN
       KGDSI(11)=KGDSI11
    ENDIF
    IF(KGDSO(1).EQ.203) THEN
       KGDSO(11)=KGDSO11
    ENDIF
  END SUBROUTINE IPOLATEV_GRIB1


end module ipolatev_mod

