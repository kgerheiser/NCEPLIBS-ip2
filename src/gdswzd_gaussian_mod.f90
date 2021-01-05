MODULE GDSWZD_GAUSSIAN_MOD
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod 
  !$$$  MODULE DOCUMENTATION BLOCK
  !
  ! MODULE:  GDSWZD_GAUSSIAN_MOD  GDS WIZARD MODULE FOR GAUSSIAN CYLINDRICAL
  !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
  !
  ! ABSTRACT: - CONVERT FROM EARTH TO GRID COORDINATES OR VICE VERSA.
  !           - COMPUTE VECTOR ROTATION SINES AND COSINES.
  !           - COMPUTE MAP JACOBIANS.
  !           - COMPUTE GRID BOX AREA.
  !
  ! PROGRAM HISTORY LOG:
  !   2015-01-21  GAYNO   INITIAL VERSION FROM A MERGER OF
  !                       ROUTINES GDSWIZ04 AND GDSWZD04.
  !   2015-09-17  GAYNO   RENAME MODULE AS "GDSWZD_GAUSSIAN_MOD"
  !   2018-07-20  WESLEY  ADD THREADS.
  !
  ! USAGE:  "USE GDSWZD_GAUSSIAN_MOD"  THEN CALL THE PUBLIC DRIVER
  !         ROUTINE "GDSWZD_GAUSSIAN".
  !
  ! ATTRIBUTES:
  !   LANGUAGE: FORTRAN 90
  !
  !$$$
  !
  IMPLICIT NONE

  PRIVATE

  PUBLIC                         :: GDSWZD_GAUSSIAN, ip_gaussian_grid

  REAL,            PARAMETER     :: PI=3.14159265358979
  REAL,            PARAMETER     :: DPR=180./PI

  INTEGER                        :: J1, JH

  REAL,            ALLOCATABLE   :: BLAT(:)
  REAL                           :: DLON, RERTH
  REAL,            ALLOCATABLE   :: YLAT_ROW(:)

  type, extends(ip_grid) :: ip_gaussian_grid
     integer :: jh
     real :: dlon, rlat1, rlon1, rlon2, hi
     integer :: jg, jscan
   contains
     procedure :: init_grib1
     procedure :: init_grib2
  end type ip_gaussian_grid

CONTAINS

  subroutine init_grib1(self, g1_desc)
    class(ip_gaussian_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6
      self%eccen_squared = 0.0
      
      self%IM=KGDS(2)
      self%JM=KGDS(3)
      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      self%JG=KGDS(10)*2
      ISCAN=MOD(KGDS(11)/128,2)
      self%JSCAN=MOD(KGDS(11)/64,2)
      self%HI=(-1.)**ISCAN
      self%JH=(-1)**self%JSCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
    end associate
  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_gaussian_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)
      
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)
      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6
      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      self%RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
      self%JG=IGDTMPL(18)*2
      ISCAN=MOD(IGDTMPL(19)/128,2)
      self%JSCAN=MOD(IGDTMPL(19)/64,2)
      self%HI=(-1.)**ISCAN
      self%JH=(-1)**self%JSCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
    end associate
    
  end subroutine init_grib2

  SUBROUTINE GDSWZD_GAUSSIAN(grid,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GDSWZD_GAUSSIAN   GDS WIZARD FOR GAUSSIAN CYLINDRICAL
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB 2 GRID DEFINITION
    !           TEMPLATE (PASSED IN INTEGER FORM AS DECODED BY THE
    !           NCEP G2 LIBRARY) AND RETURNS ONE OF THE FOLLOWING:
    !             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
    !             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
    !           WORKS FOR GAUSSIAN CYLINDRICAL PROJECTIONS.
    !           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
    !           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
    !           OUTPUT ELEMENTS ARE SET TO FILL VALUES.
    !           THE ACTUAL NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
    !           OPTIONALLY, THE VECTOR ROTATIONS, THE MAP JACOBIANS AND
    !           THE GRID BOX AREAS MAY BE RETURNED AS WELL.  TO COMPUTE
    !           THE VECTOR ROTATIONS, THE OPTIONAL ARGUMENTS 'SROT' AND 'CROT'
    !           MUST BE PRESENT.  TO COMPUTE THE MAP JACOBIANS, THE
    !           OPTIONAL ARGUMENTS 'XLON', 'XLAT', 'YLON', 'YLAT' MUST BE PRESENT.
    !           TO COMPUTE THE GRID BOX AREAS, THE OPTIONAL ARGUMENT
    !           'AREA' MUST BE PRESENT.
    !
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !   97-10-20  IREDELL  INCLUDE MAP OPTIONS
    ! 1999-04-08  IREDELL  USE SUBROUTINE SPLAT
    ! 2001-06-18  IREDELL  CORRECT AREA COMPUTATION
    ! 2012-08-01  GAYNO    CORRECT AREA COMPUTATION AT POLE.
    !                      CORRECT YLAT COMPUTATION.
    ! 2015-01-21  GAYNO    MERGER OF GDSWIZ04 AND GDSWZD04.  MAKE
    !                      CROT,SORT,XLON,XLAT,YLON,YLAT AND AREA
    !                      OPTIONAL ARGUMENTS.  MAKE PART OF A MODULE.
    !                      MOVE VECTOR ROTATION, MAP JACOBIAN AND GRID
    !                      BOX AREA COMPUTATIONS TO SEPARATE SUBROUTINES.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
    !                      RENAME AS "GDSWZD_GAUSSIAN".
    ! 2018-07-20  WESLEY   ADD THREADS.
    !
    ! USAGE:    CALL GDSWZD_GAUSSIAN(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL,
    !    &                           XPTS,YPTS,RLON,RLAT,NRET,
    !    &                           CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                MUST BE "40" FOR GAUSSIAN GRIDS.
    !     IGDTMPL  - INTEGER (IGDTLEN) GRID DEFINITION TEMPLATE ARRAY.
    !                CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE FOR SECTION
    !                THREE:
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
    !     IGDTLEN  - INTEGER NUMBER OF ELEMENTS (19) OF THE GRID DEFINITION
    !                TEMPLATE ARRAY.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IOPT     - INTEGER OPTION FLAG
    !                (+1 TO COMPUTE EARTH COORDS OF SELECTED GRID COORDS)
    !                (-1 TO COMPUTE GRID COORDS OF SELECTED EARTH COORDS)
    !     NPTS     - INTEGER MAXIMUM NUMBER OF COORDINATES
    !     FILL     - REAL FILL VALUE TO SET INVALID OUTPUT DATA
    !                (MUST BE IMPOSSIBLE VALUE; SUGGESTED VALUE: -9999.)
    !     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT>0
    !     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT>0
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT<0
    !                (ACCEPTABLE RANGE: -360. TO 360.)
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT<0
    !                (ACCEPTABLE RANGE: -90. TO 90.)
    !
    !   OUTPUT ARGUMENT LIST:
    !     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT<0
    !     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT<0
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT>0
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT>0
    !     NRET     - INTEGER NUMBER OF VALID POINTS COMPUTED
    !     CROT     - REAL, OPTIONAL (NPTS) CLOCKWISE VECTOR ROTATION COSINES
    !     SROT     - REAL, OPTIONAL (NPTS) CLOCKWISE VECTOR ROTATION SINES
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !     XLON     - REAL, OPTIONAL (NPTS) DX/DLON IN 1/DEGREES
    !     XLAT     - REAL, OPTIONAL (NPTS) DX/DLAT IN 1/DEGREES
    !     YLON     - REAL, OPTIONAL (NPTS) DY/DLON IN 1/DEGREES
    !     YLAT     - REAL, OPTIONAL (NPTS) DY/DLAT IN 1/DEGREES
    !     AREA     - REAL, OPTIONAL (NPTS) AREA WEIGHTS IN M**2
    !
    ! EXTERNAL SUBPROGRAMS CALLED:
    !   SPLAT      COMPUTE LATITUDE FUNCTIONS
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    class(ip_gaussian_grid), intent(in) :: grid
    INTEGER,         INTENT(IN   ) :: IOPT, NPTS
    INTEGER,         INTENT(  OUT) :: NRET
    !
    REAL,            INTENT(IN   ) :: FILL
    REAL,            INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,            INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL,  INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL,  INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL,  INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                        :: JSCAN, IM, JM
    INTEGER                        :: J, JA, JG
    INTEGER                        :: ISCALE, N
    !
    LOGICAL                        :: LROT, LMAP, LAREA
    !
    REAL,            ALLOCATABLE   :: ALAT(:), ALAT_JSCAN(:)
    REAL,            ALLOCATABLE   :: ALAT_TEMP(:),BLAT_TEMP(:)
    REAL                           :: HI, RLATA, RLATB, RLAT1, RLON1, RLON2
    REAL                           :: XMAX, XMIN, YMAX, YMIN, YPTSA, YPTSB
    REAL                           :: WB
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    IF(PRESENT(CROT).AND.PRESENT(SROT))THEN
       LROT=.TRUE.
    ELSE
       LROT=.FALSE.
    ENDIF
    IF(PRESENT(XLON).AND.PRESENT(XLAT).AND.PRESENT(YLON).AND.PRESENT(YLAT))THEN
       LMAP=.TRUE.
    ELSE
       LMAP=.FALSE.
    ENDIF
    IF(PRESENT(AREA))THEN
       LAREA=.TRUE.
    ELSE
       LAREA=.FALSE.
    ENDIF
    
    IM=grid%im
    JM=grid%jm
   
    RLAT1=grid%rlat1
    RLON1=grid%rlon1
    RLON2=grid%rlon2
    
    JG=grid%jg
    JSCAN=grid%jscan
    HI=grid%hi
    
    JH=grid%jh
    DLON=grid%dlon
    rerth = grid%rerth
    
    ALLOCATE(ALAT_TEMP(JG))
    ALLOCATE(BLAT_TEMP(JG))
    CALL SPLAT(4,JG,ALAT_TEMP,BLAT_TEMP)
    ALLOCATE(ALAT(0:JG+1))
    ALLOCATE(BLAT(0:JG+1))
    !$OMP PARALLEL DO PRIVATE(JA) SCHEDULE(STATIC)
    DO JA=1,JG
       ALAT(JA)=DPR*ASIN(ALAT_TEMP(JA))
       BLAT(JA)=BLAT_TEMP(JA)
    ENDDO
    !$OMP END PARALLEL DO
    DEALLOCATE(ALAT_TEMP,BLAT_TEMP)
    ALAT(0)=180.-ALAT(1)
    ALAT(JG+1)=-ALAT(0)
    BLAT(0)=-BLAT(1)
    BLAT(JG+1)=BLAT(0)
    J1=1
    DO WHILE(J1.LT.JG.AND.RLAT1.LT.(ALAT(J1)+ALAT(J1+1))/2)
       J1=J1+1
    ENDDO
    IF(LMAP)THEN
       ALLOCATE(ALAT_JSCAN(JG))
       DO JA=1,JG
          ALAT_JSCAN(J1+JH*(JA-1))=ALAT(JA)
       ENDDO
       ALLOCATE(YLAT_ROW(0:JG+1))
       DO JA=2,(JG-1)
          YLAT_ROW(JA)=2.0/(ALAT_JSCAN(JA+1)-ALAT_JSCAN(JA-1))
       ENDDO
       YLAT_ROW(1)=1.0/(ALAT_JSCAN(2)-ALAT_JSCAN(1))
       YLAT_ROW(0)=YLAT_ROW(1)
       YLAT_ROW(JG)=1.0/(ALAT_JSCAN(JG)-ALAT_JSCAN(JG-1))
       YLAT_ROW(JG+1)=YLAT_ROW(JG)
       DEALLOCATE(ALAT_JSCAN)
    ENDIF
    XMIN=0
    XMAX=IM+1
    IF(IM.EQ.NINT(360/ABS(DLON))) XMAX=IM+2
    YMIN=0.5
    YMAX=JM+0.5
    NRET=0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
       !$OMP PARALLEL DO PRIVATE(N,J,WB,RLATA,RLATB) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             RLON(N)=MOD(RLON1+DLON*(XPTS(N)-1)+3600,360.)
             J=YPTS(N)
             WB=YPTS(N)-J
             RLATA=ALAT(J1+JH*(J-1))
             RLATB=ALAT(J1+JH*J)
             RLAT(N)=RLATA+WB*(RLATB-RLATA)
             NRET=NRET+1
             IF(LROT) CALL GAUSSIAN_VECT_ROT(CROT(N),SROT(N))
             IF(LMAP) CALL GAUSSIAN_MAP_JACOB(YPTS(N),&
                  XLON(N),XLAT(N),YLON(N),YLAT(N))
             IF(LAREA) CALL GAUSSIAN_GRID_AREA(YPTS(N),AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N,JA,YPTSA, YPTSB, WB) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          XPTS(N)=FILL
          YPTS(N)=FILL
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
             XPTS(N)=1+HI*MOD(HI*(RLON(N)-RLON1)+3600,360.)/DLON
             JA=MIN(INT((JG+1)/180.*(90-RLAT(N))),JG)
             IF(RLAT(N).GT.ALAT(JA)) JA=MAX(JA-2,0)
             IF(RLAT(N).LT.ALAT(JA+1)) JA=MIN(JA+2,JG)
             IF(RLAT(N).GT.ALAT(JA)) JA=JA-1
             IF(RLAT(N).LT.ALAT(JA+1)) JA=JA+1
             YPTSA=1+JH*(JA-J1)
             YPTSB=1+JH*(JA+1-J1)
             WB=(ALAT(JA)-RLAT(N))/(ALAT(JA)-ALAT(JA+1))
             YPTS(N)=YPTSA+WB*(YPTSB-YPTSA)
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT) CALL GAUSSIAN_VECT_ROT(CROT(N),SROT(N))
                IF(LMAP) CALL GAUSSIAN_MAP_JACOB(YPTS(N), &
                     XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL GAUSSIAN_GRID_AREA(YPTS(N),AREA(N))
             ELSE
                XPTS(N)=FILL
                YPTS(N)=FILL
             ENDIF
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    DEALLOCATE(ALAT, BLAT)
    IF (ALLOCATED(YLAT_ROW)) DEALLOCATE(YLAT_ROW)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GDSWZD_GAUSSIAN
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE GAUSSIAN_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GAUSSIAN_ERROR   ERROR HANDLER
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
    !           A "FILL" VALUE TO THE OUTPUT FIELDS.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:    CALL GAUSSIAN_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !
    !   INPUT ARGUMENT LIST:
    !     IOPT     - INTEGER OPTION FLAG
    !                (+1 TO COMPUTE EARTH COORDS OF SELECTED GRID COORDS)
    !                (-1 TO COMPUTE GRID COORDS OF SELECTED EARTH COORDS)
    !     NPTS     - INTEGER MAXIMUM NUMBER OF COORDINATES
    !     FILL     - REAL FILL VALUE TO SET INVALID OUTPUT DATA
    !                (MUST BE IMPOSSIBLE VALUE; SUGGESTED VALUE: -9999.)
    !   OUTPUT ARGUMENT LIST:
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT<0
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT<0
    !     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT>0
    !     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT>0
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: IOPT, NPTS
    !
    REAL,    INTENT(IN   ) :: FILL
    REAL,    INTENT(  OUT) :: RLAT(NPTS),RLON(NPTS)
    REAL,    INTENT(  OUT) :: XPTS(NPTS),YPTS(NPTS)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(IOPT>=0) THEN
       RLON=FILL
       RLAT=FILL
    ENDIF
    IF(IOPT<=0) THEN
       XPTS=FILL
       YPTS=FILL
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GAUSSIAN_ERROR
  !
  SUBROUTINE GAUSSIAN_VECT_ROT(CROT,SROT)
    !
    ! SUBPROGRAM:  GAUSSIAN_VECT_ROT   VECTOR ROTATION FIELDS FOR
    !                                  GAUSSIAN CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE VECTOR ROTATION SINES AND
    !           COSINES FOR A GAUSSIAN CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "GAUSSIAN_VECT_ROT"
    !
    ! USAGE:    CALL GAUSSIAN_VECT_ROT(CROT,SROT)
    !
    !   INPUT ARGUMENT LIST:
    !     NONE
    !
    !   OUTPUT ARGUMENT LIST:
    !     CROT     - CLOCKWISE VECTOR ROTATION COSINES (REAL)
    !     SROT     - CLOCKWISE VECTOR ROTATION SINES (REAL)
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    !
    IMPLICIT NONE

    REAL,                INTENT(  OUT) :: CROT, SROT

    CROT=1.0
    SROT=0.0

  END SUBROUTINE GAUSSIAN_VECT_ROT
  !
  SUBROUTINE GAUSSIAN_MAP_JACOB(YPTS, XLON, XLAT, YLON, YLAT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GAUSSIAN_MAP_JACOB  MAP JACOBIANS FOR
    !                                  GAUSSIAN CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE MAP JACOBIANS FOR
    !           A GAUSSIAN CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "GAUSSIAN_MAP_JACOB"
    !
    ! USAGE:  CALL GAUSSIAN_MAP_JACOB(YPTS,XLON,XLAT,YLON,YLAT)
    !
    !   INPUT ARGUMENT LIST:
    !     YPTS     - Y-INDEX OF GRID POINT (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     XLON     - DX/DLON IN 1/DEGREES (REAL)
    !     XLAT     - DX/DLAT IN 1/DEGREES (REAL)
    !     YLON     - DY/DLON IN 1/DEGREES (REAL)
    !     YLAT     - DY/DLAT IN 1/DEGREES (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$

    IMPLICIT NONE

    REAL,                INTENT(IN   ) :: YPTS
    REAL,                INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    XLON=1/DLON
    XLAT=0.
    YLON=0.
    YLAT=YLAT_ROW(NINT(YPTS))

  END SUBROUTINE GAUSSIAN_MAP_JACOB
  !
  SUBROUTINE GAUSSIAN_GRID_AREA(YPTS,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GAUSSIAN_GRID_AREA  GRID BOX AREA FOR
    !                                  GAUSSIAN CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE GRID BOX AREA FOR
    !           A GAUSSIAN CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "GAUSSIAN_GRID_AREA"
    !
    ! USAGE:  CALL GAUSSIAN_GRID_AREA(YPTS,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     YPTS     - Y-INDEX OF GRID POINT (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     AREA     - AREA WEIGHTS IN M**2 (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE

    REAL,            INTENT(IN   ) :: YPTS
    REAL,            INTENT(  OUT) :: AREA

    INTEGER                        :: J

    REAL                           :: WB, WLAT, WLATA, WLATB

    J = YPTS
    WB=YPTS-J
    WLATA=BLAT(J1+JH*(J-1))
    WLATB=BLAT(J1+JH*J)
    WLAT=WLATA+WB*(WLATB-WLATA)
    AREA=RERTH**2*WLAT*DLON/DPR

  END SUBROUTINE GAUSSIAN_GRID_AREA

END MODULE GDSWZD_GAUSSIAN_MOD
