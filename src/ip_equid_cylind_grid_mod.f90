module ip_equid_cylind_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  implicit none

  private
  public :: ip_equid_cylind_grid

  type, extends(ip_grid) :: ip_equid_cylind_grid
     real :: hi, rlat1, rlon1, rlat2, rlon2
     real :: dlat, dlon
   contains
     procedure :: init_grib1
     procedure :: init_grib2
     procedure :: gdswzd => gdswzd_equid_cylind
  end type ip_equid_cylind_grid

  REAL                    :: DLAT ! GRID RESOLUTION IN DEGREES N/S DIRECTION
  REAL                    :: DLON ! GRID RESOLUTION IN DEGREES E/W DIRECTION
  REAL                    :: RERTH

contains

  subroutine init_grib1(self, g1_desc)
    class(ip_equid_cylind_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan

    associate(kgds => g1_desc%gds)
      self%IM=KGDS(2)
      self%JM=KGDS(3)
      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLAT2=KGDS(7)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      ISCAN=MOD(KGDS(11)/128,2)
      self%HI=(-1.)**ISCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DLAT=(self%RLAT2-self%RLAT1)/(self%JM-1)

      ! defaults
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%kscan = 0

      self%iwrap = nint(360/abs(self%dlon))

      if(self%im < self%iwrap) self%iwrap=0
      self%jwrap1 = 0
      self%jwrap2 = 0
      if(self%iwrap > 0 .and. mod(self%iwrap,2) == 0) then
         if(abs(self%rlat1) > 90-0.25*self%dlat) then
            self%jwrap1 = 2
         elseif(abs(self%rlat1) > 90-0.75*self%dlat) then
            self%jwrap1 = 1
         endif
         if(abs(self%rlat2) > 90-0.25*self%dlat) then
            self%jwrap2 = 2 * self%jm
         elseif(abs(self%rlat2) > 90-0.75*self%dlat) then
            self%jwrap2 = 2 * self%jm+1
         endif
      endif

      self%rerth = 6.3712E6
      self%eccen_squared = 0.0
    end associate

  end subroutine init_grib1


  subroutine init_grib2(self, g2_desc)
    class(ip_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)
      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6
      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      self%RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
      self%RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
      ISCAN=MOD(IGDTMPL(19)/128,2)
      self%HI=(-1.)**ISCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DLAT=(self%RLAT2-self%RLAT1)/(self%JM-1)

      self%nscan = MOD(IGDTMPL(19)/32,2)
      self%kscan = 0
      self%iwrap = NINT(360/ABS(self%DLON))

      if(self%im.lt.self%iwrap) self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0

      if(self%im < self%iwrap) self%iwrap=0
      self%jwrap1 = 0
      self%jwrap2 = 0
      if(self%iwrap > 0 .and. mod(self%iwrap,2) == 0) then
         if(abs(self%rlat1) > 90-0.25*self%dlat) then
            self%jwrap1 = 2
         elseif(abs(self%rlat1) > 90-0.75*self%dlat) then
            self%jwrap1 = 1
         endif
         if(abs(self%rlat2) > 90-0.25*self%dlat) then
            self%jwrap2 = 2 * self%jm
         elseif(abs(self%rlat2) > 90-0.75*self%dlat) then
            self%jwrap2 = 2 * self%jm+1
         endif
      endif

      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

    end associate
  end subroutine init_grib2


  SUBROUTINE GDSWZD_EQUID_CYLIND(self,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET,  &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GDSWZD_EQUID_CYLIND   GDS WIZARD FOR EQUIDISTANT CYLINDRICAL
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB 2 GRID DEFINITION
    !           TEMPLATE (PASSED IN INTEGER FORM AS DECODED BY THE
    !           NCEP G2 LIBRARY) AND RETURNS ONE OF THE FOLLOWING:
    !             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
    !             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
    !           WORKS FOR EQUIDISTANT CYLINDRICAL PROJECTIONS.
    !           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
    !           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
    !           OUTPUT ELEMENTS ARE SET TO FILL VALUES.  THE ACTUAL
    !           NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
    !           OPTIONALLY, THE VECTOR ROTATIONS, MAP JACOBIANS AND
    !           GRID BOX AREAS MAY BE RETURNED.  TO COMPUTE THE VECTOR
    !           ROTATIONS, THE OPTIONAL ARGUMENTS 'SROT' AND 'CROT'
    !           MUST BE PRESENT.  TO COMPUTE THE MAP JACOBIANS, THE
    !           OPTIONAL ARGUMENTS 'XLON', 'XLAT', 'YLON', 'YLAT' MUST 
    !           BE PRESENT. TO COMPUTE THE GRID BOX AREAS, THE OPTIONAL
    !           ARGUMENT 'AREA' MUST BE PRESENT.
    !
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !   97-10-20  IREDELL  INCLUDE MAP OPTIONS
    ! 2015-01-21  GAYNO    MERGER OF GDSWIZ00 AND GDSWZD00.  MAKE
    !                      CROT,SORT,XLON,XLAT,YLON,YLAT AND AREA
    !                      OPTIONAL ARGUMENTS.  MAKE PART OF A MODULE. 
    !                      MOVE VECTOR ROTATION, MAP JACOBIAN AND GRID
    !                      BOX AREA COMPUTATIONS TO SEPARATE SUBROUTINES.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
    !                      RENAME ROUTINE AS "GDSWZD_EQUID_CYLIND".
    ! 2018-07-20  WESLEY   ADD THREADS.
    !
    ! USAGE:    CALL GDSWZD_EQUID_CYLIND(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,
    !    &                  FILL,XPTS,YPTS,RLON,RLAT,NRET,
    !    &                  CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    ! 
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                MUST BE "0" FOR EQUIDISTANT CYLINDRICAL 
    !                PROJECTIONS.
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
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    class(ip_equid_cylind_grid), intent(in) :: self
    INTEGER,             INTENT(IN   ) :: IOPT, NPTS
    INTEGER,             INTENT(  OUT) :: NRET
    !
    REAL,                INTENT(IN   ) :: FILL
    REAL,                INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,                INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL,      INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL,      INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL,      INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                            :: IM, JM, N
    !
    LOGICAL                            :: LROT, LMAP, LAREA
    !
    REAL                               :: HI, RLAT1, RLON1, RLAT2, RLON2
    REAL                               :: XMAX, XMIN, YMAX, YMIN
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL

    IM=self%im
    JM=self%jm

    RLAT1=self%rlat1
    RLON1=self%rlon1
    RLAT2=self%rlat2
    RLON2=self%rlon2

    HI=self%hi

    rerth = self%rerth
    dlat = self%dlat
    dlon = self%dlon

    XMIN=0
    XMAX=IM+1
    IF(IM.EQ.NINT(360/ABS(DLON))) XMAX=IM+2
    YMIN=0
    YMAX=JM+1
    NRET=0
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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             RLON(N)=MOD(RLON1+DLON*(XPTS(N)-1)+3600,360.)
             RLAT(N)=MIN(MAX(RLAT1+DLAT*(YPTS(N)-1),-90.),90.)
             NRET=NRET+1
             IF(LROT)  CALL EQUID_CYLIND_VECT_ROT(CROT(N),SROT(N))
             IF(LMAP)  CALL EQUID_CYLIND_MAP_JACOB(XLON(N),XLAT(N),YLON(N),YLAT(N))
             IF(LAREA) CALL EQUID_CYLIND_GRID_AREA(RLAT(N),AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
             XPTS(N)=1+HI*MOD(HI*(RLON(N)-RLON1)+3600,360.)/DLON
             YPTS(N)=1+(RLAT(N)-RLAT1)/DLAT
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND.  &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT)  CALL EQUID_CYLIND_VECT_ROT(CROT(N),SROT(N))
                IF(LMAP)  CALL EQUID_CYLIND_MAP_JACOB(XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL EQUID_CYLIND_GRID_AREA(RLAT(N),AREA(N))
             ELSE
                XPTS(N)=FILL
                YPTS(N)=FILL
             ENDIF
          ELSE
             XPTS(N)=FILL
             YPTS(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GDSWZD_EQUID_CYLIND
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  EQUID_CYLIND_ERROR   ERROR HANDLER
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
    !           A "FILL" VALUE TO THE OUTPUT FIELDS.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:    CALL EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
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
  END SUBROUTINE EQUID_CYLIND_ERROR
  !
  SUBROUTINE EQUID_CYLIND_VECT_ROT(CROT,SROT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  EQUID_CYLIND_VECT_ROT   VECTOR ROTATION FIELDS FOR 
    !                                      EQUIDISTANT CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE VECTOR ROTATION SINES AND
    !           COSINES FOR AN EQUIDISTANT CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "EQUID_CYLIND_VECT_ROT"
    !
    ! USAGE:    CALL EQUID_CYLIND_VECT_ROT(CROT,SROT)
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
    IMPLICIT NONE

    REAL,                INTENT(  OUT) :: CROT, SROT

    CROT=1.0
    SROT=0.0

  END SUBROUTINE EQUID_CYLIND_VECT_ROT
  !
  SUBROUTINE EQUID_CYLIND_MAP_JACOB(XLON,XLAT,YLON,YLAT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  EQUID_CYLIND_MAP_JACOB  MAP JACOBIANS FOR 
    !                                      EQUIDISTANT CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE MAP JACOBIANS FOR
    !           AN EQUIDISTANT CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "EQUID_CYLIND_MAP_JACOB"
    !
    ! USAGE:  CALL EQUID_CYLIND_MAP_JACOB(XLON,XLAT,YLON,YLAT)
    ! 
    !   INPUT ARGUMENT LIST:
    !     NONE
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

    REAL,                INTENT(  OUT) :: XLON,XLAT,YLON,YLAT

    XLON=1.0/DLON
    XLAT=0.
    YLON=0.
    YLAT=1.0/DLAT

  END SUBROUTINE EQUID_CYLIND_MAP_JACOB
  !
  SUBROUTINE EQUID_CYLIND_GRID_AREA(RLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  EQUID_CYLIND_GRID_AREA  GRID BOX AREA FOR 
    !                                      EQUIDISTANT CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE GRID BOX AREA FOR
    !           AN EQUIDISTANT CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "EQUID_CYLIND_GRID_AREA"
    !
    ! USAGE:  CALL EQUID_CYLIND_GRID_AREA(RLAT,AREA)
    ! 
    !   INPUT ARGUMENT LIST:
    !     RLAT     - LATITUDE OF GRID POINT IN DEGREES (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     AREA     - AREA WEIGHTS IN M**2 (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE

    REAL,                INTENT(IN   ) :: RLAT
    REAL,                INTENT(  OUT) :: AREA

    REAL,                PARAMETER     :: PI=3.14159265358979
    REAL,                PARAMETER     :: DPR=180./PI

    REAL                               :: DSLAT, RLATU, RLATD

    RLATU=MIN(MAX(RLAT+DLAT/2,-90.),90.)
    RLATD=MIN(MAX(RLAT-DLAT/2,-90.),90.)
    DSLAT=SIN(RLATU/DPR)-SIN(RLATD/DPR)
    AREA=RERTH**2*ABS(DSLAT*DLON)/DPR

  END SUBROUTINE EQUID_CYLIND_GRID_AREA

end module ip_equid_cylind_grid_mod

