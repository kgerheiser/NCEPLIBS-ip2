module ip_polar_stereo_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use constants_mod, only: DPR, PI, PI2, PI4
  use earth_radius_mod
  implicit none

  private
  public :: ip_polar_stereo_grid

  type, extends(ip_grid) :: ip_polar_stereo_grid
     logical :: elliptical
     real :: rlat1, rlon1, orient, h, dxs, dys, slatr
     integer :: irot
   contains
     procedure :: init_grib1
     procedure :: init_grib2
     procedure :: gdswzd => gdswzd_polar_stereo
  end type ip_polar_stereo_grid

  INTEGER                         :: IROT

  REAL                            :: DE2, DXS, DYS
  REAL                            :: E2, RERTH , H, ORIENT

CONTAINS

  subroutine init_grib1(self, g1_desc)
    class(ip_polar_stereo_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    REAL, PARAMETER :: SLAT=60.0  ! standard latitude according grib1 standard

    real :: dx, dy, hi, hj
    integer :: iproj, iscan, jscan

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6
      self%eccen_squared = 0.00669437999013 !wgs84 datum
      self%ELLIPTICAL=MOD(KGDS(6)/64,2).EQ.1

      self%IM=KGDS(2)
      self%JM=KGDS(3)

      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3

      self%IROT=MOD(KGDS(6)/8,2)

      self%SLATR=SLAT/DPR

      self%ORIENT=KGDS(7)*1.E-3

      DX=KGDS(8)
      DY=KGDS(9)

      IPROJ=MOD(KGDS(10)/128,2)
      ISCAN=MOD(KGDS(11)/128,2)
      JSCAN=MOD(KGDS(11)/64,2)

      self%H=(-1.)**IPROJ
      HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)

      self%DXS=DX*HI
      self%DYS=DY*HJ

      self%iwrap= 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0
    end associate

  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_polar_stereo_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    real :: slat, dx, dy, hi, hj
    integer :: iproj, iscan, jscan

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

      self%ELLIPTICAL = self%eccen_squared > 0.0

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%RLAT1=FLOAT(IGDTMPL(10))*1.E-6
      self%RLON1=FLOAT(IGDTMPL(11))*1.E-6

      self%IROT=MOD(IGDTMPL(12)/8,2)

      SLAT=FLOAT(ABS(IGDTMPL(13)))*1.E-6
      self%SLATR=SLAT/DPR

      self%ORIENT=FLOAT(IGDTMPL(14))*1.E-6

      DX=FLOAT(IGDTMPL(15))*1.E-3
      DY=FLOAT(IGDTMPL(16))*1.E-3

      IPROJ=MOD(IGDTMPL(17)/128,2)
      ISCAN=MOD(IGDTMPL(18)/128,2)
      JSCAN=MOD(IGDTMPL(18)/64,2)

      self%H=(-1.)**IPROJ
      HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)

      self%DXS=DX*HI
      self%DYS=DY*HJ

      self%nscan = mod(igdtmpl(18) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
    end associate
  end subroutine init_grib2


  SUBROUTINE GDSWZD_POLAR_STEREO(self,IOPT,NPTS, &
       FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GDSWZD_POLAR_STEREO   GDS WIZARD FOR POLAR STEREOGRAPHIC 
    !                                    AZIMUTHAL
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB 2 GRID DEFINITION
    !           TEMPLATE (PASSED IN INTEGER FORM AS DECODED BY THE
    !           NCEP G2 LIBRARY) AND RETURNS ONE OF THE FOLLOWING:
    !             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
    !             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
    !           WORKS FOR POLAR STEREOGRAPHIC AZIMUTHAL PROJECTIONS.
    !           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
    !           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
    !           OUTPUT ELEMENTS ARE SET TO FILL VALUES.
    !           THE ACTUAL NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
    !           OPTIONALLY, THE VECTOR ROTATIONS, MAP JACOBIANS, AND
    !           GRID BOX AREAS MAY BE RETURNED AS WELL.  ROUTINE WORKS
    !           FOR BOTH SPHERICAL AND ELLIPTICAL EARTHS WITH THE
    !           EXCEPTION OF THE MAP JACOBIANS AND GRID BOX AREAS, WHICH
    !           ARE ONLY COMPUTED FOR SPHERICAL EARTHS.  TO COMPUTE
    !           THE VECTOR ROTATIONS, THE OPTIONAL ARGUMENTS 'SROT' AND 'CROT'
    !           MUST BE PRESENT.  TO COMPUTE THE MAP JACOBIANS, THE
    !           OPTIONAL ARGUMENTS 'XLON', 'XLAT', 'YLON', 'YLAT' MUST BE PRESENT.
    !           TO COMPUTE THE GRID BOX AREAS, THE OPTIONAL ARGUMENT
    !           'AREA' MUST BE PRESENT.
    !
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !   97-10-20  IREDELL  INCLUDE MAP OPTIONS
    !   09-05-13  GAYNO    ENSURE AREA ALWAYS POSITIVE
    ! 2015-01-21  GAYNO    MERGER OF GDSWIZ05 AND GDSWZD05.  MAKE
    !                      CROT,SORT,XLON,XLAT,YLON,YLAT AND AREA
    !                      OPTIONAL ARGUMENTS.  MAKE PART OF A MODULE.
    !                      MOVE VECTOR ROTATION, MAP JACOBIAN AND GRID
    !                      BOX AREA COMPUTATIONS TO SEPARATE SUBROUTINES.
    !                      INCLUDE OPTION FOR ELLIPTICAL EARTHS.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
    !                      RENAME ROUTINE AS "GDSWZD_POLAR_STEREO".
    ! 2018-07-20  WESLEY   ADD THREADING.
    !
    ! USAGE:   CALL GDSWZD_POLAR_STEREO(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,
    !    &                              FILL,XPTS,YPTS,RLON,RLAT,NRET,
    !    &                              CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                MUST BE "20" FOR POLAR STEREOGRAPHIC GRIDS.
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
    !     IGDTLEN  - INTEGER NUMBER OF ELEMENTS (18) OF THE GRID DEFINITION
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
    !                (PROPORTIONAL TO THE SQUARE OF THE MAP FACTOR)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !

    class(ip_polar_stereo_grid), intent(in) :: self
    INTEGER,          INTENT(IN   ) :: IOPT, NPTS
    INTEGER,          INTENT(  OUT) :: NRET
    !
    REAL,             INTENT(IN   ) :: FILL
    REAL,             INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,             INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL,   INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL,   INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL,   INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                         :: IM, JM
    INTEGER                         :: ITER, N
    !
    LOGICAL                         :: ELLIPTICAL, LROT, LMAP, LAREA
    !
    REAL                            :: ALAT, ALAT1, ALONG, DIFF
    REAL                            :: DI, DJ, DE
    REAL                            :: DR, E, E_OVER_2
    REAL                            :: MC, SLAT, SLATR
    REAL                            :: RLAT1, RLON1, RHO, T, TC
    REAL                            :: XMAX, XMIN, YMAX, YMIN
    REAL                            :: XP, YP, DR2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL

    elliptical = self%elliptical
    IM=self%im
    JM=self%jm

    RLAT1=self%rlat1
    RLON1=self%rlon1

    IROT=self%irot
    SLATR=self%slatr
    ORIENT=self%orient

    H=self%h
    DXS=self%dxs
    DYS=self%dys

    rerth = self%rerth
    e2 = self%eccen_squared
    !
    ! FIND X/Y OF POLE
    IF (.NOT.ELLIPTICAL) THEN
       DE=(1.+SIN(SLATR))*RERTH
       DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
       XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/DXS
       YP=1+COS((RLON1-ORIENT)/DPR)*DR/DYS
       DE2=DE**2
    ELSE
       E=SQRT(E2)
       E_OVER_2=E*0.5
       ALAT=H*RLAT1/DPR
       ALONG = (RLON1-ORIENT)/DPR
       T=TAN(PI4-ALAT/2.)/((1.-E*SIN(ALAT))/  &
            (1.+E*SIN(ALAT)))**(E_OVER_2)
       TC=TAN(PI4-SLATR/2.)/((1.-E*SIN(SLATR))/  &
            (1.+E*SIN(SLATR)))**(E_OVER_2)
       MC=COS(SLATR)/SQRT(1.0-E2*(SIN(SLATR)**2))
       RHO=RERTH*MC*T/TC
       YP = 1.0 + RHO*COS(H*ALONG)/DYS
       XP = 1.0 - RHO*SIN(H*ALONG)/DXS
    ENDIF ! ELLIPTICAL
    XMIN=0
    XMAX=IM+1
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
       IF(.NOT.ELLIPTICAL)THEN
          !$OMP PARALLEL DO PRIVATE(N,DI,DJ,DR2) REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                DI=(XPTS(N)-XP)*DXS
                DJ=(YPTS(N)-YP)*DYS
                DR2=DI**2+DJ**2
                IF(DR2.LT.DE2*1.E-6) THEN
                   RLON(N)=0.
                   RLAT(N)=H*90.
                ELSE
                   RLON(N)=MOD(ORIENT+H*DPR*ATAN2(DI,-DJ)+3600,360.)
                   RLAT(N)=H*DPR*ASIN((DE2-DR2)/(DE2+DR2))
                ENDIF
                NRET=NRET+1
                IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
                IF(LMAP) CALL POLAR_STEREO_MAP_JACOB(RLON(N),RLAT(N),DR2, &
                     XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL POLAR_STEREO_GRID_AREA(RLAT(N),DR2,AREA(N))
             ELSE
                RLON(N)=FILL
                RLAT(N)=FILL
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ELSE ! ELLIPTICAL
          !$OMP PARALLEL DO PRIVATE(N,DI,DJ,RHO,T,ALONG,ALAT1,ALAT,DIFF) &
          !$OMP& REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND.  &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                DI=(XPTS(N)-XP)*DXS
                DJ=(YPTS(N)-YP)*DYS
                RHO=SQRT(DI*DI+DJ*DJ)
                T=(RHO*TC)/(RERTH*MC)
                IF(ABS(YPTS(N)-YP)<0.01)THEN
                   IF(DI>0.0) ALONG=ORIENT+H*90.0
                   IF(DI<=0.0) ALONG=ORIENT-H*90.0
                ELSE
                   ALONG=ORIENT+H*ATAN(DI/(-DJ))*DPR
                   IF(DJ>0) ALONG=ALONG+180.
                END IF
                ALAT1=PI2-2.0*ATAN(T)
                DO ITER=1,10
                   ALAT = PI2 - 2.0*ATAN(T*(((1.0-E*SIN(ALAT1))/  &
                        (1.0+E*SIN(ALAT1)))**(E_OVER_2)))
                   DIFF = ABS(ALAT-ALAT1)*DPR
                   IF (DIFF < 0.000001) EXIT
                   ALAT1=ALAT
                ENDDO
                RLAT(N)=H*ALAT*DPR
                RLON(N)=ALONG
                IF(RLON(N)<0.0) RLON(N)=RLON(N)+360.
                IF(RLON(N)>360.0) RLON(N)=RLON(N)-360.0
                NRET=NRET+1
                IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
             ELSE
                RLON(N)=FILL
                RLAT(N)=FILL
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF ! ELLIPTICAL
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       IF(.NOT.ELLIPTICAL)THEN
          !$OMP PARALLEL DO PRIVATE(N,DR,DR2) REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90.AND. &
                  H*RLAT(N).NE.-90) THEN
                DR=DE*TAN((90-H*RLAT(N))/2/DPR)
                DR2=DR**2
                XPTS(N)=XP+H*SIN((RLON(N)-ORIENT)/DPR)*DR/DXS
                YPTS(N)=YP-COS((RLON(N)-ORIENT)/DPR)*DR/DYS
                IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                     YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                   NRET=NRET+1
                   IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
                   IF(LMAP) CALL POLAR_STEREO_MAP_JACOB(RLON(N),RLAT(N),DR2, &
                        XLON(N),XLAT(N),YLON(N),YLAT(N))
                   IF(LAREA) CALL POLAR_STEREO_GRID_AREA(RLAT(N),DR2,AREA(N))
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
       ELSE  ! ELLIPTICAL CASE
          !$OMP PARALLEL DO PRIVATE(N,ALAT,ALONG,T,RHO) REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90.AND.  &
                  H*RLAT(N).NE.-90) THEN
                ALAT = H*RLAT(N)/DPR
                ALONG = (RLON(N)-ORIENT)/DPR
                T=TAN(PI4-ALAT*0.5)/((1.-E*SIN(ALAT))/  &
                     (1.+E*SIN(ALAT)))**(E_OVER_2)
                RHO=RERTH*MC*T/TC
                XPTS(N)= XP + RHO*SIN(H*ALONG) / DXS
                YPTS(N)= YP - RHO*COS(H*ALONG) / DYS
                IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND.  &
                     YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                   NRET=NRET+1
                   IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
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
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GDSWZD_POLAR_STEREO
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE POLAR_STEREO_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLAR_STEREO_ERROR    ERROR HANDLER
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
    !           A "FILL" VALUE TO THE OUTPUT FIELDS.

    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:    CALL POLAR_STEREO_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
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
  END SUBROUTINE POLAR_STEREO_ERROR
  !
  SUBROUTINE POLAR_STEREO_VECT_ROT(RLON, CROT, SROT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLAR_STEREO_VECT_ROT   VECTOR ROTATION FIELDS FOR
    !                                      POLAR STEREOGRAPHIC GRIDS.
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE VECTOR ROTATION SINES AND
    !           COSINES FOR A POLAR STEREOGRAPHIC AZIMUTHAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "POLAR_STEREO_VECT_ROT"
    !
    ! USAGE:    CALL POLAR_STEREO_VECT_ROT(RLON,CROT,SROT)
    !
    !   INPUT ARGUMENT LIST:
    !     RLON     - GRID POINT LONGITUDE IN DEGREES (REAL)
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

    REAL,             INTENT(IN   ) :: RLON
    REAL,             INTENT(  OUT) :: CROT, SROT

    IF(IROT.EQ.1) THEN
       CROT=H*COS((RLON-ORIENT)/DPR)
       SROT=SIN((RLON-ORIENT)/DPR)
    ELSE
       CROT=1.
       SROT=0.
    ENDIF

  END SUBROUTINE POLAR_STEREO_VECT_ROT
  !
  SUBROUTINE POLAR_STEREO_MAP_JACOB(RLON,RLAT,DR2,XLON,XLAT,YLON,YLAT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLAR_STEREO_MAP_JACOB  MAP JACOBIANS FOR
    !                                      POLAR STEREOGRAPHIC GRIDS.
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE MAP JACOBIANS FOR
    !           A POLAR STEREOGRAPHIC AZIMUTHAL GRID (SPHERICAL
    !           EARTH).
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "POLAR_STEREO_MAP_JACOB"
    ! 2018-07-20  WESLEY   PASS IN DR2 FOR THREADING.
    !
    ! USAGE:  CALL POLAR_STEREO_MAP_JACOB(RLON,RLAT,DR2,XLON,XLAT,YLON,YLAT)
    !
    !   INPUT ARGUMENT LIST:
    !     RLON     - LONGITUDE IN DEGREES (REAL)
    !     RLAT     - LATITUDE IN DEGREES (REAL)
    !     DR2      - SQUARED DISTANCE FROM POLE (REAL)
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
    !
    IMPLICIT NONE

    REAL,             INTENT(IN   ) :: RLON, RLAT, DR2
    REAL,             INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    REAL                            :: CLAT, DE, DR

    IF(DR2.LT.DE2*1.E-6) THEN
       DE=SQRT(DE2)
       XLON=0.
       XLAT=-SIN((RLON-ORIENT)/DPR)/DPR*DE/DXS/2
       YLON=0.
       YLAT=H*COS((RLON-ORIENT)/DPR)/DPR*DE/DYS/2
    ELSE
       DR=SQRT(DR2)
       CLAT=COS(RLAT/DPR)
       XLON=H*COS((RLON-ORIENT)/DPR)/DPR*DR/DXS
       XLAT=-SIN((RLON-ORIENT)/DPR)/DPR*DR/DXS/CLAT
       YLON=SIN((RLON-ORIENT)/DPR)/DPR*DR/DYS
       YLAT=H*COS((RLON-ORIENT)/DPR)/DPR*DR/DYS/CLAT
    ENDIF

  END SUBROUTINE POLAR_STEREO_MAP_JACOB
  !
  SUBROUTINE POLAR_STEREO_GRID_AREA(RLAT, DR2, AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLAR_STEREO_GRID_AREA  GRID BOX AREA FOR
    !                                      POLAR STEREOGRAPHIC GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE GRID BOX AREA FOR
    !           A POLAR STEREOGRAPHIC AZIMUTHAL GRID (SPHERICAL
    !           EARTH).
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "POLAR_STEREO_GRID_AREA".
    ! 2018-07-20  WESLEY   PASS IN DR2 FOR THREADING.
    !
    ! USAGE:  CALL POLAR_STEREO_GRID_AREA(RLAT,DR2,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     RLAT     - LATITUDE OF GRID POINT IN DEGREES (REAL)
    !     DR2      - SQUARED DISTANCE FROM POLE (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     AREA     - AREA WEIGHTS IN M**2 (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    !
    IMPLICIT NONE

    REAL,             INTENT(IN   ) :: RLAT, DR2
    REAL,             INTENT(  OUT) :: AREA

    REAL                            :: CLAT

    IF(DR2.LT.DE2*1.E-6) THEN
       AREA=RERTH**2*ABS(DXS)*ABS(DYS)*4/DE2
    ELSE
       CLAT=COS(RLAT/DPR)
       AREA=RERTH**2*CLAT**2*ABS(DXS)*ABS(DYS)/DR2
    ENDIF

  END SUBROUTINE POLAR_STEREO_GRID_AREA


end module ip_polar_stereo_grid_mod

