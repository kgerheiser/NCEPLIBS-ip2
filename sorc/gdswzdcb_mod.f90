 MODULE GDSWZDCB_MOD
!$$$  MODULE DOCUMENTATION BLOCK
!
! MODULE:  GDSWZDCB_MOD  GDS WIZARD MODULE FOR ROTATED EQUIDISTANT 
!                        CYLINDRICAL GRIDS ("E" STAGGER).
!   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
!
! ABSTRACT: - CONVERT FROM EARTH TO GRID COORDINATES OR VICE VERSA.
!           - COMPUTE VECTOR ROTATION SINES AND COSINES.
!           - COMPUTE MAP JACOBIANS.
!           - COMPUTE GRID BOX AREA.
!
! PROGRAM HISTORY LOG:
!   2015-01-21  GAYNO   INITIAL VERSION FROM A MERGER OF
!                       ROUTINES GDSWIZCB AND GDSWZDCB.
!
! USAGE:  "USE GDSWZDCB_MOD"  THEN CALL THE PUBLIC DRIVER
!         ROUTINE "GDSWZDCB".
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
!
 IMPLICIT NONE

 PRIVATE

 PUBLIC                          :: GDSWZDCB

 INTEGER,         PARAMETER      :: KD=SELECTED_REAL_KIND(15,45)

 REAL(KIND=KD),   PARAMETER      :: PI=3.14159265358979_KD
 REAL(KIND=KD),   PARAMETER      :: DPR=180._KD/PI

 INTEGER                         :: IROT

 REAL(KIND=KD)                   :: CLAT, CLAT0, CLATR
 REAL(KIND=KD)                   :: CLON, DLATS, DLONS, RERTH
 REAL(KIND=KD)                   :: RLON0, SLAT, SLAT0, SLATR

 CONTAINS

 SUBROUTINE GDSWZDCB(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
                     XPTS,YPTS,RLON,RLAT,NRET, &
                     CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GDSWZDCB   GDS WIZARD FOR ROTATED EQUIDISTANT CYLINDRICAL
!   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
!
! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB 2 GRID DEFINITION
!           TEMPLATE (PASSED IN INTEGER FORM AS DECODED BY THE
!           NCEP G2 LIBRARY) AND RETURNS ONE OF THE FOLLOWING:
!             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
!             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
!           WORKS FOR E-STAGGERED ROTATED EQUIDISTANT CYLINDRICAL PROJECTIONS.
!           THE SCAN MODE (SECTION 3, OCTET 72, BITS 5-6) DETERMINE
!           WHETHER THIS IS AN "H" OR "V" GRID.  IF THE SELECTED
!           COORDINATES ARE MORE THAN ONE GRIDPOINT BEYOND THE
!           EDGES OF THE GRID DOMAIN, THEN THE RELEVANT OUTPUT ELEMENTS
!           ARE SET TO FILL VALUES. THE ACTUAL NUMBER OF VALID POINTS
!           COMPUTED IS RETURNED TOO. OPTIONALLY, THE VECTOR ROTATIONS,
!           THE MAP JACOBIANS AND THE GRID BOX AREAS MAY BE RETURNED.
!           TO COMPUTE THE VECTOR ROTATIONS, THE OPTIONAL ARGUMENTS
!           'SROT' AND 'CROT'  MUST BE PRESENT.  TO COMPUTE THE MAP
!           JACOBIANS, THE OPTIONAL ARGUMENTS 'XLON', 'XLAT', 
!           'YLON', 'YLAT' MUST BE PRESENT. TO COMPUTE THE GRID BOX
!           AREAS, THE OPTIONAL ARGUMENT 'AREA' MUST BE PRESENT.
!
! PROGRAM HISTORY LOG:
!   96-04-10  IREDELL
!   97-10-20  IREDELL  INCLUDE MAP OPTIONS
!   98-08-19  BALDWIN  MODIFY GDSWZDC9 FOR TYPE 203 ETA GRIDS
! 2003-06-11  IREDEL   INCREASE PRECISION
! 2012-08-02  GAYNO    INCREASE XMAX SO ON-GRID POINTS ARE NOT
!                      TAGGED AS OFF-GRID.
! 2015-01-21  GAYNO    MERGER OF GDSWIZCB AND GDSWZDCB.  MAKE
!                      CROT,SORT,XLON,XLAT,YLON,YLAT AND AREA
!                      OPTIONAL ARGUMENTS.  MAKE PART OF A MODULE.
!                      MOVE VECTOR ROTATION, MAP JACOBIAN AND GRID
!                      BOX AREA COMPUTATIONS TO SEPARATE SUBROUTINES.
! 2015-07-14  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
!                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
!
! USAGE:    CALL GDSWZDCB(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, 
!    &                    XPTS,YPTS,RLON,RLAT,NRET,
!    &                    CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
!
!   INPUT ARGUMENT LIST:
!     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
!                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
!                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
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
!                 (18): J-DIRECTION INCREMENT, OCTETS 68-71
!                 (19): SCANNING MODE, OCTET 72
!                 (20): LATITUDE OF SOUTHERN POLE OF PROJECTION, 
!                       OCTETS 73-76
!                 (21): LONGITUDE OF SOUTHERN POLE OF PROJECTION,
!                       OCTETS 77-80
!                 (22): ANGLE OF ROTATION OF PROJECTION, OCTS 81-84
!     IGDTLEN  - INTEGER NUMBER OF ELEMENTS (22) OF THE GRID DEFINITION
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
 INTEGER,         INTENT(IN   ) :: IGDTNUM, IGDTLEN
 INTEGER,         INTENT(IN   ) :: IGDTMPL(IGDTLEN)
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
 INTEGER                        :: IM, JM, IS1, N
 INTEGER                        :: ISCALE, ISCAN, KSCAN
 INTEGER                        :: I_OFFSET_ODD, I_OFFSET_EVEN
!
 LOGICAL                        :: LROT, LMAP, LAREA
!
 REAL(KIND=KD)                  :: RLAT0
 REAL(KIND=KD)                  :: CLONR
 REAL(KIND=KD)                  :: RLATR,RLONR
 REAL                           :: DUM1, DUM2, HI, HS
 REAL                           :: XMAX, XMIN, YMAX, YMIN, XPTF, YPTF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 IF(PRESENT(CROT)) CROT=FILL
 IF(PRESENT(SROT)) SROT=FILL
 IF(PRESENT(XLON)) XLON=FILL
 IF(PRESENT(XLAT)) XLAT=FILL
 IF(PRESENT(YLON)) YLON=FILL
 IF(PRESENT(YLAT)) YLAT=FILL
 IF(PRESENT(AREA)) AREA=FILL
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! IS THIS A ROTATED LAT/LON GRID?
 IF(IGDTNUM/=1)THEN
   CALL GDSWZDCB_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
   RETURN
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! IS THE EARTH RADIUS DEFINED?
 CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,DUM1,DUM2)
 RERTH=DUM1
 IF(RERTH<0.)THEN
   CALL GDSWZDCB_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
   RETURN
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ROUTINE ONLY WORKS FOR "E"-STAGGER GRIDS.
!   "V" GRID WHEN BIT 5 IS '1' AND BIT 6 IS '0'.
!   "H" GRID WHEN BIT 5 IS '0' AND BIT 6 IS '1'.
 I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
 I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
 IF(I_OFFSET_ODD==I_OFFSET_EVEN) THEN
   CALL GDSWZDCB_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
   RETURN
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ISCALE=IGDTMPL(10)*IGDTMPL(11)
 IF(ISCALE==0) ISCALE=1E6
 RLAT0=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
 RLAT0=RLAT0+90.0_KD
 RLON0=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)
 IROT=MOD(IGDTMPL(14)/8,2)
 IM=IGDTMPL(8)*2-1
 JM=IGDTMPL(9)
 KSCAN=I_OFFSET_ODD
 ISCAN=MOD(IGDTMPL(19)/128,2)
 HI=(-1.)**ISCAN
 SLAT0=SIN(RLAT0/DPR)
 CLAT0=COS(RLAT0/DPR)
 DLATS=FLOAT(IGDTMPL(18))/FLOAT(ISCALE)
 DLONS=FLOAT(IGDTMPL(17))/FLOAT(ISCALE)
! THE GRIB2 CONVENTION FOR "I" RESOLUTION IS TWICE WHAT THIS ROUTINE ASSUMES.
 DLONS=DLONS*0.5_KD
 IF(KSCAN.EQ.0) THEN
   IS1=(JM+1)/2
 ELSE
   IS1=JM/2
 ENDIF
 XMIN=0
 XMAX=IM+2
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
   DO N=1,NPTS
     XPTF=YPTS(N)+(XPTS(N)-IS1)
     YPTF=YPTS(N)-(XPTS(N)-IS1)+KSCAN
     IF(XPTF.GE.XMIN.AND.XPTF.LE.XMAX.AND. &
        YPTF.GE.YMIN.AND.YPTF.LE.YMAX) THEN
       HS=HI*SIGN(1.,XPTF-(IM+1)/2)
       RLONR=(XPTF-(IM+1)/2)*DLONS
       RLATR=(YPTF-(JM+1)/2)*DLATS
       CLONR=COS(RLONR/DPR)
       SLATR=SIN(RLATR/DPR)
       CLATR=COS(RLATR/DPR)
       SLAT=CLAT0*SLATR+SLAT0*CLATR*CLONR
       IF(SLAT.LE.-1) THEN
         CLAT=0.
         CLON=COS(RLON0/DPR)
         RLON(N)=0
         RLAT(N)=-90
       ELSEIF(SLAT.GE.1) THEN
         CLAT=0.
         CLON=COS(RLON0/DPR)
         RLON(N)=0
         RLAT(N)=90
       ELSE
         CLAT=SQRT(1-SLAT**2)
         CLON=(CLAT0*CLATR*CLONR-SLAT0*SLATR)/CLAT
         CLON=MIN(MAX(CLON,-1._KD),1._KD)
         RLON(N)=MOD(RLON0+HS*DPR*ACOS(CLON)+3600,360._KD)
         RLAT(N)=DPR*ASIN(SLAT)
       ENDIF
       NRET=NRET+1
       IF(LROT) CALL GDSWZDCB_VECT_ROT(RLON(N),CROT(N),SROT(N))
       IF(LMAP) CALL GDSWZDCB_MAP_JACOB(FILL, RLON(N), &
                                        XLON(N),XLAT(N),YLON(N),YLAT(N))
       IF(LAREA) CALL GDSWZDCB_GRID_AREA(FILL, AREA(N))
     ELSE
       RLON(N)=FILL
       RLAT(N)=FILL
     ENDIF
   ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
 ELSEIF(IOPT.EQ.-1) THEN
   DO N=1,NPTS
     IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
       HS=SIGN(1._KD,MOD(RLON(N)-RLON0+180+3600,360._KD)-180)
       CLON=COS((RLON(N)-RLON0)/DPR)
       SLAT=SIN(RLAT(N)/DPR)
       CLAT=COS(RLAT(N)/DPR)
       SLATR=CLAT0*SLAT-SLAT0*CLAT*CLON
       IF(SLATR.LE.-1) THEN
         CLATR=0.
         RLONR=0
         RLATR=-90
       ELSEIF(SLATR.GE.1) THEN
         CLATR=0.
         RLONR=0
         RLATR=90
       ELSE
         CLATR=SQRT(1-SLATR**2)
         CLONR=(CLAT0*CLAT*CLON+SLAT0*SLAT)/CLATR
         CLONR=MIN(MAX(CLONR,-1._KD),1._KD)
         RLONR=HS*DPR*ACOS(CLONR)
         RLATR=DPR*ASIN(SLATR)
       ENDIF
       XPTF=(IM+1)/2+RLONR/DLONS
       YPTF=(JM+1)/2+RLATR/DLATS
       IF(XPTF.GE.XMIN.AND.XPTF.LE.XMAX.AND. &
          YPTF.GE.YMIN.AND.YPTF.LE.YMAX) THEN
         XPTS(N)=IS1+(XPTF-(YPTF-KSCAN))/2
         YPTS(N)=(XPTF+(YPTF-KSCAN))/2
         NRET=NRET+1
         IF(LROT) CALL GDSWZDCB_VECT_ROT(RLON(N),CROT(N),SROT(N))
         IF(LMAP) CALL GDSWZDCB_MAP_JACOB(FILL, RLON(N), &
                                          XLON(N),XLAT(N),YLON(N),YLAT(N))
         IF(LAREA) CALL GDSWZDCB_GRID_AREA(FILL, AREA(N))
       ELSE
         XPTS(N)=FILL
         YPTS(N)=FILL
       ENDIF
     ELSE
       XPTS(N)=FILL
       YPTS(N)=FILL
     ENDIF
   ENDDO
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 END SUBROUTINE GDSWZDCB
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 SUBROUTINE GDSWZDCB_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GDSWZDCB_ERROR   GDSWZDCB ERROR HANDLER
!   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
!
! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
!           A "FILL" VALUE TO THE OUTPUT FIELDS.

! PROGRAM HISTORY LOG:
! 2015-07-13  GAYNO     INITIAL VERSION
!
! USAGE:    CALL GDSWZDCB_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
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
 END SUBROUTINE GDSWZDCB_ERROR
!
 SUBROUTINE GDSWZDCB_VECT_ROT(RLON, CROT, SROT)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GDSWZDCB_VECT_ROT   VECTOR ROTATION FIELDS FOR
!                                  ROTATED EQUIDISTANT CYLINDRICAL
!                                  GRIDS - "E" STAGGER.
!
!   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
!
! ABSTRACT: THIS SUBPROGRAM COMPUTES THE VECTOR ROTATION SINES AND
!           COSINES FOR A ROTATED EQUIDISTANT CYLINDRICAL GRID -
!           "E" STAGGER.
!
! PROGRAM HISTORY LOG:
! 2015-01-21  GAYNO    INITIAL VERSION
!
! USAGE:    CALL GDSWZDCB_VECT_ROT(RLON, CROT, SROT)
!
!   INPUT ARGUMENT LIST:
!     RLON     - LONGITUDE IN DEGREES (REAL)
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

 REAL         ,    INTENT(IN   ) :: RLON
 REAL         ,    INTENT(  OUT) :: CROT, SROT

 REAL(KIND=KD)                   :: SLON

 IF(IROT.EQ.1) THEN
   IF(CLATR.LE.0) THEN
     CROT=-SIGN(1._KD,SLATR*SLAT0)
     SROT=0.
   ELSE
     SLON=SIN((RLON-RLON0)/DPR)
     CROT=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
     SROT=SLAT0*SLON/CLATR
   ENDIF
 ELSE
   CROT=1.
   SROT=0.
 ENDIF

 END SUBROUTINE GDSWZDCB_VECT_ROT
!
 SUBROUTINE GDSWZDCB_MAP_JACOB(FILL, RLON, &
                               XLON, XLAT, YLON, YLAT)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GDSWZDCB_MAP_JACOB  MAP JACOBIANS FOR
!                                  ROTATED EQUIDISTANT CYLINDRICAL
!                                  GRIDS - "E" STAGGER.
!
!   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
!
! ABSTRACT: THIS SUBPROGRAM COMPUTES THE MAP JACOBIANS FOR
!           A ROTATED EQUIDISTANT CYLINDRICAL GRID - 
!           "E" STAGGER.
!
! PROGRAM HISTORY LOG:
! 2015-01-21  GAYNO    INITIAL VERSION
!
! USAGE:  CALL GDSWZDCB_MAP_JACOB(FILL,RLON,XLON,XLAT,YLON,YLAT)
!
!   INPUT ARGUMENT LIST:
!     FILL     - FILL VALUE FOR UNDEFINED POINTS (REAL)
!     RLON     - LONGITUDE IN DEGREES (REAL)
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

 REAL         ,    INTENT(IN   ) :: FILL, RLON
 REAL         ,    INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

 REAL(KIND=KD)                   :: SLON, TERM1, TERM2
 REAL(KIND=KD)                   :: XLATF, XLONF, YLATF, YLONF

 IF(CLATR.LE.0._KD) THEN
   XLON=FILL
   XLAT=FILL
   YLON=FILL
   YLAT=FILL
 ELSE
   SLON=SIN((RLON-RLON0)/DPR)
   TERM1=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
   TERM2=SLAT0*SLON/CLATR
   XLONF=TERM1*CLAT/(DLONS*CLATR)
   XLATF=-TERM2/(DLONS*CLATR)
   YLONF=TERM2*CLAT/DLATS
   YLATF=TERM1/DLATS
   XLON=XLONF-YLONF
   XLAT=XLATF-YLATF
   YLON=XLONF+YLONF
   YLAT=XLATF+YLATF
 ENDIF

 END SUBROUTINE GDSWZDCB_MAP_JACOB
!
 SUBROUTINE GDSWZDCB_GRID_AREA(FILL, AREA)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GDSWZDCB_GRID_AREA  GRID BOX AREA FOR
!                                  ROTATED EQUIDISTANT CYLINDRICAL
!                                  GRIDS - "E" STAGGER.
!
!   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
!
! ABSTRACT: THIS SUBPROGRAM COMPUTES THE GRID BOX AREA FOR
!           A ROTATED EQUIDISTANT CYLINDRICAL GRID -
!           "E" STAGGER.
!
! PROGRAM HISTORY LOG:
! 2015-01-21  GAYNO    INITIAL VERSION
!
! USAGE:  CALL GDSWZDCB_GRID_AREA(FILL,AREA)
!
!   INPUT ARGUMENT LIST:
!     FILL     - FILL VALUE FOR UNDEFINED POINTS (REAL)
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

 REAL,             INTENT(IN   ) :: FILL
 REAL,             INTENT(  OUT) :: AREA

 IF(CLATR.LE.0._KD) THEN
   AREA=FILL
 ELSE
   AREA=RERTH**2*CLATR*DLATS*DLONS*2/DPR**2
 ENDIF

 END SUBROUTINE GDSWZDCB_GRID_AREA

 END MODULE GDSWZDCB_MOD
