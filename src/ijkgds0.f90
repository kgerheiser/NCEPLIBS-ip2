module ijkgds0_mod
  use ip_grid_descriptor_mod
  implicit none
contains

  SUBROUTINE IJKGDS0(grid_desc,IJKGDSA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  IJKGDS0    SET UP PARAMETERS FOR FUNCTION IJKGDS1
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB 2 GRID DEFINITION TEMPLATE
    !           AND RETURNS A NAVIGATION PARAMETER ARRAY TO ALLOW FUNCTION
    !           IJKGDS1 TO DECODE THE FIELD POSITION FOR A GIVEN GRID POINT.
    !
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !   97-03-11  IREDELL  ALLOWED HEMISPHERIC GRIDS TO WRAP OVER ONE POLE
    !   98-07-13  BALDWIN  ADD 2D STAGGERED ETA GRID INDEXING (203)
    ! 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO
    ! 2015-01-27  GAYNO    REMOVE REFERENCES TO OBSOLETE NCEP GRIDS 201
    !                      AND 202.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
    !
    ! USAGE:    CALL IJKGDS0(IGDTNUM,IGDTMPL,IGDTLEN,IJKGDSA)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE:
    !                  00 - EQUIDISTANT CYLINDRICAL
    !                  01 - ROTATED EQUIDISTANT CYLINDRICAL.  "E"
    !                       AND NON-"E" STAGGERED
    !                  10 - MERCATOR CYCLINDRICAL
    !                  20 - POLAR STEREOGRAPHIC AZIMUTHAL
    !                  30 - LAMBERT CONFORMAL CONICAL
    !                  40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
    !     IGDTMPL  - INTEGER (IGDTLEN) GRID DEFINITION TEMPLATE ARRAY.
    !                CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE
    !                (SECTION 3 INFO):
    !                ALL PROJECTIONS:
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
    !     IGDTLEN  - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !
    !   OUTPUT ARGUMENT LIST:
    !     IJKGDSA  - INTEGER (20) NAVIGATION PARAMETER ARRAY
    !                IJKGDSA(1) IS NUMBER OF X POINTS
    !                IJKGDSA(2) IS NUMBER OF Y POINTS
    !                IJKGDSA(3) IS X WRAPAROUND INCREMENT
    !                           (0 IF NO WRAPAROUND)
    !                IJKGDSA(4) IS Y WRAPAROUND LOWER PIVOT POINT
    !                           (0 IF NO WRAPAROUND)
    !                IJKGDSA(5) IS Y WRAPAROUND UPPER PIVOT POINT
    !                           (0 IF NO WRAPAROUND)
    !                IJKGDSA(6) IS SCANNING MODE
    !                           (0 IF X FIRST THEN Y; 1 IF Y FIRST THEN X;
    !                            3 IF STAGGERED DIAGONAL LIKE ARAKAWA "E")
    !                IJKGDSA(7) IS MASS/WIND FLAG FOR STAGGERED DIAGONAL
    !                           (0 IF MASS; 1 IF WIND)
    !                IJKGDSA(8:20) ARE UNUSED AT THE MOMENT
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    class(ip_grid_descriptor), intent(in) :: grid_desc
    INTEGER,             INTENT(  OUT) :: IJKGDSA(20)
    !
    INTEGER                            :: IM, JM, IWRAP, JG
    INTEGER                            :: I_OFFSET_ODD, I_OFFSET_EVEN
    INTEGER                            :: ISCAN, KSCAN, NSCAN
    INTEGER                            :: JWRAP1, JWRAP2, ISCALE
    !
    REAL                               :: DLAT, DLON
    REAL                               :: RLAT1, RLAT2
    REAL                               :: RLON1, RLON2

    select type(grid_desc)
    type is(grib1_descriptor)
       error stop "not implemented yet"
    type is(grib2_descriptor)
       associate(igdtnum => grid_desc%gdt_num, igdtmpl => grid_desc%gdt_tmpl)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !  SET EXCEPTIONAL VALUES FOR LAT/LON PROJECTION
         IF(IGDTNUM.EQ.0) THEN
            IM=IGDTMPL(8)
            JM=IGDTMPL(9)
            NSCAN=MOD(IGDTMPL(19)/32,2)
            KSCAN=0
            ISCALE=IGDTMPL(10)*IGDTMPL(11)
            IF(ISCALE==0) ISCALE=10**6
            RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
            RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
            ISCAN=MOD(IGDTMPL(19)/128,2)
            IF(ISCAN.EQ.0) THEN
               DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
            ELSE
               DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
            ENDIF
            IWRAP=NINT(360/ABS(DLON))
            IF(IM.LT.IWRAP) IWRAP=0
            JWRAP1=0
            JWRAP2=0
            IF(IWRAP.GT.0.AND.MOD(IWRAP,2).EQ.0) THEN
               RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
               RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
               DLAT=ABS(RLAT2-RLAT1)/(JM-1)
               IF(ABS(RLAT1).GT.90-0.25*DLAT) THEN
                  JWRAP1=2
               ELSEIF(ABS(RLAT1).GT.90-0.75*DLAT) THEN
                  JWRAP1=1
               ENDIF
               IF(ABS(RLAT2).GT.90-0.25*DLAT) THEN
                  JWRAP2=2*JM
               ELSEIF(ABS(RLAT2).GT.90-0.75*DLAT) THEN
                  JWRAP2=2*JM+1
               ENDIF
            ENDIF
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  SET EXCEPTIONAL VALUES FOR ROTATED LAT/LON
         ELSEIF(IGDTNUM.EQ.1) THEN
            I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
            I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
            IM=IGDTMPL(8)
            JM=IGDTMPL(9)
            IWRAP=0
            JWRAP1=0
            JWRAP2=0
            KSCAN=0
            NSCAN=MOD(IGDTMPL(19)/32,2)
            IF(I_OFFSET_ODD/=I_OFFSET_EVEN)THEN
               KSCAN=I_OFFSET_ODD
               NSCAN=3
            ENDIF
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  SET EXCEPTIONAL VALUES FOR MERCATOR.
         ELSEIF(IGDTNUM.EQ.10) THEN
            IM=IGDTMPL(8)
            JM=IGDTMPL(9)
            RLON1=FLOAT(IGDTMPL(11))*1.0E-6
            RLON2=FLOAT(IGDTMPL(15))*1.0E-6
            ISCAN=MOD(IGDTMPL(16)/128,2)
            IF(ISCAN.EQ.0) THEN
               DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
            ELSE
               DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
            ENDIF
            IWRAP=NINT(360/ABS(DLON))
            IF(IM.LT.IWRAP) IWRAP=0
            JWRAP1=0
            JWRAP2=0
            KSCAN=0
            NSCAN=MOD(IGDTMPL(16)/32,2)
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  SET EXCEPTIONAL VALUES FOR POLAR STEREOGRAPHIC
         ELSEIF(IGDTNUM.EQ.20) THEN
            IM=IGDTMPL(8)
            JM=IGDTMPL(9)
            NSCAN=MOD(IGDTMPL(18)/32,2)
            IWRAP=0
            JWRAP1=0
            JWRAP2=0
            KSCAN=0
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  SET EXCEPTIONAL VALUES FOR LAMBERT CONFORMAL.
         ELSEIF(IGDTNUM.EQ.30) THEN
            IM=IGDTMPL(8)
            JM=IGDTMPL(9)
            NSCAN=MOD(IGDTMPL(18)/32,2)
            IWRAP=0
            JWRAP1=0
            JWRAP2=0
            KSCAN=0
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  SET EXCEPTIONAL VALUES FOR GAUSSIAN LAT/LON
         ELSEIF(IGDTNUM.EQ.40) THEN
            IM=IGDTMPL(8)
            JM=IGDTMPL(9)
            ISCALE=IGDTMPL(10)*IGDTMPL(11)
            IF(ISCALE==0) ISCALE=10**6
            RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
            RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
            ISCAN=MOD(IGDTMPL(19)/128,2)
            IF(ISCAN.EQ.0) THEN
               DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
            ELSE
               DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
            ENDIF
            IWRAP=NINT(360/ABS(DLON))
            IF(IM.LT.IWRAP) IWRAP=0
            JWRAP1=0
            JWRAP2=0
            IF(IWRAP.GT.0.AND.MOD(IWRAP,2).EQ.0) THEN
               JG=IGDTMPL(18)*2
               IF(JM.EQ.JG) THEN
                  JWRAP1=1
                  JWRAP2=2*JM+1
               ENDIF
            ENDIF
            NSCAN=MOD(IGDTMPL(19)/32,2)
            KSCAN=0
         ENDIF
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !  FILL NAVIGATION PARAMETER ARRAY
         IJKGDSA(1)=IM
         IJKGDSA(2)=JM
         IJKGDSA(3)=IWRAP
         IJKGDSA(4)=JWRAP1
         IJKGDSA(5)=JWRAP2
         IJKGDSA(6)=NSCAN
         IJKGDSA(7)=KSCAN
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       end associate
       class default
       error stop "descriptor class not recognized"
    end select
  end subroutine IJKGDS0


end module ijkgds0_mod
