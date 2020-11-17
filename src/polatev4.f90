module polatev4_mod
  use ip_grid_descriptor_mod
  use ijkgds0_mod, only: IJKGDS0
  implicit none

  private
  public :: polatev4

contains

  SUBROUTINE POLATEV4(IPOPT,input_desc, output_desc, &
       MI,MO,KM,IBI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATEV4   INTERPOLATE VECTOR FIELDS (SPECTRAL)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS SPECTRAL INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
    !           IT REQUIRES THAT THE INPUT FIELDS BE UNIFORMLY GLOBAL.
    !           OPTIONS ALLOW CHOICES BETWEEN TRIANGULAR SHAPE (IPOPT(1)=0)
    !           AND RHOMBOIDAL SHAPE (IPOPT(1)=1) WHICH HAS NO DEFAULT;
    !           A SECOND OPTION IS THE TRUNCATION (IPOPT(2)) WHICH DEFAULTS 
    !           TO A SENSIBLE TRUNCATION FOR THE INPUT GRID (IF OPT(2)=-1).
    !           NOTE THAT IF THE OUTPUT GRID IS NOT FOUND IN A SPECIAL LIST,
    !           THEN THE TRANSFORM BACK TO GRID IS NOT VERY FAST.
    !           THIS SPECIAL LIST CONTAINS GLOBAL CYLINDRICAL GRIDS,
    !           POLAR STEREOGRAPHIC GRIDS CENTERED AT THE POLE
    !           AND MERCATOR GRIDS. ONLY HORIZONTAL INTERPOLATION
    !           IS PERFORMED.
    !
    !           THE INPUT AND OUTPUT GRIDS ARE DEFINED BY THEIR GRIB 2 GRID
    !           DEFINITION TEMPLATE AS DECODED BY THE NCEP G2 LIBRARY. THE
    !           CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
    !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
    !           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
    !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
    !             (IGDTNUMO  =01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
    !                             NON-"E" STAGGERED
    !             (IGDTNUMO  =10) MERCATOR CYLINDRICAL
    !             (IGDTNUMO  =20) POLAR STEREOGRAPHIC AZIMUTHAL
    !             (IGDTNUMO  =30) LAMBERT CONFORMAL CONICAL
    !             (IGDTNUMI/O=40) GAUSSIAN CYLINDRICAL
    !
    !           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
    !           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
    !           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
    !           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
    !           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
    !
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
    !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
    !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
    !           IF IGDTNUMO<0, IN WHICH CASE THE NUMBER OF POINTS
    !           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT 
    !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !
    !           OUTPUT BITMAPS WILL ONLY BE CREATED WHEN THE OUTPUT GRID
    !           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
    !           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    ! 2001-06-18  IREDELL  IMPROVE DETECTION OF SPECIAL FAST TRANSFORM
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
    !                      ROUTINE GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATEV4(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,UI,VI, &
    !                         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1)=0 FOR TRIANGULAR, IPOPT(1)=1 FOR RHOMBOIDAL;
    !                IPOPT(2) IS TRUNCATION NUMBER
    !                (DEFAULTS TO SENSIBLE IF IPOPT(2)=-1).
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
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATEV FOR COMPLETE DEFINITION.
    !     IGDTLENI - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IGDTNUMO - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE. IGDTNUMO<0
    !                MEANS INTERPOLATE TO RANDOM STATION POINTS.
    !                OTHERWISE, SAME DEFINITION AS "IGDTNUMI".
    !     IGDTMPLO - INTEGER (IGDTLENO) GRID DEFINITION TEMPLATE ARRAY -
    !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATEV FOR COMPLETE DEFINITION.
    !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
    !     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
    !     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS (MUST BE ALL 0)
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
    !                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
    !                3    UNRECOGNIZED OUTPUT GRID
    !                41   INVALID NONGLOBAL INPUT GRID
    !                42   INVALID SPECTRAL METHOD PARAMETERS
    !
    ! SUBPROGRAMS CALLED:
    !   EARTH_RADIUS DETERMINE SIZE/SHAPE OF EARTH
    !   GDSWZD       GRID DESCRIPTION SECTION WIZARD
    !   SPTRUNV      SPECTRALLY TRUNCATE GRIDDED VECTOR FIELDS
    !   SPTRUNSV     SPECTRALLY INTERPOLATE VECTORS TO POLAR STEREO.
    !   SPTRUNMV     SPECTRALLY INTERPOLATE VECTORS TO MERCATOR
    !   SPTRUNGV     SPECTRALLY INTERPOLATE VECTORS TO STATIONS
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    !
    USE GDSWZD_MOD
    !
    IMPLICIT NONE
    !
    class(ip_grid_descriptor), intent(in) :: input_desc, output_desc
    INTEGER,          INTENT(IN   ) :: IPOPT(20), IBI(KM)
    INTEGER,          INTENT(IN   ) :: KM, MI, MO
    INTEGER,          INTENT(  OUT) :: IRET, IBO(KM), NO
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,             INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: CROT(MO),SROT(MO)
    !
    REAL,                 PARAMETER :: FILL=-9999.
    REAL,                 PARAMETER :: PI=3.14159265358979
    REAL,                 PARAMETER :: DPR=180./PI
    !
    INTEGER                         :: IDRTO, IROMB, ISKIPI, ISPEC
    INTEGER                         :: IDRTI, IMAXI, JMAXI, IM, JM
    INTEGER                         :: IPRIME, IG, IMO, JMO, IGO, JGO
    INTEGER                         :: ISCAN, JSCAN, NSCAN
    INTEGER                         :: ISCANO, JSCANO, NSCANO
    INTEGER                         :: ISCALE, IP, IPROJ, JSKIPI, JG
    INTEGER                         :: K, MAXWV, N, NI, NJ, NPS
    !
    REAL                            :: DLAT, DLON, DLATO, DLONO, DE, DR, DY
    REAL                            :: DUM, E2, H, HI, HJ
    REAL                            :: ORIENT, RERTH, SLAT
    REAL                            :: RLAT1, RLON1, RLAT2, RLON2, RLATI
    REAL                            :: UROT, VROT, UO2(MO,KM),VO2(MO,KM)
    REAL                            :: XMESH, X, XP, YP, XPTS(MO),YPTS(MO)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(output_desc%grid_number.GE.0) THEN
       CALL GDSWZD(output_desc, 0,MO,FILL,XPTS,YPTS, &
            RLON,RLAT,NO,CROT,SROT)
       IF(NO.EQ.0) IRET=3
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  AFFIRM APPROPRIATE INPUT GRID
    !    LAT/LON OR GAUSSIAN
    !    NO BITMAPS
    !    FULL ZONAL COVERAGE
    !    FULL MERIDIONAL COVERAGE
    
    select type(input_desc)
    type is(grib1_descriptor)
    type is(grib2_descriptor)
       select type(output_desc)
       type is(grib1_descriptor)
       type is(grib2_descriptor)
          associate(igdtnumi => input_desc%gdt_num, igdtmpli => input_desc%gdt_tmpl, &
               igdtnumo => output_desc%gdt_num, igdtmplo => output_desc%gdt_tmpl, &
               igdtleno => output_desc%gdt_len)
            IDRTI=IGDTNUMI
            IF(IDRTI==40) IDRTI=4
            IF(IDRTI==0.OR.IDRTI==4)THEN
               IM=IGDTMPLI(8)
               JM=IGDTMPLI(9)
               ISCALE=IGDTMPLI(10)*IGDTMPLI(11)
               IF(ISCALE==0) ISCALE=10**6
               RLON1=FLOAT(IGDTMPLI(13))/FLOAT(ISCALE)
               RLON2=FLOAT(IGDTMPLI(16))/FLOAT(ISCALE)
               ISCAN=MOD(IGDTMPLI(19)/128,2)
               JSCAN=MOD(IGDTMPLI(19)/64,2)
               NSCAN=MOD(IGDTMPLI(19)/32,2)
            ELSE
               IRET=41
            ENDIF
            DO K=1,KM
               IF(IBI(K).NE.0) IRET=41
            ENDDO
            IF(IRET.EQ.0) THEN
               IF(ISCAN.EQ.0) THEN
                  DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
               ELSE
                  DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
               ENDIF
               IG=NINT(360/ABS(DLON))
               IPRIME=1+MOD(-NINT(RLON1/DLON)+IG,IG)
               IMAXI=IG
               JMAXI=JM
               IF(MOD(IG,2).NE.0.OR.IM.LT.IG) IRET=41
            ENDIF
            IF(IRET.EQ.0.AND.IDRTI.EQ.0) THEN
               ISCALE=IGDTMPLI(10)*IGDTMPLI(11)
               IF(ISCALE==0) ISCALE=10**6
               RLAT1=FLOAT(IGDTMPLI(12))/FLOAT(ISCALE)
               RLAT2=FLOAT(IGDTMPLI(15))/FLOAT(ISCALE)
               DLAT=(RLAT2-RLAT1)/(JM-1)
               JG=NINT(180/ABS(DLAT))
               IF(JM.EQ.JG) IDRTI=256
               IF(JM.NE.JG.AND.JM.NE.JG+1) IRET=41
            ELSEIF(IRET.EQ.0.AND.IDRTI.EQ.4) THEN
               JG=IGDTMPLI(18)*2
               IF(JM.NE.JG) IRET=41
            ENDIF
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  SET PARAMETERS
            IF(IRET.EQ.0) THEN
               IROMB=IPOPT(1)
               MAXWV=IPOPT(2)
               IF(MAXWV.EQ.-1) THEN
                  IF(IROMB.EQ.0.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)
                  IF(IROMB.EQ.1.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)/2
                  IF(IROMB.EQ.0.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/2
                  IF(IROMB.EQ.1.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/4
                  IF(IROMB.EQ.0.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/2
                  IF(IROMB.EQ.1.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/4
               ENDIF
               IF((IROMB.NE.0.AND.IROMB.NE.1).OR.MAXWV.LT.0) IRET=42
            ENDIF
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  INTERPOLATE
            IF(IRET.EQ.0) THEN
               IF(NSCAN.EQ.0) THEN
                  ISKIPI=1
                  JSKIPI=IM
               ELSE
                  ISKIPI=JM
                  JSKIPI=1
               ENDIF
               IF(ISCAN.EQ.1) ISKIPI=-ISKIPI
               IF(JSCAN.EQ.0) JSKIPI=-JSKIPI
               ISPEC=0
               !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
               IF((IGDTNUMO.EQ.0.OR.IGDTNUMO.EQ.40).AND. &
                    MOD(IGDTMPLO(8),2).EQ.0.AND.IGDTMPLO(13).EQ.0.AND. &
                    IGDTMPLO(19).EQ.0) THEN
                  IDRTO=IGDTNUMO
                  IF(IDRTO==40)IDRTO=4
                  IMO=IGDTMPLO(8)
                  JMO=IGDTMPLO(9)
                  ISCALE=IGDTMPLO(10)*IGDTMPLO(11)
                  IF(ISCALE==0) ISCALE=10**6
                  RLON2=FLOAT(IGDTMPLO(16))/FLOAT(ISCALE)
                  DLONO=(MOD(RLON2-1+3600,360.)+1)/(IMO-1)
                  IGO=NINT(360/ABS(DLONO))
                  IF(IMO.EQ.IGO.AND.IDRTO.EQ.0) THEN
                     RLAT1=FLOAT(IGDTMPLO(12))/FLOAT(ISCALE)
                     RLAT2=FLOAT(IGDTMPLO(15))/FLOAT(ISCALE)
                     DLAT=(RLAT2-RLAT1)/(JMO-1)
                     JGO=NINT(180/ABS(DLAT))
                     IF(JMO.EQ.JGO) IDRTO=256
                     IF(JMO.EQ.JGO.OR.JMO.EQ.JGO+1) ISPEC=1
                  ELSEIF(IMO.EQ.IGO.AND.IDRTO.EQ.4) THEN
                     JGO=IGDTMPLO(18)*2
                     IF(JMO.EQ.JGO) ISPEC=1
                  ENDIF
                  IF(ISPEC.EQ.1) THEN
                     CALL SPTRUNV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMO,JMO, &
                          KM,IPRIME,ISKIPI,JSKIPI,MI,0,0,MO,0,UI,VI, &
                          .TRUE.,UO,VO,.FALSE.,DUM,DUM,.FALSE.,DUM,DUM)
                  ENDIF
                  !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
               ELSEIF(IGDTNUMO.EQ.20.AND. &
                    IGDTMPLO(8).EQ.IGDTMPLO(9).AND.MOD(IGDTMPLO(8),2).EQ.1.AND. &
                    IGDTMPLO(15).EQ.IGDTMPLO(16).AND.IGDTMPLO(18).EQ.64.AND. &
                    MOD(IGDTMPLO(12)/8,2).EQ.1) THEN
                  NPS=IGDTMPLO(8)
                  RLAT1=FLOAT(IGDTMPLO(10))*1.E-6
                  RLON1=FLOAT(IGDTMPLO(11))*1.E-6
                  ORIENT=FLOAT(IGDTMPLO(14))*1.E-6
                  XMESH=FLOAT(IGDTMPLO(15))*1.E-3
                  IPROJ=MOD(IGDTMPLO(17)/128,2)
                  IP=(NPS+1)/2
                  H=(-1.)**IPROJ
                  SLAT=FLOAT(ABS(IGDTMPLO(13)))*1.E-6
                  CALL EARTH_RADIUS(IGDTMPLO,IGDTLENO,RERTH,E2)
                  DE=(1.+SIN(SLAT/DPR))*RERTH
                  DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
                  XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/XMESH
                  YP=1+COS((RLON1-ORIENT)/DPR)*DR/XMESH
                  IF(NINT(XP).EQ.IP.AND.NINT(YP).EQ.IP) THEN
                     IF(IPROJ.EQ.0) THEN
                        CALL SPTRUNSV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                             IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                             SLAT,XMESH,ORIENT,UI,VI,.TRUE.,UO,VO,UO2,VO2, &
                             .FALSE.,DUM,DUM,DUM,DUM, &
                             .FALSE.,DUM,DUM,DUM,DUM)
                     ELSE
                        CALL SPTRUNSV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                             IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                             SLAT,XMESH,ORIENT,UI,VI,.TRUE.,UO2,VO2,UO,VO, &
                             .FALSE.,DUM,DUM,DUM,DUM, &
                             .FALSE.,DUM,DUM,DUM,DUM)
                     ENDIF
                     ISPEC=1
                  ENDIF
                  !  SPECIAL CASE OF MERCATOR GRID
               ELSEIF(IGDTNUMO.EQ.10) THEN
                  NI=IGDTMPLO(8)
                  NJ=IGDTMPLO(9)
                  RLAT1=FLOAT(IGDTMPLO(10))*1.0E-6
                  RLON1=FLOAT(IGDTMPLO(11))*1.0E-6
                  RLON2=FLOAT(IGDTMPLO(15))*1.0E-6
                  RLATI=FLOAT(IGDTMPLO(13))*1.0E-6
                  ISCANO=MOD(IGDTMPLO(16)/128,2)
                  JSCANO=MOD(IGDTMPLO(16)/64,2)
                  NSCANO=MOD(IGDTMPLO(16)/32,2)
                  DY=FLOAT(IGDTMPLO(19))*1.0E-3
                  HI=(-1.)**ISCANO
                  HJ=(-1.)**(1-JSCANO)
                  CALL EARTH_RADIUS(IGDTMPLO,IGDTLENO,RERTH,E2)
                  DLONO=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(NI-1)
                  DLATO=HJ*DY/(RERTH*COS(RLATI/DPR))*DPR
                  IF(NSCANO.EQ.0) THEN
                     CALL SPTRUNMV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NI,NJ, &
                          IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                          RLAT1,RLON1,DLATO,DLONO,UI,VI, &
                          .TRUE.,UO,VO,.FALSE.,DUM,DUM,.FALSE.,DUM,DUM)
                     ISPEC=1
                  ENDIF
               ENDIF
               !  GENERAL SLOW CASE
               IF(ISPEC.EQ.0) THEN
                  CALL SPTRUNGV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
                       IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON, &
                       UI,VI,.TRUE.,UO,VO,.FALSE.,X,X,.FALSE.,X,X)
                  DO K=1,KM
                     IBO(K)=0
                     DO N=1,NO
                        LO(N,K)=.TRUE.
                        UROT=CROT(N)*UO(N,K)-SROT(N)*VO(N,K)
                        VROT=SROT(N)*UO(N,K)+CROT(N)*VO(N,K)
                        UO(N,K)=UROT
                        VO(N,K)=VROT
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               DO K=1,KM
                  IBO(K)=1
                  DO N=1,NO
                     LO(N,K)=.FALSE.
                     UO(N,K)=0.
                     VO(N,K)=0.
                  ENDDO
               ENDDO
            ENDIF
          end associate
       class default
          print *, "unknown descriptor class"
          error stop
       end select
       class default
       print *, "unknown descriptor class"
       error stop
    end select

  END SUBROUTINE POLATEV4
end module polatev4_mod
