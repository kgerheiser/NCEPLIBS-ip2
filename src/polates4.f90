module polates4_mod
  use earth_radius_mod
  implicit none

  private
  public :: polates4

contains

  SUBROUTINE POLATES4(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       IGDTNUMO,IGDTMPLO,IGDTLENO, &
       MI,MO,KM,IBI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATES4   INTERPOLATE SCALAR FIELDS (SPECTRAL)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS SPECTRAL INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR SCALAR FIELDS.
    !           IT REQUIRES THAT THE INPUT FIELDS BE UNIFORMLY GLOBAL.
    !           OPTIONS ALLOW CHOICES BETWEEN TRIANGULAR SHAPE (IPOPT(1)=0)
    !           AND RHOMBOIDAL SHAPE (IPOPT(1)=1) WHICH HAS NO DEFAULT;
    !           A SECOND OPTION IS THE TRUNCATION (IPOPT(2)) WHICH DEFAULTS 
    !           TO A SENSIBLE TRUNCATION FOR THE INPUT GRID (IF OPT(2)=-1).
    !           NOTE THAT IF THE OUTPUT GRID IS NOT FOUND IN A SPECIAL LIST,
    !           THEN THE TRANSFORM BACK TO GRID IS NOT VERY FAST.
    !           THIS SPECIAL LIST CONTAINS GLOBL CYLINDRICAL GRIDS,
    !           POLAR STEREOGRAPHIC GRIDS CENTERED AT THE POLE
    !           AND MERCATOR GRIDS. ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
    !
    !           THE CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
    !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
    !           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
    !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
    !             (IGDTNUMO  =01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
    !                             NON-"E" STAGGERED
    !             (IGDTNUMO  =10) MERCATOR CYLINDRICAL
    !             (IGDTNUMO  =20) POLAR STEREOGRAPHIC AZIMUTHAL
    !             (IGDTNUMO  =30) LAMBERT CONFORMAL CONICAL
    !             (IGDTNUMI/O=40) GAUSSIAN CYLINDRICAL
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
    !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
    !           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
    !           IF IGDTNUMO<0, IN WHICH CASE THE NUMBER OF POINTS
    !           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
    !           OUTPUT BITMAPS WILL NOT BE CREATED.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    ! 2001-06-18  IREDELL  IMPROVE DETECTION OF SPECIAL FAST TRANSFORM
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
    !                      VERSION OF GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATES4(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,GI, &
    !                         NO,RLAT,RLON,IBO,LO,GO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1)=0 FOR TRIANGULAR, IPOPT(1)=1 FOR RHOMBOIDAL;
    !                IPOPT(2) IS TRUNCATION NUMBER
    !                (DEFAULTS TO SENSIBLE IF IPOPT(2)=-1).
    !     IGDTNUMI - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                  00 - EQUIDISTANT CYLINDRICAL
    !                  01 - ROTATED EQUIDISTANT CYLINDRICAL.  "E"
    !                       AND NON-"E" STAGGERED
    !                  10 - MERCATOR CYCLINDRICAL
    !                  20 - POLAR STEREOGRAPHIC AZIMUTHAL
    !                  30 - LAMBERT CONFORMAL CONICAL
    !                  40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
    !     IGDTMPLI - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY -
    !                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE:
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATES FOR COMPLETE DEFINITION.
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
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATES FOR COMPLETE DEFINITION.
    !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
    !     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
    !     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS (MUST BE ALL 0)
    !     GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO<0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO<0)
    !
    !   OUTPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO>=0)
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
    !     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
    !     GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
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
    !   SPTRUN       SPECTRALLY TRUNCATE GRIDDED SCALAR FIELDS
    !   SPTRUNS      SPECTRALLY INTERPOLATE SCALARS TO POLAR STEREO.
    !   SPTRUNM      SPECTRALLY INTERPOLATE SCALARS TO MERCATOR
    !   SPTRUNG      SPECTRALLY INTERPOLATE SCALARS TO STATIONS
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    USE GDSWZD_MOD_ip2
    !
    IMPLICIT NONE
    !
    INTEGER,          INTENT(IN   ) :: IGDTNUMI, IGDTLENI
    INTEGER,          INTENT(IN   ) :: IGDTMPLI(IGDTLENI)
    INTEGER,          INTENT(IN   ) :: IGDTNUMO, IGDTLENO
    INTEGER,          INTENT(IN   ) :: IGDTMPLO(IGDTLENO)
    INTEGER,          INTENT(IN   ) :: IPOPT(20)
    INTEGER,          INTENT(IN   ) :: MI, MO
    INTEGER,          INTENT(IN   ) :: IBI(KM), KM
    INTEGER,          INTENT(  OUT) :: IBO(KM), IRET, NO
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: GI(MI,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: GO(MO,KM)
    !
    REAL,             PARAMETER     :: FILL=-9999.
    REAL,             PARAMETER     :: PI=3.14159265358979
    REAL,             PARAMETER     :: DPR=180./PI
    !
    INTEGER                         :: IDRTI, IDRTO, IG, JG, IM, JM
    INTEGER                         :: IGO, JGO, IMO, JMO
    INTEGER                         :: ISCAN, JSCAN, NSCAN
    INTEGER                         :: ISCANO, JSCANO, NSCANO
    INTEGER                         :: ISKIPI, JSKIPI, ISCALE
    INTEGER                         :: IMAXI, JMAXI, ISPEC
    INTEGER                         :: IP, IPRIME, IPROJ, IROMB, K
    INTEGER                         :: MAXWV, N, NI, NJ, NPS
    !
    REAL                            :: DE, DR, DY
    REAL                            :: DLAT, DLON, DLATO, DLONO
    REAL                            :: GO2(MO,KM), H, HI, HJ
    REAL                            :: ORIENT, SLAT, RERTH, E2
    REAL                            :: RLAT1, RLON1, RLAT2, RLON2, RLATI
    REAL                            :: XMESH, XP, YP
    REAL                            :: XPTS(MO), YPTS(MO)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(IGDTNUMO.GE.0) THEN
       CALL GDSWZD(IGDTNUMO,IGDTMPLO,IGDTLENO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) IRET=3
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  AFFIRM APPROPRIATE INPUT GRID
    !    LAT/LON OR GAUSSIAN
    !    NO BITMAPS
    !    FULL ZONAL COVERAGE
    !    FULL MERIDIONAL COVERAGE
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
            MOD(IGDTMPLO(8),2).EQ.0.AND.IGDTMPLO(13).EQ.0.AND.IGDTMPLO(19).EQ.0) THEN
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
             CALL SPTRUN(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMO,JMO, &
                  KM,IPRIME,ISKIPI,JSKIPI,MI,0,0,MO,0,GI,GO)
          ENDIF
          !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
       ELSEIF(IGDTNUMO.EQ.20.AND. &
            IGDTMPLO(8).EQ.IGDTMPLO(9).AND.MOD(IGDTMPLO(8),2).EQ.1.AND. &
            IGDTMPLO(15).EQ.IGDTMPLO(16).AND.IGDTMPLO(18).EQ.64) THEN
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
                CALL SPTRUNS(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     SLAT,XMESH,ORIENT,GI,GO,GO2)
             ELSE
                CALL SPTRUNS(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     SLAT,XMESH,ORIENT,GI,GO2,GO)
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
             CALL SPTRUNM(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NI,NJ, &
                  IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                  RLAT1,RLON1,DLATO,DLONO,GI,GO)
             ISPEC=1
          ENDIF
       ENDIF
       !  GENERAL SLOW CASE
       IF(ISPEC.EQ.0) THEN
          CALL SPTRUNG(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
               IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON,GI,GO)
       ENDIF
       DO K=1,KM
          IBO(K)=0
          DO N=1,NO
             LO(N,K)=.TRUE.
          ENDDO
       ENDDO
    ELSE
       DO K=1,KM
          IBO(K)=1
          DO N=1,NO
             LO(N,K)=.FALSE.
             GO(N,K)=0.
          ENDDO
       ENDDO
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE POLATES4

end module polates4_mod
