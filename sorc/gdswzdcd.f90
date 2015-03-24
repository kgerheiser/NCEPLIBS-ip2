 SUBROUTINE GDSWZDCD(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
                     XPTS,YPTS,RLON,RLAT,NRET, &
                     LROT,CROT,SROT,LMAP,XLON,XLAT,YLON,YLAT,AREA)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GDSWZDCD   GDS WIZARD FOR ROTATED EQUIDISTANT CYLINDRICAL
!   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2007-NOV-15
!
! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB GRID DESCRIPTION SECTION
!           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63)
!           AND RETURNS ONE OF THE FOLLOWING:
!             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
!             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
!           FOR NON-"E" STAGGERED ROTATED EQUIDISTANT CYLINDRICAL PROJECTIONS.
!           (MASS OR VELOCITY POINTS.)
!           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
!           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
!           OUTPUT ELEMENTS ARE SET TO FILL VALUES.
!           THE ACTUAL NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
!
! PROGRAM HISTORY LOG:
! 2010-JAN-15  GAYNO     BASED ON ROUTINES GDSWZDCB AND GDSWZDCA
!
! USAGE:    CALL GDSWZDCD(KGDS,IOPT,NPTS,FILL,XPTS,YPTS,RLON,RLAT,NRET,
!     &                   LROT,CROT,SROT,LMAP,XLON,XLAT,YLON,YLAT,AREA)
!
!   INPUT ARGUMENT LIST:
!     KGDS     - INTEGER (200) GDS PARAMETERS AS DECODED BY W3FI63
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
!     LROT     - INTEGER FLAG TO RETURN VECTOR ROTATIONS IF 1
!     LMAP     - INTEGER FLAG TO RETURN MAP JACOBIANS IF 1
!
!   OUTPUT ARGUMENT LIST:
!     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT<0
!     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT<0
!     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT>0
!     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT>0
!     NRET     - INTEGER NUMBER OF VALID POINTS COMPUTED
!     CROT     - REAL (NPTS) CLOCKWISE VECTOR ROTATION COSINES IF LROT=1
!     SROT     - REAL (NPTS) CLOCKWISE VECTOR ROTATION SINES IF LROT=1
!                (UGRID=CROT*UEARTH-SROT*VEARTH;
!                 VGRID=SROT*UEARTH+CROT*VEARTH)
!     XLON     - REAL (NPTS) DX/DLON IN 1/DEGREES IF LMAP=1
!     XLAT     - REAL (NPTS) DX/DLAT IN 1/DEGREES IF LMAP=1
!     YLON     - REAL (NPTS) DY/DLON IN 1/DEGREES IF LMAP=1
!     YLAT     - REAL (NPTS) DY/DLAT IN 1/DEGREES IF LMAP=1
!     AREA     - REAL (NPTS) AREA WEIGHTS IN M**2 IF LMAP=1
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
 IMPLICIT NONE
!
 INTEGER,                 PARAMETER     :: KD=SELECTED_REAL_KIND(15,45)
!
 INTEGER,                 INTENT(IN   ) :: IGDTNUM, IGDTLEN
 INTEGER(KIND=4),         INTENT(IN   ) :: IGDTMPL(IGDTLEN)
 INTEGER,                 INTENT(IN   ) :: IOPT,NPTS,LROT,LMAP
 INTEGER,                 INTENT(  OUT) :: NRET
!
 REAL,                    INTENT(IN   ) :: FILL
 REAL,                    INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
 REAL,                    INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
 REAL,                    INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
 REAL,                    INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
 REAL,                    INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
!
 REAL(KIND=KD),           PARAMETER     :: PI=3.14159265358979_KD
 REAL(KIND=KD),           PARAMETER     :: DPR=180._KD/PI
!
 INTEGER                                :: IROT,IM,JM,N,ISCALE
 INTEGER                                :: I_OFFSET_ODD,I_OFFSET_EVEN,J_OFFSET
!
 REAL(KIND=KD)                          :: HS,HS2,RERTH
 REAL(KIND=KD)                          :: RLAT1,RLON1,RLAT0,RLON0,RLAT2,RLON2
 REAL(KIND=KD)                          :: SLAT1,CLAT1,SLAT0,CLAT0
 REAL(KIND=KD)                          :: SLON,SLAT2,CLAT2,CLON2
 REAL(KIND=KD)                          :: CLON1,SLATR,CLATR,CLONR
 REAL(KIND=KD)                          :: RLATR,RLONR,DLATS,DLONS
 REAL                                   :: DUM1,DUM2
 REAL(KIND=KD)                          :: SLAT,CLAT,CLON,WBD,SBD,NBD,EBD,TERM1,TERM2
 REAL                                   :: XMIN,XMAX,YMIN,YMAX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! IS THIS A ROTATED LAT/LON GRID?
 IF(IGDTNUM/=1)THEN
   CALL GDSWZDCD_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
   RETURN
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! IS THE EARTH RADIUS DEFINED?
 CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,DUM1,DUM2)
 RERTH=DUM1
 IF(RERTH<0.)THEN
   CALL GDSWZDCD_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
   RETURN
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  IS THIS AN "E"-STAGGER GRID?  ROUTINE CAN'T PROCESS THOSE.
 I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
 I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
 J_OFFSET=MOD(IGDTMPL(19)/2,2)
 IF(I_OFFSET_ODD/=I_OFFSET_EVEN) THEN
   CALL GDSWZDCD_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
   RETURN
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ISCALE=IGDTMPL(10)*IGDTMPL(11)
 IF(ISCALE==0) ISCALE=1E6
 RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
 RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
 RLAT0=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
 RLAT0=RLAT0+90.0_KD
 RLON0=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)
 RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
 RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
 IROT=MOD(IGDTMPL(14)/8,2)
 IM=IGDTMPL(8)
 JM=IGDTMPL(9)
 SLAT1=SIN(RLAT1/DPR)
 CLAT1=COS(RLAT1/DPR)
 SLAT0=SIN(RLAT0/DPR)
 CLAT0=COS(RLAT0/DPR)
 HS=SIGN(1._KD,MOD(RLON1-RLON0+180+3600,360._KD)-180)
 CLON1=COS((RLON1-RLON0)/DPR)
 SLATR=CLAT0*SLAT1-SLAT0*CLAT1*CLON1
 CLATR=SQRT(1-SLATR**2)
 CLONR=(CLAT0*CLAT1*CLON1+SLAT0*SLAT1)/CLATR
 RLATR=DPR*ASIN(SLATR)
 RLONR=HS*DPR*ACOS(CLONR)
 WBD=RLONR
 SBD=RLATR
 SLAT2=SIN(RLAT2/DPR)
 CLAT2=COS(RLAT2/DPR)
 HS2=SIGN(1._KD,MOD(RLON2-RLON0+180+3600,360._KD)-180)
 CLON2=COS((RLON2-RLON0)/DPR)
 SLATR=CLAT0*SLAT2-SLAT0*CLAT2*CLON2
 CLATR=SQRT(1-SLATR**2)
 CLONR=(CLAT0*CLAT2*CLON2+SLAT0*SLAT2)/CLATR
 NBD=DPR*ASIN(SLATR)
 EBD=HS2*DPR*ACOS(CLONR)
 DLATS=(NBD-SBD)/FLOAT(JM-1)
 DLONS=(EBD-WBD)/FLOAT(IM-1)
 IF(I_OFFSET_ODD==1) WBD=WBD+(0.5_KD*DLONS)
 IF(J_OFFSET==1) SBD=SBD+(0.5_KD*DLATS)
 XMIN=0
 XMAX=IM+1
 YMIN=0
 YMAX=JM+1
 NRET=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
 IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
   DO N=1,NPTS
     IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
        YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
       RLONR=WBD+(XPTS(N)-1._KD)*DLONS
       RLATR=SBD+(YPTS(N)-1._KD)*DLATS
       IF(RLONR <= 0._KD) THEN
         HS=-1.0_KD
       ELSE
         HS=1.0_KD
       ENDIF
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
       IF(LROT.EQ.1) THEN
         IF(IROT.EQ.1) THEN
           IF(CLATR.LE.0) THEN
             CROT(N)=-SIGN(1._KD,SLATR*SLAT0)
             SROT(N)=0
           ELSE
             SLON=SIN((RLON(N)-RLON0)/DPR)
             CROT(N)=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
             SROT(N)=SLAT0*SLON/CLATR
           ENDIF
         ELSE
           CROT(N)=1
           SROT(N)=0
         ENDIF
       ENDIF
       IF(LMAP.EQ.1) THEN
         IF(CLATR.LE.0) THEN
           XLON(N)=FILL
           XLAT(N)=FILL
           YLON(N)=FILL
           YLAT(N)=FILL
           AREA(N)=FILL
         ELSE
           SLON=SIN((RLON(N)-RLON0)/DPR)
           TERM1=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
           TERM2=SLAT0*SLON/CLATR
           XLON(N)=TERM1*CLAT/(DLONS*CLATR)
           XLAT(N)=-TERM2/(DLONS*CLATR)
           YLON(N)=TERM2*CLAT/DLATS
           YLAT(N)=TERM1/DLATS
           AREA(N)=2._KD*(RERTH**2)*CLATR*(DLONS/DPR)* &
                         SIN(0.5_KD*DLATS/DPR)
         ENDIF
       ENDIF
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
       XPTS(N)=(RLONR-WBD)/DLONS+1._KD
       YPTS(N)=(RLATR-SBD)/DLATS+1._KD
       IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
          YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
         NRET=NRET+1
         IF(LROT.EQ.1) THEN
           IF(IROT.EQ.1) THEN
             IF(CLATR.LE.0) THEN
               CROT(N)=-SIGN(1._KD,SLATR*SLAT0)
               SROT(N)=0
             ELSE
               SLON=SIN((RLON(N)-RLON0)/DPR)
               CROT(N)=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
               SROT(N)=SLAT0*SLON/CLATR
             ENDIF
           ELSE
             CROT(N)=1
             SROT(N)=0
           ENDIF
         ENDIF
         IF(LMAP.EQ.1)THEN
           IF(CLATR.LE.0) THEN
             XLON(N)=FILL
             XLAT(N)=FILL
             YLON(N)=FILL
             YLAT(N)=FILL
             AREA(N)=FILL
           ELSE
             SLON=SIN((RLON(N)-RLON0)/DPR)
             TERM1=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
             TERM2=SLAT0*SLON/CLATR
             XLON(N)=TERM1*CLAT/(DLONS*CLATR)
             XLAT(N)=-TERM2/(DLONS*CLATR)
             YLON(N)=TERM2*CLAT/DLATS
             YLAT(N)=TERM1/DLATS
             AREA(N)=2._KD*(RERTH**2)*CLATR*(DLONS/DPR)* &
                           SIN(0.5_KD*DLATS/DPR)
           ENDIF
         ENDIF
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
 END SUBROUTINE GDSWZDCD
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 SUBROUTINE GDSWZDCD_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)

 IMPLICIT NONE

 INTEGER, INTENT(IN   ) :: IOPT, NPTS

 REAL,    INTENT(IN   ) :: FILL
 REAL,    INTENT(  OUT) :: RLAT(NPTS),RLON(NPTS)
 REAL,    INTENT(  OUT) :: XPTS(NPTS),YPTS(NPTS)

 IF(IOPT>=0) THEN
   RLON=FILL
   RLAT=FILL
 ENDIF
 IF(IOPT<=0) THEN
   XPTS=FILL
   YPTS=FILL
 ENDIF

 END SUBROUTINE GDSWZDCD_ERROR
