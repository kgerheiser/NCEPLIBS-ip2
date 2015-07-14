 SUBROUTINE POLATES3(IPOPT,GDTNUMI,GDTMPLI,GDTLENI, &
                     GDTNUMO,GDTMPLO,GDTLENO, &
                     MI,MO,KM,IBI,LI,GI, &
                     NO,RLAT,RLON,IBO,LO,GO,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  POLATES3   INTERPOLATE SCALAR FIELDS (BUDGET)
!   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
!
! ABSTRACT: THIS SUBPROGRAM PERFORMS BUDGET INTERPOLATION
!           FROM ANY GRID TO ANY GRID (OR TO RANDOM STATION
!           POINTS) FOR SCALAR FIELDS.  THE ALGORITHM
!           SIMPLY COMPUTES (WEIGHTED) AVERAGES OF
!           BILINEARLY INTERPOLATED POINTS ARRANGED IN A SQUARE BOX
!           CENTERED AROUND EACH OUTPUT GRID POINT AND STRETCHING
!           NEARLY HALFWAY TO EACH OF THE NEIGHBORING GRID POINTS.
!           OPTIONS ALLOW CHOICES OF NUMBER OF POINTS IN EACH RADIUS
!           FROM THE CENTER POINT (IPOPT(1)) WHICH DEFAULTS TO 2
!           (IF IPOPT(1)=-1) MEANING THAT 25 POINTS WILL BE AVERAGED;
!           FURTHER OPTIONS ARE THE RESPECTIVE WEIGHTS FOR THE RADIUS
!           POINTS STARTING AT THE CENTER POINT (IPOPT(2:2+IPOPT(1))
!           WHICH DEFAULTS TO ALL 1 (IF IPOPT(1)=-1 OR IPOPT(2)=-1).
!           A SPECIAL INTERPOLATION IS DONE IF IPOPT(2)=-2.
!           IN THIS CASE, THE BOXES STRETCH NEARLY ALL THE WAY TO
!           EACH OF THE NEIGHBORING GRID POINTS AND THE WEIGHTS
!           ARE THE ADJOINT OF THE BILINEAR INTERPOLATION WEIGHTS.
!           THIS CASE GIVES QUASI-SECOND-ORDER BUDGET INTERPOLATION.
!           ANOTHER OPTION IS THE MINIMUM PERCENTAGE FOR MASK,
!           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
!           (IPOPT(3+IPOPT(1)) WHICH DEFAULTS TO 50 (IF -1).
!           IN CASES WHERE THERE IS NO OR INSUFFICIENT VALID INPUT DATA,
!           THE USER MAY CHOOSE TO SEARCH FOR THE NEAREST VALID DATA. 
!           THIS IS INVOKED BY SETTING IPOPT(20) TO THE WIDTH OF 
!           THE SEARCH SQUARE. THE DEFAULT IS 1 (NO SEARCH).  SQUARES ARE
!           SEARCHED FOR VALID DATA IN A SPIRAL PATTERN
!           STARTING FROM THE CENTER.  NO SEARCHING IS DONE WHERE
!           THE OUTPUT GRID IS OUTSIDE THE INPUT GRID.
!           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
!
!           THE CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
!           "GDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
!           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
!             (GDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
!             (GDTNUMI/O=01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
!                            NON-"E" STAGGERED
!             (GDTNUMI/O=10) MERCATOR CYLINDRICAL
!             (GDTNUMI/O=20) POLAR STEREOGRAPHIC AZIMUTHAL
!             (GDTNUMI/O=30) LAMBERT CONFORMAL CONICAL
!             (GDTNUMI/O=40) GAUSSIAN CYLINDRICAL
!           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
!           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
!           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
!           IF GDTNUMO=GDTNUMO-255, IN WHICH CASE THE NUMBER OF POINTS
!           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
!
!           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
!           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
!           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
!           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
!        
! PROGRAM HISTORY LOG:
!   96-04-10  IREDELL
! 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
! 1999-04-08  IREDELL  ADDED BILINEAR OPTION IPOPT(2)=-2
! 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
! 2006-01-04  GAYNO    ADDED OPTION TO DO SUBSECTION OF OUTPUT GRID.
!                      ADDED SPIRAL SEARCH OPTION.
! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
!                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
!
! USAGE:    CALL POLATES3(IPOPT,GDTNUMI,GDTMPLI,GDTLENI, &
!                         GDTNUMO,GDTMPLO,GDTLENO, &
!                         MI,MO,KM,IBI,LI,GI, &
!                         NO,RLAT,RLON,IBO,LO,GO,IRET)
!
!   INPUT ARGUMENT LIST:
!     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
!                IPOPT(1) IS NUMBER OF RADIUS POINTS
!                (DEFAULTS TO 2 IF IPOPT(1)=-1);
!                IPOPT(2:2+IPOPT(1)) ARE RESPECTIVE WEIGHTS
!                (DEFAULTS TO ALL 1 IF IPOPT(1)=-1 OR IPOPT(2)=-1).
!                IPOPT(3+IPOPT(1)) IS MINIMUM PERCENTAGE FOR MASK
!                (DEFAULTS TO 50 IF IPOPT(3+IPOPT(1)=-1)
!     GDTNUMI  - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
!                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
!                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
!     GDTMPLI  - INTEGER (GDTLENI) GRID DEFINITION TEMPLATE ARRAY -
!                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
!                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
!     GDTLENI  - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
!                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
!                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
!     GDTNUMO  - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
!                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
!                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE. GDTNUMO=GDTNUM-255
!                MEANS INTERPOLATE TO RANDOM STATION POINTS.
!     GDTMPLO  - INTEGER (GDTLENO) GRID DEFINITION TEMPLATE ARRAY -
!                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
!                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
!     GDTLENO  - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
!                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
!                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
!     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
!                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
!     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
!                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
!     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
!     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
!     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
!     GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
!     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF GDTNUMO=GDTNUMO-255)
!     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF GDTNUMO=GDTNUMO-255)
!
!   OUTPUT ARGUMENT LIST:
!     NO       - INTEGER NUMBER OF OUTPUT POINTS 
!     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (ONLY IF GDTNUMO>=0)
!     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (ONLY IF GDTNUMO>=0)
!     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
!     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
!     GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
!     IRET     - INTEGER RETURN CODE
!                0    SUCCESSFUL INTERPOLATION
!                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
!                3    UNRECOGNIZED OUTPUT GRID
!                32   INVALID BUDGET METHOD PARAMETERS
!
! SUBPROGRAMS CALLED:
!   GDSWZD       GRID DESCRIPTION SECTION WIZARD
!   IJKGDS0      SET UP PARAMETERS FOR IJKGDS1
!   IJKGDS1      RETURN FIELD POSITION FOR A GIVEN GRID POINT
!   POLFIXS      MAKE MULTIPLE POLE SCALAR VALUES CONSISTENT
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
 IMPLICIT NONE
!
 INTEGER,        INTENT(IN   ) :: GDTNUMI, GDTLENI
 INTEGER(KIND=4),INTENT(IN   ) :: GDTMPLI(GDTLENI)
 INTEGER,        INTENT(IN   ) :: GDTNUMO, GDTLENO
 INTEGER(KIND=4),INTENT(IN   ) :: GDTMPLO(GDTLENO)
 INTEGER,    INTENT(IN   )     :: IBI(KM), IPOPT(20)
 INTEGER,    INTENT(IN   )     :: KM, MI, MO
 INTEGER,    INTENT(  OUT)     :: IBO(KM), IRET, NO
!
 LOGICAL*1,  INTENT(IN   )     :: LI(MI,KM)
 LOGICAL*1,  INTENT(  OUT)     :: LO(MO,KM)
!
 REAL,       INTENT(IN   )     :: GI(MI,KM)
 REAL,       INTENT(INOUT)     :: RLAT(MO), RLON(MO)
 REAL,       INTENT(  OUT)     :: GO(MO,KM)
!
 REAL,       PARAMETER         :: FILL=-9999.
!
 INTEGER                       :: IJKGDS1, I1, J1, I2, J2, IB, JB
 INTEGER                       :: IJKGDSA(20), IX, JX, IXS, JXS
 INTEGER                       :: K, KXS, KXT, GDTNUMO2
 INTEGER                       :: LB, LSW, MP, MSPIRAL, MX
 INTEGER                       :: N, NB, NB1, NB2, NB3, NB4, NV, NX
 INTEGER                       :: N11(MO),N21(MO),N12(MO),N22(MO)
!
 REAL,    ALLOCATABLE          :: DUM1(:),DUM2(:),DUM3(:),DUM4(:)
 REAL,    ALLOCATABLE          :: DUM5(:),DUM6(:),DUM7(:)
 REAL                          :: GB, LAT(1), LON(1)
 REAL                          :: PMP, RB2, RLOB(MO), RLAB(MO), WB
 REAL                          :: W11(MO), W21(MO), W12(MO), W22(MO)
 REAL                          :: WO(MO,KM), XF, YF, XI, YI, XX, YY
 REAL                          :: XPTS(MO),YPTS(MO),XPTB(MO),YPTB(MO)
 REAL                          :: XXX(1), YYY(1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
!  DO SUBSECTION OF GRID IF KGDSO(1) IS SUBTRACTED FROM 255.
 IRET=0
 IF(GDTNUMO.GE.0) THEN
   ALLOCATE(DUM1(MO))
   ALLOCATE(DUM2(MO))
   ALLOCATE(DUM3(MO))
   ALLOCATE(DUM4(MO))
   ALLOCATE(DUM5(MO))
   ALLOCATE(DUM6(MO))
   ALLOCATE(DUM7(MO))
   CALL GDSWZD(GDTNUMO,GDTMPLO,GDTLENO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT, &
               NO,0,DUM1,DUM2,0,DUM3,DUM4,DUM5,DUM6,DUM7)
   DEALLOCATE(DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7)
   IF(NO.EQ.0) IRET=3
 ELSE
   GDTNUMO2=255+GDTNUMO
   ALLOCATE(DUM1(MO))
   ALLOCATE(DUM2(MO))
   ALLOCATE(DUM3(MO))
   ALLOCATE(DUM4(MO))
   ALLOCATE(DUM5(MO))
   ALLOCATE(DUM6(MO))
   ALLOCATE(DUM7(MO))
   CALL GDSWZD(GDTNUMO2,GDTMPLO,GDTLENO,-1,MO,FILL,XPTS,YPTS,RLON,RLAT, &
               NO,0,DUM1,DUM2,0,DUM3,DUM4,DUM5,DUM6,DUM7)
   DEALLOCATE(DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7)
   IF(NO.EQ.0) IRET=3
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SET PARAMETERS
 IF(IPOPT(1).GT.16) IRET=32  
 MSPIRAL=MAX(IPOPT(20),1)
 NB1=IPOPT(1)
 IF(NB1.EQ.-1) NB1=2
 IF(IRET.EQ.0.AND.NB1.LT.0) IRET=32
 LSW=1
 IF(IPOPT(2).EQ.-2) LSW=2
 IF(IPOPT(1).EQ.-1.OR.IPOPT(2).EQ.-1) LSW=0
 IF(IRET.EQ.0.AND.LSW.EQ.1.AND.NB1.GT.15) IRET=32
 MP=IPOPT(3+IPOPT(1))
 IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
 IF(MP.LT.0.OR.MP.GT.100) IRET=32
 PMP=MP*0.01
 IF(IRET.EQ.0) THEN
   NB2=2*NB1+1
   RB2=1./NB2
   NB3=NB2*NB2
   NB4=NB3
   IF(LSW.EQ.2) THEN
     RB2=1./(NB1+1)
     NB4=(NB1+1)**4
   ELSEIF(LSW.EQ.1) THEN
     NB4=IPOPT(2)
     DO IB=1,NB1
       NB4=NB4+8*IB*IPOPT(2+IB)
     ENDDO
   ENDIF
 ELSE
   NB3=0
   NB4=1
 ENDIF
 DO K=1,KM
   DO N=1,NO
     GO(N,K)=0.
     WO(N,K)=0.
   ENDDO
 ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX
 CALL IJKGDS0(GDTNUMI,GDTMPLI,GDTLENI,IJKGDSA)
 DO NB=1,NB3
!  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS
   JB=(NB-1)/NB2-NB1
   IB=NB-(JB+NB1)*NB2-NB1-1
   LB=MAX(ABS(IB),ABS(JB))
   WB=1
   IF(LSW.EQ.2) THEN
     WB=(NB1+1-ABS(IB))*(NB1+1-ABS(JB))
   ELSEIF(LSW.EQ.1) THEN
     WB=IPOPT(2+LB)
   ENDIF
   IF(WB.NE.0) THEN
     DO N=1,NO
       XPTB(N)=XPTS(N)+IB*RB2
       YPTB(N)=YPTS(N)+JB*RB2
     ENDDO
     ALLOCATE(DUM1(NO))
     ALLOCATE(DUM2(NO))
     ALLOCATE(DUM3(NO))
     ALLOCATE(DUM4(NO))
     ALLOCATE(DUM5(NO))
     ALLOCATE(DUM6(NO))
     ALLOCATE(DUM7(NO))
     CALL GDSWZD(GDTNUMO,GDTMPLO,GDTLENO, 1,NO,FILL,XPTB,YPTB,RLOB,RLAB, &
                 NV,0,DUM1,DUM2,0,DUM3,DUM4,DUM5,DUM6,DUM7)
     CALL GDSWZD(GDTNUMI,GDTMPLI,GDTLENI,-1,NO,FILL,XPTB,YPTB,RLOB,RLAB, &
                 NV,0,DUM1,DUM2,0,DUM3,DUM4,DUM5,DUM6,DUM7)
     DEALLOCATE(DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7)
     IF(IRET.EQ.0.AND.NV.EQ.0.AND.LB.EQ.0) IRET=2
     DO N=1,NO
       XI=XPTB(N)
       YI=YPTB(N)
       IF(XI.NE.FILL.AND.YI.NE.FILL) THEN
         I1=XI
         I2=I1+1
         J1=YI
         J2=J1+1
         XF=XI-I1
         YF=YI-J1
         N11(N)=IJKGDS1(I1,J1,IJKGDSA)
         N21(N)=IJKGDS1(I2,J1,IJKGDSA)
         N12(N)=IJKGDS1(I1,J2,IJKGDSA)
         N22(N)=IJKGDS1(I2,J2,IJKGDSA)
         IF(MIN(N11(N),N21(N),N12(N),N22(N)).GT.0) THEN
           W11(N)=(1-XF)*(1-YF)
           W21(N)=XF*(1-YF)
           W12(N)=(1-XF)*YF
           W22(N)=XF*YF
         ELSE
           N11(N)=0
           N21(N)=0
           N12(N)=0
           N22(N)=0
         ENDIF
       ELSE
         N11(N)=0
         N21(N)=0
         N12(N)=0
         N22(N)=0
       ENDIF
     ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  INTERPOLATE WITH OR WITHOUT BITMAPS
     DO K=1,KM
       DO N=1,NO
         IF(N11(N).GT.0) THEN
           IF(IBI(K).EQ.0) THEN
             GB=W11(N)*GI(N11(N),K)+W21(N)*GI(N21(N),K) &
                +W12(N)*GI(N12(N),K)+W22(N)*GI(N22(N),K)
             GO(N,K)=GO(N,K)+WB*GB
             WO(N,K)=WO(N,K)+WB
           ELSE
             IF(LI(N11(N),K)) THEN
               GO(N,K)=GO(N,K)+WB*W11(N)*GI(N11(N),K)
               WO(N,K)=WO(N,K)+WB*W11(N)
             ENDIF
             IF(LI(N21(N),K)) THEN
               GO(N,K)=GO(N,K)+WB*W21(N)*GI(N21(N),K)
               WO(N,K)=WO(N,K)+WB*W21(N)
             ENDIF
             IF(LI(N12(N),K)) THEN
               GO(N,K)=GO(N,K)+WB*W12(N)*GI(N12(N),K)
               WO(N,K)=WO(N,K)+WB*W12(N)
             ENDIF
             IF(LI(N22(N),K)) THEN
               GO(N,K)=GO(N,K)+WB*W22(N)*GI(N22(N),K)
               WO(N,K)=WO(N,K)+WB*W22(N)
             ENDIF
           ENDIF
         ENDIF
       ENDDO
     ENDDO
   ENDIF
 ENDDO   ! sub-grid points
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE OUTPUT BITMAPS AND FIELDS
 KM_LOOP : DO K=1,KM
   IBO(K)=IBI(K)
   N_LOOP : DO N=1,NO
     LO(N,K)=WO(N,K).GE.PMP*NB4
     IF(LO(N,K)) THEN
       GO(N,K)=GO(N,K)/WO(N,K)
     ELSEIF (MSPIRAL.GT.1) THEN
       ALLOCATE(DUM1(1))
       ALLOCATE(DUM2(1))
       ALLOCATE(DUM3(1))
       ALLOCATE(DUM4(1))
       ALLOCATE(DUM5(1))
       ALLOCATE(DUM6(1))
       ALLOCATE(DUM7(1))
       LAT(1)=RLAT(N)
       LON(1)=RLON(N)
       CALL GDSWZD(GDTNUMI,GDTMPLI,GDTLENI,-1,1,FILL,XXX,YYY,LON,LAT, &
                   NV,0,DUM1,DUM2,0,DUM3,DUM4,DUM5,DUM6,DUM7)
       DEALLOCATE(DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7)
       XX=XXX(1)
       YY=YYY(1)
       IF(NV.EQ.1)THEN
         I1=NINT(XX)
         J1=NINT(YY)
         IXS=SIGN(1.,XX-I1)
         JXS=SIGN(1.,YY-J1)
         SPIRAL_LOOP : DO MX=2,MSPIRAL**2
           KXS=SQRT(4*MX-2.5)
           KXT=MX-(KXS**2/4+1)
           SELECT CASE(MOD(KXS,4))
           CASE(1)
             IX=I1-IXS*(KXS/4-KXT)
             JX=J1-JXS*KXS/4
           CASE(2)
             IX=I1+IXS*(1+KXS/4)
             JX=J1-JXS*(KXS/4-KXT)
           CASE(3)
             IX=I1+IXS*(1+KXS/4-KXT)
             JX=J1+JXS*(1+KXS/4)
           CASE DEFAULT
             IX=I1-IXS*KXS/4
             JX=J1+JXS*(KXS/4-KXT)
           END SELECT
           NX=IJKGDS1(IX,JX,IJKGDSA)
           IF(NX.GT.0.)THEN
             IF(LI(NX,K).OR.IBI(K).EQ.0) THEN
               GO(N,K)=GI(NX,K)
               LO(N,K)=.TRUE.
               CYCLE N_LOOP
             ENDIF
           ENDIF
         ENDDO SPIRAL_LOOP
         IBO(K)=1
         GO(N,K)=0.
       ELSE
         IBO(K)=1
         GO(N,K)=0.
       ENDIF
     ELSE  ! no spiral search option
       IBO(K)=1
       GO(N,K)=0.
     ENDIF
   ENDDO N_LOOP
 ENDDO KM_LOOP
 IF(GDTNUMO.EQ.0) CALL POLFIXS(NO,MO,KM,RLAT,RLON,IBO,LO,GO)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 END SUBROUTINE POLATES3
