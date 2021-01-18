module polatev3_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  use ip_grid_descriptor_mod
  use ip_grid_factory_mod
  implicit none

  private
  public :: interpolate_budget_vector

contains

  SUBROUTINE interpolate_budget_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATEV3   INTERPOLATE VECTOR FIELDS (BUDGET)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS BUDGET INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
    !           THE ALGORITHM SIMPLY COMPUTES (WEIGHTED) AVERAGES
    !           OF BILINEARLY INTERPOLATED POINTS ARRANGED IN A SQUARE BOX
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
    !           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
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
    !           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
    !
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS AND
    !           THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
    !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS. ON THE OTHER
    !           THE OUTPUT CAN BE A SET OF STATION POINTS IF
    !           IGDTNUMO=IGDTNUMO-255, IN WHICH CASE THE NUMBER OF POINTS
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
    ! 2002-01-17  IREDELL  SAVE DATA FROM LAST CALL FOR OPTIMIZATION
    ! 2006-01-05  GAYNO    ADDED OPTION TO TO DO SUBSECTION OF OUTPUT GRID.
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
    !                      ROUTINE GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATEV3(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                     IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                     MI,MO,KM,IBI,LI,UI,VI, &
    !                     NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1) IS NUMBER OF RADIUS POINTS
    !                (DEFAULTS TO 2 IF IPOPT(1)=-1);
    !                IPOPT(2:2+IPOPT(1)) ARE RESPECTIVE WEIGHTS
    !                (DEFAULTS TO ALL 1 IF IPOPT(1)=-1 OR IPOPT(2)=-1).
    !                IPOPT(3+IPOPT(1)) IS MINIMUM PERCENTAGE FOR MASK
    !                (DEFAULTS TO 50 IF IPOPT(3+IPOPT(1)=-1)
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
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATEV FOR COMPLETE DEFINITION.
    !     IGDTLENI - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IGDTNUMO - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE. IGDTNUMO=IGDTNUM-255
    !                MEANS INTERPOLATE TO RANDOM STATION POINTS.
    !                OTHERWISE, SAME DEFINITION AS "IGDTNUMI".
    !     IGDTMPLO - INTEGER (IGDTLENO) GRID DEFINITION TEMPLATE ARRAY -
    !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
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
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
    !     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
    !     UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
    !     VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
    !     RLAT     - REAL (MO) INPUT LATITUDES IN DEGREES (IGDTNUMO<0)
    !     RLON     - REAL (MO) INPUT LONGITUDES IN DEGREES (IGDTNUMO<0)
    !
    !   OUTPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IGDTNUMO>=0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IGDTNUMO>=0)
    !     CROT     - REAL (MO) VECTOR ROTATION COSINES
    !     SROT     - REAL (MO) VECTOR ROTATION SINES
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
    !                32   INVALID BUDGET METHOD PARAMETERS
    !
    ! SUBPROGRAMS CALLED:
    !   CHECK_GRIDS3V CHECK IF GRID SPECS HAVE CHANGED
    !   GDSWZD        GRID DESCRIPTION SECTION WIZARD
    !   IJKGDS0       SET UP PARAMETERS FOR IJKGDS1
    !   IJKGDS1       RETURN FIELD POSITION FOR A GIVEN GRID POINT
    !   MOVECT        MOVE A VECTOR ALONG A GREAT CIRCLE
    !   POLFIXV       MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,          INTENT(IN   ) :: IPOPT(20), IBI(KM)
    INTEGER,          INTENT(IN   ) :: KM, MI, MO
    INTEGER,          INTENT(  OUT) :: IRET, NO, IBO(KM)
    !
    LOGICAL*1,        INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    REAL,             INTENT(  OUT) :: CROT(MO),SROT(MO)
    !
    REAL,             PARAMETER     :: FILL=-9999.
    !
    INTEGER                         :: IGDTNUMO2
    INTEGER                         :: I1,I2,J1,J2,IB,JB,LSW,MP
    INTEGER,          SAVE          :: MIX=-1
    INTEGER                         :: K,LB,N,NB,NB1,NB2,NB3,NB4,NV
    INTEGER                         :: N11(MO),N21(MO),N12(MO),N22(MO)
    !
    LOGICAL                         :: SAME_GRID
    !
    REAL                            :: CM11,SM11,CM12,SM12
    REAL                            :: CM21,SM21,CM22,SM22
    REAL                            :: PMP,RB2
    REAL                            :: C11(MO),C21(MO),C12(MO),C22(MO)
    REAL                            :: S11(MO),S21(MO),S12(MO),S22(MO)
    REAL                            :: W11(MO),W21(MO),W12(MO),W22(MO)
    REAL                            :: UB,VB,WB,UROT,VROT
    REAL                            :: U11,V11,U21,V21,U12,V12,U22,V22
    REAL                            :: WI1,WJ1,WI2,WJ2
    REAL                            :: WO(MO,KM),XI,YI
    REAL                            :: XPTS(MO),YPTS(MO)
    REAL                            :: XPTB(MO),YPTB(MO),RLOB(MO),RLAB(MO)
    REAL,         ALLOCATABLE, SAVE :: CROI(:),SROI(:)
    REAL,         ALLOCATABLE, SAVE :: XPTI(:),YPTI(:),RLOI(:),RLAI(:)

    class(ip_grid_descriptor), allocatable :: desc_out_subgrid
    class(ip_grid), allocatable :: subgrid

    IRET=0

    select type(grid_out)
    type is(ip_station_points_grid)
       desc_out_subgrid = grid_out%descriptor
       desc_out_subgrid%grid_num = 255 + grid_out%descriptor%grid_num

       subgrid = init_grid(desc_out_subgrid)
       CALL GDSWZD(subgrid,-1,MO,FILL,XPTS,YPTS, &
            RLON,RLAT,NO,CROT,SROT)
       IF(NO.EQ.0) IRET=3
       class default
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS, &
            RLON,RLAT,NO,CROT,SROT)
    end select

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    !CALL CHECK_GRIDS3V(IGDTNUMI,IGDTMPLI,IGDTLENI,SAME_GRID)
    same_grid = .false.
    IF(.NOT.SAME_GRID) THEN
       IF(MIX.NE.MI) THEN
          IF(MIX.GE.0) DEALLOCATE(XPTI,YPTI,RLOI,RLAI,CROI,SROI)
          ALLOCATE(XPTI(MI),YPTI(MI),RLOI(MI),RLAI(MI),CROI(MI),SROI(MI))
          MIX=MI
       ENDIF
       CALL GDSWZD(grid_in, 0,MI,FILL,XPTI,YPTI, &
            RLOI,RLAI,NV,CROI,SROI)
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
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
          UO(N,K)=0
          VO(N,K)=0
          WO(N,K)=0.
       ENDDO
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX
    DO NB=1,NB3
       !  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS AND ROTATIONS
       JB=(NB-1)/NB2-NB1
       IB=NB-(JB+NB1)*NB2-NB1-1
       LB=MAX(ABS(IB),ABS(JB))
       WB=1
       IF(IPOPT(2).EQ.-2) THEN
          WB=(NB1+1-ABS(IB))*(NB1+1-ABS(JB))
       ELSEIF(IPOPT(2).NE.-1) THEN
          WB=IPOPT(2+LB)
       ENDIF
       IF(WB.NE.0) THEN
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NO
             XPTB(N)=XPTS(N)+IB*RB2
             YPTB(N)=YPTS(N)+JB*RB2
          ENDDO
          !$OMP END PARALLEL DO
          CALL GDSWZD(grid_out, 1,NO,FILL,XPTB,YPTB, &
               RLOB,RLAB,NV)
          CALL GDSWZD(grid_in,-1,NO,FILL,XPTB,YPTB, &
               RLOB,RLAB,NV)
          IF(IRET.EQ.0.AND.NV.EQ.0.AND.LB.EQ.0) IRET=2
          !$OMP PARALLEL DO PRIVATE(N,XI,YI,I1,I2,WI1,WI2,J1,J2,WJ1,WJ2,CM11,CM21,CM12,CM22,SM11,SM21,SM12,SM22) &
          !$OMP SCHEDULE(STATIC)
          DO N=1,NO
             XI=XPTB(N)
             YI=YPTB(N)
             IF(XI.NE.FILL.AND.YI.NE.FILL) THEN
                I1=XI
                I2=I1+1
                WI2=XI-I1
                WI1=1-WI2
                J1=YI
                J2=J1+1
                WJ2=YI-J1
                WJ1=1-WJ2
                N11(N) = grid_in%field_pos(i1,j1)
                N21(N) = grid_in%field_pos(i2, j1)
                N12(N) = grid_in%field_pos(i1, j2)
                N22(N) = grid_in%field_pos(i2, j2)
                IF(MIN(N11(N),N21(N),N12(N),N22(N)).GT.0) THEN
                   W11(N)=WI1*WJ1
                   W21(N)=WI2*WJ1
                   W12(N)=WI1*WJ2
                   W22(N)=WI2*WJ2
                   CALL MOVECT(RLAI(N11(N)),RLOI(N11(N)),RLAT(N),RLON(N),CM11,SM11)
                   CALL MOVECT(RLAI(N21(N)),RLOI(N21(N)),RLAT(N),RLON(N),CM21,SM21)
                   CALL MOVECT(RLAI(N12(N)),RLOI(N12(N)),RLAT(N),RLON(N),CM12,SM12)
                   CALL MOVECT(RLAI(N22(N)),RLOI(N22(N)),RLAT(N),RLON(N),CM22,SM22)
                   C11(N)=CM11*CROI(N11(N))+SM11*SROI(N11(N))
                   S11(N)=SM11*CROI(N11(N))-CM11*SROI(N11(N))
                   C21(N)=CM21*CROI(N21(N))+SM21*SROI(N21(N))
                   S21(N)=SM21*CROI(N21(N))-CM21*SROI(N21(N))
                   C12(N)=CM12*CROI(N12(N))+SM12*SROI(N12(N))
                   S12(N)=SM12*CROI(N12(N))-CM12*SROI(N12(N))
                   C22(N)=CM22*CROI(N22(N))+SM22*SROI(N22(N))
                   S22(N)=SM22*CROI(N22(N))-CM22*SROI(N22(N))
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
          !$OMP END PARALLEL DO 
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !  INTERPOLATE WITH OR WITHOUT BITMAPS
          !  KM IS OFTEN 1 .. DO NO PUT OMP PARALLEL DO HERE
          DO K=1,KM
             !$OMP PARALLEL DO PRIVATE(N,U11,U12,U21,U22,UB,V11,V12,V21,V22,VB) SCHEDULE(STATIC)
             DO N=1,NO
                IF(N11(N).GT.0) THEN
                   IF(IBI(K).EQ.0) THEN
                      U11=C11(N)*UI(N11(N),K)-S11(N)*VI(N11(N),K)
                      V11=S11(N)*UI(N11(N),K)+C11(N)*VI(N11(N),K)
                      U21=C21(N)*UI(N21(N),K)-S21(N)*VI(N21(N),K)
                      V21=S21(N)*UI(N21(N),K)+C21(N)*VI(N21(N),K)
                      U12=C12(N)*UI(N12(N),K)-S12(N)*VI(N12(N),K)
                      V12=S12(N)*UI(N12(N),K)+C12(N)*VI(N12(N),K)
                      U22=C22(N)*UI(N22(N),K)-S22(N)*VI(N22(N),K)
                      V22=S22(N)*UI(N22(N),K)+C22(N)*VI(N22(N),K)
                      UB=W11(N)*U11+W21(N)*U21+W12(N)*U12+W22(N)*U22
                      VB=W11(N)*V11+W21(N)*V21+W12(N)*V12+W22(N)*V22
                      UO(N,K)=UO(N,K)+WB*UB
                      VO(N,K)=VO(N,K)+WB*VB
                      WO(N,K)=WO(N,K)+WB
                   ELSE
                      IF(LI(N11(N),K)) THEN
                         U11=C11(N)*UI(N11(N),K)-S11(N)*VI(N11(N),K)
                         V11=S11(N)*UI(N11(N),K)+C11(N)*VI(N11(N),K)
                         UO(N,K)=UO(N,K)+WB*W11(N)*U11
                         VO(N,K)=VO(N,K)+WB*W11(N)*V11
                         WO(N,K)=WO(N,K)+WB*W11(N)
                      ENDIF
                      IF(LI(N21(N),K)) THEN
                         U21=C21(N)*UI(N21(N),K)-S21(N)*VI(N21(N),K)
                         V21=S21(N)*UI(N21(N),K)+C21(N)*VI(N21(N),K)
                         UO(N,K)=UO(N,K)+WB*W21(N)*U21
                         VO(N,K)=VO(N,K)+WB*W21(N)*V21
                         WO(N,K)=WO(N,K)+WB*W21(N)
                      ENDIF
                      IF(LI(N12(N),K)) THEN
                         U12=C12(N)*UI(N12(N),K)-S12(N)*VI(N12(N),K)
                         V12=S12(N)*UI(N12(N),K)+C12(N)*VI(N12(N),K)
                         UO(N,K)=UO(N,K)+WB*W12(N)*U12
                         VO(N,K)=VO(N,K)+WB*W12(N)*V12
                         WO(N,K)=WO(N,K)+WB*W12(N)
                      ENDIF
                      IF(LI(N22(N),K)) THEN
                         U22=C22(N)*UI(N22(N),K)-S22(N)*VI(N22(N),K)
                         V22=S22(N)*UI(N22(N),K)+C22(N)*VI(N22(N),K)
                         UO(N,K)=UO(N,K)+WB*W22(N)*U22
                         VO(N,K)=VO(N,K)+WB*W22(N)*V22
                         WO(N,K)=WO(N,K)+WB*W22(N)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
             !$OMP END PARALLEL DO
          ENDDO
       ENDIF
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE OUTPUT BITMAPS AND FIELDS
    ! KM is often 1, do not put OMP PARALLEL here
    DO K=1,KM
       IBO(K)=IBI(K)
       !$OMP PARALLEL DO PRIVATE(N,UROT,VROT) SCHEDULE(STATIC)
       DO N=1,NO
          LO(N,K)=WO(N,K).GE.PMP*NB4
          IF(LO(N,K)) THEN
             UO(N,K)=UO(N,K)/WO(N,K)
             VO(N,K)=VO(N,K)/WO(N,K)
             UROT=CROT(N)*UO(N,K)-SROT(N)*VO(N,K)
             VROT=SROT(N)*UO(N,K)+CROT(N)*VO(N,K)
             UO(N,K)=UROT
             VO(N,K)=VROT
          ELSE
             IBO(K)=1
             UO(N,K)=0.
             VO(N,K)=0.
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
    ENDDO

    select type(grid_out)
    type is(ip_equid_cylind_grid)
       CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE INTERPOLATE_BUDGET_VECTOR
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE CHECK_GRIDS3V(IGDTNUM,IGDTMPL,IGDTLEN,SAME_GRID)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  CHECK_GRIDS3V   CHECK GRID INFORMATION
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: DETERMINE WHETHER THE GRID SPECS
    !           HAVE CHANGED.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:  CALL CHECK_GRIDS3V(IGDTNUM,IGDTMPL,IGDTLEN,SAME_GRID)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !               CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !               NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IGDTMPL - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY.
    !               CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !               OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IGDTLEN - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !               TEMPLATE ARRAY.  CORRESPONDS TO THE GFLD%IGDTLEN
    !               COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !
    !   OUTPUT ARGUMENT LIST:
    !     SAME_GRID  - WHEN TRUE, THE GRID HAS NOT CHANGED BETWEEN CALLS.
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    INTEGER,        INTENT(IN   ) :: IGDTNUM, IGDTLEN
    INTEGER,        INTENT(IN   ) :: IGDTMPL(IGDTLEN)
    !
    LOGICAL,        INTENT(  OUT) :: SAME_GRID
    !
    INTEGER, SAVE                 :: IGDTNUM_SAVE=-9999
    INTEGER, SAVE                 :: IGDTLEN_SAVE=-9999
    INTEGER, SAVE                 :: IGDTMPL_SAVE(1000)=-9999
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    SAME_GRID=.FALSE.
    IF(IGDTNUM==IGDTNUM_SAVE)THEN
       IF(IGDTLEN==IGDTLEN_SAVE)THEN
          IF(ALL(IGDTMPL==IGDTMPL_SAVE(1:IGDTLEN)))THEN
             SAME_GRID=.TRUE.
          ENDIF
       ENDIF
    ENDIF
    !
    IGDTNUM_SAVE=IGDTNUM
    IGDTLEN_SAVE=IGDTLEN
    IGDTMPL_SAVE(1:IGDTLEN)=IGDTMPL
    IGDTMPL_SAVE(IGDTLEN+1:1000)=-9999
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE CHECK_GRIDS3V

end module polatev3_mod
