module polatev0_mod
  use ijkgds_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  implicit none

  private
  public :: interpolate_bilinear_vector
  

  ! interface polatev0
  !    module procedure polatev0_grib1
  !    module procedure polatev0_grib2
  ! end interface polatev0

  class(ip_grid), allocatable :: prev_grid_in, prev_grid_out

  INTEGER,                    SAVE  :: NOX=-1,IRETX=-1
  INTEGER,        ALLOCATABLE,SAVE  :: NXY(:,:,:)
  REAL,           ALLOCATABLE,SAVE  :: RLATX(:),RLONX(:)
  REAL,           ALLOCATABLE,SAVE  :: CROTX(:),SROTX(:)
  REAL,           ALLOCATABLE,SAVE  :: WXY(:,:,:),CXY(:,:,:),SXY(:,:,:)

contains

  SUBROUTINE interpolate_bilinear_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATEV0   INTERPOLATE VECTOR FIELDS (BILINEAR)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS BILINEAR INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
    !           OPTIONS ALLOW VARYING THE MINIMUM PERCENTAGE FOR MASK,
    !           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
    !           (IPOPT(1)) WHICH DEFAULTS TO 50 (IF IPOPT(1)=-1).
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
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
    !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
    !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !           ON THE OTHER HAND, THE DATA MAY BE INTERPOLATED TO A SET OF
    !           STATION POINTS IF IGDTNUMO<0, IN WHICH CASE THE NUMBER
    !           OF POINTS AND THEIR LATITUDES AND LONGITUDES MUST BE
    !           INPUT ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !
    !           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
    !           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
    !           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
    !           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    ! 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
    ! 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
    ! 2002-01-17  IREDELL  SAVE DATA FROM LAST CALL FOR OPTIMIZATION
    ! 2007-05-22  IREDELL  EXTRAPOLATE UP TO HALF A GRID CELL
    ! 2007-10-30  IREDELL  SAVE WEIGHTS AND THREAD FOR PERFORMANCE
    ! 2012-06-26  GAYNO    FIX OUT-OF-BOUNDS ERROR.  SEE NCEPLIBS
    !                      TICKET #9.
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH MERGED VERSION
    !                      OF GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATEV0(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,UI,VI, &
    !                         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1) IS MINIMUM PERCENTAGE FOR MASK
    !                (DEFAULTS TO 50 IF IPOPT(1)=-1)
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
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO<0)
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
    !
    ! SUBPROGRAMS CALLED:
    !   GDSWZD        GRID DESCRIPTION SECTION WIZARD
    !   IJKGDS0       SET UP PARAMETERS FOR IJKGDS1
    !   IJKGDS1       RETURN FIELD POSITION FOR A GIVEN GRID POINT
    !   MOVECT        MOVE A VECTOR ALONG A GREAT CIRCLE
    !   POLFIXV       MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
    !   CHECK_GRIDS0V DETERMINE IF INPUT OR OUTPUT GRIDS HAVE CHANGED
    !                 BETWEEN CALLS TO THIS ROUTINE.
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,            INTENT(IN   ) :: IPOPT(20),IBI(KM),MI,MO,KM
    INTEGER,            INTENT(INOUT) :: NO
    INTEGER,            INTENT(  OUT) :: IRET, IBO(KM)
    !
    LOGICAL*1,          INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,          INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,               INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,               INTENT(INOUT) :: RLAT(MO),RLON(MO),CROT(MO),SROT(MO)
    REAL,               INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    !
    REAL,               PARAMETER     :: FILL=-9999.
    !
    INTEGER                           :: IJX(2),IJY(2)
    INTEGER                           :: MP,N,I,J,K,NK,NV
    !
    LOGICAL                           :: SAME_GRIDI, SAME_GRIDO
    !
    REAL                              :: CM,SM,UROT,VROT
    REAL                              :: PMP,XIJ,YIJ,XF,YF,U,V,W
    REAL                              :: XPTS(MO),YPTS(MO)
    REAL                              :: WX(2),WY(2)
    REAL                              :: XPTI(MI),YPTI(MI)
    REAL                              :: RLOI(MI),RLAI(MI)
    REAL                              :: CROI(MI),SROI(MI)

    logical :: to_station_points
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IRET=0
    MP=IPOPT(1)
    IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
    IF(MP.LT.0.OR.MP.GT.100) IRET=32
    PMP=MP*0.01

    if (.not. allocated(prev_grid_in) .or. .not. allocated(prev_grid_out)) then
       allocate(prev_grid_in, source = grid_in)
       allocate(prev_grid_out, source = grid_out)

       same_gridi = .false.
       same_grido = .false.
    else
       same_gridi = grid_in == prev_grid_in
       same_grido = grid_out == prev_grid_out

       if (.not. same_gridi .or. .not. same_grido) then
          deallocate(prev_grid_in)
          deallocate(prev_grid_out)

          allocate(prev_grid_in, source = grid_in)
          allocate(prev_grid_out, source = grid_out)
       end if
    end if

    select type(grid_out)
    type is(ip_station_points_grid)
       to_station_points = .true.
       class default
       to_station_points = .false.
    end select
   
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SAVE OR SKIP WEIGHT COMPUTATION
    IF(IRET.EQ.0.AND.(to_station_points.OR..NOT.SAME_GRIDI.OR..NOT.SAME_GRIDO))THEN
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
       IF(.not. to_station_points) THEN
          CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT, &
               NO,CROT,SROT)
          IF(NO.EQ.0) IRET=3
       ENDIF
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  LOCATE INPUT POINTS
       CALL GDSWZD(grid_in,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
       IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
       CALL GDSWZD(grid_in, 0,MI,FILL,XPTI,YPTI,RLOI,RLAI,NV,&
            CROI,SROI)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  ALLOCATE AND SAVE GRID DATA
       IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,CROTX,SROTX,NXY,WXY,CXY,SXY)
          ALLOCATE(RLATX(NO),RLONX(NO),CROTX(NO),SROTX(NO), &
               NXY(2,2,NO),WXY(2,2,NO),CXY(2,2,NO),SXY(2,2,NO))
          NOX=NO
       ENDIF
       IRETX=IRET
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE WEIGHTS
       IF(IRET.EQ.0) THEN
          !CALL IJKGDS0(IGDTNUMI,IGDTMPLI,IGDTLENI,IJKGDSA)
          !$OMP PARALLEL DO PRIVATE(N,XIJ,YIJ,IJX,IJY,XF,YF,J,I,WX,WY,CM,SM) SCHEDULE(STATIC)
          DO N=1,NO
             RLONX(N)=RLON(N)
             RLATX(N)=RLAT(N)
             CROTX(N)=CROT(N)
             SROTX(N)=SROT(N)
             XIJ=XPTS(N)
             YIJ=YPTS(N)
             IF(XIJ.NE.FILL.AND.YIJ.NE.FILL) THEN
                IJX(1:2)=FLOOR(XIJ)+(/0,1/)
                IJY(1:2)=FLOOR(YIJ)+(/0,1/)
                XF=XIJ-IJX(1)
                YF=YIJ-IJY(1)
                WX(1)=(1-XF)
                WX(2)=XF
                WY(1)=(1-YF)
                WY(2)=YF
                DO J=1,2
                   DO I=1,2
                      !NXY(I,J,N)=IJKGDS1(IJX(I),IJY(J),IJKGDSA)
                      nxy(i, j, n) = grid_in%field_pos(ijx(i), ijy(j))
                      WXY(I,J,N)=WX(I)*WY(J)
                      IF(NXY(I,J,N).GT.0) THEN
                         CALL MOVECT(RLAI(NXY(I,J,N)),RLOI(NXY(I,J,N)), &
                              RLAT(N),RLON(N),CM,SM)
                         CXY(I,J,N)=CM*CROI(NXY(I,J,N))+SM*SROI(NXY(I,J,N))
                         SXY(I,J,N)=SM*CROI(NXY(I,J,N))-CM*SROI(NXY(I,J,N))
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                NXY(:,:,N)=0
             ENDIF
          ENDDO
       ENDIF  ! IS IRET 0?
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE OVER ALL FIELDS
    IF(IRET.EQ.0.AND.IRETX.EQ.0) THEN
       IF(.not. to_station_points) THEN
          NO=NOX
          DO N=1,NO
             RLON(N)=RLONX(N)
             RLAT(N)=RLATX(N)
             CROT(N)=CROTX(N)
             SROT(N)=SROTX(N)
          ENDDO
       ENDIF
       !$OMP PARALLEL DO PRIVATE(NK,K,N,U,V,W,UROT,VROT,J,I) SCHEDULE(STATIC)
       DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          U=0
          V=0
          W=0
          DO J=1,2
             DO I=1,2
                IF(NXY(I,J,N).GT.0) THEN
                   IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K)) THEN
                      UROT=CXY(I,J,N)*UI(NXY(I,J,N),K)-SXY(I,J,N)*VI(NXY(I,J,N),K)
                      VROT=SXY(I,J,N)*UI(NXY(I,J,N),K)+CXY(I,J,N)*VI(NXY(I,J,N),K)
                      U=U+WXY(I,J,N)*UROT
                      V=V+WXY(I,J,N)*VROT
                      W=W+WXY(I,J,N)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          LO(N,K)=W.GE.PMP
          IF(LO(N,K)) THEN
             UROT=CROT(N)*U-SROT(N)*V
             VROT=SROT(N)*U+CROT(N)*V
             UO(N,K)=UROT/W
             VO(N,K)=VROT/W
          ELSE
             UO(N,K)=0.
             VO(N,K)=0.
          ENDIF
       ENDDO  ! NK LOOP
       DO K=1,KM
          IBO(K)=IBI(K)
          IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
       ENDDO
       
       select type(grid_out)
       type is(ip_equid_cylind_grid)
          CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
       end select

    ELSE
       IF(IRET.EQ.0) IRET=IRETX
       IF(.not. to_station_points) NO=0
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE INTERPOLATE_BILINEAR_VECTOR

  ! SUBROUTINE POLATEV0_grib2(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
  !      IGDTNUMO,IGDTMPLO,IGDTLENO, &
  !      MI,MO,KM,IBI,LI,UI,VI, &
  !      NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
  !   !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !   !
  !   ! SUBPROGRAM:  POLATEV0   INTERPOLATE VECTOR FIELDS (BILINEAR)
  !   !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
  !   !
  !   ! ABSTRACT: THIS SUBPROGRAM PERFORMS BILINEAR INTERPOLATION
  !   !           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
  !   !           OPTIONS ALLOW VARYING THE MINIMUM PERCENTAGE FOR MASK,
  !   !           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
  !   !           (IPOPT(1)) WHICH DEFAULTS TO 50 (IF IPOPT(1)=-1).
  !   !           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
  !   !
  !   !           THE INPUT AND OUTPUT GRIDS ARE DEFINED BY THEIR GRIB 2 GRID
  !   !           DEFINITION TEMPLATE AS DECODED BY THE NCEP G2 LIBRARY.  THE
  !   !           CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
  !   !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
  !   !           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
  !   !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
  !   !             (IGDTNUMI/O=01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
  !   !                             NON-"E" STAGGERED
  !   !             (IGDTNUMI/O=10) MERCATOR CYLINDRICAL
  !   !             (IGDTNUMI/O=20) POLAR STEREOGRAPHIC AZIMUTHAL
  !   !             (IGDTNUMI/O=30) LAMBERT CONFORMAL CONICAL
  !   !             (IGDTNUMI/O=40) GAUSSIAN CYLINDRICAL
  !   !
  !   !           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
  !   !           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
  !   !           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
  !   !           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
  !   !           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
  !   !
  !   !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
  !   !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
  !   !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  !   !           ON THE OTHER HAND, THE DATA MAY BE INTERPOLATED TO A SET OF
  !   !           STATION POINTS IF IGDTNUMO<0, IN WHICH CASE THE NUMBER
  !   !           OF POINTS AND THEIR LATITUDES AND LONGITUDES MUST BE
  !   !           INPUT ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  !   !
  !   !           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
  !   !           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
  !   !           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
  !   !           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
  !   !        
  !   ! PROGRAM HISTORY LOG:
  !   !   96-04-10  IREDELL
  !   ! 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
  !   ! 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
  !   ! 2002-01-17  IREDELL  SAVE DATA FROM LAST CALL FOR OPTIMIZATION
  !   ! 2007-05-22  IREDELL  EXTRAPOLATE UP TO HALF A GRID CELL
  !   ! 2007-10-30  IREDELL  SAVE WEIGHTS AND THREAD FOR PERFORMANCE
  !   ! 2012-06-26  GAYNO    FIX OUT-OF-BOUNDS ERROR.  SEE NCEPLIBS
  !   !                      TICKET #9.
  !   ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH MERGED VERSION
  !   !                      OF GDSWZD.
  !   ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
  !   !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
  !   !
  !   ! USAGE:    CALL POLATEV0(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
  !   !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
  !   !                         MI,MO,KM,IBI,LI,UI,VI, &
  !   !                         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
  !   !
  !   !   INPUT ARGUMENT LIST:
  !   !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
  !   !                IPOPT(1) IS MINIMUM PERCENTAGE FOR MASK
  !   !                (DEFAULTS TO 50 IF IPOPT(1)=-1)
  !   !     IGDTNUMI - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
  !   !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
  !   !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE:
  !   !                  00 - EQUIDISTANT CYLINDRICAL
  !   !                  01 - ROTATED EQUIDISTANT CYLINDRICAL.  "E"
  !   !                       AND NON-"E" STAGGERED
  !   !                  10 - MERCATOR CYCLINDRICAL
  !   !                  20 - POLAR STEREOGRAPHIC AZIMUTHAL
  !   !                  30 - LAMBERT CONFORMAL CONICAL
  !   !                  40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
  !   !     IGDTMPLI - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY -
  !   !                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
  !   !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE
  !   !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
  !   !                IPOLATEV FOR COMPLETE DEFINITION.
  !   !     IGDTLENI - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
  !   !                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
  !   !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !     IGDTNUMO - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
  !   !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
  !   !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE. IGDTNUMO<0
  !   !                MEANS INTERPOLATE TO RANDOM STATION POINTS.
  !   !                OTHERWISE, SAME DEFINITION AS "IGDTNUMI".
  !   !     IGDTMPLO - INTEGER (IGDTLENO) GRID DEFINITION TEMPLATE ARRAY -
  !   !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
  !   !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
  !   !                IPOLATEV FOR COMPLETE DEFINITION.
  !   !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
  !   !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
  !   !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
  !   !                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
  !   !     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
  !   !                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
  !   !     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
  !   !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
  !   !     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
  !   !     UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
  !   !     VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
  !   !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO<0)
  !   !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO<0)
  !   !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO<0)
  !   !     CROT     - REAL (MO) VECTOR ROTATION COSINES (IF IGDTNUMO<0)
  !   !     SROT     - REAL (MO) VECTOR ROTATION SINES (IF IGDTNUMO<0)
  !   !                (UGRID=CROT*UEARTH-SROT*VEARTH;
  !   !                 VGRID=SROT*UEARTH+CROT*VEARTH)
  !   !
  !   !   OUTPUT ARGUMENT LIST:
  !   !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO>=0)
  !   !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO>=0)
  !   !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO>=0)
  !   !     CROT     - REAL (MO) VECTOR ROTATION COSINES (IF IGDTNUMO>=0)
  !   !     SROT     - REAL (MO) VECTOR ROTATION SINES (IF IGDTNUMO>=0)
  !   !                (UGRID=CROT*UEARTH-SROT*VEARTH;
  !   !                 VGRID=SROT*UEARTH+CROT*VEARTH)
  !   !     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
  !   !     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
  !   !     UO       - REAL (MO,KM) OUTPUT U-COMPONENT FIELDS INTERPOLATED
  !   !     VO       - REAL (MO,KM) OUTPUT V-COMPONENT FIELDS INTERPOLATED
  !   !     IRET     - INTEGER RETURN CODE
  !   !                0    SUCCESSFUL INTERPOLATION
  !   !                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
  !   !                3    UNRECOGNIZED OUTPUT GRID
  !   !
  !   ! SUBPROGRAMS CALLED:
  !   !   GDSWZD        GRID DESCRIPTION SECTION WIZARD
  !   !   IJKGDS0       SET UP PARAMETERS FOR IJKGDS1
  !   !   IJKGDS1       RETURN FIELD POSITION FOR A GIVEN GRID POINT
  !   !   MOVECT        MOVE A VECTOR ALONG A GREAT CIRCLE
  !   !   POLFIXV       MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
  !   !   CHECK_GRIDS0V DETERMINE IF INPUT OR OUTPUT GRIDS HAVE CHANGED
  !   !                 BETWEEN CALLS TO THIS ROUTINE.
  !   !
  !   ! ATTRIBUTES:
  !   !   LANGUAGE: FORTRAN 90
  !   !
  !   !$$$
  !   INTEGER,            INTENT(IN   ) :: IPOPT(20),IBI(KM),MI,MO,KM
  !   INTEGER,            INTENT(IN   ) :: IGDTNUMI, IGDTLENI
  !   INTEGER,            INTENT(IN   ) :: IGDTMPLI(IGDTLENI)
  !   INTEGER,            INTENT(IN   ) :: IGDTNUMO, IGDTLENO
  !   INTEGER,            INTENT(IN   ) :: IGDTMPLO(IGDTLENO)
  !   INTEGER,            INTENT(INOUT) :: NO
  !   INTEGER,            INTENT(  OUT) :: IRET, IBO(KM)
  !   !
  !   LOGICAL*1,          INTENT(IN   ) :: LI(MI,KM)
  !   LOGICAL*1,          INTENT(  OUT) :: LO(MO,KM)
  !   !
  !   REAL,               INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
  !   REAL,               INTENT(INOUT) :: RLAT(MO),RLON(MO),CROT(MO),SROT(MO)
  !   REAL,               INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
  !   !
  !   REAL,               PARAMETER     :: FILL=-9999.
  !   !
  !   INTEGER                           :: IJX(2),IJY(2),IJKGDSA(20)
  !   INTEGER                           :: MP,N,I,J,K,NK,NV,IJKGDS1
  !   INTEGER,                    SAVE  :: NOX=-1,IRETX=-1
  !   INTEGER,        ALLOCATABLE,SAVE  :: NXY(:,:,:)
  !   !
  !   LOGICAL                           :: SAME_GRIDI, SAME_GRIDO
  !   !
  !   REAL                              :: CM,SM,UROT,VROT
  !   REAL                              :: PMP,XIJ,YIJ,XF,YF,U,V,W
  !   REAL                              :: XPTS(MO),YPTS(MO)
  !   REAL                              :: WX(2),WY(2)
  !   REAL                              :: XPTI(MI),YPTI(MI)
  !   REAL                              :: RLOI(MI),RLAI(MI)
  !   REAL                              :: CROI(MI),SROI(MI)
  !   REAL,           ALLOCATABLE,SAVE  :: RLATX(:),RLONX(:)
  !   REAL,           ALLOCATABLE,SAVE  :: CROTX(:),SROTX(:)
  !   REAL,           ALLOCATABLE,SAVE  :: WXY(:,:,:),CXY(:,:,:),SXY(:,:,:)
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   !  SET PARAMETERS
  !   IRET=0
  !   MP=IPOPT(1)
  !   IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
  !   IF(MP.LT.0.OR.MP.GT.100) IRET=32
  !   PMP=MP*0.01
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   CALL CHECK_GRIDS0V(IGDTNUMI,IGDTMPLI,IGDTLENI, &
  !        IGDTNUMO,IGDTMPLO,IGDTLENO, &
  !        SAME_GRIDI,SAME_GRIDO)
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   !  SAVE OR SKIP WEIGHT COMPUTATION
  !   IF(IRET.EQ.0.AND.(IGDTNUMO.LT.0.OR..NOT.SAME_GRIDI.OR..NOT.SAME_GRIDO))THEN
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
  !      IF(IGDTNUMO.GE.0) THEN
  !         CALL GDSWZD(IGDTNUMO,IGDTMPLO,IGDTLENO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT, &
  !              NO,CROT,SROT)
  !         IF(NO.EQ.0) IRET=3
  !      ENDIF
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  LOCATE INPUT POINTS
  !      CALL GDSWZD(IGDTNUMI,IGDTMPLI,IGDTLENI,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
  !      IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
  !      CALL GDSWZD(IGDTNUMI,IGDTMPLI,IGDTLENI, 0,MI,FILL,XPTI,YPTI,RLOI,RLAI,NV,&
  !           CROI,SROI)
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  ALLOCATE AND SAVE GRID DATA
  !      IF(NOX.NE.NO) THEN
  !         IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,CROTX,SROTX,NXY,WXY,CXY,SXY)
  !         ALLOCATE(RLATX(NO),RLONX(NO),CROTX(NO),SROTX(NO), &
  !              NXY(2,2,NO),WXY(2,2,NO),CXY(2,2,NO),SXY(2,2,NO))
  !         NOX=NO
  !      ENDIF
  !      IRETX=IRET
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  COMPUTE WEIGHTS
  !      IF(IRET.EQ.0) THEN
  !         CALL IJKGDS0(IGDTNUMI,IGDTMPLI,IGDTLENI,IJKGDSA)
  !         !$OMP PARALLEL DO PRIVATE(N,XIJ,YIJ,IJX,IJY,XF,YF,J,I,WX,WY,CM,SM) SCHEDULE(STATIC)
  !         DO N=1,NO
  !            RLONX(N)=RLON(N)
  !            RLATX(N)=RLAT(N)
  !            CROTX(N)=CROT(N)
  !            SROTX(N)=SROT(N)
  !            XIJ=XPTS(N)
  !            YIJ=YPTS(N)
  !            IF(XIJ.NE.FILL.AND.YIJ.NE.FILL) THEN
  !               IJX(1:2)=FLOOR(XIJ)+(/0,1/)
  !               IJY(1:2)=FLOOR(YIJ)+(/0,1/)
  !               XF=XIJ-IJX(1)
  !               YF=YIJ-IJY(1)
  !               WX(1)=(1-XF)
  !               WX(2)=XF
  !               WY(1)=(1-YF)
  !               WY(2)=YF
  !               DO J=1,2
  !                  DO I=1,2
  !                     NXY(I,J,N)=IJKGDS1(IJX(I),IJY(J),IJKGDSA)
  !                     WXY(I,J,N)=WX(I)*WY(J)
  !                     IF(NXY(I,J,N).GT.0) THEN
  !                        CALL MOVECT(RLAI(NXY(I,J,N)),RLOI(NXY(I,J,N)), &
  !                             RLAT(N),RLON(N),CM,SM)
  !                        CXY(I,J,N)=CM*CROI(NXY(I,J,N))+SM*SROI(NXY(I,J,N))
  !                        SXY(I,J,N)=SM*CROI(NXY(I,J,N))-CM*SROI(NXY(I,J,N))
  !                     ENDIF
  !                  ENDDO
  !               ENDDO
  !            ELSE
  !               NXY(:,:,N)=0
  !            ENDIF
  !         ENDDO
  !      ENDIF  ! IS IRET 0?
  !   ENDIF
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   !  INTERPOLATE OVER ALL FIELDS
  !   IF(IRET.EQ.0.AND.IRETX.EQ.0) THEN
  !      IF(IGDTNUMO.GE.0) THEN
  !         NO=NOX
  !         DO N=1,NO
  !            RLON(N)=RLONX(N)
  !            RLAT(N)=RLATX(N)
  !            CROT(N)=CROTX(N)
  !            SROT(N)=SROTX(N)
  !         ENDDO
  !      ENDIF
  !      !$OMP PARALLEL DO PRIVATE(NK,K,N,U,V,W,UROT,VROT,J,I) SCHEDULE(STATIC)
  !      DO NK=1,NO*KM
  !         K=(NK-1)/NO+1
  !         N=NK-NO*(K-1)
  !         U=0
  !         V=0
  !         W=0
  !         DO J=1,2
  !            DO I=1,2
  !               IF(NXY(I,J,N).GT.0) THEN
  !                  IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K)) THEN
  !                     UROT=CXY(I,J,N)*UI(NXY(I,J,N),K)-SXY(I,J,N)*VI(NXY(I,J,N),K)
  !                     VROT=SXY(I,J,N)*UI(NXY(I,J,N),K)+CXY(I,J,N)*VI(NXY(I,J,N),K)
  !                     U=U+WXY(I,J,N)*UROT
  !                     V=V+WXY(I,J,N)*VROT
  !                     W=W+WXY(I,J,N)
  !                  ENDIF
  !               ENDIF
  !            ENDDO
  !         ENDDO
  !         LO(N,K)=W.GE.PMP
  !         IF(LO(N,K)) THEN
  !            UROT=CROT(N)*U-SROT(N)*V
  !            VROT=SROT(N)*U+CROT(N)*V
  !            UO(N,K)=UROT/W
  !            VO(N,K)=VROT/W
  !         ELSE
  !            UO(N,K)=0.
  !            VO(N,K)=0.
  !         ENDIF
  !      ENDDO  ! NK LOOP
  !      DO K=1,KM
  !         IBO(K)=IBI(K)
  !         IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
  !      ENDDO
  !      IF(IGDTNUMO.EQ.0) CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   ELSE
  !      IF(IRET.EQ.0) IRET=IRETX
  !      IF(IGDTNUMO.GE.0) NO=0
  !   ENDIF
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! END SUBROUTINE POLATEV0_GRIB2
  ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! SUBROUTINE CHECK_GRIDS0V(IGDTNUMI,IGDTMPLI,IGDTLENI, &
  !      IGDTNUMO,IGDTMPLO,IGDTLENO, &
  !      SAME_GRIDI, SAME_GRIDO)
  !   !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !   !
  !   ! SUBPROGRAM:  CHECK_GRIDS0V   CHECK GRID INFORMATION
  !   !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
  !   !
  !   ! ABSTRACT: DETERMINE WHETHER THE INPUT OR OUTPUT GRID SPECS
  !   !           HAVE CHANGED.
  !   !
  !   ! PROGRAM HISTORY LOG:
  !   ! 2015-07-13  GAYNO     INITIAL VERSION
  !   !
  !   ! USAGE:  CALL CHECK_GRIDS0V(IGDTNUMI,IGDTMPLI,IGDTLENI,IGDTNUMO,IGDTMPLO, &
  !   !                            IGDTLENO, SAME_GRIDI, SAME_GRIDO)
  !   !
  !   !   INPUT ARGUMENT LIST:
  !   !     IGDTNUMI - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
  !   !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
  !   !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !     IGDTMPLI - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY -
  !   !                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
  !   !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !     IGDTLENI - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
  !   !                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
  !   !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !     IGDTNUMO - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
  !   !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
  !   !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !     IGDTMPLO - INTEGER (GDTLENO) GRID DEFINITION TEMPLATE ARRAY -
  !   !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
  !   !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
  !   !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
  !   !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
  !   !
  !   !   OUTPUT ARGUMENT LIST:
  !   !     SAME_GRIDI  - WHEN TRUE, THE INPUT GRID HAS NOT CHANGED BETWEEN CALLS.
  !   !     SAME_GRIDO  - WHEN TRUE, THE OUTPUT GRID HAS NOT CHANGED BETWEEN CALLS.
  !   !
  !   ! ATTRIBUTES:
  !   !   LANGUAGE: FORTRAN 90
  !   !
  !   !$$$
  !   IMPLICIT NONE
  !   !
  !   INTEGER,        INTENT(IN   ) :: IGDTNUMI, IGDTLENI
  !   INTEGER,        INTENT(IN   ) :: IGDTMPLI(IGDTLENI)
  !   INTEGER,        INTENT(IN   ) :: IGDTNUMO, IGDTLENO
  !   INTEGER,        INTENT(IN   ) :: IGDTMPLO(IGDTLENO)
  !   !
  !   LOGICAL,        INTENT(  OUT) :: SAME_GRIDI, SAME_GRIDO
  !   !
  !   INTEGER, SAVE                 :: IGDTNUMI_SAVE=-9999
  !   INTEGER, SAVE                 :: IGDTLENI_SAVE=-9999
  !   INTEGER, SAVE                 :: IGDTMPLI_SAVE(1000)=-9999
  !   INTEGER, SAVE                 :: IGDTNUMO_SAVE=-9999
  !   INTEGER, SAVE                 :: IGDTLENO_SAVE=-9999
  !   INTEGER, SAVE                 :: IGDTMPLO_SAVE(1000)=-9999
  !   !
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   SAME_GRIDI=.FALSE.
  !   IF(IGDTNUMI==IGDTNUMI_SAVE)THEN
  !      IF(IGDTLENI==IGDTLENI_SAVE)THEN
  !         IF(ALL(IGDTMPLI==IGDTMPLI_SAVE(1:IGDTLENI)))THEN
  !            SAME_GRIDI=.TRUE.
  !         ENDIF
  !      ENDIF
  !   ENDIF
  !   !
  !   IGDTNUMI_SAVE=IGDTNUMI
  !   IGDTLENI_SAVE=IGDTLENI
  !   IGDTMPLI_SAVE(1:IGDTLENI)=IGDTMPLI
  !   IGDTMPLI_SAVE(IGDTLENI+1:1000)=-9999
  !   !
  !   SAME_GRIDO=.FALSE.
  !   IF(IGDTNUMO==IGDTNUMO_SAVE)THEN
  !      IF(IGDTLENO==IGDTLENO_SAVE)THEN
  !         IF(ALL(IGDTMPLO==IGDTMPLO_SAVE(1:IGDTLENO)))THEN
  !            SAME_GRIDO=.TRUE.
  !         ENDIF
  !      ENDIF
  !   ENDIF
  !   !
  !   IGDTNUMO_SAVE=IGDTNUMO
  !   IGDTLENO_SAVE=IGDTLENO
  !   IGDTMPLO_SAVE(1:IGDTLENO)=IGDTMPLO
  !   IGDTMPLO_SAVE(IGDTLENO+1:1000)=-9999
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! END SUBROUTINE CHECK_GRIDS0V


  ! !> @file
  ! !! INTERPOLATE VECTOR FIELDS (BILINEAR)
  ! !! @author IREDELL @date 96-04-10
  ! !
  ! !> THIS SUBPROGRAM PERFORMS BILINEAR INTERPOLATION
  ! !!           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
  ! !!           OPTIONS ALLOW VARYING THE MINIMUM PERCENTAGE FOR MASK,
  ! !!           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
  ! !!           (IPOPT(1)) WHICH DEFAULTS TO 50 (IF IPOPT(1)=-1).
  ! !!           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
  ! !!           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
  ! !!           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
  ! !!           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
  ! !!             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
  ! !!             (KGDS(1)=001) MERCATOR CYLINDRICAL
  ! !!             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
  ! !!             (KGDS(1)=004) GAUSSIAN CYLINDRICAL (SPECTRAL NATIVE)
  ! !!             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
  ! !!             (KGDS(1)=203) ROTATED EQUIDISTANT CYLINDRICAL (E-STAGGER)
  ! !!             (KGDS(1)=205) ROTATED EQUIDISTANT CYLINDRICAL (B-STAGGER)
  ! !!           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
  ! !!           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
  ! !!           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
  ! !!           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
  ! !!           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
  ! !!           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
  ! !!           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
  ! !!           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
  ! !!           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  ! !!           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
  ! !!           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
  ! !!           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT 
  ! !!           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  ! !!           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
  ! !!           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
  ! !!           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
  ! !!           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
  ! !!        
  ! !! PROGRAM HISTORY LOG:
  ! !! -  96-04-10  IREDELL
  ! !! - 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
  ! !! - 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
  ! !! - 2002-01-17  IREDELL  SAVE DATA FROM LAST CALL FOR OPTIMIZATION
  ! !! - 2007-05-22  IREDELL  EXTRAPOLATE UP TO HALF A GRID CELL
  ! !! - 2007-10-30  IREDELL  SAVE WEIGHTS AND THREAD FOR PERFORMANCE
  ! !! - 2012-06-26  GAYNO    FIX OUT-OF-BOUNDS ERROR.  SEE NCEPLIBS
  ! !!                      TICKET #9.
  ! !! - 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH MERGED VERSION
  ! !!                      OF GDSWZD.
  ! !!
  ! !! @param IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
  ! !!                IPOPT(1) IS MINIMUM PERCENTAGE FOR MASK
  ! !!                (DEFAULTS TO 50 IF IPOPT(1)=-1)
  ! !! @param KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
  ! !! @param KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
  ! !!                (KGDSO(1)<0 IMPLIES RANDOM STATION POINTS)
  ! !! @param MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
  ! !!                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
  ! !! @param MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
  ! !!                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
  ! !! @param KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
  ! !! @param IBI      - INTEGER (KM) INPUT BITMAP FLAGS
  ! !! @param LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
  ! !! @param UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
  ! !! @param VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
  ! !! @param[out] NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)<0)
  ! !! @param[out] RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
  ! !! @param[out] RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
  ! !! @param[out] CROT     - REAL (NO) VECTOR ROTATION COSINES (IF KGDSO(1)<0)
  ! !! @param[out] SROT     - REAL (NO) VECTOR ROTATION SINES (IF KGDSO(1)<0)
  ! !!                (UGRID=CROT*UEARTH-SROT*VEARTH;
  ! !!                 VGRID=SROT*UEARTH+CROT*VEARTH)
  ! !! @param[out] IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
  ! !! @param[out] LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
  ! !! @param[out] UO       - REAL (MO,KM) OUTPUT U-COMPONENT FIELDS INTERPOLATED
  ! !! @param[out] VO       - REAL (MO,KM) OUTPUT V-COMPONENT FIELDS INTERPOLATED
  ! !! @param[out] IRET     - INTEGER RETURN CODE
  ! !!                0    SUCCESSFUL INTERPOLATION
  ! !!                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
  ! !!                3    UNRECOGNIZED OUTPUT GRID
  ! !!
  ! !! SUBPROGRAMS CALLED:
  ! !! -  GDSWZD       GRID DESCRIPTION SECTION WIZARD
  ! !! -  IJKGDS0      SET UP PARAMETERS FOR IJKGDS1
  ! !! -  (IJKGDS1)    RETURN FIELD POSITION FOR A GIVEN GRID POINT
  ! !! -  (MOVECT)     MOVE A VECTOR ALONG A GREAT CIRCLE
  ! !! -  POLFIXV      MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
  ! !!
  ! SUBROUTINE POLATEV0_grib1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI, &
  !      NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
  !   INTEGER,            INTENT(IN   ):: IPOPT(20),IBI(KM),MI,MO,KM
  !   INTEGER,            INTENT(IN   ):: KGDSI(200),KGDSO(200)
  !   INTEGER,            INTENT(INOUT):: NO
  !   INTEGER,            INTENT(  OUT):: IRET, IBO(KM)
  !   !
  !   LOGICAL*1,          INTENT(IN   ):: LI(MI,KM)
  !   LOGICAL*1,          INTENT(  OUT):: LO(MO,KM)
  !   !
  !   REAL,               INTENT(IN   ):: UI(MI,KM),VI(MI,KM)
  !   REAL,               INTENT(INOUT):: RLAT(MO),RLON(MO),CROT(MO),SROT(MO)
  !   REAL,               INTENT(  OUT):: UO(MO,KM),VO(MO,KM)
  !   !
  !   REAL,               PARAMETER    :: FILL=-9999.
  !   !
  !   INTEGER                          :: IJX(2),IJY(2),IJKGDSA(20)
  !   INTEGER                          :: MP,N,I,J,K,NK,NV,IJKGDS1
  !   INTEGER,                    SAVE :: KGDSIX(200)=-1,KGDSOX(200)=-1
  !   INTEGER,                    SAVE :: NOX=-1,IRETX=-1
  !   INTEGER,        ALLOCATABLE,SAVE :: NXY(:,:,:)
  !   !
  !   REAL                             :: CM,SM,UROT,VROT
  !   REAL,           ALLOCATABLE      :: DUM1(:),DUM2(:)
  !   REAL                             :: PMP,XIJ,YIJ,XF,YF,U,V,W
  !   REAL                             :: XPTS(MO),YPTS(MO)
  !   REAL                             :: WX(2),WY(2)
  !   REAL                             :: XPTI(MI),YPTI(MI)
  !   REAL                             :: RLOI(MI),RLAI(MI)
  !   REAL                             :: CROI(MI),SROI(MI)
  !   REAL,           ALLOCATABLE,SAVE :: RLATX(:),RLONX(:)
  !   REAL,           ALLOCATABLE,SAVE :: CROTX(:),SROTX(:)
  !   REAL,           ALLOCATABLE,SAVE :: WXY(:,:,:),CXY(:,:,:),SXY(:,:,:)
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   !  SET PARAMETERS
  !   IRET=0
  !   MP=IPOPT(1)
  !   IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
  !   IF(MP.LT.0.OR.MP.GT.100) IRET=32
  !   PMP=MP*0.01
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   !  SAVE OR SKIP WEIGHT COMPUTATION
  !   IF(IRET.EQ.0.AND.(KGDSO(1).LT.0.OR.ANY(KGDSI.NE.KGDSIX).OR.ANY(KGDSO.NE.KGDSOX))) THEN
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
  !      IF(KGDSO(1).GE.0) THEN
  !         CALL GDSWZD(KGDSO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO,CROT,SROT)
  !         IF(NO.EQ.0) IRET=3
  !      ENDIF
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  LOCATE INPUT POINTS
  !      ALLOCATE(DUM1(NO))
  !      ALLOCATE(DUM2(NO))
  !      CALL GDSWZD(KGDSI,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
  !      DEALLOCATE(DUM1,DUM2)
  !      IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
  !      CALL GDSWZD(KGDSI, 0,MI,FILL,XPTI,YPTI,RLOI,RLAI,NV,CROI,SROI)
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  ALLOCATE AND SAVE GRID DATA
  !      KGDSIX=KGDSI
  !      KGDSOX=KGDSO
  !      IF(NOX.NE.NO) THEN
  !         IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,CROTX,SROTX,NXY,WXY,CXY,SXY)
  !         ALLOCATE(RLATX(NO),RLONX(NO),CROTX(NO),SROTX(NO), &
  !              NXY(2,2,NO),WXY(2,2,NO),CXY(2,2,NO),SXY(2,2,NO))
  !         NOX=NO
  !      ENDIF
  !      IRETX=IRET
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !      !  COMPUTE WEIGHTS
  !      IF(IRET.EQ.0) THEN
  !         CALL IJKGDS0(KGDSI,IJKGDSA)
  !         !$OMP PARALLEL DO PRIVATE(N,XIJ,YIJ,IJX,IJY,XF,YF,J,I,WX,WY,CM,SM)
  !         DO N=1,NO
  !            RLONX(N)=RLON(N)
  !            RLATX(N)=RLAT(N)
  !            CROTX(N)=CROT(N)
  !            SROTX(N)=SROT(N)
  !            XIJ=XPTS(N)
  !            YIJ=YPTS(N)
  !            IF(XIJ.NE.FILL.AND.YIJ.NE.FILL) THEN
  !               IJX(1:2)=FLOOR(XIJ)+(/0,1/)
  !               IJY(1:2)=FLOOR(YIJ)+(/0,1/)
  !               XF=XIJ-IJX(1)
  !               YF=YIJ-IJY(1)
  !               WX(1)=(1-XF)
  !               WX(2)=XF
  !               WY(1)=(1-YF)
  !               WY(2)=YF
  !               DO J=1,2
  !                  DO I=1,2
  !                     NXY(I,J,N)=IJKGDS1(IJX(I),IJY(J),IJKGDSA)
  !                     WXY(I,J,N)=WX(I)*WY(J)
  !                     IF(NXY(I,J,N).GT.0) THEN
  !                        CALL MOVECT(RLAI(NXY(I,J,N)),RLOI(NXY(I,J,N)), &
  !                             RLAT(N),RLON(N),CM,SM)
  !                        CXY(I,J,N)=CM*CROI(NXY(I,J,N))+SM*SROI(NXY(I,J,N))
  !                        SXY(I,J,N)=SM*CROI(NXY(I,J,N))-CM*SROI(NXY(I,J,N))
  !                     ENDIF
  !                  ENDDO
  !               ENDDO
  !            ELSE
  !               NXY(:,:,N)=0
  !            ENDIF
  !         ENDDO
  !      ENDIF  ! IS IRET 0?
  !   ENDIF
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   !  INTERPOLATE OVER ALL FIELDS
  !   IF(IRET.EQ.0.AND.IRETX.EQ.0) THEN
  !      IF(KGDSO(1).GE.0) THEN
  !         NO=NOX
  !         DO N=1,NO
  !            RLON(N)=RLONX(N)
  !            RLAT(N)=RLATX(N)
  !            CROT(N)=CROTX(N)
  !            SROT(N)=SROTX(N)
  !         ENDDO
  !      ENDIF
  !      !$OMP PARALLEL DO PRIVATE(NK,K,N,U,V,W,UROT,VROT,J,I)
  !      DO NK=1,NO*KM
  !         K=(NK-1)/NO+1
  !         N=NK-NO*(K-1)
  !         U=0
  !         V=0
  !         W=0
  !         DO J=1,2
  !            DO I=1,2
  !               IF(NXY(I,J,N).GT.0) THEN
  !                  IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K)) THEN
  !                     UROT=CXY(I,J,N)*UI(NXY(I,J,N),K)-SXY(I,J,N)*VI(NXY(I,J,N),K)
  !                     VROT=SXY(I,J,N)*UI(NXY(I,J,N),K)+CXY(I,J,N)*VI(NXY(I,J,N),K)
  !                     U=U+WXY(I,J,N)*UROT
  !                     V=V+WXY(I,J,N)*VROT
  !                     W=W+WXY(I,J,N)
  !                  ENDIF
  !               ENDIF
  !            ENDDO
  !         ENDDO
  !         LO(N,K)=W.GE.PMP
  !         IF(LO(N,K)) THEN
  !            UROT=CROT(N)*U-SROT(N)*V
  !            VROT=SROT(N)*U+CROT(N)*V
  !            UO(N,K)=UROT/W
  !            VO(N,K)=VROT/W
  !         ELSE
  !            UO(N,K)=0.
  !            VO(N,K)=0.
  !         ENDIF
  !      ENDDO  ! NK LOOP
  !      DO K=1,KM
  !         IBO(K)=IBI(K)
  !         IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
  !      ENDDO
  !      IF(KGDSO(1).EQ.0) CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
  !      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   ELSE
  !      IF(IRET.EQ.0) IRET=IRETX
  !      IF(KGDSO(1).GE.0) NO=0
  !   ENDIF
  !   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! END SUBROUTINE POLATEV0_GRIB1



end module polatev0_mod
