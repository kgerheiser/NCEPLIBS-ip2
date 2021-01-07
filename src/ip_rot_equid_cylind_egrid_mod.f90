module ip_rot_equid_cylind_egrid_mod
  use iso_fortran_env, only: real64
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_rot_equid_cylind_egrid

  integer, parameter :: kd = real64

  type, extends(ip_grid) :: ip_rot_equid_cylind_egrid
     real(kd) :: rlon0, rlon1, rlat1, clat0, slat0
     real(kd) :: dlats, dlons, hi
     integer :: irot
   contains
     procedure :: init_grib1
     procedure :: init_grib2
  end type ip_rot_equid_cylind_egrid


contains

  subroutine init_grib1(self, g1_desc)
    class(ip_rot_equid_cylind_egrid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscale, iscan
    real(kd) :: rlat0

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6_KD
      self%eccen_squared = 0.0

      self%IM=KGDS(2)*2-1
      self%JM=KGDS(3)

      self%RLAT1=KGDS(4)*1.E-3_KD
      self%RLON1=KGDS(5)*1.E-3_KD
      RLAT0=KGDS(7)*1.E-3_KD
      self%RLON0=KGDS(8)*1.E-3_KD

      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      self%IROT=MOD(KGDS(6)/8,2)
      self%KSCAN=MOD(KGDS(11)/256,2)
      ISCAN=MOD(KGDS(11)/128,2)

      self%HI=(-1.)**ISCAN      
    end associate

  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_rot_equid_cylind_egrid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan
    real(kd) :: rlat0
    integer :: i_offset_odd!, i_offset_even

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,self%rerth,self%eccen_squared)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ROUTINE ONLY WORKS FOR "E"-STAGGER GRIDS.
      !   "V" GRID WHEN BIT 5 IS '1' AND BIT 6 IS '0'.
      !   "H" GRID WHEN BIT 5 IS '0' AND BIT 6 IS '1'.
      ! I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
      ! I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
      ! IF(I_OFFSET_ODD==I_OFFSET_EVEN) THEN
      !    CALL ROT_EQUID_CYLIND_EGRID_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
      !    RETURN
      ! ENDIF
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      self%IM=IGDTMPL(8)*2-1
      self%JM=IGDTMPL(9)

      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6

      self%RLON0=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)
      self%DLATS=FLOAT(IGDTMPL(18))/FLOAT(ISCALE)
      ! THE GRIB2 CONVENTION FOR "I" RESOLUTION IS TWICE WHAT THIS ROUTINE ASSUMES.
      self%DLONS=FLOAT(IGDTMPL(17))/FLOAT(ISCALE) * 0.5_KD

      self%IROT=MOD(IGDTMPL(14)/8,2)

      I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
      self%KSCAN=I_OFFSET_ODD
      ISCAN=MOD(IGDTMPL(19)/128,2)

      self%HI=(-1.)**ISCAN

      RLAT0=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
      RLAT0=RLAT0+90.0_KD

      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
    end associate
  end subroutine init_grib2
end module ip_rot_equid_cylind_egrid_mod


