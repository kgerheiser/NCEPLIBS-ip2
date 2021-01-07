module ip_polar_stereo_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use constants_mod, only: DPR, PI
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
  end type ip_polar_stereo_grid

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
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
    end associate
  end subroutine init_grib2


end module ip_polar_stereo_grid_mod

