module ip_mercator_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_mercator_grid

  type, extends(ip_grid) :: ip_mercator_grid
     real :: rlat1, rlon1, rlon2, rlati, hi, dlon, dphi
   contains
     procedure :: init_grib1
     procedure :: init_grib2
  end type ip_mercator_grid

CONTAINS

  subroutine init_grib1(self, g1_desc)
    class(ip_mercator_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan, jscan
    real :: dy, hj

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6
      self%eccen_squared = 0.0

      self%IM=KGDS(2)
      self%JM=KGDS(3)

      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      self%RLATI=KGDS(9)*1.E-3

      ISCAN=MOD(KGDS(11)/128,2)
      JSCAN=MOD(KGDS(11)/64,2)

      DY=KGDS(13)
      self%HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DPHI=HJ*DY/(self%RERTH*COS(self%RLATI/DPR))

      ! defaults
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%kscan = 0

      self%iwrap = nint(360 / abs(self%dlon))
      if (self%im < self%iwrap) self%iwrap = 0
    end associate

  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_mercator_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscan, jscan
    real :: hj, dy

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)

      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%RLAT1=FLOAT(IGDTMPL(10))*1.0E-6
      self%RLON1=FLOAT(IGDTMPL(11))*1.0E-6
      self%RLON2=FLOAT(IGDTMPL(15))*1.0E-6
      self%RLATI=FLOAT(IGDTMPL(13))*1.0E-6

      ISCAN=MOD(IGDTMPL(16)/128,2)
      JSCAN=MOD(IGDTMPL(16)/64,2)

      DY=FLOAT(IGDTMPL(19))*1.0E-3
      self%HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DPHI=HJ*DY/(self%RERTH*COS(self%RLATI/DPR))

      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
      self%nscan=mod(igdtmpl(16) / 32,2)

      self%iwrap = nint(360 / abs(self%dlon))
      if(self%im < self%iwrap) self%iwrap = 0

    end associate
  end subroutine init_grib2

end module ip_mercator_grid_mod

