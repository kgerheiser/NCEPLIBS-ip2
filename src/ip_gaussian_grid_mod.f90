module ip_gaussian_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  implicit none

  private
  public :: ip_gaussian_grid

  type, extends(ip_grid) :: ip_gaussian_grid
     integer :: jh
     real :: dlon, rlat1, rlon1, rlon2, hi
     integer :: jg, jscan
   contains
     procedure :: init_grib1
     procedure :: init_grib2
  end type ip_gaussian_grid

contains

  subroutine init_grib1(self, g1_desc)
    class(ip_gaussian_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6
      self%eccen_squared = 0.0

      self%IM=KGDS(2)
      self%JM=KGDS(3)
      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      self%JG=KGDS(10)*2
      ISCAN=MOD(KGDS(11)/128,2)
      self%JSCAN=MOD(KGDS(11)/64,2)
      self%HI=(-1.)**ISCAN
      self%JH=(-1)**self%JSCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
    end associate
  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_gaussian_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)
      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6
      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      self%RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
      self%JG=IGDTMPL(18)*2
      ISCAN=MOD(IGDTMPL(19)/128,2)
      self%JSCAN=MOD(IGDTMPL(19)/64,2)
      self%HI=(-1.)**ISCAN
      self%JH=(-1)**self%JSCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
    end associate

  end subroutine init_grib2

end module ip_gaussian_grid_mod

