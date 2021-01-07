module ip_equid_cylind_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  implicit none

  private
  public :: ip_equid_cylind_grid

  type, extends(ip_grid) :: ip_equid_cylind_grid
     real :: hi, rlat1, rlon1, rlat2, rlon2
     real :: dlat, dlon
   contains
     procedure :: init_grib1
     procedure :: init_grib2
  end type ip_equid_cylind_grid

contains

  subroutine init_grib1(self, g1_desc)
    class(ip_equid_cylind_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan

    associate(kgds => g1_desc%gds)
      self%IM=KGDS(2)
      self%JM=KGDS(3)
      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLAT2=KGDS(7)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      ISCAN=MOD(KGDS(11)/128,2)
      self%HI=(-1.)**ISCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DLAT=(self%RLAT2-self%RLAT1)/(self%JM-1)

      ! defaults
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%kscan = 0

      self%iwrap = nint(360/abs(self%dlon))

      if(self%im < self%iwrap) self%iwrap=0
      self%jwrap1 = 0
      self%jwrap2 = 0
      if(self%iwrap > 0 .and. mod(self%iwrap,2) == 0) then
         if(abs(self%rlat1) > 90-0.25*self%dlat) then
            self%jwrap1 = 2
         elseif(abs(self%rlat1) > 90-0.75*self%dlat) then
            self%jwrap1 = 1
         endif
         if(abs(self%rlat2) > 90-0.25*self%dlat) then
            self%jwrap2 = 2 * self%jm
         elseif(abs(self%rlat2) > 90-0.75*self%dlat) then
            self%jwrap2 = 2 * self%jm+1
         endif
      endif

      self%rerth = 6.3712E6
      self%eccen_squared = 0.0
    end associate

  end subroutine init_grib1
  

  subroutine init_grib2(self, g2_desc)
    class(ip_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)
      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6
      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      self%RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
      self%RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
      ISCAN=MOD(IGDTMPL(19)/128,2)
      self%HI=(-1.)**ISCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DLAT=(self%RLAT2-self%RLAT1)/(self%JM-1)

      self%nscan = MOD(IGDTMPL(19)/32,2)
      self%kscan = 0
      self%iwrap = NINT(360/ABS(self%DLON))

      if(self%im.lt.self%iwrap) self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0

      if(self%im < self%iwrap) self%iwrap=0
      self%jwrap1 = 0
      self%jwrap2 = 0
      if(self%iwrap > 0 .and. mod(self%iwrap,2) == 0) then
         if(abs(self%rlat1) > 90-0.25*self%dlat) then
            self%jwrap1 = 2
         elseif(abs(self%rlat1) > 90-0.75*self%dlat) then
            self%jwrap1 = 1
         endif
         if(abs(self%rlat2) > 90-0.25*self%dlat) then
            self%jwrap2 = 2 * self%jm
         elseif(abs(self%rlat2) > 90-0.75*self%dlat) then
            self%jwrap2 = 2 * self%jm+1
         endif
      endif

      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

    end associate
  end subroutine init_grib2

end module ip_equid_cylind_grid_mod

