module ip_rot_equid_cylind_grid_mod
  use iso_fortran_env, only: real64
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_rot_equid_cylind_grid

  integer, parameter :: kd = real64

  type, extends(ip_grid) :: ip_rot_equid_cylind_grid
     real(kd) :: clat0, dlats, dlons, rlon0, slat0, wbd, sbd
     integer :: irot
   contains
     procedure :: init_grib1
     procedure :: init_grib2
  end type ip_rot_equid_cylind_grid


CONTAINS

  subroutine init_grib1(self, g1_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
    real(kd) :: hs, hs2, slat1, slat2, slatr, clon1, clon2, clat1, clat2, clatr, clonr, rlonr, rlatr
    integer :: iscale

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6_KD
      self%eccen_squared = 0d0

      RLAT1=KGDS(4)*1.E-3_KD
      RLON1=KGDS(5)*1.E-3_KD
      RLAT0=KGDS(7)*1.E-3_KD
      self%RLON0=KGDS(8)*1.E-3_KD
      RLAT2=KGDS(12)*1.E-3_KD
      RLON2=KGDS(13)*1.E-3_KD

      self%IROT=MOD(KGDS(6)/8,2)
      self%IM=KGDS(2)
      self%JM=KGDS(3)

      SLAT1=SIN(RLAT1/DPR)
      CLAT1=COS(RLAT1/DPR)
      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      HS=SIGN(1._KD,MOD(RLON1-self%RLON0+180+3600,360._KD)-180)
      CLON1=COS((RLON1-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT1-self%SLAT0*CLAT1*CLON1
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT1*CLON1+self%SLAT0*SLAT1)/CLATR
      RLATR=DPR*ASIN(SLATR)
      RLONR=HS*DPR*ACOS(CLONR)

      self%WBD=RLONR
      self%SBD=RLATR
      SLAT2=SIN(RLAT2/DPR)
      CLAT2=COS(RLAT2/DPR)
      HS2=SIGN(1._KD,MOD(RLON2-self%RLON0+180+3600,360._KD)-180)
      CLON2=COS((RLON2-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT2-self%SLAT0*CLAT2*CLON2
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT2*CLON2+self%SLAT0*SLAT2)/CLATR
      NBD=DPR*ASIN(SLATR)
      EBD=HS2*DPR*ACOS(CLONR)
      self%DLATS=(NBD-self%SBD)/FLOAT(self%JM-1)
      self%DLONS=(EBD-self%WBD)/FLOAT(self%IM-1)
    end associate

  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
    integer :: iscale
    integer :: i_offset_odd, i_offset_even, j_offset

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)

      CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,self%rerth,self%eccen_squared)

      I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
      I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
      J_OFFSET=MOD(IGDTMPL(19)/2,2)

      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6

      RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      RLAT0=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
      RLAT0=RLAT0+90.0_KD

      self%RLON0=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)

      RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
      RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)

      self%IROT=MOD(IGDTMPL(14)/8,2)
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      self%WBD=RLON1
      IF (self%WBD > 180.0) self%WBD = self%WBD - 360.0
      self%SBD=RLAT1

      NBD=RLAT2
      EBD=RLON2

      self%DLATS=(NBD-self%SBD)/FLOAT(self%JM-1)
      self%DLONS=(EBD-self%WBD)/FLOAT(self%IM-1)

      IF(I_OFFSET_ODD==1) self%WBD=self%WBD+(0.5_KD*self%DLONS)
      IF(J_OFFSET==1) self%SBD=self%SBD+(0.5_KD*self%DLATS)

    end associate
  end subroutine init_grib2

end module ip_rot_equid_cylind_grid_mod

