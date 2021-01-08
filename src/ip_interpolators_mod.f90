module ip_interpolators_mod
  ! re-export specific interpolation routines
  use bilinear_interpolator_scalar_mod
  use bicubic_interpolator_scalar_mod
  use neighbor_interpolator_scalar_mod
  use budget_interpolator_scalar_mod
  use spectral_interpolator_scalar_mod
  use neighbor_budget_interpolator_scalar_mod

  implicit none
end module ip_interpolators_mod

