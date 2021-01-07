module ip_station_points_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  implicit none

  ! Not really a grid
  private
  public :: ip_station_points_grid

  type, extends(ip_grid) :: ip_station_points_grid
   contains
     procedure :: init_grib1
     procedure :: init_grib2
  end type ip_station_points_grid

contains

  subroutine init_grib1(self, g1_desc)
    class(ip_station_points_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    
  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_station_points_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc
  end subroutine init_grib2
  
end module ip_station_points_grid_mod


  
