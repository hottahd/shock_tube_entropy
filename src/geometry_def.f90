!=====================================================|
module geometry_def
!=====================================================|
  integer, parameter :: nx = 400
  integer, parameter :: margin = 2
  integer, parameter :: nxg = nx + 2*margin
  integer, parameter :: mtype = 8

  ! geometry
  real(KIND(0.d0)), save :: xmax, xmin, dx
  real(KIND(0.d0)), dimension(nxg), save :: x  ! geometry
  real(KIND(0.d0)), dimension(nxg,mtype), save :: qq, qqm, qqp ! variable
  real(KIND(0.d0)), dimension(nxg), save :: cc ! characteristic velocity

  ! time variable
  integer, save :: ns, nd
  real(KIND(0.d0)), save :: time,dt

  ! data directory
  character (len = 10),save :: data_dir
end module geometry_def
 
