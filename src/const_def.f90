module const_def

  ! ratio of heat capacities
  !real(KIND(0.d0)), parameter :: gm = 5.d0/3.d0
  !real(KIND(0.d0)), parameter :: gm = 1.4d0
  real(KIND(0.d0)), parameter :: gm = 2.d0
  real(KIND(0.d0)), parameter :: rr = 1.d0 ! gas constant
  ! heat capacity at constant volume
  real(KIND(0.d0)), parameter :: cv = rr/(gm-1.d0)
  ! heat capacity at constant pressure
  real(KIND(0.d0)), parameter :: cp = gm*cv
  real(KIND(0.d0)), parameter :: pi = 3.14159265359d0
  real(KIND(0.d0)), parameter :: pii4 = 1.d0/4.d0/pi
  real(KIND(0.d0)), parameter :: pii8 = 1.d0/8.d0/pi

end module const_def
