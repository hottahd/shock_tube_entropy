!=====================================================|
subroutine cfl
!=====================================================|
  use geometry_def, only: dt,dx,margin,nxg,qq,cc
  use const_def, only: cv,gm,pi
  implicit none
  
  integer :: i
  real(KIND(0.d0)) :: pr ! gas pressure
  real(KIND(0.d0)) :: cs ! speed of sound
  real(KIND(0.d0)) :: vv ! total fluid velocity
  real(KIND(0.d0)) :: ca ! Alfven velocity
  real(KIND(0.d0)) :: safety = 0.5d0
!-----------------------------------------------------|  
  dt = 1.d20
  
  do i = 1,nxg
     pr = qq(i,1)**gm*exp(qq(i,8)/cv)
     cs = sqrt(gm*pr/qq(i,1))
     vv = sqrt( qq(i,2)**2 + qq(i,3)**2 + qq(i,4)**2)
     ca = sqrt((qq(i,5)**2 + qq(i,6)**2 + qq(i,7)**2)/4.d0/pi/qq(i,1))

     cc(i) = cs + vv + ca
     dt = min(safety*dx/cc(i),dt)
     cc(i) = 1.d0*cc(i)
  enddo

  return
end subroutine cfl
