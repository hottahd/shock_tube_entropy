!=====================================================|
subroutine model
!=====================================================|
  use geometry_def, only: xmax,xmin,dx,x,qq &
       & ,nx,nxg,margin
  use const_def, only: cv,gm,pi
  implicit none
  integer :: i

  real(KIND(0.d0)) :: xm,dd,tmp
  real(KIND(0.d0)) :: pr ! pressure
  real(KIND(0.d0)) :: rol,ror,vxl,vxr,prl,prr,byl,byr
!-----------------------------------------------------|

  xmax = 1.d0
  xmin = 0.d0

  dx = (xmax - xmin)/real(nx)
  x(1 + margin) = xmin + 0.5d0*dx
  do i = margin,1,-1
     x(i) = x(i+1) - dx
  enddo

  do i = 2 + margin,nxg
     x(i) = x(i-1) + dx
  enddo

  ! sound wave
  xm = 0.5d0*(xmax + xmin)
  dd = 0.05d0
  do i = 1,nxg
     qq(i,1) = 1.d0 + 0.01d0*exp(-((x(i)-xm)/dd)**2) ! density
     qq(i,2) = 0.d0 ! vx
     qq(i,3) = 0.d0 ! vy
     qq(i,4) = 0.d0 ! vz
     qq(i,5) = 1.d0 ! bx
     qq(i,6) = 1.d0 ! by
     qq(i,7) = 1.d0 ! bz
     pr = qq(i,1)**gm/gm
     qq(i,8) = cv*log(pr/qq(i,1)**gm) ! entropy
  enddo
  
  rol = 1.d0
  ror = 0.125d0
  vxl = 0.75d0
  vxr = 0.d0
  prl = 1.d0
  prr = 0.1d0
  ! hydrodynamic shock tube
  do i = 1,nxg
!     if( x(i) < 0.5d0) then
!        qq(i,1) = 1.d0
!        pr = 1.d0
!        qq(i,2) = 0.75d0
!     else
!        qq(i,1) = 0.125d0
!        pr = 0.1d0
!        qq(i,2) = 0d0
!     endif
     tmp = 0.5d0 + 0.5d0*tanh((x(i) - 0.5)/(1.d0*dx))
     qq(i,1) = rol + (ror - rol)*tmp
     qq(i,2) = vxl + (vxr - vxl)*tmp
     pr      = prl + (prr - prl)*tmp
     
     qq(i,3) = 0.d0
     qq(i,4) = 0.d0
     qq(i,5) = 0.d0
     qq(i,6) = 0.d0
     qq(i,7) = 0.d0
     qq(i,8) = cv*log(pr/qq(i,1)**gm) ! entropy
  enddo

!  goto 1000
  
  ! mhd shock tube (Brio and Wu, 1988)
  rol = 1.d0
  ror = 0.125d0
  prl = 1.d0
  prr = 0.1d0
  byl = 1.d0*sqrt(4.d0*pi)
  byr = -1.d0*sqrt(4.d0*pi)
  do i = 1,nxg
!     if( x(i) < 0.5d0) then
!        qq(i,1) = 1.d0
!        pr = 1.d0
!        qq(i,6) = 1.d0*sqrt(4.d0*pi) ! by
!     else
!        qq(i,1) = 0.125d0
!        pr = 0.1d0
!        qq(i,6) = -1.d0*sqrt(4.d0*pi) ! by
!     endif
     tmp = 0.5d0 + 0.5d0*tanh((x(i) - 0.5)/(1.d0*dx))
     qq(i,1) = rol + (ror - rol)*tmp
     pr      = prl + (prr - prl)*tmp
     qq(i,6) = byl + (byr - byl)*tmp

     
     qq(i,2) = 0.d0 ! vx
     qq(i,3) = 0.d0 ! vy
     qq(i,4) = 0.d0 ! vy
     qq(i,5) = 0.75d0*sqrt(4.d0*pi)  !bx
     qq(i,7) = 0.d0  !bz
     qq(i,8) = cv*log(pr/qq(i,1)**gm) ! entropy
  enddo

1000 continue
  
  return
end subroutine model
