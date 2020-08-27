!=====================================================|
subroutine sc4rk4
!=====================================================|
  use geometry_def, only: qq,qqm,qqp,dx,dt,nxg,margin,mtype
  use const_def, only: gm,cv,pii4,pii8
  implicit none

  integer :: i,m,mm
  real(KIND(0.d0)) :: ratio
  ! total pressure (gas + magnetic)
  real(KIND(0.d0)), dimension(nxg) :: pt 
  real(KIND(0.d0)) :: dro,drx,dry,drz,dbx,dby,dbz,dse
!-----------------------------------------------------|  

  ! Runge-Kutta loop
  do mm = 1,4
     if(mm == 1) ratio = 1.d0/4.d0
     if(mm == 2) ratio = 1.d0/3.d0
     if(mm == 3) ratio = 1.d0/2.d0
     if(mm == 4) ratio = 1.d0

     if(mm == 1) then
        do i = 1,nxg
        do m = 1,mtype
           qqp(i,m) = qq(i,m)
        enddo
        enddo
     else
        do i = 1,nxg
        do m = 1,mtype
           qqp(i,m) = qqm(i,m)
        enddo
        enddo
     endif

     do i = 1,nxg
        pt(i) = qqp(i,1)**gm*exp(qqp(i,8)/cv) &
        & + (qqp(i,5)**2 + qqp(i,6)**2 + qqp(i,7)**2)*pii8
     enddo
     
     do i = 1+margin,nxg-margin
        dro = -( &
             & -      qqp(i+2,1)*qqp(i+2,2) &
             & + 8.d0*qqp(i+1,1)*qqp(i+1,2) &
             & - 8.d0*qqp(i-1,1)*qqp(i-1,2) &
             & +      qqp(i-2,1)*qqp(i-2,2) &
             & )/dx/12.d0*dt
        
        drx = -( &
             & -      qqp(i+2,1)*qqp(i+2,2)*qqp(i+2,2) &
             & + 8.d0*qqp(i+1,1)*qqp(i+1,2)*qqp(i+1,2) &
             & - 8.d0*qqp(i-1,1)*qqp(i-1,2)*qqp(i-1,2) &
             & +      qqp(i-2,1)*qqp(i-2,2)*qqp(i-2,2) &
             & -      pt(i+2) &
             & + 8.d0*pt(i+1) &
             & - 8.d0*pt(i-1) &
             & +      pt(i-2) &
             & - pii4*( &
             & -      qqp(i+2,5)*qqp(i+2,5) &
             & + 8.d0*qqp(i+1,5)*qqp(i+1,5) &
             & - 8.d0*qqp(i-1,5)*qqp(i-1,5) &
             & +      qqp(i-2,5)*qqp(i-2,5) &
             & ) &
             & )/dx/12.d0*dt

        
        dry = -( &
             & -      qqp(i+2,1)*qqp(i+2,2)*qqp(i+2,3) &
             & + 8.d0*qqp(i+1,1)*qqp(i+1,2)*qqp(i+1,3) &
             & - 8.d0*qqp(i-1,1)*qqp(i-1,2)*qqp(i-1,3) &
             & +      qqp(i-2,1)*qqp(i-2,2)*qqp(i-2,3) &
             & - pii4*( &
             & -      qqp(i+2,5)*qqp(i+2,6) &
             & + 8.d0*qqp(i+1,5)*qqp(i+1,6) &
             & - 8.d0*qqp(i-1,5)*qqp(i-1,6) &
             & +      qqp(i-2,5)*qqp(i-2,6) &
             & ) &
             & )/dx/12.d0*dt

        drz = -( &
             & -      qqp(i+2,1)*qqp(i+2,2)*qqp(i+2,4) &
             & + 8.d0*qqp(i+1,1)*qqp(i+1,2)*qqp(i+1,4) &
             & - 8.d0*qqp(i-1,1)*qqp(i-1,2)*qqp(i-1,4) &
             & +      qqp(i-2,1)*qqp(i-2,2)*qqp(i-2,4) &
             & - pii4*( & 
             & -      qqp(i+2,5)*qqp(i+2,7) &
             & + 8.d0*qqp(i+1,5)*qqp(i+1,7) &
             & - 8.d0*qqp(i-1,5)*qqp(i-1,7) &
             & +      qqp(i-2,5)*qqp(i-2,7) &
             & ) &
             & )/dx/12.d0*dt
        
        dbx = 0.d0
        
        dby = + ( &
             & - ( &
             & -      qqp(i+2,2)*qqp(i+2,6) &
             & + 8.d0*qqp(i+1,2)*qqp(i+1,6) &
             & - 8.d0*qqp(i-1,2)*qqp(i-1,6) &
             & +      qqp(i-2,2)*qqp(i-2,6) &
             & ) &
             & + ( &
             & -      qqp(i+2,3)*qqp(i+2,5) &
             & + 8.d0*qqp(i+1,3)*qqp(i+1,5) &
             & - 8.d0*qqp(i-1,3)*qqp(i-1,5) &
             & +      qqp(i-2,3)*qqp(i-2,5) &
             & ) &
             & )/dx/12.d0*dt

        dbz = + ( &
             & - ( &
             & -      qqp(i+2,2)*qqp(i+2,7) &
             & + 8.d0*qqp(i+1,2)*qqp(i+1,7) &
             & - 8.d0*qqp(i-1,2)*qqp(i-1,7) &
             & +      qqp(i-2,2)*qqp(i-2,7) &
             & ) &
             & + ( &
             & -      qqp(i+2,4)*qqp(i+2,5) &
             & + 8.d0*qqp(i+1,4)*qqp(i+1,5) &
             & - 8.d0*qqp(i-1,4)*qqp(i-1,5) &
             & +      qqp(i-2,4)*qqp(i-2,5) &
             & ) &
             & )/dx/12.d0*dt
        
        dse = -qqp(i,2)*( &
             & -      qqp(i+2,8) &
             & + 8.d0*qqp(i+1,8) &
             & - 8.d0*qqp(i-1,8) &
             & +      qqp(i-2,8) &
             & )/dx/12.d0*dt

        
        qqm(i,1) = qq(i,1) + ratio*dro
        qqm(i,2) = (qq(i,1)*qq(i,2) + ratio*drx)/qqm(i,1)
        qqm(i,3) = (qq(i,1)*qq(i,3) + ratio*dry)/qqm(i,1)
        qqm(i,4) = (qq(i,1)*qq(i,4) + ratio*drz)/qqm(i,1)
        qqm(i,5) = qq(i,5) + ratio*dbx
        qqm(i,6) = qq(i,6) + ratio*dby
        qqm(i,7) = qq(i,7) + ratio*dbz
        qqm(i,8) = qq(i,8) + ratio*dse
     enddo
     call bc(qqm)
  enddo
  
  return
end subroutine sc4rk4
