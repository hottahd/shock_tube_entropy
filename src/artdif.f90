!=====================================================|
subroutine artdif
!=====================================================|
  use geometry_def, only: qqm,qqp,nxg,dx,dt,mtype,cc
  use const_def, only: cv,cp,gm,rr,pii4
  implicit none

  integer :: i,m
  real(KIND(0.d0)), parameter :: fh = 0.d0, ep = 1.d0
  real(KIND(0.d0)), dimension(nxg,mtype) :: fx,dui
  real(KIND(0.d0)) :: dul,dur,duc,mup,mlo
  real(KIND(0.d0)) :: ul,ur,du,ra,pp
  real(KIND(0.d0)) :: dro,drx,dry,drz,dbx,dby,dbz,dse

  real(KIND(0.d0)), dimension(nxg) :: pr,te,hh ! pressure, temperature, enthalpy
!-----------------------------------------------------|  

  do i = 1,nxg
     qqp(i,1) = qqm(i,1)
     qqp(i,2) = qqm(i,2)
     qqp(i,3) = qqm(i,3)
     qqp(i,4) = qqm(i,4)
     qqp(i,5) = qqm(i,5)
     qqp(i,6) = qqm(i,6)
     qqp(i,7) = qqm(i,7)
     qqp(i,8) = qqm(i,8)

     pr(i) = qqp(i,1)**gm*exp(qqp(i,8)/cv)
     te(i) = pr(i)/rr/qqp(i,1)
     hh(i) = cp*te(i) + 0.5d0*(qqp(i,2)**2 + qqp(i,3)**2 + qqp(i,4)**2)
  enddo


  ! calculating artificial diffusivity flux
  do i = 2,nxg-1
  do m = 1,mtype
     dul = qqp(i  ,m) - qqp(i-1,m)
     dur = qqp(i+1,m) - qqp(i  ,m)
     duc = 0.5d0*(qqp(i+1,m) - qqp(i-1,m))

     mup = max(ep*dul,ep*dur,duc)
     mlo = min(ep*dul,ep*dur,duc)

     dui(i,m) = min(0.d0,mup) + max(0.d0,mlo)
  enddo
  enddo

  do i = 2,nxg-2
  do m = 1,mtype
     ul = qqp(i  ,m) + 0.5d0*dui(i  ,m)
     ur = qqp(i+1,m) - 0.5d0*dui(i+1,m)
     du = qqp(i+1,m) - qqp(i,m)
     du = sign(1.d0,du)*max(abs(du),1.d-10)
     ra = min((ur-ul)/du,1.d0)
     pp = (0.5d0 + sign(0.5d0,ra))*max(0.d0,1.d0 + fh*(ra-1.d0))
     fx(i,m) = -0.5d0*(ur - ul)*pp*max(cc(i),cc(i+1))
     !fx(i,m) = -0.5d0*(ur - ul)*pp*0.5d0*(cc(i)+cc(i+1))
  enddo  
  enddo

  do i = 2,nxg-2
     fx(i,8) = &
          & + 0.5d0*(qqp(i,1)*te(i) + qqp(i+1,1)*te(i+1))*fx(i,8) &
          & + 0.5d0*(qqp(i,1)*qqp(i,2) + qqp(i+1,1)*qqp(i+1,2))*fx(i,2) &
          & + 0.5d0*(qqp(i,1)*qqp(i,3) + qqp(i+1,1)*qqp(i+1,3))*fx(i,3) &
          & + 0.5d0*(qqp(i,1)*qqp(i,4) + qqp(i+1,1)*qqp(i+1,4))*fx(i,4) &
          & + 0.5d0*(qqp(i,5) + qqp(i+1,5))*fx(i,5)*pii4 &
          & + 0.5d0*(qqp(i,6) + qqp(i+1,6))*fx(i,6)*pii4 &
          & + 0.5d0*(qqp(i,7) + qqp(i+1,7))*fx(i,7)*pii4 
     
     fx(i,2) = + 0.5d0*(qqp(i,1) + qqp(i+1,1))*fx(i,2) &
          &    + 0.5d0*(qqp(i,2) + qqp(i+1,2))*fx(i,1)

     fx(i,3) = + 0.5d0*(qqp(i,1) + qqp(i+1,1))*fx(i,3) &
          &    + 0.5d0*(qqp(i,3) + qqp(i+1,3))*fx(i,1)

     fx(i,4) = + 0.5d0*(qqp(i,1) + qqp(i+1,1))*fx(i,4) &
          &    + 0.5d0*(qqp(i,4) + qqp(i+1,4))*fx(i,1)
  enddo
                    
  do i = 3,nxg-2
     dro = -(fx(i,1) - fx(i-1,1))/dx*dt
     drx = -(fx(i,2) - fx(i-1,2))/dx*dt
     dry = -(fx(i,3) - fx(i-1,3))/dx*dt
     drz = -(fx(i,4) - fx(i-1,4))/dx*dt
     dbx = -(fx(i,5) - fx(i-1,5))/dx*dt
     dby = -(fx(i,6) - fx(i-1,6))/dx*dt
     dbz = -(fx(i,7) - fx(i-1,7))/dx*dt
     dse =  ( &
          & - (fx(i,8) - fx(i-1,8))/dx*dt &
          & - (hh(i) - qqp(i,2)**2 - qqp(i,3)**2 - qqp(i,4)**2)*dro &
          & - (qqp(i,2)*drx + qqp(i,3)*dry + qqp(i,4)*drz) &
          & - (qqp(i,5)*dbx + qqp(i,6)*dby + qqp(i,7)*dbz)*pii4 &
          & )
     qqm(i,1) = qqm(i,1) + dro
     qqm(i,2) = (qqp(i,1)*qqm(i,2) + drx)/qqm(i,1)
     qqm(i,3) = (qqp(i,1)*qqm(i,3) + dry)/qqm(i,1)
     qqm(i,4) = (qqp(i,1)*qqm(i,4) + drz)/qqm(i,1)
     qqm(i,5) = qqm(i,5) + dbx
     qqm(i,6) = qqm(i,6) + dby
     qqm(i,7) = qqm(i,7) + dbz
     qqm(i,8) = qqm(i,8) + dse/qqp(i,1)/te(i)
  enddo

  call bc(qqm)
  return
end subroutine artdif
