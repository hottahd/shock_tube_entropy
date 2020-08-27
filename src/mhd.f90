!=====================================================|
subroutine mhd
!=====================================================|
  use geometry_def, only: nxg,qq,qqm,mtype
  implicit none

  integer :: i,m
!-----------------------------------------------------|

  call cfl
  call sc4rk4
  call artdif

  

  do i = 1,nxg
  do m = 1,mtype
     qq(i,m) = qqm(i,m)
  enddo
  enddo

  return
end subroutine mhd
