!=====================================================|
subroutine bc(qq)
!=====================================================|
  use geometry_def, only: nxg,margin,mtype
  implicit none

  integer :: i,m
  real(KIND(0.d0)), dimension(nxg,mtype), intent(inout) :: qq
!-----------------------------------------------------|  

  do i = 1,margin
  do m = 1,mtype
     ! periodic boundary
     !qq(i,m) = qq(nxg - 2*margin + i,m)
     qq(i,m) = qq(2*margin-i+1,m)
  enddo
  enddo

  do i = 1,margin
  do m = 1,mtype
     ! periodic boundary
     !qq(nxg-i+1,m) = qq(2*margin - i + 1 ,m)
     qq(nxg-i+1,m) = qq(nxg - 2*margin + i,m)
  enddo
  enddo
  
  return
end subroutine bc
