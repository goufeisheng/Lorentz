subroutine fitfunc ( ndata, ncoeff, coeff, fvec, iflag )
  use fit_mod
  implicit none
 interface
   function lrzfunc ( x, ilor, coeff)
     real*8,intent(in) :: x
     integer,intent(in) ::  ilor
     real*8,dimension(:),intent(in) :: coeff
     real*8 :: lrzfunc
   end function  
 end interface


  integer,intent(in) ::  ncoeff,ndata
  real*8,intent(in) :: coeff(ncoeff)
  real*8,dimension(ndata),intent(out) ::  fvec
  integer ( kind = 4 ) iflag

  !working var
  integer :: i,j
  real*8 ::  rtmp
  fvec(1:ndata)=0d0
  do i=1,nlorentz
    !fvec(1:ndata) =fvec(1:ndata)+coeff((i-1)*2+1)/((datax(1:ndata)-(i*3.0-3.0))**2+coeff((i-1)*2+2)**2)
    !fvec(1:ndata) =fvec(1:ndata)+coeff((i-1)*nlrzparam+1)/((datax(1:ndata)-(i*3.0-3.0))**2+9.0)
    !fvec(1:ndata) =fvec(1:ndata)+coeff((i-1)*3+1)/((datax(1:ndata)-coeff((i-1)*3+2))**2+coeff((i-1)*3+3)**2)
    do j=1,ndata
      fvec(j)=fvec(j)+lrzfunc(datax(j),i,coeff) 
      !call lrzfunc(datax(j),i,coeff,rtmp,0)
      !fvec(i)=fvec(i)+rtmp 
    enddo
  enddo 
  !coeff(1) + coeff(2) * datax(1:ndata)**coeff(3) +coeff(4)*datax(1:ndata)
 
  do i =1, ndata 
      if(fvec(i) .lt. 0.0d0) then 
          fvec(i) =  fvec(i) - 0.01
      endif
  enddo

  fvec(1:ndata) = fvec(1:ndata)- datay(1:ndata)

end subroutine fitfunc
