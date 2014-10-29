subroutine readdata(filename,range1,range2,ndata,nlrzin)
  use fit_mod
  implicit none
  integer :: nlrzin
  character*100,intent(in) :: filename
  real*8,intent(in) :: range1,range2
  !real*8,dimension(:),intent(out) :: datax,datay
  integer,intent(out) :: ndata

  integer :: i,j,ios
  real*8  :: tmp,tmp1,tmp2


  open(unit=19,file=filename)
  rewind 19
  ndata=0
  do
    read(19,*,iostat=ios),tmp,tmp2
    if(ios/=0)then
      exit
    endif
    if(tmp>=range1 .and. tmp <= range2) ndata=ndata+1
  enddo
  print*,"   lrz.x: data points: ", ndata
 !read data
  rewind(19)
  i=1
  do
    read(19,*,iostat=ios),tmp,tmp2
    if(tmp>=range1 .and. tmp <= range2)then
      datax(i)=tmp
      datay(i)=tmp2
      i=i+1
    endif
    if(ios/=0)then
      exit
    endif
  enddo
  close(19)

  tmp = 0.0d0 
  do i = 1, ndata
       tmp = tmp + dabs(datay(i))
       !Y=DABS(X)    倍精度实数绝对值    REAL*8    REAL*8  
  enddo
  tmp = tmp/dble(ndata*10)
  !Y=DBLE(X)    转换为倍精度实数    ALL    REAL*8  
  tmp1 = 1.0d-50
  tmp2 = 1.0d50
  do i = 1, ndata 
       if(datax(i) .gt. tmp1 .and. dabs(datay(i)) .gt. tmp) then
       !GT是大于号(>)
            tmp1 = datax(i)
       endif 
       if(datax(i) .lt. tmp2 .and. dabs(datay(i)) .gt. tmp) then
       !LT是小于号(<)
            tmp2 = datax(i)
       endif
  enddo

  tmp = (tmp1 - tmp2)/dble(nlrzin/2-1) 

  do i = 1, nlrzin/2
      do j = 1, 2
           centre(i+nlrzin/2*(j-1)) = tmp2 + (dble(i-1)*tmp) 
           width(i+nlrzin/2*(j-1)) = dble(j) * (tmp/2.0d0)
      enddo
  enddo


end subroutine readdata
