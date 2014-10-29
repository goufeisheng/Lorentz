
subroutine outputfit(filename,lorentz,nlorentz,maxpeak,range1,range2,fitmode, & 
                     lscreen)
  character*30 :: fitmode
  character*100,intent(in) :: filename
  real*8,intent(in) :: range1,range2
  real*8,dimension(maxpeak,3),intent(in) :: lorentz
  integer,intent(in) :: nlorentz
  integer,dimension(nlorentz),intent(in) :: lscreen
  complex*16 :: ctmp3 

  integer :: i,ios
  real*8  :: x,y,yfit,ysum
  real*8  :: w2
  !print*,"-------------------------------------"
  !print*,"----------------outputfit-----------------"
  !print*,nlorentz
  !print*,lorentz(1,1:3)
  open(unit=19,file='hybrd_heom_in.data')
  open(unit=90,file="hybrd_heom_fitted.data")
  open(unit=99,file='se_heom_fitted.data')  
  rewind(19)
  rewind(90)
  rewind(99)
  !count number of data
 
  do
    read(19,*,iostat=ios),x,y
    if(ios/=0)then
      exit
    endif
      yfit=0.0
      ctmp3=(0.d0, 0.d0)
      do i=1,nlorentz
        if ( lscreen(i)==0) then
          w2=(lorentz(i,3)**2)
          if(fitmode=='Aw')then
            yfit=yfit+(lorentz(i,1)**2)*w2/(((x-lorentz(i,2))**2)+w2)
          ctmp3 = ctmp3 + (lorentz(i,1)**2)*dabs(lorentz(i,3))/dcmplx(x-lorentz(i,2), dabs(lorentz(i,3)))
          else if(fitmode=='w')then
            yfit=yfit+(lorentz(i,1))*w2/(((x-lorentz(i,2))**2)+w2)
!#          ctmp3 = ctmp3 + lorentz(i,1)*dabs(lorentz(i,3))/dcmplx(x-lorentz(i,2), dabs(lorentz(i,3)))
          ctmp3 = ctmp3 + lorentz(i,1)*lorentz(i,3)/dcmplx((x-lorentz(i,2)), lorentz(i,3))
          else 
            yfit=yfit+(lorentz(i,1))/(((x-lorentz(i,2))**2)+w2)
           ctmp3 = ctmp3 + lorentz(i,1)*dabs(lorentz(i,3))/dcmplx(x-lorentz(i,2), dabs(lorentz(i,3)))
          endif
        endif
      enddo
      write(90,518),x,yfit !,yfit-y
      write(99,518),x, dble(ctmp3), dimag(ctmp3)
  enddo
  close(90)
  close(99)
  close(19)
  open(unit=100, file='fitted_long.data')
  rewind(100)
do x=-10.0,10.0,0.01
      yfit=0.0
      do i=1,nlorentz
        if ( lscreen(i)==0) then
          w2=(lorentz(i,3)**2)
          if(fitmode=='Aw')then
            yfit=yfit+(lorentz(i,1)**2)*w2/(((x-lorentz(i,2))**2)+w2)
          else if(fitmode=='w')then
            yfit=yfit+(lorentz(i,1))*w2/(((x-lorentz(i,2))**2)+w2)
          else
            yfit=yfit+(lorentz(i,1))/(((x-lorentz(i,2))**2)+w2)
          endif
        endif
      enddo
      write(100,518),x,yfit !,yfit-y
  enddo
518 format(5(2x, e18.6e3))
end subroutine outputfit


