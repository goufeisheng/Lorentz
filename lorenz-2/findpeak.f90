subroutine findpeak(datax,datay,ndata,peaks,npeak)
  !record non-smooth points, peaks centre&width
  implicit none
  real*8,dimension(ndata),intent(in)     :: datax,datay
  real*8,dimension(npeak,3),intent(out)   :: peaks
  integer,intent(in)  :: ndata
  integer,intent(inout) :: npeak
  real*8  :: thres0
  real*8  :: thres1
  real*8  :: thres2
  real*8  :: thres3
  real*8  :: thresx
  !working variables
  integer :: i,j,k,ipeak,sharpstate,flag
  real*8  :: d1,d2,d3
  real*8  :: slope,slopeold,datamax
  character*5 :: state 
  real*8  :: upx1,downx1,downx2,downy1,maxx1,maxy1
  real*8  :: updx,downdx,maxdx

  state='0'
  npeak=0
  datamax=0  
  do i=1,ndata
    datamax=max(datamax,datay(i))
  enddo
  thres0=datamax/100
  thres1=datamax/100
  thres2=datamax/500
  thres3=datamax/100
  thresx=0.01
!  print*,"datamax:",datamax
!  print*,"threshold:",thres0,thres1,thres2,thres3
!  print*,ndata,datay(1),datax(1)
  ipeak=0
  do i=1+1,ndata-1
    !--------------------------
    ! find max
    !--------------------------
    !slopeold=slope
    slope=datay(i+1)-datay(i-1)
    if(state=='0')then
      if(slope>thres0)then
        state='up'
        upx1=datax(i)
      else if(slope< -thres0)then
        state='down'
        downx1=datax(i)
        downy1=datay(i)
      endif
    else if(state=='up')then
      flag=0
      if(slope<-thres1)then
        state='down'
        updx=datax(i)-upx1
        downy1=datay(i)
        flag=1
        if(.true.)then
          maxdx=0.0
          ipeak=ipeak+1
          peaks(ipeak,2)=datax(i)
          peaks(ipeak,3)=(datax(i)-datax(i-1))*3
          peaks(ipeak,1)=datay(i)*(peaks(ipeak,3)**2)
        endif
      else if(abs(slope)<thres2)then
        state='max'
        updx=datax(i)-upx1
        maxx1=datax(i)
        maxy1=datay(i)
        flag=1
      endif
      if(flag==1)then  
        if(updx>thresx)then
          ipeak=ipeak+1
          peaks(ipeak,2)=upx1+updx/4
          peaks(ipeak,3)=updx/6
          peaks(ipeak,1)=datay(i)/4*(peaks(ipeak,3)**2)
          ipeak=ipeak+1
          peaks(ipeak,2)=upx1+updx/5
          peaks(ipeak,3)=updx/6
          peaks(ipeak,1)=datay(i)/4*(peaks(ipeak,3)**2)
        endif
      endif
    else if(state=='down')then
      if(slope>thres1)then
        state='up'
        downdx=datax(i)-downx1
      else if(abs(slope)<thres2)then
        state='0'
        downdx=datax(i)-downx1
        downx2=datax(i)
        if(downdx>thresx)then
          ipeak=ipeak+1
          peaks(ipeak,2)=downx2-downdx/4
          peaks(ipeak,3)=downdx/6
          peaks(ipeak,1)=downy1/4*(peaks(ipeak,3)**2)
          ipeak=ipeak+1
          peaks(ipeak,2)=downx2-downdx/5
          peaks(ipeak,3)=downdx/6
          peaks(ipeak,1)=downy1/4*(peaks(ipeak,3)**2)
        endif
      endif
    else if(state=='max')then
      flag=0
      if(abs(slope)>thres2)then
        if(slope>0)then
          state='up'
          flag=1
          maxdx=datax(i)-maxx1
        else
          state='down'
          downx1=datax(i)
          downy1=datay(i)
          flag=1
        endif
        if(flag==1)then  
          maxdx=datax(i)-maxx1
          ipeak=ipeak+1
          peaks(ipeak,2)=maxx1+maxdx/2
          peaks(ipeak,3)=maxdx/2
          peaks(ipeak,1)=maxy1*(peaks(ipeak,3)**2)
        endif
      endif
    endif
    !1st order derivatives
    !d1=datay(i+1)-datay(i-1)
    !d2=datay(i+2)-datay(i)
    !d3=datay(i)-datay(i-2)
    !print*,npeak
    !npeak=npeak+1
    !peaks(npeak)=datax(i)
    !width(npeak)=datay(i)
  enddo
  npeak=ipeak
!  open(unit=90,file='opeaks.data')
!  do i=1,ipeak
!    write(90,*),peaks(i,1),peaks(i,2),peaks(i,3)
!  enddo
!  close(90)
end subroutine findpeak
