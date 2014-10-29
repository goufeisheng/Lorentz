program main
  use fit_mod
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! lorentz fitting program
! input: 
!    filename       :  name of file containing fitting data!!!
!    range1, range2 :  (low and up) range of used data
!    outrange1,outrange2: range of output data
!    nlrz
!    nlrzin
!    nlrzout
!    nlrzp
!    inmode
!    fitmode
!    tol           : tolerence in FCN fitting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    function lrzfunc ( x, ilor, coeff)
      real*8,intent(in) :: x
      integer,intent(in) ::  ilor
      real*8,dimension(:),intent(in) :: coeff
      real*8 :: lrzfunc
    end function  
  end interface
  !
  character*100 :: filename
  integer :: maxdata=1000000
  integer :: maxpeak=1000
!  integer,dimension(:),allocatable :: lscreen
  integer :: ndata, nscreen
  real*8,dimension(:),allocatable :: coeff
  integer :: ncoeff,i,j,k,il,ios
  real*8  :: step,x,y,tmp,tmp2
  !-----input parameters
  character*30 :: inmode 
  !character*30 :: fitmode 
  real*8  :: range1,range2,outrange1,outrange2
  integer :: nlrz,nlrzin,nlrzout,nlrzp,ierror
  real*8  :: tol,linewidth
  real*8  :: errmax,errsum,rtmp

  namelist /general/ filename,range1,range2,tol, linewidth

  namelist /lrz/ nlrz,nlrzin,nlrzout,inmode,fitmode

  namelist /output/ outrange1,outrange2

  !default values:
  nlrzin=0
  nlrzout=4
  range1=-2.0
  range2=2.0
  nlrzp=2
  outrange1=-999.0
  outrange2= 999.0
  inmode='uniform'
  !
  rewind 5
  read(5, nml=general, end=110)
  110 continue
  rewind 5
  read(5, nml=lrz, end=120)
  120 continue
  rewind 5
  read(5, nml=output, end=130)
  130 continue

  allocate(lorentz(maxpeak,3))
  allocate(datax(maxdata),datay(maxdata))
  allocate(centre(nlrzin), width(nlrzin))

  call readdata(filename,range1,range2,ndata,nlrzin)

  if(inmode=='auto')then
    print*,"   lrz.x: run in auto mode."
    nlorentz=maxpeak
    call findpeak(datax,datay,ndata,lorentz,nlorentz)
    ncoeff=nlorentz*3
    allocate(coeff(ncoeff))
    !nlrzin=nlorentz
    !nlrzout=0
    print*,"   lrz.x: identified data peaks:", nlorentz,ncoeff
  else
    nlorentz=nlrzin+nlrzout
    ncoeff=nlorentz*3
    allocate(coeff(ncoeff))
    !------------------------------------------
    ! Set Peak Centre & Width
    !------------------------------------------
    !outside range:
    step = (range2-range1)/2
    if(nlrzout>=2)then
      lorentz(1,2) = range1 - step
      lorentz(1,3) = range2-range1
      lorentz(2,2) = range2 + step
      lorentz(2,3) = range2-range1
    endif
    if(nlrzout>=4)then
      lorentz(3,2) = range1 - step/5
      lorentz(3,3) = (range2-range1)/5
      lorentz(4,2) = range2 + step/5
      lorentz(4,3) = (range2-range1)/5
    endif
    !
    !inside range
    step = (range2-range1)
    do i=nlrzout+1,nlorentz
      j=i-nlrzout
      if(inmode=='uniform')then
        lorentz(i,2) = range1+i*step/(nlrzin+1)
        if(mod(i,2)==0)then
          lorentz(i,3) = 5
        else
          lorentz(i,3) = 2
        endif
      else if(inmode=='tb')then
        x=(j*1.0-1-(nlrzin-1)/2.0)*1.0
        print*,x,atan(x)
        lorentz(i,2) = (atan(x))*1.1*2/3.14159265*2.0d0
        lorentz(i,3) = 1.0*exp(-abs(x))*2.0d0
        if(fitmode=='Aw')then
          lorentz(i,1)=sqrt(datay(i))/(lorentz(i,3)**2)
        else if(fitmode=='A')then
          lorentz(i,1)=sqrt(datay(i))
        else if(fitmode=='w')then
          lorentz(i,1)=datay(i)/(lorentz(i,3)**2)
        else     
          lorentz(i,1)=datay(i)
        endif 
      elseif(inmode=='fixed')then       
        lorentz(i,2) = centre(j)
        lorentz(i,3) = width(j)

        lorentz(i,1) = datay(i)
      else
        lorentz(i,2) = 0.0d0
        lorentz(i,3) = 1.0d0
      endif
    enddo
  endif
  !-----------------------------------------------
  print*,"   lrz.x: initial guess:"
  do i=1,nlorentz
      print*,"        " ,i,lorentz(i,2),lorentz(i,3)
  enddo
  !

  !

  !----------------------------------
  !set initial guess for test07
  !----------------------------------
  do i=1,nlorentz
      coeff(i)=lorentz(i,1)
      coeff(nlorentz+i)=lorentz(i,2)
      coeff(nlorentz*2+i)=lorentz(i,3)
  enddo

!  call timestamp ()
  print*,"   lrz.x: start fitting, please wait....."
  call test07 (datax,datay,ndata,coeff,ncoeff,tol,errmax)
!  call timestamp ()
  !---------
  !output
  !---------
!  print*,"c="
!  do i=1,ncoeff
!    print*,i,coeff(i)
!  enddo
  print*,"   lrz.x:  errmax=",errmax
  do i=1,nlorentz
      lorentz(i,1)=coeff(i)
      lorentz(i,2)=coeff(nlorentz+i)
      lorentz(i,3)=coeff(nlorentz*2+i)
  enddo

  do i=1,nlorentz-1
     do j=i+1, nlorentz
        if ( dabs(lorentz(i,2)-lorentz(j,2)) .le. 1.d-3 .and. dabs(lorentz(i,3)-lorentz(j,3)).le. 1.d-3 )then
           if( fitmode=='Aw' .or. fitmode=='A') then
               lorentz(i,1) = dsqrt(lorentz(i,1)**2 + lorentz(j,1)**2)
           else
               lorentz(i,1) = lorentz(i,1) + lorentz(j,1)
           endif
        lorentz(j,1) = 0.00
        endif
     enddo
  enddo
          

  allocate(lscreen(nlorentz))
  do i=1,nlorentz  
    lscreen(i)=0
  enddo
  nscreen=0
  do i=1,nlorentz
   if( fitmode=='Aw' .or. fitmode=='A') then
     if ( lorentz(i,1)**2/linewidth .le. 1.d-3 .or. dabs(lorentz(i,2)) .ge. range2) then
       lscreen(i)=1
       nscreen=nscreen+lscreen(i)
     endif
   else
     if ( dabs(lorentz(i,1)/linewidth) .le. 1.d-3 .or. dabs(lorentz(i,2)) .ge. range2) then
       lscreen(i)=1
       nscreen=nscreen+lscreen(i)   
     endif
   endif
enddo
  open(unit=91,file="eachlrz.data")
  do j=1,nlorentz
    if ( lscreen(j)==0) then
      do x=-10.0, 10.0, 0.01
        write(91,*),x,lrzfunc(x,j,coeff)
      enddo
        write(91,*)
        write(91,*)
    endif
  enddo
  close(91)

  call outputfit(filename,lorentz,nlorentz,maxpeak,outrange1,outrange2,fitmode, &
          lscreen)
  do i=1,nlorentz
     lorentz(i,3)=dabs(lorentz(i,3))
  enddo
  open(unit=92,file='peaks.data')
  print*,"   lrz.x: used peak number:", nlorentz-nscreen
  print*,"          A(i)        width(i)         center(i)"
  do i=1,nlorentz
    if ( lscreen(i)==0)then
      if(fitmode=='Aw' .or. fitmode=='A')then
        write(92,518),lorentz(i,1)**2/linewidth,lorentz(i,3), lorentz(i,2)
        print*,"      ",lorentz(i,1)**2/linewidth,lorentz(i,3), lorentz(i,2)
      else 
        write(92,518),lorentz(i,1), lorentz(i,3), lorentz(i,2)
        print*,"      ",lorentz(i,1), lorentz(i,3), lorentz(i,2)
      endif
    endif
  enddo
  close(92)
  print*,"   lrz.x: unused peak number:", nscreen
  print*,"          A(i)        width(i)         center(i)"
  do i=1,nlorentz
    if ( lscreen(i)==1)then
      if(fitmode=='Aw' .or. fitmode=='A')then
        print*,"      ",lorentz(i,1)**2/linewidth,lorentz(i,3), lorentz(i,2)
      else 
        print*,"      ",lorentz(i,1), lorentz(i,3), lorentz(i,2)
      endif
    endif
  enddo
  !
  write ( *, '(a)' ) '   lrz.x: end of lorentz fitting. Check INFO!'
518 format( 3(e12.5e2, 2x))
end program

subroutine test07 (datax,datay,ndata,coeff,ncoeff,tol,errmax)

  !*****************************************************************************80
  !
  !! TEST07 tests LMDIF1.
  !
  !  Discussion:
  !
  !    LMDIF1 solves M nonlinear equations in N unknowns, where M is greater
  !    than N.  It is similar to test02, except that the functional fit is
  !    nonlinear this time, of the form
  !
  !      y = a + b * x**c,
  !
  !    with x and y data, and a, b and c unknown.
  !
  !    This problem is set up so that the data is exactly fit by by
  !    a=1, b=3, c=2.  Normally, the data would only be approximately
  !    fit by the best possible solution.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 December 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  implicit none
  real*8,dimension(ndata),intent(in)  :: datax,datay
  integer,intent(in) :: ndata,ncoeff
  real*8,intent(in) :: tol
  real*8,dimension(ncoeff),intent(out) :: coeff
  real*8,intent(out) :: errmax
  external   fitfunc 

  !working var
  integer :: i,info,iflag
  real*8,dimension(:),allocatable :: fvec
  real*8 :: errsum
!  print*,"N:",ndata,ncoeff,tol
  allocate(fvec(ndata))
!  print*,"------------------------------"
  !print*,"initial err:"
  call fitfunc (ndata, ncoeff, coeff, fvec, iflag )

  !x(1:3) = (/ 0.0D+00, 5.0D+00, 1.3D+00 /)
  !call r8vec_print ( n, x, '  X:' )
  iflag = 1

  !call r8vec_print ( m, fvec, '  F(X):' )

  call lmdif1 ( fitfunc, ndata, ncoeff, coeff, fvec, tol, info )

  print*,""
  print*,"   lrz.x:  Returned INFO = ", info
  !call r8vec_print ( n, x, '  X:' )
  !call r8vec_print ( m, fvec, '  F(X):' )

  errsum = 0.0
  errmax = 0.0
!  open(unit=90,file="result.data")
  do i=1,ndata
    errsum=errsum+abs(fvec(i))
    errmax=max(errmax,abs(fvec(i)))
    !write(90,*),datax(i),fvec(i)+datay(i),abs(fvec(i))
  enddo
!  close(90)
end subroutine test07

