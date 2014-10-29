function lrzfunc ( x, ilor, coeff)
  use fit_mod
  implicit none
  real*8,intent(in) :: x
  integer,intent(in) ::  ilor
  real*8,dimension(:),intent(in) :: coeff
  real*8 :: value
  real*8 :: lrzfunc
  !working var
  integer :: i,j
  real*8 :: w2
  value=0.0
    w2=coeff(nlorentz*2+ilor)**2
    if(fitmode=='Aw')then
      value=(coeff(ilor)**2)*w2/((x-coeff(nlorentz+ilor))**2+w2)
    else if(fitmode=='w')then
      value=(coeff(ilor))*w2/((x-coeff(nlorentz+ilor))**2+w2)
    else if(fitmode=='A')then
      value=(coeff(ilor)**2)/((x-coeff(nlorentz+ilor))**2+w2)
    else
      value=(coeff(ilor))/((x-coeff(nlorentz+ilor))**2+w2)
    endif
  lrzfunc=value
end function lrzfunc
