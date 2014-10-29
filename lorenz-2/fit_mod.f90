module fit_mod
  real*8,dimension(:),allocatable,save :: datax,datay
  !real*8,dimension(:,:),allocatable,save :: coeff
  real*8,dimension(:,:),allocatable,save :: lorentz
  real*8,dimension(:),allocatable :: centre, width
  integer,dimension(:),allocatable,save :: lscreen
 integer,save :: nlorentz
  character*30 :: fitmode
end module fit_mod
