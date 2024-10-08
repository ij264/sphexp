program sphexp_chi_comp

  use nrtype
  use module_util
  use module_sphexp
  use module_function
  implicit none


  character(len=10) :: string
  character(len=256) :: data_file
  character(len=256) :: coef_file
  integer(i4b), parameter :: io1 = 7

  integer(i4b) :: narg,ierr,ndata,lmax,ncoef,idata

  real(dp) :: chi,err
  real(dp), dimension(:), allocatable :: lata,lona,funa,sigmaa

  real(dp), dimension(:), allocatable :: flm,flm_var
  real(dp), dimension(:,:), allocatable :: a

  !---------------------------------------!
  !      read in the input parameters     !
  !---------------------------------------!

  narg = command_argument_count()
  
  if(narg == 0 .or. narg .ne. 2) then
     print *, '======================================================'
     print *, ' code to compute least squares spherical harmonic '
     print *, ' coefficients of a given data set'
     print *, '======================================================'
     print *, ''
     print *, ' inputs: [data file] [coef file]'
     print *, ''
     print *, ' further details:'
     print *, ''
     print *, ' -- data file should have the format: lat lon fun sigma'
     print *, '    angles in degrees, and the variance is optional'
     print *, ''
     print *, ' -- coef file contains a list of spherical harmonic    '
     print *, '    coefficients in the format produced by sphexp_inv  '
     print *, '======================================================'
     stop
  end if

  call get_command_argument(1,data_file)
  
  call get_command_argument(2,coef_file)
  

  !---------------------------------------------!
  !             read in the data file           !
  !---------------------------------------------!
  call read_data_SH(io1,data_file,ndata,lata,lona,funa,sigmaa)
  
  
  !---------------------------------------------!
  !         read in the coefficient file        !
  !---------------------------------------------!
  call read_coef_SH(io1,coef_file,lmax,flm,flm_var)
  ncoef = (lmax+1)*(lmax+1)

  
  !---------------------------------------------!
  !     build the spherical harmonic matrix     !
  !---------------------------------------------!
  allocate(a(ndata,ncoef))
  call build_ylm_matrix_SH(0,lmax,ndata,lata,lona,a)


  !---------------------------------------------!
  !              compute chi squared            !
  !---------------------------------------------!
  chi = 0.0_dp
  do idata = 1,ndata     
     err = dot_product(a(idata,:),flm(:))
     err = err-funa(idata)
     chi = chi + err**2/sigmaa(idata)**2
  end do

  print *, ' chi squared over ndata = ',chi/ndata

  open(io1,file = trim(data_file)//'.chi')
  write(io1,*) chi,ndata,ncoef
  close(io1)

  
end program sphexp_chi_comp
