program sphexp_data_power_spectrum_recovery

  use nrtype
  use module_util
  use module_ice5g
  use module_sphexp
  use ran_state
  implicit none


  character(len=10) :: string
  character(len=256) :: data_in
  character(len=256) :: data_out
  character(len=256) :: coef_in
  integer(i4b), parameter :: io1 = 7
  integer(i4b) :: lmax,ilat,ilon,nsample,isample,msec, & 
       narg,l,m,ll,mm,ilm,ndata,ierr,idata,lpeak,nth,nph, & 
       ith,iph,nd,lmin,ilm1,ilm2,lc,mc
  integer(i4b), dimension(8) :: time_info

  real(dp) :: rdate,rfilt,lat,lon,sigma,ran,fun,ran2, & 
              sigma_ran,sigma0,theta,phi,am,dth,dph,lambda,amp,fac
  real(dp), dimension(:), allocatable :: flm,flm_var
  real(dp), dimension(:), allocatable :: lata,lona,funa,sigmaa

  real(dp), dimension(:,:), allocatable :: a
  
  complex(dpc), dimension(:), allocatable :: flm_c

  ! read in the input parameters
  narg = command_argument_count()
  
  if(narg == 0 .or. narg .ne. 3) then

     print *, '======================================================'
     print *, ' code to sample a spherical harmonic grid at a range'
     print *, ' of input data locations and optionally add random '
     print *, ' errors to the results normalised by the uncertainty'
     print *, ' of the input data points '
     print *, '======================================================'
     print *, ''
     print *, ' inputs: [data file] [coef in] [ierr]'
     print *, ''
     print *, ' further details:'
     print *, ''
     print *, ' -- data file should have the format:'
     print *, '    lat lon dummy_value uncertainty_at_that_datapoint dummy_number'
     print *, '    dummy numbers have no effect (relict)'
     print *, '    uncertainty is used to weight the misfit function'
     print *, ''
     print *, ' -- coef in == input s.h. coefficients'
     print *, ''
     print *, ' -- ierr == 0 or 1 [add random noise to input s.h.]'
     print *, ''
     print *, '======================================================'
     stop


     stop
  end if
 
  call get_command_argument(1,data_in)

  call get_command_argument(2,coef_in)

  call get_command_argument(3,string)
  read(string,*) ierr


  ! open and read the coefficient file
  call read_coef_SH(io1,coef_in,lmax,flm,flm_var)

  

  data_out = trim(data_in)//'.id'

  ! set up the spherical harmonic expansions
  call set_SHTOOLS(lmax)

  call  read_data_SH(io1,data_in,ndata,lata,lona,funa,sigmaa)



  open(io1,file=trim(data_out))


  ! write out the value of the synthetic spherical harmonic at each location within the input data file 
  do idata = 1,ndata

     theta = (90.0_dp-lata(idata))*deg2rad
     phi   = lona(idata)*deg2rad
     call sphsum_real(lmax,theta,phi,flm,fun)

     call random_number(ran)
     ran = ran-0.5_dp

     ! add additional noise to the synthetic scaled by the data uncertainties
     if(ierr == 1) then
        sigma = sigmaa(idata)
     else if(ierr == 0) then
        sigma = 0.0_dp
     else
        stop 'bad value of ierr'
     end if

     write(io1,*) lata(idata),lona(idata),fun+ran*sigma,sigmaa(idata)
     

  end do

  close(io1)

end program sphexp_data_power_spectrum_recovery
