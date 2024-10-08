program sphexp_data_idlm

  use nrtype
  use module_util
  use module_ice5g
  use module_sphexp
  use ran_state
  implicit none


  character(len=10) :: string
  character(len=256) :: data_in
  character(len=256) :: data_out
  integer(i4b), parameter :: io1 = 7
  integer(i4b) :: lmax,ilat,ilon,nsample,isample,msec, & 
       narg,l,m,ll,mm,ilm,ndata,ierr,idata,lpeak,nth,nph, & 
       ith,iph,nd,lmin,ilm1,ilm2,lc,mc
  integer(i4b), dimension(8) :: time_info

  real(dp) :: rdate,rfilt,lat,lon,sigma,ran,fun,ran2, & 
              sigma_ran,sigma0,theta,phi,am,dth,dph,lambda,amp,fac
  real(dp), dimension(:), allocatable :: flm
  real(dp), dimension(:), allocatable :: lata,lona,funa,sigmaa

  real(dp), dimension(:,:), allocatable :: a
  
  complex(dpc), dimension(:), allocatable :: flm_c

  ! read in the input parameters
  narg = command_argument_count()
  
  if(narg == 0 .or. narg .ne. 5) then
     print *, ' inputs: [data in] [l] [m] [ierr] [amp]'

     print *, '======================================================'
     print *, ' code to compute least squares spherical harmonic '
     print *, ' coefficients of a single input spherical harmonic '
     print *, ' sampled at locations contained within a data file. '
     print *, ' Errors can be read from the data file and additional '
     print *, ' uncertainty scaled by these errors can be added to '
     print *, ' the spherical harmonic chosen to simulate noise '
     print *, '======================================================'
     print *, ''
     print *, ' inputs: [data file] [l] [m] [ierr] [amp]'
     print *, ''
     print *, ' further details:'
     print *, ''
     print *, ' -- data file should have the format:'
     print *, '    lat lon dummy_value uncertainty_at_that_datapoint dummy_number'
     print *, '    dummy numbers have no effect (relict)'
     print *, '    uncertainty is used to weight the misfit function'
     print *, '    set uncertainty to same value to equal weight all points'
     print *, ''
     print *, ' -- l == degree of the input spherical harmonic'
     print *, ''
     print *, ' -- m == order of the input spherical harmonic'
     print *, ''
     print *, ' -- ierr == 0 or 1 [add random noise to input s.h.]'
     print *, ''
     print *, ' -- amp == amplitude of the input coefficient'
     print *, ''
     print *, '======================================================'
     stop


     stop
  end if
 
  call get_command_argument(1,data_in)

  call get_command_argument(2,string)
  read(string,*) lc

  call get_command_argument(3,string)
  read(string,*) mc

  call get_command_argument(4,string)
  read(string,*) ierr

  call get_command_argument(5,string)
  read(string,*) amp



  

  data_out = trim(data_in)//'.id'

  ! set up the spherical harmonic expansions
  lmax = lc
  call set_SHTOOLS(lc)

  allocate(flm(ncoef_real_SH))

  open(io1,file=trim(data_out)//'.coef')

  ! set all but the coefficient being tested to 0
  ilm = 0
  do l = 0,lmax
     do m = -l,l
        ilm = ilm+1
        if(l == lc .and. m == mc) then
           flm(ilm) = amp
        else
           flm(ilm) = 0.0_dp
        end if
        write(io1,*) l,m,flm(ilm),0
     end do
  end do

  close(io1)




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

end program sphexp_data_idlm
