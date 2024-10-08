program sphexp_global_grid

  use nrtype
  use module_util
  use module_ice5g
  use module_sphexp
  use ran_state
  implicit none


  character(len=10) :: string
  character(len=256) :: coef_file
  integer(i4b), parameter :: io1 = 7
  integer(i4b) :: lmax,ilat,ilon,nsample,isample,msec, & 
       narg,l,m,ll,mm,ilm,ndata,ierr,idata,lpeak,nth,nph, & 
       ith,iph,nd,lmin,ilm1,ilm2,lc,mc
  
  real(dp) :: rdate,rfilt,lat,lon,sigma,ran,fun,ran2, & 
              sigma_ran,sigma0,theta,phi,am,dth,dph,lambda,amp,fac
  real(dp), dimension(:), allocatable :: flm,flm_var
  

  
  

  ! read in the input parameters
  narg = command_argument_count()
  
  if(narg == 0 .or. narg .ne. 1)then
     print *, ' inputs: [coef in]'

     print *, '======================================================'
     print *, ' code to read in spherical harmonic coefficients and'
     print *, ' generate a global grid sampled every 1 degree'
     print *, '======================================================'
     print *, ''
     print *, ' inputs: [coef file]'
     print *, ''
     print *, ' input format of Davids other codes'
     print *, ''
     print *, ' outputs: [data file]'
     print *, ''
     print *, ' output format: lat lon value'
     print *, ''
     print *, '======================================================'
     stop
  end if
 
  call get_command_argument(1,coef_file)

  ! open and read the coefficient file
  call read_coef_SH(io1,coef_file,lmax,flm,flm_var)



  ! set up the spherical harmonic expansions (generate some of the factors called later)
  call set_SHTOOLS(lmax)




  ! generate a global grid of the input spherical harmonic

  open(io1,file=trim(coef_file)//'.global_interp')

  nth = 181
  nph = 360
  dth = pi_d/(nth-1)
  dph = twopi_d/(nph)

  write(io1,*) nth,nph

  do ith = 1,nth
     theta = (ith-1)*dth     
     do iph = 1,nph
        phi = (iph-1)*dph
        call sphsum_real(lmax,theta,phi,flm,fun)
        lat = 90.0_dp-theta*rad2deg
        lon = phi*rad2deg
        write(io1,*) lat,lon,fun        
     end do
  end do

  close(io1)



end program sphexp_global_grid
