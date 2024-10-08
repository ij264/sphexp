program sphexp_backus_gilbert
  
  use nrtype
  use module_util
  use module_sphexp
  use module_function
  implicit none


  character(len=10) :: string
  character(len=256) :: data_file
  integer(i4b), parameter :: io1 = 7
  integer(i4b), parameter :: etol = 6
  integer(i4b) :: narg,ndata,ierr,lmax,ios,idata,ncoef, & 
                  icoef,info,iph,ith,irec,nth,nph,l,jcoef,i, & 
                  lout,ncoef_out,lmin,ilm,m
  integer(i4b), dimension(:), allocatable :: ipiv
  

  real(dp) :: theta,phi,nu,err,lat,lon,dth,dph,fun,tol,rvar,amp
  real(dp), dimension(:), allocatable :: lata,lona,funa,sigmaa,lambdaa
  real(dp), dimension(:), allocatable :: flm
  real(dp), dimension(:), allocatable :: chi
  real(dp), dimension(:), allocatable :: var


  complex(dpc), dimension(:), allocatable :: flm_c

  !---------------------------------------!
  !      read in the input parameters     !
  !---------------------------------------!

  narg = command_argument_count()
  
  if(narg == 0 .or. narg .ne. 4) then
     print *, '======================================================'
     print *, ' code to compute least squares spherical harmonic '
     print *, ' coefficients of a given data set'
     print *, '======================================================'
     print *, ''
     print *, ' inputs: [data file] [lmin] [nu] [lout]'
     print *, ''
     print *, ' further details:'
     print *, ''
     print *, ' -- data file should have the format: lat lon fun sigma'
     print *, '    angles in degrees, and the variance is optional'
     print *, ''
     print *, ' -- nu == scaling factor for the variance damping'
     print *, ''
     print *, ' -- lout == range of coefficients to estimate'
     print *, '======================================================'
     stop
  end if


  call get_command_argument(1,data_file)

  call get_command_argument(2,string)
  read(string,*) lmin
  
  call get_command_argument(3,string)
  read(string,*) nu
  
  call get_command_argument(4,string)
  read(string,*) lout



  if(lout < lmin) stop ' lout < lmin'


  !---------------------------------------------!
  !             read in the data file           !
  !---------------------------------------------!

  call read_data_SH_av(io1,data_file,ndata,lata,lona,funa,sigmaa,lambdaa)

  lmax = sqrt(real(etol))/minval(lambdaa)

  print *, 'lmax for data = ',lmax

  if(lout > lmax) lout = lmax

  !---------------------------------------------!
  !       solve the Backus-Gilbert Problem      !
  !---------------------------------------------!

  ncoef_out = (lout+1)*(lout+1)
  allocate(flm(ncoef_out))
  allocate(chi(ncoef_out))
  allocate(var(ncoef_out))
  
  call  backus_gilbert_estaimte_SH(ndata,lmin,lmax,lata,lona,funa, & 
                                   sigmaa,lambdaa,nu,lout,flm,chi,var)
  



  open(io1,file=trim(data_file)//'.coef.out')
  ilm = 0
  do l = 0,lout
     do m = -l,l
        ilm = ilm+1
        write(io1,*) l,m,flm(ilm),chi(ilm),var(ilm)
     end do
  end do  
  close(io1)



  open(io1,file=trim(data_file)//'.coef.power')
  ilm = 0
  do l = 0,lout
     
     amp = 0.0_dp
     rvar =  0.0_dp
     do m = -l,l
        
        ilm = ilm+1
        amp = amp + flm(ilm)*flm(ilm)
        rvar = rvar + 4.0_dp*flm(ilm)*flm(ilm)*var(ilm)

     end do

     write(io1,*) l,amp,sqrt(rvar)
     
  end do




  close(io1)


 

  !-----------------------------------------!
  !  expand the solution on a global grid   !
  !-----------------------------------------!

  call  sph_polar_grid(lout)  
  allocate(flm_c(ncoef_sph))
  call real_to_complex_SH(lout,flm,flm_c)

  open(io1,file=trim(data_file)//'.global_interp')

  
  nth = nth_sph*2
  nph = nph_sph*2
  dth = pi_d/(nth-1)
  dph = twopi_d/(nph-1)

  write(io1,*) nth,nph

  do ith = 1,nth
     theta = (ith-1)*dth     
     do iph = 1,nph
        phi = (iph-1)*dph
        call sphsum(theta,phi,flm_c,fun)
        lat = 90.0_dp-theta*rad2deg
        lon = phi*rad2deg
        write(io1,*) lat,lon,fun        
     end do
  end do

  close(io1)


  

end program sphexp_backus_gilbert
