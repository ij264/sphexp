program sphexp_data_id

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
       ith,iph,nd,lmin,ilm1,ilm2
  integer(i4b), dimension(8) :: time_info

  real(dp) :: rdate,rfilt,lat,lon,sigma,ran,fun,ran2, & 
              sigma_ran,sigma0,theta,phi,am,dth,dph,lambda,amp,fac
  real(dp), dimension(:), allocatable :: flm
  real(dp), dimension(:), allocatable :: lata,lona,funa,sigmaa

  real(dp), dimension(:,:), allocatable :: a
  
  complex(dpc), dimension(:), allocatable :: flm_c

  ! read in the input parameters
  narg = command_argument_count()
  
  if(narg == 0 .or. narg .ne. 7)then
     print *, ' inputs: [data in] [lmin] [lmax] [lpeak] [sigma] [lambda] [nd]'
     stop
  end if
 
  call get_command_argument(1,data_in)

  call get_command_argument(2,string)
  read(string,*) lmin

  call get_command_argument(3,string)
  read(string,*) lmax

  call get_command_argument(4,string)
  read(string,*) lpeak

  call get_command_argument(5,string)
  read(string,*) sigma

  call get_command_argument(6,string)
  read(string,*) lambda

  call get_command_argument(7,string)
  read(string,*) nd

  data_out = trim(data_in)//'.id'

  ! set up the spherical harmonic expansions
  call set_SHTOOLS(lmax)

  allocate(flm(ncoef_real_SH))

  open(io1,file=trim(data_out)//'.coef')

  am = lpeak/sqrt2

  ilm = 0
  do l = 0,lmax
     do m = -l,l
        ilm = ilm+1
        if(l >= lmin) then
           call random_number(ran)        
           ran = ran-0.5_dp
           flm(ilm) = ran*l**2*exp(-0.5_dp*l*l/(am*am))/am**3  & 
                + 0.1*ran*lpeak**2*exp(-0.5_dp*lpeak*lpeak/(am*am))/am**3
        else
           flm(ilm) = 0.0_dp
        end if
        write(io1,*) flm(ilm)
     end do
  end do

  close(io1)


  open(io1,file=trim(data_out)//'.coef.power.in')

  am = lpeak/sqrt2

  ilm = 0
  do l = 0,lmax
     amp = 0.0_dp
     do m = -l,l
        ilm = ilm+1
        call random_number(ran)        
        ran = ran-0.5_dp
        flm(ilm) = ran*l**2*exp(-0.5_dp*l*l/(am*am))/am**3  & 
             + 0.1*ran*lpeak**2*exp(-0.5_dp*lpeak*lpeak/(am*am))/am**3
        amp = amp + flm(ilm)*flm(ilm)
     end do
     write(io1,*) l,amp,0.0_dp
  end do

  close(io1)



  call  read_data_SH(io1,data_in,ndata,lata,lona,funa,sigmaa)

  open(io1,file=trim(data_out))


  ilm2 = 0
  do l = 0,lmax
     ilm1 = ilm2+1
     ilm2 = ilm1+l
     fac = lambda/6371.0_dp
     fac = fac*fac
     fac = fac*l*(l+1)
     fac = exp(-fac)
     flm(ilm1:ilm2) = flm(ilm1:ilm2)*fac
  end do

  
  do idata = 1,ndata,nd

     theta = (90.0_dp-lata(idata))*deg2rad
     phi   = lona(idata)*deg2rad
     call sphsum_real(lmax,theta,phi,flm,fun)

     call random_number(ran)
     ran = ran-0.5_dp
     
     write(io1,*) lata(idata),lona(idata),fun+ran*sigma,sigma,lambda
     

  end do

  close(io1)



  !-----------------------------------------!
  !  expand the solution on a global grid   !
  !-----------------------------------------!

  call  sph_polar_grid(lmax)  
  allocate(flm_c(ncoef_sph))
  call real_to_complex_SH(lmax,flm,flm_c)

  open(io1,file=trim(data_out)//'.global_interp.syn')

  
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





end program sphexp_data_id
