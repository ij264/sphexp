program sphexp_inv
  
  use nrtype
  use module_util
  use module_sphexp
  use module_function
  implicit none


  character(len=10) :: string
  character(len=256) :: data_file
  integer(i4b), parameter :: io1 = 7
  integer(i4b) :: narg,ndata,ierr,lmax,ios,idata,ncoef, & 
                  icoef,info,iph,ith,irec,nth,nph,l,jcoef, & 
                  i,lmin,ncoef_max,ncoef_min,m,ilm,il,mp,  &
                  ilm1,ilm2,l2,l1,m1,m2,lmax_exp
                 
  integer(i4b), dimension(:), allocatable :: ipiv
  

  real(dp) :: theta,phi,nu1,nu2,err,chi,lat,lon,dth, & 
              dph,fun,tol,amp,var,fun2
  real(dp), dimension(:), allocatable :: lata,lona,funa,sigmaa,funa_syn
  real(dp), dimension(:), allocatable :: flm
  real(dp), dimension(:,:), allocatable :: resmat
  real(dp), dimension(:,:), allocatable :: covm

  complex(dpc), dimension(:), allocatable :: flm_c

  !---------------------------------------!
  !      read in the input parameters     !
  !---------------------------------------!

  narg = command_argument_count()
  
  if(narg == 0 .or. narg .ne. 6) then
     print *, '======================================================'
     print *, ' code to compute least squares spherical harmonic '
     print *, ' coefficients of a given data set and write out the'
     print *, ' coefficients, power spectra, final grid, final grid '
     print *, ' expanded to a different maximum degree, misfit '
     print *, ' values, resolution and covariance matrices, and '
     print *, ' fits to the input data at the same locations. '
     print *, '======================================================'
     print *, ''
     print *, ' inputs: [data file] [lmin] [lmax] [nu1] [nu2] [lmax_exp]'
     print *, ''
     print *, ' further details:'
     print *, ''
     print *, ' -- data file should have the format: lat lon fun sigma'
     print *, '    angles in degrees, and the variance is optional'
     print *, ''
     print *, ' -- lmin == minimum spherical harmonic degree in expansion'
     print *, ''
     print *, ' -- lmax == maximum spherical harmonic degree in expansion'
     print *, ''
     print *, ' -- nu1 == scaling factor for norm damping'
     print *, ''
     print *, ' -- nu2 == scaling factor for gradient damping'
     print *, ''
     print *, ' -- lmax_exp == maximum degree for smaller expansion'
     print *, ''
     print *, '======================================================'
     stop
  end if


  call get_command_argument(1,data_file)

  call get_command_argument(2,string)
  read(string,*) lmin
  
  call get_command_argument(3,string)
  read(string,*) lmax
  
  call get_command_argument(4,string)
  read(string,*) nu1

  
  call get_command_argument(5,string)
  read(string,*) nu2

  call get_command_argument(6,string)
  read(string,*) lmax_exp


  !---------------------------------------------!
  !             read in the data file           !
  !---------------------------------------------!

  call read_data_SH(io1,data_file,ndata,lata,lona,funa,sigmaa)

  allocate(funa_syn(ndata))
  
  !---------------------------------------------!
  !       solve the least squares problem       !
  !---------------------------------------------!
  ncoef_max = (lmax+1)*(lmax+1)
  ncoef_min = lmin*lmin
  ncoef = ncoef_max-ncoef_min
  allocate(flm(ncoef_max))  
  allocate(resmat(ncoef,ncoef))
  allocate(covm(ncoef,ncoef))


  if(ncoef > ndata) stop ' more parameters than data'

  call least_squares_fit_SH(ndata,lmin,lmax,lata,lona,funa,sigmaa,nu1, & 
       nu2,flm,chi,resmat = resmat,covm = covm, funa_syn = funa_syn)
  print *, ' chi squared over ndata = ',chi/ndata

  ! write out the misfit
  open(io1,file = trim(data_file)//'.chi')
  write(io1,*) chi,ndata,ncoef
  close(io1)


  ! write out the synthetic data
  open(io1,file=trim(data_file)//'.syn')
  do idata = 1,ndata
     write(io1,*) lata(idata),lona(idata),funa_syn(idata),sigmaa(idata)   
  end do
  close(io1)



  ! write out spherical harmonic coefficients and standard deviation  
  open(io1,file=trim(data_file)//'.coef.out')
  icoef = 0
  do l = 0,lmax
     do m = -l,l
        icoef = icoef+1
        write(io1,*) l,m,flm(icoef),sqrt(covm(icoef,icoef))
     end  do
  end do
  close(io1)





  open(io1,file=trim(data_file)//'.coef.power')
  ilm = 0
  do l = 0,lmax


     ! compute the variance of the spectral power
     var = 0.0_dp
     ilm1 = ilm     
     do m = -l,l

        ilm1 = ilm1+1

        icoef = ilm1 - ncoef_min

        ilm2 = ilm
        do mp = -l,l

           ilm2 = ilm2+1
           
           jcoef = ilm2 - ncoef_min

           if(l > lmin) then
              var = var + 4.0_dp*flm(ilm1)*flm(ilm2)*covm(icoef,jcoef)
           end if

        end do

     end do

     
     ! compute the spectral power
     amp = 0.0_dp
     do m = -l,l
        ilm = ilm+1
        icoef = ilm - ncoef_min
        amp = amp + flm(ilm)*flm(ilm)
     end do



     
     ! print out l, spectral ampliture, and its standard deviation
     write(io1,*) l,amp,sqrt(var)
     

  end do



  close(io1)



  !-----------------------------------------!
  !  expand the solution on a global grid   !
  !-----------------------------------------!

  call  sph_polar_grid(lmax)  
  allocate(flm_c(ncoef_sph))
  call real_to_complex_SH(lmax,flm,flm_c)

  open(io1,file=trim(data_file)//'.global_interp')

  
  nth = 181
  nph = 361
  dth = pi_d/(nth-1)
  dph = twopi_d/(nph-1)

  write(io1,*) nth,nph

  ! writes out lat, lon, value from full expansion, value from smaller expansion

  do ith = 1,nth
     theta = (ith-1)*dth     
     do iph = 1,nph
        phi = (iph-1)*dph
        call sphsum(theta,phi,flm_c,fun)
        call sphsum(theta,phi,flm_c,fun2,lmax_exp)
        lat = 90.0_dp-theta*rad2deg
        lon = phi*rad2deg
        write(io1,*) lat,lon,fun,fun2  
     end do
  end do

  close(io1)


  ! writing out the resolution matrix
  open(io1,file=trim(data_file)//'.resmat')

  ilm1 = 0
  icoef = 0
  ! do the first loop of (l,m)
  do l1 = lmin,lmax
     do m1 = -l1,l1
        ilm1 = ilm1+1
        icoef = icoef+1

        ! do second loop over (l',m')
        ilm2 = 0        
        jcoef = 0
        do l2 = lmin,lmax
           do m2 = -l2,l2
              ilm2 = ilm2+1
              jcoef = jcoef+1
              write(io1,*) icoef,jcoef,l1,m1,l2,m2,resmat(ilm1,ilm2)
           end do
        end do
        ! end second loop

     end do
  end do
  ! end first loop

  close(io1)
  

  ! writing out the model covariance matrix
  open(io1,file=trim(data_file)//'.covm')

  ilm1 = 0
  icoef = 0
  ! do the first loop of (l,m)
  do l1 = lmin,lmax
     do m1 = -l1,l1
        ilm1 = ilm1+1
        icoef = icoef+1
              
        ! do second loop over (l',m')
        ilm2 = 0        
        jcoef = 0
        do l2 = lmin,lmax
           do m2 = -l2,l2
              ilm2 = ilm2+1
              jcoef = jcoef+1
              write(io1,*) icoef,jcoef,l1,m1,l2,m2,covm(ilm1,ilm2)
           end do
        end do
        ! end second loop

     end do
  end do
  ! end first loop


  close(io1)
  



end program sphexp_inv
