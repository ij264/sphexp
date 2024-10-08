program sphexp_data

  use nrtype
  use module_util
  use module_ice5g
  use module_sphexp
  use ran_state
  implicit none


  character(len=10) :: string
  integer(i4b), parameter :: io1 = 7
  integer(i4b) :: lmax,ilat,ilon,nsample,isample,msec,narg
  integer(i4b), dimension(8) :: time_info

  real(dp) :: rdate,rfilt,lat,lon,sigma,ran,fun,ran2, & 
              sigma_ran,sigma0
  real(dp), dimension(:,:), allocatable :: sl_grid

  complex(dpc), dimension(:), allocatable :: sl_lm


  ! read in the input parameters
  narg = command_argument_count()
  
  if(narg == 0) then
     print *, ' inputs: [lmax] [nsample] [sigma]'
     stop
  end if

  if((narg .ne. 3)) then
     print *, 'wrong number of arguments'
      print *, ' inputs: [lmax] [nsample] [sigma]'
     stop     
  end if
 
  call get_command_argument(1,string)
  read(string,*) lmax
  

  call get_command_argument(2,string)
  read(string,*) nsample


  call get_command_argument(3,string)
  read(string,*) sigma0


  ! set up the spherical harmonic expansions
  call set_lat_lon
  call set_SHTOOLS(lmax)

  ! expand the sea level
  rfilt = 0.0_dp
  rdate = 00.0_dp
  call sea_level_ice5g(io1,rfilt,rdate,sl_lm)


 
  
  ! calculate the sea level on an equally spaced grid
  call set_flm_SH(sl_lm)
  call MakeGridDH_wrapper(sl_grid)

  ! write out the full data file
  open(io1,file='sphexp.data.full')
  write(io1,*) nlat_SH-1,nlon_SH-2
  do ilat = 1,nlat_SH-1
     lat = lat_SH(ilat)
     do ilon = 1,nlon_SH-2
        lon = lon_SH(ilon)
        write(io1,*) lat,lon,sl_grid(ilat,ilon)
     end do
  end do
  close(io1)


  call  sph_polar_grid(lmax)


  ! set up random number generator
  call date_and_time( values = time_info )
  msec = 1000 * time_info(7) + time_info(8)
  call ran_seed(sequence=msec)

  ! set variance of the data
  sigma_ran = 0.1_dp

  ! randomly sample the sea level to produce the data file 
  ! to be used in our test
  open(io1,file='sphexp.data')
  isample = 0
  do isample = 1,nsample

     ! get the position
     call random_number(ran)
     ilat = nlat_SH*ran
     if(ilat == 0) ilat = 1
     call random_number(ran)
     ilon = nlon_SH*ran    
     if(ilon == 0) ilon = 1
     if(ilon > nlon_SH-2) ilon = nlon_SH-2
     lat = lat_SH(ilat)
     lon = lon_SH(ilon)


     ! get the error term
     call gasdev(ran)
     ran = ran-0.5_dp
     sigma = sigma0*(1.0_dp+ran)

     ! set the data point
     call gasdev(ran)
     ran = ran-0.5_dp
     fun = sl_grid(ilat,ilon) + sigma*ran

     ! write out the data point
     if(sigma0 > 0.0_dp) then
        write(io1,*) lat,lon,fun,sigma
     else
        write(io1,*) lat,lon,fun      
     end if

  end do
  close(io1)



end program sphexp_data
