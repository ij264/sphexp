module module_ice5g

  use nrtype
  implicit none





  integer(i4b), parameter, private :: nlat_ice5g = 180
  integer(i4b), parameter, private :: nlon_ice5g = 360

  real(dp), dimension(nlat_ice5g), save, private :: lat_ice5g
  real(dp), dimension(nlon_ice5g), save, private :: lon_ice5g


  type ice5g_slice
     real(dp) :: time
     real(dp), dimension(nlat_ice5g,nlon_ice5g) :: orog
     real(dp), dimension(nlat_ice5g,nlon_ice5g) :: iceh
     real(dp), dimension(nlat_ice5g,nlon_ice5g) :: icem     
     integer(i4b) :: lmax
     complex(dpc), dimension(:), allocatable :: orog_lm
     complex(dpc), dimension(:), allocatable :: iceh_lm
     complex(dpc), dimension(:), allocatable :: icem_lm
  end type ice5g_slice

  ! dates of the ice5g time-slices
  integer(i4b), parameter :: ndate_ice5g = 39
  real(dp), dimension(ndate_ice5g) :: rdate_ice5g = (/00.0_dp, & 
                                                      00.5_dp, &
                                                      01.0_dp, &
                                                      01.5_dp, & 
                                                      02.0_dp, & 
                                                      02.5_dp, & 
                                                      03.0_dp, &
                                                      03.5_dp, & 
                                                      04.0_dp, &
                                                      04.5_dp, &
                                                      05.0_dp, &
                                                      05.5_dp, & 
                                                      06.0_dp, &
                                                      06.5_dp, &
                                                      07.0_dp, &
                                                      07.5_dp, &
                                                      08.0_dp, &
                                                      08.5_dp, & 
                                                      09.0_dp, &
                                                      09.5_dp, &
                                                      10.0_dp, & 
                                                      10.5_dp, &
                                                      11.0_dp, & 
                                                      11.5_dp, & 
                                                      12.0_dp, & 
                                                      12.5_dp, &
                                                      13.0_dp, & 
                                                      13.5_dp, & 
                                                      14.0_dp, & 
                                                      14.5_dp, & 
                                                      15.0_dp, & 
                                                      15.5_dp, & 
                                                      16.0_dp, & 
                                                      16.5_dp, & 
                                                      17.0_dp, & 
                                                      18.0_dp, & 
                                                      19.0_dp, & 
                                                      20.0_dp, &
                                                      21.0_dp/)


  character(len=4), dimension(ndate_ice5g) :: sdate_ice5g = (/'00.0', & 
                                                              '00.5', &
                                                              '01.0', &
                                                              '01.5', & 
                                                              '02.0', & 
                                                              '02.5', & 
                                                              '03.0', &
                                                              '03.5', & 
                                                              '04.0', &
                                                              '04.5', &
                                                              '05.0', &
                                                              '05.5', & 
                                                              '06.0', &
                                                              '06.5', &
                                                              '07.0', &
                                                              '07.5', &
                                                              '08.0', &
                                                              '08.5', & 
                                                              '09.0', &
                                                              '09.5', &
                                                              '10.0', & 
                                                              '10.5', &
                                                              '11.0', & 
                                                              '11.5', & 
                                                              '12.0', & 
                                                              '12.5', &
                                                              '13.0', & 
                                                              '13.5', & 
                                                              '14.0', & 
                                                              '14.5', & 
                                                              '15.0', & 
                                                              '15.5', & 
                                                              '16.0', & 
                                                              '16.5', & 
                                                              '17.0', & 
                                                              '18.0', & 
                                                              '19.0', & 
                                                              '20.0', &
                                                              '21.0'/)

  

contains


  subroutine sea_level_ice5g(io,rfilt,rdate,sl_lm)
    use nrtype
    use module_sphexp
    integer(i4b), intent(in) :: io
    real(dp), intent(in) :: rfilt
    real(dp), intent(in) :: rdate
    complex(dpc), dimension(:), allocatable, intent(inout) :: sl_lm

    integer(i4b) :: idate1,idate2,idate,ilat,ilon
    real(dp) :: time,time1,time2,lat,lon,sl1,sl2,fac
    type(ice5g_slice) :: slice1
    type(ice5g_slice) :: slice2

    ! work out the time
    time = t_present-rdate

    ! work out where this date lies in the ice5g lists
    if(rdate <= rdate_ice5g(1)) then
       idate1 = 1
       idate2 = 1
    else if(rdate >= rdate_ice5g(ndate_ice5g)) then
       idate1 = ndate_ice5g
       idate2 = ndate_ice5g
    else       
       do idate = 1,ndate_ice5g-1
          if(rdate_ice5g(idate) < rdate .and.  rdate <= rdate_ice5g(idate+1)) then
             idate1 = idate
             idate2 = idate+1
          end if
       end do
    end if

    if(idate1 == idate2) then

       ! read in the one time slice 
       call read_time_slice(io,sdate_ice5g(idate1),slice1)
       
       ! set the initial ice height
       do ilat = 1,nlat_SH
          lat = lat_SH(ilat)
          do ilon = 1,nlon_SH
             lon = lon_SH(ilon)
             fgrid_SH(ilat,ilon) = -ice5g_eval(slice1%orog,lat,lon)
          end do
       end do
       
       ! expand in spherical harmonics
       call SHExpandDH_wrapper(sl_lm,rfilt = rfilt) 


    else

       ! read in the two time slices
       call read_time_slice(io,sdate_ice5g(idate1),slice1)
       call read_time_slice(io,sdate_ice5g(idate2),slice2)

       time1 = slice1%time
       time2 = slice2%time
       fac = (time-time2)/(time1-time2)
       
       ! set the initial ice height
       do ilat = 1,nlat_SH
          lat = lat_SH(ilat)
          do ilon = 1,nlon_SH
             lon = lon_SH(ilon)
             sl1 = -ice5g_eval(slice1%orog,lat,lon)
             sl2 = -ice5g_eval(slice2%orog,lat,lon)
             fgrid_SH(ilat,ilon) = sl2 + fac*(sl1-sl2)
          end do
       end do
       
       ! expand in spherical harmonics
       call SHExpandDH_wrapper(sl_lm,rfilt = rfilt) 

    end if


    return
  end subroutine sea_level_ice5g



  subroutine ice5g_sphexp(lmax,slice)
    use nrtype
    use module_sphexp
    implicit none
    integer(i4b), intent(in) :: lmax
    type(ice5g_slice), intent(inout) :: slice


    integer(i4b) :: irec
    real(dp) :: lat,lon
    real(dp), dimension(:), allocatable :: f_sph


    ! build up orog array
    allocate(f_sph(nrec_sph))
    do irec = 1,nrec_sph       
       lat = 90.0_dp-tha_sph(irec)*rad2deg
       lon = pha_sph(irec)*rad2deg
       f_sph(irec) = ice5g_eval(slice%orog,lat,lon)       
    end do
    
    ! calculate orog coefficients
    call sphexp_cal(f_sph,slice%orog_lm)
    
    ! build up iceh array
    do irec = 1,nrec_sph       
       lat = 90.0_dp-tha_sph(irec)*rad2deg
       lon = pha_sph(irec)*rad2deg
       f_sph(irec) = ice5g_eval(slice%iceh,lat,lon)       
    end do
    
    ! calculate orog coefficients
    call sphexp_cal(f_sph,slice%iceh_lm)


    ! build up icem array
    do irec = 1,nrec_sph       
       lat = 90.0_dp-tha_sph(irec)*rad2deg
       lon = pha_sph(irec)*rad2deg
       f_sph(irec) = ice5g_eval(slice%icem,lat,lon)       
    end do
    
    ! calculate orog coefficients
    call sphexp_cal(f_sph,slice%icem_lm)


    return
  end subroutine ice5g_sphexp
    
  

  
  subroutine set_lat_lon
    use nrtype
    implicit none

    integer(i4b) :: ilat,ilon
    
    do ilat = 1,nlat_ice5g
       lat_ice5g(ilat) = -90.5_dp + ilat
    end do
    do ilon = 1,nlon_ice5g
       lon_ice5g(ilon) = -1.0_dp + ilon
    end do

    return
  end subroutine set_lat_lon


    
  subroutine read_time_slice(io1,date,slice)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: io1
    character(len=4), intent(in) :: date
    type(ice5g_slice), intent(inout) :: slice
    
    character(len=256) :: file_pref
    character(len=256) :: file

    integer(i4b) :: ilat,ilon,k

    ! set slice time (relative to present
    read(date,*) slice%time
    slice%time = t_present-slice%time*1000.0_dp*yr2sec

    file_pref = '/home/da380/raid/dta/ice5g/ascii/ice5g_v1.2_'//date 

    file = trim(file_pref)//'k_1deg.orog.ascii'
    open(io1,file=trim(file),action='read')
    k = 0
    do ilat = 1,nlat_ice5g
       do ilon = 1,nlon_ice5g
          k = k+1
          read(io1,*) slice%orog(ilat,ilon)
       end do
    end do
    close(io1)
      
    file = trim(file_pref)//'k_1deg.sftgit.ascii'
    open(io1,file=trim(file),action='read')
    k = 0
    do ilat = 1,nlat_ice5g
       do ilon = 1,nlon_ice5g
          k = k+1
          read(io1,*) slice%iceh(ilat,ilon)
       end do
    end do
    close(io1)
    
    file = trim(file_pref)//'k_1deg.sftgif.ascii'
    open(io1,file=trim(file),action='read')
    k = 0
    do ilat = 1,nlat_ice5g
       do ilon = 1,nlon_ice5g
          k = k+1
          read(io1,*) slice%icem(ilat,ilon)
          slice%icem(ilat,ilon) = slice%icem(ilat,ilon)/100.0_dp
       end do
    end do
    close(io1)
    

    return
  end subroutine read_time_slice

  
  function ice5g_eval(data,lat,lon)
    use nrtype
    implicit none
    real(dp) :: ice5g_eval
    real(dp), dimension(nlat_ice5g,nlon_ice5g), intent(in) :: data
    real(dp), intent(in) :: lat
    real(dp), intent(in) :: lon
    
    integer(i4b) :: ilat1,ilat2,ilon1,ilon2

    real(dp) :: lat1,lat2,lon1,lon2,dlat,dlon, & 
         rlat,rlon,latf,lonf,f,f1,f2,f11,f12,f21,f22


    ! set latitude range
    lat1 = lat_ice5g(1)
    lat2 = lat_ice5g(nlat_ice5g)
    dlat = lat_ice5g(2)-lat_ice5g(1)

    ! set longitude range
    lon1 = lon_ice5g(1)
    lon2 = lon_ice5g(nlon_ice5g)
    dlon = lon_ice5g(2)-lon_ice5g(1)


    ! adjust latitude 
    if(lat >= -90.0_dp .and. lat < lat1) then
       rlat = -89.5_dp
    else if(lat > lat2 .and. lat <= 90.0_dp) then
       rlat = 89.5_dp
    else
       rlat = lat
    end if

    ! adjust longitude
    if(lon < 0.0_dp) then
       rlon = 360.0_dp+lon
    else
       rlon = lon
    end if

    ! find latitude indices
    if(rlat == lat1) then
       ilat1 = 1
       ilat2 = 1
       latf  = 0.0_dp
    elseif(rlat == lat2) then
       ilat1 = nlat_ice5g
       ilat2 = nlat_ice5g
       latf  = 0.0_dp
    else
       ilat1 = floor((rlat-lat1)/dlat)+1
       ilat2 = ilat1+1
       latf  = (rlat-lat_ice5g(ilat1))/dlat
    end if

    ! find longitude indices
    if(rlon >= lon2) then
       ilon1 = nlon_ice5g
       ilon2 = 1
       lonf  = (rlon-lon_ice5g(ilon1))/dlon
    else
       ilon1 = floor((rlon-lon1)/dlon)+1
       ilon2 = ilon1+1
       lonf = (rlon-lon_ice5g(ilon1))/dlon
    end if


    ! get function values at nodes
    f11 = data(ilat1,ilon1)
    f21 = data(ilat2,ilon1)
    f12 = data(ilat1,ilon2)
    f22 = data(ilat2,ilon2)


    ! perform latitude interpolations
    f1 = f11 + (f21-f11)*latf
    f2 = f12 + (f22-f12)*latf

    ! perform longitude interpolation
    f = f1 + (f2-f1)*lonf

    ice5g_eval = 0.5_dp*f


    ! perform longitude interpolations
    f1 = f11 + (f12-f11)*lonf
    f2 = f21 + (f22-f21)*lonf

    ! perform longitude interpolation
    f = f1 + (f2-f1)*latf

    ice5g_eval = ice5g_eval + 0.5_dp*f



    return
  end function ice5g_eval


  subroutine deallocate_slice(slice)
    use nrtype
    implicit none
    type(ice5g_slice), intent(inout) :: slice
    
    if(allocated(slice%orog_lm)) deallocate(slice%orog_lm)
    if(allocated(slice%iceh_lm)) deallocate(slice%iceh_lm)
    if(allocated(slice%icem_lm)) deallocate(slice%icem_lm)
    

    return
  end subroutine deallocate_slice




end module module_ice5g
