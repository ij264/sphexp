module module_maps

  use nrtype
  implicit none

  
  real(dp), parameter :: r_nic  = 1221500.0_dp
  real(dp), parameter :: r_noc  = 3480000.0_dp
  real(dp), parameter :: r_moho = 6251000.0_dp

contains


  subroutine write_coast(io1,proj,lateq,loneq,azi,lat1_in,lat2_in, & 
       lon1_in,lon2_in,out_file,rotb_in)
    use nrtype
    implicit none

    integer(i4b), intent(in) :: io1
    character(len=*), intent(in) :: proj
    real(dp), intent(in) :: lateq
    real(dp), intent(in) :: loneq
    real(dp), intent(in) :: azi
    real(dp), intent(in) :: lat1_in
    real(dp), intent(in) :: lat2_in
    real(dp), intent(in) :: lon1_in
    real(dp), intent(in) :: lon2_in
    character(len=*), intent(in) :: out_file
    integer(i4b), intent(in), optional :: rotb_in

    logical(lgt), dimension(:), allocatable :: ifplot
    integer(i4b) :: ncoast,ios,icoast,nl,il,rotb
    real(dp) :: lat,lon,x,y,dist,lon0,dlat,dlon, & 
                lon1,lon2,lat1,lat2,lat0
    real(dp), dimension(:), allocatable :: xa,ya

    if(present(rotb_in)) then
       rotb = rotb_in
    else
       rotb = 0
    end if

    !----------------------------!
    !     read in coastlines     !
    !----------------------------!
    open(io1,file='/home/da380/raid/dta/coast.dat',action='read')
    
    ncoast = 0
    do 
       read(io1,*,iostat = ios) lat,lon
       if(ios /= 0) exit
       ncoast = ncoast+1
    end do
    rewind(io1)

    lat1 = lat1_in
    lat2 = lat2_in
    lon1 = lon1_in
    lon2 = lon2_in
    if(lon2 ==  180.0_dp) lon2 =  179.9_dp
    if(lon1 == -180.0_dp) lon1 = -179.9_dp
    
    allocate(xa(ncoast),ya(ncoast),ifplot(ncoast))
    ifplot = .true.
    ncoast = 0
    do 
       if(ncoast > 0) lon0 = lon
       if(ncoast > 0) lat0 = lat
       read(io1,*,iostat = ios) lat,lon
       if(ios /= 0) exit          
       ncoast = ncoast+1
       if(lat /= 999) then
          if(lat < lat1 .or. lat > lat2) then
             ifplot(ncoast) = .false.
          end if
          if(lon < lon1 .or. lon > lon2) then
             ifplot(ncoast) = .false.
          end if
          call rotate_to_equator(lat,lon,lateq,loneq,azi,lat,lon)
          if(ncoast == 1 .or. lon0 == 999) lon0 = lon
          if(ncoast == 1 .or. lat0 == 999) lat0 = lat
          if(abs(lon0-lon) < 350.0_dp) then 
             if(trim(proj) == 'ait') then
                call pr_aitoff(lat*deg2rad,lon*deg2rad,x,y)
             else if(trim(proj) == 'geo') then
                x = lon
                y = lat
             else if(trim(proj) == 'hammer') then
                call pr_hammer_aitoff(lat*deg2rad,lon*deg2rad,x,y)
             end if
             xa(ncoast) = x
             ya(ncoast) = y
          else
             xa(ncoast) = 999
             ya(ncoast) = 999
          end if
       else
          xa(ncoast) = 999
          ya(ncoast) = 999
       end if
    end do
    close(io1)


    open(io1,file=trim(out_file),form='formatted')
    do icoast = 1,ncoast
       if(ifplot(icoast)) then
          if(xa(icoast) /= 999) then
             write(io1,*) xa(icoast),ya(icoast)
          else
             write(io1,*) 'NaN      NaN'
          end if
       else
          write(io1,*) 'NaN      NaN'
       end if
    end do
    close(io1)

    
    nl = 180
    dlon = (lon2-lon1)/(nl-1)
    dlat = (lat2-lat1)/(nl-1)
    open(io1,file=trim(out_file)//'.border',form='formatted')
    do il = 1,nl
       lat = lat1+(il-1)*dlat
       lon = lon1
       if(rotb == 1) call rotate_to_equator(lat,lon,lateq,loneq,azi,lat,lon)
       if(trim(proj) == 'ait') then
          call pr_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       else if(trim(proj) == 'geo') then
          x = lon
          y = lat
       else if(trim(proj) == 'hammer') then
          call pr_hammer_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       end if
       write(io1,*) x,y
    end do
    do il = 1,nl
       lat = lat2
       lon = lon1+(il-1)*dlon
       if(rotb == 1) call rotate_to_equator(lat,lon,lateq,loneq,azi,lat,lon)
       if(trim(proj) == 'ait') then
          call pr_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       else if(trim(proj) == 'geo') then
          x = lon
          y = lat
       else if(trim(proj) == 'hammer') then
          call pr_hammer_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       end if
       write(io1,*) x,y
    end do
    do il = 1,nl
       lat = lat2-(il-1)*dlat
       lon = lon2
       if(rotb == 1) call rotate_to_equator(lat,lon,lateq,loneq,azi,lat,lon)
       if(trim(proj) == 'ait') then
          call pr_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       else if(trim(proj) == 'geo') then
          x = lon
          y = lat
       else if(trim(proj) == 'hammer') then
          call pr_hammer_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       end if
       write(io1,*) x,y
    end do
    do il = 1,nl
       lat = lat1
       lon = lon2-(il-1)*dlon
       if(rotb == 1) call rotate_to_equator(lat,lon,lateq,loneq,azi,lat,lon)
       if(trim(proj) == 'ait') then
          call pr_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       else if(trim(proj) == 'geo') then
          x = lon
          y = lat
       else if(trim(proj) == 'hammer') then
          call pr_hammer_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       end if
       write(io1,*) x,y
    end do
    close(io1)

    
    return
  end subroutine write_coast
  


  subroutine write_arc(io1,proj,lateq,loneq,azi,lat1_in,lat2_in, & 
       lon1_in,lon2_in,ddel_in,out_file,rotb_in)
    use nrtype
    implicit none  

    integer(i4b), intent(in) :: io1
    character(len=*), intent(in) :: proj
    real(dp), intent(in) :: lateq
    real(dp), intent(in) :: loneq
    real(dp), intent(in) :: azi
    real(dp), intent(in) :: lat1_in
    real(dp), intent(in) :: lat2_in
    real(dp), intent(in) :: lon1_in
    real(dp), intent(in) :: lon2_in
    real(dp), intent(in) :: ddel_in
    character(len=*), intent(in) :: out_file
    integer(i4b), intent(in), optional :: rotb_in

    integer(i4b) :: idel,ndel
    real(dp) :: th1,th2,ph1,ph2,ddel,lat,lon, & 
         x,y,th,ph,lon0,delta,del

    ! convert to radians
    th1 = (90.0_dp-lat1_in)*deg2rad
    th2 = (90.0_dp-lat2_in)*deg2rad
    ph1 = lon1_in*deg2rad
    ph2 = lon2_in*deg2rad
    ddel = ddel_in*deg2rad
    
    ! calculate arc angle
    call delta_cal(th1,ph1,th2,ph2,delta)
    
    if(ddel == 0.0_dp) then
       ddel = 1.0_dp*deg2rad
       ndel = twopi_d/ddel+1
       ddel = twopi_d/(ndel-1)
    else
       ndel = delta/ddel+1
       ddel = delta/(ndel-1)
    end if
    
    open(io1,file=trim(out_file),form='formatted',action='write')
    
    do idel = 1,ndel
       del = (idel-1)*ddel
       call great_circle(th1,ph1,th2,ph2,del,th,ph)
       lat = 90.0_dp-th*rad2deg
       lon = ph*rad2deg
       if(lon > 180.0_dp) lon = lon-360.0_dp
       call rotate_to_equator(lat,lon,lateq,loneq,azi,lat,lon)
       if(idel == 1) lon0 = 0.0_dp
       if(abs(lon-lon0) > 350.0_dp) then
          write(io1,*) 'NaN  NaN'
       end if
       if(trim(proj) == 'ait') then
          call pr_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       else if(trim(proj) == 'geo') then
          x = lon
          y = lat
       else if(trim(proj) == 'hammer') then
          call pr_hammer_aitoff(lat*deg2rad,lon*deg2rad,x,y)
       end if
       write(io1,*) x,y
       lon0 = lon
    end do
    close(io1)
    

    return
  end subroutine write_arc


  subroutine write_slice_border(io1,proj,lat1_in,lat2_in, & 
       lon1_in,lon2_in,ddel_in,r1_in,r2_in,out_file)
    use nrtype
!    use module_model
    implicit none  

    integer(i4b), intent(in) :: io1
    character(len=*), intent(in) :: proj
    real(dp), intent(in) :: lat1_in
    real(dp), intent(in) :: lat2_in
    real(dp), intent(in) :: lon1_in
    real(dp), intent(in) :: lon2_in
    real(dp), intent(in) :: ddel_in
    real(dp), intent(in) :: r1_in
    real(dp), intent(in) :: r2_in
    character(len=*), intent(in) :: out_file

    integer(i4b) :: idel,ndel,nr,ir,ivert,irplot
    real(dp) :: th1,th2,ph1,ph2,ddel,lat,lon, & 
         x,y,z,th,ph,lon0,delta,del,r1,r2,rr,dr

    real(dp), dimension(3) :: e1,e2,e3,f1,f2,f3

    ! convert to radians
    th1 = (90.0_dp-lat1_in)*deg2rad
    th2 = (90.0_dp-lat2_in)*deg2rad
    ph1 = lon1_in*deg2rad
    ph2 = lon2_in*deg2rad
    ddel = ddel_in*deg2rad


    ! convert to km
    r1 = r1_in*1000.0_dp
    r2 = r2_in*1000.0_dp
    
    ! calculate arc angle
    call delta_cal(th1,ph1,th2,ph2,delta)
    
    if(ddel == 0.0_dp) then
       ddel = 1.0_dp*deg2rad
       ndel = twopi_d/ddel+1
       ddel = twopi_d/(ndel-1)
       ivert = 0
    else
       ndel = delta/ddel+1
       ddel = delta/(ndel-1)
       ivert = 1
    end if
    
    open(io1,file=trim(out_file),form='formatted',action='write')
    

    !------------------------------------!
    !       write the radial lines       !
    !------------------------------------!


    if(ivert == 1) then

       del = 0.0_dp
1000   continue
       nr = 2
       dr = (r2-r1)/(nr-1)
       call great_circle(th1,ph1,th2,ph2,del,th,ph)      
       
       do ir = 1,nr

          rr = r1+(ir-1)*dr
          
          if(trim(proj) == 'geo') then
            
             ! set the local co-ordinates
             x = ph*rad2deg
             y = 90.0_dp-th*rad2deg
             z = rr
             if(x > 180.0_dp) then
                x = x-360.0_dp
             end if
             
          else if(trim(proj) == 'cart') then
             
             ! get local basis vectors
             call sph_unvec(th,ph,e1,e2,e3)
             
             ! set local co-ordinates
             x = rr*e1(1)
             y = rr*e1(2)
             z = rr*e1(3)
             
          else if(trim(proj) == 'proj') then
             
             ! set local co-ordinates
             x = rr*sin(del-0.5_dp*delta)
             y = rr*cos(del-0.5_dp*delta)
             z = 0.0_dp
             
          else
             
             stop 'write_slice_border: bad proj'
             
          end if
          
          write(io1,*) x,y,z

       end do


       write(io1,*) 'NaN     NaN      NaN'

       if(del /= delta) then
          del = delta
          goto 1000
       end if

    end if


    !---------------------------------------------!
    !            write the circlar arcs           !
    !---------------------------------------------!

    ! set the radius
    irplot = 1
    rr = r1

1001 continue 

    do idel = 1,ndel
       del = (idel-1)*ddel
       call great_circle(th1,ph1,th2,ph2,del,th,ph)              

       if(trim(proj) == 'geo') then
          
          ! set the local co-ordinates
          x = ph*rad2deg
          y = 90.0_dp-th*rad2deg
          z = rr
          if(x > 180.0_dp) then
             x = x-360.0_dp
          end if
          
       else if(trim(proj) == 'cart') then
          
          ! get local basis vectors
          call sph_unvec(th,ph,e1,e2,e3)
          
          ! set local co-ordinates
          x = rr*e1(1)
          y = rr*e1(2)
          z = rr*e1(3)
             
       else if(trim(proj) == 'proj') then
          
          ! set local co-ordinates
          x = rr*sin(del-0.5_dp*delta)
          y = rr*cos(del-0.5_dp*delta)
          z = 0.0_dp
          
       else
          
          stop 'write_slice_border: bad proj'
          
       end if

       write(io1,*) x,y,z

    end do

    write(io1,*) 'NaN     NaN      NaN'


    if(irplot == 1) then
       rr = r_nic
       if(rr < r1 .or. rr > r2) then
          irplot = irplot+1
       else
          irplot = irplot+1
          goto 1001
       end if
    end if


    if(irplot == 2) then
       rr = r_noc
       if(rr < r1 .or. rr > r2) then
          irplot = irplot+1
       else
          irplot = irplot+1
          goto 1001
       end if
    end if

    if(irplot == 3) then
       rr = 6371000.0_dp-670000.0_dp
       if(rr < r1 .or. rr > r2) then
          irplot = irplot+1
       else
          irplot = irplot+1
          goto 1001
       end if
    end if

    if(irplot == 4) then
       rr = r_moho
       if(rr < r1 .or. rr > r2) then
          irplot = irplot+1
       else
          irplot = irplot+1
          goto 1001
       end if
    end if


    if(irplot == 5) then
       irplot = irplot+1
       rr = r2       
       goto 1001
    end if


    ! close the border file
    close(io1)
    

    return
  end subroutine write_slice_border



  !===========================================================!
  !===========================================================!
  !     routines for calculation of distances on sphere       !
  !===========================================================!
  !===========================================================!


  subroutine angles(lats,lons,latr,lonr,delta,azep,azst,corr)
    ! this routine computes the epicentral angle and azimuth 
    ! from a given source to a given receiver. the inputs are
    ! the longitude and latitude of the two locations in degress.
    ! note that the answer is returned in radians, and that the 
    ! azimuth is measures anti-clockwise from south.
    
    ! if the optional input 'corr' is set to .true. that an 
    ! ellipticity correction is made 

    use nrtype
    implicit none
    real(dp), intent(in) :: latr,lonr,lats,lons
    real(dp), intent(out) :: delta,azep,azst
    logical(lgt), intent(in), optional :: corr
    integer(i4b) :: i,nr
    real(dp) :: conv,sgn
    real(dp) :: cltr,clts,clnr,clns,sltr,slnr,slts,slns,saz,sazr
    real(dp), dimension(3) :: rs,rr,ts,ps,ks,kr,tr,pr
      
    ! ellipticity correction factor = 1 - 2 / 298.3
    real(dp), parameter :: flat = 0.9932953402614817_dp


    conv=pi_d/180.0_dp
    
    ! convert angles to radians
    cltr=(90.0_dp-latr)*conv
    clts=(90.0_dp-lats)*conv
    clnr=lonr*conv
    clns=lons*conv

    
    ! compute the required sines and cosines
    sltr=sin(cltr); slts=sin(clts); cltr=cos(cltr); clts=cos(clts)
    slnr=sin(clnr); slns=sin(clns); clnr=cos(clnr); clns=cos(clns)
    
    ! compute the required unit vectors
    rs(1)=slts*clns; rs(2)=slts*slns; rs(3)=clts
    rr(1)=sltr*clnr; rr(2)=sltr*slnr; rr(3)=cltr
    ts(1)=clts*clns; ts(2)=clts*slns; ts(3)=-slts
    tr(1)=cltr*clnr; tr(2)=cltr*slnr; tr(3)=-sltr
    ps(1)=-slns;     ps(2)=clns;      ps(3)=0;
    pr(1)=-slnr;     pr(2)=clnr;      pr(3)=0;


    ! compute cos(delta)
    delta=dot_product(rr,rs)      
    ! compute the tangent vector to the source-recevier
    ! great circle a the source
    if(abs(delta) < 1.0_dp) then
       ks=(rr-delta*rs)/sqrt(1-delta*delta)
    else
       ks=0.0_dp
    end if
    if(abs(delta) < 1.0_dp) then
       kr=(delta*rr-rs)/sqrt(1-delta*delta)
    else
       kr=0.0_dp
    end if
    !  compute the sine and consine of the azimuth
    azep=dot_product(ks,ts)
    saz=dot_product(ks,ps)
    azst=dot_product(kr,tr)
    sazr=dot_product(kr,pr)
    ! check the deltas
    if(delta > 1.0_dp) delta=1.0_dp
    ! check the azimuths
    if(abs(azep) > 1.0_dp) azep=sign(1.0_dp,azep)
    if(abs(azst) > 1.0_dp) azst=sign(1.0_dp,azst)
    ! compute delta
    delta=acos(delta)
    ! compute azimuth picking the right quadrant
    if(saz >= 0.0_dp) then
       azep=acos(azep)
    else
       azep=acos(azep)
       azep=twopi_d-azep
    end if
    if(sazr >= 0.0_dp) then
       azst=acos(azst)
    else
       azst=acos(azst)
       azst=twopi_d-azst
    end if
    

    return
  end subroutine angles


  subroutine delaz(eplati,eplong,stlati,stlong,delta,azep,azst,ellip)
    use nrtype
    implicit none
    ! in/out
    real(dp), intent(in) :: eplati,eplong,stlati,stlong
    real(dp), intent(out) :: delta,azep,azst
    logical(lgt), intent(in), optional :: ellip
    
    ! parameters
    real(dp), parameter :: geoco  = 0.993277_dp
    real(dp), parameter :: hpi    = 0.5_dp*pi_d
    real(dp), parameter :: rad    = 0.0174533_dp
    real(dp), parameter :: reprad = 57.29578_dp
    real(dp), parameter :: delsmall = 0.00001_dp
    
    ! local variables
    logical(lgt) :: ellip_do
    real(dp) :: el,stl,elon,slon,as,bs,cs,ds,a,b,c,d, & 
         cdel,delt,sdel,caze,aze,azs,dif,cazs,eplat,stlat
    
    
    if(present(ellip)) then
       ellip_do = ellip
    else
       ellip_do = .true.
    end if
    
    eplat = eplati
    stlat = stlati
    
      
    if(eplat == 90.0_dp) eplat = 89.999_dp
    if(stlat == 90.0_dp) stlat = 89.999_dp  !station at N pole

    if(eplat == -90.0_dp) eplat = -89.999_dp
    if(stlat == -90.0_dp) stlat = -89.999_dp  !station at S pole
      
    if(eplat == stlat .and. eplong == stlong) then
       delta = 0.0_dp
       azep  = 0.0_dp
       azst  = 0.0_dp
       return
    end if

      
    if(ellip_do) then
       el  = atan(geoco*tan(eplat*rad))
    else
       el = eplat*rad
    end if
    el  = hpi-el
    if(ellip_do) then
       stl = atan(geoco*tan(stlat*rad))
    else
       stl = stlat*rad
    end if
    
    stl = hpi-stl
    
    elon = eplong*rad
    slon = stlong*rad
    
    as = cos(stl)
    bs = sin(stl)
    cs = cos(slon)
    ds = sin(slon)
    
    a = cos(el)
    b = sin(el)
    c = cos(elon)
    d = sin(elon)
    
    cdel = a*as+b*bs*(c*cs+d*ds)
    if(abs(cdel) > 1.) cdel = sign(1.0_dp,cdel)
    delt  = acos(cdel)
    delta = delt
    

    sdel = sin(delt)
    caze = (as-a*cdel)/(sdel*b)
    if(abs(caze) > 1.0_dp) caze = sign(1.0_dp,caze)
    aze = acos(caze)
    
    if(bs > 0.0_dp) cazs = (a-as*cdel)/(bs*sdel)
    if(bs == 0.0_dp) cazs = sign(1.0_dp,cazs)
    if(abs(cazs) > 1.0_dp) cazs = sign(1.0_dp,cazs)
    azs = acos(cazs)
    dif = ds*c-cs*d
    if(dif < 0.0_dp) aze = twopi_d-aze
    azep = aze
    if(dif > 0.0_dp) azs = twopi_d-azs
    azst = azs
    

    return
  end subroutine delaz



  
  subroutine delta_cal(th1,ph1,th2,ph2,delta)
    use nrtype
    implicit none
    real(dp), intent(in) :: th1
    real(dp), intent(in) :: ph1
    real(dp), intent(in) :: th2
    real(dp), intent(in) :: ph2
    real(dp), intent(out) :: delta


    real(dp), dimension(3) :: rr1,rr2

    call sph_2_cart(1.0_dp,th1,ph1,rr1(1),rr1(2),rr1(3))
    call sph_2_cart(1.0_dp,th2,ph2,rr2(1),rr2(2),rr2(3))
    
    delta = dot_product(rr1,rr2)
    delta = acos(delta)

    return
  end subroutine delta_cal



  
  !===============================================================!
  !===============================================================!
  !             routines for calculation of rotations             !
  !===============================================================!
  !===============================================================!


  subroutine rotmat(alpha,beta,gamma,rmat)
    ! given the Euler angles (using the convention of Edmonds)
    ! this routine returns the rotation matrix which takes 
    ! the unprimed co-ordinates (x,y,z) into the primed 
    ! co-ordinates (x',y',z') according to 
    !
    ! (x',y',z') = rmat*(x,y,z)
    !
    ! the formula used is taken from (C.242) of Dahlen & Tromp (1998)
    use nrtype
    implicit none
    real(dp), intent(in) :: alpha,beta,gamma
    real(dp), dimension(3,3) :: rmat
    
    real(dp) :: sa,sb,sg,ca,cb,cg
    
    
    sa = sin(alpha)
    sb = sin(beta)
    sg = sin(gamma)
    
    ca = cos(alpha)
    cb = cos(beta)
    cg = cos(gamma)
    
    rmat(1,1) =  ca*cb*cg-sa*sg
    rmat(1,2) =  sa*cb*cg+ca*sg
    rmat(1,3) = -sb*cg
    rmat(2,1) = -ca*cb*sg-sa*cg
    rmat(2,2) = -sa*cb*sg+ca*cg
    rmat(2,3) =  sb*sg
    rmat(3,1) =  ca*sb
    rmat(3,2) =  sa*sb
    rmat(3,3) =  cb

    return
  end subroutine rotmat


  subroutine rotate_to_pole(lat,lon,alpha,beta,gamma)
    ! given the lat and lon of a point in the earth model, 
    ! this routine returns the Euler angles that rotates this point 
    ! the north pole of the new co-ordinate system and aligns the 
    ! great circle between the point and the geographic north pole
    ! with the prime merdian of the new co-ordinate system
    use nrtype
    implicit none
    real(dp), intent(in) :: lat,lon
    real(dp), intent(out) :: alpha,beta,gamma


    alpha = lon*deg2rad
    beta  = (90.0_dp-lat)*deg2rad
    gamma = 0.0_dp
    

      
    return
  end subroutine rotate_to_pole


    
  subroutine rotate_great_circle(lat1,lon1,lat2,lon2,alpha, & 
       beta,gamma,delout,azout)
    ! given the latitude and longitude of two points, this routine returns 
    ! the Euler angles that produce a new co-ordinate system in which 
    ! (lat1,lon1) lies at the north pole and (lat2,lon2) lies on the 
    ! prime meridian.
    use nrtype      
    implicit none
    real(dp), intent(in) :: lat1,lon1,lat2,lon2
    real(dp), intent(out) :: alpha,beta,gamma
    real(dp), intent(out), optional :: delout
    real(dp), intent(out), optional :: azout
    
    real(dp) :: delta,azst,azep,azimuth
    
    ! work out epicentral angle and azimuth from point 1 to point 2
    call delaz(lat1,lon1,lat2,lon2,delta,azep,azst,.false.)
    azimuth = pi_d-azep
    

    alpha = lon1*deg2rad
    beta  = (90.0_dp-lat1)*deg2rad
    gamma = azimuth
    
    if(present(delout)) delout = delta
    if(present(azout)) azout = azout
    
    return
  end subroutine rotate_great_circle  


  subroutine rotate_to_equator(lat,lon,lateq,loneq,azi,latr,lonr,rmt)
    use nrtype
    implicit none
    real(dp), intent(in) :: lat
    real(dp), intent(in) :: lon
    real(dp), intent(in) :: lateq
    real(dp), intent(in) :: loneq
    real(dp), intent(in) :: azi
    real(dp), intent(out) :: latr
    real(dp), intent(out) :: lonr
    real(dp), dimension(3,3), intent(out), optional :: rmt
    
    real(dp) :: r,th,ph,xi
    real(dp), dimension(3) :: x
    real(dp), dimension(3,3) :: rot

    r = 1.0_dp
    th = (90.0_dp-lat)*deg2rad
    ph = lon*deg2rad
    call sph_2_cart(1.0_dp,th,ph,x(1),x(2),x(3))
    
    ! perform the first rotation
    xi = loneq*deg2rad
    rot(:,:) =  0.0_dp
    rot(1,1) =  cos(xi) 
    rot(1,2) =  sin(xi)
    rot(2,1) = -sin(xi)
    rot(2,2) =  cos(xi)
    rot(3,3) =  1.0_dp
    x = matmul(rot,x)
    if(present(rmt)) rmt = rot


    ! perform the second rotation
    xi = lateq*deg2rad
    rot(:,:) =  0.0_dp
    rot(1,1) =  cos(xi) 
    rot(1,3) =  sin(xi)
    rot(3,1) = -sin(xi)
    rot(3,3) =  cos(xi)
    rot(2,2) =  1.0_dp    
    x = matmul(rot,x)
    if(present(rmt)) rmt = matmul(rot,rmt)


    ! perform the third rotation
    xi = azi*deg2rad
    rot(:,:) =  0.0_dp
    rot(1,1) =  1.0_dp
    rot(2,2) =  cos(xi)
    rot(2,3) =  sin(xi)
    rot(3,2) =  -sin(xi)
    rot(3,3) =  cos(xi)
    x = matmul(rot,x)
    if(present(rmt)) rmt = matmul(rot,rmt)

    ! transform back to spherical co-ordinates
    call cart_2_sph(x(1),x(2),x(3),r,th,ph)

    latr = 90.0_dp-th*rad2deg
    lonr = ph*rad2deg
    if(lonr > 180.0_dp) lonr = lonr-360.0_dp


    return
  end subroutine rotate_to_equator



  !===========================================================!
  !===========================================================!
  !      routines for calculation of great circular arcs      !
  !===========================================================!
  !===========================================================!


  subroutine great_circle(th1,ph1,th2,ph2,delta,th,ph)
    use nrtype 
    implicit none
    real(dp), intent(in) :: th1
    real(dp), intent(in) :: ph1
    real(dp), intent(in) :: th2
    real(dp), intent(in) :: ph2
    real(dp), intent(in) :: delta
    real(dp), intent(out) :: th
    real(dp), intent(out) :: ph


    real(dp) :: cdelta,sdelta,bdelta,cbdelta,sbdelta,r
    real(dp), dimension(3) :: rr1,rr2,rr

    call sph_2_cart(1.0_dp,th1,ph1,rr1(1),rr1(2),rr1(3))
    call sph_2_cart(1.0_dp,th2,ph2,rr2(1),rr2(2),rr2(3))

    cbdelta = dot_product(rr1,rr2)
    bdelta = acos(cbdelta)
    sbdelta = sin(bdelta)

    sdelta = sin(delta)
    cdelta = cos(delta)

    rr = cdelta*rr1 + sdelta*(rr2-cbdelta*rr1)/sbdelta
    call cart_2_sph(rr(1),rr(2),rr(3),r,th,ph)

    return
  end subroutine great_circle


  subroutine great_circle_bv(th1,ph1,th2,ph2,delta,th,ph,e1,e2,e3)
    use nrtype 
    implicit none
    real(dp), intent(in) :: th1
    real(dp), intent(in) :: ph1
    real(dp), intent(in) :: th2
    real(dp), intent(in) :: ph2
    real(dp), intent(in) :: delta
    real(dp), intent(out) :: th
    real(dp), intent(out) :: ph
    real(dp), dimension(3), intent(out) :: e1
    real(dp), dimension(3), intent(out) :: e2
    real(dp), dimension(3), intent(out) :: e3


    real(dp) :: cdelta,sdelta,bdelta,cbdelta,sbdelta,r
    real(dp), dimension(3) :: rr1,rr2,rr


    call sph_2_cart(1.0_dp,th1,ph1,rr1(1),rr1(2),rr1(3))
    call sph_2_cart(1.0_dp,th2,ph2,rr2(1),rr2(2),rr2(3))

    cbdelta = dot_product(rr1,rr2)
    bdelta = acos(cbdelta)
    sbdelta = sin(bdelta)

    sdelta = sin(delta)
    cdelta = cos(delta)
    
    ! compute the basis vectors
    e2 =  cdelta*rr1 + sdelta*(rr2-cbdelta*rr1)/sbdelta
    e1 = -sdelta*rr1 + cdelta*(rr2-cbdelta*rr1)/sbdelta
    e3(1) = rr1(2)*rr2(3)-rr1(3)*rr2(2)
    e3(2) = rr1(3)*rr2(1)-rr1(1)*rr2(3)
    e3(3) = rr1(1)*rr2(2)-rr1(2)*rr2(1)
    e3 = e3/sbdelta

    call cart_2_sph(e2(1),e2(2),e2(3),r,th,ph)


    return
  end subroutine great_circle_bv


  !=============================================================!
  !=============================================================!
  !  spherical polar to cartesian co-ordinate transformations   !
  !=============================================================!
  !=============================================================!

  subroutine sph_2_cart(r,theta,phi,x,y,z)
    use nrtype
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(in) :: theta
    real(dp), intent(in) :: phi
    real(dp), intent(out) :: x
    real(dp), intent(out) :: y
    real(dp), intent(out) :: z
    
    real(dp) :: cth,sth,cph,sph

    cth = cos(theta)
    sth = sin(theta)
    cph = cos(phi)
    sph = sin(phi)

    x = sth*cph
    y = sth*sph
    z = cth

    return
  end subroutine sph_2_cart


  subroutine cart_2_sph(x,y,z,r,theta,phi)
    use nrtype
    implicit none
    real(dp), intent(in) :: x,y,z
    real(dp), intent(out) :: r,theta,phi

    r=x*x+y*y+z*z
    r=sqrt(r)

    if(z == 0.0_dp) then
       theta=0.5_dp*pi_d
    else
       theta=z/r
       theta=acos(theta)
    end if

    phi = atan2(y,x)
    
       
    return
  end subroutine cart_2_sph


  subroutine sph_unvec(th,ph,rv,thv,phv)
    use nrtype
    implicit none
    real(dp), intent(in) :: th
    real(dp), intent(in) :: ph
    real(dp), dimension(3), intent(out) :: rv
    real(dp), dimension(3), intent(out) :: thv
    real(dp), dimension(3), intent(out) :: phv

    real(dp) :: sth,cth,sph,cph

    cth = cos(th)
    sth = sin(th)
    sph = sin(ph)
    cph = cos(ph)

    rv(1)  =  sth*cph
    rv(2)  =  sth*sph
    rv(3)  =  cth
    thv(1) =  cth*cph
    thv(2) =  cth*sph
    thv(3) = -sth
    phv(1) = -sph
    phv(2) =  cph
    phv(3) =  0.0_dp


    return
  end subroutine sph_unvec


  !=====================================================================!
  !=====================================================================!
  !                      map projection routines                        !
  !=====================================================================!
  !=====================================================================!

    
  subroutine pr_aitoff(lat,lon,x,y,e1,e2)
    ! routines performs Aitoff projection 
    ! input angles are in radians:
    !  -pi/2 <= lat <= pi/2
    !  -pi <= lon <= pi
    use nrtype
    implicit none
    real(dp), intent(in) :: lat
    real(dp), intent(in) :: lon
    real(dp), intent(out) :: x
    real(dp), intent(out) :: y
    real(dp), dimension(2), intent(out), optional  :: e1
    real(dp), dimension(2), intent(out), optional  :: e2

    
    real(dp) :: fac

    fac = cos(lat)*cos(0.5_dp*lon)
    fac = acos(fac)
    if(fac == 0.0_dp) then
       fac = 1.0_dp
    else
       fac = sin(fac)/fac
    end if
    fac = 1.0_dp/fac
    
    x = 2.0_dp*fac*cos(lat)*sin(0.5_dp*lon)
    y = fac*sin(lat)
    



    if(present(e1) .and. present(e2)) then
    end if

    return
  end subroutine pr_aitoff



    
  subroutine pr_hammer_aitoff(lat,lon,x,y,e1,e2)
    ! routines performs Aitoff projection 
    ! input angles are in radians:
    !  -pi/2 <= lat <= pi/2
    !  -pi <= lon <= pi
    use nrtype
    implicit none
    real(dp), intent(in) :: lat
    real(dp), intent(in) :: lon
    real(dp), intent(out) :: x
    real(dp), intent(out) :: y
    real(dp), dimension(2), intent(out), optional  :: e1
    real(dp), dimension(2), intent(out), optional  :: e2



    x = 2.0_dp*cos(lat)*sin(0.5_dp*lon) & 
         /sqrt(1+cos(lat)*cos(0.5_dp*lon))
    y = sin(lat)/sqrt(1+cos(lat)*cos(0.5_dp*lon))


    if(present(e1) .and. present(e2)) then
       e1(1) = sin(lat)*sin(0.5_dp*lon)*(2.0_dp+cos(lat)*cos(0.5_dp*lon)) & 
            / (1+cos(lat)*cos(0.5_dp*lon))**(1.5_dp)
       e1(2) = 0.25_dp*cos(lat)*(4.0_dp*cos(0.5_dp*lon) + cos(lat)*(3.0_dp & 
            + cos(lon)))/(1+cos(lat)*cos(0.5_dp*lon))**(1.5_dp)
       
       e2(1) = 0.25_dp*(4.0_dp*cos(lat) + (3.0_dp+cos(2.0_dp*lat)) & 
            *cos(0.5_dp*lon))/(1+cos(lat)*cos(0.5_dp*lon))**(1.5_dp)
       
       e2(2) = 0.25_dp*cos(lat)*sin(lat)*sin(0.5_dp*lon) & 
            / (1+cos(lat)*cos(0.5_dp*lon))**(1.5_dp)
    end if

    return
  end subroutine pr_hammer_aitoff

  



end module module_maps
