module module_crust


  use nrtype
  use module_util
  implicit none

  !-------------------------------------------------------!
  ! definitions of the arrays
  !
  !
  ! elev(:,:)       = top of the water
  ! amapl(1,:,:)    = bottom of the water
  ! amapl(2,:,:)    = bottom of ice
  ! amapl(3,:,:)    = bottom of soft sediments
  ! amapl(4,:,:)    = bottom of hard sediments
  ! amapl(5,:,:)    = bottom of upper crust
  ! amapl(6,:,:)    = bottom of middle crust
  ! amapt(:,:)      = bottom of the lower crust = MOHO
  ! amapth(:,:)     = total crustal thickness
  ! amaps(:,:)      = sediment thickness
  ! amapi(:,:)      = ice thickness
  ! amapvp(1,:,:)   = vp, water layer
  ! amapvp(2,:,:)   = vp, ice layer
  ! amapvp(3,:,:)   = vp, soft sediments 
  ! amapvp(4,:,:)   = vp, hard sediments 
  ! amapvp(5,:,:)   = vp, upper crust
  ! amapvp(6,:,:)   = vp, middle crust
  ! amapvp(7,:,:)   = vp, lower crust
  ! amapvp(8,:,:)   = vp, mantle below the moho
  !
  ! amapvs and amaprho are as for amapvp, but for s-velocity 
  ! and density respectively.  

  

  integer(i4b), parameter, private :: ntyp = 360
  integer(i4b), parameter, private :: nla = 90
  integer(i4b), parameter, private :: nlo = 180


  integer(i4b), dimension(nlo), private :: ilon
  integer(i4b), dimension(nla), private :: il
  real(dp), dimension(nlo,nla), private :: elev
  real(dp), dimension(nlo,nla), private :: amapt
  real(dp), dimension(nlo,nla), private :: amapth
  real(dp), dimension(nlo,nla), private :: amaps
  real(dp), dimension(nlo,nla), private :: amapi
  real(dp), dimension(6,nlo,nla), private :: amapl
  real(dp), dimension(8,nlo,nla), private :: amapvp
  real(dp), dimension(8,nlo,nla), private :: amaprho
  real(dp), dimension(8,nlo,nla), private :: amapvs


   contains


     subroutine crust_eval(lat,lon,felev,famapl,famapt, & 
          famapth,famaps,famapi,famapvp,famapvs,famaprho)
       ! given the latitude (lat) and longitude (lon), 
       ! this routine returns the values of the crustal 
       ! parameters using bilinear interpolation from the 
       ! CRUST 2.0 grid
       use nrtype
       implicit none
       real(dp), intent(in) :: lat,lon
       real(dp), intent(out), optional :: felev       
       real(dp), dimension(6),intent(out), optional :: famapl
       real(dp), intent(out), optional :: famapt
       real(dp), intent(out), optional :: famapth
       real(dp), intent(out), optional :: famaps
       real(dp), intent(out), optional :: famapi
       real(dp), dimension(8), intent(out), optional :: famapvp
       real(dp), dimension(8), intent(out), optional :: famapvs
       real(dp), dimension(8), intent(out), optional :: famaprho
       

       integer(i4b) :: jlat1,jlon1,jlat2,jlon2,k
       real(dp) :: x1,x2,y1,y2,f11,f12,f21,f22,lont

       ! work out which cell the point lies in
       
       if(lon < 0.0_dp) then
          lont = lon+360.0_dp
       else
          lont = lon
       end if

       if(lat > -88.0_dp) then
          jlat1 = floor((90.0_dp-lat)/2.0_dp)+1
          jlat2 = jlat1+1
          x1 = il(jlat1)
          x2 = il(jlat2)
       else
          jlat1 = nla
          jlat2 = nla
          x1 = -88.0_dp
          x2 = -90.0_dp
       end if

       if(lont < 358.0_dp) then
          jlon1 = floor(lont/2.0_dp)+1
          jlon2 = jlon1+1
          y1 = ilon(jlon1)
          y2 = ilon(jlon2)
       else
          jlon1 = nlo
          jlon2 = 1
          y1 = 358.0_dp
          y2 = 360.0_dp
       end if



       if(present(felev)) then
          f11 = elev(jlon1,jlat1)
          f12 = elev(jlon2,jlat1)
          f21 = elev(jlon1,jlat2)
          f22 = elev(jlon2,jlat2)
          felev = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont)          
       end if

       
       if(present(famapl)) then
          do k = 1,6
             f11 = amapl(k,jlon1,jlat1)
             f12 = amapl(k,jlon2,jlat1)
             f21 = amapl(k,jlon1,jlat2)
             f22 = amapl(k,jlon2,jlat2)
             famapl(k) = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont) 
          end do
       end if

       if(present(famapt)) then
          f11 = amapt(jlon1,jlat1)
          f12 = amapt(jlon2,jlat1)
          f21 = amapt(jlon1,jlat2)
          f22 = amapt(jlon2,jlat2)
          famapt = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont)          
       end if

       if(present(famapth)) then
          f11 = amapth(jlon1,jlat1)
          f12 = amapth(jlon2,jlat1)
          f21 = amapth(jlon1,jlat2)
          f22 = amapth(jlon2,jlat2)
          famapth = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont)          
       end if

       if(present(famaps)) then
          f11 = amaps(jlon1,jlat1)
          f12 = amaps(jlon2,jlat1)
          f21 = amaps(jlon1,jlat2)
          f22 = amaps(jlon2,jlat2)
          famaps = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont)          
       end if

       if(present(famapi)) then
          f11 = amapi(jlon1,jlat1)
          f12 = amapi(jlon2,jlat1)
          f21 = amapi(jlon1,jlat2)
          f22 = amapi(jlon2,jlat2)
          famapi = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont)          
       end if

       if(present(famapvp)) then
          do k = 1,8
             f11 = amapvp(k,jlon1,jlat1)
             f12 = amapvp(k,jlon2,jlat1)
             f21 = amapvp(k,jlon1,jlat2)
             f22 = amapvp(k,jlon2,jlat2)
             famapvp(k) = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont) 
          end do
       end if


       if(present(famapvs)) then
          do k = 1,8
             f11 = amapvs(k,jlon1,jlat1)
             f12 = amapvs(k,jlon2,jlat1)
             f21 = amapvs(k,jlon1,jlat2)
             f22 = amapvs(k,jlon2,jlat2)
             famapvs(k) = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont) 
          end do
       end if


       if(present(famaprho)) then
          do k = 1,8
             f11 = amaprho(k,jlon1,jlat1)
             f12 = amaprho(k,jlon2,jlat1)
             f21 = amaprho(k,jlon1,jlat2)
             f22 = amaprho(k,jlon2,jlat2)
             famaprho(k) = bilinear(x1,x2,y1,y2,f11,f12,f21,f22,lat,lont) 
          end do
       end if



       return
     end subroutine crust_eval





    subroutine get_crust
      use nrtype
      implicit none

      integer(i4b) :: i,j,k,ilat,ilong      
      integer(i4b), dimension(nlo) :: itmp
      real(dp), dimension(nlo,nla) :: atmp

      ! get the arrays
      call readCN2_t7_extract(ilon,il,elev,amapt,amapth,amaps,amapi,amapl)
      call readCN2_7vr_extract(amapvp,amaprho,amapvs)
      
      
      ! re-order so that longitude ranges from 0 to 360 instead of -180 to 180
      itmp(1:nlo/2) = ilon(nlo/2+1:nlo)
      itmp(nlo/2+1:nlo) = ilon(1:nlo/2)+360
      ilon = itmp
      
      ! perform the same re-ording for the other arrays
      atmp(1:nlo/2,:) = elev(nlo/2+1:nlo,:)
      atmp(nlo/2+1:nlo,:) = elev(1:nlo/2,:)
      elev = atmp
      
      atmp(1:nlo/2,:) = amapt(nlo/2+1:nlo,:)
      atmp(nlo/2+1:nlo,:) = amapt(1:nlo/2,:)
      amapt = atmp
      
      atmp(1:nlo/2,:) = amapth(nlo/2+1:nlo,:)
      atmp(nlo/2+1:nlo,:) = amapth(1:nlo/2,:)
      amapth = atmp
      
      atmp(1:nlo/2,:) = amaps(nlo/2+1:nlo,:)
      atmp(nlo/2+1:nlo,:) = amaps(1:nlo/2,:)
      amaps = atmp
      
      atmp(1:nlo/2,:) = amapi(nlo/2+1:nlo,:)
      atmp(nlo/2+1:nlo,:) = amapi(1:nlo/2,:)
      amapi = atmp
      
      do i = 1,6
         atmp(1:nlo/2,:) = amapl(i,nlo/2+1:nlo,:)
         atmp(nlo/2+1:nlo,:) = amapl(i,1:nlo/2,:)
         amapl(i,:,:) = atmp     
      end do
      
      
      do i = 1,8
         atmp(1:nlo/2,:) = amapvp(i,nlo/2+1:nlo,:)
         atmp(nlo/2+1:nlo,:) = amapvp(i,1:nlo/2,:)
         amapvp(i,:,:) = atmp     
         atmp(1:nlo/2,:) = amaprho(i,nlo/2+1:nlo,:)
         atmp(nlo/2+1:nlo,:) = amaprho(i,1:nlo/2,:)
         amaprho(i,:,:) = atmp 
         atmp(1:nlo/2,:) = amapvs(i,nlo/2+1:nlo,:)
         atmp(nlo/2+1:nlo,:) = amapvs(i,1:nlo/2,:)
         amapvs(i,:,:) = atmp 
      end do
      
      
      return
    end subroutine get_crust


    subroutine readCN2_t7_extract(ilon,il,elev,amapt,amapth,amaps,amapi,amapl)
      use nrtype
      implicit none
    
      integer(i4b), dimension(nlo), intent(out) :: ilon
      integer(i4b), dimension(nla), intent(out) :: il
      real(dp), dimension(nlo,nla), intent(out) :: elev
      real(dp), dimension(nlo,nla), intent(out) :: amapt
      real(dp), dimension(nlo,nla), intent(out) :: amapth
      real(dp), dimension(nlo,nla), intent(out) :: amaps
      real(dp), dimension(nlo,nla), intent(out) :: amapi
      real(dp), dimension(6,nlo,nla), intent(out) :: amapl



      character(len=2), dimension(ntyp) :: ctype
      character(len=1) :: dum
      character(len=5) :: dum0
      character(len=506) :: line
      character(len=2), dimension(nlo) :: types
      
      integer(i4b) :: i,j,k,l,m,ilat
      
      real(dp) :: aux,felev,fthick,flons
      real(dp), dimension(ntyp,7) :: flev
      real(dp), dimension(7) :: flev2

      
      
      ! read in elevation map      
      open(3,file='/home/da380/raid/dta/crust2/CNelevatio2.txt')
      open(19,file='map2.topo')
      read(3,*) ilon
      do j=1,nla
         read(3,*) il(j),(elev(i,j),i=1,nlo) 
      end do
      close(3)
      close(19)
      
      open(2,file='/home/da380/raid/dta/crust2/CNtype2_key.txt')
      open(7,file='/home/da380/raid/dta/crust2/CNtype2.txt')
            
      ! read in key for crust types      
      read(2,890)dum
      do  i=1,ntyp
         read(2,899)ctype(i)
         read(2,891)dum
         read(2,899)line
         read(line,*)(flev(i,l),l=1,7),dum0,fthick
         ! switch layer 1 and 2
         ! new layer 1 is water
         ! new layer 2 is ice
         aux=flev(i,2)
         flev(i,2)=flev(i,1)
         flev(i,1)=aux
         felev=0.
         do l=1,7
            felev=felev+flev(i,l)
         end do
         if(felev /= fthick) then
            print*,' read error in thickness of crust'
            print 905,' type: ',i,ctype(i)
            print*,felev,fthick 
            stop
         end if
      end do
      close(2)
      
      
    ! read CNtype file
      read(7,*)flons
      do  j=1,nla
         read(7,901) ilat,types
         outer:do   i=1,nlo
            do  l=1,ntyp
               if(types(i) == ctype(l)) then
                  do  m=1,7
                     flev2(m)=flev(l,m)
                  end do
                  amaps(i,j)=flev2(3)+flev2(4)
                  amapi(i,j)=flev2(2)
                  if(elev(i,j) < 0.0_dp)then
                     flev2(1) = -elev(i,j)/1000.0_dp
                     elev(i,j) = 0.0_dp      
                  end if
                  do k=2,7
                     flev2(k)=flev2(k)+flev2(k-1)
                  end do
                  amapt(i,j)=elev(i,j)/1000.-flev2(7) 
                  amapth(i,j)=flev2(7) 
                  ! uncomment next line if thickness without water is desired
                  if(elev(i,j) == 0.0_dp) amapth(i,j)=amapth(i,j)-flev2(1)
                  do  k=1,6
                     amapl(k,i,j)=elev(i,j)/1000.-flev2(k)
                  end do
                  cycle outer
               endif
            end do
         end do outer
      end do

      
      amapt = amapt*1000.0_dp      
      amapth = amapth*1000.0_dp
      amaps = amaps*1000.0_dp
      amapi = amapi*1000.0_dp
      amapl = amapl*1000.0_dp

890   format(////a)
891   format(//a)
899   format(a)
901   format(i4,1x,180(2x,a2,1x))
902   format(30(8e15.5/))
905   format(a,i6,1x,a2)

      close(7)

      return            
    end subroutine readCN2_t7_extract

    
    subroutine readCN2_t7

    use nrtype
    implicit none
    

    
    character(len=2), dimension(ntyp) :: ctype
    character(len=1) :: dum
    character(len=5) :: dum0
    character(len=506) :: line
    character(len=2), dimension(nlo) :: types
    
    integer(i4b) :: i,j,k,l,m,ilat
    integer(i4b), dimension(nlo) :: ilon
    integer(i4b), dimension(nla) :: il      

    real(dp) :: aux,felev,fthick,flons
    real(dp), dimension(nlo) :: amapt
    real(dp), dimension(nlo) :: amapth
    real(dp), dimension(nlo) :: amaps
    real(dp), dimension(nlo) :: amapi
    real(dp), dimension(6,nlo) :: amapl
    real(dp), dimension(ntyp,7) :: flev
    real(dp), dimension(7) :: flev2
    real(dp), dimension(nlo,nla) :: elev
 
    
    ! read in elevation map

    open(3,file='/home//da380/raid/dta/crust2/CNelevatio2.txt')
    open(19,file='map2.topo')
    read(3,*) ilon
    do j=1,nla
       read(3,*) il(j),(elev(i,j),i=1,nlo) 
       write(19,902)(elev(i,j)/1000.,i=nla+1,nlo), &
            (elev(i,j)/1000.,i=1,nla)
    end do
    close(3)
    close(19)

    open(2,file='/home/da380/raid/dta/crust2/CNtype2_key.txt')
    open(7,file='/home/da380/raid/dta/crust2/CNtype2.txt')
    open(8,file='map2.t7')
    open(16,file='map2.t0')
    open(9,file='map2.t1')
    open(10,file='map2.t2')
    open(11,file='map2.t3')
    open(12,file='map2.t4')
    open(13,file='map2.t5')
    open(14,file='map2.t6')
    open(15,file='map2.thick')
    open(17,file='map2.sed')
    open(18,file='map2.ice')


    ! read in key for crust types

    read(2,890)dum
    do  i=1,ntyp
       read(2,899)ctype(i)
       print 899,ctype(i)
       read(2,891)dum
       read(2,899)line
       read(line,*)(flev(i,l),l=1,7),dum0,fthick
       ! switch layer 1 and 2
       ! new layer 1 is water
       ! new layer 2 is ice
       aux=flev(i,2)
       flev(i,2)=flev(i,1)
       flev(i,1)=aux
       felev=0.
       do l=1,7
          felev=felev+flev(i,l)
       end do
       if(felev /= fthick) then
          print*,' read error in thickness of crust'
          print 905,' type: ',i,ctype(i)
          print*,felev,fthick 
          stop
       end if
    end do
    print*,'reading types completed'
    close(2)


    ! read CNtype file
    read(7,*)flons
    do  j=1,nla
       read(7,901) ilat,types
       print*,ilat
       outer:do   i=1,nlo
          do  l=1,ntyp
             if(types(i) == ctype(l)) then
                do  m=1,7
                   flev2(m)=flev(l,m)
                end do
                amaps(i)=flev2(3)+flev2(4)
                amapi(i)=flev2(2)
                if(elev(i,j) < 0.0_dp)then
                   flev2(1) = -elev(i,j)/1000.0_dp
                   elev(i,j) = 0.0_dp      
                end if
                do k=2,7
                   flev2(k)=flev2(k)+flev2(k-1)
                end do
                amapt(i)=elev(i,j)/1000.-flev2(7) 
                amapth(i)=flev2(7) 
                ! uncomment next line if thickness without water is desired
                if(elev(i,j) == 0.0_dp) amapth(i)=amapth(i)-flev2(1)
                do  k=1,6
                   amapl(k,i)=elev(i,j)/1000.-flev2(k)
                end do
                cycle outer
             endif
          end do
          print*,' crust type code not found: ',types(i)
          print*,' latitude: ',ilat,' long index: ',i
       end do outer
       write(8,902) (amapt(l),l=nla+1,nlo),(amapt(l),l=1,nla)
       write(15,902) (amapth(l),l=nla+1,nlo),(amapth(l),l=1,nla)
       write(16,902)(elev(l,j)/1000.,l=nla+1,nlo),(elev(l,j)/1000.,l=1,nla)
       write(17,902)(amaps(l),l=nla+1,nlo),(amaps(l),l=1,nla)
       write(18,902) (amapi(l),l=nla+1,nlo),(amapi(l),l=1,nla)
       do  k=1,6
          write(8+k,902) (amapl(k,l),l=nla+1,nlo),(amapl(k,l),l=1,nla)
       end do
    end do
    
890 format(////a)
891 format(//a)
899 format(a)
901 format(i4,1x,180(2x,a2,1x))
902 format(30(8e15.5/))
905 format(a,i6,1x,a2)


    close(7)
    close(8)
    close(16)
    close(9)
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(17)
    close(18)

    return
  end subroutine readCN2_t7
  
  
  subroutine readCN2_7vr_extract(amapvp,amaprho,amapvs)
    use nrtype
    implicit none

    real(dp), dimension(8,nlo,nla), intent(out) :: amapvp,amaprho,amapvs

    character(len=2), dimension(ntyp) :: ctype
    character(len=1) :: dum
    character(len=506) :: line
    character(len=2), dimension(nlo) :: types

    integer(i4b) :: i,j,k,l,ilat
    
    real(dp) :: aux,flons
    real(dp), dimension(ntyp,8) :: fvel,fvels,frho



    open(2,file='/home/da380/raid/dta/crust2/CNtype2_key.txt')
    open(7,file='/home/da380/raid/dta/crust2/CNtype2.txt')


    ! read in key for crust types
    
    read(2,890)dum
    do i=1,ntyp
       read(2,899)ctype(i)
       read(2,899)line
       read(line,*)(fvel(i,l),l=1,8)
       read(2,899)line
       read(line,*)(fvels(i,l),l=1,8)
       read(2,899)line
       read(line,*)(frho(i,l),l=1,8)
       read(2,899)dum
       ! flip layers
       aux=fvel(i,1)
       fvel(i,1)=fvel(i,2)
       fvel(i,2)=aux
       aux=fvels(i,1)
       fvels(i,1)=fvels(i,2)
       fvels(i,2)=aux
       aux=frho(i,1)
       frho(i,1)=frho(i,2)
       frho(i,2)=aux
    end do

    ! read CNtype file
    read(7,899)line
    read(line,*)flons
    do j=1,nla
       read(7,901)ilat,types       
       inner: do i=1,nlo
          do l=1,ntyp
             if(types(i).eq.ctype(l))then
                do  k=1,8
                   amapvp(k,i,j)=fvel(l,k)
                   amapvs(k,i,j)=fvels(l,k)
                   amaprho(k,i,j)=frho(l,k)
                end do
                cycle inner
             endif
          end do
       end do inner
    end do

    amapvp = amapvp*1000.0_dp
    amaprho = amaprho*1000.0_dp
    amapvs = amapvs*1000.0_dp

    
890 format(////a)
899 format(a)
901 format(i4,1x,180(2x,a2,1x))
902 format(30(8e15.5/))
    
    close(7)
    close(2)
    
    return
  end subroutine readCN2_7vr_extract
  

  subroutine readCN2_7vr
    use nrtype
    implicit none

    character(len=2), dimension(ntyp) :: ctype
    character(len=1) :: dum
    character(len=506) :: line
    character(len=2), dimension(nlo) :: types

    integer(i4b) :: i,j,k,l,ilat
    
    real(dp) :: aux,flons
    real(dp), dimension(ntyp,8) :: fvel,fvels,frho
    real(dp), dimension(8,nlo,nla) :: amapvp,amaprho,amapvs


    open(2,file='/home/da380/raid/dta/crust2/CNtype2_key.txt')
    open(7,file='/home/da380/raid/dta/crust2/CNtype2.txt')


    ! read in key for crust types
    
    read(2,890)dum
    do i=1,ntyp
       read(2,899)ctype(i)
       print 899,ctype(i)
       read(2,899)line
       read(line,*)(fvel(i,l),l=1,8)
       read(2,899)line
       read(line,*)(fvels(i,l),l=1,8)
       read(2,899)line
       read(line,*)(frho(i,l),l=1,8)
       read(2,899)dum
       ! flip layers
       aux=fvel(i,1)
       fvel(i,1)=fvel(i,2)
       fvel(i,2)=aux
       aux=fvels(i,1)
       fvels(i,1)=fvels(i,2)
       fvels(i,2)=aux
       aux=frho(i,1)
       frho(i,1)=frho(i,2)
       frho(i,2)=aux
    end do
    
    ! read CNtype file
    read(7,899)line
    read(line,*)flons
    do j=1,nla
       read(7,901)ilat,types
       print*,ilat
       inner: do i=1,nlo
          do l=1,ntyp
             if(types(i).eq.ctype(l))then
                do  k=1,8
                   amapvp(k,i,j)=fvel(l,k)
                   amapvs(k,i,j)=fvels(l,k)
                   amaprho(k,i,j)=frho(l,k)
                end do
                cycle inner
             endif
          end do
          print*,' crust type code not found: ',types(i)
          print*,' latitude: ',ilat,' long index: ',i
       end do inner
    end do
    do k=1,8
       write(dum,'(i1)')k
       open(8+k,file='map2.vp'//dum)
       open(18+k,file='map2.vs'//dum)
       open(28+k,file='map2.rho'//dum)
       do i=1,nla
          write(8+k,902)(amapvp(k,l,i),l=nla+1,nlo), & 
               (amapvp(k,l,i),l=1,nla)
          write(18+k,902)(amapvs(k,l,i),l=nla+1,nlo), & 
               (amapvs(k,l,i),l=1,nla)
          write(28+k,902)(amaprho(k,l,i),l=nla+1,nlo), & 
               (amaprho(k,l,i),l=1,nla)
       end do
       close(8+k)
       close(18+k)
       close(28+k)
    end do
    
890 format(////a)
899 format(a)
901 format(i4,1x,180(2x,a2,1x))
902 format(30(8e15.5/))
    
    close(7)
    close(2)
    return
  end subroutine readCN2_7vr






  
  subroutine crustal_average(lat,lon,dr,drho,dvp,dvs)
    use nrtype
    implicit none
    
    real(dp), intent(in) :: lat,lon
    
    real(dp), dimension(4), intent(out) :: dr
    real(dp), dimension(3), intent(out) :: drho
    real(dp), dimension(3), intent(out) :: dvp
    real(dp), dimension(3), intent(out) :: dvs
    
    
    real(dp), parameter :: r0 = 6371000.0_dp
    real(dp), parameter :: rsol = 6368000.0_dp
    real(dp), parameter :: rmidc = 6356000.0_dp
    real(dp), parameter :: rmoho = 6346600.0_dp
    
    real(dp), parameter :: rho_ocean  = 1020.0_dp
    real(dp), parameter :: rho_ucrust = 2600.0_dp
    real(dp), parameter :: rho_lcrust = 2900.0_dp
    real(dp), parameter :: rho_smoho  = 3380.75_dp
    
    real(dp), parameter :: vp_ocean  = 1450.0_dp
    real(dp), parameter :: vp_ucrust = 5800.0_dp
    real(dp), parameter :: vp_lcrust = 6800.0_dp
    real(dp), parameter :: vp_smoho  = 8022.06_dp
    
    real(dp), parameter :: vs_ocean  = 0.0_dp
    real(dp), parameter :: vs_ucrust = 3200.0_dp
    real(dp), parameter :: vs_lcrust = 3900.0_dp
    real(dp), parameter :: vs_smoho  = 4396.02_dp
    
    
    real(dp) :: felev,famapt,famapth,famaps,famapi,tmp1,tmp2
    real(dp), dimension(6) :: famapl
    real(dp), dimension(8) :: famapvp
    real(dp), dimension(8) :: famapvs
    real(dp), dimension(8) :: famaprho
    

    real(dp) :: thick_ocean,thick_ucrust,thick_lcrust
    
    ! get the local crustal structure
    call crust_eval(lat,lon,felev=felev,        & 
         famapt  = famapt,   &
         famapth = famapth,  &
         famaps  = famaps,   &
         famapi  = famapi,   &
         famapl = famapl,    & 
         famapvp = famapvp,   & 
         famapvs = famapvs,   & 
         famaprho = famaprho)  
    
    

    ! get the boundary perturbations
    dr(1) = felev/r0
    dr(2) = ((r0-rsol)+famapl(1))/rsol
    dr(3) = ((r0-rmidc)+famapl(5))/rmidc
    dr(4) = ((r0-rmoho)+famapt)/rmoho
    

    
    ! parameter perturbations in the upper crust
    
    tmp1 =  famaprho(2)*(famapl(1)-famapl(2)) & 
         +famaprho(3)*(famapl(2)-famapl(3)) & 
         +famaprho(4)*(famapl(3)-famapl(4)) & 
         +famaprho(5)*(famapl(4)-famapl(5))
    tmp2 = famapl(1)-famapl(5)
    drho(1) = (tmp1/tmp2-rho_ucrust)/rho_ucrust
    
    tmp1 =  famapvp(2)*(famapl(1)-famapl(2)) & 
         +famapvp(3)*(famapl(2)-famapl(3)) & 
         +famapvp(4)*(famapl(3)-famapl(4)) & 
         +famapvp(5)*(famapl(4)-famapl(5))
    tmp2 = famapl(1)-famapl(5)
    dvp(1) = (tmp1/tmp2-vp_ucrust)/vp_ucrust
    
    
    tmp1 =  famapvs(2)*(famapl(1)-famapl(2)) & 
         +famapvs(3)*(famapl(2)-famapl(3)) & 
         +famapvs(4)*(famapl(3)-famapl(4)) & 
         +famapvs(5)*(famapl(4)-famapl(5))
    tmp2 = famapl(1)-famapl(5)
    dvs(1) = (tmp1/tmp2-vs_ucrust)/vs_ucrust
    
    
    ! parameter perturbations in the lower crust
    
    tmp1 =  famaprho(6)*(famapl(5)-famapl(6)) & 
         +famaprho(7)*(famapl(6)-famapt)
    tmp2 = famapl(5)-famapt
    drho(2) = (tmp1/tmp2-rho_lcrust)/rho_lcrust
    
    tmp1 =  famapvp(6)*(famapl(5)-famapl(6)) & 
         +famapvp(7)*(famapl(6)-famapt)
    tmp2 = famapl(5)-famapt
    dvp(2) = (tmp1/tmp2-vp_lcrust)/vp_lcrust
    
    tmp1 =  famapvs(6)*(famapl(5)-famapl(6)) & 
         +famapvs(7)*(famapl(6)-famapt)
    tmp2 = famapl(5)-famapt
    dvs(2) = (tmp1/tmp2-vs_lcrust)/vs_lcrust
    
    
    ! parameter perturbation below the moho
    drho(3) = (famaprho(8)-rho_smoho)/rho_smoho
    dvp(3)  = (famapvp(8)-vp_smoho)/vp_smoho
    dvs(3)  = (famapvs(8)-vs_smoho)/vs_smoho

      

    return
  end subroutine crustal_average




  subroutine crustal_average_iso(lat,lon,dr,drho,dvp,dvs)
    use nrtype
    implicit none
    
    real(dp), intent(in) :: lat,lon
    
    real(dp), dimension(4), intent(out) :: dr
    real(dp), dimension(3), intent(out) :: drho
    real(dp), dimension(3), intent(out) :: dvp
    real(dp), dimension(3), intent(out) :: dvs
    
    
    real(dp), parameter :: r0 = 6371000.0_dp
    real(dp), parameter :: rsol = 6368000.0_dp
    real(dp), parameter :: rmidc = 6356000.0_dp
    real(dp), parameter :: rmoho = 6346600.0_dp
    
    real(dp), parameter :: rho_ocean  = 1020.0_dp
    real(dp), parameter :: rho_ucrust = 2600.0_dp
    real(dp), parameter :: rho_lcrust = 2900.0_dp
    real(dp), parameter :: rho_smoho  = 3380.75_dp
    
    real(dp), parameter :: vp_ocean  = 1450.0_dp
    real(dp), parameter :: vp_ucrust = 5800.0_dp
    real(dp), parameter :: vp_lcrust = 6800.0_dp
    real(dp), parameter :: vp_smoho  = 8022.06_dp
    
    real(dp), parameter :: vs_ocean  = 0.0_dp
    real(dp), parameter :: vs_ucrust = 3200.0_dp
    real(dp), parameter :: vs_lcrust = 3900.0_dp
    real(dp), parameter :: vs_smoho  = 4396.02_dp
    
    
    real(dp) :: felev,famapt,famapth,famaps,famapi,tmp1,tmp2
    real(dp), dimension(6) :: famapl
    real(dp), dimension(8) :: famapvp
    real(dp), dimension(8) :: famapvs
    real(dp), dimension(8) :: famaprho
    

    real(dp) :: thick_ocean,thick_ucrust,thick_lcrust
    real(dp) :: prem_mass,rcomp,loc_mass,loc_cmass,loc_wmass
    real(dp), parameter :: dcomp = 100.0_dp ! chosen compensation
                                            ! depth in km
    real(dp) :: cdepth
    

    ! compute the mass of prem crust down to the 
    ! given compensation depth
    rcomp = r0-dcomp*1000.0_dp
    prem_mass = rho_ocean*(r0-rsol)        &
                +rho_ucrust*(rsol-rmidc)   &
                +rho_lcrust*(rmidc-rmoho)  & 
                +rho_smoho*(rmoho-rcomp)

    ! get the local crustal structure
    call crust_eval(lat,lon,felev=felev,        & 
         famapt  = famapt,   &
         famapth = famapth,  &
         famaps  = famaps,   &
         famapi  = famapi,   &
         famapl = famapl,    & 
         famapvp = famapvp,   & 
         famapvs = famapvs,   & 
         famaprho = famaprho)  


    
    

    ! get the boundary perturbations
    dr(1) = felev/r0
    dr(2) = ((r0-rsol)+famapl(1))/rsol
    dr(3) = ((r0-rmidc)+famapl(5))/rmidc
    dr(4) = ((r0-rmoho)+famapt)/rmoho
    
    cdepth = -1000.0_dp*dcomp

    call  crustal_balance(felev,famapl,famapt,famaprho,cdepth)    


    
    ! parameter perturbations in the upper crust
    
    tmp1 =  famaprho(2)*(famapl(1)-famapl(2)) & 
         +famaprho(3)*(famapl(2)-famapl(3)) & 
         +famaprho(4)*(famapl(3)-famapl(4)) & 
         +famaprho(5)*(famapl(4)-famapl(5))
    tmp2 = famapl(1)-famapl(5)
    drho(1) = (tmp1/tmp2-rho_ucrust)/rho_ucrust
    
    tmp1 =  famapvp(2)*(famapl(1)-famapl(2)) & 
         +famapvp(3)*(famapl(2)-famapl(3)) & 
         +famapvp(4)*(famapl(3)-famapl(4)) & 
         +famapvp(5)*(famapl(4)-famapl(5))
    tmp2 = famapl(1)-famapl(5)
    dvp(1) = (tmp1/tmp2-vp_ucrust)/vp_ucrust
    
    
    tmp1 =  famapvs(2)*(famapl(1)-famapl(2)) & 
         +famapvs(3)*(famapl(2)-famapl(3)) & 
         +famapvs(4)*(famapl(3)-famapl(4)) & 
         +famapvs(5)*(famapl(4)-famapl(5))
    tmp2 = famapl(1)-famapl(5)
    dvs(1) = (tmp1/tmp2-vs_ucrust)/vs_ucrust
    
    
    ! parameter perturbations in the lower crust
    
    tmp1 =  famaprho(6)*(famapl(5)-famapl(6)) & 
         +famaprho(7)*(famapl(6)-famapt)

    tmp2 = famapl(5)-famapt
    drho(2) = (tmp1/tmp2-rho_lcrust)/rho_lcrust
    
    tmp1 =  famapvp(6)*(famapl(5)-famapl(6)) & 
         +famapvp(7)*(famapl(6)-famapt)
    tmp2 = famapl(5)-famapt
    dvp(2) = (tmp1/tmp2-vp_lcrust)/vp_lcrust
    
    tmp1 =  famapvs(6)*(famapl(5)-famapl(6)) & 
         +famapvs(7)*(famapl(6)-famapt)
    tmp2 = famapl(5)-famapt
    dvs(2) = (tmp1/tmp2-vs_lcrust)/vs_lcrust
    
    
    ! parameter perturbation below the moho
    drho(3) = (famaprho(8)-rho_smoho)/rho_smoho
    dvp(3)  = (famapvp(8)-vp_smoho)/vp_smoho
    dvs(3)  = (famapvs(8)-vs_smoho)/vs_smoho






    return
  end subroutine crustal_average_iso
  

    

  subroutine crustal_balance(felev,famapl,famapt,famaprho,cdepth)
    use nrtype
    implicit none
    real(dp), intent(in) :: felev,famapt,cdepth
    real(dp), dimension(6), intent(in) :: famapl
    real(dp), dimension(8), intent(inout) :: famaprho
    integer(i4b) :: i
    
    
    real(dp), parameter :: r0 = 6371000.0_dp
    real(dp), parameter :: rsol = 6368000.0_dp
    real(dp), parameter :: rmidc = 6356000.0_dp
    real(dp), parameter :: rmoho = 6346600.0_dp
    
    real(dp), parameter :: rho_ocean  = 1020.0_dp
    real(dp), parameter :: rho_ucrust = 2600.0_dp
    real(dp), parameter :: rho_lcrust = 2900.0_dp
    real(dp), parameter :: rho_smoho  = 3380.75_dp
    

    real(dp) :: hw,hi,hss,hhs,huc,hmc,hlc,hcd
    real(dp) :: hwp,hucp,hlcp,hcdp,rcd
    real(dp) :: rhow,rhoi,rhoss,rhohs,rhouc,rhomc,rholc
    real(dp) :: mc,m1,m2,m3,m4,mp
    real(dp) :: fac
    
    ! work out the thicknesses of the different layers
    hw = felev-famapl(1)
    hi = famapl(1)-famapl(2)
    hss = famapl(2)-famapl(3)
    hhs = famapl(3)-famapl(4)
    huc = famapl(4)-famapl(5)
    hmc = famapl(5)-famapl(6)
    hlc = famapl(6)-famapt
    hcd = famapt-cdepth
      
    ! get the densities for the different layers
    rhow   = famaprho(1)
    rhoi   = famaprho(2)
    rhoss  = famaprho(3)
    rhohs  = famaprho(4)
    rhouc  = famaprho(5)
    rhomc  = famaprho(6)
    rholc  = famaprho(7)

    ! work out the mass of the various layers of the crust
    m1 = rhow*hw+rhoi*hi
    m2 = rhoss*hss+rhohs*hhs+rhouc*huc 
    m3 = rhomc*hmc+rholc*hlc
    m4 = rho_smoho*hcd
    
    ! work out the mass of the prem crust
    rcd = r0+cdepth
    hwp = r0-rsol
    hucp = rsol-rmidc
    hlcp = rmidc-rmoho
    hcdp = rmoho-rcd
    mp = rho_ocean*hwp+rho_ucrust*hucp & 
         +rho_lcrust*hlcp + rho_smoho*hcdp
    
      
    ! calculate the density scaling factor
    fac = (mp-m1-m4)/(m2+m3)


    ! scale the density in the crust
    famaprho(3:7) = fac*famaprho(3:7)

    return
  end subroutine crustal_balance



  subroutine crustal_mass(felev,famapl,famapt,famaprho,cdepth,cmass)
    use nrtype
    implicit none
    real(dp), intent(in) :: felev,famapt,cdepth
    real(dp), dimension(6), intent(in) :: famapl
    real(dp), dimension(8), intent(in) :: famaprho
    real(dp), intent(out) :: cmass
    integer(i4b) :: i

    
    real(dp), parameter :: r0 = 6371000.0_dp
    real(dp), parameter :: rsol = 6368000.0_dp
    real(dp), parameter :: rmidc = 6356000.0_dp
    real(dp), parameter :: rmoho = 6346600.0_dp
      
    real(dp), parameter :: rho_ocean  = 1020.0_dp
    real(dp), parameter :: rho_ucrust = 2600.0_dp
    real(dp), parameter :: rho_lcrust = 2900.0_dp
    real(dp), parameter :: rho_smoho  = 3380.75_dp


    real(dp) :: hw,hi,hss,hhs,huc,hmc,hlc,hcd
    real(dp) :: hwp,hucp,hlcp,hcdp,rcd
    real(dp) :: rhow,rhoi,rhoss,rhohs,rhouc,rhomc,rholc
    real(dp) :: mc,m1,m2,m3,m4,mp
    real(dp) :: fac

    ! work out the thicknesses of the different layers
    hw = felev-famapl(1)
    hi = famapl(1)-famapl(2)
    hss = famapl(2)-famapl(3)
    hhs = famapl(3)-famapl(4)
    huc = famapl(4)-famapl(5)
    hmc = famapl(5)-famapl(6)
    hlc = famapl(6)-famapt
    hcd = famapt-cdepth

    ! get the densities for the different layers
    rhow   = famaprho(1)
    rhoi   = famaprho(2)
    rhoss  = famaprho(3)
    rhohs  = famaprho(4)
    rhouc  = famaprho(5)
    rhomc  = famaprho(6)
    rholc  = famaprho(7)

    ! work out the mass of the various layers of the crust
    m1 = rhow*hw+rhoi*hi
    m2 = rhoss*hss+rhohs*hhs+rhouc*huc 
    m3 = rhomc*hmc+rholc*hlc
    m4 = rho_smoho*hcd
    
    ! sum the masses
    cmass = m1+m2+m3+m4


    return
  end subroutine crustal_mass
  
  



  subroutine get_topo(lat,lon,topo)
    use nrtype
    implicit none
    
    real(dp), intent(in) :: lat,lon    
    real(dp), intent(out) :: topo  
    real(dp), dimension(6) :: famapl

    ! get the local crustal structure
    call crust_eval(lat,lon,famapl=famapl)        
    topo = famapl(1)

    return
  end subroutine get_topo



  subroutine ocean_function(lat,lon,ofun)
    use nrtype
    implicit none
    
    real(dp), intent(in) :: lat,lon    
    real(dp), intent(out) :: ofun
      

    real(dp), dimension(6) :: famapl
    real(dp) :: topo

    ! get the local crustal structure
    call crust_eval(lat,lon,famapl=famapl)        
    topo = famapl(1)
    if(topo > 0.0_dp) then
       ofun = 0.0_dp
    else
      ofun = 1.0_dp
    end if

    return
  end subroutine ocean_function


  function topo_filter(topo,t11,t12)
    ! imputs in m
    use nrtype
    implicit none
    real(dp) :: topo_filter
    real(dp), intent(in) :: topo,t11,t12
    real(dp) :: tmp

    if(topo < t11) then
       topo_filter = 0.0_dp
    elseif(topo >= t11 .and. topo < t12) then
       tmp = pi_d*(topo-t11)/(t12-t11)
       topo_filter = 0.5_dp*(1.0_dp-cos(tmp))
    elseif(topo >= t12) then
       topo_filter = 1.0_dp
    end if
    return    
  end function topo_filter





end module module_crust
