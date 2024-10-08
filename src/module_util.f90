module module_util

  interface my_reallocate
     module procedure my_reallocate_iv,my_reallocate_ia,my_reallocate_rv, &
          my_reallocate_dv,my_reallocate_cv,my_reallocate_zv
  end interface

  interface poly_interp
     module procedure poly_interp_r, poly_interp_z
  end interface


  interface bilinear
     module procedure bilinear_d
  end interface

  
  interface get_float
     module procedure get_float_sp, get_float_dp
  end interface


  interface get_floats
     module procedure get_floats_sp, get_floats_dp
  end interface


  interface my_dot
     module procedure my_dot_c,my_dot_z
  end interface



contains
    

    subroutine my_reallocate_iv(a,n)
      use nrtype
      implicit none
      integer(i4b), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      integer(i4b), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_iv


    subroutine my_reallocate_ia(a,n,m)
      use nrtype
      implicit none
      integer(i4b), dimension(:,:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n,m
      integer(i4b), dimension(:,:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store,m_tmp,m_store
      n_tmp=size(a,1); m_tmp=size(a,2)
      n_store=min(n,n_tmp); m_store=min(m,m_tmp)      
      allocate(a_tmp(n_store,m_store))
      a_tmp(1:n_store,1:m_store)=a(1:n_store,1:m_store)
      deallocate(a); allocate(a(n,m))
      a(1:n_store,1:m_store)=a_tmp(1:n_store,1:m_store)
      deallocate(a_tmp)
      return
    end subroutine my_reallocate_ia



    subroutine my_reallocate_rv(a,n)
      use nrtype
      implicit none
      real(sp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      real(sp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_rv


    subroutine my_reallocate_dv(a,n)
      use nrtype
      implicit none
      real(dp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      real(dp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_dv
      

    subroutine my_reallocate_cv(a,n)
      use nrtype
      implicit none
      complex(sp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      complex(sp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_cv


    subroutine my_reallocate_zv(a,n)
      use nrtype
      implicit none
      complex(dpc), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      complex(dpc), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_zv


    subroutine poly_interp_r(x1,f1,fp1,x2,f2,fp2,x,f,fp)
      use nrtype
      implicit none
      real(dp), intent(in) :: x,x1,x2
      real(dp), intent(in) :: f1,fp1,f2,fp2
      real(dp), intent(out) :: f,fp
      real(dp) :: c2,c3,delta,deltai,deltai2

      if(x == x1) then
        f=f1
        fp=fp1
        return
      else if(x == x2) then
         f=f2
         fp=fp2
         return
      end if

      delta=x2-x1

      if(delta == 0) stop 'poly_interp: delta = 0'
      
      deltai=1.0_dp/delta
      deltai2=deltai*deltai

      c2=3.0_dp*deltai2*(f2-f1)-deltai*(fp2+2.0_dp*fp1)
      c3=deltai2*(fp2+fp1)-2.0_dp*deltai*deltai2*(f2-f1)

      delta=x-x1
      deltai=delta*delta
      deltai2=deltai*delta

      f=f1+fp1*delta+c2*deltai+c3*deltai2
      fp=fp1+2.0_dp*c2*delta+3.0_dp*c3*deltai


      return
    end subroutine poly_interp_r




    subroutine poly_interp_z(x1,f1,fp1,x2,f2,fp2,x,f,fp)
      use nrtype
      implicit none
      real(dp), intent(in) :: x,x1,x2
      complex(dpc), intent(in) :: f1,fp1,f2,fp2
      complex(dpc), intent(out) :: f,fp
      real(dp) :: delta,deltai,deltai2
      complex(dpc) :: c2,c3

      if(x == x1) then
        f=f1
        fp=fp1
        return
      else if(x == x2) then
         f=f2
         fp=fp2
         return
      end if

      delta=x2-x1

      if(delta == 0) stop 'poly_interp: delta = 0'
      
      deltai=1.0_dp/delta
      deltai2=deltai*deltai

      c2=3.0_dp*deltai2*(f2-f1)-deltai*(fp2+2.0_dp*fp1)
      c3=deltai2*(fp2+fp1)-2.0_dp*deltai*deltai2*(f2-f1)

      delta=x-x1
      deltai=delta*delta
      deltai2=deltai*delta

      f=f1+fp1*delta+c2*deltai+c3*deltai2
      fp=fp1+2.0_dp*c2*delta+3.0_dp*c3*deltai


      return
    end subroutine poly_interp_z


    function bilinear_d(x1,x2,y1,y2,f11,f12,f21,f22,x,y)
      use nrtype
      implicit none
      real(dp) :: bilinear_d
      real(dp), intent(in) :: x1,x2,y1,y2,f11,f12,f21,f22,x,y
      

      bilinear_d = (f11*(x2-x)*(y2-y) & 
                   +f12*(x2-x)*(y-y1) & 
                   +f21*(x-x1)*(y2-y) &
                   +f22*(x-x1)*(y-y1))/((x2-x1)*(y2-y1))         

      return
    end function bilinear_d


    subroutine get_string(tag,string)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      character(len=*), intent(out) :: string

      write(6,'("'//tag//'")',advance='no')
      read(5,*) string
      string=adjustl(string)
      return
    end subroutine get_string

    
    subroutine get_integer(tag,int)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      integer(i4b), intent(out) :: int
      character(len=256) :: string,form,sls
      integer(i4b) :: sl,ios
      
      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         write(sls,*) sl
         sls=adjustl(sls)
         sl=len_trim(sls)
         form='(i'//sls(1:sl)//')'
         read(string,form,iostat=ios) int
         if(ios == 0) exit
         print *, ' input must be an integer'
      end do
            
      return
    end subroutine get_integer


    subroutine get_float_dp(tag,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      real(dp), intent(out) :: fpn
      character(len=256) :: string,form,sls,ips
      integer(i4b) :: sl,ios,ip

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         ip=index(string,'.')
         if(ip == 0) ip=sl
         write(sls,*) sl
         write(ips,*) sl-ip
         sls=adjustl(sls)
         sl=len_trim(sls)
         ips=adjustl(ips)
         ip=len_trim(ips)
         form='(f'//sls(1:sl)//'.'//ips(1:ip)//')'
         read(string,form,iostat=ios) fpn
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_float_dp

    subroutine get_float_sp(tag,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      real(sp), intent(out) :: fpn
      character(len=256) :: string,form,sls,ips
      integer(i4b) :: sl,ios,ip

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         ip=index(string,'.')
         if(ip == 0) ip=sl
         write(sls,*) sl
         write(ips,*) sl-ip
         sls=adjustl(sls)
         sl=len_trim(sls)
         ips=adjustl(ips)
         ip=len_trim(ips)
         form='(f'//sls(1:sl)//'.'//ips(1:ip)//')'
         read(string,form,iostat=ios) fpn
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_float_sp



    subroutine get_floats_dp(tag,nd,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      integer(i4b), intent(in) :: nd
      real(dp), dimension(nd), intent(out) :: fpn      
      character(len=256) :: form,sls,ips
      character(len=256), dimension(nd) :: strings
      integer(i4b) :: sl,ios,ip,id

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) strings(1:nd)
         do id = 1,nd
            read(strings(id),*,iostat=ios) fpn(id)
         end do
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_floats_dp


    subroutine get_floats_sp(tag,nd,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      integer(i4b), intent(in) :: nd
      real(sp), dimension(nd), intent(out) :: fpn
      character(len=256) :: string,form,sls,ips
      integer(i4b) :: sl,ios,ip

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         read(string,*,iostat=ios) fpn(:)
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_floats_sp



    function string_length(string,sign)
      use nrtype
      implicit none
      integer(i4b) :: string_length
      integer(i4b), intent(in) :: sign
      character (len=*), intent(in) :: string
      select case(sign)
      case(1)
         string_length=0
         do 
            if(string(string_length+sign:string_length+sign) == ' ') exit
            string_length=string_length+sign
         end do
      case(-1)
         string_length=len(string)
         do 
            if(string(string_length:string_length) == ' ') exit
            string_length=string_length+sign
         end do
     case default
        stop 'bad input of sign to string_length'
     end select
     return
   end function string_length


    subroutine bound(xx,x,i1,i2,err)
      ! Given a monotonically increasing array, xx, and a value
      ! x, returns the indices i1 and i2 such that
      ! xx(i1) <= x <= xx(i2). The routine checks that x is in 
      ! range, and if not returns err=.true.
      use nrtype
      implicit none
      real(dp), dimension(:), intent(in) :: xx
      real(dp), intent(in) :: x
      integer(i4b), intent(out) :: i1,i2
      logical(lgt), intent(out) :: err
      integer(i4b) :: i,n
      real(dp), parameter :: tol=1.0e-5_dp
      real(dp) :: rtol
      
      err=.false.
      n=size(xx)
      rtol=tol*xx(n)
      if(x > xx(n) .and. x < xx(1)) then
         err=.true.
         return
      end if

      if(abs(x-xx(n)) < rtol) then
         i2=n
         i1=n-1
         return
      end if
      if(abs(x-xx(1)) < rtol) then
         i2=2
         i1=1
         return
      end if
      do i=n-1,1,-1
         if(x <= xx(i+1) .and. x >= xx(i)) then
            i2=i+1
            i1=i
            return
         end if
      end do
      err=.true.
      return
    end subroutine bound



    function my_dot_c(x,y)
      use nrtype
      implicit none
      complex(spc) :: my_dot_c
     complex(spc), dimension(:), intent(in) :: x,y
      integer(i4b) :: i,n
      
      my_dot_c = 0.0_sp
     n = size(x,1)
      do i = 1,n
         my_dot_c = my_dot_c + x(i)*y(i)
      end do
      
      return
      
    end function my_dot_c


    function my_dot_z(x,y)
      use nrtype
      implicit none
      complex(dpc) :: my_dot_z
      complex(dpc), dimension(:), intent(in) :: x,y
      integer(i4b) :: i,n
      
      my_dot_z = 0.0_dp
      n = size(x,1)
      do i = 1,n
         my_dot_z = my_dot_z + x(i)*y(i)
      end do
     return
      
    end function my_dot_z




    subroutine string_cat_int(stri,i,stro)
      use nrtype
      implicit none
      
      character(len=*), intent(in) :: stri
      integer(i4b), intent(in) :: i
      character(len=*), intent(out) :: stro

      character(len=10) :: stmp
      integer(i4b) :: j

      write(stmp,'(i10)') i
      if(i /= 0) then
         j = floor(log10(real(i)))
      else
         j = 1
      end if
      stro = trim(stri)//stmp(10-j:10)
      
      return
    end subroutine string_cat_int



    !========================================================!
    !                       sorting routines                 !
    !========================================================!


    subroutine reduce_array(nrec,xa,indx,nxd,xd,rnkx,dxm_in)
      use nrtype
      implicit none
      integer(i4b), intent(in) :: nrec
      real(dp), dimension(nrec), intent(in) :: xa
      integer(i4b), dimension(nrec), intent(out) :: indx
      integer(i4b), intent(out) :: nxd
      real(dp), dimension(:), allocatable, intent(inout) :: xd
      integer(i4b), dimension(nrec), intent(out) :: rnkx
      real(dp), intent(in), optional :: dxm_in

      integer(i4b) :: irec
      real(dp) :: dxm

      if(present(dxm_in)) then
         dxm = dxm_in
      else
         dxm = 0.0_dp
      end if

      ! index the array
      call indexx_sp(real(xa),indx)

      ! work out how many different points are present
      nxd = 1
      do irec = 2,nrec
         if(abs(xa(indx(irec)) - xa(indx(irec-1))) > dxm) then
            nxd = nxd+1            
         end if
      end do
      ! store the different reciever depths
      allocate(xd(nxd))
      nxd = 1
      xd(1) = xa(indx(1))
      rnkx(indx(1)) = 1
      do irec = 2,nrec
         if(abs(xa(indx(irec)) - xa(indx(irec-1))) > dxm) then
            nxd = nxd+1            
            xd(nxd) = xa(indx(irec))
         end if
         rnkx(indx(irec)) = nxd
      end do


      return
    end subroutine reduce_array

    SUBROUTINE indexx_sp(arr,index)
      USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: arr
      INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
      INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
      REAL(SP) :: a
      INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
      INTEGER(I4B), DIMENSION(NSTACK) :: istack
      n=assert_eq(size(index),size(arr),'indexx_sp')
      index=arth(1,1,n)
      jstack=0
      l=1
      r=n
      do
         if (r-l < NN) then
            do j=l+1,r
               indext=index(j)
               a=arr(indext)
               do i=j-1,l,-1
                  if (arr(index(i)) <= a) exit
                  index(i+1)=index(i)
               end do
               index(i+1)=indext
            end do
            if (jstack == 0) RETURN
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         else
            k=(l+r)/2
            call swap(index(k),index(l+1))
            call icomp_xchg(index(l),index(r))
            call icomp_xchg(index(l+1),index(r))
            call icomp_xchg(index(l),index(l+1))
            i=l+1
            j=r
            indext=index(l+1)
            a=arr(indext)
            do
               do
                  i=i+1
                  if (arr(index(i)) >= a) exit
               end do
               do
                  j=j-1
                  if (arr(index(j)) <= a) exit
               end do
               if (j < i) exit
               call swap(index(i),index(j))
            end do
            index(l+1)=index(j)
            index(j)=indext
            jstack=jstack+2
            if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
            if (r-i+1 >= j-l) then
               istack(jstack)=r
               istack(jstack-1)=i
               r=j-1
            else
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            end if
         end if
      end do
    CONTAINS
      !BL
      SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
         INTEGER(I4B) :: swp
         if (arr(j) < arr(i)) then
            swp=i
            i=j
            j=swp
         end if
       END SUBROUTINE icomp_xchg
     END SUBROUTINE indexx_sp
     
     SUBROUTINE indexx_i4b(iarr,index)
       USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
       IMPLICIT NONE
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
       INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
       INTEGER(I4B) :: a
       INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
       INTEGER(I4B), DIMENSION(NSTACK) :: istack
       n=assert_eq(size(index),size(iarr),'indexx_sp')
       index=arth(1,1,n)
       jstack=0
       l=1
       r=n
       do
          if (r-l < NN) then
             do j=l+1,r
                indext=index(j)
                a=iarr(indext)
                do i=j-1,l,-1
                   if (iarr(index(i)) <= a) exit
                   index(i+1)=index(i)
                end do
                index(i+1)=indext
             end do
             if (jstack == 0) RETURN
             r=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+r)/2
             call swap(index(k),index(l+1))
             call icomp_xchg(index(l),index(r))
             call icomp_xchg(index(l+1),index(r))
             call icomp_xchg(index(l),index(l+1))
             i=l+1
             j=r
             indext=index(l+1)
             a=iarr(indext)
             do
                do
                   i=i+1
                   if (iarr(index(i)) >= a) exit
                end do
                do
                   j=j-1
                   if (iarr(index(j)) <= a) exit
                end do
                if (j < i) exit
                call swap(index(i),index(j))
             end do
             index(l+1)=index(j)
             index(j)=indext
             jstack=jstack+2
             if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
             if (r-i+1 >= j-l) then
                istack(jstack)=r
                istack(jstack-1)=i
                r=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do
     CONTAINS
!BL
       SUBROUTINE icomp_xchg(i,j)
         INTEGER(I4B), INTENT(INOUT) :: i,j
         INTEGER(I4B) :: swp
         if (iarr(j) < iarr(i)) then
            swp=i
            i=j
            j=swp
         end if
       END SUBROUTINE icomp_xchg
     END SUBROUTINE indexx_i4b


     SUBROUTINE sort(arr)
       USE nrtype; USE nrutil, ONLY : swap
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
       INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
       REAL(DP) :: a
       INTEGER(I4B) :: n,k,i,j,jstack,l,r
       INTEGER(I4B), DIMENSION(NSTACK) :: istack
       n=size(arr)
       jstack=0
       l=1
       r=n
       do
          if (r-l < NN) then
             do j=l+1,r
                a=arr(j)
                do i=j-1,l,-1
                   if (arr(i) <= a) exit
                   arr(i+1)=arr(i)
                end do
                arr(i+1)=a
             end do
             if (jstack == 0) RETURN
             r=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+r)/2
             call swap(arr(k),arr(l+1))
             call swap(arr(l),arr(r),arr(l)>arr(r))
             call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
             call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
             i=l+1
             j=r
             a=arr(l+1)
             do
                do
                   i=i+1
                   if (arr(i) >= a) exit
                end do
                do
                   j=j-1
                   if (arr(j) <= a) exit
                end do
                if (j < i) exit
                call swap(arr(i),arr(j))
             end do
             arr(l+1)=arr(j)
             arr(j)=a
             jstack=jstack+2
             if (jstack > NSTACK) stop 'sort: NSTACK too small'
             if (r-i+1 >= j-l) then
                istack(jstack)=r
                istack(jstack-1)=i
                r=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             end if
          end if
       end do
     END SUBROUTINE sort


     
   end module module_util
