module module_function


  contains



    subroutine ylm_real(lmax,theta,phi,ylm)
      use nrtype
      implicit none
      ! given a position (theta,phi) on the unit sphere, and a maximum
      ! spherical harmonic degree lmax, this routine returns a vector
      ! of the real spherical harmonics at this point ordered as follows:
      !
      ! ylm(1) = Y_{0,0}
      ! ylm(2) = Y_{1,-1}
      ! ylm(3) = Y_{1,0}
      ! ylm(4) = Y_{1,1}
      ! ylm(5) = Y_{2,-2} ...
      ! 
      ! note these functions are the fully normalizes real spherical 
      ! harmonics as defined in Appendix B6 of Dahlen & Tromp 1998.
      
      integer(i4b), intent(in) :: lmax
      real(dp), intent(in) :: theta
      real(dp), intent(in) :: phi
      real(dp), dimension(:), intent(out) :: ylm

      real(dp), dimension(lmax+1) :: x,xp,xc


      integer(i4b) :: l,m,ilm,im
      real(dp) :: cmp,smp

      ilm = 0
      do l = 0,lmax

         ! calculate the Legendre functions for this l
         call legendre(theta,l,l,x(1:l+1),xp(1:l+1),xc(1:l+1))
      

         do m = -l,l
            ilm = ilm+1

            if(m < 0)  then
               cmp = cos(m*phi)
               ylm(ilm) = sqrt2*x(1-m)*cmp
            else if(m == 0) then
               ylm(ilm) = x(1)
            else if(m > 0) then
               smp = sin(m*phi)
               ylm(ilm) = sqrt2*x(1+m)*smp
            end if

         end do

      end do

      return
    end subroutine ylm_real

    subroutine glegendre(theta,l,x0,x1m,x1p,x2m,x2p)
      ! this routine returns the generalized Legendre polynomials
      ! X^{N}_{lm} for a given l in the range  0 <= m <= l and 
      ! -2 <= N <= 2 
      ! the outputs are given as linear poloidal and toroidal 
      ! combinations of the functions as follows:
      !
      !  x0 = X^{0}_{lm}
      !
      !  x1m = 0.5*zeta*(X^{-1}_{lm}-X^{+1}_{lm})
      !  x1p = 0.5*zeta*(X^{-1}_{lm}+X^{+1}_{lm})
      !
      !  x2m = 0.5*zeta*sqrt(zeta^2-2)*(X^{-2}_{lm}-X^{+2}_{lm})
      !  x2p = 0.5*zeta*sqrt(zeta^2-2)*(X^{-2}_{lm}+X^{+2}_{lm})
      !
      ! where zeta^2 = l*(l+1)
      !
      ! note that the returned functions are fully normalized
      ! such that 
      !
      ! Y^{N}_{lm}(\theta,\phi) = X^{N}_{lm}(\theta) \ee^{\ii m \phi}
      ! 
      ! is the fully normalized generalized spherical harmonic
      ! as defined in Dalhen \& Tromp (1998) 
      !
      ! the routine is based on a called to John Woodhouse's
      ! code rotmx2 

      use nrtype
      implicit none
      real(dp), intent(in) :: theta
      integer(i4b), intent(in) :: l
      real(dp), dimension(:), intent(out) :: x0
      real(dp), dimension(:), intent(out) :: x1m
      real(dp), dimension(:), intent(out) :: x1p
      real(dp), dimension(:), intent(out) :: x2m
      real(dp), dimension(:), intent(out) :: x2p

      integer(i4b) :: m,im
      real(dp) :: zeta,szeta2m2,nu,fac1,fac2
      real(dp), dimension(5,2*l+1) :: dnm

      ! calculate the required d-functions
      dnm = 0.0_dp
      if(l == 0) then
         call ROTMX2(0,l,theta,dnm(3:3,:),1,2*l+1)
      else if(l == 1) then
         call ROTMX2(1,l,theta,dnm(2:4,:),3,2*l+1)
      else if(l >= 2) then
         call ROTMX2(2,l,theta,dnm,5,2*l+1)
      end if

      
      ! calculate some parameters
      zeta = l*(l+1)
      if(l /= 0.0_dp) then
         szeta2m2 = sqrt(zeta-2.0_dp)
      else
         szeta2m2 = 0.0_dp
      end if
      zeta = sqrt(zeta)
      nu = sqrt((2*l+1)/fourpi_d)
      fac1 = 0.5_dp*nu*zeta
      fac2 = 0.5_dp*nu*zeta*szeta2m2



      do m = 0,l
         im = m+l+1
         x0(m+1) = nu*dnm(3,im)
         x1m(m+1) = fac1*(dnm(2,im)-dnm(4,im))
         x1p(m+1) = fac1*(dnm(2,im)+dnm(4,im))
         x2m(m+1) = fac2*(dnm(1,im)-dnm(5,im))
         x2p(m+1) = fac2*(dnm(1,im)+dnm(5,im))
      end do


      return
    end subroutine glegendre


    subroutine legendre(theta,l,m,x,xp,xc)
      
      ! This routine computes the associated legendre 
      ! polynomials P_{lm}(cos(theta)) along with
      ! their angular derivative and  products with
      ! cosec(theta). This routine is a f90 translation
      ! of the routine legndr by John Woodhouse      
      ! the output arrays should have dimension at least as 
      ! big as m+1
      
      ! inputs: 
      !          theta = angular argument in radians
      !          l     = angular degree of the Legendre functions
      !          m     = maximum angular order to calculate
      !
      ! outputs:
      !          x  = array of the Legendre functions; the array
      !               order is such that x(1) = P_{l0}, x(2) = P_{l1},
      !               etc upto x(m+1) = P_{lm}. 
      !          xp = theta derivatives of the Legendre functions
      !          xc = legendre functions times by cosec(theta).
      
      ! the legendre functions returned are such that the fully normalized 
      ! spherical harmonic functions Y_{lm} are given by
      ! 
      ! Y_{lm}(theta,phi) = P_{lm}(cos(theta))*exp(i*m*phi),
      !
      ! The P_{lm}(cos(theta)) calculated are equal to the
      ! X_{lm}(theta) described in appendix B of Dahlen & Tromp (1998).
      
      
      use nrtype
      implicit none
      
      ! inputs:
      real(dp), intent(in) :: theta
      integer(i4b), intent(in) :: l,m
      ! outputs:
      real(dp), dimension(:), intent(out) :: x,xp,xc
      
      ! local variables
      integer(i4b) :: lp1,mp1,i,k
      real(dp) :: sum,th,ct,st,fct,sfl3,compar,dsfl3, &
           cot,cosec,x1,x2,x3,f1,f2,xm,small,lsign
      
      sum = 0.0_dp
      lp1 = l+1
      mp1 = m+1
      
      th = theta
      ct = cos(th)
      st = sin(th)
      
      fct    = sqrt(real(2*l+1)/fourpi_d)
      sfl3   = sqrt(real(l*(l+1)))
      compar = real(2*l+1)/fourpi_d
      dsfl3  = sfl3
      small  = 1.0e-16_dp*compar
      
      x      = 0.0_dp
      xp     = 0.0_dp
      xc = 0.0_dp
      
      if(l <= 1 .or. abs(theta) <= 1.0e-5_dp) then
         x(1) = fct
         if(l == 0) return     
         x(1)  = ct*fct
         x(2)  = -0.5_dp*st*fct*dsfl3
         xp(1) = -0.5_dp*st*fct*dsfl3*dsfl3
         xp(2)  = -0.5_dp*ct*fct*dsfl3
         if(abs(theta) <  1.0e-5_dp) xc(2) = xp(2)
         if(abs(theta) >= 1.0e-5_dp) xc(2) = x(2)/st
         return
      end if

      if(abs(pi_d-theta) <= 1.0e-4_dp) then
         
         lsign = (-1)**l
         x(1) = -lsign*fct*ct
         x(2) = lsign*0.5_dp*fct*dsfl3*st
         xp(1) = lsign*0.5_dp*fct*dsfl3**2*st
         xp(2) = lsign*0.5_dp*fct*dsfl3*ct
         xc(2) = -xp(2)

         return
      end if
      
      x1 = 1.0_dp
      x2 = ct
      
      do i = 2,l
         x3 = (real(2*i-1)*ct*x2-real(i-1)*x1)/real(i)
         x1 = x2
         x2 = x3
      end do
      
      cot   = ct/st;
      cosec = 1.0_dp/st
      
      x3 = x2*fct
      x2 = real(l)*(x1-ct*x2)*fct/st

      x(1) = x3
      x(2) = x2
      sum  = x3*x3
      
      xp(1) = -x2
      xp(2) = real(l*(l+1))*x3-cot*x2
      
      x(2)      = -x(2)/sfl3
      xc(2) = x(2)*cosec
      xp(2)     = -xp(2)/sfl3
      
      sum = sum+2.0_dp*x(2)*x(2)
      if(sum-compar > small) return

      x1 =  x3
      x2 = -x2/sqrt(real(l*(l+1)))
      do i = 3,mp1
         k   = i-1
         f1  = sqrt(real(l+i-1)*(l-i+2))
         f2  = sqrt(real(l+i-2)*(l-i+3))
         xm  = k
         x3  = -(2.0_dp*cot*(xm-1.0_dp)*x2+f2*x1)/f1
         sum = sum+2.0_dp*x3*x3
         if(sum-compar > small .and. i /= lp1) return
         x(i)      = x3
         xc(i) = x(i)*cosec
         x1        = x2
         xp(i)     = -(f1*x2+xm*cot*x3)
         x2        = x3
      end do
      
      return
    end subroutine legendre



  subroutine sphsum2(th,ph,lmax,flm,f)
    use nrtype
    implicit none
    real(dp), intent(in) :: th,ph
    integer(i4b), intent(in) :: lmax
    complex(dpc), dimension(0:lmax,0:lmax) :: flm
    real(dp), intent(out) :: f
    
    integer(i4b) :: l,m
    real(dp) :: fac
    real(dp), dimension(lmax+1) :: x,xp,xcosec

    
    f = 0.0_dp
    do l = 0,lmax
       call legendre(th,l,l,x,xp,xcosec)
       do m = 0,l
          if(m == 0) then
             fac = 1.0_dp
          else
             fac = 2.0_dp
          end if
          f = f +fac*real(flm(l,m)*x(m+1)*exp(ii*m*ph))
       end do
       
    end do     

    return
  end subroutine sphsum2


  function erf(x)

    ! # MS Fortran
    ! Error function from Numerical Recipes.
    ! erf(x) = 1 - erfc(x)
    use nrtype
    
    implicit none
    
    real(dp) :: erf
    real(dp), intent(in) :: x
    real(dp) :: dumerfc
    real(dp)::  t, z


    z = abs(x)
    t = 1.0 / ( 1.0 + 0.5 * z )

    dumerfc = t * exp(-z * z - 1.26551223 + t *	     &
              ( 1.00002368 + t * ( 0.37409196 + t *  &
              ( 0.09678418 + t * (-0.18628806 + t *  &
              ( 0.27886807 + t * (-1.13520398 + t *  &
              ( 1.48851587 + t * (-0.82215223 + t *  &
              0.17087277 )))))))))

    if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc
    
    erf = 1.0 - dumerfc
    
  end function erf
   


  
end module module_function
