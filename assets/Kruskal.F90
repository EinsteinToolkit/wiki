! $Header$

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



module kruskal
  use cctk
  use constants
  use tensor
  implicit none
  
  integer, parameter :: slicing_Kruskal        = 1
  integer, parameter :: slicing_Schwarzschild  = 2
  integer, parameter :: slicing_Wald_Iyer      = 3
  
contains
  
  function decode_slicing (slicing) result (slicing_type)
    DECLARE_CCTK_FUNCTIONS
    
    CCTK_POINTER_TO_CONST, intent(in) :: slicing
    integer :: slicing_type
    
    if (CCTK_EQUALS (slicing, "Kruskal")) then
       slicing_type = slicing_Kruskal
    else if (CCTK_EQUALS (slicing, "Schwarzschild")) then
       slicing_type = slicing_Schwarzschild
    else if (CCTK_EQUALS (slicing, "Wald-Iyer")) then
       slicing_type = slicing_Wald_Iyer
    else
       call CCTK_WARN (0, "internal error")
    end if
  end function decode_slicing
  
  
  
  subroutine calc_kruskal (xx, tt, alfa, beta, gama, uu, vv, sr, st, &
       slicing_type, allow_singular_slicing)
    
    DECLARE_CCTK_PARAMETERS
    
    CCTK_REAL, intent(in)  :: xx(3), tt
    CCTK_REAL, intent(out) :: alfa, beta(3), gama(3,3)
    CCTK_REAL, intent(out) :: uu, vv, sr, st
    integer,   intent(in)  :: slicing_type
    logical,   intent(in)  :: allow_singular_slicing
    
    CCTK_REAL :: rr, rr01, rr02
    CCTK_REAL :: ww
    CCTK_REAL :: guu, gvv, gqq
    CCTK_REAL :: duudww, dwwdrr, duudrr, duudtt, dvvdrr, dvvdtt
    CCTK_REAL :: grr, grt, gtt
    CCTK_REAL :: betal(3), detgama, gamau(3,3)
    
    integer   :: a, b
  
  
  
    interface
       ! Solve y = (x-1) exp (x) for x
       CCTK_REAL function Kruskal_Solve (y)
         implicit none
         CCTK_REAL y
       end function Kruskal_Solve
    end interface
    
    
    
    ! Choose a nonzero square of the coordinate radius
    rr = sqrt (sum (xx**2))
    rr02 = rr**2 + eps**2
    rr01 = sqrt (rr02)
    
    
    
    ! u is the Kruskal radial coordinate
    ! Choose it such that u in [-inf; +inf]
    ! u := sqrt(w exp (w)) - sqrt(1/w exp (1/w))
    ! with w := r / 2M
    ww = rr01 / (2*mass)
    if (ww > 100 .or. 1/ww > 100) goto 9999
    uu = sqrt(ww * exp (ww)) - sqrt(1/ww * exp (1/ww))
    
    ! v is the Kruskal time coordinate
    select case (slicing_type)
    case (slicing_Kruskal)
       vv = tt
    case (slicing_Schwarzschild)
       vv = uu * tanh (tt / (4 * mass))
    case (slicing_Wald_Iyer)
       vv = tt * xx(3) / rr01
    case default
       call CCTK_WARN (0, "internal error")
    end select
    
    
    
    ! r is the Schwarzschild radial coordinate
    ! t is the Schwarzschild time
    
    ! The singularities are at r=0
    ! The horizons are at r=2M, which are u=+v and u=-v
    
    ! v = (+/-) sqrt(|r/2M - 1|) exp(r/4M) sinh(t/4M)
    ! u = (+/-) sqrt(|r/2M - 1|) exp(r/4M) cosh(t/4M)
    
    ! u^2 - v^2 =: exp(r/2M) (r/2M - 1)
    ! v / u = tanh(t/4M)
    
    sr = 2*mass * Kruskal_Solve (uu**2 - vv**2)
    if (abs (uu) < 1.0d-8) then
       ! TODO: this is not good, use an analytic expression for v/u instead
       st = 0
    else
       st = 4*mass * atanh (vv / uu)
    end if
    
    
    
    ! The Kruskal metric is
    ! ds^2 = 32 M^3/r exp(-r/2M) (- dv^2 + du^2) + r^2 dOmega^2
    guu = 32 * mass**3 / sr * exp(- sr / (2*mass))
    gvv = - guu
    gqq = sr**2
    
    
    
    ! Convert to coordinates
    duudww = (ww * sqrt(ww * exp (ww)) + sqrt(1/ww * exp (1/ww))) &
         &   * (ww + 1) / (2 * ww**2)
    dwwdrr = 1 / (2*mass)
    duudrr = duudww * dwwdrr
    duudtt = 0
    
    select case (slicing_type)
    case (slicing_Kruskal)
       ! vv = tt
       dvvdrr = 0
       dvvdtt = 1
    case (slicing_Schwarzschild)
       ! vv = uu * tanh (tt / (4 * mass))
       if (allow_singular_slicing) then
          dvvdrr = duudrr * tanh (tt / (4 * mass))
          dvvdtt = uu / (4*mass * cosh (tt / (4 * mass)) ** 2)
       else
          ! vv = uu * tanh (tt0 / (4 * mass)) + (tt - tt0)
          dvvdrr = duudrr * tanh (tt / (4 * mass))
          dvvdtt = 1
       end if
    case (slicing_Wald_Iyer)
       ! vv = tt * xx(3) / rr01
       dvvdrr = 0
       dvvdtt = xx(3) / rr01
#warning "TODO: dv/dq is missing"
    case default
       call CCTK_WARN (0, "internal error")
    end select
    
    
    
    grr = duudrr * duudrr * guu + dvvdrr * dvvdrr * gvv
    grt = duudrr * duudtt * guu + dvvdrr * dvvdtt * gvv
    gtt = duudtt * duudtt * guu + dvvdtt * dvvdtt * gvv
    
    
    
    ! Calculate metric
#warning "TODO: beta^2 term is missing"
    alfa = sign (sqrt(- gtt), dvvdtt)
    
    forall (a=1:3)
       betal(a) = grt * (xx(a) / rr01)
    end forall
    
    forall (a=1:3, b=1:3)
       gama(a,b) = gqq / rr02 * delta3(a,b) &
            &      + (grr - gqq / rr02) * xx(a) * xx(b) / rr02
    end forall
    
    call calc_det (gama, detgama)
    call calc_inv (gama, detgama, gamau)
    
    beta = matmul (gamau, betal)
    
    ! All is fine, skip error handling
    goto 1000
    
    
    
    ! Error handling
9999 continue
    
    ! Invent arbitrary data
    alfa = 1
    beta = 0
    gama = delta3
    
    uu = 0
    vv = 0
    sr = 0
    st = 0
    
    
    
    ! Continue whether there was an error or not
1000 continue
  end subroutine calc_kruskal
  
  elemental function atanh (x) result (r)
    CCTK_REAL, intent(in) :: x
    CCTK_REAL :: r
    r = log ((1 + x) / (1 - x)) / 2
  end function atanh
  
end module kruskal
