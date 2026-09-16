! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine Kruskal_InitialData (CCTK_ARGUMENTS)
  use adm_metric
  use cctk
  use kruskal
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL :: xx(3), tt
  CCTK_REAL :: alfa, beta(3), gama(3,3)
  CCTK_REAL :: uu, vv, sr, st
  CCTK_REAL :: xx1(3), tt1
  CCTK_REAL :: alfa1, beta1(3), gama1(3,3)
  CCTK_REAL :: xx2(3), tt2
  CCTK_REAL :: alfa2, beta2(3), gama2(3,3)
  CCTK_REAL :: dtalfa, dtbeta(3), dtgama(3,3)
  CCTK_REAL :: dalfa(3), dbeta(3,3), dgama(3,3,3)
  CCTK_REAL :: kk(3,3)
  
  integer   :: i, j, k
  integer   :: d
  
  logical   :: want_lapse, want_shift, want_data
  integer   :: slicing_type
  
  
  
  want_lapse = CCTK_EQUALS (initial_lapse, "Kruskal")
  want_shift = CCTK_EQUALS (initial_shift, "Kruskal")
  want_data  = CCTK_EQUALS (initial_data,  "Kruskal")
  
  
  
  slicing_type = decode_slicing (slicing)
  
  
  
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           
           ! Get the current position
           xx(1) = x(i,j,k)
           xx(2) = y(i,j,k)
           xx(3) = z(i,j,k)
           tt = cctk_time
           
           
           
           ! Calculate spatial derivatives
           do d = 1, 3
              xx1 = xx
              xx2 = xx
              xx1(d) = xx1(d) - delta
              xx2(d) = xx2(d) + delta
              call calc_kruskal &
                   (xx1, tt, alfa1, beta1, gama1, uu, vv, sr, st, &
                   slicing_type, .false.)
              call calc_kruskal &
                   (xx2, tt, alfa2, beta2, gama2, uu, vv, sr, st, &
                   slicing_type, .false.)
              dalfa(d)     = (alfa2 - alfa1) / (2 * delta)
              dbeta(:,d)   = (beta2 - beta1) / (2 * delta)
              dgama(:,:,d) = (gama2 - gama1) / (2 * delta)
           end do
           
           ! Calculate temporal derivatives
           tt1 = tt - delta
           tt2 = tt + delta
           call calc_kruskal (xx, tt1, alfa1, beta1, gama1, uu, vv, sr, st, &
                slicing_type, .false.)
           call calc_kruskal (xx, tt2, alfa2, beta2, gama2, uu, vv, sr, st, &
                slicing_type, .false.)
           dtalfa = (alfa2 - alfa1) / (2 * delta)
           dtbeta = (beta2 - beta1) / (2 * delta)
           dtgama = (gama2 - gama1) / (2 * delta)
           
           call calc_kruskal (xx, tt, alfa, beta, gama, uu, vv, sr, st, &
                slicing_type, .false.)
           
           call calc_extcurv (gama, dgama, dtgama, alfa, beta, dbeta, kk)
           
           ! Re-evaluate the metric with the correct slicing
           call calc_kruskal (xx, tt, alfa, beta, gama, uu, vv, sr, st, &
                slicing_type, .true.)
           
           
           
           ! Store initial data
           
           if (want_data) then
              gxx(i,j,k) = gama(1,1)
              gxy(i,j,k) = gama(1,2)
              gxz(i,j,k) = gama(1,3)
              gyy(i,j,k) = gama(2,2)
              gyz(i,j,k) = gama(2,3)
              gzz(i,j,k) = gama(3,3)
              
              kxx(i,j,k) = kk(1,1)
              kxy(i,j,k) = kk(1,2)
              kxz(i,j,k) = kk(1,3)
              kyy(i,j,k) = kk(2,2)
              kyz(i,j,k) = kk(2,3)
              kzz(i,j,k) = kk(3,3)
           end if
           
           if (want_lapse) then
              alp(i,j,k) = alfa
           end if
           
           if (want_shift) then
              betax(i,j,k) = beta(1)
              betay(i,j,k) = beta(2)
              betaz(i,j,k) = beta(3)
           end if
           
           if (calculate_Schwarzschild_radius /= 0) then
              schw_radius(i,j,k) = sr
              schw_time  (i,j,k) = st
              kruskal_u  (i,j,k) = uu
              kruskal_v  (i,j,k) = vv
           end if
           
        end do
     end do
  end do
  
  
  
  if (CCTK_EQUALS (metric_type, "physical")) then
     ! do nothing
  else if (CCTK_EQUALS (metric_type, "static conformal")) then
     if (CCTK_EQUALS (conformal_storage, "factor")) then
        psi = 1
        conformal_state = 1
     else if (CCTK_EQUALS (conformal_storage, "factor+derivs")) then
        psi = 1
        psix = 0
        psiy = 0
        psiz = 0
        conformal_state = 2
     else if (CCTK_EQUALS (conformal_storage, "factor+derivs+2nd derivs")) then
        psi = 1
        psix = 0
        psiy = 0
        psiz = 0
        psixx = 0
        psixy = 0
        psixz = 0
        psiyy = 0
        psiyz = 0
        psizz = 0
        conformal_state = 3
     else
        call CCTK_WARN (0, "Unknown value of StaticConformal::conformal_storage")
     end if
  else
     call CCTK_WARN (0, "Unknown value of ADMBase::metric_type")
  end if
  
end subroutine Kruskal_InitialData
