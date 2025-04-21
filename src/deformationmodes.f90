! N. M. Rathmann <rathmann@nbi.ku.dk>, 2023-

! Idealized kinematic modes of deformation 

!----------
! Notation:
!----------
! ugrad_ij := d(u_i)/d(x_j) 
!----------

module deformationmodes 

    use header
    use tensorproducts

    implicit none 
        
contains      

    ! General
    
    function F_to_strain(F) result (strain)
        implicit none
        real(kind=dp), intent(in) :: F(3,3)      ! deformation gradient
        real(kind=dp)             :: strain(3,3) ! strain tensor 
        strain = symmetricpart(F) - identity
    end
    
    subroutine ugrad_to_D_and_W(ugrad, D,W)
        implicit none
        real(kind=dp), intent(in)  :: ugrad(3,3)
        real(kind=dp), intent(out) :: D(3,3), W(3,3)
        D = symmetricpart(ugrad)
        W = antisymmetricpart(ugrad)
    end
    
    ! Pure shear

    function pureshear_r(tau, time) result (r)
        implicit none
        real(kind=dp), intent(in) :: tau, time
        real(kind=dp)             :: r
        r = exp(-time/tau) ! scaling parameter
    end
    
    function pureshear_strainii_to_t(strainii, tau) result (time)
        ! Time taken to reach strainii strain for ax=i:
        !   strainii = Fii - 1 = r - 1 = exp(-t/tau) - 1 => exp(-t/tau) = strainii + 1
        implicit none
        real(kind=dp), intent(in) :: strainii,tau
        real(kind=dp)             :: time
        time = -tau*log(strainii+1.0d0)
    end
    
    function pureshear_F(ax,q,tau, time) result (F)
        implicit none
        integer, intent(in) :: ax
        real(kind=dp), intent(in) :: q,tau,time
        real(kind=dp)             :: r, F(3,3), qp,qm
        r = pureshear_r(tau,time)
        qp = -(1+q)/2.0d0
        qm = -(1-q)/2.0d0
        if (ax == 0) F = diag([r,     r**qp, r**qm]) ! compression along x
        if (ax == 1) F = diag([r**qm, r,     r**qp]) ! compression along y
        if (ax == 2) F = diag([r**qp, r**qm, r])     ! compression along z
    end
       
    function pureshear_ugrad(ax,q,tau) result(ugrad)
        implicit none
        integer, intent(in)       :: ax
        real(kind=dp), intent(in) :: q,tau
        real(kind=dp)             :: ugrad(3,3), qp,qm
        qp = -(1+q)/2.0d0
        qm = -(1-q)/2.0d0
        if (ax == 0) ugrad = -1/tau * diag([1.0d0, qp,    qm])    ! compression along x
        if (ax == 1) ugrad = -1/tau * diag([qm,    1.0d0, qp])    ! compression along y
        if (ax == 2) ugrad = -1/tau * diag([qp,    qm,    1.0d0]) ! compression along z
    end

    ! Simple shear

    function simpleshear_gamma(tau, time) result (gam)
        implicit none
        real(kind=dp), intent(in) :: tau, time
        real(kind=dp)             :: gam
        gam = atan(time/tau) ! shear angle
    end
    
    function simpleshear_gamma_to_t(gam, tau) result (time)
        ! Time taken to reach shear angle gam
        implicit none
        real(kind=dp), intent(in) :: gam,tau
        real(kind=dp)             :: time
        time = tau*tan(gam)
    end
    
    function simpleshear_F(plane,tau, time) result (F)
        implicit none
        integer, intent(in)       :: plane
        real(kind=dp), intent(in) :: tau,time
        real(kind=dp)             :: gam, F(3,3)
        gam = simpleshear_gamma(tau, time)
        F = identity
        if (plane == 0) F(2,3) = tan(gam) ! yz shear (du_y/dz)
        if (plane == 1) F(1,3) = tan(gam) ! xz shear (du_x/dz)
        if (plane == 2) F(1,2) = tan(gam) ! xy shear (du_x/dy)
    end
       
    function simpleshear_ugrad(plane,tau) result(ugrad)
        implicit none
        integer, intent(in)       :: plane
        real(kind=dp), intent(in) :: tau
        real(kind=dp)             :: ugrad(3,3)
        ugrad = 0.0d0
        if (plane == 0) ugrad(2,3) = 1/tau ! yz shear
        if (plane == 1) ugrad(1,3) = 1/tau ! xz shear
        if (plane == 2) ugrad(1,2) = 1/tau ! xy shear
    end

end module deformationmodes
