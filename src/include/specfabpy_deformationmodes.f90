! General 

function F_to_strain(F) result (strain)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: F(3,3)      ! deformation gradient
    real(kind=dp)             :: strain(3,3) ! strain tensor 
    strain = F_to_strain__sf(F)
end

subroutine ugrad_to_D_and_W(ugrad, D,W)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in)  :: ugrad(3,3)
    real(kind=dp), intent(out) :: D(3,3), W(3,3)
    call ugrad_to_D_and_W__sf(ugrad, D,W)
end

! Pure shear

function pureshear_r(tau, time) result (r)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: tau, time
    real(kind=dp)             :: r
    r = pureshear_r__sf(tau, time)
end

function pureshear_strainii_to_t(strainii, tau) result (time)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: strainii, tau
    real(kind=dp)             :: time
    time = pureshear_strainii_to_t__sf(strainii, tau)
end

function pureshear_F(ax,q,tau, time) result (F)
    use specfabpy_const
    implicit none
    integer, intent(in) :: ax
    real(kind=dp), intent(in) :: q,tau,time
    real(kind=dp)             :: F(3,3)
    F = pureshear_F__sf(ax,q,tau, time)
end
   
function pureshear_ugrad(ax,q,tau) result(ugrad)
    use specfabpy_const
    implicit none
    integer, intent(in)       :: ax
    real(kind=dp), intent(in) :: q,tau
    real(kind=dp)             :: ugrad(3,3)
    ugrad = pureshear_ugrad__sf(ax,q,tau)
end

! Simple shear 

function simpleshear_gamma(tau, time) result (gam)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: tau, time
    real(kind=dp)             :: gam
    gam = simpleshear_gamma__sf(tau, time)
end

function simpleshear_gamma_to_t(gam, tau) result (time)
    ! Time taken to reach shear angle gam
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: gam,tau
    real(kind=dp)             :: time
    time = simpleshear_gamma_to_t__sf(gam, tau)
end

function simpleshear_F(plane,tau, time) result (F)
    use specfabpy_const
    implicit none
    integer, intent(in)       :: plane
    real(kind=dp), intent(in) :: tau,time
    real(kind=dp)             :: F(3,3)
    F = simpleshear_F__sf(plane,tau, time)
end
   
function simpleshear_ugrad(plane,tau) result(ugrad)
    use specfabpy_const
    implicit none
    integer, intent(in)       :: plane
    real(kind=dp), intent(in) :: tau
    real(kind=dp)             :: ugrad(3,3)
    ugrad = simpleshear_ugrad__sf(plane,tau)
end

