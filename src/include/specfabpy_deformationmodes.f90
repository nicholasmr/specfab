! General 

function F_to_strain(F) result (strain)
    implicit none
    real(kind=dp), intent(in) :: F(3,3)      ! deformation gradient
    real(kind=dp)             :: strain(3,3) ! strain tensor 
    strain = F_to_strain__sf(F)
end

subroutine ugrad_to_D_and_W(ugrad, D,W)
    implicit none
    real(kind=dp), intent(in)  :: ugrad(3,3)
    real(kind=dp), intent(out) :: D(3,3), W(3,3)
    call ugrad_to_D_and_W__sf(ugrad, D,W)
end

! Pure shear

function pureshear_b(T, time) result (b)
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: T, time
    real(kind=dp)             :: b
    b = pureshear_b__sf(T, time)
end

function pureshear_strainii_to_t(strainii, T) result (time)
    ! Time taken to reach strainii strain for ax=i:
    !   strainii = Fii - 1 = b^(-1) - 1 = exp(-t/T) - 1 => exp(-t/T) = strainii + 1
    use specfabpy_const
    implicit none
    real(kind=dp), intent(in) :: strainii,T
    real(kind=dp)             :: time
    time = pureshear_strainii_to_t__sf(strainii, T)
end

function pureshear_F(ax,r,T, time) result (F)
    use specfabpy_const
    implicit none
    integer, intent(in) :: ax
    real(kind=dp), intent(in) :: r,T,time
    real(kind=dp)             :: F(3,3)
    F = pureshear_F__sf(ax,r,T, time)
end
   
function pureshear_ugrad(ax,r,T) result(ugrad)
    use specfabpy_const
    implicit none
    integer, intent(in)       :: ax
    real(kind=dp), intent(in) :: r,T
    real(kind=dp)             :: ugrad(3,3)
    ugrad = pureshear_ugrad__sf(ax,r,T)
end

! Simple shear 

function simpleshear_gamma(T, time) result (gam)
    implicit none
    real(kind=dp), intent(in) :: T, time
    real(kind=dp)             :: gam
    gam = simpleshear_gamma__sf(T, time)
end

function simpleshear_gamma_to_t(gam, T) result (time)
    ! Time taken to reach shear angle gam
    implicit none
    real(kind=dp), intent(in) :: gam,T
    real(kind=dp)             :: time
    time = simpleshear_gamma_to_t__sf(gam, T)
end

function simpleshear_F(plane,T, time) result (F)
    implicit none
    integer, intent(in)       :: plane
    real(kind=dp), intent(in) :: T,time
    real(kind=dp)             :: F(3,3)
    F = simpleshear_F__sf(plane,T, time)
end
   
function simpleshear_ugrad(plane,T) result(ugrad)
    implicit none
    integer, intent(in)       :: plane
    real(kind=dp), intent(in) :: T
    real(kind=dp)             :: ugrad(3,3)
    ugrad = simpleshear_ugrad__sf(plane,T)
end

