! N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

! Idealized modes of deformation

module deformationmodes 

    use tensorproducts

    implicit none 

    integer, parameter, private       :: dp = 8 ! Default precision
    real(kind=dp), parameter, private :: identity(3,3) = reshape([1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0], [3,3])
        
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

    function pureshear_b(T, time) result (b)
        implicit none
        real(kind=dp), intent(in) :: T, time
        real(kind=dp)             :: b
        b = exp(time/T) ! scaling parameter
    end
    
    function pureshear_strainii_to_t(strainii, T) result (time)
        ! Time taken to reach strainii strain for ax=i:
        !   strainii = Fii - 1 = b^(-1) - 1 = exp(-t/T) - 1 => exp(-t/T) = strainii + 1
        implicit none
        real(kind=dp), intent(in) :: strainii,T
        real(kind=dp)             :: time
        time = -T*log(strainii+1.0d0)
    end
    
    function pureshear_F(ax,r,T, time) result (F)
        implicit none
        integer, intent(in) :: ax
        real(kind=dp), intent(in) :: r,T,time
        real(kind=dp)             :: b, F(3,3)
        b = pureshear_b(T,time)
        if (ax == 0) F = diag([b**(-1), b**((1+r)/2), b**((1-r)/2)]) ! compression along x
        if (ax == 1) F = diag([b**((1-r)/2), b**(-1), b**((1+r)/2)]) ! compression along y
        if (ax == 2) F = diag([b**((1+r)/2), b**((1-r)/2), b**(-1)]) ! compression along z
    end
       
    function pureshear_ugrad(ax,r,T) result(ugrad)
        implicit none
        integer, intent(in)       :: ax
        real(kind=dp), intent(in) :: r,T
        real(kind=dp)             :: ugrad(3,3)
        if (ax == 0) ugrad = 1/T*diag([-1.0d0, (1+r)/2, (1-r)/2]) ! compression along x
        if (ax == 1) ugrad = 1/T*diag([(1-r)/2, -1.0d0, (1+r)/2]) ! compression along y
        if (ax == 2) ugrad = 1/T*diag([(1+r)/2, (1-r)/2, -1.0d0]) ! compression along z
    end


    ! Simple shear

    function simpleshear_gamma(T, time) result (gam)
        implicit none
        real(kind=dp), intent(in) :: T, time
        real(kind=dp)             :: gam
        gam = atan(time/T) ! shear angle
    end
    
    function simpleshear_gamma_to_t(gam, T) result (time)
        ! Time taken to reach shear angle gam
        implicit none
        real(kind=dp), intent(in) :: gam,T
        real(kind=dp)             :: time
        time = T*tan(gam)
    end
    
    function simpleshear_F(plane,T, time) result (F)
        implicit none
        integer, intent(in)       :: plane
        real(kind=dp), intent(in) :: T,time
        real(kind=dp)             :: gam, F(3,3)
        gam = simpleshear_gamma(T, time)
        F = identity
        if (plane == 0) F(2,3) = tan(gam) ! y--z shear
        if (plane == 1) F(1,3) = tan(gam) ! x--z shear
        if (plane == 2) F(1,2) = tan(gam) ! x--y shear
    end
       
    function simpleshear_ugrad(plane,T) result(ugrad)
        implicit none
        integer, intent(in)       :: plane
        real(kind=dp), intent(in) :: T
        real(kind=dp)             :: ugrad(3,3)
        ugrad = 0.0d0
        if (plane == 0) ugrad(2,3) = 1/T ! y--z shear
        if (plane == 1) ugrad(1,3) = 1/T ! x--z shear
        if (plane == 2) ugrad(1,2) = 1/T ! x--y shear
    end

end module deformationmodes
