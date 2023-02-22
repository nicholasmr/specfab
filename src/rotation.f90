! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

! Routines for rotation ODFs 

module rotation  

    use moments 
    use dynamics ! for nlm_lenvec, Lcap__max, etc.
    
    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real,    parameter, private :: Pi = 3.141592653589793
    
    integer, private :: ll ! Loop indices
        
contains      

    function rotate_nlm4(nlm4, theta,phi) result (nlm4_rot)
    
        implicit none

        complex(kind=dp), intent(in) :: nlm4(15) ! nlm truncated at L=4
        real(kind=dp), intent(in)    :: theta, phi 
        complex(kind=dp)             :: nlm4_rot(15) !, n2m_rot(5), n4m_rot(9)
        complex(kind=dp)             :: D2mn(-2:2,-2:2), D4mn(-4:4,-4:4)
    
        include "include/Dlmn_L4.f90"
        nlm4_rot(I_l0) = nlm4(I_l0)
        nlm4_rot(I_l2:(I_l4-1)) = matmul(D2mn, nlm4(I_l2:(I_l4-1)))
        nlm4_rot(I_l4:(I_l6-1)) = matmul(D4mn, nlm4(I_l4:(I_l6-1)))
    end
    
    function rotate_nlm(nlm, theta,phi) result (nlm_rot)
    
        implicit none

        complex(kind=dp), intent(in) :: nlm(:) ! nlm truncated at L=12
        real(kind=dp), intent(in)    :: theta, phi 
        complex(kind=dp)             :: nlm_rot(size(nlm))
        integer, parameter           :: Lmax_wignerD = 12
        integer                      :: nlm_len
        complex(kind=dp)             :: D2mn(-2:2,-2:2), D4mn(-4:4,-4:4), D6mn(-6:6,-6:6), D8mn(-8:8,-8:8), D10mn(-10:10,-10:10), D12mn(-12:12,-12:12) !, D14mn(-14:14,-14:14)
        integer, parameter           :: nlm_lenvec(0:Lcap__max) = [((ll+1)*(ll+2)/2, ll=0, Lcap__max, 1)] ! nlm length for a given Lcap
    
        nlm_len = size(nlm) ! for convenience
    
        include "include/Dlmn.f90"
        nlm_rot(I_l0) = nlm(I_l0)
        if (nlm_len .ge. nlm_lenvec(2))  nlm_rot(I_l2:(I_l4-1))   = matmul(D2mn,  nlm(I_l2:(I_l4-1)))
        if (nlm_len .ge. nlm_lenvec(4))  nlm_rot(I_l4:(I_l6-1))   = matmul(D4mn,  nlm(I_l4:(I_l6-1)))
        if (nlm_len .ge. nlm_lenvec(6))  nlm_rot(I_l6:(I_l8-1))   = matmul(D6mn,  nlm(I_l6:(I_l8-1)))
        if (nlm_len .ge. nlm_lenvec(8))  nlm_rot(I_l8:(I_l10-1))  = matmul(D8mn,  nlm(I_l8:(I_l10-1)))
        if (nlm_len .ge. nlm_lenvec(10)) nlm_rot(I_l10:(I_l12-1)) = matmul(D10mn, nlm(I_l10:(I_l12-1)))
        if (nlm_len .ge. nlm_lenvec(12)) nlm_rot(I_l12:(I_l14-1)) = matmul(D12mn, nlm(I_l12:(I_l14-1)))
!       if (nlm_len .ge. nlm_lenvec(6)) nlm_rot(I_l14:(I_l16-1)) = matmul(D14mn, nlm(I_l14:(I_l16-1)))

!        if (nlm_len .gt. nlm_lenvec(Lmax_wignerD)) then
!            print *, 'specfab error: L > L_max for rotating nlm with wigner D matrices'
!        end if
    end
    
    function rotate_vector(v, theta, phi) result(w)
    
        ! theta is angle rotating z->x, i.e. around y
        ! phi   is angle rotating x->y, i.e. around z
        
        implicit none
        
        real(kind=dp), intent(in) :: v(3), theta, phi
        real(kind=dp) :: w(3), Rz(3,3), Ry(3,3)
        
        Rz = reshape([+cos(phi),+sin(phi),0.0d0, & ! column major, so this row is Rz(:,1)
                      -sin(phi),+cos(phi),0.0d0, &
                          0.0d0,0.0d0,1.0d0],[3,3])
                          
        Ry = reshape([ +cos(theta),0.0d0,-sin(theta), & ! column major, so this row is Ry(:,1)
                            0.0d0,1.0d0,0.0d0, &
                       +sin(theta),0.0d0,+cos(theta)],[3,3])
                       
        w = matmul(Ry,matmul(Rz,v))
    end

end module rotation

