! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

! Routines for rotating distributions 

module rotation  

    use header
    use moments 
    
    implicit none 
        
contains      

    function rotate_nlm4(nlm4, theta,phi) result (nlm4_rot)
    
        implicit none

        complex(kind=dp), intent(in) :: nlm4(nlm_lenvec(4)) ! nlm truncated at L=4
        real(kind=dp), intent(in)    :: theta, phi 
        complex(kind=dp)             :: nlm4_rot(nlm_lenvec(4)) !, n2m_rot(5), n4m_rot(9)
        complex(kind=dp)             :: D2mn(-2:2,-2:2), D4mn(-4:4,-4:4)
    
        include "include/Dlmn_L4.f90"
        nlm4_rot(L0rng) = nlm4(L0rng)
        nlm4_rot(L2rng) = matmul(D2mn, nlm4(L2rng))
        nlm4_rot(L4rng) = matmul(D4mn, nlm4(L4rng))
    end
    
    function rotate_nlm(nlm, theta,phi) result (nlm_rot)
    
        implicit none

        complex(kind=dp), intent(in) :: nlm(:) ! nlm truncated at L=12
        real(kind=dp), intent(in)    :: theta, phi 
        complex(kind=dp)             :: nlm_rot(size(nlm))
        integer, parameter           :: Lmax_wignerD = 12
        integer                      :: nlm_len
        complex(kind=dp)             :: D2mn(-2:2,-2:2), D4mn(-4:4,-4:4), D6mn(-6:6,-6:6), D8mn(-8:8,-8:8), D10mn(-10:10,-10:10), D12mn(-12:12,-12:12) !, D14mn(-14:14,-14:14)
    
        nlm_rot(:) = 0.0d0
        nlm_len = size(nlm) ! for convenience
    
        include "include/Dlmn.f90"
        nlm_rot(L0rng) = nlm(L0rng)
        if (nlm_len .ge. nlm_lenvec(2))  nlm_rot(L2rng)  = matmul(D2mn,  nlm(L2rng))
        if (nlm_len .ge. nlm_lenvec(4))  nlm_rot(L4rng)  = matmul(D4mn,  nlm(L4rng))
        if (nlm_len .ge. nlm_lenvec(6))  nlm_rot(L6rng)  = matmul(D6mn,  nlm(L6rng))
        if (nlm_len .ge. nlm_lenvec(8))  nlm_rot(L8rng)  = matmul(D8mn,  nlm(L8rng))
        if (nlm_len .ge. nlm_lenvec(10)) nlm_rot(L10rng) = matmul(D10mn, nlm(L10rng))
        if (nlm_len .ge. nlm_lenvec(12)) nlm_rot(L12rng) = matmul(D12mn, nlm(L12rng))

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

