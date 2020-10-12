! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

! The structure and declarations in this file might seem odd, but allow a simple wrapping of specfab.f90 into Python module using f2py

module specfabpy
    
    use specfab    

    implicit none
    
    integer, parameter :: num_dofs = nlm_len
    integer, parameter :: Ltrunc = Lcap
    
contains

    subroutine init(nlmlen, lm_list) 
        implicit none
        integer, intent(in) :: nlmlen
        integer, intent(out) :: lm_list(2,nlmlen)
        call initspecfab()
        lm_list = lm 
    end
    
    function get_dndt_ij(nlmlen, eps,omg)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        integer, intent(in) :: nlmlen
        real(kind=dp), intent(in) :: eps(3,3), omg(3,3) ! strain-rate (eps) and spin (omg)
        complex(kind=dp) :: get_dndt_ij(nlmlen,nlmlen)
        get_dndt_ij = dndt_ij(eps,omg)
    end
    
    subroutine get_eigenframe(nlm, e1,e2,e3, eigvals)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
!        integer, intent(in) :: nlmlen
!        complex(kind=dp), intent(in)  :: nlm(nlmlen)
        complex(kind=dp), intent(in)  :: nlm(:) ! f2py does not need to know the dimension
        real(kind=dp), dimension(3), intent(out) :: e1,e2,e3, eigvals
        call eigenframe(nlm, e1,e2,e3, eigvals)
    end
    
    subroutine get_pqframe(nlm, p23,p12,p13, q23,q12,q13)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), dimension(3), intent(out) :: p23,p12,p13, q23,q12,q13
        call pqframe(nlm, p23,p12,p13, q23,q12,q13)
    end
    
    function get_Eeiej(nlm, Ecc,Eca,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: Ecc, Eca
        integer, intent(in) :: nprime
        real(kind=dp) :: get_Eeiej(3,3)
        get_Eeiej = Eeiej(nlm, Ecc,Eca,nprime)
    end

    function get_Epijqij(nlm, Ecc,Eca,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: Ecc, Eca
        integer, intent(in) :: nprime    
        real(kind=dp) :: get_Epijqij(3)
        get_Epijqij = Epijqij(nlm, Ecc,Eca,nprime)
    end
    
    function get_nu(eps)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        real(kind=dp), intent(in) :: eps(3,3)
        real(kind=dp) :: get_nu
        get_nu = f_nu(eps)
    end
        
end module specfabpy 

