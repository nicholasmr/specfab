! N. M. Rathmann <rathmann@nbi.ku.dk>, 2023-

! Idealized CPO states

module idealstate 

    use header
    use rotation
    use frames

    implicit none 
    
contains

    function nlm_idealz(colat, Lmax) 
    
        ! delta function (if colat=0) or ring delta function (if colat>0) aligned with the vertical direction z (and antipodally symmetric)
    
        implicit none
        
        real(kind=dp), intent(in) :: colat ! colat of delta
        integer, intent(in)       :: Lmax
        complex(kind=dp)          :: nlm_idealz( nlm_lenvec(Lmax) )
        integer, parameter        :: N = 10
        real(kind=dp)             :: Yl0_list(0:N)
        
        nlm_idealz(:) = 0.0d0
        Yl0_list = [Y00(colat), Y20(colat), Y40(colat), Y60(colat), Y80(colat), Y100(colat), Y120(colat), Y140(colat), Y160(colat), Y180(colat), Y200(colat)]
        do ii=0,min(N,int(Lmax/2))
            nlm_idealz(ILm(ii*2)) = Yl0_list(ii)
        end do
    end

    function nlm_ideal(m, colat, Lmax) result(nlm_rot)
    
        implicit none
        
        real(kind=dp), intent(in) :: m(3), colat ! symmetry axis, circle colat w.r.t. m
        integer, intent(in)       :: Lmax
        complex(kind=dp)          :: nlm(nlm_lenvec(Lmax)), nlm_rot(nlm_lenvec(Lmax))
        integer, parameter        :: N = 10
        real(kind=dp)             :: theta, phi
        
        nlm(:) = nlm_idealz(colat, Lmax)
        theta = acos(m(3))
        phi = atan2(m(2),m(1))
        nlm_rot = rotate_nlm(rotate_nlm(nlm, theta,0.0d0), 0.0d0,phi)
    end
    
    function nlm_isvalid(nhat20, nhat40) result(isvalid)
        
        ! Test whether normalized, low-order nlm components are within eigenvalue bounds
        
        implicit none
        
        complex(kind=dp), intent(in) :: nhat20(:), nhat40(:) 
!        integer                      :: N = 
        logical                      :: isvalid(size(nhat20)), isvalid_a2, idvalid_a4
        
        complex(kind=dp)             :: nlm(nlm_lenvec(4))
        real(kind=dp)                :: ei(3,3), a2_eigvals(3)
        real(kind=dp)                :: Q1(3,3),Q2(3,3),Q3(3,3),Q4(3,3),Q5(3,3),Q6(3,3), a4_eigvals(6)

        do ii=1,size(nhat20)
            nlm(:) = 0.0d0
            nlm(1) = 1.0d0
            nlm(ILm(2)) = nhat20(ii)
            nlm(ILm(4)) = nhat40(ii)
            
            call eig(nlm, ei, a2_eigvals)
            isvalid_a2 = (maxval(a2_eigvals) < 1.0d0) .and. (minval(a2_eigvals) > 0.0d0) ! a2 eigenvalues are within bounds
            
            call a4_eigentensors(nlm, Q1,Q2,Q3,Q4,Q5,Q6, a4_eigvals)  
            idvalid_a4 = (maxval(a4_eigvals) < 1.0d0) .and. (minval(a4_eigvals) > 0.0d0) ! a4 eigenvalues are within bounds
            
            isvalid(ii) = isvalid_a2 .and. idvalid_a4
        end do                
    end 

    function Sl(nlm, l)
    
        ! Angular power spectrum S(l) at single degree l
        
        implicit none

        complex(kind=dp), intent(in) :: nlm(:)
        integer, intent(in)          :: l
        real(kind=dp)                :: Sl
        complex(kind=dp)             :: nlm_sub(2*l+1)
        
        nlm_sub(:) = nlm(IL(l):(IL(l+2)-1)) ! [n_l^-l, ... n_l^l]
        Sl = 1.0d0/(2*l+1) * sum(abs(nlm_sub)**2) 
    end
    
    subroutine powerspectrum(nlm,Lmax, S,l)
    
        ! Normalized angular power spectrum S(l)/S(0) for all degrees l <= Lmax
        
        implicit none

        complex(kind=dp), intent(in) :: nlm(:)
        integer, intent(in)          :: Lmax
        real(kind=dp), intent(out)   :: S(Lmax/2+1) ! S(l)
        integer, intent(out)         :: l(Lmax/2+1) ! list of l values (e.g. 0,2,4 for L=4)
        
        l(:) = [(ll, ll=0, Lmax, 2)] ! (0, 2, 4, ..., L)
        S(:) = [(Sl(nlm,l(ii)), ii=1, size(l))]
        S(:) = S(:)/S(1) ! normalize
    end

    function pfJ(nlm, Lmax) 
    
        ! Pole figure J index, truncated at Lmax
        
        implicit none

        real(kind=dp)                :: pfJ
        complex(kind=dp), intent(in) :: nlm(:)
        integer, intent(in)          :: Lmax
        integer                      :: l(Lmax/2+1) ! list of l values (e.g. 0,2,4 for L=4)
        
        l(:) = [(ll, ll=0, Lmax, 2)] ! (0, 2, 4, ..., L)
        pfJ  = 4*Pi * sum( [( (2*l(ii)+1) * Sl(nlm,l(ii)), ii=1, size(l))] )
    end
    
    ! SH functions for m=0
    include "include/Yl0-functions.f90"
    
end module idealstate
