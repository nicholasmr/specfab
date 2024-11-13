! N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

! Idealized CPO states

module idealstate 

    use header
    use rotation
    use frames

    implicit none 
    
contains

    function nlm_ideal(m, colat, L) result(nlm_rot)
    
        implicit none
        
        real(kind=dp), intent(in) :: m(3), colat ! symmetry axis, circle colat w.r.t. m
        integer, intent(in)       :: L
        complex(kind=dp)          :: nlm(nlm_lenvec(L)), nlm_rot(nlm_lenvec(L))
        integer, parameter        :: N = 6
        real(kind=dp)             :: Yl0_list(0:N), theta, phi
        
        nlm(:) = 0.0d0
        Yl0_list = [Y00(colat), Y20(colat), Y40(colat), Y60(colat), Y80(colat), Y100(colat), Y120(colat)]
        
        do ii=0,min(N,int(L/2))
!            print *, ii, ILm(ii*2)
            nlm(ILm(ii*2)) = Yl0_list(ii)
        end do

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
        real(kind=dp)                :: e1(3),e2(3),e3(3), a2_eigvals(3)
        real(kind=dp)                :: Q1(3,3),Q2(3,3),Q3(3,3),Q4(3,3),Q5(3,3),Q6(3,3), a4_eigvals(6)

        do ii=1,size(nhat20)
            nlm(:) = 0.0d0
            nlm(1) = 1.0d0
            nlm(ILm(2)) = nhat20(ii)
            nlm(ILm(4)) = nhat40(ii)
            
            call frame(nlm, 'e', e1,e2,e3, a2_eigvals)
            isvalid_a2 = (maxval(a2_eigvals) < 1.0d0) .and. (minval(a2_eigvals) > 0.0d0) ! a2 eigenvalues are within bounds
            
            call a4_eigentensors(nlm, Q1,Q2,Q3,Q4,Q5,Q6, a4_eigvals)  
            idvalid_a4 = (maxval(a4_eigvals) < 1.0d0) .and. (minval(a4_eigvals) > 0.0d0) ! a4 eigenvalues are within bounds
            
            isvalid(ii) = isvalid_a2 .and. idvalid_a4
        end do                
    end 
    
    ! SH functions for m=0
    
    function Y00(theta) result (Y)
        implicit none
        real(kind=dp), intent(in) :: theta ! latitude
        real(kind=dp)             :: Y
        Y = 1/2.0d0 * sqrt(1/Pi)
    end
    
    function Y20(theta) result (Y)
        implicit none
        real(kind=dp), intent(in) :: theta ! latitude
        real(kind=dp)             :: Y
        Y = 1/4.0d0 * sqrt(5/Pi) * (3*cos(theta)**2 - 1)
    end
    
    function Y40(theta) result (Y)
        implicit none
        real(kind=dp), intent(in) :: theta ! latitude
        real(kind=dp)             :: Y
        Y = 3/16.0d0 * sqrt(1/Pi) * (35*cos(theta)**4 - 30*cos(theta)**2 + 3)
    end
    
    function Y60(theta) result (Y)
        implicit none
        real(kind=dp), intent(in) :: theta ! latitude
        real(kind=dp)             :: Y
        Y = 1/32.0d0 * sqrt(13/Pi) * (231*cos(theta)**6 - 315*cos(theta)**4 + 105*cos(theta)**2 - 5)
    end

    function Y80(theta) result (Y)
        implicit none
        real(kind=dp), intent(in) :: theta ! latitude
        real(kind=dp)             :: Y
        Y = 1/256.0d0 * sqrt(17/Pi) * (6435*cos(theta)**8 - 12012*cos(theta)**6 + 6930*cos(theta)**4 - 1260*cos(theta)**2 + 35)
    end
    
    function Y100(theta) result (Y)
        implicit none
        real(kind=dp), intent(in) :: theta ! latitude
        real(kind=dp)             :: Y
        Y = 1/512.0d0 * sqrt(21/Pi) * (46189*cos(theta)**10 - 109395*cos(theta)**8 + 90090*cos(theta)**6 - 30030*cos(theta)**4 + 3465*cos(theta)**2 - 63)
    end
       
    function Y120(theta) result (Y)
        implicit none
        real(kind=dp), intent(in) :: theta ! latitude
        real(kind=dp)             :: Y
        Y = 1/2048.0d0 * sqrt(25/Pi) * (676039*cos(theta)**12 - 1939938*cos(theta)**10 + 2078505*cos(theta)**8 - 1021020*cos(theta)**6 + 225225*cos(theta)**4 - 18018*cos(theta)**2 + 231)
    end

!    function Y140(theta) result (Y)
!        implicit none
!        real(kind=dp), intent(in) :: theta ! latitude
!        real(kind=dp)             :: Y
!        Y = 1/4096.0d0 * sqrt(29/Pi) * (5014575*cos(theta)**14 - 16900975*cos(theta)**12 + 22309287*cos(theta)**10 - 14549535*cos(theta)**8 + 4849845*cos(theta)**6 - 765765*cos(theta)**4 + 45045*cos(theta)**2 - 429)
!    end
        
end module idealstate
