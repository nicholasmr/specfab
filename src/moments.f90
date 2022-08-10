! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2022

! Calculations of normalized ODF moments (structure tensors) given the spectral expansion coefficieints of an ODF.

module moments 

    use mandel

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real, parameter, private    :: Pi = 3.141592653589793
    integer, parameter, private :: identity(3,3)  = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    integer, parameter          :: I_l0=1, I_l2=I_l0+1, I_l4=I_l2+(2*2+1), I_l6=I_l4+(2*4+1), I_l8=I_l6+(2*6+1), I_l10=I_l8+(2*8+1) ! Indices for extracting l=0,2,4,6,8 coefs of nlm
    integer, parameter          :: I_l12=I_l10+(2*10+1), I_l14=I_l12+(2*12+1), I_l16=I_l14+(2*14+1), I_l18=I_l16+(2*16+1), I_l20=I_l18+(2*18+1)

contains      

    !---------------------------------
    ! Conversion between nlm and <c^k> 
    !---------------------------------
              
    function a2(nlm) 
        ! a^(2) := <c^2> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a2(3,3)
        complex(kind=dp)             :: n2m(-2:2) 
        n2m = nlm(I_l2:(I_l4-1))
        a2 = f_ev_c2(nlm(1),n2m)
    end

    function a4(nlm) 
        ! a^(4) := <c^4> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a4(3,3,3,3)
        complex(kind=dp)             :: n2m(-2:2), n4m(-4:4)
        n2m = nlm(I_l2:(I_l4-1))
        n4m = nlm(I_l4:(I_l6-1))
        a4 = f_ev_c4(nlm(1),n2m,n4m)
    end
    
    function a2_to_nlm(a2) result(nlm)
        ! Get n_2^m from a^(2)
        implicit none
        real(kind=dp), intent(in) :: a2(3,3)
        complex(kind=dp)          :: nlm(1+5)
        nlm = (0.0, 0.0) ! init
        include "include/a2_to_nlm__body.f90"
    end

    function a4_to_nlm(a4) result(nlm)
        ! Get n_2^m and n_4^m from a^(4)
        implicit none
        real(kind=dp), intent(in) :: a4(3,3,3,3)
        complex(kind=dp)          :: nlm(1+5+9)
        real(kind=dp)             :: a4Mandel(6,6)        
        nlm = 0.0 ! init
        a4Mandel = a4_to_mat(a4) ! 6x6 Mandel matrix of a4
        include "include/a4_to_nlm__body.f90"
    end

    !---------------------------------
    ! All moments <c^k> (k = 2,4,6,8)
    !---------------------------------
    
    subroutine f_ev_ck(nlm, opt, ev_c2,ev_c4,ev_c6,ev_c8)
        
        ! "ev_ck" are the structure tensors <c^k> := a^(k) for a given psi(theta,phi) (grain orientation distribution)
        !    prescribed in terms of the array "nlm" of spectral expansion coefficients.
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        character*1, intent(in)      :: opt 
        real(kind=dp), intent(inout) :: ev_c2(3,3),ev_c4(3,3,3,3),ev_c6(3,3,3,3, 3,3),ev_c8(3,3,3,3, 3,3,3,3)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
        
        n00 = nlm(1)
        n2m = nlm(I_l2:(I_l4-1))
        n4m = nlm(I_l4:(I_l6-1))
        ev_c2 = f_ev_c2(n00,n2m)
        ev_c4 = f_ev_c4(n00,n2m,n4m)
        
        if (opt == 'f') then
            ! Full calculation?
            n6m = nlm(I_l6:(I_l8-1))
            n8m = nlm(I_l8:(I_l10-1))
            ev_c6 = f_ev_c6(n00,n2m,n4m,n6m)
            ev_c8 = f_ev_c8(n00,n2m,n4m,n6m,n8m)
        else
            ! Reduced calculation? Higher-order contributions not needed by caller.
            ev_c6 = 0.0d0
            ev_c8 = 0.0d0
        end if
    end
    
    subroutine f_ev_ck_Mandel(nlm, ev_c2_Mandel, ev_c4_Mandel)
            
        ! a2 and a4 in Mandel's notation. 
        ! Mandel's notation should be used whereever possible for performance.
            
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(inout) :: ev_c2_Mandel(6),ev_c4_Mandel(6,6)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4)
        
        n00 = nlm(1)
        n2m = nlm(I_l2:(I_l4-1))
        n4m = nlm(I_l4:(I_l6-1))
        ev_c2_Mandel = mat_to_vec(f_ev_c2(n00,n2m)) ! 6x1 Mandel vector of a2
        ev_c4_Mandel = f_ev_c4_Mandel(n00,n2m,n4m)  ! 6x6 Mandel matrix of a4
    end

    !---------------------------------
    ! Indidual <c^k> moments
    !---------------------------------

    function f_ev_c0(n00) result(ev)
        ! Integral over orientation distribution = total number of c-axes (grains) = N
        complex(kind=dp), intent(in) :: n00
        real(kind=dp) :: ev
        ev = REAL(sqrt(4*Pi)*n00)
    end

    function f_ev_c2(n00,n2m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2)
        complex(kind=dp) :: n2mhat(0:2)
!        real(kind=dp) :: k = 0.0
        real(kind=dp) :: ev(3,3)
        ev = 0.0
        include "include/ev_c2__body.f90"
!        ev = ev * k/f_ev_c0(n00) ! not needed, <c^2> already normalized
    end
    
    function f_ev_c4(n00,n2m,n4m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp) :: k = 0.0, ev(3,3, 3,3)
        ev = 0.0
        include "include/ev_c4__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c4_Mandel(n00,n2m,n4m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp) :: k = 0.0, ev(6,6)
        include "include/ev_c4_Mandel__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c6(n00,n2m,n4m,n6m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        real(kind=dp) :: k = 0.0, ev(3,3, 3,3, 3,3)
        ev = 0.0
        include "include/ev_c6__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c8(n00,n2m,n4m,n6m,n8m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
        real(kind=dp) :: k = 0.0, ev(3,3, 3,3, 3,3, 3,3)
        ev = 0.0
        include "include/ev_c8__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    !---------------------------------
    ! Fabric eigen frames
    !---------------------------------

    subroutine frame(nlm, ftype, e1,e2,e3, eigvals)

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        character*1, intent(in)      :: ftype ! 'x','e','p' (cartensian frame, a2 eigen frame, 45deg-rotated eigen frame)
        integer, parameter           :: n = 3
        real(kind=dp), intent(out)   :: e1(n),e2(n),e3(n), eigvals(3)
        real(kind=dp)                :: p(n),q(n)
        ! If eigen frame
        integer            :: inf
        integer, parameter :: l=3*3-1
        real(kind=dp)      :: e_ij(n,n), work(l)
        
        eigvals = 0.0
        
        ! Cartesian frame
        if (ftype == 'x') then
            e1 = [1,0,0] 
            e2 = [0,1,0] 
            e3 = [0,0,1] 
            
        else
            ! a2 eigen frame    
            e_ij = a2(nlm)
            call dsyev('V','U',n,e_ij,n,eigvals,work,l,inf)
            e1 = e_ij(:,3)
            e2 = e_ij(:,2)
            e3 = e_ij(:,1)
            eigvals = [eigvals(3),eigvals(2),eigvals(1)] ! Largest first
            
            ! Rotated a2 eigen frame
            if (ftype == 'p') then
                p = (e1+e2)/sqrt(2.0) 
                q = (e1-e2)/sqrt(2.0) 
                e1 = p
                e2 = q
                ! cross product
                e3(1) = p(2) * q(3) - p(3) * q(2)
                e3(2) = p(3) * q(1) - p(1) * q(3)
                e3(3) = p(1) * q(2) - p(2) * q(1)
            end if     
               
        end if
    end

    subroutine a4_eigentensors(nlm, Q1,Q2,Q3,Q4,Q5,Q6, eigvals6)  

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: M(6,6)
        real(kind=dp), intent(out)   :: Q1(3,3),Q2(3,3),Q3(3,3),Q4(3,3),Q5(3,3),Q6(3,3)
        
        integer, parameter :: n = 6
        real(kind=dp), intent(out) :: eigvals6(n)
        integer            :: inf
        integer, parameter :: l6=n*n-1
        real(kind=dp)      :: work6(l6)
        
        M = a4_to_mat(a4(nlm)) ! 6x6
        call dsyev('V','U',n,M,n,eigvals6,work6,l6,inf)
        Q1 = vec_to_mat(M(:,1)) 
        Q2 = vec_to_mat(M(:,2)) 
        Q3 = vec_to_mat(M(:,3)) 
        Q4 = vec_to_mat(M(:,4))
        Q5 = vec_to_mat(M(:,5))
        Q6 = vec_to_mat(M(:,6))
    end

end module moments
