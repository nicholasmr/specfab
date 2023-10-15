! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2023

! Normalized vector moments (structure tensors) given the spectral expansion coefficients of CPO.

module moments 

    use header
    use mandel
    use tensorproducts

    implicit none 
    
    real(kind=dp), parameter :: a2_iso(3,3)     = identity/3
    real(kind=dp), parameter :: a4_iso(3,3,3,3) = reshape([ ( ( ( ( (identity(ii,jj)*identity(kk,ll)+identity(ii,kk)*identity(jj,ll)+identity(ii,ll)*identity(jj,kk))/15, ii=1,3), jj=1,3), kk=1,3), ll=1,3) ], [3,3,3,3]) 
    
contains      

    !---------------------------------
    ! Conversion between nlm and <c^k> 
    !---------------------------------
              
    function a2(nlm) 
        ! a^(2) := <c^2> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a2(3,3)

        a2 = f_ev_c2(nlm(L0rng),nlm(L2rng))
    end

    function a4(nlm) 
        ! a^(4) := <c^4> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a4(3,3,3,3)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)

        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        a4 = f_ev_c4(n00,n2m,n4m)
    end
    
    function a6(nlm) 
        ! a^(6) := <c^6> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a6(3,3,3,3,3,3)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        a6 = f_ev_c6(n00,n2m,n4m,n6m)
    end
    
    function a2_to_nlm(a2) result(nlm)
        ! Get n_2^m from a^(2)
        implicit none
        real(kind=dp), intent(in) :: a2(3,3)
        complex(kind=dp)          :: nlm(nlm_lenvec(2)) 
        nlm = (0.0d0, 0.0d0) ! init
        include "include/a2_to_nlm__body.f90"
    end

    function a4_to_nlm(a4) result(nlm)
        ! Get n_2^m and n_4^m from a^(4)
        implicit none
        real(kind=dp), intent(in) :: a4(3,3,3,3)
        complex(kind=dp)          :: nlm(nlm_lenvec(4)) 
        real(kind=dp)             :: a4Mandel(6,6)        
        nlm = 0.0d0 ! init
        a4Mandel = a4_to_mat(a4) ! 6x6 Mandel matrix of a4
        include "include/a4_to_nlm__body.f90"
    end

    function a6_to_nlm(a6) result(nlm)
        ! Get n_2^m and n_4^m from a^(4)
        implicit none
        real(kind=dp), intent(in) :: a6(3,3,3,3,3,3)
        complex(kind=dp)          :: nlm(nlm_lenvec(6)) 
        nlm = 0.0d0 ! init
        include "include/a6_to_nlm__body.f90"
    end

    !---------------------------------
    ! All moments <c^k> (k = 2,4,6,8)
    !---------------------------------
    
    subroutine f_ev_ck(nlm, opt, ev_c2,ev_c4,ev_c6,ev_c8)
        
        ! "ev_ck" are the structure tensors <c^k> := a^(k) for a given n(theta,phi) prescribed in terms of "nlm" expansion coefficients.
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        character*1, intent(in)      :: opt 
        real(kind=dp), intent(inout) :: ev_c2(3,3),ev_c4(3,3,3,3),ev_c6(3,3,3,3, 3,3),ev_c8(3,3,3,3, 3,3,3,3)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
        
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        ev_c2 = f_ev_c2(n00,n2m)
        ev_c4 = f_ev_c4(n00,n2m,n4m)
        
        if (opt == 'f') then
            ! Full calculation?
            n8m = nlm(L8rng)
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
        ! Mandel's notation should be used wherever possible for performance.
            
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(out)   :: ev_c2_Mandel(6),ev_c4_Mandel(6,6)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        ev_c2_Mandel = mat_to_vec(f_ev_c2(n00,n2m)) ! 6x1 Mandel vector of a2
        ev_c4_Mandel = f_ev_c4_Mandel(n00,n2m,n4m)  ! 6x6 Mandel matrix of a4
    end

    !---------------------------------
    ! Individual <c^k> moments
    !---------------------------------

    function f_ev_c0(n00) result(ev)
        ! Integral over orientation distribution = total number of c-axes = N
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
        ev = 0.0d0
        include "include/ev_c2__body.f90"
!        ev = ev * k/f_ev_c0(n00) ! not needed, <c^2> already normalized
    end
    
    function f_ev_c4(n00,n2m,n4m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp) :: k = 0.0d0, ev(3,3, 3,3)
        ev = 0.0d0
        include "include/ev_c4__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c4_Mandel(n00,n2m,n4m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp) :: k = 0.0d0, ev(6,6)
        ev = 0.0d0
        include "include/ev_c4_Mandel__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c6(n00,n2m,n4m,n6m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        real(kind=dp) :: k = 0.0d0, ev(3,3, 3,3, 3,3)
        ev = 0.0d0
        include "include/ev_c6__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c8(n00,n2m,n4m,n6m,n8m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
        real(kind=dp) :: k = 0.0d0, ev(3,3, 3,3, 3,3, 3,3)
        ev = 0.0d0
        include "include/ev_c8__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end

    !---------------------------------
    ! Orthogonal distribution
    !---------------------------------
    
    function a2_orth(blm,nlm) result(ev)
    
        ! Calculates <v^2> where v = b x n

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)
        
        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_v2__body.f90"
        ev = ev * k/norm
    end
    
    function a4_orth(blm,nlm) result(ev)
    
        ! Calculates <v^4> where v = b x n

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3,3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)

        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_v4__body.f90"
        ev = ev * k/norm 
    end
  
    function a4_joint(blm,nlm) result(ev)
    
        ! Calculates < b^2 n^2 xi(r,r') > where kernel xi=sin^2(alpha) where alpha = angle between n and b vectors

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3,3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)

        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_c2b2__body.f90"
        ev = ev * k/norm 
    end
    
    function a4_jointcross(blm,nlm) result(ev)
    
        ! Calculates <n^2 v^2> where v = b x n

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3,3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)

        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_c2v2__body.f90"
        ev = ev * k/norm 
    end
    
    !---------------------------------
    ! Fabric eigen frames
    !---------------------------------

    subroutine frame(nlm, ftype, e1,e2,e3, eigvals)

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        character*1, intent(in)      :: ftype ! 'x','e','p' (Cartesian frame, a2 eigenframe, 45deg rotated eigenframe)
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
    
    !---------------------------------
    ! AUX
    !---------------------------------
    
    subroutine decompose_nlm(nlm, n00,n2m,n4m,n6m)

        implicit none
        
        complex(kind=dp), intent(in)  :: nlm(:)
        complex(kind=dp), intent(out) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)

        n4m = 0.0d0
        n6m = 0.0d0
        
        n00 = nlm(L0rng)
        n2m = nlm(L2rng)
        if (size(nlm) >= nlm_lenvec(4)) n4m = nlm(L4rng)
        if (size(nlm) >= nlm_lenvec(6)) n6m = nlm(L6rng)
    end
    
    subroutine ai_orthotropic(nlm_1,nlm_2,nlm_3, a2_i, a4_ii, a4_jk)
    
        implicit none
        
        complex(kind=dp), intent(in) :: nlm_1(:), nlm_2(:), nlm_3(:)
        real(kind=dp), intent(out)   :: a2_i(3, 3,3), a4_ii(3, 3,3,3,3), a4_jk(3, 3,3,3,3) ! a4_jk = <m_j m_j m_k m_k>

        a2_i(1, :,:)      = a2(nlm_1) ! <m1'^2>
        a2_i(2, :,:)      = a2(nlm_2) ! <m2'^2>
        a4_ii(1, :,:,:,:) = a4(nlm_1) ! <m1'^4>
        a4_ii(2, :,:,:,:) = a4(nlm_2) ! <m2'^4>
        a4_jk(3, :,:,:,:) = a4_joint(nlm_1,nlm_2) ! <m1'^2 m2'^2>

        ! Is nlm_3 given? Use it
        if (real(nlm_3(1)) > 1e-8) then
            a2_i(3, :,:)      = a2(nlm_3) ! <m3'^2>
            a4_ii(3, :,:,:,:) = a4(nlm_3) ! <m3'^4>
            a4_jk(1, :,:,:,:) = a4_joint(nlm_2,nlm_3) ! <m2'^2 m3'^2>
            a4_jk(2, :,:,:,:) = a4_joint(nlm_1,nlm_3) ! <m1'^2 m3'^2>
            
        ! else derive nlm_3 from nlm_1, nlm_2
        else
            a2_i(3, :,:)      = a2_orth(nlm_1, nlm_2) ! <m3'^2> 
            a4_ii(3, :,:,:,:) = a4_orth(nlm_1, nlm_2) ! <m3'^4>
            a4_jk(1, :,:,:,:) = a4_jointcross(nlm_2,nlm_1) ! <m2'^2 m3'^2> 
            a4_jk(2, :,:,:,:) = a4_jointcross(nlm_1,nlm_2) ! <m1'^2 m3'^2> 
        end if
    end
    
    subroutine ai_orthotropic_discrete(mi, a2_i, a4_ii, a4_jk)

        implicit none
        
        real(kind=dp), intent(in)  :: mi(:,:,:) ! (3,3,N) = (m'_i, xyz comp., grain no.) 
        real(kind=dp), intent(out) :: a2_i(3, 3,3), a4_ii(3, 3,3,3,3), a4_jk(3, 3,3,3,3)
        real(kind=dp) :: m1(3), m2(3), m3(3)
        integer       :: ii, nn, N

        a2_i  = 0.0d0
        a4_ii = 0.0d0
        a4_jk = 0.0d0
        
        N = size(mi, dim=3) ! number of grains
        
        do nn=1,N
        
            do ii=1,3
                m1 = mi(ii,:,nn) ! use m1 as dummy
                a2_i(ii, :,:)      = a2_i(ii, :,:)      + outerprod(m1,m1)/N
                a4_ii(ii, :,:,:,:) = a4_ii(ii, :,:,:,:) + outerprod4(m1,m1,m1,m1)/N
            end do 
            
            m1 = mi(1,:,nn)
            m2 = mi(2,:,nn)
            m3 = mi(3,:,nn)
            
            a4_jk(1, :,:,:,:) = a4_jk(1, :,:,:,:) + outerprod4(m2,m2, m3,m3)/N
            a4_jk(2, :,:,:,:) = a4_jk(2, :,:,:,:) + outerprod4(m1,m1, m3,m3)/N
            a4_jk(3, :,:,:,:) = a4_jk(3, :,:,:,:) + outerprod4(m1,m1, m2,m2)/N
        end do
    end
    
end module moments
