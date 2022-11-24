! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2022

! For converting between full (nlm) and reduced (rnlm) forms of nlm by noting that real-valued ODFs fulfill: n_l^{-m} = (-1)^m conj(n_l^m) 

module reducedform

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
!    integer, private :: ii, jj, kk, ll, mm ! Loop indices
    integer, private :: jj, kk, ll, mm ! Loop indices (ii not used)
    integer, parameter, private :: Lcap__max = 30 ! Hard limit (should be the same value as in dynamics.f90)
    
    !----------
    
    ! Set by init()
    integer, private :: Lcap ! Truncation "L" of expansion series 
    integer, private :: nlm_len
    integer          :: rnlm_len 
    integer, private :: rnlm_pos_len 
    integer, private :: rnlm_neg_len 

    ! Max lengths of nlm and rnlm (all,pos,neg) arrays
    integer, parameter, private :: nlm_lenvec(0:Lcap__max)      = [( (ll+1)*(ll+2)/2, ll=0, Lcap__max, 1)] ! nlm length for a given Lcap: sum_{l=0}^L (2*l+1) **for even l** = sum_{l=0}^{L/2} (4*l+1) = (L+1)*(L+2)/2
    integer, parameter          :: rnlm_all_lenvec(0:Lcap__max) = [( (ll+2)**2/4,     ll=0, Lcap__max, 1)] ! nlm length for a given Lcap: sum_{l=0}^L l+1     **for even l** = sum_{l=0}^{L/2} (2*l+1) = (L+2)^2/4
    integer, parameter, private :: rnlm_pos_lenvec(0:Lcap__max) = [( ll*(ll+2)/4,     ll=0, Lcap__max, 1)] ! nlm length for a given Lcap: sum_{l=0}^L l       **for even l** = sum_{l=0}^{L/2} (2*l)   = L*(L+2)/4
    integer, parameter, private :: rnlm_neg_lenvec(0:Lcap__max) = rnlm_pos_lenvec
    
    integer, parameter, private :: nlm_lenmax      = nlm_lenvec(Lcap__max)
    integer, parameter          :: rnlm_lenmax     = rnlm_all_lenvec(Lcap__max)
    integer, parameter, private :: rnlm_pos_lenmax = rnlm_pos_lenvec(Lcap__max)
    integer, parameter, private :: rnlm_neg_lenmax = rnlm_neg_lenvec(Lcap__max)
    
    ! Full --> reduced
    integer, parameter, private :: I_all(rnlm_lenmax)     = [( (kk+nlm_lenvec(ll)-ll,   kk=0, ll,   1), ll=0, Lcap__max, 2)] ! rnlm     = nlm(I_all(1:rnlm_all_lenvec(L))) = (n_0^0, n_2^0, n_2^1, n_2^2, n_4^0, ... ) 
    integer, parameter, private :: I_pos(rnlm_pos_lenmax) = [( (kk+nlm_lenvec(ll)-ll,   kk=1, ll,   1), ll=0, Lcap__max, 2)] ! rnlm_pos = nlm(I_pos(1:rnlm_pos_lenvec(L))) = (n_2^1, n_2^2, n_4^1, n_4^2, n_4^3, n_4^4, ... ) 
    integer, parameter, private :: I_neg(rnlm_neg_lenmax) = [( (kk+nlm_lenvec(ll)-2*ll, kk=0, ll-1, 1), ll=0, Lcap__max, 2)] ! rnlm_neg = nlm(I_neg(1:rnlm_neg_lenvec(L))) = (n_2^-2, n_2^-1, n_4^-4, n_4^-3, n_4^-2, n_4^-1, n_6^-6, ... ) 

    ! Reduced --> full
    integer, parameter, private :: I_full(nlm_lenmax) = [( (kk+rnlm_all_lenvec(ll)-1*ll, kk=ll, 1, -1), (kk+rnlm_all_lenvec(ll)-1*ll, kk=0, ll, 1), ll=0, Lcap__max, 2)] ! nlm(1:nlm_len) = Mdiag(1:nlm_len) * rnlm(I_full(1:nlm_len(L))) ...and then nlm(I_neg(1:rnlm_neg_len(L))) = conj( nlm(I_neg(1:rnlm_neg_len(L))) )
    integer, parameter, private :: Mdiag(nlm_lenmax)  = [( ((-1)**(mm), mm=ll, 1, -1), (1, kk=0, ll, 1), ll=0, Lcap__max, 2)] ! For adjusting odd negative m values with a factor of -1 (i.e. the (-1)^m factor in n_l^{-m} = (-1)^m conj(n_l^m))

!    ! rnlm indices
!    integer, parameter :: rI_pos(rnlm_pos_lenmax) = [( (kk+rnlm_all_lenvec(ll), kk=-ll+1, 0, 1), ll=0, Lcap__max, 2)] ! rnlm(rI_pos(1:rnlm_pos_len)) = (n_2^1, n_2^2, n_4^1, n_4^2, n_4^3, n_4^4, ...)

    ! (l,m) list
!    integer, parameter, private :: lm(2,nlm_lenmax)   = reshape([( (ll,mm, mm=-ll,ll), ll=0,  Lcap__max,2)], [2,nlm_lenmax])  ! These are the (l,m) values corresponding to the coefficients in "nlm". (copied from specfab.f90). Not used, so outcommented to avoid compiler warnings.
    integer, parameter          :: rlm(2,rnlm_lenmax) = reshape([( (ll,mm, mm=0,ll),   ll=0,  Lcap__max,2)], [2,rnlm_lenmax]) ! These are the (l,m) values corresponding to the coefficients in rnlm.

contains      

    !---------------------------------
    ! Init
    !---------------------------------

    subroutine initreduced(Lcap_) 

        ! Needs to be called once before using the module routines.

        implicit none    
        integer, intent(in) :: Lcap_ ! Truncation "Lcap"
        
        Lcap = Lcap_ ! Save internal copy
        
        nlm_len      = nlm_lenvec(Lcap)
        rnlm_len     = rnlm_all_lenvec(Lcap)
        rnlm_pos_len = rnlm_pos_lenvec(Lcap)
        rnlm_neg_len = rnlm_neg_lenvec(Lcap)
        
    end

    !---------------------------------
    ! Convert fabric dynamics operators to work on reduced form of nlm 
    !---------------------------------

    subroutine reduce_M(M, Mrr, Mri, Mir, Mii)
        
        ! *** Reduce full fabric matrix M to Mrr, Mri, Mii, Mir such that: 
        ! 
        !   d(rnlm)/dt =                  matmul(Mrr,real(rnlm)) + matmul(Mri,aimag(rnlm))  
        !                + complex(0,1)*( matmul(Mir,real(rnlm)) + matmul(Mii,aimag(rnlm)) ) 
        ! 
        ! where rnlm is the reduced state vector
        
        implicit none

        complex(kind=dp), intent(in)                             :: M(nlm_len,nlm_len)
        real(kind=dp), intent(out), dimension(rnlm_len,rnlm_len) :: Mrr, Mri, Mir, Mii 
        complex(kind=dp)                                         :: Z(rnlm_len,2)
        real(kind=dp), dimension(rnlm_len)                       :: Vp,Vn, Qp,Qn
        integer                                                  :: s, m_jj, jj_p,jj_n
        
        Mrr = 0.0d0
        Mri = 0.0d0
        Mir = 0.0d0
        Mii = 0.0d0

        do jj = 1, rnlm_len

            m_jj = rlm(2,jj) ! +m value               
            jj_p = I_all(jj)         ! +m entry in full array
            jj_n = jj_p - (2*m_jj+0) ! -m entry in full array

            Z = M(I_all(1:rnlm_len), [jj_p,jj_n])
            Vp = real( Z(:,1))
            Qp = aimag(Z(:,1))
            Vn = real( Z(:,2))
            Qn = aimag(Z(:,2))
            
            s = -1 ! sign if m is odd
            if (MOD(m_jj,2) .eq. 0) s = +1 ! sign if m is even
            if (jj_p .eq. jj_n) s = 0 ! prevent counting m=0 modes twice

            Mrr(:,jj) = Mrr(:,jj) + Vp + s*Vn
            Mri(:,jj) = Mri(:,jj) - Qp + s*Qn
            Mir(:,jj) = Mir(:,jj) + Qp + s*Qn
            Mii(:,jj) = Mii(:,jj) + Vp - s*Vn
        end do

    end

!   *** DEPRECATED, USE reduce_M() ***
!    function dndt_to_drndt(dndt) result (drndt)
!        
!        ! *** Convert dndt_ij (full matrix) to drndt_ij (reduced matrix) such that: d/dt (rnlm)_i = drndt_ij (rnlm)_j

!        ! NOTICE: x-y rotation will make rnlm_j complex valued, in which case drndt_ij does not work as intended (i.e. the full problem involving dndt_ij and nlm_j must be solved).

!        implicit none

!        complex(kind=dp), intent(in) :: dndt(nlm_len,nlm_len)
!        real(kind=dp)                :: drndt(rnlm_len,rnlm_len)
!        integer                      :: kk

!        drndt = 0.0d0

!        do ii = 1, rnlm_len
!            do jj = 1, nlm_len
!                
!                kk = I_full(jj)
!                
!                if (lm(2,jj) >= 0) then ! m >= 0
!                    drndt(ii,kk) = drndt(ii,kk) + real(dndt(I_all(ii), jj))
!                else ! m < 0
!                    if (MOD(lm(2,jj),2) .eq. 0) then ! m even
!                        drndt(ii,kk) = drndt(ii,kk) + real(dndt(I_all(ii), jj))
!                    else ! m odd
!                        drndt(ii,kk) = drndt(ii,kk) - real(dndt(I_all(ii), jj))
!                    end if
!                end if
!                
!            end do
!        end do
!    end

    !---------------------------------
    ! Convert between reduced and full forms
    !---------------------------------

    function rnlm_to_nlm(rnlm) result(nlm)

        ! *** Convert reduced (rnlm) to full (nlm)

        implicit none

        complex(kind=dp), intent(in) :: rnlm(rnlm_len)
        complex(kind=dp)             :: nlm(nlm_len)

        nlm(1:nlm_len) = Mdiag(1:nlm_len) * rnlm(I_full(1:nlm_len)) 
        nlm(I_neg(1:rnlm_neg_len)) = conjg( nlm(I_neg(1:rnlm_neg_len)) )
    end

    function nlm_to_rnlm(nlm) result(rnlm)

        ! *** Convert full (nlm) to reduced (rnlm) 

        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        complex(kind=dp)             :: rnlm(rnlm_len)
        
        rnlm = nlm(I_all(1:rnlm_len))
    end

    !---------------------------------
    ! For validation purposes only
    !---------------------------------

    function nlm_to_rnlm_pos(nlm) result(rnlm_pos)

        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        complex(kind=dp)             :: rnlm_pos(rnlm_pos_len)
        
        rnlm_pos = nlm(I_pos(1:rnlm_pos_len))
    end

    function nlm_to_rnlm_neg(nlm) result(rnlm_neg)

        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        complex(kind=dp)             :: rnlm_neg(rnlm_neg_len)
        
        rnlm_neg = nlm(I_neg(1:rnlm_neg_len))
    end

end module reducedform 
