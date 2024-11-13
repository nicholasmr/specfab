! N. M. Rathmann <rathmann@nbi.ku.dk>, 2021

program demo

    use specfab
    
    implicit none

!    integer, parameter :: dp = 8
!    integer            :: ii ! loop counter
    
    integer, parameter :: L = 2
    
    ! don't use names nlm_len or rnlm_len (used by specfab)   
    integer, parameter :: nlmlen  = (L+1)*(L+2)/2 
    integer, parameter :: rnlmlen = (L+2)**2/4    
    
    complex(kind=dp)   :: nlm0(nlmlen)   = 9999.0d0
    complex(kind=dp)   :: nlm1(nlmlen)   = 9999.0d0
    complex(kind=dp)   :: nlm2(nlmlen)   = 9999.0d0
    complex(kind=dp)   :: rnlm0(rnlmlen) = 9999.0d0
    complex(kind=dp)   :: rnlm1(rnlmlen) = 9999.0d0

!    complex(kind=dp)   :: nlm_re_cmp(nlmlen) = 9999.0d0

    real(kind=dp)    :: eps(3,3), omg(3,3), ugrad(3,3)
    complex(kind=dp) :: M(nlmlen,nlmlen) !, drndt(rnlmlen,rnlmlen)
    real(kind=dp), dimension(rnlmlen,rnlmlen) :: Mrr,Mii,Mri,Mir

    call initreduced(L) 
    
    print *, 'nlm_len, nlm_reduced_len = ', nlmlen, rnlmlen

    ! Fill nlm with some random values
    do ii = 1, nlmlen ! loop over unique values  
        ! insert some arbitrary values
        if (lm(2,ii) .eq. 0) then
            nlm0(ii) = complex(ii, 0) ! m=0 components are always real
        else
            nlm0(ii) = complex(ii, 10*ii) 
        end if
    end do            

    ! *** Test 1: extract the parts _all, _pos, _neg 

    print *, '------------------------------'    
    print *, '*** Test extraction of rnlm_pos, rnlm_pos and rnlm_neg'
    print *, '------------------------------'  
    print *, 'nlm0     ', nlm0
    print *, 'rnlm0_all', nlm_to_rnlm(nlm0)
    print *, 'rnlm0_pos', nlm_to_rnlm_pos(nlm0)
    print *, 'rnlm0_neg', nlm_to_rnlm_neg(nlm0)
    print *, ' '
    
    ! *** Test 2: Conversion between nlm (full) <--> rnlm (reduced)
    
    rnlm0 = nlm_to_rnlm(nlm0)
    nlm1  = rnlm_to_nlm(rnlm0)
    rnlm1 = nlm_to_rnlm(nlm1)
    nlm2  = rnlm_to_nlm(rnlm1)

    print *, '------------------------------'    
    print *, '*** Test conversion: real parts'
    print *, '------------------------------'    
    print *, 'real nlm0  = ', real(nlm0)
    print *, 'real rnlm0 = ', real(rnlm0)
    print *, 'real nlm1  = ', real(nlm1)
    print *, 'real rnlm1 = ', real(rnlm1)
    print *, 'real nlm2  = ', real(nlm2)
    print *, '------------------------------'    
    print *, '*** Test conversion: imaginary parts'
    print *, '------------------------------'    
    print *, 'imag nlm0  = ', aimag(nlm0)
    print *, 'imag rnlm0 = ', aimag(rnlm0)
    print *, 'imag nlm1  = ', aimag(nlm1)
    print *, 'imag rnlm1 = ', aimag(rnlm1)
    print *, 'imag nlm2  = ', aimag(nlm2)
    print *, '------------------------------'    
    print *, '... last line (nlm2) should equal the third line (nlm1) for both real and imaginary parts. If so, the reduced (rnlm) <--> full (nlm) conversion works.'
    print *, ' '

    ! *** Test 3: M --> M* decomposition
    ! Test using lattice rotation
    
    print *, '------------------------------'    
    print *, '*** Test M --> (Mrr,Mri,Mir,Mii) decomposition'
    print *, '------------------------------'  
    
    ! Some random strain-rate tensor
    ugrad = reshape([0.75,0.,0., 0.,0.25,0., 0.,0.,-1.], [3, 3]) ! vertical compression
    ugrad = ugrad + reshape([0.,0.,-1., 0.,0.,0., +1.,0.,0.], [3, 3]) ! x-z shear
    ugrad = ugrad + reshape([0.,-1.,0., +1.,0.,0., 0.,0.,0.], [3, 3]) ! x-y shear
    ugrad = ugrad + reshape([0.,0.,0., 0.,0.,-1., 0.,+1.,0.], [3, 3]) ! y-z shear
    eps = (ugrad+transpose(ugrad))/2 ! strain-rate
    omg = (ugrad-transpose(ugrad))/2 ! spin
    
    call initspecfab(L) 
    M  = M_LROT(eps, omg, 1d0, 0d0) 
    call reduce_M(M, Mrr,Mri,Mir,Mii) ! Any x-y rotation will make nlm (rnlm) complex valued, all four Mrr,Mri,Mir,Mii must be considere (else Mrr suffices)
   
    nlm0  = nlm2 ! Test fabric. Must be has correctly conjugate components for -m modes, hence we use nlm2
    rnlm0 = nlm_to_rnlm(nlm0) ! reduced form
    print *, 'nlm0 and rnlm0=nlm_to_rnlm(nlm0) used for test are:'
    print *, nlm0
    print *, rnlm0
    print *, 'Testing... '
    nlm1  = matmul(M, nlm0) ! true d(nlm)/dt
    rnlm1 = matmul(Mrr,real(rnlm0)) + matmul(Mri,aimag(rnlm0))  + complex(0,1)*( matmul(Mir,real(rnlm0)) + matmul(Mii,aimag(rnlm0)) ) ! d(rnlm)/dt
    print *, '             matmul(M, nlm0)                  = ', nlm1
    print *, 'rnlm_to_nlm( matmul(Mrr, real(rnlm0)) + ... ) = ', rnlm_to_nlm(rnlm1)
    print *, 'Difference  = ', nlm1 - rnlm_to_nlm(rnlm1)
    
end program

