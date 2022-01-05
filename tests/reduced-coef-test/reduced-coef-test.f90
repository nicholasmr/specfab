! N. M. Rathmann <rathmann@nbi.ku.dk>, 2021

program demo

    use specfab
    
    implicit none

    integer, parameter :: dp = 8
    integer            :: ii ! loop counter
    
    integer, parameter :: L = 2
    
    ! don't use names nlm_len or rnlm_len (used by specfab)   
    integer, parameter :: nlmlen  = (L+1)*(L+2)/2 
    integer, parameter :: rnlmlen = (L+2)**2/4    
    
    complex(kind=dp)   :: nlm0(nlmlen)   = 9999.0d0
    complex(kind=dp)   :: nlm1(nlmlen)   = 9999.0d0
    complex(kind=dp)   :: nlm2(nlmlen)   = 9999.0d0
    complex(kind=dp)   :: rnlm0(rnlmlen) = 9999.0d0
    complex(kind=dp)   :: rnlm1(rnlmlen) = 9999.0d0

    complex(kind=dp)   :: nlm_re_cmp(nlmlen) = 9999.0d0

    real(kind=dp)    :: eps(3,3), omg(3,3), ugrad(3,3)
    complex(kind=dp) :: dndt(nlmlen,nlmlen), drndt(rnlmlen,rnlmlen)

    call initreduced(L) 
    
    print *, 'nlm_len, nlm_reduced_len = ', nlmlen, rnlmlen

    ! Fill nlm with some random values
    do ii = 1, nlmlen ! loop over unique values  
        nlm0(ii) = complex(ii, 10*ii) ! insert some arbitrary values
        nlm_re_cmp(ii) = complex(ii, 0) ! used to test drndt below
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

    ! *** Test 3: dndt_ij --> drndt_ij
    ! Test using lattice rotation
    
    print *, '------------------------------'    
    print *, '*** Test dndt_ij <--> drndt_ij (assumes the operand rnlm is real valued)'
    print *, '------------------------------'  
    
    ! Some random strain-rate tensor
    ugrad = reshape([1.,0.,0., 0.,0.,0., 0.,0.,-1.], [3, 3]) + reshape([0.,0.,-1., 0.,0.,0., 1.,0.,0.], [3, 3]) ! NOTE x-y rotation will make rnlm complex valued, in which case drndt does not work as intended (i.e. the full dndt must be used).
    eps = (ugrad+transpose(ugrad))/2 ! strain-rate
    omg = (ugrad-transpose(ugrad))/2 ! spin
    
    call initspecfab(L) 
    dndt  = dndt_ij_LATROT(eps,omg, 0*eps,0d0,0d0,0d0, 1d0) 
    drndt = dndt_to_drndt(dndt)
    
    nlm_re_cmp = rnlm_to_nlm(nlm_to_rnlm(nlm_re_cmp)) ! the random (real-valued) (this is needed to that -m components are equal to the (negative) +m components.
    print *, 'nlm_re_cmp and nlm_to_rnlm(nlm_re_cmp) used for test are (MUST be real valued):'
    print *, nlm_re_cmp
    print *, nlm_to_rnlm(nlm_re_cmp)
    print *, 'Testing... '
    print *, 'nlm_to_rnlm(matmul(dndt, nlm_re_cmp))  = ', nlm_to_rnlm(matmul(dndt, nlm_re_cmp))
    print *, 'matmul(drndt, nlm_to_rnlm(nlm_re_cmp)) = ', matmul(drndt, nlm_to_rnlm(nlm_re_cmp))
    
end program

