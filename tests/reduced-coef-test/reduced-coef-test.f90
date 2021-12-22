! N. M. Rathmann <rathmann@nbi.ku.dk>, 2021

program demo

    use specfab    
    
    implicit none

    integer, parameter :: dp = 8
    integer            :: ii ! loop counter
    integer, parameter :: L = 4
    complex(kind=dp)   :: nlm0((L+1)*(L+2)/2)  = 9999.0d0
    complex(kind=dp)   :: nlm1((L+1)*(L+2)/2)  = 9999.0d0
    complex(kind=dp)   :: nlm2((L+1)*(L+2)/2)  = 9999.0d0
    complex(kind=dp)   :: nlm_red0((L+2)**2/4) = 9999.0d0
    complex(kind=dp)   :: nlm_red1((L+2)**2/4) = 9999.0d0

    call initspecfab(L) 
    
    print *, 'nlm_len ?= nlm_reduced_len + nlm_reduced_neg_len :: each are: ', nlm_len, nlm_reduced_len, nlm_reduced_neg_len

    ! Fill nlm with some random values
    do ii = 1, nlm_len ! loop over unique values  
        nlm0(ii) = complex(ii, 10*ii) ! insert some arbritrary values
    end do            

    ! Reduced (unique) part of nlm
    nlm_red0 = nlm_reduced(nlm0)
    nlm1     = nlm_full(nlm_red0)
    nlm_red1 = nlm_reduced(nlm1)
    nlm2     = nlm_full(nlm_red1)

    print *, '------------------------------'    
    print *, 'Comparing real parts...'
    print *, '------------------------------'    
    print *, 'real nlm                        = ', real(nlm0)
    print *, 'real nlm_reduced(nlm)           = ', real(nlm_red0)
    print *, 'real nlm_full(nlm_reduced(nlm)) = ', real(nlm1)
    print *, 'real ...reduced again           = ', real(nlm_red1)
    print *, 'real ...full again              = ', real(nlm2)
    print *, '------------------------------'    
    print *, 'Comparing imaginary parts...'
    print *, '------------------------------'    
    print *, 'imag nlm                        = ', aimag(nlm0)
    print *, 'imag nlm_reduced(nlm)           = ', aimag(nlm_red0)
    print *, 'imag nlm_full(nlm_reduced(nlm)) = ', aimag(nlm1)
    print *, 'imag ...reduced again           = ', aimag(nlm_red1)
    print *, 'imag ...full again              = ', aimag(nlm2)
    print *, '------------------------------'    
    print *, 'Last line (nlm2) should equal the third line (nlm1) for both real and imaginary parts. If so, the reduced <--> full nlm conversion works.'

end program

