! N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2022

program demo

    use specfab    
    
    implicit none

    integer, parameter :: dp = 8
    real(kind=dp)      :: A, fa,fb!, m(3),t(3)
    integer            :: n
    
    real(kind=dp), dimension(3), parameter :: x1 = [1,0,0], x2 = [0,1,0], x3 = [0,0,1] ! x,y,z dir.

    call initspecfab(4) ! Not going to use specfab, but is good practice to initialize it even when using static routines (as is the case here).
   
    A = 3.5d0 ! Rate factor (arbitrary value)
    n = 3     ! Flow law exponent
   
    ! rho -> rho_ice
    fa = 2.0d0
    fb = 0.12d0
   
    ! tho -> rho_snow
    fa = 20.0d0
    fb = 6.0d0
   
    print *, '--------------------------------------------------------'    
    print *, 'Test self consistency of forward and inverse rheologies (verify that tau0 = tau(eps(tau0)) for different stress tensors)'
    print *, '--------------------------------------------------------'
    call sig_of_eps_of_sig(tau_vv_cmp(x1), A,n, fa,fb) 
    call sig_of_eps_of_sig(tau_vv_cmp(x2), A,n, fa,fb) 
    call sig_of_eps_of_sig(tau_vv_cmp(x3), A,n, fa,fb) 
    call sig_of_eps_of_sig(tau_vw(x1,x2), A,n, fa,fb) 
    call sig_of_eps_of_sig(tau_vw(x1,x3), A,n, fa,fb) 
    call sig_of_eps_of_sig(tau_vw(x2,x3), A,n, fa,fb) 

contains

subroutine sig_of_eps_of_sig(sig_in, A,n, fa,fb) 

    implicit none

    integer, parameter        :: dp = 8
    real(kind=dp), intent(in) :: sig_in(3,3), A, fa,fb
    integer, intent(in)       :: n
    real(kind=dp)             :: eps(3,3), sig(3,3)

    eps = rheo_fwd__isotropic_compressible(sig_in, A,n, fa,fb)
    sig = rheo_rev__isotropic_compressible(eps,    A,n, fa,fb)
    
    print *, 'sig0             = ', sig_in
    print *, 'sig(eps(sig0))   = ', sig
    print *, 'norm2(abs(diff)) = ', norm2(abs(sig-sig_in))
    print *, '----------------'
        
end
    
function tau_vv_cmp(v) 
    ! v--v compression/extension 
    implicit none
    real(kind=dp), intent(in) :: v(3)
    integer, parameter :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=dp) :: tau_vv_cmp(3,3)
    tau_vv_cmp = identity/2.25d0 - outerprod(v,v) ! This version does (on purpose) not have a vanishing trace!
end

end program

