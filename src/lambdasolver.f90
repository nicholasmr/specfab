! N. M. Rathmann <rathmann@nbi.ku.dk>, 2021

! Module for solving for lambda'_i given nglen, E_11, E_22, E_33 
! Required for Martin's forward orthotropic rheology

module lambdasolver

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    integer, parameter :: NMAX = 3, LWA = 3+(15*NMAX+3*NMAX*NMAX)/2 

    real(kind=dp), parameter :: TOL = 1d-6 

    ! Common block 
    integer :: nglen_saved
    real(kind=dp) :: Eii_saved(3)

contains

    function lambdaprime(nglen, Eii) result(X)

        implicit none

        ! Input
        integer :: nglen
        real(kind=dp) :: Eii(3)

        ! For DNQFJ
        integer :: IOPT(5)
        real(kind=dp) :: FVEC(NMAX), WA(LWA), X(NMAX)
!        real(kind=dp) :: FNORM
       
        ! Constant parameters are shared with (passed to) DNQFJ using the common block
        nglen_saved = nglen 
        Eii_saved = Eii

        ! Numerically calcualte Jacobian matrix
        IOPT(4) = 1
        IOPT(5) = 0
        
        ! Initial guess
        X = 1.0D0

        call DNQSOL(DNQFJ, NMAX, X, FVEC, TOL, IOPT, WA, LWA) 

        ! Debug
        !FNORM = DNRM2(N,FVEC,1)
        !print *, X
        !print *, 'Termination status: IOPT(1),IOPT(2),IOPT(3) = ', IOPT(1), IOPT(2), IOPT(3)
        !print *, 'FNORM: ', FNORM    

    end
    
        
    subroutine DNQFJ(N, X, FVEC, FJAC, IFLAG)

        implicit none

        integer       :: IFLAG, N
        real(kind=dp) :: FJAC(N,N), FVEC(N), X(N)
        real(kind=dp) :: expo, f0, f1

        expo = (nglen_saved-1)/2.0d0

        f0 = 3.0d0/8 
        f1 = 3.0d0/16

        if (IFLAG .eq. 1) then
            FVEC(1) = f0 * (f1)**(expo) * (x(2)+x(3))*(x(2)**2+x(2)*x(3)+x(3)**2)**(expo) - Eii_saved(1)
            FVEC(2) = f0 * (f1)**(expo) * (x(3)+x(1))*(x(3)**2+x(3)*x(1)+x(1)**2)**(expo) - Eii_saved(2)
            FVEC(3) = f0 * (f1)**(expo) * (x(1)+x(2))*(x(1)**2+x(1)*x(2)+x(2)**2)**(expo) - Eii_saved(3)
        endif
        
        return
    end

end module lambdasolver
