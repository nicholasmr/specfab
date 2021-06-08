! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

program demo

    use specfab    
    use netcdf
    
    implicit none

    integer, parameter :: dp = 8

    ! Numerics
    real, parameter    :: dt = 0.02 ! Time-step size
    integer, parameter :: Nt = 230  ! Number of time steps
    integer            :: Lcap = 12     ! Expansion series truncation
    real(kind=dp)      :: nu0  = 5.0d-3 ! Regularization magnitude calibrated for demo with L=12
    
    ! Constants and argv strings    
    integer :: ii,tt ! loop vars
    character(len=5) :: arg_exp ! experiment type (see below)

    ! Fabric state and evolution
    complex(kind=dp), allocatable :: nlm(:), dndt(:,:), dndt_ROT(:,:), dndt_REG(:,:) ! Series expansion coefs and evolution matrix
    real(kind=dp) :: ugrad(3,3), eps(3,3), omg(3,3) ! Large-scale deformation

    ! For dumping state to netCDF
    complex(kind=dp), allocatable   :: nlm_save(:,:)
    real(kind=dp)                   :: eigvals_save(3,Nt), Eeiej_lin_save(3,3,Nt), Eeiej_nlin_save(3,3,Nt), Epijqij_lin_save(3,Nt), Epijqij_nlin_save(3,Nt)
    real(kind=dp), dimension(3,Nt)  :: e1_save,e2_save,e3_save, p23_save,p12_save,p13_save, q23_save,q12_save,q13_save
    character(len=30) :: fname_sol
    integer :: ncid, c_did, time_did, eig_did, dim_did, pair_did ! Dimension IDs
    integer :: id_cre,id_cim,id_lm, id_eig, id_e1,id_e2,id_e3, id_p23,id_p12,id_p13, id_q23,id_q12,id_q13 ! Var IDs
    integer :: id_Eeiej_lin, id_Eeiej_nlin, id_Epijqij_lin, id_Epijqij_nlin ! Var IDs

    if (command_argument_count() .ne. 1) then
        print *,'usage: ./demo ugrad'
        call exit(0)
    end if

    !-------------------------------------------------------------------
    ! Velocity gradient tensor
    !-------------------------------------------------------------------

    call get_command_argument(1, arg_exp)
    select case (arg_exp)

        ! RECALL COLUMN FIRST IN FORTRAN
        
        ! Uniaxial compression (uc) and uniaxial extension (ue)
        case ('uc_xx', 'ue_xx')
            ugrad = reshape([-1.0,0.,0., 0.,0.5,0., 0.,0.,0.5], [3, 3])
        case ('uc_yy', 'ue_yy')
            ugrad = reshape([0.5,0.,0., 0.,-1.0,0., 0.,0.,0.5], [3, 3])
        case ('uc_zz', 'ue_zz')
            ugrad = reshape([0.5,0.,0., 0.,0.5,0., 0.,0.,-1.], [3, 3])

        ! Confined compression (cc)
        case ('cc_zx')
            ugrad = reshape([1.,0.,0., 0.,0.,0., 0.,0.,-1.], [3, 3]) ! confined in y
        case ('cc_zy')
            ugrad = reshape([0.,0.,0., 0.,1.,0., 0.,0.,-1.], [3, 3]) ! confined in x
        case ('cc_yx')
            ugrad = reshape([1.,0.,0., 0.,-1.,0., 0.,0.,0.], [3, 3]) ! confined in z
                    
        ! Simple shear (ss)
        case ('ss_xz')
            ugrad = reshape([0.,0.,0., 0.,0.,0., 1.,0.,0.], [3, 3]) 
        case ('ss_xy')
            ugrad = reshape([0.,0.,0., 1.,0.,0., 0.,0.,0.], [3, 3]) 
        case ('ss_yz')
            ugrad = reshape([0.,0.,0., 0.,0.,0., 0.,1.,0.], [3, 3]) 

        ! Ridig rotation (rr)
        case ('rr_xz')
            ugrad = reshape([0.,0.,-1., 0.,0.,0., +1.,0.,0.], [3, 3])
        case ('rr_xy')
            ugrad = reshape([0.,-1.,0., +1.,0.,0., 0.,0.,0.], [3, 3]) 
        case ('rr_yz')
            ugrad = reshape([0.,0.,0., 0.,0.,-1., 0.,+1.,0.], [3, 3]) 

        case default
            print *,'argv error: valid "ugrad" are "uc_ij", "ue_ij", "cc_ij", "ss_ij", "rr_ij" where i,j are any of x,y,z'
            call exit(0)
    end select

    select case (arg_exp)
        case ('ue_xx','ue_yy','ue_zz')
            ugrad = -1*ugrad
    end select

    eps = (ugrad+transpose(ugrad))/2 ! strain-rate
    omg = (ugrad-transpose(ugrad))/2 ! spin
            
    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

    call initspecfab(Lcap) ! nlm_len is now defined (number of expansion coeffcients, i.e. #DOFs)
    allocate(nlm(nlm_len))
    nlm = [(0,ii=1,nlm_len)] ! Expansion coefs "n_l^m" are saved in the 1D array "nlm". Corresponding (l,m) values for the i'th coef (i.e. nlm(i)) are (l,m) = (lm(1,i),lm(2,i))
    allocate(nlm_save(nlm_len,Nt))
    allocate(dndt(nlm_len,nlm_len))
    allocate(dndt_ROT(nlm_len,nlm_len))
    allocate(dndt_REG(nlm_len,nlm_len))    

    
    select case (arg_exp)
        case ('rr_xz','rr_xy','rr_yz')
            ! Init with some anisotropy
            nlm(1) = sqrt(cmplx(1/2.,0))
            nlm(3) = sqrt(cmplx(1/2.,0))
        case default
            ! Init with isotropy
            nlm(1) = (1,0)
            nlm(1) = nlm(1)/f_ev_c0(nlm(1)) ! Normalize
    end select
    
    write(*,"(A13,I4,A5,F12.10,A4,I2,A10,I3,A1)") 'Numerics: Nt=', Nt, ', dt=', dt, ', L=', Lcap, ' (nlm_len=',nlm_len,')'

    !-------------------------------------------------------------------
    ! Integrate
    !-------------------------------------------------------------------

    call savestate(nlm, 1) ! Save initial state    
    dndt_ROT = dndt_ij_LATROT(eps,omg, 0*eps,0d0,0d0,0d0, 1d0) ! Assume constant strain-rate and spin with Taylor style plastic spin for lattice rotation (beta=1).
    dndt_REG = f_nu_eps(nu0, eps) * dndt_ij_REG() ! Regularization: nu * (reg. mag.) *  reg. matrix
            
    do tt = 2, Nt
!        write(*,"(A9,I3)") '*** Step ', tt
        dndt = dndt_ROT + dndt_REG 
        nlm = nlm + dt * matmul(dndt, nlm) ! Spectral coefficients evolve by a linear transformation
        call savestate(nlm, tt)
    end do
    
    !-------------------------------------------------------------------
    ! Dump solution to netCDF
    !-------------------------------------------------------------------
    
    write (fname_sol,"('solutions/ROT_',A5,'.nc')") arg_exp
    call check( nf90_create(fname_sol, NF90_CLOBBER, ncid) )
    
    call check(nf90_put_att(ncid,NF90_GLOBAL, "tsteps", Nt))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "dt",     dt))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "nu",     f_nu_eps(nu0, eps) ))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "L",      Lcap))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "ugrad",  reshape(ugrad, [size(ugrad)]) ))
    
    call check(nf90_put_att(ncid,NF90_GLOBAL, "Eca_opt_lin", Eca_opt_lin))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "Ecc_opt_lin", Ecc_opt_lin))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "Eca_opt_nlin", Eca_opt_nlin))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "Ecc_opt_nlin", Ecc_opt_nlin))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "alpha_opt_lin",  alpha_opt_lin))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "alpha_opt_nlin", alpha_opt_nlin))
    
    call check( nf90_def_dim(ncid, "DOF",    nlm_len,   c_did) )
    call check( nf90_def_dim(ncid, "tstep",  Nt,        time_did) )
    call check( nf90_def_dim(ncid, "eigval", 3,         eig_did) )
    call check( nf90_def_dim(ncid, "dim",    3,         dim_did) )
    call check( nf90_def_dim(ncid, "pair",   2,         pair_did) )

    call check( nf90_def_var(ncid, "lm",    NF90_INT,    [pair_did, c_did], id_lm) )
    call check( nf90_def_var(ncid, "c_re",  NF90_DOUBLE, [c_did,   time_did], id_cre) )
    call check( nf90_def_var(ncid, "c_im",  NF90_DOUBLE, [c_did,   time_did], id_cim) )
    
    call check( nf90_def_var(ncid, "eigvals", NF90_DOUBLE, [eig_did, time_did], id_eig) )
    call check( nf90_def_var(ncid, "Eeiej_nlin",   NF90_DOUBLE, [dim_did,dim_did, time_did], id_Eeiej_nlin) )
    call check( nf90_def_var(ncid, "Eeiej_lin",    NF90_DOUBLE, [dim_did,dim_did, time_did], id_Eeiej_lin) )
    call check( nf90_def_var(ncid, "e1",      NF90_DOUBLE, [dim_did, time_did], id_e1) )
    call check( nf90_def_var(ncid, "e2",      NF90_DOUBLE, [dim_did, time_did], id_e2) )
    call check( nf90_def_var(ncid, "e3",      NF90_DOUBLE, [dim_did, time_did], id_e3) )
    call check( nf90_def_var(ncid, "Epijqij_nlin", NF90_DOUBLE, [dim_did, time_did], id_Epijqij_nlin) )
    call check( nf90_def_var(ncid, "Epijqij_lin",  NF90_DOUBLE, [dim_did, time_did], id_Epijqij_lin) )
    call check( nf90_def_var(ncid, "p23",     NF90_DOUBLE, [dim_did, time_did], id_p23) )
    call check( nf90_def_var(ncid, "p12",     NF90_DOUBLE, [dim_did, time_did], id_p12) )
    call check( nf90_def_var(ncid, "p13",     NF90_DOUBLE, [dim_did, time_did], id_p13) )
    call check( nf90_def_var(ncid, "q23",     NF90_DOUBLE, [dim_did, time_did], id_q23) )
    call check( nf90_def_var(ncid, "q12",     NF90_DOUBLE, [dim_did, time_did], id_q12) )
    call check( nf90_def_var(ncid, "q13",     NF90_DOUBLE, [dim_did, time_did], id_q13) )
    
    call check( nf90_enddef(ncid) )
    
    call check( nf90_put_var(ncid, id_cre,   real(nlm_save)) )
    call check( nf90_put_var(ncid, id_cim,   aimag(nlm_save)) )
    call check( nf90_put_var(ncid, id_lm,    lm(:,1:nlm_len)) )
    
    call check( nf90_put_var(ncid, id_eig,  eigvals_save) )
    call check( nf90_put_var(ncid, id_Eeiej_lin,  Eeiej_lin_save) )
    call check( nf90_put_var(ncid, id_Eeiej_nlin, Eeiej_nlin_save) )
    call check( nf90_put_var(ncid, id_e1,   e1_save) )
    call check( nf90_put_var(ncid, id_e2,   e2_save) )
    call check( nf90_put_var(ncid, id_e3,   e3_save) )
    call check( nf90_put_var(ncid, id_Epijqij_lin,  Epijqij_lin_save) )
    call check( nf90_put_var(ncid, id_Epijqij_nlin, Epijqij_nlin_save) )
    call check( nf90_put_var(ncid, id_p23,  p23_save) )
    call check( nf90_put_var(ncid, id_p12,  p12_save) )
    call check( nf90_put_var(ncid, id_p13,  p13_save) )
    call check( nf90_put_var(ncid, id_q23,  q23_save) )
    call check( nf90_put_var(ncid, id_q12,  q12_save) )
    call check( nf90_put_var(ncid, id_q13,  q13_save) )
    
    call check( nf90_close(ncid) )

    print *, 'Solution dumped in ', fname_sol
    print *, "Plot result:"
    write(*,"(A25,A5)") "python3 plot_demo_LATROT.py ", arg_exp

contains

    subroutine savestate(nlm, tt)
    
        implicit none
    
        complex(kind=dp), intent(in) :: nlm(nlm_len)
        integer, intent(in) :: tt
    
        nlm_save(:,tt)  = nlm

        call eigenframe(nlm, e1_save(:,tt),e2_save(:,tt),e3_save(:,tt), eigvals_save(:,tt)) 
        call pqframe(nlm, p23_save(:,tt),p12_save(:,tt),p13_save(:,tt), q23_save(:,tt),q12_save(:,tt),q13_save(:,tt))

        ! Linear (n'=1) mixed Taylor--Sachs enhancements        
        Eeiej_lin_save(:,:,tt) = Eeiej(  nlm, Ecc_opt_lin, Eca_opt_lin, alpha_opt_lin, 1)
        Epijqij_lin_save(:,tt) = Epijqij(nlm, Ecc_opt_lin, Eca_opt_lin, alpha_opt_lin, 1)
        
        ! Nonlinear (n'=3) Sachs enhancements
        Eeiej_nlin_save(:,:,tt) = Eeiej(  nlm, Ecc_opt_nlin, Eca_opt_nlin, alpha_opt_nlin, 3)
        Epijqij_nlin_save(:,tt) = Epijqij(nlm, Ecc_opt_nlin, Eca_opt_nlin, alpha_opt_nlin, 3)

    end
    
    subroutine check(status)
        implicit none
        integer, intent (in) :: status
        if(status /= nf90_noerr) then 
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check  
    
end program

