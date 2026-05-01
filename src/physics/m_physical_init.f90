!==================================================================================================
!> @fichier m_physical_init.f90
!> @auteur  Nathan ZIMNIAK
!> @date    10-11-2025
!> @brief   Paramètres physiques du modèle MHD.
!==================================================================================================
module m_physical_init
    use m_numerical_init, only: dp, ip

    implicit none
    private
    public :: compute_radial_exponents, compute_midplane_parameters, initialize_solutions

    ! Paramètres d'entrée du modèle.
    type, public :: input_parameters
        real(dp) :: xi
        real(dp) :: ep
        real(dp) :: alpham
        real(dp) :: chim
        real(dp) :: Pm
        real(dp) :: alphap
        real(dp) :: mu
        real(dp) :: p
    end type input_parameters

    ! Exposants radiaux (r^alpha_i).
    type, public :: radial_exponents
        real(dp) :: be
        real(dp) :: alpha0
        real(dp) :: alpha1
        real(dp) :: alpha2
        real(dp) :: alpha3
        real(dp) :: alpha4
    end type radial_exponents

    ! Grandeurs définies au plan médian du disque.
    type, public :: midplane_parameters
        real(dp) :: q
        real(dp) :: alphav
        real(dp) :: ms
        real(dp) :: delta
        real(dp) :: Rm
        real(dp) :: Rmt
        real(dp) :: Ga
        real(dp) :: Lambda
        real(dp) :: gamma1
        real(dp) :: ddmt
        real(dp) :: ddmv
        real(dp) :: ddmp
        real(dp) :: ddmb
    end type midplane_parameters

contains

    !----------------------------------------------------------------------------------------------
    !> Calcule les exposants radiaux.
    !>
    !> @param[in]  input_param Paramètres d'entrée.
    !> @param[out] radial_exps Exposants radiaux.
    !----------------------------------------------------------------------------------------------
    pure subroutine compute_radial_exponents(input_param, radial_exps)
        type(input_parameters), intent(in) :: input_param
        type(radial_exponents), intent(out) :: radial_exps

        ! Alias locaux pour lisibilité.
        associate(alpha0 => radial_exps%alpha0, &
                  alpha1 => radial_exps%alpha1, &
                  alpha2 => radial_exps%alpha2, &
                  alpha3 => radial_exps%alpha3, &
                  alpha4 => radial_exps%alpha4, &
                  be     => radial_exps%be,     &
                  xi     => input_param%xi)

        be     = 0.75_dp + 0.5_dp*xi
        alpha4 = 2.0_dp*be - 3.0_dp
        alpha1 = 0.5_dp*(alpha4 + 1.0_dp)
        alpha2 = be - 2.0_dp - alpha4/2.0_dp
        alpha3 = alpha2
        alpha0 = alpha4 - 1.0_dp
        end associate
    end subroutine compute_radial_exponents

    !----------------------------------------------------------------------------------------------
    !> Calcule les grandeurs au plan médian du disque.
    !>
    !> @param[in]  Input_param    Paramètres d'entrée.
    !> @param[in]  Radial_exps    Exposants radiaux.
    !> @param[out] Midplane_params Grandeurs au plan médian.
    !----------------------------------------------------------------------------------------------
    pure subroutine compute_midplane_parameters(input_param, radial_exps, midplane_params)
        type(input_parameters), intent(in) :: input_param
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(out) :: midplane_params

        ! Alias locaux pour lisibilité.
        associate(ddmv   => midplane_params%ddmv,   &
                  ddmb   => midplane_params%ddmb,   &
                  ddmp   => midplane_params%ddmp,   &
                  ddmt   => midplane_params%ddmt,   &
                  gamma1 => midplane_params%gamma1, &
                  ms     => midplane_params%ms,     &
                  delta  => midplane_params%delta,  &
                  alphav => midplane_params%alphav, &
                  q      => midplane_params%q,      &
                  Ga     => midplane_params%Ga,     &
                  Lambda => midplane_params%Lambda, &
                  Rm     => midplane_params%Rm,     &
                  Rmt    => midplane_params%Rmt,    &
                  mu     => input_param%mu,        &
                  p      => input_param%p,         &
                  ep     => input_param%ep,        &
                  alpham => input_param%alpham,    &
                  chim   => input_param%chim,      &
                  Pm     => input_param%Pm,        &
                  alphap => input_param%alphap,    &
                  alpha0 => radial_exps%alpha0,    &
                  alpha2 => radial_exps%alpha2)

        ddmv   = - 2.0_dp
        ddmb   = - 2.0_dp
        ddmp   = - 2.0_dp
        ddmt   = - 2.0_dp
        gamma1 = 1.0_dp
        ms     = alpham*p*sqrt(mu)
        delta  = sqrt(1.0_dp + alpha0*ep*ep*(1.0_dp + alphap*sqrt(mu)) - mu*p*ep + alpha2*ms*ms*ep*ep)
        alphav = alpham*Pm*sqrt(mu)
        q      = 0.5_dp*delta*(ms - alphav*ep)/mu
        Ga     = 3.0_dp*sqrt(mu)*chim/(alpham*(ms - alphav*ep))
        Rm     = p/ep
        Rmt    = Rm*chim
        if (alphav <= 1e-10_dp) then
            Lambda = huge(1.0_dp)
        else
            Lambda = 2.0_dp*q*mu/(alphav*ep*delta)
        end if
        end associate
    end subroutine compute_midplane_parameters

    !-----------------------------------------------------------------
    !> Initialise les solutions au plan médian du disque.
    !>
    !> @param[in]  input_param Paramètres d'entrée.
    !> @param[out] y           Vecteur des solutions.
    !> @param[out] dydx        Vecteur des dérivées.
    !-----------------------------------------------------------------
    pure subroutine initialize_solutions(input_param, y, dydx)
        type(input_parameters), intent(in) :: input_param
        real(dp), intent(out) :: y(:)
        real(dp), intent(out) :: dydx(:)

        associate(xi => input_param%xi)

        y(1)  = 0.0_dp      ! Bphi
        y(2)  = 1.0_dp      ! vr
        y(3)  = 0.0_dp      ! vz
        y(4)  = log(1.0_dp) ! rho
        y(5)  = 1.0_dp      ! Omega
        y(6)  = 1.0_dp      ! Psi
        y(7)  = 1.0_dp      ! T
        y(8)  = -1.0_dp     ! Bphi'
        y(9)  = 0.0_dp      ! Psi'
        y(10) = 1.0_dp      ! P

        dydx(1)  = -1.0_dp     ! Bphi'
        dydx(2)  = 0.0_dp      ! vr'
        dydx(3)  = xi - 1.0_dp ! vz'
        dydx(4)  = 0.0_dp      ! rho'
        dydx(5)  = 0.0_dp      ! Omega'
        dydx(6)  = 0.0_dp      ! Psi'
        dydx(7)  = 0.0_dp      ! T'
        dydx(8)  = 0.0_dp      ! Bphi''
        dydx(9)  = 1.0_dp      ! Psi''
        dydx(10) = 0.0_dp      ! P'
        end associate
    end subroutine initialize_solutions
end module m_physical_init