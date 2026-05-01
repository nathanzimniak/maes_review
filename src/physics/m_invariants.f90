!==================================================================================================
!> @fichier  m_invariants.f90
!> @auteur   Nathan ZIMNIAK
!> @date     10-11-2025
!> @brief    Calcul des invariants MHD.
!==================================================================================================
module m_invariants
    use m_numerical_init, only: dp
    use m_physical_init,  only: input_parameters, radial_exponents, midplane_parameters

    implicit none
    private
    public :: compute_invariants

contains

    !----------------------------------------------------------------------------------------------
    !> Calcule les invariants MHD.
    !>
    !> @param[in]  x               Position courante.
    !> @param[in]  y               Vecteur des solutions.
    !> @param[in]  input_params    Paramètres d'entrée du modèle.
    !> @param[in]  radial_exps     Exposants radiaux.
    !> @param[in]  midplane_params Paramètres physiques au plan médian.
    !> @param[out] mhd_invariants  Invariants MHD.
    !----------------------------------------------------------------------------------------------
    pure function compute_invariants(x, y, input_params, radial_exps, midplane_params) result(mhd_invariants)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params

        real(dp) :: mhd_invariants(4), r_psi, d_psi, uz, bz, g_a, &
                    Epol, Egrav, Emag, Eth, omega_star, lambdaBP, &
                    kappaBP, E

        ! Alias locaux pour lisibilité.
        associate(ep     => input_params%ep,       &
                  be     => radial_exps%be,        &
                  mu     => input_params%mu,       &
                  alpha4 => radial_exps%alpha4,    &
                  delta  => midplane_params%delta, &
                  ms     => midplane_params%ms,    &
                  q      => midplane_params%q)

        ! Grandeurs intermédiaires nécessaires au calcul des invariants.
        r_psi = y(6)**(-1.0_dp/be)
        d_psi = be*y(6)*y(2)/(y(3) + x*y(2))
        uz    = ep*y(3)
        bz    = y(6) - x*d_psi/be

        ! Évaluation des invariants MHD.
        omega_star = (delta*y(5) - q*ms*ep*y(1)*uz/bz)*r_psi**(-1.5_dp)
        lambdaBP   = (delta*y(5) - mu*q*ep*y(1)*bz/(ms*uz*exp(y(4))))*r_psi**(0.5_dp)
        kappaBP    = (exp(y(4))*uz/bz)*r_psi**(alpha4/2.0_dp)*ms/(mu*ep)
        g_a        = 1.0_dp - delta*y(5)*r_psi**(-1.5_dp)/omega_star
        Epol       = 0.50_dp*ms*ms*ep*ep*(y(2)*y(2) + uz*uz)/r_psi
        Egrav      = - 1.0_dp/(sqrt(1.0_dp + x*x*ep*ep)*r_psi)
        Emag       = - 0.5_dp*omega_star*omega_star*r_psi*r_psi*(1.0_dp - g_a*g_a)
        Eth        = (5.0_dp/2.0_dp)*ep*ep*y(7)/r_psi
        E          = Epol + Eth + Egrav + Emag + lambdaBP*omega_star

        ! Stockage des invariants.
        mhd_invariants(1) = omega_star
        mhd_invariants(2) = lambdaBP
        mhd_invariants(3) = kappaBP
        mhd_invariants(4) = E
        end associate
    end function compute_invariants
end module m_invariants