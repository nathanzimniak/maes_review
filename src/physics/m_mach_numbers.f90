!==================================================================================================
!> @fichier  m_mach_numbers.f90
!> @auteur   Nathan ZIMNIAK
!> @date     10-11-2025
!> @brief    Calcul des nombres de Mach caractéristiques du modèle MHD.
!==================================================================================================
module m_mach_numbers
    use m_numerical_init, only: dp, ip, IDEAL_LIN, IDEAL_LOG, integration_history
    use m_physical_init,  only: input_parameters, radial_exponents, midplane_parameters
    use m_mhd,            only: rhs_res, rhs_id

    implicit none
    private
    public :: compute_mach_id, compute_mach_sm, compute_mach_a, compute_mach_fm

contains

    !----------------------------------------------------------------------------------------------
    !> Calcule le nombre de Mach "idéal".
    !>
    !> @param[in]    x                Position courante.
    !> @param[in]    y                Vecteur des solutions.
    !> @param[in]    input_params     Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps      Exposants radiaux.
    !> @param[in]    midplane_params  Paramètres physiques au plan médian.
    !> @param[inout] integration_hist Historique d'intégration.
    !> @param[out]   mach_id          Nombre de Mach "idéal".
    !----------------------------------------------------------------------------------------------
    pure subroutine compute_mach_id(x, y, input_params, radial_exps, midplane_params, integration_hist, mach_id)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params
        type(integration_history), intent(inout) :: integration_hist
        real(dp), intent(out) :: mach_id

        real(dp) :: dydx(size(y)), mv, Turb, Gao, mach_id_p, mach_id_t, &
                    mach_id_v, gamma, dmach_id_p, dmach_id_t, dmach_id_v, y4loc

        ! Alias locaux pour lisibilité.
        associate(ep             => input_params%ep,                 &
                  be             => radial_exps%be,                  &
                  ddmv           => midplane_params%ddmv,            &
                  gamma1         => midplane_params%gamma1,          &
                  ms             => midplane_params%ms,              &
                  delta          => midplane_params%delta,           &
                  q              => midplane_params%q,               &
                  Lambda         => midplane_params%Lambda,          &
                  Rmt            => midplane_params%Rmt,             &
                  x_last         => integration_hist%x_last,         &
                  mach_id_p_last => integration_hist%mach_id_p_last, &
                  mach_id_t_last => integration_hist%mach_id_t_last, &
                  mach_id_v_last => integration_hist%mach_id_v_last)

        ! Grandeurs intermédiaires nécessaires au calcul du nombre de Mach "idéal".
        call rhs_res(x, y, input_params, radial_exps, midplane_params, dydx)
        y4loc = exp(y(4))
        mv    = exp(0.5_dp*ddmv*x*x)
        gamma = gamma1
        Turb  = y4loc*mv
        Gao   = Rmt*delta/(q*ms)
        mach_id_p = 1.0_dp - abs(y(2)*y(6) - (y(3) + x*y(2))*y(9)/be)
        mach_id_t = 1.0_dp - abs(Gao*(- 1.5_dp*y(9)*y(5)/be - y(6)*dydx(5)) + Rmt*ep*ep*(be*y(1)*y(2) + (y(3) + x*y(2))*(dydx(1) - y(1)*dydx(4))))
        mach_id_v = 1.0_dp - abs(Turb/(1.0_dp + Lambda))
        dmach_id_p = 1.0_dp - abs((mach_id_p - mach_id_p_last)/(x - x_last))
        dmach_id_t = 1.0_dp - abs((mach_id_t - mach_id_t_last)/(x - x_last))
        dmach_id_v = 1.0_dp - abs((mach_id_v - mach_id_v_last)/(x - x_last))
        mach_id_p_last = mach_id_p
        mach_id_t_last = mach_id_t
        mach_id_v_last = mach_id_v
        x_last = x

        ! Nombre de Mach "idéal".
        mach_id = minval([mach_id_p, mach_id_t, mach_id_v, dmach_id_p, dmach_id_t, dmach_id_v])
        end associate
    end subroutine compute_mach_id

    !----------------------------------------------------------------------------------------------
    !> Calcule le nombre de Mach magnétosonique lent (SM).
    !>
    !> @param[in]    x                Position courante.
    !> @param[in]    y                Vecteur des solutions.
    !> @param[in]    input_params     Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps      Exposants radiaux.
    !> @param[in]    midplane_params  Paramètres physiques au plan médian.
    !> @param[in]    int_state        État d'intégration courant.
    !> @param[out]   mach_sm          Nombre de Mach SM.
    !----------------------------------------------------------------------------------------------
    pure subroutine compute_mach_sm(x, y, input_params, radial_exps, midplane_params, int_state, mach_sm)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params
        integer(ip), intent(in) :: int_state
        real(dp), intent(out) :: mach_sm

        real(dp) :: dydx(size(y)), mb, gamma, V, Vs2, Va2, Van2, Vsm2, y4loc

        ! Alias locaux pour lisibilité.
        associate(ep      => input_params%ep,        &
                  alphap  => input_params%alphap,    &
                  mu      => input_params%mu,        &
                  be      => radial_exps%be,         &
                  ddmb    => midplane_params%ddmb,   &
                  gamma1  => midplane_params%gamma1, &
                  ms      => midplane_params%ms,     &
                  q       => midplane_params%q)

        ! Grandeurs intermédiaires nécessaires au calcul du nombre de Mach SM.
        call rhs_id(x, y, int_state, input_params, radial_exps, midplane_params, dydx)
        y4loc = exp(y(4))
        mb    = exp(0.5_dp*ddmb*x*x)
        gamma = gamma1
        V     = y(3) + x*y(2)
        Vs2   = gamma*y(7)*(1.0_dp + alphap*sqrt(mu)*mb)
        Va2   = mu*((y(6) - x*dydx(6)/be)*(y(6) - x*dydx(6)/be) + (dydx(6)/(be*ep))*(dydx(6)/(be*ep)) + (q*y(1))*(q*y(1)))/y4loc
        Van2  = mu*y(6)*y(6)/(y4loc*(1.0_dp + ep*ep*x*x))
        Vsm2  = 0.5_dp*(Vs2 + Va2 - sqrt((Vs2+Va2)*(Vs2+Va2) - 4.0_dp*Vs2*Van2))
        
        ! Nombre de Mach SM.
        mach_sm = sqrt((ms*ep*V)*(ms*ep*V)/(Vsm2*(1.0_dp + ep*ep*x*x)))
        end associate
    end subroutine compute_mach_sm

    !----------------------------------------------------------------------------------------------
    !> Calcule le nombre de Mach d'Alfvén (A).
    !>
    !> @param[in]    x                Position courante.
    !> @param[in]    y                Vecteur des solutions.
    !> @param[in]    input_params     Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps      Exposants radiaux.
    !> @param[in]    midplane_params  Paramètres physiques au plan médian.
    !> @param[in]    int_state        État d'intégration courant.
    !> @param[out]   mach_a           Nombre de Mach A.
    !----------------------------------------------------------------------------------------------
    pure subroutine compute_mach_a(x, y, input_params, radial_exps, midplane_params, int_state, mach_a)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params
        integer(ip), intent(in) :: int_state
        real(dp), intent(out) :: mach_a

        real(dp) :: dydx(size(y)), mb, gamma, V, Vs2, Va2, Van2, y4loc

        ! Alias locaux pour lisibilité.
        associate(ep     => input_params%ep,        &
                  alphap => input_params%alphap,    &
                  mu     => input_params%mu,        &
                  be     => radial_exps%be,         &
                  ddmb   => midplane_params%ddmb,   &
                  gamma1 => midplane_params%gamma1, &
                  ms     => midplane_params%ms,     &
                  q      => midplane_params%q)

        ! Grandeurs intermédiaires nécessaires au calcul du nombre de Mach A.
        call rhs_id(x, y, int_state, input_params, radial_exps, midplane_params, dydx)
        y4loc = exp(y(4))
        mb    = exp(0.5_dp*ddmb*x*x)
        gamma = gamma1
        V     = y(3) + x*y(2)
        Vs2   = gamma*y(7)*(1.0_dp + alphap*sqrt(mu)*mb)
        Va2   = mu*((y(6) - x*dydx(6)/be)*(y(6) - x*dydx(6)/be) + (dydx(6)/(be*ep))*(dydx(6)/(be*ep)) + (q*y(1))*(q*y(1)))/y4loc
        Van2  = mu*y(6)*y(6)/(y4loc*(1.0_dp + ep*ep*x*x))

        ! Nombre de Mach A.
        mach_a = sqrt((ms*ep*V)*(ms*ep*V)/(Van2*(1.0_dp + ep*ep*x*x)))
        end associate
    end subroutine compute_mach_a

    !----------------------------------------------------------------------------------------------
    !> Calcule le nombre de Mach magnétosonique rapide (FM).
    !>
    !> @param[in]    x                Position courante.
    !> @param[in]    y                Vecteur des solutions.
    !> @param[in]    input_params     Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps      Exposants radiaux.
    !> @param[in]    midplane_params  Paramètres physiques au plan médian.
    !> @param[in]    int_state        État d'intégration courant.
    !> @param[out]   mach_fm          Nombre de Mach FM.
    !----------------------------------------------------------------------------------------------
    pure subroutine compute_mach_fm(x, y, input_params, radial_exps, midplane_params, int_state, mach_fm)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params
        integer(ip), intent(in) :: int_state
        real(dp), intent(out) :: mach_fm

        real(dp) :: dydx(size(y)), mb, gamma, V, Vs2, Va2, Van2, Vfm2, xloc, y1loc, y4loc, y6loc, dydx6loc

        ! Alias locaux pour lisibilité.
        associate(ep     => input_params%ep,        &
                  alphap => input_params%alphap,    &
                  mu     => input_params%mu,        &
                  be     => radial_exps%be,         &
                  ddmb   => midplane_params%ddmb,   &
                  gamma1 => midplane_params%gamma1, &
                  ms     => midplane_params%ms,     &
                  q      => midplane_params%q)

        ! Grandeurs intermédiaires nécessaires au calcul du nombre de Mach FM.
        call rhs_id(x, y, int_state, input_params, radial_exps, midplane_params, dydx)

        select case (int_state)
        case (IDEAL_LIN)
            xloc     = x
            y1loc    = y(1)
            y4loc    = exp(y(4))
            y6loc    = y(6)
            dydx6loc = dydx(6)
        case (IDEAL_LOG)
            xloc     = exp(x)
            y1loc    = - exp(y(1))
            y4loc    = exp(y(4))
            y6loc    = exp(y(6))
            dydx6loc = (y6loc/xloc)*dydx(6)
        case default; error stop "État d'intégration invalide."
        end select

        mb    = exp(0.5_dp*ddmb*xloc*xloc)
        gamma = gamma1
        V     = y(3) + xloc*y(2)
        Vs2   = gamma*y(7)*(1.0_dp + alphap*sqrt(mu)*mb)
        Va2   = mu*((y6loc - xloc*dydx6loc/be)*(y6loc - xloc*dydx6loc/be) + (dydx6loc/(be*ep))*(dydx6loc/(be*ep)) + (q*y1loc)*(q*y1loc))/y4loc
        Van2  = mu*y6loc*y6loc/(y4loc*(1.0_dp + ep*ep*xloc*xloc))
        Vfm2  = 0.5_dp*(Vs2 + Va2 + sqrt((Vs2+Va2)*(Vs2+Va2) - 4.0_dp*Vs2*Van2))

        ! Nombre de Mach FM.
        mach_fm = sqrt((ms*ep*V)*(ms*ep*V)/(Vfm2*(1.0_dp + ep*ep*xloc*xloc)))
        end associate
    end subroutine compute_mach_fm
end module m_mach_numbers