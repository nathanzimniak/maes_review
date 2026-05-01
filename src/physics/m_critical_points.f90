!==================================================================================================
!> @fichier  m_critical_points.f90
!> @auteur   Nathan ZIMNIAK
!> @date     10-11-2025
!> @brief    Traitement des points critiques du modèle MHD.
!==================================================================================================
module m_critical_points
    use, intrinsic :: ieee_arithmetic
    use m_numerical_init, only: dp, ip, X_STEP_DEFAULT, JUMP_ERROR, NO_ERROR, integration_history
    use m_physical_init,  only: input_parameters, radial_exponents, midplane_parameters
    use m_mhd,            only: rhs_id
    use m_invariants,     only: compute_invariants
    use m_mach_numbers,   only: compute_mach_sm, compute_mach_a, compute_mach_fm

    implicit none
    private
    public :: cross_midplane, cross_sm_point, cross_a_point

contains

    !----------------------------------------------------------------------------------------------
    !> Développement limité près du plan médian du disque.
    !>
    !> @param[in]    x               Position courante.
    !> @param[inout] y               Vecteur des solutions.
    !> @param[in]    input_params    Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps     Exposants radiaux.
    !> @param[in]    midplane_params Paramètres physiques au plan médian.
    !----------------------------------------------------------------------------------------------
    pure subroutine cross_midplane(x, y, input_params, radial_exps, midplane_params)
        real(dp), intent(in) :: x
        real(dp), intent(inout) :: y(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params

        real(dp) :: A, B, C, D, E, F, b8, a85,ddf8, aT4, ddf5, G,  &
                    b4, dHeat, a42, ddf0, ddpsi, df3, ddf2, dddf3, &
                    ddf4, d4psi

        ! Alias locaux pour lisibilité.
        associate(xi     => input_params%xi,        &
                  ep     => input_params%ep,        &
                  alphap => input_params%alphap,    &
                  mu     => input_params%mu,        &
                  be     => radial_exps%be,         &
                  alpha0 => radial_exps%alpha0,     &
                  alpha2 => radial_exps%alpha2,     &
                  alpha3 => radial_exps%alpha3,     &
                  alpha4 => radial_exps%alpha4,     &
                  ddmv   => midplane_params%ddmv,   &
                  ddmb   => midplane_params%ddmb,   &
                  ddmp   => midplane_params%ddmp,   &
                  ddmt   => midplane_params%ddmt,   &
                  gamma1 => midplane_params%gamma1, &
                  ms     => midplane_params%ms,     &
                  delta  => midplane_params%delta,  &
                  q      => midplane_params%q,      &
                  Ga     => midplane_params%Ga,     &
                  Lambda => midplane_params%Lambda, &
                  Rm     => midplane_params%Rm,     &
                  Rmt    => midplane_params%Rmt)

        ! Coefficients du développement limité au voisinage du plan médian.
        ddpsi = - be*ep*ep*(Rm + be - 2.0_dp)
        df3   = xi - 1.0_dp
        ddf0  = (-1.0_dp - mu*q*q + mu*Rm*ddpsi/be + ms*ms*ep*ep*df3*(alpha3 &
                - xi) - alphap*sqrt(mu)*ddmb)/(1.0_dp + alphap*sqrt(mu))
        dHeat = 0.0_dp
        b4    = ddf0/gamma1 + ddpsi*(1.0_dp + alpha4*(gamma1 - 1.0_dp))/be - (gamma1 &
                - 1.0_dp)*ep*dHeat/gamma1
        a42   = 0.0_dp
        G     = Ga*Rm*ep*(1.0_dp + (be - 2.0_dp)/Rm)
        b8    = ddmt + G - Rmt*ep*ep*(be + xi) + ep*ep*(2.0_dp - be)*(3.5_dp &
                - be)
        a85   = - delta*Rmt/(q*ms)
        A     = ms*ms*ep*ep*(-alpha2*a42 + 2.0_dp*(xi - alpha2)) + delta*delta*a42 &
                - a42 + mu*Rm*ep*ep
        B     = 2.0_dp*delta*delta
        E     = b4*(alpha2*ms*ms*ep*ep - delta*delta + 1.0_dp) - 3.0_dp*ep*ep &
                + ep*ep*(alpha0 - 2.0_dp)*ddf0 + 2.0_dp*(be - 2.0_dp)*mu*q*q*ep*ep &
                + mu*Rm*ep*ep*(ddmp - 2.0_dp*ddpsi*(be - 1.0_dp - xi)/be) &
                + ep*ep*alphap*sqrt(mu)*((alpha0 - 1.0_dp)*ddmb + (alpha0 - 2.0_dp)*ddf0)
        aT4   = 1.0_dp
        C     = 1.0_dp + a42 - aT4*a42/(1.0_dp + Lambda)
        D     = 1.0_dp - 4.0_dp*xi + a85*Lambda/(1.0_dp + Lambda)
        F     = - b4 - Lambda*(b8 + ddpsi*(be - 2.0_dp)/be)/(1.0_dp + Lambda) + (ddmv &
                + aT4*b4)/(1.0_dp + Lambda)
        ddf5  = (A*F-E*C)/(A*D-B*C)
        ddf8  = b8 + a85*ddf5
        ddf2  = (E*D-F*B)/(A*D-B*C)
        ddf4  = b4 + a42*ddf2
        dddf3 = ddf2*(df3 - 2.0_dp) - 2.0_dp*xi*ddf4
        d4psi = - be*Rm*ep*ep*(ddf2 - ddmp - 3.0_dp*ddpsi*(1.0_dp - 1.0_dp/be)) &
                + ep*ep*ddpsi*(6.0_dp*be - be*be - 8.0_dp)
      
        ! Reconstruction des solutions à partir du développement limité.
        y(1)  = - x + x*x*x * ddf8 / 6.0_dp      ! Bphi
        y(2)  = 1.0_dp + 0.5_dp*x*x * ddf2       ! vr
        y(3)  = x*df3 + x*x*x * dddf3 / 6.0_dp   ! vz
        y(4)  = log(1.0_dp + 0.5_dp*x*x * ddf4)  ! rho
        y(5)  = 1.0_dp + 0.5_dp*x*x * ddf5       ! Omega
        y(6)  = 1.0_dp + 0.5_dp*x*x * ddpsi      ! Psi
        y(8)  = - 1.0_dp + 0.5_dp*x*x * ddf8     ! Bphi'
        y(9)  = x*ddpsi + x*x*x * d4psi / 6.0_dp ! Psi'
        y(10) = 1.0_dp + 0.5_dp*x*x * ddf0       ! P
        y(7)  = y(10) / exp(y(4))                ! T
        end associate
    end subroutine cross_midplane

    !----------------------------------------------------------------------------------------------
    !> Franchit le point magnétosonique lent (SM) par extrapolation.
    !>
    !> @param[inout] x                 Position courante.
    !> @param[inout] y                 Vecteur des solutions.
    !> @param[out]   dydx              Vecteur des dérivées.
    !> @param[in]    input_params      Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps       Exposants radiaux.
    !> @param[in]    midplane_params   Paramètres physiques au plan médian.
    !> @param[in]    integration_hist  Historique d'intégration.
    !> @param[in]    int_state         État d'intégration.
    !> @param[out]   mach_sm           Nombre de Mach SM.
    !> @param[out]   error_state       Code d'erreur.
    !----------------------------------------------------------------------------------------------
    pure subroutine cross_sm_point(x, y, dydx, input_params, radial_exps, midplane_params, integration_hist, int_state, mach_sm, error_state)
        real(dp), intent(inout) :: x
        real(dp), intent(inout) :: y(:)
        real(dp), intent(out) :: dydx(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params
        type(integration_history), intent(in) :: integration_hist
        integer(ip), intent(in) :: int_state
        real(dp), intent(inout) :: mach_sm
        integer(ip), intent(out) :: error_state

        integer(ip) :: k0, k1, k2
        real(dp) :: mhd_invariants(4), omega_star, lambdaBP,       &
                    kappaBP, r_psi0, r_psi1, r_psi2, tan_theta0,   &
                    tan_theta1, tan_theta2, mach_a20, mach_a21,    &
                    mach_a22, dr_psi, dmach_a2, dtan_theta, r_psi, &
                    tan_theta, mach_a2, x_jump, br, bz, uz

        ! Alias locaux pour lisibilité.
        associate(ep      => input_params%ep,       &
                  mu      => input_params%mu,       &
                  be      => radial_exps%be,        &
                  alpha4  => radial_exps%alpha4,    &
                  delta   => midplane_params%delta, &
                  ms      => midplane_params%ms,    &
                  q       => midplane_params%q,     &
                  i_id    => integration_hist%i_id, &
                  x_hist  => integration_hist%x,    &
                  y_hist  => integration_hist%y)

        ! Vérification de la profondeur d’historique nécessaire à l’extrapolation.
        if (i_id <= 3_ip) then
            error_state = JUMP_ERROR
            return
        end if

        ! Sélection des points d’historique utilisés pour l’extrapolation.
        k0 = i_id - 1_ip
        k1 = i_id - 2_ip
        k2 = 1_ip
        x_jump = X_STEP_DEFAULT

        ! Évaluation des invariants utilisés pour reconstruire la solution après le saut.
        mhd_invariants = compute_invariants(x_hist(k2), y_hist(:,k2), input_params, radial_exps, midplane_params)
        omega_star = mhd_invariants(1)
        lambdaBP   = mhd_invariants(2)
        kappaBP    = mhd_invariants(3)

        do
            ! Grandeurs "pilotes" en amont du point critique.
            r_psi0     = y_hist(6,k0)**(-1.0_dp/be)
            r_psi1     = y_hist(6,k1)**(-1.0_dp/be)
            r_psi2     = y_hist(6,k2)**(-1.0_dp/be)

            tan_theta0 = - y_hist(2,k0)/(ep*y_hist(3,k0))
            tan_theta1 = - y_hist(2,k1)/(ep*y_hist(3,k1))
            tan_theta2 = - y_hist(2,k2)/(ep*y_hist(3,k2))

            mach_a20   = kappaBP*kappaBP*mu*ep*ep/(exp(y_hist(4,k0))*r_psi0**alpha4)
            mach_a21   = kappaBP*kappaBP*mu*ep*ep/(exp(y_hist(4,k1))*r_psi1**alpha4)
            mach_a22   = kappaBP*kappaBP*mu*ep*ep/(exp(y_hist(4,k2))*r_psi2**alpha4)

            ! Extrapolation des grandeurs pilotes en aval du point critique.
            dr_psi     = (r_psi1 - r_psi2)/(x_hist(k1)-x_hist(k2))
            dmach_a2   = (mach_a21 - mach_a22)/(x_hist(k1)-x_hist(k2))
            dtan_theta = (tan_theta1 - tan_theta2)/(x_hist(k1)-x_hist(k2))

            r_psi     = r_psi0 + x_jump*dr_psi
            tan_theta = tan_theta0 + x_jump*dtan_theta
            mach_a2   = mach_a20 + x_jump*dmach_a2

            ! Reconstruction des solutions extrapolées à partir des invariants MHD et des grandeurs pilotes.
            x = x_hist(k0) + x_jump
            y(7) = y_hist(7,k0)*r_psi/r_psi0
            y(6) = r_psi**(-be)
            y(4) = log(kappaBP*kappaBP*r_psi**(-alpha4)*mu*ep*ep/mach_a2)
            y(8) = exp(y(4))*y(7)

            br = y(6)*tan_theta/(1.0_dp - x*ep*tan_theta)
            bz = y(6) + ep*x*br
            uz = mu*ep*kappaBP*bz*r_psi**(-alpha4/2.0_dp)/(exp(y(4))*ms)

            y(3) = uz/ep
            y(2) = - uz*tan_theta
            y(1) = - kappaBP*(omega_star/r_psi**(-1.5_dp) - lambdaBP/r_psi**(0.5_dp))/(q*(mach_a2 - 1.0_dp)*r_psi**(alpha4/2.0_dp))
            y(5) = (omega_star/r_psi**(-1.5_dp) + q*ms*ep*y(1)*uz/bz)/delta

            ! Contrôle de validité du point extrapolé.
            call compute_mach_sm(x, y, input_params, radial_exps, midplane_params, int_state, mach_sm)
            call rhs_id(x, y, int_state, input_params, radial_exps, midplane_params, dydx)

            if ((mach_sm < 1.0_dp) .or. (dydx(4) >= 0.0_dp) .or. (dydx(5) <= 0.0_dp) .or. any(ieee_is_nan(y)) .or. (x > x_hist(k0)+0.2_dp)) then
                k1  = k1 - 1_ip
                x_jump = x_jump + 2.0_dp*X_STEP_DEFAULT
                if (k1 == k2) then
                    error_state = JUMP_ERROR
                    return
                end if
            else
                error_state = NO_ERROR
                exit
            end if
        end do
        end associate
    end subroutine cross_sm_point

    !----------------------------------------------------------------------------------------------
    !> Franchit le point d'Alfvén (A) par extrapolation.
    !>
    !> @param[inout] x                 Position courante.
    !> @param[inout] y                 Vecteur des solutions.
    !> @param[out]   dydx              Vecteur des dérivées.
    !> @param[in]    input_params      Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps       Exposants radiaux.
    !> @param[in]    midplane_params   Paramètres physiques au plan médian.
    !> @param[in]    integration_hist  Historique d'intégration.
    !> @param[in]    int_state         État d'intégration.
    !> @param[out]   mach_a            Nombre de Mach A.
    !> @param[out]   error_state       Code d'erreur.
    !----------------------------------------------------------------------------------------------
    pure subroutine cross_a_point(x, y, dydx, input_params, radial_exps, midplane_params, integration_hist, int_state, mach_a, error_state)
        real(dp), intent(inout) :: x
        real(dp), intent(inout) :: y(:)
        real(dp), intent(out) :: dydx(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params
        type(integration_history), intent(in) :: integration_hist
        integer(ip), intent(in) :: int_state
        real(dp), intent(inout) :: mach_a
        integer(ip), intent(out) :: error_state

        integer(ip) :: k0, k1, k2
        real(dp) :: mhd_invariants(4), omega_star, lambdaBP,       &
                    kappaBP, r_psi0, r_psi1, r_psi2, g_a0, g_a1,   &
                    g_a2, mach_a20, mach_a21, mach_a22, dr_psi,    &
                    dmach_a2, dg_a, r_psi, g_a, mach_a2, x_jump,   &
                    br, bz, uz, dT, T0,T1,T2, r_a, dd, T, Epol,    &
                    Egrav, Emag, Eth,varE, tan_theta, up, bp

        ! Alias locaux pour lisibilité.
        associate(ep      => input_params%ep,       &
                  mu      => input_params%mu,       &
                  be      => radial_exps%be,        &
                  alpha4  => radial_exps%alpha4,    &
                  delta   => midplane_params%delta, &
                  ms      => midplane_params%ms,    &
                  q       => midplane_params%q,     &
                  i_id    => integration_hist%i_id, &
                  x_hist  => integration_hist%x,    &
                  y_hist  => integration_hist%y)

        ! Vérification de la profondeur d’historique nécessaire à l’extrapolation.
        if (i_id <= 10_ip) then
            error_state = JUMP_ERROR
            return
        end if

        ! Sélection des points d’historique utilisés pour l’extrapolation.
        k0 = i_id - 1_ip
        k1 = i_id - 2_ip
        k2 = i_id - 10_ip
        x_jump = X_STEP_DEFAULT

        ! Évaluation des invariants utilisés pour reconstruire la solution après le saut.
        mhd_invariants = compute_invariants(x_hist(k2), y_hist(:,k2), input_params, radial_exps, midplane_params)
        omega_star = mhd_invariants(1)
        lambdaBP   = mhd_invariants(2)
        kappaBP    = mhd_invariants(3)

        do
            r_psi = y_hist(6,k2)**(-1.0_dp/be)
            uz    = ep*y_hist(3,k2)
            g_a   = 1.0_dp - delta*y_hist(5,k2)*r_psi**(-1.5_dp)/omega_star
            Epol  = 0.5_dp*ms*ms*ep*ep*(y_hist(2,k2)*y_hist(2,k2) + uz*uz)/r_psi
            Egrav = - 1.0_dp/(sqrt(1.0_dp + x_hist(k2)*x_hist(k2)*ep*ep)*r_psi)
            Emag  = - 0.5_dp*omega_star*omega_star*r_psi*r_psi*(1.0_dp - g_a*g_a)
            Eth   = (5.0_dp/2.0_dp)*ep*ep*y_hist(7,k2)/r_psi
            varE  = Epol + Eth + Egrav + Emag

            ! Grandeurs "pilotes" en amont du point critique.
            r_psi0     = y_hist(6,k0)**(-1.0_dp/be)
            r_psi1     = y_hist(6,k1)**(-1.0_dp/be)
            r_psi2     = y_hist(6,k2)**(-1.0_dp/be)

            mach_a20   = kappaBP*kappaBP*mu*ep*ep/(exp(y_hist(4,k0))*r_psi0**alpha4)
            mach_a21   = kappaBP*kappaBP*mu*ep*ep/(exp(y_hist(4,k1))*r_psi1**alpha4)
            mach_a22   = kappaBP*kappaBP*mu*ep*ep/(exp(y_hist(4,k2))*r_psi2**alpha4)

            g_a0 = mach_a20*(1.0_dp - lambdaBP/(omega_star*r_psi0*r_psi0))/(mach_a20-1.0_dp)
            g_a1 = mach_a21*(1.0_dp - lambdaBP/(omega_star*r_psi1*r_psi1))/(mach_a21-1.0_dp)
            g_a2 = mach_a22*(1.0_dp - lambdaBP/(omega_star*r_psi2*r_psi2))/(mach_a22-1.0_dp)

            T0 = y_hist(7,k0)/r_psi0
            T1 = y_hist(7,k1)/r_psi1
            T2 = y_hist(7,k2)/r_psi2

            ! Extrapolation des grandeurs pilotes en aval du point critique.
            dr_psi   = (r_psi1 - r_psi2)/(x_hist(k1) - x_hist(k2))
            dmach_a2 = (mach_a21 - mach_a22)/(x_hist(k1) - x_hist(k2))
            dg_a     = (g_a1 - g_a2)/(x_hist(k1) - x_hist(k2))
            dT       = (T1 - T2)/(x_hist(k1) - x_hist(k2))

            ! Estimation du saut nécessaire pour dépasser le point d’Alfvén.
            mach_a2 = 1.0_dp + 0.1_dp
            x_jump = (mach_a2 - mach_a20)/dmach_a2

            mach_a2 = mach_a20 + x_jump*dmach_a2
            g_a     = g_a0 + x_jump*dg_a
            T       = T0 + x_jump*dT

            r_a = sqrt(lambdaBP/omega_star)
            dd  = 0.1_dp*g_a/(1.0_dp + 0.1_dp*(1.0_dp - g_a))
            r_psi = r_a*sqrt(1.0_dp + dd)

            ! Reconstruction des solutions extrapolées à partir des invariants MHD et des grandeurs pilotes.
            x = x_hist(k0) + x_jump
            y(7) = T*r_psi
            y(6) = r_psi**(-be)
            y(4) = log(kappaBP*kappaBP*r_psi**(-alpha4)*mu*ep*ep/mach_a2)
            y(8) = exp(y(4))*y(7)

            Egrav = - 1.0_dp/sqrt(1.0_dp + x*x*ep*ep)
            Eth   = (5.0_dp/2.0_dp)*ep*ep*y(7)
            Emag  = - 0.5_dp*omega_star*omega_star*r_psi*r_psi*r_psi*(1.0_dp - g_a*g_a)
            Epol  = varE*r_psi - Egrav - Eth - Emag
            up    = sqrt(2.0_dp*Epol)
            bp    = up*kappaBP*r_psi**(-alpha4/2.0_dp)/mach_a2
            br    = (sqrt((1.0_dp+x*x*ep*ep)*bp*bp-y(6)*y(6)) - x*ep*y(6))/(1.0_dp + x*x*ep*ep)
            bz    = y(6) + x*ep*br
            uz    = mu*ep*kappaBP*bz*r_psi**(-alpha4/2.0_dp)/(exp(y(4))*ms)
            tan_theta = br/bz

            y(2) = - uz*tan_theta
            y(3) = uz/ep
            y(1) = - kappaBP*(omega_star/r_psi**(-1.5_dp) - lambdaBP/r_psi**(0.5_dp))/(q*(mach_a2 - 1.0_dp)*r_psi**(alpha4/2.0_dp))
            y(5) = (omega_star/r_psi**(-1.5_dp) + q*ms*ep*y(1)*uz/bz)/delta

            ! Contrôle de validité du point extrapolé.
            call compute_mach_a(x, y, input_params, radial_exps, midplane_params, int_state, mach_a)
            call rhs_id(x, y, int_state, input_params, radial_exps, midplane_params, dydx)

            if ((mach_a < 1.0_dp) .or. any(ieee_is_nan(y)) .or. (x > x_hist(k0)+1.0_dp)) then
                k1  = k1 - 1_ip
                if ((k1 == k2)) then
                    error_state = JUMP_ERROR
                    return
                end if
            else
                error_state = NO_ERROR
                exit
            end if
        end do
        end associate
    end subroutine cross_a_point
end module m_critical_points