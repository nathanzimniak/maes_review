!==================================================================================================
!> @fichier  m_mhd.f90
!> @auteur   Nathan ZIMNIAK
!> @date     10-11-2025
!> @brief    Calcul des dérivées du système MHD (résistif et idéal).
!==================================================================================================
module m_mhd
    use m_numerical_init, only: dp, ip, IDEAL_LIN, IDEAL_LOG
    use m_physical_init,  only: input_parameters, radial_exponents, midplane_parameters

    implicit none
    private
    public :: rhs_res, rhs_id

contains

    !----------------------------------------------------------------------------------------------
    !> Calcule les dérivées des solutions (MHD résistive).
    !>
    !> @param[in]  x               Position courante.
    !> @param[in]  y               Vecteur des solutions.
    !> @param[in]  input_params    Paramètres d'entrée du modèle.
    !> @param[in]  radial_exps     Exposants radiaux.
    !> @param[in]  midplane_params Grandeurs physiques au plan médian.
    !> @param[out] dydx            Vecteur des dérivées.
    !----------------------------------------------------------------------------------------------
    pure subroutine rhs_res(x, y, input_params, radial_exps, midplane_params, dydx)
        real(dp), intent(in)  :: x
        real(dp), intent(in)  :: y(:)
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in)  :: midplane_params
        real(dp), intent(out) :: dydx(:)

        real(dp) :: mp, mt, mv, mb, dmt, dmb, V, d2psi, Cte, b4, a43,  &
                    a42, a33, b3, a22, b2, a73, a72, A, B, C, D, E, F, &
                    Deter, Gao, Turb, b7, gamma, y4loc

        ! Alias locaux pour lisibilité.
        associate(xi     => input_params%xi,        &
                  ep     => input_params%ep,        &
                  mu     => input_params%mu,        &
                  alphap => input_params%alphap,    &
                  be     => radial_exps%be,         &
                  alpha0 => radial_exps%alpha0,     &
                  alpha1 => radial_exps%alpha1,     &
                  alpha2 => radial_exps%alpha2,     &
                  alpha3 => radial_exps%alpha3,     &
                  ddmv   => midplane_params%ddmv,   &
                  ddmb   => midplane_params%ddmb,   &
                  ddmp   => midplane_params%ddmp,   &
                  ddmt   => midplane_params%ddmt,   &
                  gamma1 => midplane_params%gamma1, &
                  ms     => midplane_params%ms,     &
                  delta  => midplane_params%delta,  &
                  q      => midplane_params%q,      &
                  Lambda => midplane_params%Lambda, &
                  Rm     => midplane_params%Rm,     &
                  Rmt    => midplane_params%Rmt)

        ! Grandeurs intermédiaires nécessaires à l’évaluation du système.
        y4loc = exp(y(4))
        mp    = exp(0.5_dp*ddmp*x*x)
        mt    = exp(0.5_dp*ddmt*x*x)
        mv    = exp(0.5_dp*ddmv*x*x)
        mb    = exp(0.5_dp*ddmb*x*x)
        dmt   = ddmt*x*mt
        dmb   = ddmb*x*mb
        Turb  = y4loc*mv
        V     = y(3) + x*y(2)
        b4    = (xi - 1.0_dp)*y(2)/V
        b7    = y(10)*b4 - y(9)*y(10)/(be*y(6))
        gamma = gamma1
        d2psi = (y(2)*y(6) - V*y(9)/be)/mp
        Cte   = Lambda/(1.0_dp + Lambda)
        a43   = - 1.0_dp/V
        a42   = x*a43
        a33   = ms*ms*ep*ep*V
        b3    = alpha3*ms*ms*ep*ep*y(2)*y(3) - x*(1.0_dp + ep*ep*x*x)**(-1.5_dp) - b7/y4loc - q*q*mu*y(1)*y(8)/y4loc + Rm*mu*y(9)*d2psi/(y4loc*be)
        a22   = a33
        b2    = alpha2*ms*ms*ep*ep*y(2)*y(2) - delta*delta*y(5)*y(5) + (1.0_dp + ep*ep*x*x)**(-1.5_dp) + ep*ep*(alpha0*y(10) - x*b7)/y4loc + q*q*mu*ep*ep*y(1)*(alpha1*y(1) - x*y(8))/y4loc - Rm*mu*ep*ep*d2psi*(y(6)- x*y(9)/be)/y4loc
        a73   = - gamma*y(10)/V
        a72   = x*a73
        A     = a33 + a73*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc
        D     = a22 + x*ep*ep*a72*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc
        B     = - a72*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc
        C     = - a73*x*ep*ep*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc
        E     = b3 - alphap*sqrt(mu)*(dmb*y(10) + b7*mb)/y4loc
        F     = b2 + alphap*sqrt(mu)*ep*ep*((alpha0*mb - x*dmb)*y(10) - b7*x*mb)/y4loc
        Deter = A*D - B*C
        Gao   = Rmt*delta/(q*ms)

        ! Assemblage final du second membre résistif.
        dydx(1)  = y(8)
        dydx(2)  = (A*F + E*C)/Deter
        dydx(3)  = (E*D + B*F)/Deter
        dydx(4)  = b4 + a42*dydx(2) + a43*dydx(3)
        dydx(5)  = (y(2)*y(5) + (Cte*(y(6)*y(8) - alpha1*y(9)*y(1)/be) - Turb/(1.0_dp + Lambda))/y4loc)/(2.0_dp*V)
        dydx(6)  = y(9)
        dydx(8)  = (- Gao*(1.5_dp*y(9)*y(5)/be + y(6)*dydx(5)) - dmt*(y(8)*(1.0_dp + ep*ep*x*x) - alpha1*ep*ep*x*y(1)) - ep*ep*mt*(alpha1*(be - 2.5_dp)*y(1) - x*y(8)*(2.0_dp*be - 4.5_dp)) + Rmt*ep*ep*(be*y(1)*y(2) + V*(y(8) - y(1)*dydx(4))))/(mt*(1.0_dp + ep*ep*x*x))
        dydx(9)  = (- be*Rm*ep*ep*d2psi + ep*ep*(be*(2.0_dp - be)*y(6) + (2.0_dp*be-3.0_dp)*x*y(9)))/(1.0_dp + ep*ep*x*x)
        dydx(10) = b7 + a72*dydx(2) + a73*dydx(3)
        dydx(7)  = dydx(10)/y4loc - dydx(4)*y(7)
        end associate
    end subroutine rhs_res

    !----------------------------------------------------------------------------------------------
    !> Calcule les dérivées des solutions (MHD idéale).
    !>
    !> @param[in]  x               Position courante.
    !> @param[in]  y               Vecteur des solutions.
    !> @param[in]  int_state       État d'intégration.
    !> @param[in]  input_params    Paramètres d'entrée du modèle.
    !> @param[in]  radial_exps      Exposants radiaux.
    !> @param[in]  midplane_params Grandeurs physiques au plan médian.
    !> @param[out] dydx            Vecteur des dérivées.
    !----------------------------------------------------------------------------------------------
    pure subroutine rhs_id(x, y, int_state, input_params, radial_exps, midplane_params, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        integer(ip), intent(in)  :: int_state
        type(input_parameters), intent(in)  :: input_params
        type(radial_exponents), intent(in)  :: radial_exps
        type(midplane_parameters), intent(in)  :: midplane_params
        real(dp), intent(out) :: dydx(:)

        real(dp) :: mb, dmb, V, Cte, b4, a43, a42, a33, b3, a22, b2,&
                  a73, a72, A, B, C, D, E, F, Deter, b7, gamma, b5, &
                  a51, b1, a14, a15, alpha14, beta1, b0, a03, a02,  &
                  a31, a30, a21, a20, QQ, xloc, y1loc, y4loc, y6loc,&
                  dydx6loc

        ! Alias locaux pour lisibilité.
        associate(xi        => input_params%xi,        &
                  ep        => input_params%ep,        &
                  alphap    => input_params%alphap,    &
                  mu        => input_params%mu,        &
                  be        => radial_exps%be,         &
                  alpha0    => radial_exps%alpha0,     &
                  alpha1    => radial_exps%alpha1,     &
                  alpha2    => radial_exps%alpha2,     &
                  alpha3    => radial_exps%alpha3,     &
                  ddmb      => midplane_params%ddmb,   &
                  gamma1    => midplane_params%gamma1, &
                  ms        => midplane_params%ms,     &
                  delta     => midplane_params%delta,  &
                  q         => midplane_params%q,      &
                  Lambda    => midplane_params%Lambda)

        ! Grandeurs intermédiaires nécessaires à l’évaluation du système.
        select case (int_state)
        case (IDEAL_LIN)
            xloc  = x
            y1loc = y(1)
            y4loc = exp(y(4))
            y6loc = y(6)
        case (IDEAL_LOG)
            xloc  = exp(x)
            y1loc = -exp(y(1))
            y4loc = exp(y(4))
            y6loc = exp(y(6))
        case default; error stop
        end select

        mb    = exp(0.5_dp*ddmb*xloc*xloc)
        dmb   = ddmb*xloc*mb
        V       = y(3) + xloc*y(2)
        b7    = y(8)*y(2)*(xi - 2.0_dp)/V
        gamma = gamma1

        select case (int_state)
        case (IDEAL_LIN)
            dydx(6) = be*y(6)*y(2)/V
            dydx6loc = dydx(6)
        case (IDEAL_LOG)
            dydx(6)  = be*xloc*y(2)/V
            dydx6loc = (y6loc/xloc)*dydx(6)
        case default
        end select

        Cte     = Lambda/(1.0_dp + Lambda)
        b5      = 0.5_dp*y(2)*y(5)/V - 0.5_dp*alpha1*Cte*dydx6loc*y1loc/(be*y4loc*V)
        a51     = 0.5_dp*y6loc*Cte/(V*y4loc)
        b1      = 1.5_dp*delta*dydx6loc*y(5)/(be*q*ms*ep*ep*V) - y1loc*dydx6loc/y6loc
        a14     = y1loc
        a15     = y6loc*delta/(q*ms*ep*ep*V)
        alpha14 = a14/(1.0_dp - a15*a51)
        beta1   = (b1 + a15*b5)/(1.0_dp - a15*a51)
        b0      = dydx6loc*dydx6loc*(1.0_dp + ep*ep*xloc*xloc)*(be - 1.0_dp)/(be*y6loc) - ep*ep*(be*(2.0_dp - be)*y6loc + (2.0_dp*be - 3.0_dp)*xloc*dydx6loc)
        a03     = - dydx6loc*(1.0_dp + ep*ep*xloc*xloc)/V
        a02     = - y(3)*a03/y(2)
        b4      = (xi - 1.0_dp)*y(2)/V
        a43     = - 1.0_dp/V
        a42     = xloc*a43
        a73     = - gamma*y(8)/V
        a72     = xloc*a73
        a33     = ms*ms*ep*ep*V
        a31     = - q*q*mu*y1loc/y4loc
        a30     = - mu*dydx6loc/(y4loc*be*be*ep*ep)
        b3      = alpha3*ms*ms*ep*ep*y(2)*y(3) - xloc*(1.0_dp + ep*ep*xloc*xloc)**(-1.5_dp)- b7/y4loc
        a22     = a33
        a21     = xloc*ep*ep*a31
        a20     = xloc*ep*ep*a30 + mu*y6loc/(y4loc*be)
        b2      = alpha2*ms*ms*ep*ep*y(2)*y(2) - delta*delta*y(5)*y(5) + (1.0_dp + ep*ep*xloc*xloc)**(-1.5_dp) + ep*ep*(alpha0*y(8) - xloc*b7)/y4loc + q*q*mu*ep*ep*alpha1*y1loc*y1loc/y4loc
        QQ      = a31*alpha14
        A       = a33 + a73*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc - a30*a03 - QQ*a43
        D       = a22 + xloc*ep*ep*a72*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc - a20*a02 - QQ*xloc*ep*ep*a42
        B       = - a72*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc + a30*a02 + QQ*a42
        C       = - a73*xloc*ep*ep*(1.0_dp + alphap*sqrt(mu)*mb)/y4loc + a20*a03 + QQ*xloc*ep*ep*a43
        Deter   = A*D - B*C
        E       = b3 - alphap*sqrt(mu)*(dmb*y(8) + b7*mb)/y4loc + a30*b0 + a31*beta1 + QQ*b4
        F       = b2 + alphap*sqrt(mu)*ep*ep*((alpha0*mb - xloc*dmb)*y(8) - b7*xloc*mb)/y4loc + a20*b0 + a21*beta1 + QQ*xloc*ep*ep*b4

        if (int_state == IDEAL_LOG) then
            E = xloc * E
            F = xloc * F
        end if

        ! Assemblage final du second membre idéal.
        dydx(2) = (A*F + E*C)/Deter
        dydx(3) = (E*D + B*F)/Deter

        if (int_state == IDEAL_LIN) then
            dydx(4) = b4 + a42*dydx(2) + a43*dydx(3)
            dydx(1) = beta1 + alpha14*dydx(4)
            dydx(5) = b5 + a51*dydx(1)
            dydx(8) = b7 + a72*dydx(2) + a73*dydx(3)
        else
            dydx(4) = xloc*b4 + a42*dydx(2) + a43*dydx(3)
            dydx(1) = (xloc*beta1 + alpha14*dydx(4)) / y1loc
            dydx(5) =  xloc*b5 + a51*dydx(1)*y1loc
            dydx(8) =  xloc*b7 + a72*dydx(2) + a73*dydx(3)
        end if

        dydx(7) = dydx(8)/y4loc - dydx(4)*y(7)
        end associate
    end subroutine rhs_id
end module m_mhd