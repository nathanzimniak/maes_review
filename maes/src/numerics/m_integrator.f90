!==================================================================================================
!> @fichier  m_integrator.f90
!> @auteur   Nathan ZIMNIAK
!> @date     10-11-2025
!> @brief    Interface d'intégration des équations MHD à l'aide du solveur LSODA (ODEPACK).
!==================================================================================================
module m_integrator
    use, intrinsic :: ieee_arithmetic
    use m_numerical_init, only: dp, ip, RESISTIVE, IDEAL_LIN, IDEAL_LOG, INTEGRATION_ERROR, NO_ERROR
    use m_physical_init,  only: input_parameters, radial_exponents, midplane_parameters
    use m_mhd,            only: rhs_res, rhs_id
    use odepack_mod,      only: lsoda_class

    implicit none
    private
    public :: integrate

contains

    !----------------------------------------------------------------------------------------------
    !> Intègre le système MHD depuis la position courante jusqu'à une position cible.
    !>
    !> @param[inout] x               Position courante.
    !> @param[inout] y               Vecteur des olutions
    !> @param[in]    xtarget         Position cible.
    !> @param[in]    input_params    Paramètres d'entrée du modèle.
    !> @param[in]    radial_exps     Exposants radiaux.
    !> @param[in]    midplane_params Grandeurs physiques au plan médian.
    !> @param[in]    int_state       État d'intégration courant.
    !----------------------------------------------------------------------------------------------
    subroutine integrate(x, y, x_target, input_params, radial_exps, midplane_params, int_state, error_state)
        real(dp), intent(inout) :: x
        real(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: x_target
        type(input_parameters), intent(in) :: input_params
        type(radial_exponents), intent(in) :: radial_exps
        type(midplane_parameters), intent(in) :: midplane_params
        integer(ip), intent(in) :: int_state
        integer(ip), intent(out) :: error_state

        integer :: istate, itask
        real(dp) :: rtol, atol(1)
        type(lsoda_class) :: ls

        ! Tolérances relatives et absolues utilisées par LSODA.
        rtol   = 1.0e-10_dp
        atol   = 1.0e-10_dp

        ! Configuration du mode d'intégration LSODA.
        istate = 1
        itask  = 1

        ! Initialisation du solveur avec le wrapper local du second membre.
        call ls%initialize(rhs, size(y), iprint=0, istate=istate)
        if (istate < 0) then
            error_state = INTEGRATION_ERROR
            return
        end if

        ! Intégration de x jusqu'à x_target.
        call ls%integrate(y, x, x_target, rtol, atol, itask, istate)
        if ((istate < 0) .or. any(ieee_is_nan(y))) then
            error_state = INTEGRATION_ERROR
        else
            error_state = NO_ERROR
        end if

    contains

        !----------------------------------------------------------------------------------------------
        !> Évalue le second membre demandé par LSODA.
        !>
        !> @param[inout] self Instance interne du solveur LSODA.
        !> @param[in]    neq  Taille du système.
        !> @param[in]    x    Position courante?
        !> @param[in]    y    Vecteur des solutions.
        !> @param[out]   dydx Dérivées des solutions.
        !> @param[out]   ierr Code d’erreur.
        !----------------------------------------------------------------------------------------------
        subroutine rhs(self, neq, x, y, dydx, ierr)
            class(lsoda_class), intent(inout) :: self
            integer, intent(in) :: neq
            real(dp), intent(in) :: x
            real(dp), intent(in) :: y(neq)
            real(dp), intent(out) :: dydx(neq)
            integer, intent(out) :: ierr

            ! Sélection du modèle physique selon le régime d'intégration.
            select case (int_state)
            case (RESISTIVE)
                call rhs_res(x, y, input_params, radial_exps, midplane_params, dydx)
            case (IDEAL_LIN, IDEAL_LOG)
                call rhs_id (x, y, int_state, input_params, radial_exps, midplane_params, dydx)
            end select

            ierr = 0
        end subroutine rhs
    end subroutine integrate
end module m_integrator
