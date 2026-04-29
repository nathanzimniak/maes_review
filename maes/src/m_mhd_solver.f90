!==================================================================================================
!> @fichier m_mhd_solver.f90
!> @auteur  Nathan ZIMNIAK
!> @date    10-11-2025
!> @brief   Intégrateur MHD 1D le long de la variable auto-similaire x=z/(epsilon*r), du plan
!>          médian jusqu’à l’infini, avec franchissement des points critiques (SM, A, FM) et
!>          bascule du régime résistif au régime idéal.
!==================================================================================================
module m_mhd_solver
    use m_numerical_init,  only: dp, ip, X_START, MAX_STEP, X_STEP_DEFAULT, X_STEP_FAR, &
                                 X_STEP_INFINITY, X_STEP_CRIT, X_FAR, X_INFINITY,       &
                                 MACH_SM_THR_STEP, MACH_A_THR_STEP, MACH_FM_THR_STEP,   &
                                 RESISTIVE, IDEAL_LIN, IDEAL_LOG, COMPLETED,            &
                                 SUB_SM, SUPER_SM, SUPER_A,                             &
                                 NO_ERROR, LOOP_ERROR,                                  &
                                 MACH_ID_THR, MACH_SM_THR, MACH_A_THR, MACH_FM_THR,     &
                                 N_HIST, integration_history
    use m_physical_init,   only: input_parameters, radial_exponents, midplane_parameters, &
                                 compute_radial_exponents, compute_midplane_parameters,   &
                                 initialize_solutions
    use m_io,              only: generate_id, create_file, write_solutions
    use m_mhd,             only: rhs_res, rhs_id
    use m_integrator,      only: integrate
    use m_mach_numbers,    only: compute_mach_id, compute_mach_sm, compute_mach_a, compute_mach_fm
    use m_critical_points, only: cross_midplane, cross_sm_point, cross_a_point

    implicit none
    private
    public :: compute_solution

contains

    !----------------------------------------------------------------------------------------------
    !> Intègre les équations MHD depuis le plan médian jusqu'à l'infini.
    !>
    !> @param[in]  input_params  paramètres d'entrée.
    !> @param[out] int_state     état d'intégration.
    !> @param[out] crit_state    état critique.
    !> @param[out] error_state   code d'erreur.
    !> @param[in]  write_enabled écriture des résultats.
    !----------------------------------------------------------------------------------------------
    subroutine compute_solution(input_params, int_state, crit_state, error_state, write_enabled)
        type(input_parameters), intent(in) :: input_params
        integer(ip), intent(out) :: int_state
        integer(ip), intent(out) :: crit_state
        integer(ip), intent(out) :: error_state
        logical, intent(in) :: write_enabled

        real(dp) :: x                ! Position courante le long de x.
        real(dp) :: x_step           ! Pas d'intégration courant.
        real(dp) :: x_target         ! Position cible pour l'étape d'intégration.
        real(dp) :: x_switch_lin_log ! Position de transition vers la formulation logarithmique.
        real(dp) :: mach_id          ! Nombre de Mach idéal.
        real(dp) :: mach_sm          ! Nombre de Mach SM.
        real(dp) :: mach_a           ! Nombre de Mach A.
        real(dp) :: mach_fm          ! Nombre de Mach FM.
        real(dp), allocatable :: y(:)    ! Vecteur des solutions.
        real(dp), allocatable :: dydx(:) ! Vecteur des dérivées des solutions.
        integer :: out_unit ! Unité logique pour les sorties fichier.
        character(len=11) :: id ! Identifiant unique de la solution.
        type(radial_exponents)    :: radial_exps      ! Exposants radiaux.
        type(midplane_parameters) :: midplane_params  ! Grandeurs au plan médian.
        type(integration_history) :: integration_hist ! Historique des dernières étapes d'intégration.

        !------------------------------!
        !-- Initialisation numérique --!
        !------------------------------!

        ! États initiaux de l’intégration numérique.
        error_state = NO_ERROR
        int_state = RESISTIVE
        crit_state = SUB_SM

        ! Pas initial et seuil de basculement vers la formulation logarithmique.
        x_step = X_STEP_DEFAULT
        x_switch_lin_log = huge(1.0_dp)

        ! Alias locaux pour lisibilité.
        associate(i_res => integration_hist%i_res, i_id => integration_hist%i_id)

        ! Allocation mémoire des tableaux.
        allocate(y(10))
        allocate(dydx(10))

        ! Création du fichier de sortie.
        if (write_enabled) then
            id = generate_id(input_params)
            call create_file(id, out_unit)
        end if

        !-----------------------------!
        !-- Initialisation physique --!
        !-----------------------------!
        
        ! Calcul des exposants radiaux (be, alpha0 etc.).
        call compute_radial_exponents(input_params, radial_exps)
        
        ! Calcul des grandeurs au plan médian (mu, p etc.).
        call compute_midplane_parameters(input_params, radial_exps, midplane_params)

        ! Initialisation des solutions (calcul de y en x = 0).
        x = 0.0_dp
        call initialize_solutions(input_params, y, dydx)
        
        ! Écriture des solutions.
        if (write_enabled)  call write_solutions(out_unit, x, y, dydx, int_state)
        
        ! Traitement du point critique x = 0 par développement limité.
        x = X_START
        call cross_midplane(x, y, input_params, radial_exps, midplane_params)
        
        ! Mise à jour des dérivées au nouveau point.
        call rhs_res(x, y, input_params, radial_exps, midplane_params, dydx)
        
        ! Écriture des solutions.
        if (write_enabled) call write_solutions(out_unit, x, y, dydx, int_state)

        !-------------------------------------!
        !-- Intégration du système résistif --!
        !-------------------------------------!
        
        do while ((int_state == RESISTIVE) .and. (error_state == NO_ERROR))
            ! Incrémente le nombre d'itérations en régime résistif.
            i_res = i_res + 1_ip

            ! Intégration du système d’EDO de x vers x_target.
            x_target = x + x_step
            call integrate(x, y, x_target, input_params, radial_exps, midplane_params, int_state, error_state)
            if (error_state /= NO_ERROR) exit

            ! Mise à jour des dérivées au nouveau point.
            call rhs_res(x, y, input_params, radial_exps, midplane_params, dydx)

            ! Écriture des solutions.
            if (write_enabled) call write_solutions(out_unit, x, y, dydx, int_state)

            ! Test du critère de transition vers le régime idéal.
            call compute_mach_id(x, y, input_params, radial_exps, midplane_params, integration_hist, mach_id)
            if (mach_id >= MACH_ID_THR) then
                ! Changement de régime.
                y    = [y([1,2,3,4,5,6,7,10])]
                dydx = [dydx([1,2,3,4,5,6,7,10])]
                int_state = IDEAL_LIN
                allocate(integration_hist%x(N_HIST))
                allocate(integration_hist%y(size(y),N_HIST))
            end if

            ! Arrêt de sécurité en cas de non-convergence de la phase résistive.
            if (i_res >= MAX_STEP) error_state = LOOP_ERROR
        enddo

        !---------------------------------------------!
        !-- Intégration du système idéal (linéaire) --!
        !---------------------------------------------!

        do while ((int_state /= COMPLETED) .and. (error_state == NO_ERROR))
            ! Incrémente le nombre d'itérations en régime idéal.
            i_id = i_id + 1_ip

            ! Sauvegarde du point courant pour les extrapolations critiques.
            integration_hist%x(i_id)    = x
            integration_hist%y(:, i_id) = y

            ! Ajustement du pas d'intégration selon le prochain point critique à franchir.
            select case (crit_state)
            case (SUB_SM)
                call compute_mach_sm(x, y, input_params, radial_exps, midplane_params, int_state, mach_sm)
                x_step = merge(X_STEP_CRIT, X_STEP_DEFAULT, mach_sm >= MACH_SM_THR_STEP)

            case (SUPER_SM)
                call compute_mach_a(x, y, input_params, radial_exps, midplane_params, int_state, mach_a)
                x_step = merge(X_STEP_CRIT, X_STEP_DEFAULT, mach_a >= MACH_A_THR_STEP)

            case (SUPER_A)
                call compute_mach_fm(x, y, input_params, radial_exps, midplane_params, int_state, mach_fm)
                select case (int_state)
                case (IDEAL_LIN)
                    x_step = merge(X_STEP_CRIT, X_STEP_DEFAULT, mach_fm >= MACH_FM_THR_STEP)
                case (IDEAL_LOG)
                    if (mach_fm >= MACH_FM_THR_STEP) then; x_step = X_STEP_CRIT
                    else if (x < X_FAR)              then; x_step = X_STEP_DEFAULT
                    else if (x < X_INFINITY)         then; x_step = X_STEP_FAR
                    else                                 ; x_step = X_STEP_INFINITY
                    end if
                end select
            end select

            ! Intégration du système d’EDO de x vers x_target.
            x_target = x + x_step
            call integrate(x, y, x_target, input_params, radial_exps, midplane_params, int_state, error_state)
            if (error_state /= NO_ERROR) exit

            ! Mise à jour des dérivées au nouveau point.
            call rhs_id(x, y, int_state, input_params, radial_exps, midplane_params, dydx)

            ! Écriture des solutions.
            if (write_enabled) call write_solutions(out_unit, x, y, dydx, int_state)

            ! Traitement des points critiques x = x_SM, x_A, x_FM par extrapolation.
            select case (crit_state)
            case (SUB_SM)
                ! Test du critère de transition vers le régime SM.
                call compute_mach_sm(x, y, input_params, radial_exps, midplane_params, int_state, mach_sm)
                if ((mach_sm >= MACH_SM_THR) .or. (dydx(5) <= 0.0_dp) .or. (dydx(4) >= 0.0_dp) .or. (-dydx(6) <= 0)) then
                    !write(*,*) "-------------------------"
                    !write(*,*) "Début du saut en x = ", x
                    !write(*,*) "dydx(4) = ", dydx(4)
                    !write(*,*) "dydx(5) = ", dydx(5)
                    !write(*,*) "Br = ", -dydx(6)
                    !write(*,*) "Bphi = ", y(1)
                    !write(*,*) "Msm2 = ", mach_sm**2
                    !write(*,*) "-------------------------"
                    ! Franchissement du point SM.
                    call cross_sm_point(x, y, dydx, input_params, radial_exps, midplane_params, integration_hist, int_state, mach_sm, error_state)
                    if (error_state /= NO_ERROR) exit
                    crit_state = SUPER_SM
                    ! Écriture des solutions.
                    if (write_enabled) call write_solutions(out_unit, x, y, dydx, int_state)
                end if

            case (SUPER_SM)
                ! Test du critère de transition vers le régime A.
                call compute_mach_a(x, y, input_params, radial_exps, midplane_params, int_state, mach_a)
                if ((mach_a >= MACH_A_THR) .or. (dydx(5) <= 0.0_dp) .or. (dydx(4) >= 0.0_dp) .or. (-dydx(6) <= 0)) then
                    ! Franchissement du point A.
                    call cross_a_point(x, y, dydx, input_params, radial_exps, midplane_params, integration_hist, int_state, mach_a, error_state)
                    if (error_state /= NO_ERROR) exit
                    crit_state = SUPER_A
                    ! Écriture des solutions.
                    if (write_enabled) call write_solutions(out_unit, x, y, dydx, int_state)
                    ! Prépare la transition vers la résolution logarithmique.
                    x_switch_lin_log = x + 10*X_STEP_DEFAULT
                end if

            case (SUPER_A)
                ! Test du critère de transition vers le régime FM.
                call compute_mach_fm(x, y, input_params, radial_exps, midplane_params, int_state, mach_fm)
                if (mach_fm >= MACH_FM_THR) then
                    int_state = COMPLETED
                endif
                ! Vérification du basculement en résolution logarithmique.
                if ((int_state /= IDEAL_LOG) .and. (x >= x_switch_lin_log)) then
                    int_state = IDEAL_LOG
                    x    = log(x)
                    y(1) = log(-y(1))
                    y(6) = log(y(6))
                endif
            end select

            ! Arrêt de sécurité en cas de non-convergence de la phase idéale.
            if (i_id >= MAX_STEP) error_state = LOOP_ERROR
        enddo

        ! Libération mémoire.
        deallocate(y)
        deallocate(dydx)
        if (allocated(integration_hist%x)) deallocate(integration_hist%x)
        if (allocated(integration_hist%y)) deallocate(integration_hist%y)
        if (write_enabled) close(out_unit)
        end associate
    end subroutine compute_solution
end module m_mhd_solver