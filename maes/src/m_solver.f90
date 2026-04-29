!==================================================================================================
!> @fichier m_solver.f90
!> @auteur  Nathan ZIMNIAK
!> @date    10-11-2025
!> @brief   Balayage paramétrique du modèle MHD auto-similaire.
!==================================================================================================
module m_solver
    use m_numerical_init, only: dp, ip
    use m_physical_init,  only: input_parameters
    use m_io,             only: write_result_hdf5
    use m_mhd_solver,     only: compute_solution
    use omp_lib

    implicit none
    private
    public :: run

contains

    !----------------------------------------------------------------------------------------------
    !> Résout le modèle MHD sur une grille de paramètres en parallèle.
    !>
    !> @param[in] input_params_grid tableaux des paramètres d'entrée.
    !----------------------------------------------------------------------------------------------
    subroutine run(input_params_grid)
        type(input_parameters), intent(in) :: input_params_grid(:)

        integer(ip) :: num_points   ! Nombre total de points à explorer.
        integer(ip) :: i            ! Indice de boucle.
        integer(ip) :: current_step ! Nombre de points explorés.
        integer(ip), allocatable :: error_states(:) ! Codes d’erreur par cas.
        integer(ip), allocatable :: int_states(:)   ! États d’intégration finaux.
        integer(ip), allocatable :: crit_states(:)  ! États critiques atteints.
        real(dp) :: start_wtime   ! Temps initial.
        real(dp) :: current_wtime ! Temps courant.
        real(dp) :: elapsed_wtime ! Temps écoulé.
        real(dp) :: progress      ! Fraction de progression globale.

        ! Démarrage du chronomètre.
        start_wtime = omp_get_wtime()

        ! Nombre total de points dans la grille paramétrique.
        num_points = size(input_params_grid, kind=ip)

        ! Allocation mémoire des tableaux.
        allocate(error_states(num_points))
        allocate(int_states(num_points))
        allocate(crit_states(num_points))

        ! Compteur partagé du nombre de points explorés.
        current_step = 0_ip

        ! Distribution de la boucle sur les threads OpenMP.
        !$omp parallel do default(none) &
        !$omp private(i, current_wtime, elapsed_wtime, progress) &
        !$omp shared(input_params_grid, int_states, crit_states, error_states, num_points, current_step, start_wtime)
        do i = 1, num_points
            ! Résolution des équations MHD pour un jeu de paramètres donné (point de la grille).
            call compute_solution(input_params_grid(i), int_states(i), crit_states(i), error_states(i), write_enabled = .false.)

            ! Temps réel écoulé et progression du calcul.
            current_wtime = omp_get_wtime()
            elapsed_wtime = current_wtime - start_wtime
            progress = real(current_step, dp)/real(num_points, dp)

            ! Mise à jour atomique du compteur de progression.
            !$omp atomic
            current_step = current_step + 1_ip
            !$omp end atomic

            ! Affichage de l’itération courante, de la progression, et du temps écoulé.
            !$omp critical
            write(*,'(A,I9,A,F5.1,A,A,F8.1,A)')           &
                  "Iteration ", current_step,             &
                  " | Progress ", 100.0_dp*progress, "%", &
                  " | Elapsed time ", elapsed_wtime, " s"
            !$omp end critical

        end do
        !$omp end parallel do

        ! Écriture des résultats.
        call write_result_hdf5("results.h5", input_params_grid, int_states, crit_states, error_states)

        ! Libération mémoire.
        deallocate(error_states, int_states, crit_states)
    end subroutine run
end module m_solver
