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
    use mpi_f08
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
        integer(ip) :: local_step   ! Nombre local de points explorés par le processus courant.
        integer :: mpi_size ! Nombre total de processus MPI.
        integer :: mpi_rank ! Rang du processus MPI courant.
        integer(ip), allocatable :: error_states(:) ! Codes d’erreur.
        integer(ip), allocatable :: int_states(:)   ! États d’intégration.
        integer(ip), allocatable :: crit_states(:)  ! États critiques.
        real(dp) :: start_wtime   ! Temps initial.
        real(dp) :: current_wtime ! Temps courant.
        real(dp) :: elapsed_wtime ! Temps écoulé.
        real(dp) :: progress      ! Fraction de progression globale.

        ! Synchronisation de tous les processus avant le démarrage du chronomètre.
        call MPI_Barrier(MPI_COMM_WORLD)

        ! Démarrage du chronomètre.
        start_wtime = MPI_Wtime()

        ! Récupération du nombre total de processus MPI et du rang du processus courant.
        call MPI_Comm_size(MPI_COMM_WORLD, mpi_size)
        call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank)

        ! Nombre total de points dans la grille paramétrique.
        num_points = size(input_params_grid, kind=ip)

        ! Allocation mémoire des tableaux.
        allocate(error_states(num_points))
        allocate(int_states(num_points))
        allocate(crit_states(num_points))

        ! Initialisation des tableaux de résultats.
        error_states = 0_ip
        int_states   = 0_ip
        crit_states  = 0_ip

        ! Compteur partagé du nombre de points explorés.
        current_step = 0_ip

        ! Distribution cyclique des points entre processus MPI et parallélisation sur les threads OpenMP.
        !$omp parallel do default(none) &
        !$omp private(i, current_wtime, elapsed_wtime, progress, local_step) &
        !$omp shared(input_params_grid, int_states, crit_states, error_states, num_points, current_step, start_wtime, mpi_rank, mpi_size)
        do i = mpi_rank + 1, num_points, mpi_size
            ! Résolution des équations MHD pour un jeu de paramètres donné (point de la grille).
            call compute_solution(input_params_grid(i), int_states(i), crit_states(i), error_states(i), write_enabled = .false.)

            ! Mise à jour atomique du compteur de progression.
            !$omp atomic capture
            current_step = current_step + 1_ip
            local_step = current_step
            !$omp end atomic

            ! Affichage périodique de la progression estimée depuis le rank 0.
            if (mpi_rank == 0 .and. (mod(local_step, 50_ip/mpi_size) == 0_ip .or. local_step*mpi_size >= num_points)) then
                current_wtime = MPI_Wtime()
                elapsed_wtime = current_wtime - start_wtime
                progress = min(100.0_dp, 100.0_dp*real(local_step*mpi_size, dp)/real(num_points, dp))
            
                !$omp critical
                write(*,'(A,F5.1,A,A,F8.1,A)')    &
                      "Progress ", progress, "%", &
                      " | Elapsed time ", elapsed_wtime, " s"
                !$omp end critical
            end if
        end do
        !$omp end parallel do

        ! Reconstruction des tableaux de résultats complets sur tous les processus.
        call MPI_Allreduce(MPI_IN_PLACE, int_states, num_points, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)
        call MPI_Allreduce(MPI_IN_PLACE, crit_states, num_points, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)
        call MPI_Allreduce(MPI_IN_PLACE, error_states, num_points, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)

        ! Écriture des résultats.
        if (mpi_rank == 0) then
            call write_result_hdf5("results.h5", input_params_grid, int_states, crit_states, error_states)
        end if

        ! Libération mémoire.
        deallocate(error_states, int_states, crit_states)
    end subroutine run
end module m_solver
