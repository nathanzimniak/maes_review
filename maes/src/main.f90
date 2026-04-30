!==================================================================================================
!> @fichier main.f90
!> @auteur  Nathan ZIMNIAK
!> @date    28-04-2025
!> @brief   Programme principal de balayage paramétrique du modèle MHD auto-similaire : génération
!>          de la grille (xi,mu,p) et exécution du solveur sur l’ensemble des points.
!==================================================================================================
program main
    use m_solver,         only: run
    use m_numerical_init, only: dp, ip
    use m_physical_init,  only: input_parameters
    use m_io,             only: read_setup
    use mpi_f08
    
    implicit none

    integer(ip) :: n_xi, n_mu, n_p, n_ep, n_alpham, n_chim, n_Pm, n_alphap, n_pts
    integer(ip) :: i, i_xi, i_mu, i_p, i_ep, i_alpham, i_chim, i_Pm, i_alphap
    integer :: mpi_rank
    real(dp), allocatable :: xi_list(:), mu_list(:), p_list(:)
    real(dp), allocatable :: ep_list(:), alpham_list(:), chim_list(:), Pm_list(:), alphap_list(:)
    type(input_parameters), allocatable :: input_parameters_grid(:)

    ! Initialisation de l'environnement d'exécution MPI.
    call MPI_Init()

    ! Récupération du rang du processus courant.
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank)

    ! Lecture du fichier setup et génération des grilles paramétriques pour le balayage.
    call read_setup("setup.in", xi_list, mu_list, p_list, ep_list, alpham_list, chim_list, Pm_list, alphap_list)

    ! Récupération du nombre de valeurs dans chaque direction de la grille.
    n_xi = size(xi_list, kind=ip)
    n_mu = size(mu_list, kind=ip)
    n_p = size(p_list, kind=ip)
    n_ep = size(ep_list, kind=ip)
    n_alpham = size(alpham_list, kind=ip)
    n_chim = size(chim_list, kind=ip)
    n_Pm = size(Pm_list, kind=ip)
    n_alphap = size(alphap_list, kind=ip)

    ! Nombre total de points explorés.
    n_pts = n_xi*n_mu*n_p*n_ep*n_alpham*n_chim*n_Pm*n_alphap
    
    ! Affichage du résumé du balayage paramétrique.
    if (mpi_rank == 0) then
        write(*,*) ""
        write(*,'(A)') "========================================================"
        write(*,'(A,I0)')                 "Nombre total de points explorés : ", n_pts
        write(*,'(A)') "--------------------------------------------------------"
        write(*,'(A,F0.6,A,F0.6,A,I0,A)')      "xi     ∈ [", xi_list(1), ", ", xi_list(n_xi), "] ; " , n_xi, " points en échelle log"
        write(*,'(A,F0.6,A,F0.6,A,I0,A)')      "mu     ∈ [", mu_list(1), ", ", mu_list(n_mu), "] ; " , n_mu, " points en échelle log"
        write(*,'(A,F0.3,A,F0.3,A,I0,A)')      "p      ∈ [", p_list(1),  ", ", p_list(n_p),   "] ; ",  n_p,  " points en échelle lin"
        write(*,'(A,*(F0.3,A))', advance='no') "ep     = [", (ep_list(i),     merge(", ", char(0)//char(0), i < n_ep),     i=1,n_ep);     write(*,'(A)') "]"
        write(*,'(A,*(F0.1,A))', advance='no') "alpham = [", (alpham_list(i), merge(", ", char(0)//char(0), i < n_alpham), i=1,n_alpham); write(*,'(A)') "]"
        write(*,'(A,*(F0.1,A))', advance='no') "chim   = [", (chim_list(i),   merge(", ", char(0)//char(0), i < n_chim),   i=1,n_chim);   write(*,'(A)') "]"
        write(*,'(A,*(F0.1,A))', advance='no') "Pm     = [", (Pm_list(i),     merge(", ", char(0)//char(0), i < n_Pm),     i=1,n_Pm);     write(*,'(A)') "]"
        write(*,'(A,*(F0.1,A))', advance='no') "alphap = [", (alphap_list(i), merge(", ", char(0)//char(0), i < n_alphap), i=1,n_alphap); write(*,'(A)') "]"
        write(*,'(A)') "========================================================"
        write(*,*) ""
    end if

    ! Allocation mémoire de la grille explorée.
    allocate(input_parameters_grid(n_pts))
    
    ! Remplissage de la grille explorée.
    i = 0_ip
    do i_ep = 1, n_ep
        do i_alpham = 1, n_alpham
            do i_chim = 1, n_chim
                do i_Pm = 1, n_Pm
                    do i_alphap = 1, n_alphap
                        do i_xi = 1, n_xi
                            do i_p = 1, n_p
                                do i_mu = 1, n_mu
                                    i = i + 1_ip
                                    input_parameters_grid(i) = input_parameters(xi     = xi_list(i_xi),         &
                                                                                ep     = ep_list(i_ep),         &
                                                                                alpham = alpham_list(i_alpham), &
                                                                                chim   = chim_list(i_chim),     &
                                                                                Pm     = Pm_list(i_Pm),         &
                                                                                alphap = alphap_list(i_alphap), &
                                                                                mu     = mu_list(i_mu),         &
                                                                                p      = p_list(i_p))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    ! Lancement du balayage paramétrique.
    call run(input_parameters_grid)

    ! Libération mémoire.
    deallocate(xi_list, mu_list, p_list, ep_list, alpham_list, chim_list, Pm_list, alphap_list, input_parameters_grid)

    ! Finalisation de l'environnement MPI et libération des ressources associées.
    call MPI_Finalize()
end program main
