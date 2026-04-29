!==================================================================================================
!> @fichier main.f90
!> @auteur  Nathan ZIMNIAK
!> @date    28-04-2025
!> @brief   Programme principal de balayage paramétrique du modèle MHD auto-similaire : génération
!>          de la grille (xi, mu, p) et exécution du solveur sur l’ensemble des points.
!==================================================================================================
program main
    use m_solver,         only: run
    use m_numerical_init, only: dp, ip
    use m_physical_init,  only: input_parameters
    
    implicit none

    integer(ip) :: n_xi, n_mu, n_p
    integer(ip) :: i, i_xi, i_mu, i_p
    real(dp) :: xi_min, xi_max
    real(dp) :: mu_min, mu_max
    real(dp) :: p_min, p_max
    real(dp), allocatable :: xi_list(:), mu_list(:), p_list(:)
    type(input_parameters), allocatable :: input_parameters_grid(:)

    ! Nombre de points dans chaque direction de la grille paramétrique.
    n_xi = 20_ip
    n_mu = 20_ip
    n_p  = 20_ip

    ! Bornes des paramètres explorés.
    xi_min = 0.000100_dp
    xi_max = 0.000200_dp
    mu_min = 0.400000_dp
    mu_max = 0.600000_dp
    p_min = 4.000_dp
    p_max = 5.000_dp

    ! Allocation mémoire des tableaux.
    allocate(xi_list(n_xi))
    allocate(mu_list(n_mu))
    allocate(p_list(n_p))
    allocate(input_parameters_grid(n_xi*n_mu*n_p))

    ! Construction des listes définissant la grille explorée.
    do i = 1, n_xi
        xi_list(i) = 10.0_dp**(log10(xi_min) + (log10(xi_max) - log10(xi_min))*real(i-1,dp)/real(n_xi-1,dp))
        xi_list(i) = real(nint(xi_list(i)*1.0e6_dp), dp)/1.0e6_dp
    end do
    do i = 1, n_mu
        mu_list(i) = 10.0_dp**(log10(mu_min) + (log10(mu_max) - log10(mu_min))*real(i-1,dp)/real(n_mu-1,dp))
        mu_list(i) = real(nint(mu_list(i)*1.0e6_dp), dp)/1.0e6_dp
    end do
    do i = 1, n_p
        p_list(i) = p_min + (p_max - p_min)*real(i-1,dp)/real(n_p-1,dp)
        p_list(i) = real(nint(p_list(i)*1.0e3_dp), dp)/1.0e3_dp
    end do
    
    ! Remplissage de la grille explorée
    i = 0_ip
    do i_p = 1, n_p
        do i_mu = 1, n_mu
            do i_xi = 1, n_xi
                i = i + 1_ip
                input_parameters_grid(i) = input_parameters(xi     = xi_list(i_xi), &
                                                            ep     = 0.100_dp,      &
                                                            alpham = 1.0_dp,        &
                                                            chim   = 1.0_dp,        &
                                                            Pm     = 1.0_dp,        &
                                                            alphap = 0.0_dp,        &
                                                            mu     = mu_list(i_mu), &
                                                            p      = p_list(i_p))
            end do
        end do
    end do

    ! Lancement du balayage paramétrique.
    call run(input_parameters_grid)

    ! Libération mémoire.
    deallocate(xi_list, mu_list, p_list, input_parameters_grid)

end program main
