!==================================================================================================
!> @fichier m_io.f90
!> @auteur  Nathan ZIMNIAK
!> @date    10-11-2025
!> @brief   Gestion des entrées/sorties.
!==================================================================================================
module m_io
    use m_numerical_init, only: dp, ip, RESISTIVE, IDEAL_LIN, IDEAL_LOG
    use m_physical_init,  only: input_parameters
    use hdf5

    implicit none
    private
    public :: generate_id, create_file, write_solutions, write_result_hdf5, read_setup

contains

    !----------------------------------------------------------------------------------------------
    !> Lit le fichier de configuration du balayage paramétrique.
    !>
    !>
    !> @param[in]  filename    Nom du fichier de configuration.
    !> @param[out] xi_list     Liste des valeurs de xi.
    !> @param[out] mu_list     Liste des valeurs de mu.
    !> @param[out] p_list      Liste des valeurs de p.
    !> @param[out] ep_list     Liste des valeurs de ep.
    !> @param[out] alpham_list Liste des valeurs de alpham.
    !> @param[out] chim_list   Liste des valeurs de chim.
    !> @param[out] Pm_list     Liste des valeurs de Pm.
    !> @param[out] alphap_list Liste des valeurs de alphap.
    !----------------------------------------------------------------------------------------------
    subroutine read_setup(filename, xi_list, mu_list, p_list, ep_list, alpham_list, chim_list, Pm_list, alphap_list)
        character(*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: xi_list(:)
        real(dp), allocatable, intent(out) :: mu_list(:)
        real(dp), allocatable, intent(out) :: p_list(:)
        real(dp), allocatable, intent(out) :: ep_list(:)
        real(dp), allocatable, intent(out) :: alpham_list(:)
        real(dp), allocatable, intent(out) :: chim_list(:)
        real(dp), allocatable, intent(out) :: Pm_list(:)
        real(dp), allocatable, intent(out) :: alphap_list(:)

        integer(ip) :: n_xi, n_mu, n_p
        integer(ip) :: i
        real(dp) :: xi_min, xi_max
        real(dp) :: mu_min, mu_max
        real(dp) :: p_min, p_max

        character(len=256) :: line, key, value
        integer :: unit, ios, pos, j, n_values

        ! Ouverture du fichier de configuration en lecture.
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "Erreur : impossible d'ouvrir le fichier setup."

        ! Lecture du fichier ligne par ligne.
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            ! Suppression des espaces initiaux.
            line = adjustl(line)

            ! Ignorer les lignes vides et les commentaires.
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "#") cycle

            ! Recherche du séparateur clé/valeur.
            pos = index(line, "=")
            if (pos == 0) cycle

            ! Extraction de la clé et de la valeur associée
            key   = adjustl(trim(line(:pos-1)))
            value = adjustl(trim(line(pos+1:)))

            ! Lecture de la valeur selon la clé rencontrée.
            select case (trim(key))
            case ("n_xi"); read(value, *) n_xi
            case ("n_mu"); read(value, *) n_mu
            case ("n_p"); read(value, *) n_p
            case ("xi_min"); read(value, *) xi_min
            case ("xi_max"); read(value, *) xi_max
            case ("mu_min"); read(value, *) mu_min
            case ("mu_max"); read(value, *) mu_max
            case ("p_min"); read(value, *) p_min
            case ("p_max"); read(value, *) p_max
            case ("ep")
                n_values = 1
                do j = 1, len_trim(value)
                    if (value(j:j) == ",") n_values = n_values + 1
                end do
                allocate(ep_list(n_values))
                read(value, *) ep_list
            case ("alpham")
                n_values = 1
                do j = 1, len_trim(value)
                    if (value(j:j) == ",") n_values = n_values + 1
                end do
                allocate(alpham_list(n_values))
                read(value, *) alpham_list
            case ("chim")
                n_values = 1
                do j = 1, len_trim(value)
                    if (value(j:j) == ",") n_values = n_values + 1
                end do
                allocate(chim_list(n_values))
                read(value, *) chim_list
            case ("Pm")
                n_values = 1
                do j = 1, len_trim(value)
                    if (value(j:j) == ",") n_values = n_values + 1
                end do
                allocate(Pm_list(n_values))
                read(value, *) Pm_list
            case ("alphap")
                n_values = 1
                do j = 1, len_trim(value)
                    if (value(j:j) == ",") n_values = n_values + 1
                end do
                allocate(alphap_list(n_values))
                read(value, *) alphap_list
            case default
                error stop "Erreur : clé inconnue dans le fichier setup."
            end select
        end do

        close(unit)
    
        ! Allocation mémoire des listes générées à partir des bornes lues.
        allocate(xi_list(n_xi))
        allocate(mu_list(n_mu))
        allocate(p_list(n_p))

        ! Construction des listes définissant la grille explorée.
        do i = 1, n_xi
            xi_list(i) = 10.0_dp**(log10(xi_min) + (log10(xi_max) - log10(xi_min))*real(i-1, dp)/real(n_xi-1, dp))
            xi_list(i) = real(nint(xi_list(i)*1.0e6_dp), dp)/1.0e6_dp
        end do
        do i = 1, n_mu
            mu_list(i) = 10.0_dp**(log10(mu_min) + (log10(mu_max) - log10(mu_min))*real(i-1, dp)/real(n_mu-1, dp))
            mu_list(i) = real(nint(mu_list(i)*1.0e6_dp), dp)/1.0e6_dp
        end do
        do i = 1, n_p
            p_list(i) = p_min + (p_max - p_min)*real(i-1, dp)/real(n_p-1, dp)
            p_list(i) = real(nint(p_list(i)*1.0e3_dp), dp)/1.0e3_dp
        end do
    end subroutine read_setup

    !----------------------------------------------------------------------------------------------
    !> Écrit le résumé du balayage paramétrique dans un fichier HDF5.
    !>
    !> @param[in] filename          Nom du fichier HDF5.
    !> @param[in] input_params_grid Grille des paramètres d'entrée.
    !> @param[in] int_states        États d'intégration finaux.
    !> @param[in] crit_states       États critiques atteints.
    !> @param[in] error_states      Codes d'erreur.
    !----------------------------------------------------------------------------------------------
    subroutine write_result_hdf5(filename, input_params_grid, int_states, crit_states, error_states)
        character(*), intent(in) :: filename
        type(input_parameters), intent(in) :: input_params_grid(:)
        integer(ip), intent(in) :: int_states(:)
        integer(ip), intent(in) :: crit_states(:)
        integer(ip), intent(in) :: error_states(:)

        integer(ip) :: i
        integer(ip) :: num_points
        integer(ip) :: n_xi, n_mu, n_p
        integer :: h5err
        integer(HID_T) :: id_file, id_dspace, id_dset, id_plist, id_aspace, id_attr, id_atype
        integer(HSIZE_T) :: strlen_h5
        integer(HSIZE_T), dimension(1) :: adims
        integer(HSIZE_T), dimension(2) :: dims
        integer(HSIZE_T), dimension(2) :: chunk_dims
        real(dp), allocatable :: buffer(:,:)
        character(len=256) :: column_names
        character(len=4096) :: setup_summary
        character(len=512) :: line

        ! Nombre total de points dans la grille paramétrique.
        num_points = size(input_params_grid, kind=ip)

        ! Allocation mémoire du tableau.
        allocate(buffer(num_points, 11))

        ! Construction du tableau de sortie.
        do i = 1, num_points
            buffer(i, 1) = input_params_grid(i)%xi
            buffer(i, 2) = input_params_grid(i)%ep
            buffer(i, 3) = input_params_grid(i)%alpham
            buffer(i, 4) = input_params_grid(i)%chim
            buffer(i, 5) = input_params_grid(i)%Pm
            buffer(i, 6) = input_params_grid(i)%alphap
            buffer(i, 7) = input_params_grid(i)%mu
            buffer(i, 8) = input_params_grid(i)%p
            buffer(i, 9) = real(int_states(i), dp)
            buffer(i,10) = real(crit_states(i), dp)
            buffer(i,11) = real(error_states(i), dp)
        end do

        ! Reconstruction des nombres de valeurs distinctes.
        n_xi = count_unique_real(buffer(:,1))
        n_mu = count_unique_real(buffer(:,7))
        n_p  = count_unique_real(buffer(:,8))

        ! Construction du résumé textuel du setup.
        setup_summary = ""

        write(line,'(A)') "========================================================"
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A,I0)') "Nombre total de points explorés : ", num_points
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A)') "--------------------------------------------------------"
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A,F0.6,A,F0.6,A,I0,A)') &
            "xi     ∈ [", minval(buffer(:,1)), ", ", maxval(buffer(:,1)), "] ; ", n_xi, " points en échelle log"
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A,F0.6,A,F0.6,A,I0,A)') &
            "mu     ∈ [", minval(buffer(:,7)), ", ", maxval(buffer(:,7)), "] ; ", n_mu, " points en échelle log"
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A,F0.3,A,F0.3,A,I0,A)') &
            "p      ∈ [", minval(buffer(:,8)), ", ", maxval(buffer(:,8)), "] ; ", n_p, " points en échelle lin"
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A,A)') "ep     = ", trim(real_list_to_string(unique_real(buffer(:,2)), 3))
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A,A)') "alpham = ", trim(real_list_to_string(unique_real(buffer(:,3)), 1))
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')
            
        write(line,'(A,A)') "chim   = ", trim(real_list_to_string(unique_real(buffer(:,4)), 1))
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')
            
        write(line,'(A,A)') "Pm     = ", trim(real_list_to_string(unique_real(buffer(:,5)), 1))
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')
            
        write(line,'(A,A)') "alphap = ", trim(real_list_to_string(unique_real(buffer(:,6)), 1))
        setup_summary = trim(setup_summary)//trim(line)//new_line('a')

        write(line,'(A)') "========================================================"
        setup_summary = trim(setup_summary)//trim(line)

        ! Initialisation de la bibliothèque HDF5.
        call h5open_f(h5err)

        ! Création du fichier HDF5 de sortie.
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, id_file, h5err)

        ! Définition de la taille du dataset "results".
        dims = (/ int(num_points, kind=HSIZE_T), int(11, kind=HSIZE_T) /)
        call h5screate_simple_f(2, dims, id_dspace, h5err)

        ! Création de la liste de propriétés du dataset.
        call h5pcreate_f(H5P_DATASET_CREATE_F, id_plist, h5err)

        ! Activation du stockage par chunks, nécessaire pour la compression.
        chunk_dims = (/ int(min(num_points, 1024_ip), kind=HSIZE_T), int(11, kind=HSIZE_T) /)
        call h5pset_chunk_f(id_plist, 2, chunk_dims, h5err)

        ! Activation de la compression gzip/deflate.
        call h5pset_deflate_f(id_plist, 4, h5err)

        ! Création du dataset principal contenant les résultats.
        call h5dcreate_f(id_file, "results", H5T_NATIVE_DOUBLE, id_dspace, &
                         id_dset, h5err, dcpl_id=id_plist)

        ! Écriture du tableau de résultats.
        call h5dwrite_f(id_dset, H5T_NATIVE_DOUBLE, buffer, dims, h5err)

        ! Attribut texte décrivant l'ordre des colonnes du dataset.
        column_names = "xi,ep,alpham,chim,Pm,alphap,mu,p,int_state,crit_state,error_state"

        ! Création du dataspace scalaire de l'attribut.
        adims(1) = 1_HSIZE_T
        call h5screate_simple_f(1, adims, id_aspace, h5err)

        ! Création et écriture de l'attribut "column_names".
        call h5tcopy_f(H5T_FORTRAN_S1, id_atype, h5err)
        strlen_h5 = int(len_trim(column_names), kind=HSIZE_T)
        call h5tset_size_f(id_atype, strlen_h5, h5err)
        call h5tset_strpad_f(id_atype, H5T_STR_NULLTERM_F, h5err)

        call h5acreate_f(id_dset, "column_names", id_atype, id_aspace, id_attr, h5err)
        call h5awrite_f(id_attr, id_atype, column_names, adims, h5err)

        call h5aclose_f(id_attr, h5err)
        call h5tclose_f(id_atype, h5err)

        ! Création et écriture de l'attribut global "setup_summary".
        call h5tcopy_f(H5T_FORTRAN_S1, id_atype, h5err)
        strlen_h5 = int(len_trim(setup_summary), kind=HSIZE_T)
        call h5tset_size_f(id_atype, strlen_h5, h5err)
        call h5tset_strpad_f(id_atype, H5T_STR_NULLTERM_F, h5err)

        call h5acreate_f(id_file, "setup_summary", id_atype, id_aspace, id_attr, h5err)
        call h5awrite_f(id_attr, id_atype, setup_summary, adims, h5err)

        call h5aclose_f(id_attr, h5err)
        call h5sclose_f(id_aspace, h5err)
        call h5tclose_f(id_atype, h5err)

        ! Fermeture des objets HDF5 principaux.
        call h5pclose_f(id_plist, h5err)
        call h5dclose_f(id_dset, h5err)
        call h5sclose_f(id_dspace, h5err)
        call h5fclose_f(id_file, h5err)
        call h5close_f(h5err)

        ! Libération mémoire.
        deallocate(buffer)
    end subroutine write_result_hdf5

    !----------------------------------------------------------------------------------------------
    !> Écrit les solutions et leurs dérivées dans un fichier.
    !>
    !> @param[in] out_unit  Unité logique du fichier de sortie.
    !> @param[in] x         Position courante.
    !> @param[in] y         Vecteur des solutions.
    !> @param[in] dydx      Vecteur des dérivées.
    !> @param[in] int_state État d'intégration courant.
    !----------------------------------------------------------------------------------------------
    subroutine write_solutions(out_unit, x, y, dydx, int_state)
        integer, intent(in) :: out_unit
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(in) :: dydx(:)
        integer(ip), intent(in) :: int_state

        real(dp) :: x_local
        real(dp) :: y1, y2, y3, y4, y5, y6, y7, y8
        real(dp) :: dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8

        select case (int_state)
        case (RESISTIVE)
            ! Écriture depuis le vecteur résistif complet.
            x_local = x
            y1 = y(1)
            y2 = y(2)
            y3 = y(3)
            y4 = exp(y(4))
            y5 = y(5)
            y6 = y(6)
            y7 = y(7)
            y8 = y(10)
            dydx1 = dydx(1)
            dydx2 = dydx(2)
            dydx3 = dydx(3)
            dydx4 = dydx(4)*exp(y(4))
            dydx5 = dydx(5)
            dydx6 = dydx(6)
            dydx7 = dydx(7)
            dydx8 = dydx(10)

        case (IDEAL_LIN)
            ! Écriture depuis le vecteur idéal en coordonnée linéaire.
            x_local = x
            y1 = y(1)
            y2 = y(2)
            y3 = y(3)
            y4 = exp(y(4))
            y5 = y(5)
            y6 = y(6)
            y7 = y(7)
            y8 = y(8)
            dydx1 = dydx(1)
            dydx2 = dydx(2)
            dydx3 = dydx(3)
            dydx4 = dydx(4)*exp(y(4))
            dydx5 = dydx(5)
            dydx6 = dydx(6)
            dydx7 = dydx(7)
            dydx8 = dydx(8)

        case (IDEAL_LOG)
            ! Reconversion des variables logarithmiques avant écriture.
            x_local = exp(x)
            y1 = -exp(y(1))
            y2 = y(2)
            y3 = y(3)
            y4 = exp(y(4))
            y5 = y(5)
            y6 = exp(y(6))
            y7 = y(7)
            y8 = y(8)
            dydx1 = dydx(1)*(-exp(y(1)))/exp(x)
            dydx2 = dydx(2)/exp(x)
            dydx3 = dydx(3)/exp(x)
            dydx4 = dydx(4)*exp(y(4))/exp(x)
            dydx5 = dydx(5)/exp(x)
            dydx6 = dydx(6)*exp(y(6))/exp(x)
            dydx7 = dydx(7)/exp(x)
            dydx8 = dydx(8)/exp(x)
        end select

        ! Écriture d'une ligne contenant les solutions et leurs dérivées au point courant.
        write(out_unit,'(*(ES30.15E3,:,","))') x_local, y1, y2, y3, y4, y5, y6, y7, y8, dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7, dydx8
    end subroutine write_solutions

    !----------------------------------------------------------------------------------------------
    !> Crée un fichier CSV pour l’écriture des solutions et de leurs dérivées.
    !>
    !> @param[in]  name     Nom du fichier.
    !> @param[out] out_unit Unité logique associée au fichier ouvert.
    !----------------------------------------------------------------------------------------------
    subroutine create_file(name, out_unit)
        character(*), intent(in) :: name
        integer, intent(out) :: out_unit

        ! Ouverture du fichier CSV en mode écriture (remplacement si existant).
        open(newunit=out_unit, file=name//'.csv', status='replace', action='write')

        ! Écriture de la ligne d'en-tête.
        write(out_unit,'(A)') "x, Bphi, vr, vz, rho, Omega, Psi, T, P, Bphi', vr', vz', rho', Omega', Psi', T', P'"
    end subroutine create_file

    !----------------------------------------------------------------------------------------------
    !> Génère un identifiant unique et déterministe pour une solution.
    !>
    !> @param[in] input_param Paramètres d'entrée du modèle.
    !> @return    id          Identifiant chaîne de longueur fixe.
    !----------------------------------------------------------------------------------------------
    pure function generate_id(input_param) result(id)
        type(input_parameters), intent(in) :: input_param

        integer(ip) :: ix, ie
        character(len=11) :: id

        ! Discrétisation des paramètres.
        ix = nint(input_param%xi * 1.0e5_dp)
        ie = nint(input_param%ep * 1.0e3_dp)

        ! Encodage en chaîne.
        write(id, '(I6.6,I5.5)') ix, ie
    end function generate_id

    !----------------------------------------------------------------------------------------------
    !> Compte le nombre de valeurs réelles distinctes dans un tableau.
    !----------------------------------------------------------------------------------------------
    pure function count_unique_real(values) result(n_unique)
        real(dp), intent(in) :: values(:)
    
        integer(ip) :: n_unique
        integer(ip) :: i, j
        logical :: found
    
        n_unique = 0_ip
    
        do i = 1, size(values, kind=ip)
            found = .false.
        
            do j = 1, i - 1_ip
                if (values(i) == values(j)) then
                    found = .true.
                    exit
                end if
            end do
        
            if (.not. found) n_unique = n_unique + 1_ip
        end do
    end function count_unique_real
    
    !----------------------------------------------------------------------------------------------
    !> Extrait les valeurs réelles distinctes d'un tableau en conservant leur ordre d'apparition.
    !----------------------------------------------------------------------------------------------
    function unique_real(values) result(unique_values)
        real(dp), intent(in) :: values(:)
        real(dp), allocatable :: unique_values(:)
    
        integer(ip) :: i, j, n_unique
        logical :: found
    
        n_unique = count_unique_real(values)
        allocate(unique_values(n_unique))
    
        n_unique = 0_ip
    
        do i = 1, size(values, kind=ip)
            found = .false.
        
            do j = 1, n_unique
                if (values(i) == unique_values(j)) then
                    found = .true.
                    exit
                end if
            end do
        
            if (.not. found) then
                n_unique = n_unique + 1_ip
                unique_values(n_unique) = values(i)
            end if
        end do
    end function unique_real

    !----------------------------------------------------------------------------------------------
    !> Convertit une liste de réels en chaîne au format [x, y, z].
    !----------------------------------------------------------------------------------------------
    function real_list_to_string(values, decimals) result(str)
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: decimals
        character(len=1024) :: str

        integer(ip) :: i
        character(len=64) :: fmt
        character(len=64) :: value_str

        write(fmt,'(A,I0,A)') "(F0.", decimals, ")"

        str = "["

        do i = 1, size(values, kind=ip)
            write(value_str, fmt) values(i)

            if (i > 1_ip) then
                str = trim(str)//", "//trim(value_str)
            else
                str = trim(str)//trim(value_str)
            end if
        end do

        str = trim(str)//"]"
    end function real_list_to_string
end module m_io