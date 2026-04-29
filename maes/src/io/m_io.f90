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
    public :: generate_id, create_file, write_solutions, write_result_hdf5

contains

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
        integer :: h5err
        integer(HID_T) :: id_file, id_dspace, id_dset, id_plist, id_aspace, id_attr, id_atype
        integer(HSIZE_T) :: strlen_h5
        integer(HSIZE_T), dimension(1) :: adims
        integer(HSIZE_T), dimension(2) :: dims
        integer(HSIZE_T), dimension(2) :: chunk_dims
        real(dp), allocatable :: buffer(:,:)
        character(len=256) :: column_names

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

        ! Création du type chaîne de caractères adapté à l'attribut.
        call h5tcopy_f(H5T_FORTRAN_S1, id_atype, h5err)
        strlen_h5 = int(len_trim(column_names), kind=HSIZE_T)
        call h5tset_size_f(id_atype, strlen_h5, h5err)
        call h5tset_strpad_f(id_atype, H5T_STR_NULLTERM_F, h5err)

        ! Création et écriture de l'attribut "column_names".
        call h5acreate_f(id_dset, "column_names", id_atype, id_aspace, id_attr, h5err)
        call h5awrite_f(id_attr, id_atype, column_names, adims, h5err)

        ! Fermeture des objets HDF5 associés à l'attribut.
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
end module m_io