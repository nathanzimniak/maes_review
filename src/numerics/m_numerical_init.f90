!==================================================================================================
!> @fichier m_numerical_init.f90
!> @auteur  Nathan ZIMNIAK
!> @date    10-11-2025
!> @brief   Paramètres numériques globaux et constantes de contrôle utilisées par le solveur MHD.
!==================================================================================================
module m_numerical_init
    use, intrinsic :: iso_fortran_env, only : real64, int32

    implicit none
    private

    ! Précision numérique.
    integer, parameter, public :: dp = real64 ! Réels double précision.
    integer, parameter, public :: ip = int32  ! Entiers 32 bits.

    ! Bornes d'intégration.
    real(dp), parameter, public :: X_START  = 0.0001_dp ! Début de l'intégration.
    integer(ip), parameter, public :: MAX_STEP = 5000_ip ! Nombre maximal de pas d'intégration.

    ! Pas d'intégration.
    real(dp), parameter, public :: X_STEP_DEFAULT  = 0.010_dp ! Pas d'intégration par défaut.
    real(dp), parameter, public :: X_STEP_FAR      = 0.100_dp ! Pas d'intégration à partir de X_FAR.
    real(dp), parameter, public :: X_STEP_INFINITY = 1.000_dp ! Pas d'intégration à partir de X_INF.
    real(dp), parameter, public :: X_STEP_CRIT     = 0.001_dp ! Pas d'intégration à partir de MACH_*_THR_STEP.

    ! Seuils en position pour l'adaptation du pas d'intégration.
    real(dp), parameter, public :: X_FAR = 10_dp
    real(dp), parameter, public :: X_INFINITY = 100_dp
    
    ! Seuils en nombre de Mach pour la réduction du pas d’intégration.
    real(dp), parameter, public :: MACH_SM_THR_STEP = 0.700_dp
    real(dp), parameter, public :: MACH_A_THR_STEP  = 0.700_dp
    real(dp), parameter, public :: MACH_FM_THR_STEP = 0.400_dp

    ! État d'intégration du système.
    integer(ip), parameter, public :: RESISTIVE = 0_ip
    integer(ip), parameter, public :: IDEAL_LIN = 1_ip
    integer(ip), parameter, public :: IDEAL_LOG = 2_ip
    integer(ip), parameter, public :: COMPLETED = 3_ip

    ! État de vitesse du système.
    integer(ip), parameter, public :: SUB_SM   = 0_ip
    integer(ip), parameter, public :: SUPER_SM = 1_ip
    integer(ip), parameter, public :: SUPER_A  = 2_ip
    integer(ip), parameter, public :: SUPER_FM = 3_ip

    ! État d'erreur du système.
    integer(ip), parameter, public :: NO_ERROR          = 0_ip
    integer(ip), parameter, public :: INTEGRATION_ERROR = 1_ip
    integer(ip), parameter, public :: LOOP_ERROR        = 2_ip
    integer(ip), parameter, public :: JUMP_ERROR        = 3_ip

    ! Seuils de passage des points critiques.
    real(dp), parameter, public :: MACH_ID_THR = 0.997_dp
    real(dp), parameter, public :: MACH_SM_THR = 1.000_dp
    real(dp), parameter, public :: MACH_A_THR  = 0.950_dp
    real(dp), parameter, public :: MACH_FM_THR = 0.900_dp

    ! Historique d'intégration.
    integer(ip), parameter, public :: N_HIST = 5000_ip ! Profondeur de l'historique.
    type, public :: integration_history
        integer(ip) :: i_res = 0_ip ! Compteur du nombre d'intégration en MHD résistive.
        integer(ip) :: i_id  = 0_ip ! Compteur du nombre d'intégration en MHD idéale.
        real(dp) :: x_last = 0.0_dp         ! Dernière valeur de la position.
        real(dp) :: mach_id_p_last = 0.0_dp ! Dernière valeur du nombre de Mach "idéal" poloidal.
        real(dp) :: mach_id_t_last = 0.0_dp ! Dernière valeur du nombre de Mach "idéal" toroidal.
        real(dp) :: mach_id_v_last = 0.0_dp ! Dernière valeur du nombre de Mach "idéal" visqueux.
        real(dp), allocatable :: x(:)   ! Dernières valeurs de la position.
        real(dp), allocatable :: y(:,:) ! Dernières valeurs des solutions.
    end type integration_history
end module m_numerical_init