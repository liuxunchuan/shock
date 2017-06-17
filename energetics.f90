MODULE MODULE_ENERGETICS
  !*****************************************************************************
  !** The module 'MODULE_ENERGETICS' contains variables and one               **
  !** subroutine related to the fluxes of energy, mass and momentum.          **
  !** All the variables are calculated here.                                  **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !-----------
  ! mass flux
  !-----------
  REAL(KIND=DP) :: Mass_flux      ! mass flux (g/s/cm2)
  REAL(KIND=DP) :: Mass_flux_init ! initial mass flux (g/s/cm2)
  REAL(KIND=DP) :: Mass_cons      ! relative mass flux conservation

  !---------------
  ! momentum flux
  !---------------
  REAL(KIND=DP) :: Momentum_flux_kin  ! kinetic momentum flux (erg/cm3)
  REAL(KIND=DP) :: Momentum_flux_the  ! thermal momentum flux (erg/cm3)
  REAL(KIND=DP) :: Momentum_flux_mag  ! magnetic momentum flux (erg/cm3)
  REAL(KIND=DP) :: Momentum_flux_vis  ! viscous momentum flux (erg/cm3)
  REAL(KIND=DP) :: Momentum_flux      ! total momentum flux (erg/cm3)
  REAL(KIND=DP) :: Momentum_flux_init ! initial total momentum flux (erg/cm3)
  REAL(KIND=DP) :: Momentum_cons      ! relative momentum flux conservation

  !-------------
  ! energy flux
  !-------------
  REAL(KIND=DP) :: Energy_flux_kin      ! kinetic Energy_flux (erg/s/cm2)
  REAL(KIND=DP) :: Energy_flux_the      ! thermal Energy_flux (erg/s/cm2)
  REAL(KIND=DP) :: Energy_flux_int      ! H2 internal Energy_flux (erg/s/cm2)
  REAL(KIND=DP) :: Energy_flux_mag      ! magnetic Energy_flux (erg/s/cm2)
  REAL(KIND=DP) :: Energy_flux_vis      ! viscous Energy_flux (erg/s/cm2)
  REAL(KIND=DP) :: Energy_flux          ! total Energy_flux (erg/s/cm2)
  REAL(KIND=DP) :: Energy_flux_init     ! initial total Energy_flux (erg/s/cm2)
  REAL(KIND=DP) :: Energy_cons          ! relative energy flux conservation
  REAL(KIND=DP) :: Heating_old = 0.0_DP ! heating (erg/s/cm3) before last call to DRIVE
  REAL(KIND=DP) :: Energy_gain = 0.0_DP ! energy gain (erg/s/cm2) of the gas due to heating

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  SUBROUTINE ENERGETIC_FLUXES
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD (twice)
    ! purpose :
    !    Calculate the fluxes of mass(g/s/cm2), momentum (erg/cm3) and energy
    !    (erg/s/cm2). When this subroutine is called for the first time, it
    !    save the initial fluxes for comparison purpose.
    !    It also checks the conservation of those fluxes.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    all the variables of this module are calculated here
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR
    USE MODULE_EVOLUTION, ONLY : Bn, Bi, Bneg, H2_int_energy
    USE MODULE_CONSTANTS, ONLY : kB, pi, mP, EVerg
    USE MODULE_CHEMICAL_SPECIES, ONLY : speci, ind_H, ind_H2, ind_He
    USE MODULE_DEBUG_JLB

    IMPLICIT NONE
    INTEGER(KIND=LONG), SAVE :: Ncall = 0 ! counts the number of call to this routine

    Ncall = Ncall + 1


    !--- mass (g/s/cm2) ---
    Mass_flux = &
         RhoN   * Vn + &
         RhoI   * Vi + &
         RhoNEG * Vi

    !--- momentum (erg/cm3) ---
    Momentum_flux_kin = &
         RhoN   * Vn**2._DP + &
         RhoI   * Vi**2._DP + &
         RhoNEG * Vi**2._DP
    Momentum_flux_the = &
         RhoN   * kB * Tn / muN + &
         RhoI   * kB * Ti / muI + &
         RhoNEG * kB * Te / muNEG
    Momentum_flux_mag = &
         (Bfield * Vs_cm/Vi)**2._DP / (8._DP*pi)
    !ccccccccccccccccccccccccccccccccccccccccc
    ! calculate the viscosity term for J-type
    !ccccccccccccccccccccccccccccccccccccccccc
    if ((shock_type == "J") .AND. (viscosity)) then
      Momentum_flux_vis = &
         (RhoN + RhoI + RhoNEG) * (XLL * grad_V)**2.0_DP
    else
      Momentum_flux_vis = 0.0_dp
    endif
    Momentum_flux = &
         Momentum_flux_kin + &
         Momentum_flux_the + &
         Momentum_flux_mag + &
         Momentum_flux_vis

    !--- energy flux (erg/s/cm2) ---
    Energy_flux_kin = &
         0.5_DP * RhoN   * Vn**3._DP + &
         0.5_DP * RhoI   * Vi**3._DP + &
         0.5_DP * RhoNEG * Vi**3._DP
    Energy_flux_the = &
         2.5_DP * RhoN   * Vn * kB * Tn / muN + &
         2.5_DP * RhoI   * Vi * kB * Ti / muI + &
         2.5_DP * RhoNEG * Vi * kB * Te / muNEG
    Energy_flux_int = &
         H2_int_energy   * Vn * kB
    Energy_flux_mag = &
         Bfield**2._DP / (4._DP*pi) * Vs_cm**2._DP / Vi
    !ccccccccccccccccccccccccccccccccccccccccc
    ! calculate the viscosity term for J-type
    !ccccccccccccccccccccccccccccccccccccccccc
    if ((shock_type == "J") .AND. (viscosity)) then
      Energy_flux_vis = &
         (RhoN + RhoI + RhoNEG) * Vn * (XLL * grad_V)**2.0_DP
    else
      Energy_flux_vis = 0.0_dp
    endif
    Energy_flux = &
         Energy_flux_kin + &
         Energy_flux_the + &
         Energy_flux_int + &
         Energy_flux_mag + &
         Energy_flux_vis

    ! if this is the first call of this routine, then save the fluxes
    IF (Ncall <= 2) THEN
       Mass_flux_init     = Mass_flux
       Momentum_flux_init = Momentum_flux
       Energy_flux_init   = Energy_flux
    END IF

    !--------------------
    ! check conservation
    !--------------------

    !--- relative mass conservation ---
    Mass_cons = (Mass_flux_init - Mass_flux) / Mass_flux_init

    !--- relative momentum conservation ---
    Momentum_cons = (Momentum_flux_init - Momentum_flux) / Momentum_flux_init

    !--- relative energy conservation ---
    Energy_gain = Energy_gain + &
         0.5_DP * dist_step * (Heating_old + (Bn + Bi + Bneg))
    Heating_old = Bn + Bi + Bneg


! 23 oct 2001 - petit bug de normalisation (GPdF + JLB)

    Energy_cons = (Energy_flux_init - Energy_flux + Energy_gain) &
                / Energy_flux_init

  END SUBROUTINE ENERGETIC_FLUXES

END MODULE MODULE_ENERGETICS
