MODULE MODULE_CONSTANTS
  !*****************************************************************************
  !** The module 'MODULE_CONSTANTS' contains mathematical and physical        **
  !** constants and unit conversion factors                                  **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !--- mathematical constants ---
  REAL(KIND=DP), PARAMETER :: pi = 3.141592653589793238_DP
  REAL(KIND=DP), PARAMETER :: Zero = 0.0_DP

  !--- physical constants ---
  REAL(KIND=DP), PARAMETER :: mP = 1.67262D-24 ! mass of the proton (g)
  REAL(KIND=DP), PARAMETER :: me = 9.10939D-28 ! mass of the electron (g)
  REAL(KIND=DP), PARAMETER :: qe = 4.80325D-10 ! charge of the electron (esu)

  REAL(KIND=DP), PARAMETER :: alpha_H  = 6.67D-25 ! polarisability of H (cm-3)
  REAL(KIND=DP), PARAMETER :: alpha_H2 = 7.70D-25 ! polarisability of H2 (cm-3)
  REAL(KIND=DP), PARAMETER :: alpha_He = 2.10D-25 ! polarisability of He (cm-3)

  REAL(KIND=DP), PARAMETER :: bohr = 5.29177D-9 ! bohr radius
  REAL(KIND=DP), PARAMETER :: kB = 1.380658D-16 ! Boltzmann's constant (erg.K-1)
  REAL(KIND=DP), PARAMETER :: R  = 8.314510D7   ! molar gas constant (erg.K-1.mol-1)

  !--- units conversion ---
  REAL(KIND=DP), PARAMETER :: EVerg   = 1.60218D-12 ! 1 eV       =  1.60218D-12 erg
  REAL(KIND=DP), PARAMETER :: amu     = 1.66054D-24 ! 1 amu      =  1.66054D-24 g
  REAL(KIND=DP), PARAMETER :: YEARsec = 3.15569D7   ! 1 year     =  3.15569D7 s
  REAL(KIND=DP), PARAMETER :: kCaleV  = 4.3363D-2   ! 1 kCal/mol =  4.3363D-2 eV
  real(KIND=DP), parameter :: AUcgs   = 6.127D-9    ! conversion between au and cgs
  REAL(KIND=DP), PARAMETER :: parsec  = 3.0857D18   ! 1 pc = 3.0857D18 cm
  REAL(KIND=DP), PARAMETER :: arc_ster= 1/206265D0  ! conversion from arcsec to steradians 

  !--- standard output : UNIX -> 6, MAC -> 9 ---
  INTEGER(KIND=LONG), PARAMETER :: screen = 6


  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_CONSTANTS
