MODULE MODULE_GAMMA

  !*****************************************************************************
  !** The module 'MODULE_GAMMA' contains gamma (ratio of specific heats at    **
  !** constant pressure and volume) and related variables.                    **
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! gamma = 5/3 for three degrees of freedom
  REAL(KIND=DP), PRIVATE, PARAMETER :: freedom_deg=3._DP ! degrees of freedom
  REAL(KIND=DP), PARAMETER :: GAMMA=(freedom_deg+2._DP)/freedom_deg
  REAL(KIND=DP), PARAMETER :: GAMMA1=1._DP/(GAMMA-1._DP)
  REAL(KIND=DP), PARAMETER :: GAMMA2=GAMMA/(GAMMA-1._DP)
  REAL(KIND=DP), PARAMETER :: GAMMA3=0.5_DP*(GAMMA+1._DP)/(GAMMA-1._DP)

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_GAMMA


MODULE MODULE_DEBUG_JLB

  !*****************************************************************************
  !** The module 'DEBUG_JLB' contains variables for checking ....
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  REAL(KIND=DP) :: tr_1, tr_2, tr_3, tr_4, tr_5, tr_6
  REAL(KIND=DP) :: tr_7, tr_8, tr_9, tr_10, tr_11, tr_12
  REAL(KIND=DP) :: tr_13, tr_14, tr_15, tr_16, tr_17, tr_18, tr_19
  REAL(KIND=DP) :: tr_20, tr_21, tr_22, tr_23, tr_24, tr_25, tr_26
  REAL(KIND=DP) :: tr_27, tr_28, tr_29, tr_30

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_DEBUG_JLB


MODULE MODULE_GRAINS

  !*****************************************************************************
  !** The module 'GRAINS' contains variables related to the dust grains       **
  !*****************************************************************************
! SC May 2006 changes made since last release 16II04
! - lay: added variable R_gr_scale12 = <R> / <R**2>**1/2 
!        added variable layer_thickness 
!(needed in evolution.f90 to compute grain cross section, 
! including contribution of mantles)
!        modified Nsites_grain = 5e4 (to have initial Nlayer compatible
!                                     with mantle density of 1 g cm-3)

  IMPLICIT NONE
  INCLUDE "precision.f90"

  !--- physical characteristics ---
  REAL(KIND=DP), PARAMETER :: Rmin_grain        = 1.0D-6 ! min. radius (cm) (or 4.99D-5)
  REAL(KIND=DP), PARAMETER :: Rmax_grain        = 3.0D-5 ! max. radius (cm) (or 5.00D-5)
  REAL(KIND=DP), PARAMETER :: Rho_GRAIN         = 3.0_DP ! volumic mass (g/cm3) 3.0_DP
  REAL(KIND=DP) :: ratio_GRAIN_gas                       ! mass ratio (grains/gas)
  REAL(KIND=DP) :: R_gr_scale                            ! <R**2> = R_gr_scale * <R**3>**2/3
  REAL(KIND=DP) :: R_gr_scale12                          ! <R> = R_gr_scale12 * <R**2>**1/2
  REAL(KIND=DP) :: Nsites_grain                          ! number of sites per grain 
  REAL(KIND=DP), PARAMETER :: dsites_grain      = 3.4D-8 ! mean distance between sites on grain
  REAL(KIND=DP), PARAMETER :: layer_thickness   = 3.4D-8 ! thickness of one ice layer, in cm

  !--- density, temperature, ...  ---
  REAL(KIND=DP) :: Rsquare_grain ! mean square radius (cm2) calculated in INITIALIZE
  REAL(KIND=DP) :: R_grain       ! mean radius (cm)

  REAL(KIND=DP) :: Dens_grain_init ! initial grain number density (cm-3) calculated in INITIALIZE
  REAL(KIND=DP) :: Dens_grain    ! grain number density (cm-3) calculated in DIFFUN
  REAL(KIND=DP) :: MD_grain! mass of grains per cm3 of gas (g/cm3)
  REAL(KIND=DP) :: Rho_GRAIN_charge ! mass density (g/cm3) of charged grains

  REAL(KIND=DP) :: Tgrain        ! grain temperature (K) read in READ_PARAMETERS
  REAL(KIND=DP) :: Nlayers_grain ! number of layers in grain mantle
  REAL(KIND=DP) :: Teff_grain    ! effective temperature for sputtering reactions

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_GRAINS


MODULE MODULE_PHYS_VAR

  !*****************************************************************************
  !** The module 'MODULE_PHYS_VAR' contains all                               **
  !** the physical variables of the model, stored in two vectors :            **
  !**    * v_variab -> contains all physical variables in the order :         **
  !**                  1) MHD variables (Nv_MHD)                              **
  !**                  2) density of each chemical specy (Nspec)              **
  !**                  3) density of H2 levels (Nlevels_H2)                   **
  !**      (the parameter in parentheses is the number of variables)          **
  !**                                                                         **
  !**    * v_lvariab  = LOG(v_variab)                                         **
  !**                   this vector is transmitted to subroutine DRIVE        **
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  INTEGER(KIND=LONG) :: Nv_MHD     ! number of MHD variables (initialized in MHD)
  INTEGER(KIND=LONG) :: d_v_var ! dimension of the two vectors (calculated in MHD)

  !------------------------------------------------------------------
  ! the two vectors are allocated and initialized in MHD
  ! their dimension is
  ! d_v_var = (Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + NpH2O_lev + NoNH3_lev & 
  !         + NpNH3_lev + NOH_lev + Natype_lev + Netype_lev)
  !------------------------------------------------------------------
  REAL(KIND=DP),DIMENSION(:), SAVE, ALLOCATABLE :: v_variab
  REAL(KIND=DP),DIMENSION(:), SAVE, ALLOCATABLE :: v_lvariab

  !-----------------------------------------------------------------------------
  ! physical variables of the model :
  ! initialized in INITIALIZE and calculated in DIFFUN (except where indicated)
  !-----------------------------------------------------------------------------

  ! neutrals
  REAL(KIND=DP) :: timeN=0.0_DP ! flow time (year) calculated in MHD
  REAL(KIND=DP) :: RhoN         ! mass density (g.cm-3)
  REAL(KIND=DP) :: DensityN     ! number density (cm-3)
  REAL(KIND=DP) :: muN          ! mean mass (g)
  REAL(KIND=DP) :: Tn           ! temperature (K), read in READ_PARAMETERS
  REAL(KIND=DP) :: nH_init      ! initial proton density (cm-3), used in READ_SPECIES
  REAL(KIND=DP) :: Vn           ! velocity (cm/s)
  REAL(KIND=DP) :: dVn          ! neutral velocity gradient
! Used only for J-type shocks:
  REAL(KIND=DP):: XLL           ! characteristic viscous length (cm)
  REAL(KIND=DP) :: grad_V       ! neutral velocity gradient also (with viscosity...)
  REAL(KIND=DP) :: save_dv      ! neutral velocity gradient also (without viscosity...)

  ! positive ions
  REAL(KIND=DP) :: timeI=0.0_DP ! flow time (year) calculated in MHD
  REAL(KIND=DP) :: RhoI         ! mass density (g.cm-3)
  REAL(KIND=DP) :: DensityI     ! number density (cm-3)
  REAL(KIND=DP) :: muI          ! mean mass (g)
  REAL(KIND=DP) :: Ti           ! temperature (K)
  REAL(KIND=DP) :: Vi           ! velocity (cm/s)
  REAL(KIND=DP) :: dVi          ! ions' velocity gradient
 
  ! negative ions
  REAL(KIND=DP) :: RhoNEG       ! mass density (g.cm-3)
  REAL(KIND=DP) :: DensityNEG   ! number density (cm-3)
  REAL(KIND=DP) :: muNEG        ! mean mass (g)

  ! electrons
  REAL(KIND=DP) :: Te             ! temperature (K) (for e- and negative ions)

  ! common
  REAL(KIND=DP) :: nH             ! density of protons nH=n(H)+2n(H2)+n(H+)
  REAL(KIND=DP) :: DeltaV         ! Vi-Vn (cm/s)
  REAL(KIND=DP) :: ABS_DeltaV     ! ABS(Vi-Vn) (cm/s)
  REAL(KIND=DP) :: Vgrad          ! abs. value of the gradient of Vn (km.s-1.cm-1), 
                                  ! calculated in MHD
  REAL(KIND=DP) :: DeltaVmin=1.D3 ! (cm/s) : if ABS_DeltaV < DeltaVmin then stop
  REAL(KIND=DP) :: Vsound         ! sound velocity (cm/s)
  REAL(KIND=DP) :: Vmagnet        ! magnetosonic velocity (cm/s)
  REAL(KIND=DP) :: dmin=1.0D13    ! used when calculating Vgrad

  !--- shock ---
  CHARACTER(len=1)   :: shock_type     ! shock type ('C' or 'J')
  INTEGER(KIND=LONG) :: Nfluids        ! number of fluids (neutrals, ions, e-)
  REAL(KIND=DP)      :: Vs_km          ! shock velocity (km/s), read in READ_PARAMETERS
  REAL(KIND=DP)      :: Vs_cm          ! shock velocity (cm/s), read in READ_PARAMETERS
  REAL(KIND=DP)      :: timeJ          ! shock age (years), read in READ_PARAMETERS
  LOGICAL            :: viscosity      ! .TRUE. if artificial viscosity is to be used
  INTEGER            :: Force_I_C = 1  ! flag : if 1, then impose DensityI = SUM(ions)
  INTEGER            :: Cool_KN = 0    ! flag : if 1, use Kaufman & Neufeld cooling rates
                                       !              for H20 and CO

  ! H2 levels and related variables
  !--- ortho/para ratio            ---
  !--- calculated in COMPUTE_OP_H2 ---
  REAL(KIND=DP)      :: op_H2         ! H2-ortho/para, initialized in READ_PARAMETERS
  INTEGER(KIND=LONG) :: NH2_lev       ! number of levels of H2 included
  INTEGER(KIND=LONG) :: NH2_lines_out ! maximum number of H2 lines in output file
  CHARACTER(LEN=4)   :: H_H2_flag     ! selects H--H2 collision data; options:
                                      !       DRF : Flower et al.
                                      !        MM : Martin & Mandy
                                      !      BOTH : DRF where possible, and MM otherwise
  INTEGER            :: iforH2 = 1    ! flag that determines how much of the H2 binding energy  
                                      ! is released as internal energy when H2 forms on grains  
                                      !   0: 1/3 of 4.4781 eV as internal energy (=> 17249 K) (Allen, 1999)
                                      !   1: proportional to Boltzmann distribution at 17249 K
                                      !   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)
                                      !   3: v = 6, J = 0,1
                                      !   4: fraction = relative populations at t, initialised as H2_lev%density
                                      !                 and changed during integration
  INTEGER            :: ikinH2 = 1    ! flag that determines how much of the H2 binding energy
                                      ! is released as kinetic energy when H2 forms on grains
                                      !   1: 0.5 * (4.4781 - internal)
                                      !   2: Inf(1.4927 eV, 4.4781 - internal)
  ! CH3OH levels and related variables
  INTEGER            :: imeth = 0     ! flag that determines whether the methanol level population
                                      ! densities and related quantities are calculated
                                      !   0: not computed
                                      !   1: are computed
  ! LVG treatment
  INTEGER            :: LVG   = 0     ! flag that determines whether the LVG approximation
                                      ! is used to determine the population densities
                                      ! and related quantities; otherwise, molecular_cooling.f90
                                      !   0: not used
                                      !   1: used

  !--- magnetic field (GAUSS), read in READ_PARAMETERS ---
  REAL(KIND=DP):: Bfield
  REAL(KIND=DP):: Bbeta

  !--- environment, read in READ_PARAMETERS ---
  REAL(KIND=DP) :: RAD  ! multipicative factor for the radiation field 
  REAL(KIND=DP) :: Av   ! extinction (magnitudes)
  REAL(KIND=DP) :: Zeta ! cosmic ray ionization rate (s-1)

  !--- distance (calculated in MHD) ---
  REAL(KIND=DP) :: distance      = 0.0_DP ! distance from start of the shock (cm)
  REAL(KIND=DP) :: dist_step = 0.0_DP     ! step (cm) between two consecutive calls to DRIVE
  !--- calculation step (calc. in MHD) ---
  INTEGER(KIND=LONG) :: counter = 0       ! counts the number of calls to DRIVE

  !-----------------------------------------------------------------
  ! indices enabling variables to be extracted from the two vectors
  ! initialized in MHD
  !-----------------------------------------------------------------
  INTEGER(KIND=LONG) :: iv_Vn, iv_Vi
  INTEGER(KIND=LONG) :: iv_RhoN, iv_RhoI, iv_RhoNEG
  INTEGER(KIND=LONG) :: iv_Tn, iv_Ti, iv_Te
  INTEGER(KIND=LONG) :: iv_DensityN, iv_DensityI
  INTEGER(KIND=LONG) :: iv_gv
  INTEGER(KIND=LONG) :: bv_speci,ev_speci
  INTEGER(KIND=LONG) :: bv_neu,ev_neu
  INTEGER(KIND=LONG) :: bv_ion,ev_ion
  INTEGER(KIND=LONG) :: bv_neg,ev_neg
  INTEGER(KIND=LONG) :: bv_gra,ev_gra
  INTEGER(KIND=LONG) :: bv_cor,ev_cor
  INTEGER(KIND=LONG) :: bv_H2_lev,ev_H2_lev
  INTEGER(KIND=LONG) :: bv_CO_lev,ev_CO_lev
  INTEGER(KIND=LONG) :: bv_SiO_lev,ev_SiO_lev
  INTEGER(KIND=LONG) :: bv_oH2O_lev,ev_oH2O_lev
  INTEGER(KIND=LONG) :: bv_pH2O_lev,ev_pH2O_lev
  INTEGER(KIND=LONG) :: bv_oNH3_lev,ev_oNH3_lev
  INTEGER(KIND=LONG) :: bv_pNH3_lev,ev_pNH3_lev
  INTEGER(KIND=LONG) :: bv_atype_lev,ev_atype_lev
  INTEGER(KIND=LONG) :: bv_etype_lev,ev_etype_lev
  INTEGER(KIND=LONG) :: bv_OH_lev,ev_OH_lev

  !--- used for outputs (used in WRITE_OUTPUTS) ---
  CHARACTER(len=2)  :: speci_out
  CHARACTER(len=7)  :: H2_out
  CHARACTER(len=10) :: line_out
  INTEGER           :: Nstep_w     ! number of steps between two outputs
  INTEGER           :: Nstep_max   ! maximum number of integration steps
  REAL(KIND=DP) :: sum_H2

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  SUBROUTINE READ_PARAMETERS

    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     reads the shock parameters from input file.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     shock_type, Vs_km, Vs_cm, Tn          (MODULE_PHYS_VAR)
    !     op_H2, Zeta, RAD, Av, Bfield          (MODULE_PHYS_VAR)
    !     Nstep_max, Nstep_w                    (MODULE_PHYS_VAR)
    !     Tgrain                                (MODULE_GRAINS)
    !     XLL, Eps_V                            (MODULE_VAR_VODE)
    !     duration_max, length_max              (MODULE_VAR_VODE)
    !---------------------------------------------------------------------------

    USE MODULE_VAR_VODE, ONLY : &
         duration_max, length_max, Eps_V
    USE MODULE_GRAINS, ONLY : Tgrain
    USE MODULE_CONSTANTS, ONLY : YEARsec, screen
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: name_file_input='input/input_mhd.in'
    INTEGER                     :: file_input
    CHARACTER(len=1) :: charact

    !--- opening file ---
    file_input=GET_FILE_NUMBER()
    OPEN(file_input,file=name_file_input,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    !--- shock parameters ---
    READ(file_input,'(A1)')charact    ! comment
    READ(file_input,'(A)') shock_type ! module MODULE_PHYS_VAR
    READ(file_input,*) Nfluids        ! module MODULE_PHYS_VAR
    READ(file_input,*) Bbeta          ! module MODULE_PHYS_VAR
    READ(file_input,*) Vs_km          ! module MODULE_PHYS_VAR
    Vs_cm=1.D5*Vs_km                  ! module MODULE_PHYS_VAR
    READ(file_input,*) DeltaVmin      ! Vi - Vn initial (in cm s-1)
    READ(file_input,*) op_H2          ! module MODULE_H2
    READ(file_input,*) Tn             ! module MODULE_PHYS_VAR
    READ(file_input,*) nH_init        ! module MODULE_PHYS_VAR
    READ(file_input,*) Tgrain         ! module MODULE_GRAINS
    READ(file_input,*) Cool_KN        ! module MODULE_MOLECULAR_COOLING

    ! compute magnetic field in microGauss: (nH in cm-3)
    Bfield = Bbeta * sqrt(nH_init) * 1.0D-6

    !--- environment (module MODULE_PHYS_VAR) ---
    READ(file_input,'(A1)')charact ! comment
    READ(file_input,*) Zeta
    READ(file_input,*) RAD
    READ(file_input,*) Av

    !--- numerical parameters (module MODULE_VAR_VODE) ---
    READ(file_input,'(A1)')charact ! comment
    READ(file_input,*) Nstep_max
    READ(file_input,*) Nstep_w
    READ(file_input,*) NH2_lev
    READ(file_input,*) NH2_lines_out
    READ(file_input,'(A4)') H_H2_flag
    READ(file_input,*) iforH2
    if (iforH2 < 0 .OR. iforH2 > 4) then
      print *, "  iforH2 =", iforH2, "  is not allowed"
      stop
    endif
    READ(file_input,*) ikinH2
    if (ikinH2 /= 1 .AND. ikinH2 /= 2) then
      print *, "  ikinH2 =", ikinH2, "  is not allowed"
      stop
    endif
    READ(file_input,*) imeth
    if (imeth.eq.0) then
      print *, "  methanol level populations are not computed  "
    endif
    READ(file_input,*) LVG
    if (LVG.eq.0) then
      print *, "  LVG approximation is not used; molecular_cooling.f90 instead  "
    endif
    READ(file_input,*) XLL
    READ(file_input,*) Eps_V
    READ(file_input,*) timeJ
    READ(file_input,*) duration_max
    READ(file_input,*) Force_I_C
    length_max=duration_max*Vs_cm*YEARsec

    !--- output specifications (module MODULE_PHYS_VAR)---
    READ(file_input,'(A1)')charact ! comment
    READ(file_input,'(A2)')speci_out
    READ(file_input,'(A7)')H2_out
    READ(file_input,'(A10)')line_out

    !--- closing file ---
    CLOSE(file_input)

    !--- verify some values ---
    if (H_H2_flag /= "DRF" .AND. H_H2_flag /= "MM" .AND. H_H2_flag /= "BOTH") then
      WRITE(screen,'(" ")')
      WRITE(screen,'("  Wrong flag for H--H2 collisions: ",A4)') H_H2_flag
      stop
    endif

    grad_V = 1.0e-2_dp * Vs_cm / XLL

    if (shock_type == "J") then
      viscosity = .TRUE.
      print *, "    grad_V =", grad_V
    endif

  END SUBROUTINE READ_PARAMETERS


END MODULE MODULE_PHYS_VAR


