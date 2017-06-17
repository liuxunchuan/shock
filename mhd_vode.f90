PROGRAM MHD

  !#############################################################################
  !#############################################################################
  !##                                                                         ##
  !##                  Main program for the MHD shock model                   ##
  !##                  ------------------------------------                   ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !## updates :                                                               ##
  !## --------                                                                ##
  !##    * Oct 2000 : David Wilgenbus                                         ##
  !##        - translation into fortran 90 and revision of mhd C-type code    ##
  !##          of G. Pineau des Forets and D. Flower                          ##
  !##                                                                         ##
  !##    * Nov 2000 : David Wilgenbus                                         ##
  !##        - correction of grain compression                                ##
  !##        - H2O and CO cooling taken from Neufeld & Kaufman (1993)         ##
  !##          (* 2011-2012 : David Flower substituted LVG calculations of    ##
  !##          level populations and line intensities, from which the         ##
  !##          corresponding cooling rates are derived)                       ##
  !##    * May/June 2001 : Jacques Le Bourlot                                 ##
  !##        - Replacement of GEAR by DVODE (Pierre Hily-Blant)               ##
  !##        - General revision (with simplifications)                        ##
  !##        - Correction of some minor bugs                                  ##
  !##        - Suppression of variables 11 et 12                              ##
  !##          (grain mass and radius)                                        ##
  !##                                                                         ##
  !##   * winter/spring 2001/2002 : DRF - GPdF - JLB                          ##
  !##        - Correction of several bugs                                     ##
  !##        - J-type shocks                                                  ##
  !##                                                                         ##
  !#############################################################################
  !#############################################################################

!-------------------------------------------------
!--- call of modules (variables, subroutines)  ---
!-------------------------------------------------

  USE MODULE_TOOLS,             ONLY : CREATE_ARCHIVE_FILE, GET_FILE_NUMBER, JACO
  USE MODULE_GRAINS,            ONLY : MD_grain, R_grain
  USE MODULE_CONSTANTS,         ONLY : YEARsec, pi, Zero, parsec, mP, kB, EVerg
  USE MODULE_INITIALIZE,        ONLY : INITIALIZE
  USE MODULE_CHEMICAL_SPECIES
  USE MODULE_H2,                ONLY : H2_lev, H2_lines, &
                                       NH2_lines, Index_VJ_H2, op_LTE
  USE MODULE_PHYS_VAR
  USE MODULE_VAR_VODE,          ONLY : duration_max, length_max, Eps_V, T0_V, H0_V, &
                                       Tout_V, MF_V, Itask_V, Istate_V, Hdone, Hnext, &
                                       liw_V, lrw_V, itol_V, iopt_V, atol_V, rtol_V, &
                                       iwork_V, rwork_V, ipar_V, rpar_V
  USE mvode
  USE MODULE_OUTPUTS
  USE MODULE_EVOLUTION
  USE MODULE_ENERGETICS,        ONLY : ENERGETIC_FLUXES
  USE MODULE_MOLECULAR_COOLING, ONLY : err_cool, f_err_cool, &
                                       n_err_cool
  USE MODULE_DEBUG_JLB
  USE MODULE_LINE_EXCIT

!--------------------------------
!--- declaration of variables ---
!--------------------------------

  IMPLICIT NONE
    INTEGER :: NCO_lev, NSiO_lev, NoH2O_lev, NpH2O_lev, NoNH3_lev, NpNH3_lev, NOH_lev &
               , Natype_lev, Netype_lev
! The values of the following PARAMETERS must be the same as those of the 
! PARAMETER (NLEVCO) in the corresponding LVG modules (CO.f, SiO.f, ... )
! See also evolution.f90 and outputs.f90
    PARAMETER (NCO_lev=41)
    PARAMETER (NSiO_lev=41)
    PARAMETER (NoH2O_lev=45)
    PARAMETER (NpH2O_lev=45)
    PARAMETER (NoNH3_lev=17)
    PARAMETER (NpNH3_lev=24)
    PARAMETER (NOH_lev=20)
    PARAMETER (Natype_lev=256)
    PARAMETER (Netype_lev=256)
! INCLUDE "precision.f90"

  INTEGER (KIND=LONG) :: i
  LOGICAL             :: its_OK = .FALSE. ! integration step is OK
  LOGICAL             :: again = .TRUE.   ! continue integration while it is .TRUE.
  LOGICAL             :: writ_step        ! if .TRUE. then write to output files
  LOGICAL             :: time_up = .FALSE.! transition to J-shock when time_up = .TRUE.
  LOGICAL             :: preJ=.TRUE.      ! if .TRUE. C-shock prior to J-discontinuity
  CHARACTER (LEN=80)  :: message_end      ! message written at the end of integration

! values saved before call to DRIVE
  REAL (KIND=DP)      :: distance_old = 0.0_DP, Tn_old = 0.0_DP
  REAL (KIND=DP)      :: Vn_old = 0.0_DP, Vi_old = 0.0_DP

  REAL (KIND=DP), DIMENSION (:), ALLOCATABLE :: dlv_var
  INTEGER :: iflag, jj

    REAL(KIND=DP) :: dens_CO_lev(1:NCO_lev),TAU_CO(NCO_lev,NCO_lev),WWT_CO, &
                     TVGRAD_CO(1:NCO_lev), &
                     TVGRAD_CO_old(1:NCO_lev),CO_lines(1:NCO_lev), &
                     dens_CO_lev_old(NCO_lev),CO_cdens(1:NCO_lev)
    REAL(KIND=DP) :: dens_SiO_lev(1:NSiO_lev),TAU_SiO(NSiO_lev,NSiO_lev),WWT_SiO, &
                     TVGRAD_SiO(1:NSiO_lev), &
                     TVGRAD_SiO_old(1:NSiO_lev),SiO_lines(1:NSiO_lev), &
                     dens_SiO_lev_old(NSiO_lev),SiO_cdens(1:NSiO_lev)
    REAL(KIND=DP) :: dens_oH2O_lev(1:NoH2O_lev),freq_oH2O(NoH2O_lev,NoH2O_lev), &
                     tau_oH2O(NoH2O_lev,NoH2O_lev),WWT_oH2O, &
                     TVGRAD_oH2O(NoH2O_lev,NoH2O_lev), &
                     TVGRAD_oH2O_old(NoH2O_lev,NoH2O_lev), &
                     oH2O_lines(NoH2O_lev,NoH2O_lev), &
                     dens_oH2O_lev_old(NoH2O_lev),oH2O_cdens(1:NoH2O_lev)
    REAL(KIND=DP) :: dens_pH2O_lev(1:NpH2O_lev),freq_pH2O(NpH2O_lev,NpH2O_lev), &
                     tau_pH2O(NpH2O_lev,NpH2O_lev),WWT_pH2O, &
                     TVGRAD_pH2O(NpH2O_lev,NpH2O_lev), &
                     TVGRAD_pH2O_old(NpH2O_lev,NpH2O_lev), &
                     pH2O_lines(NpH2O_lev,NpH2O_lev), &
                     dens_pH2O_lev_old(NpH2O_lev),pH2O_cdens(1:NpH2O_lev)
    REAL(KIND=DP) :: dens_oNH3_lev(1:NoNH3_lev),freq_oNH3(NoNH3_lev,NoNH3_lev), &
                     tau_oNH3(NoNH3_lev,NoNH3_lev),WWT_oNH3, &
                     TVGRAD_oNH3(NoNH3_lev,NoNH3_lev), &
                     TVGRAD_oNH3_old(NoNH3_lev,NoNH3_lev), &
                     oNH3_lines(NoNH3_lev,NoNH3_lev), &
                     dens_oNH3_lev_old(NoNH3_lev),oNH3_cdens(1:NoNH3_lev)
    REAL(KIND=DP) :: dens_pNH3_lev(1:NpNH3_lev),freq_pNH3(NpNH3_lev,NpNH3_lev), &
                     tau_pNH3(NpNH3_lev,NpNH3_lev),WWT_pNH3, &
                     TVGRAD_pNH3(NpNH3_lev,NpNH3_lev), &
                     TVGRAD_pNH3_old(NpNH3_lev,NpNH3_lev), &
                     pNH3_lines(NpNH3_lev,NpNH3_lev), &
                     dens_pNH3_lev_old(NpNH3_lev),pNH3_cdens(1:NpNH3_lev)
    REAL(KIND=DP) :: dens_OH_lev(1:NOH_lev),freq_OH(NOH_lev,NOH_lev), &
                     tau_OH(NOH_lev,NOH_lev),WWT_OH, &
                     TVGRAD_OH(NOH_lev,NOH_lev), &
                     TVGRAD_OH_old(NOH_lev,NOH_lev), &
                     OH_lines(NOH_lev,NOH_lev), &
                     dens_OH_lev_old(NOH_lev),OH_cdens(1:NOH_lev)
    REAL(KIND=DP) :: dens_atype_lev(1:Natype_lev),freq_atype(Natype_lev,Natype_lev) &
                    ,tau_atype(Natype_lev,Natype_lev),WWT_atype,TVGRAD_atype(Natype_lev,Natype_lev) &
                    ,TVGRAD_atype_old(Natype_lev,Natype_lev) &
                    ,atype_lines(Natype_lev,Natype_lev) &
                    ,dens_atype_lev_old(Natype_lev) &
                    ,atype_cdens(Natype_lev)
    REAL(KIND=DP) :: dens_etype_lev(1:Netype_lev),freq_etype(Netype_lev,Netype_lev) &
                    ,tau_etype(Netype_lev,Netype_lev),WWT_etype,TVGRAD_etype(Netype_lev,Netype_lev) &
                    ,TVGRAD_etype_old(Netype_lev,Netype_lev) &
                    ,etype_lines(Netype_lev,Netype_lev) &
                    ,dens_etype_lev_old(Netype_lev) &
                    ,etype_cdens(Netype_lev)
      COMMON /FCTCO/ dens_CO_lev,TAU_CO,WWT_CO,TVGRAD_CO,CO_lines,CO_cdens
      COMMON /FCTSiO/ dens_SiO_lev,TAU_SiO,WWT_SiO,TVGRAD_SiO,SiO_lines,SiO_cdens
      COMMON /FCToH2O/ dens_oH2O_lev,freq_oH2O,tau_oH2O,WWT_oH2O,TVGRAD_oH2O,oH2O_lines,oH2O_cdens
      COMMON /FCTpH2O/ dens_pH2O_lev,freq_pH2O,tau_pH2O,WWT_pH2O,TVGRAD_pH2O,pH2O_lines,pH2O_cdens
      COMMON /FCToNH3/ dens_oNH3_lev,freq_oNH3,tau_oNH3,WWT_oNH3,TVGRAD_oNH3,oNH3_lines,oNH3_cdens
      COMMON /FCTpNH3/ dens_pNH3_lev,freq_pNH3,tau_pNH3,WWT_pNH3,TVGRAD_pNH3,pNH3_lines,pNH3_cdens
      COMMON /FCTOH/ dens_OH_lev,freq_OH,tau_OH,WWT_OH,TVGRAD_OH,OH_lines,OH_cdens
      COMMON /FCTa_CH3OH/ dens_atype_lev,freq_atype,tau_atype &
                         ,WWT_atype,TVGRAD_atype,atype_lines,atype_cdens
      COMMON /FCTe_CH3OH/ dens_etype_lev,freq_etype,tau_etype &
                         ,WWT_etype,TVGRAD_etype,etype_lines,etype_cdens
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the CO lines
! and the rate of cooling per unit volume by CO
!------------------------------------------------------------------
  CO_lines(1:NCO_lev) =  0.0_DP
  CO_cdens(1:NCO_lev) =  0.0_DP
  WWT_CO =  0.0_DP
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the SiO lines
! and the rate of cooling per unit volume by SiO
!------------------------------------------------------------------
  SiO_lines(1:NSiO_lev) =  0.0_DP
  SiO_cdens(1:NSiO_lev) =  0.0_DP
  WWT_SiO =  0.0_DP
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the o-H2O lines
! and the rate of cooling per unit volume by o-H2O
!------------------------------------------------------------------
  oH2O_lines(1:NoH2O_lev,1:NoH2O_lev) =  0.0_DP
  oH2O_cdens(1:NoH2O_lev) =  0.0_DP
  WWT_oH2O =  0.0_DP
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the p-H2O lines
! and the rate of cooling per unit volume by p-H2O
!------------------------------------------------------------------
  pH2O_lines(1:NpH2O_lev,1:NpH2O_lev) =  0.0_DP
  pH2O_cdens(1:NpH2O_lev) =  0.0_DP
  WWT_pH2O =  0.0_DP
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the o-NH3 lines
! and the rate of cooling per unit volume by o-NH3
!------------------------------------------------------------------
  oNH3_lines(1:NoNH3_lev,1:NoNH3_lev) =  0.0_DP
  oNH3_cdens(1:NoNH3_lev) =  0.0_DP
  WWT_oNH3 =  0.0_DP
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the p-NH3 lines
! and the rate of cooling per unit volume by p-NH3
!------------------------------------------------------------------
  pNH3_lines(1:NpNH3_lev,1:NpNH3_lev) =  0.0_DP
  pNH3_cdens(1:NpNH3_lev) =  0.0_DP
  WWT_pNH3 =  0.0_DP
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the OH lines
! and the rate of cooling per unit volume by OH
!------------------------------------------------------------------
  OH_lines(1:NOH_lev,1:NOH_lev) =  0.0_DP
  OH_cdens(1:NOH_lev) =  0.0_DP
  WWT_OH =  0.0_DP
!------------------------------------------------------------------
! initialize the integrated intensities (K km s-1) of the A-type lines
! the column densities of the A-type levels,
! and the rate of cooling per unit volume by A-type
!------------------------------------------------------------------
  atype_lines(1:Natype_lev,1:Natype_lev) =  0.0_DP
  atype_cdens(1:Natype_lev) =  0.0_DP
  WWT_atype =  0.0_DP
!----------------------------
! initialize the integrated intensities (K km s-1) of the E-type lines,
! the column densities of the E-type levels,
! and the rate of cooling per unit volume by E-type
!------------------------------------------------------------------
  etype_lines(1:Netype_lev,1:Netype_lev) =  0.0_DP
  etype_cdens(1:Netype_lev) =  0.0_DP
  WWT_etype =  0.0_DP
! message written to screen
!----------------------------
  WRITE(*,'(80("-"))')
  WRITE(*,'(20X,"MHD SHOCK")')
  WRITE(*,'(20X,"---------")')


!-------------------------------------------------------
! initializations :
! -----------------
!    * shock parameters
!    * chemical species (+ check)
!    * chemical reactions (+ check)
!    * molecule H2 (levels, collision rates, lines)
!    * molecule SiO (collision rates, Einstein coefficients)
!-------------------------------------------------------
  CALL INITIALIZE

!---------------------------------
! create the archive-utility file
!---------------------------------
  CALL CREATE_ARCHIVE_FILE(shock_type,op_LTE,op_H2,nH,Vs_km)

!----------------------------------
! open the file for error messages
!----------------------------------
  f_err_cool = GET_FILE_NUMBER()
  OPEN(f_err_cool, file=n_err_cool, STATUS='REPLACE', &
       access='SEQUENTIAL', form='FORMATTED', action='WRITE')

!---------------------------------------------------
! compute elemental abundances (H, C, N, ...)
! result in elements%ab, saved in elements%ab_init
!---------------------------------------------------
  CALL ELEMENTAL_ABUNDANCES
  elements(:)%ab_init = elements(:)%ab


!---------------------------
! put into vector form
!---------------------------
!--- allocation of the vectors containing physical variables ---
!--- dimension = 1:d_v_var                                   ---
!--- deallocation at the end of MHD                          ---
! if (shock_type == "J") then
    Nv_MHD = 11
! else
!   Nv_MHD = 10
! endif
  d_v_var = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + NpH2O_lev + NoNH3_lev & 
            + NpNH3_lev + NOH_lev + Natype_lev + Netype_lev
            
  ALLOCATE(v_variab  (1:d_v_var)) ! physical variables (Y)
  ALLOCATE(v_lvariab (1:d_v_var)) ! logarithm of the variables
  ALLOCATE(v_dvariab (1:d_v_var)) ! value of Y*dY at the last call to DIFFUN
  ALLOCATE(v_l_var   (1:d_v_var)) ! Y at the last call to DIFFUN
  ALLOCATE(v_l_der   (1:d_v_var)) ! logarithm of the variables
  ALLOCATE(dlv_var   (1:d_v_var)) ! right-hand side
  ALLOCATE(YN        (0:Nspec))   ! change in number density (cm-3.s-1); calculated in CHEMISTRY

 999  continue
!--- initialization ---
  v_variab(:)  = Zero
  v_lvariab(:) = Zero
  v_dvariab(:) = Zero
  v_l_var(:)   = Zero
  v_l_der(:)   = Zero

!--- MHD variables, number = Nv_MHD ---
  i = 1  ; v_variab(i) = Vn       ; iv_Vn = i       ! neutrals (cm/s)
  i = 2  ; v_variab(i) = Vi       ; iv_Vi = i       ! ions (cm/s)
  i = 3  ; v_variab(i) = RhoN     ; iv_RhoN = i     ! neutrals (g.cm-3)
  i = 4  ; v_variab(i) = RhoI     ; iv_RhoI = i     ! ions (g.cm-3)
  i = 5  ; v_variab(i) = Tn       ; iv_Tn = i       ! neutrals (K)
  i = 6  ; v_variab(i) = Ti       ; iv_Ti = i       ! ions (K)
  i = 7  ; v_variab(i) = Te       ; iv_Te = i       ! electrons (K)
  i = 8  ; v_variab(i) = DensityN ; iv_DensityN = i ! neutrals (cm-3)
  i = 9  ; v_variab(i) = DensityI ; iv_DensityI = i ! ions (cm-3)
  i = 10 ; v_variab(i) = RhoNEG   ; iv_RhoNEG = i   ! negative ions (g.cm-3)
! if (shock_type == "J") then
    i = 11 ; v_variab(i) = grad_V   ; iv_gv     = i   ! - neutral velocity gradient
! endif

!--- density of species (cm-3), number = Nspec (does not include added species) ---
  bv_speci = Nv_MHD + 1
  ev_speci = bv_speci + Nspec - 1
  v_variab(bv_speci:ev_speci) = speci(1:Nspec)%density

! index of neutrals, species on grain mantles, species on cores, positive ions, negative ions
  bv_neu = b_neu + bv_speci - 1
  ev_neu = e_neu + bv_speci - 1
  bv_gra = b_gra + bv_speci - 1
  ev_gra = e_gra + bv_speci - 1
  bv_cor = b_cor + bv_speci - 1
  ev_cor = e_cor + bv_speci - 1
  bv_ion = b_ion + bv_speci - 1
  ev_ion = e_ion + bv_speci - 1
  bv_neg = b_neg + bv_speci - 1
  ev_neg = e_neg + bv_speci - 1

!  print *, " bv_neu    =", bv_neu, ", ev_neu    =", ev_neu
!  print *, " bv_gra    =", bv_gra, ", ev_gra    =", ev_gra
!  print *, " bv_cor    =", bv_cor, ", ev_cor    =", ev_cor
!  print *, " bv_ion    =", bv_ion, ", ev_ion    =", ev_ion
!  print *, " bv_neg    =", bv_neg, ", ev_neg    =", ev_neg

!--- density of H2 levels (cm-3), number = NH2_lev ---
  bv_H2_lev = Nv_MHD + Nspec + 1
  ev_H2_lev = bv_H2_lev + NH2_lev - 1
  v_variab(bv_H2_lev:ev_H2_lev) = H2_lev(1:NH2_lev)%density

!  print *, " bv_H2_lev =", bv_H2_lev, ", ev_H2_lev =", ev_H2_lev

!--- density of CO levels (cm-3), number = NCO_lev ---
  bv_CO_lev = Nv_MHD + Nspec + NH2_lev + 1
  ev_CO_lev = bv_CO_lev + NCO_lev - 1
  v_variab(bv_CO_lev:ev_CO_lev) = dens_CO_lev(1:NCO_lev)
!  print *, " bv_CO_lev =", bv_CO_lev, ", ev_CO_lev =", ev_CO_lev

!--- density of SiO levels (cm-3), number = NSiO_lev ---
  bv_SiO_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + 1
  ev_SiO_lev = bv_SiO_lev + NSiO_lev - 1
  v_variab(bv_SiO_lev:ev_SiO_lev) = dens_SiO_lev(1:NSiO_lev)
!  print *, " bv_SiO_lev =", bv_SiO_lev, ", ev_SiO_lev =", ev_SiO_lev

!--- density of o-H2O levels (cm-3), number = NoH2O_lev ---
  bv_oH2O_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + 1
  ev_oH2O_lev = bv_oH2O_lev + NoH2O_lev - 1
  v_variab(bv_oH2O_lev:ev_oH2O_lev) = dens_oH2O_lev(1:NoH2O_lev)
!  print *, " bv_oH2O_lev =", bv_oH2O_lev, ", ev_oH2O_lev =", ev_oH2O_lev

!--- density of p-H2O levels (cm-3), number = NpH2O_lev ---
  bv_pH2O_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + 1
  ev_pH2O_lev = bv_pH2O_lev + NpH2O_lev - 1
  v_variab(bv_pH2O_lev:ev_pH2O_lev) = dens_pH2O_lev(1:NpH2O_lev)
!  print *, " bv_pH2O_lev =", bv_pH2O_lev, ", ev_pH2O_lev =", ev_pH2O_lev

!--- density of o-NH3 levels (cm-3), number = NoNH3_lev ---
  bv_oNH3_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + NpH2O_lev + 1
  ev_oNH3_lev = bv_oNH3_lev + NoNH3_lev - 1
  v_variab(bv_oNH3_lev:ev_oNH3_lev) = dens_oNH3_lev(1:NoNH3_lev)
!  print *, " bv_oNH3_lev =", bv_oNH3_lev, ", ev_oNH3_lev =", ev_oNH3_lev

!--- density of p-NH3 levels (cm-3), number = NpNH3_lev ---
  bv_pNH3_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + NpH2O_lev + NoNH3_lev &
                + 1
  ev_pNH3_lev = bv_pNH3_lev + NpNH3_lev - 1
  v_variab(bv_pNH3_lev:ev_pNH3_lev) = dens_pNH3_lev(1:NpNH3_lev)
!  print *, " bv_pNH3_lev =", bv_pNH3_lev, ", ev_pNH3_lev =", ev_pNH3_lev

!--- density of OH levels (cm-3), number = NOH_lev ---
  bv_OH_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + NpH2O_lev + NoNH3_lev &
                + NpNH3_lev + 1
  ev_OH_lev = bv_OH_lev + NOH_lev - 1
  v_variab(bv_OH_lev:ev_OH_lev) = dens_OH_lev(1:NOH_lev)
!  print *, " bv_OH_lev =", bv_OH_lev, ", ev_OH_lev =", ev_OH_lev

!--- density of A-type levels (cm-3), number = Natype_lev ---
  bv_atype_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + NpH2O_lev &
                        + NoNH3_lev + NpNH3_lev + NOH_lev + 1
  ev_atype_lev = bv_atype_lev + Natype_lev - 1
  v_variab(bv_atype_lev:ev_atype_lev) = dens_atype_lev(1:Natype_lev)
!  print *, " bv_atype_lev =", bv_atype_lev, ", ev_atype_lev =", ev_atype_lev
!  print *, dens_atype_lev(1:Natype_lev)

!--- density of E-type levels (cm-3), number = Netype_lev ---
  bv_etype_lev = Nv_MHD + Nspec + NH2_lev + NCO_lev + NSiO_lev + NoH2O_lev + NpH2O_lev &
                        + NoNH3_lev + NpNH3_lev + NOH_lev + Natype_lev + 1
  ev_etype_lev = bv_etype_lev + Netype_lev - 1
  v_variab(bv_etype_lev:ev_etype_lev) = dens_etype_lev(1:Netype_lev)
!  print *, " bv_etype_lev =", bv_etype_lev, ", ev_etype_lev =", ev_etype_lev


!--- log of this vector ---
  WHERE (v_variab > 0)
     v_lvariab = LOG(v_variab)
     ELSEWHERE
     v_lvariab = minus_infinity
  END WHERE

!---------------------------------------------
! calculate initial values of mass, momentum
! and energy fluxes
!---------------------------------------------
  CALL ENERGETIC_FLUXES

!
!-- Initialization of the arguments of dvode
!
!--- skip the initializations if transition from C- to J-shock
  if(time_up) go to 1000
  liw_V = d_v_var+30                     ! size of array iwork_V
  lrw_V = 22 + 9*d_v_var + 2*d_v_var**2  ! size of array rwork_V

  allocate(iwork_V(liw_V))
  allocate(rwork_V(lrw_V))
  allocate(rtol_V(d_v_var))
  allocate(atol_V(d_v_var))
  iwork_V = 0
  rwork_V = 0.0_DP
  rtol_V = 0.0_DP
  atol_V = 0.0_DP

1000 continue

  atol_V     = 0.d-10          ! relative (not absolute) error check  
  rtol_V     = Eps_V * 10.0_DP ! error = rtol_V * abs(y(i))
  rwork_V(5) = H0_V * 1.0D-10  ! step size on first step
! rwork_V(6) = 1.0D17          ! HMAX
! rwork_V(7) = H0_V * 1.0D-10  ! HMIN

!---------------------------------------------
! write information about the current model
! (parameters, species, chemistry, H2, grains)
! in the file file_info
!---------------------------------------------
  CALL WRITE_INFO

!--------------------------------------
!---    start of the integration    ---
!---  (aborted when again=.FALSE.)  ---
!--------------------------------------
  WRITE(*,*) '--- start of the integration ---'
  DO WHILE (again)
     counter=counter+1 ! counts the number of integration steps

 !---------------------------------------
 ! save some values before call to DRIVE
 !---------------------------------------
     distance_old = distance
     Tn_old       = Tn
     Vn_old       = Vn
     Vi_old       = Vi
     speci(1:Nspec)%Dens_old = speci(1:Nspec)%density
     H2_lev(1:NH2_lev)%Dens_old = H2_lev(1:NH2_lev)%density
     H2_lines(1:NH2_lines)%emiss_old = H2_lines(1:NH2_lines)%emiss
     dens_CO_lev_old(1:NCO_lev) = dens_CO_lev(1:NCO_lev)
     dens_SiO_lev_old(1:NSiO_lev) = dens_SiO_lev(1:NSiO_lev)
     dens_oH2O_lev_old(1:NoH2O_lev) = dens_oH2O_lev(1:NoH2O_lev)
     dens_pH2O_lev_old(1:NpH2O_lev) = dens_pH2O_lev(1:NpH2O_lev)
     dens_oNH3_lev_old(1:NoNH3_lev) = dens_oNH3_lev(1:NoNH3_lev)
     dens_pNH3_lev_old(1:NpNH3_lev) = dens_pNH3_lev(1:NpNH3_lev)
     dens_OH_lev_old(1:NOH_lev) = dens_OH_lev(1:NOH_lev)
     dens_atype_lev_old(1:Natype_lev) = dens_atype_lev(1:Natype_lev)
     dens_etype_lev_old(1:Netype_lev) = dens_etype_lev(1:Netype_lev)
     TVGRAD_CO_old(1:NCO_lev) = TVGRAD_CO(1:NCO_lev)
     TVGRAD_SiO_old(1:NSiO_lev) = TVGRAD_SiO(1:NSiO_lev)
     TVGRAD_oH2O_old(1:NoH2O_lev,1:NoH2O_lev) = TVGRAD_oH2O(1:NoH2O_lev,1:NoH2O_lev)
     TVGRAD_pH2O_old(1:NpH2O_lev,1:NpH2O_lev) = TVGRAD_pH2O(1:NpH2O_lev,1:NpH2O_lev)
     TVGRAD_oNH3_old(1:NoNH3_lev,1:NoNH3_lev) = TVGRAD_oNH3(1:NoNH3_lev,1:NoNH3_lev)
     TVGRAD_pNH3_old(1:NpNH3_lev,1:NpNH3_lev) = TVGRAD_pNH3(1:NpNH3_lev,1:NpNH3_lev)
     TVGRAD_OH_old(1:NOH_lev,1:NOH_lev) = TVGRAD_OH(1:NOH_lev,1:NOH_lev)
     TVGRAD_atype_old(1:Natype_lev,1:Natype_lev) = TVGRAD_atype(1:Natype_lev,1:Natype_lev)
     TVGRAD_etype_old(1:Netype_lev,1:Netype_lev) = TVGRAD_etype(1:Netype_lev,1:Netype_lev)

 !------------------------------------------------------------
 !           numerical integration
 ! DVODE, as modified by Pierre Hily-Blant
 !
 ! remark :
 ! --------
 !         most of physical variables are calculated in
 !         DIFFUN, so they are available after call to DRIVE.
 !------------------------------------------------------------

     do while (.not. its_OK)

       call dvode(diffun, d_v_var, v_lvariab, T0_V, Tout_V, itol_V, &
                  rtol_V, atol_V, Itask_V, Istate_V, iopt_V, rwork_V, lrw_V, &
                  iwork_V, liw_V, jaco, MF_V, rpar_V, ipar_V)

       call diffun(d_v_var,Tout_V,v_lvariab,v_dvariab)  ! Refresh values of variables (PL)

       Hdone  = rwork_V(11)
       Hnext  = rwork_V(12)

       if (Istate_V == 2) then

         its_OK = .TRUE.

         call DVINDY (Tout_V, 1, rwork_V(21), d_v_var, dlv_var, iflag)

         dVn = dlv_var(1) * Vn
         dVi = dlv_var(2) * Vi

!        write (75, "(1p,300es17.8e3)") Tout_V, Hdone, (v_lvariab(jj), jj=1,d_v_var)
!        write (76, "(1p,300es18.9e3)") Tout_V, Hdone, &
!                     (dlv_var(jj), jj=1,d_v_var)
!                     (dlv_var(jj)*v_variab(jj), jj=1,d_v_var)

       else if (Istate_V == -1) then

         print *, " "
         print *, " Istate_V = ", Istate_V
         print *, " MXSTEP = ", iwork_V(6)
         print *, " "
         Tout_V = T0_V
         Hnext = Hnext * 0.1_DP
         iwork_V(6) = iwork_V(6) + 100
         Istate_V = 3

       else if (Istate_V == -2) then
         print *, " Istate_V = ", Istate_V
         print *, "  To be done"
         stop
       else if (Istate_V == -3) then
         print *, " Istate_V = ", Istate_V
         print *, "  To be done"
         stop
       else if (Istate_V == -4) then
         print *, " "
         print *, " Istate_V = ", Istate_V
         print *, " IMXER = ", iwork_V(16), v_variab(iwork_V(16))
         print *, " T = ", T0_V
         print *, " Tout = ", Tout_V
         print *, " Hnext = ", Hnext
         print *, " "
         Tout_V = T0_V
         Hnext = Hnext * 0.01_DP
         Istate_V = 2
       else if (Istate_V == -5) then
         print *, " Istate_V = ", Istate_V
         print *, "  To be done"
         stop
       else if (Istate_V == -6) then
         print *, " Istate_V = ", Istate_V
         print *, "  To be done"
         stop
       else
         print *, " Istate_V = ", Istate_V
         print *, "  Unknown error"
         stop
       end if

     end do

     its_OK = .FALSE.

     write (70,'(30ES20.11E3)') Tout_V, Hdone, &
                 Tn, Ti, Te, &
                 tr_1, tr_2, tr_3, tr_4, tr_5, &
                 tr_6, tr_7, tr_8, tr_9, tr_10, &
                 tr_11, tr_12, tr_13, tr_14, tr_15

 !-------------------------------------------
 ! calculate physical variables that are not
 ! already calculated in DIFFUN
 !-------------------------------------------
 !--- distance (cm) ---
     distance  = Tout_V
     dist_step = distance - distance_old

 !--- flow times (yr) of neutrals and ions ---
     timeN = timeN + 0.5_DP*dist_step*(1._DP/Vn_old+1._DP/Vn)/YEARsec
     timeI = timeI + 0.5_DP*dist_step*(1._DP/Vi_old+1._DP/Vi)/YEARsec

 !--------------------------------------------------------------
 ! calculate the column density (cm-2) of each chemical species
 ! and each H2 energy level.
 ! note : trapezium integration rule applied
 !--------------------------------------------------------------
 ! a) chemical species
     speci(1:Nspec)%Col_dens = speci(1:Nspec)%Col_dens + &
          0.5_DP * dist_step * (speci(1:Nspec)%Dens_old + speci(1:Nspec)%density)

 ! b) H2 levels
     H2_lev(1:NH2_lev)%Col_dens = H2_lev(1:NH2_lev)%Col_dens + &
          0.5_DP * dist_step * (H2_lev(1:NH2_lev)%Dens_old + H2_lev(1:NH2_lev)%density)


 !---------------------------------------------------
 ! Calculate the integrated intensity (erg/s/cm2/sr)
 ! of each H2 line and each fine structure line.
 !---------------------------------------------------
 ! a) H2 lines
     H2_lines(1:NH2_lines)%intensity = H2_lines(1:NH2_lines)%intensity + &
          0.5_DP * dist_step / (4._DP*pi) * &
          (H2_lines(1:NH2_lines)%emiss_old + H2_lines(1:NH2_lines)%emiss)

 ! b) integrated line intensities
     CALL LINE_INTEG

 !---------------------------------------------------------
 ! Calculate the integrated CO line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      CO_lines(1:NCO_lev) = CO_lines(1:NCO_lev) + &
          0.5_DP * dist_step * Vgrad * &
          (TVGRAD_CO_old(1:NCO_lev) + TVGRAD_CO(1:NCO_lev))
          
      CO_cdens(1:NCO_lev) = CO_cdens(1:NCO_lev) + &
          0.5_DP * dist_step * (dens_CO_lev_old(1:NCO_lev) &
          + dens_CO_lev(1:NCO_lev))
 !---------------------------------------------------------
 ! Calculate the integrated SiO line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      SiO_lines(1:NSiO_lev) = SiO_lines(1:NSiO_lev) + &
          0.5_DP * dist_step * Vgrad * &
          (TVGRAD_SiO_old(1:NSiO_lev) + TVGRAD_SiO(1:NSiO_lev))
          
      SiO_cdens(1:NSiO_lev) = SiO_cdens(1:NSiO_lev) + &
          0.5_DP * dist_step * (dens_SiO_lev_old(1:NSiO_lev) &
          + dens_SiO_lev(1:NSiO_lev))
 !---------------------------------------------------------
 ! Calculate the integrated o-H2O line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      oH2O_lines(1:NoH2O_lev,1:NoH2O_lev) = oH2O_lines(1:NoH2O_lev,1:NoH2O_lev) + &
          0.5_DP * dist_step * Vgrad * (TVGRAD_oH2O_old(1:NoH2O_lev,1:NoH2O_lev) &
          + TVGRAD_oH2O(1:NoH2O_lev,1:NoH2O_lev))
          
      oH2O_cdens(1:NoH2O_lev) = oH2O_cdens(1:NoH2O_lev) + &
          0.5_DP * dist_step * (dens_oH2O_lev_old(1:NoH2O_lev) &
          + dens_oH2O_lev(1:NoH2O_lev))
 !---------------------------------------------------------
 ! Calculate the integrated p-H2O line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      pH2O_lines(1:NpH2O_lev,1:NpH2O_lev) = pH2O_lines(1:NpH2O_lev,1:NpH2O_lev) + &
          0.5_DP * dist_step * Vgrad * (TVGRAD_pH2O_old(1:NpH2O_lev,1:NpH2O_lev) &
          + TVGRAD_pH2O(1:NpH2O_lev,1:NpH2O_lev))
          
      pH2O_cdens(1:NpH2O_lev) = pH2O_cdens(1:NpH2O_lev) + &
          0.5_DP * dist_step * (dens_pH2O_lev_old(1:NpH2O_lev) &
          + dens_pH2O_lev(1:NpH2O_lev))
 !---------------------------------------------------------
 ! Calculate the integrated o-NH3 line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      oNH3_lines(1:NoNH3_lev,1:NoNH3_lev) = oNH3_lines(1:NoNH3_lev,1:NoNH3_lev) + &
          0.5_DP * dist_step * Vgrad * (TVGRAD_oNH3_old(1:NoNH3_lev,1:NoNH3_lev) &
          + TVGRAD_oNH3(1:NoNH3_lev,1:NoNH3_lev))
          
      oNH3_cdens(1:NoNH3_lev) = oNH3_cdens(1:NoNH3_lev) + &
          0.5_DP * dist_step * (dens_oNH3_lev_old(1:NoNH3_lev) &
          + dens_oNH3_lev(1:NoNH3_lev))
 !---------------------------------------------------------
 ! Calculate the integrated p-NH3 line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      pNH3_lines(1:NpNH3_lev,1:NpNH3_lev) = pNH3_lines(1:NpNH3_lev,1:NpNH3_lev) + &
          0.5_DP * dist_step * Vgrad * (TVGRAD_pNH3_old(1:NpNH3_lev,1:NpNH3_lev) &
          + TVGRAD_pNH3(1:NpNH3_lev,1:NpNH3_lev))
          
      pNH3_cdens(1:NpNH3_lev) = pNH3_cdens(1:NpNH3_lev) + &
          0.5_DP * dist_step * (dens_pNH3_lev_old(1:NpNH3_lev) &
          + dens_pNH3_lev(1:NpNH3_lev))
 !---------------------------------------------------------
 ! Calculate the integrated OH line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      OH_lines(1:NOH_lev,1:NOH_lev) = OH_lines(1:NOH_lev,1:NOH_lev) + &
          0.5_DP * dist_step * Vgrad * (TVGRAD_OH_old(1:NOH_lev,1:NOH_lev) &
          + TVGRAD_OH(1:NOH_lev,1:NOH_lev))
          
      OH_cdens(1:NOH_lev) = OH_cdens(1:NOH_lev) + &
          0.5_DP * dist_step * (dens_OH_lev_old(1:NOH_lev) &
          + dens_OH_lev(1:NOH_lev))
 !---------------------------------------------------------
 ! Calculate the integrated A-type line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      atype_lines(1:Natype_lev,1:Natype_lev) = atype_lines(1:Natype_lev,1:Natype_lev) + &
          0.5_DP * dist_step * Vgrad * (TVGRAD_atype_old(1:Natype_lev,1:Natype_lev) &
          + TVGRAD_atype(1:Natype_lev,1:Natype_lev))
          
      atype_cdens(1:Natype_lev) = atype_cdens(1:Natype_lev) + &
          0.5_DP * dist_step * (dens_atype_lev_old(1:Natype_lev) &
          + dens_atype_lev(1:Natype_lev))
 !---------------------------------------------------------
 ! Calculate the integrated E-type line intensities (K km s-1)
 ! and level column densities (cm-2)
 !---------------------------------------------------------
      etype_lines(1:Netype_lev,1:Netype_lev) = etype_lines(1:Netype_lev,1:Netype_lev) + &
          0.5_DP * dist_step * Vgrad * (TVGRAD_etype_old(1:Netype_lev,1:Netype_lev) &
          + TVGRAD_etype(1:Netype_lev,1:Netype_lev))
          
      etype_cdens(1:Netype_lev) = etype_cdens(1:Netype_lev) + &
          0.5_DP * dist_step * (dens_etype_lev_old(1:Netype_lev) &
          + dens_etype_lev(1:Netype_lev))
 !------------------------------------------------
 ! calculate mass, momentum and energy fluxes
 ! and check the conservation of these quantities
 !------------------------------------------------
     CALL ENERGETIC_FLUXES

 !-----------------------
 ! STOP integration if :
 !-----------------------
 ! (1) Maximum number of steps has been reached
     IF (counter == Nstep_max) THEN
        again = .FALSE.
        message_end = "Maximum number of steps has been reached."
     END IF

 ! (2) Maximum evolution time has been reached
     IF (timeN >= duration_max) THEN
        again = .FALSE.
        message_end = "Maximum evolution time has been reached."
     END IF

 ! (3) Maximum shock length has been reached
     IF (distance >= length_max) THEN
        again = .FALSE.
        message_end = "Maximum shock length has been reached."
     END IF

 ! (4) drift velocity < DeltaVmin
!     IF (ABS_DeltaV < DeltaVmin * 1.0d-1) THEN
!!    IF ((Vn - Vi) < 10.0_DP) THEN
!        again = .FALSE.
!        message_end = "Equilibrium has been reached."
!     END IF

 ! (5) we reach a sonic point
!    IF (Vn < (Vsound + 1.0D+1)) THEN
!       again = .FALSE.
!       message_end = "Sonic Point!"
!    END IF

 ! (6) Temperature is below 5.0 K
 !   IF (counter > 1000 .AND. Tn < 5.0_DP) THEN
 !      again = .FALSE.
 !      message_end = "temperature lower than 5 K"
 !   END IF

 ! (7) Numerical instabilities prevent conservation of ions
     IF (abs((DensityI-SUM(v_variab(bv_ion:ev_ion)))/DensityI) > 1.0e-2_DP) THEN
        again = .FALSE.
        message_end = "Ions are NOT conserved"
     END IF

!--- write to several output files ---
     if (mod(counter-1,NStep_w) == 0) then
     CALL WRITE_OUTPUT(again)

!--- write basic information to screen ---
     CALL WRITE_SCREEN
!     write (*,'(" C  :",1P,5E10.3)') pop_cat
!     write (*,'(" N  :",1P,5E10.3)') pop_nat
!     write (*,'(" O  :",1P,5E10.3)') pop_oat
!     write (*,'(" S  :",1P,5E10.3)') pop_sat
!     write (*,'(" Si :",1P,5E10.3)') pop_siat
!     write (*,'(" C+ :",1P,5E10.3)') pop_cpl
!     write (*,'(" N+ :",1P,5E10.3)') pop_npl
!     write (*,'(" O+ :",1P,5E10.3)') pop_opl
!     write (*,'(" S+ :",1P,5E10.3)') pop_spl
!     write (*,'(" Si+:",1P,5E10.3)') pop_sipl
!     print *, Cool_n, Cool_i, Cool_e
      endif

!--- set next requested output time
     Tout_V = Tout_V + 1.5_DP * Hnext

!---------------------------------------------------------------------------------
!--- change to single fluid model when Vi and Vn have re-equalized

!   if ( ( shock_type == "J" )                     .AND. &
!        (Vn <= Vi) ) then
         if (Vn <= Vi) then
      Nfluids = 1
    endif
!---------------------------------------------------------------------------------
!!   Switch off artificial viscosity after discontinuity in J shocks.
!!   If (viscosity = .TRUE.), as initially, then the derivative of the flow
!!   velocity v is determined by integrating a second order equation for v,
!!   which arises from including the artificial viscosity. This integration proceeds
!!   through the J-shock 'discontinuity', where v is varying rapidly, up to
!!   the point at which the velocity gradient becomes equal to that calculated
!!   neglecting artificial viscosity, to within a certain tolerance (see below).
!!   viscosity is then set equal to .FALSE., and the integration is pursued neglecting
!!   the artificial viscosity terms (which have become negligible).
!!   In addition, the logarithmic derivatives are integrated beyond this point.

    if ( ( shock_type == "J" )                     .AND. &
!        ( viscosity == .TRUE. )                   .AND. &
         ( viscosity)                   .AND. &
         ( ABS((save_dv + grad_V) / grad_V) < 2.0e-3_dp ) ) then
!    Switch back to C-type and change preJ to .FALSE.
      shock_type = "C"
      preJ = .FALSE.
      viscosity = .FALSE.
   rtol_V = rtol_V * 1.d0   ! Reset to a less demanding value after the "discontinuity"
    endif
 !---------------------------------------------------------------------------------
 !--- switch to J-type if specified evolution time has been reached
     time_up = (shock_type == "C") .AND. (preJ) .AND. (timeI >= timeJ)
     if(time_up) then
      shock_type = "J"
      viscosity = .TRUE.
 !--- reset DVODE parameters
   T0_V   = distance    ! value of the integration variable
   H0_V   = XLL * 1.d-6 ! step length
   Tout_V = T0_V + H0_V ! value at which output should occur
   MF_V = 22            ! chosen numerical procedure
   Itask_V = 1          ! the task DVODE should perform
   Istate_V = 1         ! the outcome of the task (= 2 on normal return)

       go to 999
     endif
 !---------------------------------------------------------------------------------

 !--------------------------
 !--- end of integration ---
 !--------------------------
  END DO

!--- write parameters for H2 excitation diagram
  CALL WRITE_EXCIT

!--- write final abundances in input format
  CALL WRITE_SPECIES

! write messages to screen
  WRITE(*,*) message_end
  WRITE(*,'(80("-"))')
  IF (err_cool) WRITE(*,'("*** WARNING, check the file:", A)') n_err_cool

!------------------
! file closure
!------------------
  CLOSE(f_err_cool)

!----------------------
! Array deallocations
!----------------------
! arrays allocated in MHD:

  DEALLOCATE(v_variab)
  DEALLOCATE(v_lvariab)
  DEALLOCATE(v_dvariab)
  DEALLOCATE(v_l_var)
  DEALLOCATE(v_l_der)
  DEALLOCATE(YN)

! arrays allocated elsewhere: 

  DEALLOCATE(speci)
  DEALLOCATE(index_VJ_H2)
  DEALLOCATE(H2_lines)

!#############################################################################
!#############################################################################
!##                                                                         ##
!##                 END OF THE MHD SHOCK CODE                               ##
!##                                                                         ##
!#############################################################################
!#############################################################################

END PROGRAM MHD

