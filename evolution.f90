MODULE MODULE_EVOLUTION

  !*****************************************************************************
  !** The module 'MODULE_EVOLUTION' contains 3 subroutines (CHEMISTRY, SOURCE **
  !** and DIFFUN) used to calculate at each step :                            **
  !**           - the source terms (CHEMISTRY and SOURCE)                     **
  !**           - the partial derivatives for DVODE integration package.      **
  !*****************************************************************************
! changes made by SC since last release 16II04: search for "SC" in file
! 
! coding of the evolution.f90 file name:
! - kn: bug corrected in Newton solving for dVn/dz when KN=1 (2 May 2006)
!       + minimum gradient = thermal gradient = Vsound/(Z+dmin)
! - rg: use only core species to calculate mean radius of grain core
!       assuming a Matthis distribution (2 May 2006)
! - lay: add mantle layers to grain cross section assuming Nsites = cste
!        using: rsquare_grain = (r_core + Nlayer * layer_thickness)^2
! modified variable: Nsites_grain = 5e4 instead of 1e6 (variables_mhd)
! new variables needed: layer_thickness = 2e-8 cm (variables_mhd)
!                       R_gr_scale12 (variables_mhd, calculated in initialize)
! - G: species C60,C60+,C60- replaced by G,G+,G-

  IMPLICIT NONE
  INCLUDE "precision.f90"

  !-----------------------------------------------
  ! vectors containing variables and derivatives
  ! at the last call to DIFFUN
  ! dimension = 1:d_v_var
  ! allocated and initialized in MHD
  !-----------------------------------------------
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: v_dvariab ! Y * dY at the last call of DIFFUN
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: v_l_var   ! Y at the last call of DIFFUN
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: v_l_der   ! DY at the last call to DIFFUN

  REAL(KIND=DP), PARAMETER :: tiny = 1.0D-100 ! avoids overflow in DIFFUN

  !------------------------------------------------------------------
  ! source terms for each fluid
  ! (N -> neutrals, I -> ions, NEG -> negative ions)
  ! = rates of change in number, mass, momentum and energy densities
  ! calculated in CHEMISTRY, and in SOURCE
  !------------------------------------------------------------------
  REAL(KIND=DP) :: CupN, CupI, CupNEG        ! increase in number density (cm-3.s-1)
  REAL(KIND=DP) :: CdownN, CdownI, CdownNEG  ! decrease in number density (cm-3.s-1)
  REAL(KIND=DP) :: YNn, YNi, YNneg           ! change in number density (cm-3.s-1)
  REAL(KIND=DP) :: Sn, Si, Sneg, Score       ! change in mass density (g.cm-3.s-1)
  REAL(KIND=DP) :: An, Ai, Aneg              ! change in momentum density (g.cm-2.s-2)
  REAL(KIND=DP) :: Bn, Bi, Bneg              ! change in energy density (erg.cm-3.s-1)
  REAL(KIND=DP) :: B_i_n, B_neg_n, B_e_n, B_inelastic_e_n, B_grain_n
  REAL(KIND=DP) :: H2_int_energy             ! internal energy of H2 (K cm-3)

  !------------------------------
  ! source terms for each species
  ! dimension = 0:Nspec
  ! calculated in CHEMISTRY
  !------------------------------
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: YN ! change in number density (cm-3.s-1)

  PRIVATE :: tiny, CHEMISTRY, SOURCE

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS


  SUBROUTINE CHEMISTRY
    !---------------------------------------------------------------------------
    ! called by :
    !     DIFFUN
    ! purpose :
    !     calculates source terms related to chemistry. Each type of reaction
    !     has its own reaction rate (with different meanings and units for
    !     gamma, alpha and beta ).
    !     For each species :
    !         * YN : rate at which it is formed through chemical reactions
    !                [1], eqs. (14) and (20)
    !         * YS : rate of change in mass (g.cm-3.s-1)
    !         * YA : rate of change in momentum (g.cm-2.s-2)
    !         * YB : rate of change in energy (erg.cm-3.s-1)
    !     For each fluid (x= n -> neutrals, i -> ions, neg -> negative ions and e-)
    !         * CupX   : rate of increase in number density (cm-3.s-1)
    !         * CdownX : rate of decrease in number density (cm-3.s-1)
    !         * Nx     : rate of change in number density (cm-3.s-1)
    !         * Sx     : rate of change in mass density (g.cm-3.s-1)
    !         * Ax     : rate of change in momentum density (g.cm-2.s-2)
    !         * Bx     : rate of change in energy density (erg.cm-3.s-1)
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     YN
    !     CupN,   CupI,   CupNEG
    !     CdownN, CdownI, CdownNEG
    !     YNn,    YNi,    YNneg
    !     Sn,     Si,     Sneg
    !     An,     Ai,     Aneg
    !     Bn,     Bi,     Bneg
    ! references :
    !     [1] Flower et al., 1985, MNRAS 216, 775
    !     [2] Flower et al., 1986, MNRAS 218, 729
    !     [3] Pineau des Forets et al., 1986, MNRAS 220, 801
    !     [4] Viala et al., 1988, A&A 190, 215
    !---------------------------------------------------------------------------
    USE MODULE_CHEM_REACT
    USE MODULE_GRAINS
    USE MODULE_H2
!   USE MODULE_PHYS_VAR,    ONLY : Av, RAD, Zeta, Te, Tn, Vi, Vn, nH, ABS_DeltaV, NH2_lev
    USE MODULE_PHYS_VAR
    USE MODULE_CONSTANTS,   ONLY : pi, kB, me, EVerg, Zero
    USE MODULE_DEBUG_JLB
    IMPLICIT NONE

    ! index 0 -> inexistent species, Nspec_plus -> added species (e-, GRAIN, ...)
    REAL(KIND=DP), DIMENSION(0:Nspec_plus) :: up_N, up_MV, up_MV2, up_DE
    REAL(KIND=DP), DIMENSION(0:Nspec_plus) :: down_N, down_MV, down_MV2
    REAL(KIND=DP), DIMENSION(0:Nspec) :: YS, YA, YB
    INTEGER(KIND=LONG) :: i, IR1, IR2, IP1, IP2, IP3, IP4
    REAL(KIND=DP) :: massR1, massR2, densityR1, densityR2, V_R1, V_R2, Mass_prod, Nprod_m1
    REAL(KIND=DP) :: Vcm, Vcm2
    REAL(KIND=DP) :: creation_rate, destr_R1, destr_R2
    REAL(KIND=DP) :: creation_Vcm, creation_Vcm2, creation_DE
    REAL(KIND=DP) :: destr1_Vcm, destr2_Vcm, destr1_Vcm2, destr2_Vcm2
    REAL(KIND=DP) :: exp_factor, exp_factor1, exp_factor2, coeff, Tr, Ts, energy, Teff
    REAL(KIND=DP) :: react_rate, r_alph, r_beta, r_gamm
    REAL(KIND=DP) :: fracH2
    INTEGER :: lev

    REAL(KIND=DP) :: stick_H_gr = 1.0_DP

    ! initialization
    up_N   = Zero ; up_MV   = Zero ; up_MV2   = Zero ; up_DE = Zero
    down_N = Zero ; down_MV = Zero ; down_MV2 = Zero

    Sel_ch_H2(1:NH2_lev) = Zero
    Sel_ne_H2(1:NH2_lev) = Zero
    Sel_io_H2(1:NH2_lev) = Zero
    Sel_el_H2(1:NH2_lev) = Zero

    !---------------------------------------------
    ! calculates creation and destruction terms
    ! for each reaction (i=index of the reaction)
    !---------------------------------------------
    DO i=1,Nreact
       !--- index, mass, abundance, velocity of each reactant and product ---
       IR1       = react(i)%R(1)       ; IR2       = react(i)%R(2)
       IP1       = react(i)%P(1)       ; IP2       = react(i)%P(2)
       IP3       = react(i)%P(3)       ; IP4       = react(i)%P(4)
       massR1    = speci(IR1)%mass     ; massR2    = speci(IR2)%mass
       densityR1 = speci(IR1)%density  ; densityR2 = speci(IR2)%density
       V_R1      = speci(IR1)%velocity ; V_R2      = speci(IR2)%velocity
       Mass_prod = react(i)%Mass_prod  ! sum of the masses of the products
       Nprod_m1  = react(i)%Nprod_m1   ! number of products minus one

       r_alph = react(i)%alpha
       r_beta = react(i)%beta
       r_gamm = react(i)%gamma

       !--- center of mass velocity (cm/s) and its square ---
       Vcm  = (massR1*V_R1 + massR2*V_R2) / (massR1+massR2)
       Vcm2 = Vcm*Vcm

       !---------------------------------------------------------------------
       ! now calculate the reaction rate coefficient (cm3.s-1) according to 
       ! the type of the reaction (effective temperature, significance of gamma,
       ! alpha, beta, ...).
       !---------------------------------------------------------------------
       SELECT CASE (react(i)%type)

          !-----------------
          ! photo reactions
          !-----------------
       CASE ('PHOTO')
          !--- reaction rate coefficient (cm3.s-1) ---

          react_rate = r_gamm * EXP(-r_beta*AV) * RAD

          creation_rate = react_rate * densityR1 * densityR2

          !---------------------------------------
          ! cosmic ray ionization or dissociation
          !---------------------------------------
       CASE ('CR_IO')
          !--- reaction rate coefficient (cm3.s-1) ---
          ! We take into account :
          !     0) direct cosmic ray (CR) induced ionization/dissociation
          !     1) the CR induced secondary photons
          !     2) the UV photons produced by collisional excitation of H2 by electrons
          !        following shock heating
          !     3) allow for the differential compression of the charged grains and the
          !        neutrals (Vi/Vn). This factor is applicable to cases (1) and
          !        (2) above, but not to case (0). However, the direct CR induced
          !        processes are completely negligible within the shock, which is
          !        where Vi/Vn differs from 1.
          ! note that react(i)%beta is preset at a very large value (in REACTION_TYPE)
          ! for direct cosmic ray ionization : react(i)%beta = 1.0D8

! JLB + GP : 13 juin 2002  - correction de la fraction de H2 dans le gaz
! JLB     exp_factor = MIN(r_beta/Te,180._DP)
          exp_factor = r_beta / Te
          coeff = SQRT(8.0_DP*kB*Te/(pi*me)) * 1.0D-16 * (4.29D-1 + 6.14D-6*Te)
          if (r_beta < 0.99D8) then
            fracH2 = Dens_H2 / (Dens_H + Dens_H2)
          else
            fracH2 = 1.0_DP
          endif
          react_rate = Zeta * ((Tn/300._DP)**r_alph) * r_gamm &
                     + Dens_e * coeff * EXP(-exp_factor) / 0.15_DP * r_gamm * Vi / Vn
          react_rate = react_rate * fracH2

          creation_rate = react_rate * densityR1 * densityR2

          !-------------------------------------------
          ! cosmic ray induced desorption from grains
          !-------------------------------------------
       CASE ('CR_DE')
          !--- reaction rate coefficient (cm3.s-1) ---

          react_rate = r_gamm * Dens_GRAIN / Dens_ongrains * pi * Rsquare_GRAIN

          creation_rate = react_rate * densityR1 * densityR2

          !-------------------------------
          ! H2 and HD formation on grains
          !-------------------------------
       CASE ('H2_FO')
          !--- reaction rate coefficient (cm3.s-1) ---

!         react_rate = r_gamm * nH/Dens_H * ((Tn/300._DP)**r_alph)

! JLB - 21 juin 2001 - Correction du taux de formation de H2

          !--- effective temperature ---
          Teff_grain = massR1*(Vi-Vn)*(Vi-Vn)/(3.0_DP*kB) + speci(IR1)%temperature
          stick_H_gr = 1.0_DP / sqrt(1.0_DP + Teff_grain / 30.0_DP)

          !--- reaction rate coefficient (cm3.s-1) (Andersson & Wannier for H adsorption) ---
          !    or try Hollenbach & McKee (1979, ApJ Suppl, 41, 555)
          stick_H_gr = 1.0_DP / (1.D0 + 0.04_DP*(Teff_grain+Tgrain)**0.5_DP &
                                 + 2.0D-3*Teff_grain + 8.0D-6*Teff_grain*Teff_grain)
          react_rate = stick_H_gr * pi * Rsquare_GRAIN * &
               SQRT(8.0_DP*kB*Teff_grain/massR1/pi)

          creation_rate = react_rate * Dens_grain * densityR2

          if (IR2 == ind_H) then
            For_gr_H2 = creation_rate
          endif
          !----------------------------------------
          ! three-body reactions on grains surface
          !----------------------------------------
       CASE ('THREE')
          ! exclude grain reactions when calculating steady state conditions
      if (shock_type == "S") then
          creation_rate = 0._DP
                             else
          !--- effective temperature ---

          Teff_grain = massR1*(Vi-Vn)*(Vi-Vn)/(3.0_DP*kB) + speci(IR1)%temperature

          !--- reaction rate coefficient (cm3.s-1) (Andersson & Wannier for H adsorption) ---

          react_rate = r_gamm * pi * Rsquare_GRAIN * &
               SQRT(8.0_DP*kB*Teff_grain/massR1/pi)
          react_rate = react_rate / (Teff_grain/r_beta+1.0_DP)**2.0_DP
          react_rate = react_rate * Dens_GRAIN / Dens_ongrains

          creation_rate = react_rate * densityR1 * densityR2
      endif
          !------------------------------
          ! sputtering of grain mantles
          !------------------------------
       CASE ('SPUTT')
          ! exclude grain reactions when calculating steady state conditions
      if (shock_type == "S") then
          creation_rate = 0._DP
                             else
          !--- effective temperature of grain sputtering ---

          Ts = massR2*(V_R2-V_R1)*(V_R2-V_R1)/(3.0_DP*kB)
          Tr =  speci(IR2)%temperature
          Teff_grain = Ts + Tr

          !--- reaction rate coefficient (cm3.s-1) ---
          ! reaction  : X* + IR2 -> X + IR2 + GRAIN
          ! Nlayers_grain = number of layers in the mantles
          ! IF Nlayers_grain > 1 => the detachment probability is
          !                         gamma * n(X*)/n(ongrains)
          ! IF Nlayers_grain < 1 => the detachment probability is
          !                         gamma * n(X*)/(Nsites_grain*n(grains))

          Nlayers_grain = Dens_ongrains/Nsites_grain/Dens_GRAIN
          exp_factor1   = (4.0_DP*r_beta-1.5_DP*Ts)/Tr
          exp_factor2   = 4.0_DP*r_beta/Teff_grain
          exp_factor    = MAX(exp_factor1,exp_factor2)
! JLB     exp_factor    = MIN(exp_factor,180._DP)

          react_rate = 4.0_DP * r_gamm * Dens_GRAIN/Dens_ongrains * &
               pi*Rsquare_GRAIN * SQRT(8.0_DP*kB*Teff_grain/massR2/pi)
          react_rate = react_rate * (1.0_DP+2.0_DP*Teff_grain/(4.0_DP*r_beta)) * &
               EXP(-exp_factor)
          react_rate = react_rate * MIN(Nlayers_grain,1.0_DP)

          creation_rate = react_rate * densityR1 * densityR2
      endif 

          !------------------------
          ! erosion of grain cores
          !------------------------
       CASE ('EROSI')
          ! exclude grain reactions when calculating steady state conditions
      if (shock_type == "S") then
          creation_rate = 0._DP
                             else
          !--- reaction rate coefficient (cm3.s-1) ---
          react_rate = Zero
          !--- absolute value of the ion-neutral velocity drift ---

          energy = massR2 * ABS_DeltaV*ABS_DeltaV / 2.0_DP / EVerg ! energy (eV)

          IF (energy > r_alph) THEN
          react_rate = pi * Rsquare_GRAIN * ABS_DeltaV * &
               r_gamm * EXP(-r_beta/(Energy-r_alph))
          END IF

          ! the rate coefficient is proportional to the grain density
          ! and not to the core species density

! Original form:
!         react_rate = react_rate * Dens_grain / densityR1

          if (densityR1 > Dens_grain) then
            creation_rate = react_rate * Dens_grain * densityR2
          else
            creation_rate = react_rate * densityR1 * densityR2
          endif
      endif 

          !-------------------------
          ! adsorption on to grains
          !-------------------------
       CASE ('ADSOR')
          ! exclude grain reactions when calculating steady state conditions
      if (shock_type == "S") then
          creation_rate = 0._DP
                             else
          !--- effective temperature (K) ---

          Teff_grain    = massR1*(Vi-Vn)*(Vi-Vn)/(3.0_DP*kB) + speci(IR1)%temperature

          !--- reaction rate coefficient (cm3.s-1) ---
          ! here, gamma = sticking coefficient
          ! react_rate = STICK * SIGMA(GRAIN) * V(MOYENNE)

          react_rate = r_gamm * pi * Rsquare_GRAIN * &
               SQRT(8.0_DP*kB*Teff_grain/massR1/pi)
          IF(r_beta /=Zero) THEN
             ! Hollenbach & McKee (1979, ApJ Suppl, 41, 555)
             react_rate = react_rate / (1.0_DP + 0.04_DP*(Teff_grain+Tgrain)**0.5_DP + &
                  2.0D-3 * Teff_grain + 8.0D-6*Teff_grain*Teff_grain)
          END IF

          creation_rate = react_rate * densityR1 * densityR2
      endif 

          !--------------------------------------------------------------
          ! collisional dissociation of H2
          !--------------------------------------------------------------
       CASE ('DISSO')
          !--- effective temperature of the reaction                               ---
          !--- see [1], eq. (43) and [3], eqs.(2.1)-(2.3)                          ---
          !--- for reactions between members of the same fluid : Ts = 0, Teff = Tr ---

          Ts = massR1*massR2*(V_R1-V_R2)*(V_R1-V_R2)/(massR1+massR2)/3.0_DP/kB
          Tr = (massR1*speci(IR2)%temperature + massR2*speci(IR1)%temperature) / &
               (massR1+massR2)
          Teff = Ts + Tr

          !--- destruction is computed state by state                             ---
          !--- reduce dissociation energy r_beta by the excitation energy         ---
          !--- of the level

          do lev=1,NH2_lev

            !--- exponential factor ([3], PAGE 810) ---
            exp_factor1 = (r_beta - H2_lev(lev)%energy - 3.0_DP*Ts) / Tr
            exp_factor2 = (r_beta - H2_lev(lev)%energy) / Teff
            exp_factor  = MAX(exp_factor1,exp_factor2)
            exp_factor  = MAX(exp_factor1,0.0_DP)

            Sel_rx_H2(lev) = EXP(-exp_factor)

          end do

          !--- reaction rate coefficient (cm3.s-1) ---
          Sel_rx_H2 = r_gamm * Sel_rx_H2 * ((Teff/300._DP)**r_alph)

          !--- Note: we avoid dividing by n(H2) just to be able to remultiply...

          react_rate = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_rx_H2))

          creation_rate = react_rate * densityR2

          Sel_ch_H2 = Sel_ch_H2 - Sel_rx_H2 * densityR2
          if (IR2 >= b_neu .AND. IR2 <= e_neu) then
            Sel_ne_H2 = Sel_ne_H2 - Sel_rx_H2 * densityR2
          else if (IR2 >= b_ion .AND. IR2 <= e_ion) then
            Sel_io_H2 = Sel_io_H2 - Sel_rx_H2 * densityR2
          else if (IR2 == ind_e) then
            Sel_el_H2 = Sel_el_H2 - Sel_rx_H2 * densityR2
          endif

          !--------------------------------------------------------------
          ! all other reactions
          ! (including reverse reactions added in ADD_INVERSE_REACTIONS)
          !--------------------------------------------------------------
       CASE ('OTHER', 'REVER')
          !--- effective temperature of the reaction                               ---
          !--- see [1], eq. (43) and [3], eqs.(2.1)-(2.3)                          ---
          !--- for reactions between members of the same fluid : Ts = 0, Teff = Tr ---

          Ts = massR1*massR2*(V_R1-V_R2)*(V_R1-V_R2)/(massR1+massR2)/3.0_DP/kB
          Tr = (massR1*speci(IR2)%temperature + massR2*speci(IR1)%temperature) / &
               (massR1+massR2)
          Teff = Ts + Tr

          !--- exponential factor ([3], PAGE 810) ---
          exp_factor1 = (r_beta-3.0_DP*Ts) / Tr
          exp_factor2 = r_beta / Teff
          exp_factor  = MAX(exp_factor1,exp_factor2)
! JLB     exp_factor  = MIN(exp_factor,180._DP)

          !--- reaction rate coefficient (cm3.s-1) ---
          react_rate = r_gamm * EXP(-exp_factor) * ((Teff/300._DP)**r_alph)

!  include temperature dependent Coulomb enhancement factor in the cases of neutralization
!  of charged grains, G+ and G-, by electrons and positive ions, respectively
          if (IP1 == ind_G) react_rate = react_rate * (1._DP + 450._DP/Teff)

          creation_rate = react_rate * densityR1 * densityR2

       END SELECT ! end of selection for the reaction 

       !--- creation rate per unit volume (cm-3.s-1)
       !--- destruction rate (s-1) of each reactant 
!  creation computed in 'select case' above
!  already multiplied by density  
!      creation_rate = react_rate * densityR1 * densityR2
!      destr_R1      = react_rate * densityR2
!      destr_R2      = react_rate * densityR1

       destr_R1      = creation_rate
       destr_R2      = creation_rate

       !--- source terms for destruction (multiplied by mass later) ---
       destr1_Vcm   = destr_R1 * Vcm
       destr2_Vcm   = destr_R2 * Vcm
       destr1_Vcm2  = destr_R1 * Vcm2
       destr2_Vcm2  = destr_R2 * Vcm2

       !--- source terms for creation ---
       creation_Vcm  = creation_rate * Vcm
       creation_Vcm2 = creation_rate * Vcm2
       creation_DE   = creation_rate * react(i)%DE

       !--- for the reactants, increase in the rate of destruction (s-1) of : ---
       !---      down_N(IR)   -> number                                       ---
       !---      down_MV(IR)  -> momentum                                     ---
       !---      down_MV2(IR) -> energy                                       ---

       down_N(IR1)   = down_N(IR1)   + destr_R1
       down_N(IR2)   = down_N(IR2)   + destr_R2
       down_MV(IR1)  = down_MV(IR1)  + destr1_Vcm
       down_MV(IR2)  = down_MV(IR2)  + destr2_Vcm
       down_MV2(IR1) = down_MV2(IR1) + destr1_Vcm2
       down_MV2(IR2) = down_MV2(IR2) + destr2_Vcm2

       !--- for the products, increase in the rate of creation (s-1) of :    ---
       !---      up_N(IP)   -> number                                        ---
       !---      up_MV(IP)  -> momentum                                      ---
       !---      up_MV2(IP) -> kinetic energy                                ---
       !---      up_DE(IP)  -> energy defect of the reaction                 ---
       !---                    DE is distributed over the reaction products  ---
       !---                    with a weighing factor : [1], eq. (31)        ---
       !---                        =  DE/(n-1)*[1-mass(i)/MASS] if n > 1     ---
       !---                        =  DE if n=1                              ---
       !---                    where n=number of produts, i=index of product ---
       !---                    and MASS=sum of mass(i), i=1..n               ---
       !---                    By this means, energy is conserved, whatever  ---
       !---                    the number of products in the reaction.       ---

       ! first product, always present !
       up_N(IP1)   = up_N(IP1)   + creation_rate
       up_MV(IP1)  = up_MV(IP1)  + creation_Vcm
       up_MV2(IP1) = up_MV2(IP1) + creation_Vcm2

       IF (Nprod_m1 == 0) THEN
          ! allow for reactions with only one product : H2_FO, ADSOR

!  JLB - May 2001
!  only 1/3 of H2 enthalpy of formation goes to kinetic energy 
!  Maybe should be corrected elsewhere ???
!         up_DE(IP1) = up_DE(IP1) + creation_DE
!  + (24 I 2002) account for loss of kinetic energy by H attaching to grain
!    (proposed by David)

          if (react(i)%type == "H2_FO") then
            if (ikinH2 == 1) then
              up_DE(IP1) = up_DE(IP1) + 0.5_DP * (creation_DE - creation_rate * H2_int_E * kB / EVerg) &
              - creation_rate * 1.5_DP * kB * speci(IR2)%temperature / EVerg
            else
              up_DE(IP1) = up_DE(IP1) + min(creation_DE / 3.0_DP, creation_DE - creation_rate * H2_int_E * kB / EVerg) &
              - creation_rate * 1.5_DP * kB * speci(IR2)%temperature / EVerg
            endif
          else
            up_DE(IP1) = up_DE(IP1) + creation_DE
          endif
       ELSE
          ! two or more products
          up_DE(IP1) = up_DE(IP1) + creation_DE / Nprod_m1 * &
                      (1.0_DP-speci(IP1)%mass/Mass_prod)
       ENDIF

       ! second product
       IF (IP2 > 0) THEN
          up_N(IP2)   = up_N(IP2)   + creation_rate
          up_MV(IP2)  = up_MV(IP2)  + creation_Vcm
          up_MV2(IP2) = up_MV2(IP2) + creation_Vcm2
          ! in this case, the reaction has at least two products !
          up_DE(IP2)  = up_DE(IP2)  + creation_DE/Nprod_m1 * &
                       (1.0_DP-speci(IP2)%mass/Mass_prod)
       END IF
       ! third product
       IF (IP3 > 0) THEN
          up_N(IP3)   = up_N(IP3)   + creation_rate
          up_MV(IP3)  = up_MV(IP3)  + creation_Vcm
          up_MV2(IP3) = up_MV2(IP3) + creation_Vcm2
          ! in this case, the reaction has at least three products !
          up_DE(IP3)  = up_DE(IP3)  + creation_DE / Nprod_m1 * &
                       (1.0_DP-speci(IP3)%mass/Mass_prod)
       END IF
       ! fourth product
       IF (IP4 > 0) THEN
          up_N(IP4)   = up_N(IP4)   + creation_rate
          up_MV(IP4)  = up_MV(IP4)  + creation_Vcm
          up_MV2(IP4) = up_MV2(IP4) + creation_Vcm2
          ! in this case, the reaction has four products !
          up_DE(IP4)  = up_DE(IP4)  + creation_DE / Nprod_m1 * &
                       (1.0_DP-speci(IP4)%mass/Mass_prod)
       END IF

    END DO ! end of the loop on chemical reactions

    !------------------------------------------------------------------
    ! Calculate for every chemical species the change per unit volume
    ! and time in :
    !       YN(i) -> number (cm-3.s-1)
    !       YS(i) -> mass (g.cm-3.s-1)
    !       YA(i) -> momentum (g.cm-2.s-2)
    !       YB(i) -> energy (erg.cm-3.s-1)
    !----------------------------------------------------------------

    YN(1:Nspec) = up_N(1:Nspec) - down_N(1:Nspec)
    YS(1:Nspec) = YN(1:Nspec) * speci(1:Nspec)%mass
    YA(1:Nspec) = (up_MV(1:Nspec) - down_MV(1:Nspec)) * speci(1:Nspec)%mass
    YB(1:Nspec) = (up_MV2(1:Nspec) - down_MV2(1:Nspec)) * speci(1:Nspec)%mass

    !--------------------------------------------------------------------------
    ! Calculate the source terms (the changes per unit volume and time)
    ! of each fluid by summing over the source terms for the chemical species :
    !      Cup    -> increase in number (cm-3.s-1)
    !      Cdown  -> decrease in number (cm-3.s-1)
    !      YN     -> net change in number (cm-3.s-1) [1], eq. (16,17,20)
    !      S      -> mass (g.cm-3.s-1)               [1], eq. (18,19)
    !      A      -> momentum (g.cm-2.s-2)           [1], eq. (21)
    !      B      -> energy (erg.cm-3.s-1)           [1], eq. (25,26,31)
    !--------------------------------------------------------------------------

    ! neutral species
    CupN   = SUM(up_N(b_neu:e_neu))
    CdownN = SUM(down_N(b_neu:e_neu))
    YNn    = SUM(YN(b_neu:e_neu))
    Sn     = SUM(YS(b_neu:e_neu))
    An     = SUM(YA(b_neu:e_neu))
    Bn     = 0.5_DP * SUM(YB(b_neu:e_neu)) + &
             SUM(up_DE(b_neu:e_neu)) * EVerg

    ! neutrals on mantles : add these terms to the neutrals
    CupN   = CupN   + SUM(up_N(b_gra:e_gra))
    CdownN = CdownN + SUM(down_N(b_gra:e_gra))
    YNn    = YNn    + SUM(YN(b_gra:e_gra))
    Sn     = Sn     + SUM(YS(b_gra:e_gra))
    An     = An     + SUM(YA(b_gra:e_gra))
    Bn     = Bn     + 0.5_DP * SUM(YB(b_gra:e_gra)) + &
             SUM(up_DE(b_gra:e_gra)) * EVerg

    ! neutrals eroded from grain cores : add these terms to the neutrals
    CupN   = CupN   + SUM(up_N(b_cor:e_cor))
    CdownN = CdownN + SUM(down_N(b_cor:e_cor))
    YNn    = YNn    + SUM(YN(b_cor:e_cor))
    Sn     = Sn     + SUM(YS(b_cor:e_cor))
    An     = An     + SUM(YA(b_cor:e_cor))
    Bn     = Bn     + 0.5_DP * SUM(YB(b_cor:e_cor)) + &
             SUM(up_DE(b_cor:e_cor)) * EVerg

    ! positive ions
    CupI   = SUM(up_N(b_ion:e_ion))
    CdownI = SUM(down_N(b_ion:e_ion))
    YNi    = SUM(YN(b_ion:e_ion))
    Si     = SUM(YS(b_ion:e_ion))
    Ai     = SUM(YA(b_ion:e_ion))
    Bi     = 0.5_DP * SUM(YB(b_ion:e_ion)) + SUM(up_DE(b_ion:e_ion)) * EVerg

    ! negative ions
    CupNEG   = SUM(up_N(b_neg:e_neg))
    CdownNEG = SUM(down_N(b_neg:e_neg))
    YNneg    = SUM(YN(b_neg:e_neg))
    Sneg     = SUM(YS(b_neg:e_neg))
    Aneg     = SUM(YA(b_neg:e_neg))
    Bneg     = 0.5_DP * SUM(YB(b_neg:e_neg)) + SUM(up_DE(b_neg:e_neg)) * EVerg
    ! heating of the electrons by photo-ionization
    Bneg     = Bneg + up_DE(ind_e) * EVerg

    ! change in mass density of grains by erosion
    Score    = SUM(YS(b_cor:e_cor))

  END SUBROUTINE CHEMISTRY



  SUBROUTINE SOURCE

    !---------------------------------------------------------------------------
    ! called by :
    !     DIFFUN
    ! purpose :
    !     First calculates the cooling of the different fluids, then computes
    !     the source terms needed in DIFFUN: momentum and energy. These terms
    !     are updated from their values calculated in CHEMISTRY. Note that
    !     source terms for number and mass density are not modified.
    !
    !     Source takes into account the cooling due to H2 lines, fine structure
    !     lines and other molecular cooling.
    !     It includes also other physical processes :
    !        - heat transfer between neutral and charged fluids;
    !        - photo-electric effect on grains;
    !        - thermalisation with grains;
    !        - momentum and energy transfer between neutral and charged fluids
    !          through elastic scattering;
    !        - energy transfer between charged fluids through
    !          Coulomb scattering.
    ! subroutine/function needed :
    !     LINE_THERMAL_BALANCE
    !     COMPUTE_MOLECULAR_COOLING
    !     COMPUTE_H2
    ! input variables :
    ! ouput variables :
    ! results :
    !     An,  Ai,  Aneg  -> source terms for momentum (g.cm-2.s-2)
    !     Bn,  Bi,  Bneg  -> source terms for energy (erg.cm-3.s-1)
    ! references :
    !     [1]  Flower et al., 1985, MNRAS 216, 775
    !     [2]  Flower et al., 1986, MNRAS 218, 729
    !     [3]  Pineau des Forets et al., 1986, MNRAS 220, 801
    !     [4]  Monteiro & Flower, 1987, MNRAS 228, 101
    !     [5]  Black, 1987, in Interstellar processes, eds. Hollenbach & !
    !          Thronson (Reidel, Dordrecht), P.731
    !     [6]  Viala et al., 1988, A&A 190, 215
    !     [7]  Spitzer, 1962, Physics of fully ionised gases (Wiley, New York)
    !     [8]  Spitzer & Scott, 1969, AP.J. 158, 161
    !     [9]  Flower & Watt, 1984, MNRAS 209, 25
    !     [10] Roueff 1990, A&A 234, 567
    ! modification history :
    !   - nov 2000 : David Wilgenbus (suggestion of David Flower)
    !       correction of the source term for momentum transfer between grains
    !       and neutrals.
    !       "A_grain_n = RhoN * 1.1D-21 * nH * ABS_DeltaV * DeltaV"
    !       is replaced by
    !       "A_grain_n = RhoN*Dens_grain *pi*Rsquare_GRAIN * ABS_DeltaV*DeltaV"
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_GAMMA, ONLY : GAMMA1
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_H2
    USE MODULE_GRAINS, ONLY : Tgrain, Dens_grain, Rsquare_GRAIN
    USE MODULE_DEBUG_JLB
    USE MODULE_LINE_EXCIT

    IMPLICIT NONE
    REAL(KIND=DP) :: cooling_n, cooling_i, cooling_neg ! cooling rate (erg/cm3/s)
    REAL(KIND=DP) :: Teff, Ve, alpha_n
    REAL(KIND=DP) :: muN_plus_muI, muN_plus_muNEG, RhoN_times_RhoI, muN_NEG
    REAL(KIND=DP) :: Sigma_IN, Sigma_eN, Sigma_NegN, Sigma_inelastic_eN
    REAL(KIND=DP) :: A_i_n, A_e_n, A_neg_n, A_grain_n
!   REAL(KIND=DP) :: B_i_n, B_neg_n, B_e_n, B_inelastic_e_n, B_grain_n
    REAL(KIND=DP) :: B_heat_exchange_n, B_photo_grain_n, B_therm_grain_n, B_ioniz_RC_n
    REAL(KIND=DP) :: lambda, B_i_e, B_i_neg, B_heat_exchange_i
    REAL(KIND=DP) :: B_heat_exchange_neg
    REAL(KIND=DP) :: Sel_net_H2

    !----------------------------------------------------------
    ! in FINE_STRUCTURE_COOLING, cooling rates (erg/cm3/s) for
    ! fine structure lines of C+, Si+, C, Si, O, S+, and N+
    ! are calculated.
    ! results are :
    ! cooling_Cplus, cooling_Siplus, cooling_C, cooling_Si
    ! cooling_O, cooling_Splus, cooling_Nplus, cooling_Cplus_e,
    ! cooling_Siplus_e, cooling_Splus_e, cooling_Nplus_e
    !----------------------------------------------------------
!   CALL FINE_STRUCTURE_COOLING

!  JLB - september 2002
!  FINE_STRUCTURE_COOLING is obsolete and superseded by :

    call LINE_THERMAL_BALANCE

    !---------------------------------------------------------
    ! in COMPUTE_MOLECULAR_COOLING, cooling rate (erg/cm3/s)
    ! of 13CO is calculated.
    ! useful result for SOURCE is :
    !    molec_cool
    !---------------------------------------------------------
!
! 29 mai 2001 - JLB : compute molec_cool directly in DIFFUN
!
!   CALL COMPUTE_MOLECULAR_COOLING
!  DRF - november 2012
!  COMPUTE_MOLECULAR_COOLING now calculates only the 13CO cooling rate;
!  the cooling rates for OH, NH3, CO and H2O are calculated explicitly
!  when solving for their level population densities (LVG approximation)

    !-----------------------------------------------------
    ! in COMPUTE_H2, cooling rate (erg/cm3/s) for
    ! ro-vibrational lines of H2 and source terms for
    ! H2 levels are calculated.
    ! results are :
    ! Cooling_H2, YN_rovib_H2, H2_Energy
    !-----------------------------------------------------
    CALL COMPUTE_H2

    !------------------------------------------------------------------
    ! Add different contributions to the cooling of the neutral fluid
    ! cooling_n (erg/cm3/s) :
    ! (1) fine structure lines (from FINE_STRUCTURE_COOLING)
    !     remark : we assume that, in collisions with C+ or Si+, 
    !              the kinetic energy is lost by the (light) neutrals, 
    !              not by the (heavy) ions; this is why the terms  
    !              cooling_Cplus and cooling_Siplus are here and
    !              not in the cooling of the positively charged fluid.
    ! (1) - bis : Superseded by LINE_THERMAL_BALANCE
    !     For each collision, energy loss is allocated accurately to
    !     each fluid.
    ! (2) molecular (not H2) cooling (from COMPUTE_MOLECULAR_COOLING)
    ! (3) H2 rovibrational cooling (from COMPUTE_H2)
    ! (4) dissociation of H2 by collisions with H, He and H2
    !     we assume a rate 8 times smaller for H2 than for H
    !     and 10 times smaller for He
    !     remark :
    !        this process is not taken into account in CHEMISTRY as DE=0._DP
    !        for endothermic reactions (see ENERGY_DEFECT_REACTION).
    !------------------------------------------------------------------
    cooling_n = Zero
    ! (1)
    cooling_n = cooling_n - Cool_n
    ! (2)
!
! 29 mai 2001 - JLB : added directly in DIFFUN
!
!   cooling_n = cooling_n + &
!        molec_cool
    ! (3)
    cooling_n = cooling_n + &
         cooling_H2
    ! (4) to be supressed if endothermic reaction do not have DE = Zero
    ! Use consistent rate of collisional dissociation by neutrals 
    ! (56641.590 is energy of highest H2 level (in K))
    Sel_net_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_ne_H2 &
                * (56641.590-H2_lev(1:NH2_lev)%energy))

    cooling_n = cooling_n - Sel_net_H2 * kB

!   cooling_n = cooling_n + &
!        Dens_H * Dens_H2 * 7.2D-22*EXP(-52000._DP/Tn) + &
!        Dens_H2* Dens_H2 * 9.0D-23*EXP(-52000._DP/Tn)

    !------------------------------------------------------------------------
    ! Add the different contributions to the rate of radiative cooling of the
    ! negatively charged fluid cooling_neg (erg/cm3/s) through excitation by
    ! electron collisions of :
    ! (1) O and C : 3P -> 1D; N 4S -> 2D (MENDOZA, 1983, PNE SYMPOSIUM)
    ! (2) C+, Si+, S+, and N+ (from FINE_STRUCTURE_COOLING)
    ! (1) and (2) - bis : this done in LINE_THERMAL_BALANCE now
    ! (3) rovibrational excitation of H2, rotational excitation of CO
    !-------------------------------------------------------------------------
    ! (1) and (2)
    cooling_neg = Zero
    cooling_neg = cooling_neg - Cool_e
    ! (3)
    Ve = SQRT(8.0_DP*kB*Te/(pi*me)) ! thermal velocity of e- (cm/s)
    cooling_neg = cooling_neg + &
         Dens_e *  Dens_H2 * EVerg * &
         Ve * (4.0D-18*EXP(-510._DP/Te)*0.044_DP + 4.0D-17*EXP(-6330._DP/Te)*0.544_DP)
    cooling_neg = cooling_neg + &
         Dens_e * Dens_CO * EVerg * Ve * 1.0D-15*EXP(-5.0_DP/Te)*0.00043_DP

    !-------------------------------------------------------------------------
    ! Add the different contributions to the rate of radiative cooling of the
    ! positively charged fluid, cooling_i (erg/cm3/s) :
    ! (1) ro-vibrational excitation of H2 by collisions with ions
    ! (2) dissociation of H2 by collisions with ions
    !     remark :
    !        this is not taken into account in CHEMISTRY, as DE=Zero
    !        for endothermic reactions (see ENERGY_DEFECT_REACTION).
    ! (3) ionization of H2 by collisions with ions
    ! The energy loss by the positively charged fluid is proportional to
    !     mass_H2/(mass_H2 + muI); the corresponding loss by the neutral
    ! fluid is proportional to
    !     muI/(mass_H2 + muI)
    ! (4) ionization of H by collisions with ions
    ! The energy loss by the positively charged fluid is proportional to
    !     mass_H/(mass_H + muI); the corresponding loss by the neutral
    ! fluid is proportional to
    !     muI/(mass_H + muI)
    !-------------------------------------------------------------------------
    cooling_i = Zero
    ! (1)
    Teff = & ! effective temperature, used in (1), (2) and (3)
         (mass_H2*muI*ABS_DeltaV*ABS_DeltaV/3.0_DP/kB + (mass_H2*Ti + muI*Tn)) / &
         (mass_H2 + muI)
    cooling_i = cooling_i + &
         ABS_DeltaV * EVerg * Dens_H2 * DensityI * &
         (4.0D-16*EXP(-510._DP/Teff)*0.044_DP + 4.0D-17*EXP(-6330._DP/Teff)*0.544_DP)
    ! (2) to be supressed if endothermic reaction do not have DE = Zero
    Sel_net_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_io_H2 &
                * (56641.590-H2_lev(1:NH2_lev)%energy))
    cooling_i = cooling_i - Sel_net_H2 * kB
    ! (3) to be supressed if endothermic reaction do not have DE = Zero
    cooling_i = cooling_i + &
         DensityI * Dens_H2 * 1.10D-13*(Teff/300._DP)**0.5_DP* &
         EXP(-179160.0_DP/Teff)*15.43_DP * EVerg
    cooling_i = cooling_i * mass_H2/(mass_H2 + muI)
    cooling_n = cooling_n + cooling_i * muI/mass_H2
    ! (4) to be supressed if endothermic reaction do not have DE = Zero
    Teff = & ! effective temperature used in (4)
         (mass_H*muI*ABS_DeltaV*ABS_DeltaV/3.0_DP/kB + (mass_H*Ti + muI*Tn)) / &
         (mass_H + muI)
    cooling_i = cooling_i + &
         DensityI * Dens_H * 1.30D-13*(Teff/300._DP)**0.5_DP* &
         EXP(-157890.0_DP/Teff)*13.60_DP * EVerg * mass_H/(mass_H + muI)
    cooling_n = cooling_n + &
         DensityI * Dens_H * 1.30D-13*(Teff/300._DP)**0.5_DP* &
         EXP(-157890.0_DP/Teff)*13.60_DP * EVerg * muI/(mass_H + muI)

    ! (0) add contribution from LINE_THERMAL_BALANCE
    ! Added after (1)-(4) so as not to add it also to cooling_n
    cooling_i = cooling_i - Cool_i

    !--- some useful variables ---
    muN_plus_muI    = muN  + muI           ! sum of neutral and positive ion mass (g)
    muN_plus_muNEG  = muN  + muNEG         ! sum of neutral and negative ion mass (g)
    RhoN_times_RhoI = RhoN * RhoI          ! product of neutral and positive ion mass densities (g2/cm6)
    muN_NEG = muN * muNEG / muN_plus_muNEG ! reduced mass of a neutral/negative ion pair (g)

    !--- polarisability of the neutrals                    ---
    !--- = weighted average of polarisabilities of H and H2  ---
    alpha_n=(Dens_H*alpha_H + Dens_H2*alpha_H2) / (Dens_H+Dens_H2)

    !------------------------------------------------------------------
    ! calculate the cross-sections (cm-2)
    ! (1) Sigma_IN           -> elastic ions > 0  - neutral scattering
    !                           ([1], eq.(23))
    ! (2) Sigma_eN           -> elastic electron - neutral scattering
    !                           ([1], eq.(34))
    ! (3) Sigma_NegN         -> elastic ions < 0  - neutral scattering
    !                           ([1], eq.(23))
    ! (4) Sigma_inelastic_eN -> inelastic electron-neutral scattering
    !------------------------------------------------------------------
    ! (1)
    Sigma_IN = 2.41_DP * pi * qe * SQRT(muN_plus_muI*alpha_n/muN/muI)
    Sigma_IN = MAX(Sigma_IN, 1.0D-15*ABS_DeltaV)
    ! (2)
    Sigma_eN = 1.0D-15 * SQRT(8.0_DP*kB*Te/(pi*me))
    ! (3)
    Sigma_NegN = 2.41_DP * pi * qe * SQRT(alpha_n/muN_NEG)
    Sigma_NegN = MAX(Sigma_NegN, 1.0D-14*ABS_DeltaV)
    ! (4)
    Sigma_inelastic_eN = 1.0D-16 * (4.29D-1 + 6.14D-6*Te) * SQRT(8.0_DP*kB*Te/(pi*me))

    !--------------------------------------------------------------------
    ! Rates of production of momentum (g.cm-2.s-2) in the neutral fluid
    !   A_i_n     -> elastic ions > 0  - neutral scattering
    !   A_e_n     -> elastic electron - neutral scattering
    !   A_neg_n   -> elastic ions < 0  - neutral scattering
    !   A_grain_n -> elastic charged grains - neutral scattering
    !--------------------------------------------------------------------
    A_i_n     = Sigma_IN * DeltaV * RhoN_times_RhoI / muN_plus_muI
    A_e_n     = DensityN * Dens_e * me * Sigma_eN * DeltaV
    A_neg_n   = DensityN * DensityNEG * muN_NEG * Sigma_NegN * DeltaV
    A_grain_n = RhoN * Dens_grain * pi * Rsquare_GRAIN * ABS_DeltaV * DeltaV &
                * (Dens_Gminus + Dens_Gplus)/ (Dens_Gminus + Dens_Gplus + Dens_G)
!   A_grain_n = RhoN * Dens_grain * pi * Rsquare_GRAIN * ABS_DeltaV * DeltaV

    !------------------------------------------------------------------
    ! Source terms of number and mass density are calculated in
    ! CHEMISTRY and not in SOURCE.
    ! These terms are (for neutral fluid, positive fluid, negative fluid) :
    ! YNn, YNi, YNneg -> change in number density (cm-3.s-1)
    ! Sn,  Si,  Sneg  -> change in mass density (g.cm-3.s-1)
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Source terms of momentum (g.cm-2.s-2) of each fluid (update of
    ! the calculation in CHEMISTRY).
    !   An   -> neutral fluid
    !   Ai   -> positive fluid
    !   Aneg -> negative fluid
    !----------------------------------------------------------------
    An   = An   + A_i_n + A_e_n + A_neg_n + A_grain_n
    Ai   = Ai   - A_i_n
    Aneg = Aneg         - A_e_n - A_neg_n - A_grain_n

    !--------------------------------------------------------------------------
    ! Source term of energy Bn (erg.cm-3.s-1) in the neutral fluid
    ! obtained by adding :
    ! +Bn                 -> exo/endothermicity of the reactions
    !                        (from CHEMISTRY)
    ! -cooling_n          -> energy loss through cooling
    ! +B_i_n              -> elastic diffusion of the neutrals by the
    !                        positive ions ([1], eq.(32))
    ! +B_neg_n            -> elastic diffusion of the neutrals by the
    !                        negative ions ([1], eq.(32))
    ! +B_e_n              -> elastic diffusion of the neutrals by the
    !                        electrons ([1], eq.(33))
    ! +B_inelastic_e_n    -> inelastic scattering of the electrons on H2;
    !                        mean excitation energy 12 eV = 140 000 K for the
    !                        singlet excited states, 10 eV = 116 300 K for the
    !                        triplet, of which 10 - 4.48 = 5.52 eV are recovered
    !                        by the neutrals through the dissociation products;
    !                        collisional excitation of the n=2 state of H.
    ! +B_grain_n          -> elastic diffusion of the neutrals by the
    !                        charged grains
    ! +B_heat_exchange_n  -> heat exchange
    !                        [1],eq.(27) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1)
    ! +B_photo_grain_n    -> heating through the photo-electric effect on grains
    !                        [5]
    ! +B_therm_grain_n    -> energy loss/gain through thermalisation with grains
    !                        TIELENS  & HOLLENBACH (1987, AP.J. 291, 722)
    ! +B_ioniz_RC_n       -> heating through cosmic ray ionization [8]
    !                        (not included in CHEMISTRY as DE = Zero for
    !                        reactions of the type 'CR_IO')
    !--------------------------------------------------------------------------

    B_i_n = (3.0_DP*kB*(Ti-Tn) + DeltaV*(muI*Vi+muN*Vn)) * &
         RhoN_times_RhoI/muN_plus_muI * Sigma_IN/muN_plus_muI

    B_neg_n = (3.0_DP*kB*(Te-Tn) + DeltaV*(muNEG*Vi+muN*Vn)) * &
         DensityN*DensityNEG * Sigma_NegN*muN_NEG/muN_plus_muNEG

    B_e_n = (4.0_DP*kB*(Te-Tn) + muN*Vn*DeltaV) * &
         DensityN*Dens_e * me/muN*Sigma_eN

    B_inelastic_e_n = - SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_el_H2) * &
                      5.52_DP*EVerg

    B_grain_n = A_grain_n*Vi

    B_heat_exchange_n = GAMMA1*kB*(CupN*(Ti+Te)/2._DP - CdownN*Tn)

    B_photo_grain_n = 4.0D-26*nH*RAD*EXP(-2.5_DP*AV)

    B_therm_grain_n = 3.5D-34*SQRT(Tn)*(Tgrain-Tn)*nH**2.0_DP

    B_ioniz_RC_n = Zeta * EVerg * DensityN * &
         MAX(5.7_DP,32.0_DP-7.0_DP*LOG10(DensityN/Dens_e))

    !--- sum of the various contributions to Bn ---
    Bn = Bn                  &
         - cooling_n         &
         + B_i_n             &
         + B_neg_n           &
         + B_e_n             &
         + B_inelastic_e_n   &
         + B_grain_n         &
         + B_photo_grain_n   &
         + B_therm_grain_n   &
         + B_ioniz_RC_n      !&
         !+ B_heat_exchange_n ! include only Bn and not B_heat_exchange_n (March 2000)

    !--------------------------------------------------------------------------
    ! Source term of energy Bi (erg.cm-3.s-1) for the positive ions fluid
    ! obtained by adding :
    !  +Bi                -> exo/endothermicity of the reactions
    !                        (from CHEMISTRY)
    !  -cooling_i         -> energy loss through cooling
    !  -B_i_e             -> diffusion of electrons by positive ions
    !                        [7], see also [1], eq.(35)
    !  -B_i_neg           -> diffusion of negative ions by positive ions
    !                        [7]
    !  -B_i_n             -> elastic diffusion of the neutrals by the
    !                        positive ions. See the calculation for Bn above
    !  +B_heat_exchange_i -> heat exchange
    !                        [1],eq.(28) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1)
    !--------------------------------------------------------------------------
    lambda = 1.5_DP/(qe**3.0_DP) * SQRT((kB*Te)**3.0_DP/(Dens_e*pi))
    B_i_e = 4.0_DP*qe**4._DP/muI/Te * SQRT(2.0_DP*me*pi/(kB*Te)) * &
         DensityI*Dens_e*LOG(lambda)*(Ti-Te)

    lambda = 1.5_DP/(qe**3.0_DP) * SQRT((kB*Te)**3.0_DP/(DensityI*pi))
    B_i_neg = qe**4.0_DP/(muI*muNEG) * (Ti/muI+Te/muNEG)**(-1.5_DP)* &
         (32.0_DP*pi/kB)**0.5_DP*DensityI*DensityNEG*LOG(lambda)*(Ti-Te)

    B_heat_exchange_i=GAMMA1*kB*(CupI*Tn-CdownI*Ti)

    !--- sum of the different contributions to Bi ---
    Bi = Bi                   &
         - cooling_i          &
         - B_i_e              &
         - B_i_neg            &
         - B_i_n              &
         + B_heat_exchange_i

    !--------------------------------------------------------------------------
    ! CALCULATE THE SOURCE TERM OF ENERGY Bneg (ERG.CM-3.S-1) FOR THE
    ! NEGATIVE FLUID BY ADDING :
    !   B_i_e              -> diffusion of electrons by positive ions
    !                         [7], see also [1], eq.(35).
    !                         see calculation for Bi above.
    !   B_i_neg            -> diffusion of negative ions by positive ions
    !                         [7], see calculation for Bi above.
    !  -B_e_n              -> elastic diffusion of the neutrals by the
    !                         electrons ([1], eq.(33));
    !                         see calculation for Bn above.
    !  -B_inelastic_e_n    -> inelastic scattering of the electrons on H2 and H.
    !                         For H, excitation of the n=2 state :
    !                         Aggarwal et al. 1991, J. PHYS. B, 24, 1385
    !                                ionization:
    !                         Hummer, 1963, MNRAS, 125, 461.
    !                         For H2, excitation of the repulsive B triplet state
    !                         and of the B and C singlet states
    !                                ionization:
    !                         Rapp & Englander-Golden, 1965, JCP, 43, 1464.
    !  -B_grain_n          -> elastic diffusion of the neutrals by the
    !                         charged grains, see calculation for Bn above.
    !  -B_neg_n            -> elastic diffusion of the neutrals by the
    !                         negative ions ([1], eq.(32));
    !                         see calculation for Bn above.
    !  cooling_neg         -> energy loss through cooling
    !  B_heat_exchange_neg -> heat exchange
    !                        [1],eq.(28) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1)
    !--------------------------------------------------------------------------
    B_inelastic_e_n = Sigma_inelastic_eN*EXP(-1.4D5/Te)
    Sel_net_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_el_H2 &
                * (116300.0_DP-H2_lev(1:NH2_lev)%energy))
    B_inelastic_e_n = Dens_H2*Dens_e * &
                      B_inelastic_e_n * 12.0_DP * EVerg - &
                      Sel_net_H2 * kB
    B_inelastic_e_n = B_inelastic_e_n + &
         Dens_H*Dens_e*EVerg*3.8D-06*EXP(-1.1842D5/Te)*(Te/300._DP)**(-0.5_DP)
    B_inelastic_e_n = B_inelastic_e_n + &
         Dens_H*Dens_e*EVerg*13.6_DP*9.2D-10*EXP(-1.5789D5/Te)*(Te/300.0D0)**0.5_DP + &
         Dens_H2*Dens_e*EVerg*15.43_DP*1.40D-09*EXP(-1.7916D5/Te)*(Te/300._DP)**0.5_DP

    B_heat_exchange_neg=0.5_DP*me*Vi*Vi*(YNi-YNneg)-GAMMA1*kB*Te*CdownI

    !--- sum of the different contributions to Bneg ---
    Bneg = Bneg                 &
         - cooling_neg          &
         + B_i_e                &
         + B_i_neg              &
         - B_e_n                &
         - B_inelastic_e_n      &
         - B_grain_n            &
         - B_neg_n              &
         + B_heat_exchange_neg

  END SUBROUTINE SOURCE



  SUBROUTINE DIFFUN(N, Z, Y, DY)

    !---------------------------------------------------------------------------
    ! called by :
    !    STIFF
    !    PSET
    ! purpose :
    !    calculates DY, the partial derivatives of the Y vector.
    !               DY(i)/DZ= F(Y,Z), 1 <= i <= N
    !    see [1], equations (3)-(14)
    !    The source terms (change per unit volume and time) are calculated in
    !    CHEMISTRY and SOURCE.
    !    DIFFUN also extracts all physical variables from the Y vector
    !    (these variables are used in SOURCE).
    ! subroutine/function needed :
    !    COMPUTE_OP_H2
    !    CHEMISTRY
    !    SOURCE
    ! input variables :
    !    N -> dimension of Y and DY
    !    Z -> independent variable (distance)
    !    Y -> vector containing MHD variables + species + para-H2 + H2 levels
    !         + level populations of 'heavy' molecules (CO, SiO, H2O, NH3, OH, CH3OH)
    ! output variables :
    !    DY -> first order derivatives DY/dZ
    ! results :
    !    v_l_var, v_l_der ...
    ! references :
    !   [1] Flower et al., 1985, MNRAS 216, 775
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_GAMMA
    USE MODULE_GRAINS
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_H2
    USE MODULE_CONSTANTS, ONLY : kB, pi, Zero, parsec
!   USE MODULE_PHYS_VAR,  ONLY : NH2_lev
    USE MODULE_DEBUG_JLB
    USE mco
    USE msio
    USE moh2o
    USE mph2o
    USE monh3
    USE mpnh3
    USE moh
    USE m_atype
    USE m_etype
    IMPLICIT NONE
    INTEGER :: NCO_lev, NSiO_lev, NoH2O_lev, NpH2O_lev, NoNH3_lev, NpNH3_lev, NOH_lev  &
               , Natype_lev, Netype_lev &
               ,ICOUNT_CO, ICOUNT_SiO, ICOUNT_oH2O, ICOUNT_pH2O, ICOUNT_oNH3, ICOUNT_pNH3 &
               ,ICOUNT_OH, ICOUNT_atype, ICOUNT_etype
! The values of the following PARAMETERS must be the same as those of the 
! PARAMETER (NLEVCO) in the corresponding LVG modules (CO.f, SiO.f, ... )
! See also mhd_vode.f90 and outputs.f90
    PARAMETER (NCO_lev=41)
    PARAMETER (NSiO_lev=41)
    PARAMETER (NoH2O_lev=45)
    PARAMETER (NpH2O_lev=45)
    PARAMETER (NoNH3_lev=17)
    PARAMETER (NpNH3_lev=24)
    PARAMETER (NOH_lev=20)
    PARAMETER (Natype_lev=256)
    PARAMETER (Netype_lev=256)
    INTEGER(KIND=LONG), INTENT(in) :: N ! dimension of Y and DY
    REAL(KIND=DP),INTENT(in) :: Z ! useless here : the dependence in Z is implicit
!   REAL(KIND=DP),DIMENSION(N), INTENT(in) :: Y ! vector containing the log of the variables
    REAL(KIND=DP),DIMENSION(N), INTENT(inout) :: Y ! vector containing the log of the variables
    REAL(KIND=DP),DIMENSION(N), INTENT(inout) :: DY ! vector containing the derivatives
    REAL(KIND=DP) :: Bfield2, kT_muN, Vi2, Vn2, V_DensityI, Vmagson ! useful for the calculations
    INTEGER(KIND=LONG) :: b_specy ! useful for densities of chemical species
    REAL(KIND=DP) :: aux1_dvn, aux2_dvn, dvn1, dvn2, gvn1, gvn2, corrdv
    INTEGER :: ii, i
    REAL(KIND=DP) :: Rgrain_cube ! useful for computing R_square
    REAL(KIND=DP) :: XLL2, Bne, Ane, Sne, YNne, DensityNe, RhoNe
    REAL(KIND=DP) :: dens_CO_lev(1:NCO_lev), YN_rovib_CO(1:NCO_lev), sum_CO
    REAL(KIND=DP) :: dens_SiO_lev(1:NSiO_lev), YN_rovib_SiO(1:NSiO_lev), sum_SiO
    REAL(KIND=DP) :: dens_oH2O_lev(1:NoH2O_lev), YN_rovib_oH2O(1:NoH2O_lev), sum_oH2O, op_H2O
    REAL(KIND=DP) :: dens_pH2O_lev(1:NpH2O_lev), YN_rovib_pH2O(1:NpH2O_lev), sum_pH2O
    REAL(KIND=DP) :: dens_oNH3_lev(1:NoNH3_lev), YN_rovib_oNH3(1:NoNH3_lev), sum_oNH3, op_NH3
    REAL(KIND=DP) :: dens_pNH3_lev(1:NpNH3_lev), YN_rovib_pNH3(1:NpNH3_lev), sum_pNH3
    REAL(KIND=DP) :: dens_OH_lev(1:NpNH3_lev), YN_rovib_OH(1:NpNH3_lev), sum_OH
    REAL(KIND=DP) :: dens_atype_lev(1:Natype_lev), YN_rovib_atype(1:Natype_lev), sum_atype
    REAL(KIND=DP) :: dens_etype_lev(1:Netype_lev), YN_rovib_etype(1:Netype_lev), sum_etype, ae_CH3OH
    REAL(KIND=DP) :: TAU_CO(NCO_lev,NCO_lev)
    REAL(KIND=DP) :: TAU_SiO(NSiO_lev,NSiO_lev)
    REAL(KIND=DP) :: freq_oH2O(NoH2O_lev,NoH2O_lev), freq_pH2O(NpH2O_lev,NpH2O_lev) &
                   , freq_oNH3(NoNH3_lev,NoNH3_lev), freq_pNH3(NpNH3_lev,NpNH3_lev) &
                   , freq_OH(NOH_lev,NOH_lev) &
                   , freq_atype(Natype_lev,Natype_lev), freq_etype(Netype_lev,Netype_lev) &
                   , TAU_oH2O(NoH2O_lev,NoH2O_lev), TAU_pH2O(NpH2O_lev,NpH2O_lev) &
                   , TAU_oNH3(NoNH3_lev,NoNH3_lev), TAU_pNH3(NpNH3_lev,NpNH3_lev) & 
                   , TAU_OH(NOH_lev,NOH_lev) &
                   , TAU_atype(Natype_lev,Natype_lev), TAU_etype(Netype_lev,Netype_lev)
    REAL(KIND=DP) :: WWT_CO, WWT_SiO, WWT_oH2O, WWT_pH2O, WWT_oNH3, WWT_pNH3, WWT_OH, WWT_atype, WWT_etype
    REAL(KIND=DP) :: timescale_CO(1:NCO_lev),timescale_SiO(1:NSiO_lev),timescale_oH2O(1:NoH2O_lev) &
                    ,timescale_pH2O(1:NpH2O_lev),timescale_oNH3(1:NoNH3_lev) &
                    ,timescale_pNH3(1:NpNH3_lev),timescale_OH(1:NOH_lev),timescale_atype(1:Natype_lev) &
                    ,timescale_etype(1:Netype_lev),dyntime
      COMMON /FCTCO/ dens_CO_lev,TAU_CO,WWT_CO
      COMMON /FCTSiO/ dens_SiO_lev,TAU_SiO,WWT_SiO
      COMMON /FCToH2O/ dens_oH2O_lev,freq_oH2O,TAU_oH2O,WWT_oH2O
      COMMON /FCTpH2O/ dens_pH2O_lev,freq_pH2O,TAU_pH2O,WWT_pH2O
      COMMON /FCToNH3/ dens_oNH3_lev,freq_oNH3,TAU_oNH3,WWT_oNH3
      COMMON /FCTpNH3/ dens_pNH3_lev,freq_pNH3,TAU_pNH3,WWT_pNH3
      COMMON /FCTOH/ dens_OH_lev,freq_OH,TAU_OH,WWT_OH
      COMMON /FCTa_CH3OH/ dens_atype_lev,freq_atype,TAU_atype,WWT_atype
      COMMON /FCTe_CH3OH/ dens_etype_lev,freq_etype,TAU_etype,WWT_etype
      COMMON /tscales/ timescale_CO,timescale_SiO,timescale_oH2O,timescale_pH2O,timescale_oNH3 &
                      ,timescale_pNH3,timescale_OH,timescale_atype,timescale_etype,dyntime
      DATA ICOUNT_CO,ICOUNT_SiO,ICOUNT_oH2O,ICOUNT_pH2O,ICOUNT_oNH3,ICOUNT_pNH3,ICOUNT_OH &
           ,ICOUNT_atype,ICOUNT_etype / 0, 0, 0, 0, 0, 0, 0, 0, 0 /
! THE ORTHO:PARA H2O RATIO
      DATA op_H2O / 3.D0 /
! THE ORTHO:PARA NH3 RATIO
      DATA op_NH3 / 1.D0 /
! THE A:E CH3OH RATIO
      DATA ae_CH3OH / 1.D0 /

    !-------------------------------------
    ! convert from logarithmic variables
    ! (Y = v_lvariab)
    !-------------------------------------
    v_variab(1:N)=EXP(Y(1:N))
    ! if v_variab <= Zero (unphysical for these variables) for this interpolation of DVODE
    ! use the last value v_l_var
    WHERE (v_variab <= Zero)
       v_variab=v_l_var
    END WHERE
    v_dvariab(1:N) = v_variab(1:N) * DY(1:N)

    b_specy = bv_speci-1
!  JLB test - 25 Jan. 2002
    sum_H2 = SUM(v_variab(bv_H2_lev:ev_H2_lev))
    if (Z > 1.0e6) then
      Y(b_specy+ind_H2) = log(sum_H2)
    endif
    v_variab(b_specy+ind_H2) = sum_H2
!  End test

!  JLB test - 19 Feb. 2002
!  ... similar type
    if (Force_I_C == 1) then
      DensityI = SUM(v_variab(bv_ion:ev_ion))
      if (Z > 1.0e6) then
        Y(iv_DensityI) = log(DensityI)
      endif
      v_variab(iv_DensityI) = DensityI
    endif
!  End test

    !--------------------------------------------------------------------
    ! Saveguard against round-off errors in the derivation of Vi, Ti, Te
    ! for one fluid model, we have :
    !    Vn=Vi (velocities), Tn=Ti=Te (temperatures)
    ! for two fluids model, we have :
    !    Ti=Te (temperatures)
    ! Can be checked in printouts periodically performed in MHD
    !----------------------------------------------------------------------
    SELECT CASE (Nfluids)
    CASE (1)
       v_variab(iv_Vi) = v_variab(iv_Vn) ! Vi=Vn
       v_variab(iv_Ti) = v_variab(iv_Tn) ! Ti=Tn
       v_variab(iv_Te) = v_variab(iv_Tn) ! Te=Tn
    CASE (2)
       v_variab(iv_Te) = v_variab(iv_Ti) ! Te=Ti
    CASE DEFAULT ! do nothing for three-fluid model
    END SELECT

    !------------------------------------------
    !--- extracts physical variables        ---
    !--- (module MODULE_PHYS_VAR) ---
    !------------------------------------------

    !--- temperatures (K) : neutrals, ions > 0, electrons and ions < 0 ---
    Tn = v_variab(iv_Tn)
    Ti = v_variab(iv_Ti)
    Te = v_variab(iv_Te)

    !--- velocities (cm/s) : neutrals, ions, ion-neutral drift velocity ---
    Vn = v_variab(iv_Vn)
    Vi = v_variab(iv_Vi)
    DeltaV     = Vi - Vn
    ABS_DeltaV = ABS(DeltaV)

    if ((shock_type == "J") .AND. (viscosity)) then
      grad_V = v_variab(iv_gv)
    endif

    !--- number densities (cm-3) : neutrals, ions > 0, ions < 0, species on grains ---
    DensityN = v_variab(iv_DensityN)
    DensityI = v_variab(iv_DensityI)
    DensityNEG = SUM(v_variab(bv_neg:ev_neg))
    Dens_ongrains = SUM(v_variab(bv_gra:ev_gra))
    IF (Dens_ongrains < 1.0D-11) Dens_ongrains = 1.0D-11       ! pourquoi ????????????????
    Dens_cor = SUM(v_variab(bv_cor:ev_cor))

    !--- mass densities (g/cm3) : neutrals, ions > 0, ions < 0 ---
    RhoN = v_variab(iv_RhoN)
    RhoI = v_variab(iv_RhoI)
    RhoNEG = v_variab(iv_RhoNEG)

    !--- average masses (g) : neutrals, ions > 0, ions < 0 ---
    muN = RhoN / DensityN
    muI = RhoI / DensityI
    muNEG = RhoNEG / DensityNEG

    !--- sound speed and ion magnetosonic speed (cm.s-1) ---
    Vsound = SQRT(Gamma*kB*Tn/muN)
    Vmagnet = Gamma * kB * (Ti+Te) / &
         ((DensityI*muI+DensityNEG*muNEG) / (DensityI+DensityNEG)) + &
         (Bfield*Vs_cm/Vi)**2._DP / (4._DP*pi*RhoI)
    Vmagnet = SQRT(Vmagnet)
    !--- magnetosonic speed (cm.s-1) ---
    Vmagson = Gamma * kB * (Tn+Ti+Te) / &
         ((DensityN*muN+DensityI*muI+DensityNEG*muNEG) / (DensityN+DensityI+DensityNEG)) + &
         (Bfield*Vs_cm/Vi)**2._DP / (4._DP*pi*(RhoN+RhoI+RhoNEG))
    Vmagson = SQRT(Vmagson)

    !--- variables useful for calculations ---
    kT_muN = kB * Tn / muN
    Vn2 = Vn * Vn
    Vi2 = Vi * Vi
    V_DensityI = Vi * DensityI
    ! square of the magnetic field strength
    Bfield2 = (Vs_cm*Vs_cm) * (Bfield*Bfield) / (4.0_DP*pi*Vi2)
    XLL2 = XLL * XLL

    !----------------------------------------------------------------
    !--- density of some species (module MODULE_CHEMICAL_SPECIES) ---
    !---        used in SOURCE and in RSF                         ---
    !--- for e-, see later                                        ---
    !----------------------------------------------------------------
    Dens_H = v_variab(b_specy+ind_H)
      IF (ind_H==0) Dens_H = Zero            ! if the species is not present, its density is Zero
    Dens_H2 = v_variab(b_specy+ind_H2)
      IF (ind_H2==0) Dens_H2 = Zero
    Dens_He = v_variab(b_specy+ind_He)
      IF (ind_He==0) Dens_He = Zero
    Dens_O = v_variab(b_specy+ind_O)
      IF (ind_O==0) Dens_O = Zero
    Dens_Oplus = v_variab(b_specy+ind_Oplus)
      IF (ind_Oplus==0) Dens_Oplus = Zero
    Dens_N = v_variab(b_specy+ind_N)
      IF (ind_N==0) Dens_N = Zero
    Dens_C = v_variab(b_specy+ind_C)
      IF (ind_C==0) Dens_C = Zero
    Dens_S = v_variab(b_specy+ind_S)
      IF (ind_S==0) Dens_S = Zero
    Dens_Si = v_variab(b_specy+ind_Si)
      IF (ind_Si==0) Dens_Si = Zero
    Dens_SiO = v_variab(b_specy+ind_SiO)
      IF (ind_SiO==0) Dens_SiO = Zero
    Dens_H2O = v_variab(b_specy+ind_H2O)
      IF (ind_H2O==0) Dens_H2O = Zero
    Dens_OH = v_variab(b_specy+ind_OH)
      IF (ind_OH==0) Dens_OH = Zero
    Dens_CO = v_variab(b_specy+ind_CO)
      IF (ind_CO==0) Dens_CO = Zero
    Dens_NH3 = v_variab(b_specy+ind_NH3)
      IF (ind_NH3==0) Dens_NH3 = Zero
    Dens_CH3OH = v_variab(b_specy+ind_CH3OH)
      IF (ind_CH3OH==0) Dens_CH3OH = Zero
    Dens_G = v_variab(b_specy+ind_G)
      IF (ind_G==0) Dens_G = Zero
    Dens_Cplus = v_variab(b_specy+ind_Cplus)
      IF (ind_Cplus==0) Dens_Cplus = Zero
    Dens_Siplus = v_variab(b_specy+ind_Siplus)
      IF (ind_Siplus==0) Dens_Siplus = Zero
    Dens_Hplus = v_variab(b_specy+ind_Hplus)
      IF (ind_Hplus==0) Dens_Hplus = Zero
    Dens_Splus = v_variab(b_specy+ind_Splus)
      IF (ind_Splus==0) Dens_Splus = Zero
    Dens_Nplus = v_variab(b_specy+ind_Nplus)
      IF (ind_Nplus==0) Dens_Nplus = Zero
    Dens_Feplus = v_variab(b_specy+ind_Feplus)
      IF (ind_Feplus==0) Dens_Feplus = Zero
    Dens_Gplus = v_variab(b_specy+ind_Gplus)
      IF (ind_Gplus==0) Dens_Gplus = Zero
    Dens_Gminus = v_variab(b_specy+ind_Gminus)
      IF (ind_Gminus==0) Dens_Gminus = Zero

    !--- 'proton' density (cm-3) nH=n(H)+2n(H2)+n(H+) ---
    nH = Dens_H + 2._DP * Dens_H2 + Dens_Hplus

    !---------------------------------------------
    !--- chemical species (used in CHEMISTRY)  ---
    !--- density, velocity, temperature        ---
    !---------------------------------------------

    !--- ordinary species ---
    speci(1:Nspec)%density = v_variab(bv_speci:ev_speci)
 ! neutrals : V = Vn, T = Tn
    speci(b_neu:e_neu)%velocity    = Vn
    speci(b_neu:e_neu)%temperature = Tn
 ! species on grain mantles : V = Vi, T = Tgrain (not used)
    speci(b_gra:e_gra)%velocity    = Vi
    speci(b_gra:e_gra)%temperature = Tgrain
 ! species in grain cores : V = Vi, T = Tgrain (not used)
    speci(b_cor:e_cor)%velocity    = Vi
    speci(b_cor:e_cor)%temperature = Tgrain
 ! positive ions : V = Vi, T = Ti
    speci(b_ion:e_ion)%velocity    = Vi
    speci(b_ion:e_ion)%temperature = Ti
 ! negative ions : V = Vi, T = Te
    speci(b_neg:e_neg)%velocity    = Vi
    speci(b_neg:e_neg)%temperature = Te

 !--- added species ---
 ! electrons, with charge neutrality, and V = Vi, T = Te
    speci(ind_e)%density     = DensityI - DensityNEG
    Dens_e                   = speci(ind_e)%density
    speci(ind_e)%velocity    = Vi
    speci(ind_e)%temperature = Te
 ! photons, cosmic rays and secondary photons
 ! density=cst=1.0_DP, V=cst=Zero, T=cst=Zero (READ_SPECIES) => do nothing
 ! grains, compressed with ions, and V = 0, T = 0
    speci(ind_GRAIN)%density = Dens_grain_init * Vs_cm / Vi
    Dens_GRAIN = speci(ind_GRAIN)%density
 ! V=cst=Zero, T=cst=Zero (READ_SPECIES) => do nothing

    !--- update grain mass density (g/cm3), average radius (cm), square radius (cm2)
! SC apr 06: use only cores to compute grain radius, add mantles afterwards.
    MD_grain = DOT_PRODUCT(v_variab(bv_cor:ev_cor),DBLE(speci(b_cor:e_cor)%mass))
    Rgrain_cube = 3.0_DP * MD_grain / (4.0_DP * PI * Rho_GRAIN * Dens_GRAIN)

    Rsquare_GRAIN = R_gr_scale * Rgrain_cube**(2.0_DP/3.0_DP)

! number of sites per grain
    Nsites_grain = 4._DP*PI*Rsquare_grain / (dsites_grain * dsites_grain)

! SC add contribution of mantle thickness to grain cross section: 
    r_grain = R_gr_scale12 * sqrt(Rsquare_GRAIN) ! mean radius of grain cores
    Nlayers_grain = Dens_ongrains/Nsites_grain/Dens_GRAIN
    Rsquare_GRAIN = Rsquare_GRAIN + 2._DP * Nlayers_grain * layer_thickness * r_grain &
          + (Nlayers_grain * layer_thickness)**2 

    r_grain = sqrt(Rsquare_GRAIN)

    MD_grain = MD_grain  + DOT_PRODUCT(v_variab(bv_gra:ev_gra),DBLE(speci(b_gra:e_gra)%mass))

!  mass of a grain (same for neutral and charged grains)
!            speci(ind_G)%mass = MD_grain/Dens_GRAIN
!            speci(ind_Gplus)%mass = MD_grain/Dens_GRAIN
!            speci(ind_Gminus)%mass = MD_grain/Dens_GRAIN


!  mass of charged grains
    Rho_GRAIN_charge = MD_grain * (Dens_Gminus + Dens_Gplus)/ &
                       (Dens_Gminus + Dens_Gplus + Dens_G)

    !--- recompute magnetosonic speed (cm.s-1) allowing for charged grains ---
    Vmagnet = Gamma * kB * (Ti+Te) / &
         ((DensityI*muI+DensityNEG*muNEG) / (DensityI+DensityNEG)) + &
         (Bfield*Vs_cm/Vi)**2._DP / (4._DP*pi*(RhoI+RhoNEG+Rho_GRAIN_charge))
!        (Bfield*Vs_cm/Vi)**2._DP / (4._DP*pi*(RhoI+RhoNEG+MD_grain))
    Vmagnet = SQRT(Vmagnet)

    !--------------------------
    ! density of each H2 level
    !--------------------------
    H2_lev(1:NH2_lev)%density = v_variab(bv_H2_lev:ev_H2_lev)

! JLB + GF (25 VI 2002) - Update H2_lev%Form_gr if required

    if (iforH2 == 4) then
      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%density / Dens_H2
      H2_int_E = SUM(H2_lev(1:NH2_lev)%Form_gr * H2_lev(1:NH2_lev)%Energy)
    endif

! Vgrad should be larger than 1 km s-1 pc-1 = 1 / parsec = 1 / 3.0857d13
!SC March 2004: minimum gradient = thermal gradient over scale z=(distance+dmin)
!   dmin = 1.0D13

    ! initialize the level populations on the first call of CO_LVG
    IF(ICOUNT_CO.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL CO_LVG(Dens_H,Dens_H2,Dens_He,Dens_CO,Tn,op_H2,Vgrad,  &
                NCO_lev,YN_rovib_CO,ICOUNT_CO)
     v_variab(bv_CO_lev:ev_CO_lev) = dens_CO_lev(1:NCO_lev)
     Y(bv_CO_lev:ev_CO_lev) = DLOG(v_variab(bv_CO_lev:ev_CO_lev))
    ENDIF
    ! calculate the source terms for the population densities of CO
    ! and the collisional heating and cooling rates
    sum_CO = SUM(v_variab(bv_CO_lev:ev_CO_lev))
    ! Renormalize the CO level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_CO_lev:ev_CO_lev)=v_variab(bv_CO_lev:ev_CO_lev)*Dens_CO/sum_CO
    ! print*,v_variab(bv_CO_lev:ev_CO_lev)
    !--------------------------
    ! density of each CO level 
    !--------------------------
          dens_CO_lev(1:NCO_lev) = v_variab(bv_CO_lev:ev_CO_lev)

    ! initialize the level populations on the first call of SiO_LVG
    IF(ICOUNT_SiO.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL SiO_LVG(Dens_H,Dens_H2,Dens_He,Dens_SiO,Tn,op_H2,Vgrad,  &
                NSiO_lev,YN_rovib_SiO,ICOUNT_SiO)
     v_variab(bv_SiO_lev:ev_SiO_lev) = dens_SiO_lev(1:NSiO_lev)
     Y(bv_SiO_lev:ev_SiO_lev) = DLOG(v_variab(bv_SiO_lev:ev_SiO_lev))
    ENDIF
    ! calculate the source terms for the population densities of SiO
    ! and the collisional heating and cooling rates
    sum_SiO = SUM(v_variab(bv_SiO_lev:ev_SiO_lev))
    ! Renormalize the SiO level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_SiO_lev:ev_SiO_lev)=v_variab(bv_SiO_lev:ev_SiO_lev)*Dens_SiO/sum_SiO
    ! print*,v_variab(bv_SiO_lev:ev_SiO_lev)
    !--------------------------
    ! density of each SiO level 
    !--------------------------
          dens_SiO_lev(1:NSiO_lev) = v_variab(bv_SiO_lev:ev_SiO_lev)

    ! initialize the level populations on the first call of oH2O_LVG
    IF(ICOUNT_oH2O.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL oH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad,  &
                NoH2O_lev,YN_rovib_oH2O,op_H2O,ICOUNT_oH2O)
     v_variab(bv_oH2O_lev:ev_oH2O_lev) = dens_oH2O_lev(1:NoH2O_lev)
     Y(bv_oH2O_lev:ev_oH2O_lev) = DLOG(v_variab(bv_oH2O_lev:ev_oH2O_lev))
    ENDIF
    ! calculate the source terms for the population densities of o-H2O
    ! and the collisional heating and cooling rates
    sum_oH2O = SUM(v_variab(bv_oH2O_lev:ev_oH2O_lev))
    ! Renormalize the o-H2O level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_oH2O_lev:ev_oH2O_lev)=v_variab(bv_oH2O_lev:ev_oH2O_lev)*   &
                                      Dens_H2O*op_H2O/(1.+op_H2O)/sum_oH2O
    ! print*,v_variab(bv_oH2O_lev:ev_oH2O_lev)
    !--------------------------
    ! density of each oH2O level 
    !--------------------------
          dens_oH2O_lev(1:NoH2O_lev) = v_variab(bv_oH2O_lev:ev_oH2O_lev)

    ! initialize the level populations on the first call of pH2O_LVG
    IF(ICOUNT_pH2O.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL pH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad,  &
                NpH2O_lev,YN_rovib_pH2O,op_H2O,ICOUNT_pH2O)
     v_variab(bv_pH2O_lev:ev_pH2O_lev) = dens_pH2O_lev(1:NpH2O_lev)
     Y(bv_pH2O_lev:ev_pH2O_lev) = DLOG(v_variab(bv_pH2O_lev:ev_pH2O_lev))
    ENDIF
    ! calculate the source terms for the population densities of p-H2O
    ! and the collisional heating and cooling rates
    sum_pH2O = SUM(v_variab(bv_pH2O_lev:ev_pH2O_lev))
    ! Renormalize the p-H2O level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_pH2O_lev:ev_pH2O_lev)=v_variab(bv_pH2O_lev:ev_pH2O_lev)*   &
                                      Dens_H2O/(1.+op_H2O)/sum_pH2O
    ! print*,v_variab(bv_pH2O_lev:ev_pH2O_lev)
    !--------------------------
    ! density of each pH2O level 
    !--------------------------
          dens_pH2O_lev(1:NpH2O_lev) = v_variab(bv_pH2O_lev:ev_pH2O_lev)

   ! initialize the level populations on the first call of oNH3_LVG
    IF(ICOUNT_oNH3.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL oNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad,  &
                NoNH3_lev,YN_rovib_oNH3,op_NH3,ICOUNT_oNH3)
     v_variab(bv_oNH3_lev:ev_oNH3_lev) = dens_oNH3_lev(1:NoNH3_lev)
     Y(bv_oNH3_lev:ev_oNH3_lev) = DLOG(v_variab(bv_oNH3_lev:ev_oNH3_lev))
    ENDIF
    ! calculate the source terms for the population densities of o-NH3
    ! and the collisional heating and cooling rates
    sum_oNH3 = SUM(v_variab(bv_oNH3_lev:ev_oNH3_lev))
    ! Renormalize the o-NH3 level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_oNH3_lev:ev_oNH3_lev)=v_variab(bv_oNH3_lev:ev_oNH3_lev)*   &
                                      Dens_NH3*op_NH3/(1.+op_NH3)/sum_oNH3
    ! print*,v_variab(bv_oNH3_lev:ev_oNH3_lev)
    !--------------------------
    ! density of each oNH3 level 
    !--------------------------
          dens_oNH3_lev(1:NoNH3_lev) = v_variab(bv_oNH3_lev:ev_oNH3_lev)

    ! initialize the level populations on the first call of pNH3_LVG
    IF(ICOUNT_pNH3.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL pNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad,  &
                NpNH3_lev,YN_rovib_pNH3,op_NH3,ICOUNT_pNH3)
     v_variab(bv_pNH3_lev:ev_pNH3_lev) = dens_pNH3_lev(1:NpNH3_lev)
     Y(bv_pNH3_lev:ev_pNH3_lev) = DLOG(v_variab(bv_pNH3_lev:ev_pNH3_lev))
    ENDIF
    ! calculate the source terms for the population densities of p-NH3
    ! and the collisional heating and cooling rates
    sum_pNH3 = SUM(v_variab(bv_pNH3_lev:ev_pNH3_lev))
    ! Renormalize the p-NH3 level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_pNH3_lev:ev_pNH3_lev)=v_variab(bv_pNH3_lev:ev_pNH3_lev)*   &
                                      Dens_NH3/(1.+op_NH3)/sum_pNH3
    ! print*,v_variab(bv_pNH3_lev:ev_pNH3_lev)
    !--------------------------
    ! density of each pNH3 level 
    !--------------------------
          dens_pNH3_lev(1:NpNH3_lev) = v_variab(bv_pNH3_lev:ev_pNH3_lev)

    ! initialize the level populations on the first call of OH_LVG
    IF(ICOUNT_OH.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL OH_LVG(Dens_H,Dens_H2,Dens_He,Dens_OH,Tn,op_H2,Vgrad,  &
                NOH_lev,YN_rovib_OH,ICOUNT_OH)
     v_variab(bv_OH_lev:ev_OH_lev) = dens_OH_lev(1:NOH_lev)
     Y(bv_OH_lev:ev_OH_lev) = DLOG(v_variab(bv_OH_lev:ev_OH_lev))
    ENDIF
    ! calculate the source terms for the population densities of OH
    ! and the collisional heating and cooling rates
    sum_OH = SUM(v_variab(bv_OH_lev:ev_OH_lev))
    ! Renormalize the OH level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_OH_lev:ev_OH_lev)=v_variab(bv_OH_lev:ev_OH_lev)*   &
                                      Dens_OH/sum_OH
    ! print*,v_variab(bv_OH_lev:ev_OH_lev)
    !--------------------------
    ! density of each OH level 
    !--------------------------
          dens_OH_lev(1:NOH_lev) = v_variab(bv_OH_lev:ev_OH_lev)
          
    ! initialize the level populations on the first call of atype_LVG
    IF(ICOUNT_atype.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL atype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad,  &
                Natype_lev,YN_rovib_atype,ae_CH3OH,ICOUNT_atype)
     v_variab(bv_atype_lev:ev_atype_lev) = dens_atype_lev(1:Natype_lev)
     Y(bv_atype_lev:ev_atype_lev) = DLOG(v_variab(bv_atype_lev:ev_atype_lev))
    ENDIF
    ! calculate the source terms for the population densities of A-type
    ! and the collisional heating and cooling rates
    sum_atype = SUM(v_variab(bv_atype_lev:ev_atype_lev))
    ! Renormalize the CH3OH level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_atype_lev:ev_atype_lev)=v_variab(bv_atype_lev:ev_atype_lev)*   &
                                      Dens_CH3OH*ae_CH3OH/(1.+ae_CH3OH)/sum_atype
    ! print*,v_variab(bv_atype_lev:ev_atype_lev)
    !--------------------------
    ! density of each atype level 
    !--------------------------
          dens_atype_lev(1:Natype_lev) = v_variab(bv_atype_lev:ev_atype_lev)

    ! initialize the level populations on the first call of etype_LVG
    IF(ICOUNT_etype.EQ.0) THEN 
      Vgrad = Vmagson/(distance+dmin) 
    CALL etype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad,  &
                Netype_lev,YN_rovib_etype,ae_CH3OH,ICOUNT_etype)
     v_variab(bv_etype_lev:ev_etype_lev) = dens_etype_lev(1:Netype_lev)
     Y(bv_etype_lev:ev_etype_lev) = DLOG(v_variab(bv_etype_lev:ev_etype_lev))
    ENDIF
    ! calculate the source terms for the population densities of E-type
    ! and the collisional heating and cooling rates
    sum_etype = SUM(v_variab(bv_etype_lev:ev_etype_lev))
    ! Renormalize the CH3OH level population densities so that their sum
    ! is equal to the density of the corresponding chemical species
    v_variab(bv_etype_lev:ev_etype_lev)=v_variab(bv_etype_lev:ev_etype_lev)*   &
                                      Dens_CH3OH*ae_CH3OH/(1.+ae_CH3OH)/sum_etype
    ! print*,v_variab(bv_etype_lev:ev_etype_lev)
    !--------------------------
    ! density of each etype level 
    !--------------------------
          dens_etype_lev(1:Netype_lev) = v_variab(bv_etype_lev:ev_etype_lev)
                   
    !------------------------------------------------------------------
    ! In subroutine COMPUTE_OP are calculated the ortho:para H2 ratio,
    ! and the densities of ortho-H2 and para-H2 (cm-3).
    ! results : op_H2, Dens_orthoH2, Dens_paraH2
    ! (used in FINE_STRUCTURE_LINE, called in SOURCE)
    !------------------------------------------------------------------
    CALL COMPUTE_OP_H2

    !--------------------------------------------------------------------
    ! In subroutine CHEMISTRY the source terms of ([1], equations(3)-(14))
    ! are calculated as functions of the abundances of
    ! the species and the global properties of the fluids;
    ! results in :
    !    YN
    !    CupN,   CupI,   CupNEG
    !    CdownN, CdownI, CdownNEG
    !    YNn,    YNi,    YNneg
    !    Sn,     Si,     Sneg
    !    An,     Ai,     Aneg
    !    Bn,     Bi,     Bneg
    !--------------------------------------------------------------------
    CALL CHEMISTRY

    !--------------------------------------------------------------------
    ! In subroutine SOURCE, the source terms for the MHD equations are
    ! computed, depending on the results of subroutine CHEMISTRY and the
    ! function values : v_variab(i).
    ! SOURCE also computes H2 cooling and internal energy.
    ! Results in :
    !   An, Ai, Aneg
    !   Bn, Bi, Bneg
    !   H2_energy, YN_rovib_H2
    !--------------------------------------------------------------------
    CALL SOURCE

    !--- add chemical contribution YN (from CHEMISTRY) to the source of ---
    !--- internal energy of H2 H2_energy (from COMPUTE_H2 via SOURCE)   ---

    !--- Sel_tot_H2 is the rate of destruction of H2 by collisional dissociation
    Sel_tot_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_ch_H2))
    Sel_tne_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_ne_H2))
    Sel_tio_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_io_H2))
    Sel_tel_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_el_H2))

    H2_int_energy = SUM(DBLE(H2_lev(1:NH2_lev)%density * H2_lev(1:NH2_lev)%energy))
    H2_energy = H2_energy + kB * (YN(ind_H2) - Sel_tot_H2) / sum_H2 * &
         H2_int_energy &
       + kB * SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_ch_H2 * H2_lev(1:NH2_lev)%energy))

    !-----------------------------------------------------
    ! calculation of the derivatives of the MHD variables
    ! with respect to Z (distance into the cloud)
    !-----------------------------------------------------
    !--- Vn ---
!
! 29 mai 2001 - JLB : molecular cooling no longer included in Bn 
!                     but computed directly here. Old expression was:
!                     (with molec_cool in Bn, via cooling_n)
! one fluid model : Bne = Bn + Bi + Bneg

    SELECT CASE (Nfluids)
    CASE (1)
       ! one-fluid model : Bne = Bn + Bi + Bneg; Ane = 0; Sne = 0
       Bne = Bn + Bi + Bneg
       Ane = 0.0_DP
       Sne = 0.0_DP
       YNne = YNn + YNi + YNi
       DensityNe = DensityN + DensityI + DensityI
       RhoNe = RhoN + RhoI + RhoNEG
    CASE (2)
       ! two-fluid model : Bne = Bn
       Bne = Bn
       Ane = An
       Sne = Sn
       YNne = YNn
       DensityNe = DensityN
       RhoNe = RhoN
    CASE (3)
       ! three-fluid model : Bne = Bn
       Bne = Bn
       Ane = An
       Sne = Sn
       YNne = YNn
       DensityNe = DensityN
       RhoNe = RhoN
    END SELECT

!   DY(iv_Vn) = (GAMMA3*Sn*Vn2 - H2_energy - GAMMA2*An*Vn+Bn) &
!        / (RhoN*(GAMMA2*kT_muN-GAMMA1*Vn2))
!
!  molecular cooling is a direct function of DY(iv_Vn), so this is really an
!            implicit equation
!
!    CALL COMPUTE_MOLECULAR_COOLING
    aux1_dvn = GAMMA3*Sne*Vn2 - H2_energy - GAMMA2*Ane*Vn + Bne
    aux2_dvn = (DensityNe*GAMMA2*kB*Tn-RhoNe*GAMMA1*Vn2)

!  Add Magnetic term for J shocks

!   if ((shock_type == "J") .AND. (Nfluids == 1)) then
    if (Nfluids == 1) then
      aux2_dvn = aux2_dvn + GAMMA1 * Bfield2
    end if

    if (shock_type == "J") then

!SC      Vgrad = SQRT(grad_V*grad_V + 1.0D10 / (parsec*parsec)) * 1.0D-5
      Vgrad = SQRT(grad_V*grad_V + Vmagson*Vmagson/(distance+dmin)**2) * 1.0D-5 
      CALL COMPUTE_MOLECULAR_COOLING

    ! skip LVG treatment if LVG = 0
    if(LVG.ne.0) then
    ! calculate the source terms for the population densities of CO
    ! and the collisional heating and cooling rates
    CALL CO_LVG(Dens_H,Dens_H2,Dens_He,Dens_CO,Tn,op_H2,Vgrad*1.D5,  &
                NCO_lev,YN_rovib_CO,ICOUNT_CO)
      molec_cool = molec_cool + WWT_CO

    ! calculate the source terms for the population densities of SiO
    ! and the collisional heating and cooling rates
    CALL SiO_LVG(Dens_H,Dens_H2,Dens_He,Dens_SiO,Tn,op_H2,Vgrad*1.D5,  &
                NSiO_lev,YN_rovib_SiO,ICOUNT_SiO)
      molec_cool = molec_cool + WWT_SiO

    ! calculate the source terms for the population densities of o-H2O
    ! and the collisional heating and cooling rates
    CALL oH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NoH2O_lev,YN_rovib_oH2O,op_H2O,ICOUNT_oH2O)
      molec_cool = molec_cool + WWT_oH2O

    ! calculate the source terms for the population densities of p-H2O
    ! and the collisional heating and cooling rates
    CALL pH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NpH2O_lev,YN_rovib_pH2O,op_H2O,ICOUNT_pH2O)
      molec_cool = molec_cool + WWT_pH2O

    ! calculate the source terms for the population densities of o-NH3
    ! and the collisional heating and cooling rates
    CALL oNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NoNH3_lev,YN_rovib_oNH3,op_NH3,ICOUNT_oNH3)
      molec_cool = molec_cool + WWT_oNH3

    ! calculate the source terms for the population densities of p-NH3
    ! and the collisional heating and cooling rates
    CALL pNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NpNH3_lev,YN_rovib_pNH3,op_NH3,ICOUNT_pNH3)
      molec_cool = molec_cool + WWT_pNH3

    ! calculate the source terms for the population densities of OH
    ! and the collisional heating and cooling rates
    CALL OH_LVG(Dens_H,Dens_H2,Dens_He,Dens_OH,Tn,op_H2,Vgrad*1.D5,  &
                NOH_lev,YN_rovib_OH,ICOUNT_OH)
      molec_cool = molec_cool + WWT_OH

    ! skip methanol level populations if imeth = 0
    if(imeth.ne.0) then
    ! calculate the source terms for the population densities of A-type
    ! and the collisional heating and cooling rates
    CALL atype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Natype_lev,YN_rovib_atype,ae_CH3OH,ICOUNT_atype)
      molec_cool = molec_cool + WWT_atype

    ! calculate the source terms for the population densities of E-type
    ! and the collisional heating and cooling rates
    CALL etype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Netype_lev,YN_rovib_etype,ae_CH3OH,ICOUNT_etype)
      molec_cool = molec_cool + WWT_etype
    endif
    endif


      dvn2 = (aux1_dvn - molec_cool) / aux2_dvn
      DY(iv_Vn) = dvn2

    else

!      dvn1 = (aux1_dvn - molec_cool) / aux2_dvn
!      Vgrad = SQRT(dvn1*dvn1 + 1.0D10 / (parsec*parsec)) * 1.0D-5
! SC May 06
      dvn1 = v_l_der(iv_Vn)  
      Vgrad = SQRT(dvn1*dvn1 + Vmagson*Vmagson/(distance+dmin)**2) * 1.0D-5

      CALL COMPUTE_MOLECULAR_COOLING

    ! skip LVG treatment if LVG = 0
    if(LVG.ne.0) then
    ! calculate the source terms for the population densities of CO
    ! and the collisional heating and cooling rates
    CALL CO_LVG(Dens_H,Dens_H2,Dens_He,Dens_CO,Tn,op_H2,Vgrad*1.D5,  &
                NCO_lev,YN_rovib_CO,ICOUNT_CO)
      molec_cool = molec_cool + WWT_CO

    ! calculate the source terms for the population densities of SiO
    ! and the collisional heating and cooling rates
    CALL SiO_LVG(Dens_H,Dens_H2,Dens_He,Dens_SiO,Tn,op_H2,Vgrad*1.D5,  &
                NSiO_lev,YN_rovib_SiO,ICOUNT_SiO)
      molec_cool = molec_cool + WWT_SiO

    ! calculate the source terms for the population densities of o-H2O
    ! and the collisional heating and cooling rates
    CALL oH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NoH2O_lev,YN_rovib_oH2O,op_H2O,ICOUNT_oH2O)
      molec_cool = molec_cool + WWT_oH2O

    ! calculate the source terms for the population densities of p-H2O
    ! and the collisional heating and cooling rates
    CALL pH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NpH2O_lev,YN_rovib_pH2O,op_H2O,ICOUNT_pH2O)
      molec_cool = molec_cool + WWT_pH2O

    ! calculate the source terms for the population densities of o-NH3
    ! and the collisional heating and cooling rates
    CALL oNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NoNH3_lev,YN_rovib_oNH3,op_NH3,ICOUNT_oNH3)
      molec_cool = molec_cool + WWT_oNH3

    ! calculate the source terms for the population densities of p-NH3
    ! and the collisional heating and cooling rates
    CALL pNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NpNH3_lev,YN_rovib_pNH3,op_NH3,ICOUNT_pNH3)
      molec_cool = molec_cool + WWT_pNH3

    ! calculate the source terms for the population densities of OH
    ! and the collisional heating and cooling rates
    CALL OH_LVG(Dens_H,Dens_H2,Dens_He,Dens_OH,Tn,op_H2,Vgrad*1.D5,  &
                NOH_lev,YN_rovib_OH,ICOUNT_OH)
      molec_cool = molec_cool + WWT_OH

    ! skip methanol level populations if imeth = 0
    if(imeth.ne.0) then
    ! calculate the source terms for the population densities of A-type
    ! and the collisional heating and cooling rates
    CALL atype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Natype_lev,YN_rovib_atype,ae_CH3OH,ICOUNT_atype)
      molec_cool = molec_cool + WWT_atype

    ! calculate the source terms for the population densities of E-type
    ! and the collisional heating and cooling rates
    CALL etype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Netype_lev,YN_rovib_etype,ae_CH3OH,ICOUNT_etype)
      molec_cool = molec_cool + WWT_etype
    endif
    endif

      dvn2 = (aux1_dvn - molec_cool) / aux2_dvn
      gvn1 = dvn1 - dvn2

      Vgrad = SQRT(dvn2*dvn2 + Vmagson*Vmagson/(distance+dmin)**2) * 1.0D-5
      CALL COMPUTE_MOLECULAR_COOLING

    ! skip LVG treatment if LVG = 0
    if(LVG.ne.0) then
    ! calculate the source terms for the population densities of CO
    ! and the collisional heating and cooling rates
    CALL CO_LVG(Dens_H,Dens_H2,Dens_He,Dens_CO,Tn,op_H2,Vgrad*1.D5,  &
                NCO_lev,YN_rovib_CO,ICOUNT_CO)
      molec_cool = molec_cool + WWT_CO

    ! calculate the source terms for the population densities of SiO
    ! and the collisional heating and cooling rates
    CALL SiO_LVG(Dens_H,Dens_H2,Dens_He,Dens_SiO,Tn,op_H2,Vgrad*1.D5,  &
                NSiO_lev,YN_rovib_SiO,ICOUNT_SiO)
      molec_cool = molec_cool + WWT_SiO

    ! calculate the source terms for the population densities of o-H2O
    ! and the collisional heating and cooling rates
    CALL oH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NoH2O_lev,YN_rovib_oH2O,op_H2O,ICOUNT_oH2O)
      molec_cool = molec_cool + WWT_oH2O

    ! calculate the source terms for the population densities of p-H2O
    ! and the collisional heating and cooling rates
    CALL pH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NpH2O_lev,YN_rovib_pH2O,op_H2O,ICOUNT_pH2O)
      molec_cool = molec_cool + WWT_pH2O

    ! calculate the source terms for the population densities of o-NH3
    ! and the collisional heating and cooling rates
    CALL oNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NoNH3_lev,YN_rovib_oNH3,op_NH3,ICOUNT_oNH3)
      molec_cool = molec_cool + WWT_oNH3

    ! calculate the source terms for the population densities of p-NH3
    ! and the collisional heating and cooling rates
    CALL pNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NpNH3_lev,YN_rovib_pNH3,op_NH3,ICOUNT_pNH3)
      molec_cool = molec_cool + WWT_pNH3

    ! calculate the source terms for the population densities of OH
    ! and the collisional heating and cooling rates
    CALL OH_LVG(Dens_H,Dens_H2,Dens_He,Dens_OH,Tn,op_H2,Vgrad*1.D5,  &
                NOH_lev,YN_rovib_OH,ICOUNT_OH)
      molec_cool = molec_cool + WWT_OH

    ! skip methanol level populations if imeth = 0
    if(imeth.ne.0) then
    ! calculate the source terms for the population densities of A-type
    ! and the collisional heating and cooling rates
    CALL atype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Natype_lev,YN_rovib_atype,ae_CH3OH,ICOUNT_atype)
      molec_cool = molec_cool + WWT_atype

    ! calculate the source terms for the population densities of E-type
    ! and the collisional heating and cooling rates
    CALL etype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Netype_lev,YN_rovib_etype,ae_CH3OH,ICOUNT_etype)
      molec_cool = molec_cool + WWT_etype
    endif
    endif

      gvn2 = dvn2 - (aux1_dvn - molec_cool) / aux2_dvn

! SC: 2 problems with loop: gvn1 = 0 at end of first iteration, 
! & update of dvn2 incorrect: should be gvn2*(dvn1-dvn2) instead of gvn2*gvn1
!      ii = 0
!      do while ( abs(gvn1) > 1.0e-20 .AND. ii < 1000)
!        Vgrad = SQRT(dvn1*dvn1 + 1.0D10 / (parsec*parsec)) * 1.0D-5
!        CALL COMPUTE_MOLECULAR_COOLING
!        gvn2 = dvn2 - (aux1_dvn - molec_cool) / aux2_dvn
!        dvn1 = dvn2
!        dvn2 = dvn2 - gvn2 * gvn1 / (gvn1 - gvn2)
!        gvn1 = gvn2
!        ii = ii + 1
!      end do

! SC: proposed loop for Newton's method of solving for dVn/dz
      ii = 0
      do while ( abs(gvn1-gvn2) > 1.0e-6*(abs(gvn1)).AND. ii < 1000)
        corrdv = dvn2 - gvn2 * (dvn1 - dvn2) / (gvn1 - gvn2)
        dvn1 = dvn2
        gvn1 = gvn2
        dvn2 = corrdv
        Vgrad = SQRT(dvn2*dvn2 + Vmagson*Vmagson/(distance+dmin)**2) * 1.0D-5
        CALL COMPUTE_MOLECULAR_COOLING

    ! skip LVG treatment if LVG = 0
    if(LVG.ne.0) then
    ! calculate the source terms for the population densities of CO
    ! and the collisional heating and cooling rates
    CALL CO_LVG(Dens_H,Dens_H2,Dens_He,Dens_CO,Tn,op_H2,Vgrad*1.D5,  &
                NCO_lev,YN_rovib_CO,ICOUNT_CO)
      molec_cool = molec_cool + WWT_CO

    ! calculate the source terms for the population densities of SiO
    ! and the collisional heating and cooling rates
    CALL SiO_LVG(Dens_H,Dens_H2,Dens_He,Dens_SiO,Tn,op_H2,Vgrad*1.D5,  &
                NSiO_lev,YN_rovib_SiO,ICOUNT_SiO)
      molec_cool = molec_cool + WWT_SiO

    ! calculate the source terms for the population densities of o-H2O
    ! and the collisional heating and cooling rates
    CALL oH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NoH2O_lev,YN_rovib_oH2O,op_H2O,ICOUNT_oH2O)
      molec_cool = molec_cool + WWT_oH2O

    ! calculate the source terms for the population densities of p-H2O
    ! and the collisional heating and cooling rates
    CALL pH2O_LVG(Dens_H,Dens_H2,Dens_He,Dens_H2O,Tn,op_H2,Vgrad*1.D5,  &
                NpH2O_lev,YN_rovib_pH2O,op_H2O,ICOUNT_pH2O)
      molec_cool = molec_cool + WWT_pH2O

    ! calculate the source terms for the population densities of o-NH3
    ! and the collisional heating and cooling rates
    CALL oNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NoNH3_lev,YN_rovib_oNH3,op_NH3,ICOUNT_oNH3)
      molec_cool = molec_cool + WWT_oNH3

    ! calculate the source terms for the population densities of p-NH3
    ! and the collisional heating and cooling rates
    CALL pNH3_LVG(Dens_H,Dens_H2,Dens_He,Dens_NH3,Tn,op_H2,Vgrad*1.D5,  &
                NpNH3_lev,YN_rovib_pNH3,op_NH3,ICOUNT_pNH3)
      molec_cool = molec_cool + WWT_pNH3

    ! calculate the source terms for the population densities of OH
    ! and the collisional heating and cooling rates
    CALL OH_LVG(Dens_H,Dens_H2,Dens_He,Dens_OH,Tn,op_H2,Vgrad*1.D5,  &
                NOH_lev,YN_rovib_OH,ICOUNT_OH)
      molec_cool = molec_cool + WWT_OH

    ! skip methanol level populations if imeth = 0
    if(imeth.ne.0) then
    ! calculate the source terms for the population densities of A-type
    ! and the collisional heating and cooling rates
    CALL atype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Natype_lev,YN_rovib_atype,ae_CH3OH,ICOUNT_atype)
      molec_cool = molec_cool + WWT_atype

    ! calculate the source terms for the population densities of E-type
    ! and the collisional heating and cooling rates
    CALL etype_LVG(Dens_H,Dens_H2,Dens_He,Dens_CH3OH,Tn,op_H2,Vgrad*1.D5,  &
                Netype_lev,YN_rovib_etype,ae_CH3OH,ICOUNT_etype)
      molec_cool = molec_cool + WWT_etype
    endif
    endif

        gvn2 = dvn2 - (aux1_dvn - molec_cool) / aux2_dvn
        ii = ii + 1
      end do

      DY(iv_Vn) = dvn2
!SC     Vgrad = SQRT(dvn1*dvn1 + 1.0D10 / (parsec*parsec)) * 1.0D-5
      Vgrad = SQRT(dvn2*dvn2 + Vmagson*Vmagson/(distance+dmin)**2) * 1.0D-5
    endif

!  Save DY(iv_Vn) for comparison with second order term in MAIN

    save_dv = dvn2
        if ((shock_type == "J") .AND. (viscosity)) then
      DY(iv_Vn) = - grad_V
    end if

    !--- Vi, according to the number of fluids ---
    SELECT CASE (Nfluids)
    CASE (1)
       ! one-fluid model : Vi=Vn
       DY(iv_Vi) = DY(iv_Vn)
    CASE DEFAULT
       ! two- or three- fluid model : ions and electrons have the same velocity
       DY(iv_Vi) = (-GAMMA3*Sn*Vi2 + GAMMA2*An*Vi + Bi + Bneg) &
            / (GAMMA2*kB*DensityI*(Ti+Te) + GAMMA1*Bfield2 - GAMMA1*Vi2*(RhoI+RhoNEG+Rho_GRAIN_charge))
    END SELECT

    !--- RhoN ---
    DY(iv_RhoN) = (Sn-RhoN*DY(iv_Vn)) / Vn
    !--- RhoI ---
    DY(iv_RhoI) = (Si-RhoI*DY(iv_Vi)) / Vi
    !--- RHONEG ---
    DY(iv_RhoNEG) = (Sneg-RhoNEG*DY(iv_Vi)) / Vi
    
    !--- Tn ---
    
!    molec_cool = molec_cool - SUM(dens_CO_lev(1:NCO_lev) * HANDC(1:NCO_lev)) * kB

    if ((shock_type == "J") .AND. (viscosity)) then
      DY(iv_Tn) = Bne - molec_cool - Ane * Vn + 0.5_dp * Sne * Vn2 &
                - GAMMA1 * kB * Tn * YNne - H2_energy &
                - DensityNe*kB*Tn * DY(iv_Vn) &
                - RhoNe * XLL2 * DY(iv_Vn) * DY(iv_Vn) * DY(iv_Vn)
      DY(iv_Tn) = DY(iv_Tn) / (GAMMA1 * kB * Vn * DensityNe)
    else
      DY(iv_Tn) = Bne - molec_cool - Ane * Vn + 0.5_dp * Sne * Vn2 &
                - GAMMA1 * kB * Tn * YNne - H2_energy &
                - DensityNe*kB*Tn * DY(iv_Vn)
      DY(iv_Tn) = DY(iv_Tn) / (GAMMA1 * kB * Vn * DensityNe)
    endif
     !DY(iv_Tn) = 0.0_DP

    !--- Ti and Te, according to the number of fluids ---
    SELECT CASE (Nfluids)
    CASE (1)
       ! one fluid model : Tn=Ti=Te
       DY(iv_Ti) = DY(iv_Tn)
       DY(iv_Te) = DY(iv_Tn)
    CASE (2)
       ! two fluid model : Ti=Te
!      DY(iv_Ti)=-0.25_DP*Sne*Vi2-GAMMA1*YNi*kB*Ti+0.5_DP*Ane*Vi
       DY(iv_Ti)=-0.25_DP*Sn*Vi2-GAMMA1*YNi*kB*Ti+0.5_DP*An*Vi
       DY(iv_Ti)=DY(iv_Ti)+0.5_DP*(Bi+Bneg)-DensityI*kB*Ti*DY(iv_Vi)
       DY(iv_Ti)=DY(iv_Ti)/(GAMMA1*kB*V_DensityI)
       DY(iv_Te)=DY(iv_Ti)
    CASE (3)
       ! three fluid model
       DY(iv_Ti) = (0.5_DP*Si*Vi2 - GAMMA1*YNi*kB*Ti - Ai*Vi &
            +Bi - DensityI*kB*Ti*DY(iv_Vi)) / (GAMMA1*kB*V_DensityI)
       DY(iv_Te) = (0.5_DP*Sneg*Vi2 - GAMMA1*YNi*kB*Te -Aneg*Vi &
            +Bneg - DensityI*kB*Te*DY(iv_Vi)) /(GAMMA1*kB*V_DensityI)
    END SELECT

    !--- DensityN ---
!   DY(iv_DensityN)=(YNne-DensityNe*DY(iv_Vn))/Vn
    DY(iv_DensityN)=(YNn-DensityN*DY(iv_Vn))/Vn
    !--- DensityI ---
    DY(iv_DensityI)=(YNi-DensityI*DY(iv_Vi))/Vi

    DY(iv_gv) = 0.0_dp
    !--- Velocity gradient for J shocks ---
    if (shock_type == "J") then
     if (viscosity) then
      if (Nfluids == 1) then
      DY(iv_gv) = Bne - molec_cool - 0.5_dp * Sne * Vn2 - XLL2 * DY(iv_Vn) * DY(iv_Vn) * Sne &
              - GAMMA2 * kB * Tn * YNne - H2_energy - RhoNe * Vn2 * DY(iv_Vn) &
              - GAMMA2 * DensityNe * Vn * kB * DY(iv_Tn) + Bfield2 * DY(iv_Vn)
      DY(iv_gv) = DY(iv_gv) / (2.0_dp * RhoNe * XLL2 * DY(iv_Vn) * Vn)
       else
      DY(iv_gv) = Bne - molec_cool - 0.5_dp * Sne * Vn2 - XLL2 * DY(iv_Vn) * DY(iv_Vn) * Sne &
              - GAMMA2 * kB * Tn * YNne - H2_energy - RhoNe * Vn2 * DY(iv_Vn) &
              - GAMMA2 * DensityNe * Vn * kB * DY(iv_Tn)
      DY(iv_gv) = DY(iv_gv) / (2.0_dp * RhoNe * XLL2 * DY(iv_Vn) * Vn)
      end if
! iv_gv is used for -dVn/dz
      DY(iv_gv) = - DY(iv_gv)
     end if

    end if

    !----------------------------------------------------------------
    ! calculates the derivatives of the chemical species' densities
    ! with respect to Z (distance into the cloud)
    !----------------------------------------------------------------
    !--- neutrals (velocity = Vn) ---
    DY(bv_neu:ev_neu)=(YN(b_neu:e_neu) - &
         v_variab(bv_neu:ev_neu)*DY(iv_Vn) )/Vn
    !--- species on grains (velocity = Vi) ---
    DY(bv_gra:ev_gra)=(YN(b_gra:e_gra) - &
         v_variab(bv_gra:ev_gra)*DY(iv_Vi) )/Vi
    !--- species in grain cores (velocity = Vi) ---
    DY(bv_cor:ev_cor)=(YN(b_cor:e_cor) - &
         v_variab(bv_cor:ev_cor)*DY(iv_Vi) )/Vi
    !--- ions > 0 (velocity = Vi) ---
    DY(bv_ion:ev_ion)=(YN(b_ion:e_ion) - &
         v_variab(bv_ion:ev_ion)*DY(iv_Vi) )/Vi
    !--- ions < 0 (velocity = Vi) ---
    DY(bv_neg:ev_neg)=(YN(b_neg:e_neg) - &
         v_variab(bv_neg:ev_neg)*DY(iv_Vi) )/Vi

    !------------------------------------------
    ! densities of ro-vibrational levels of H2
    !------------------------------------------
    ! The contribution of chemical reactions YN(ind_H2) is included in proportion
    ! to the fractional population of the level : v_variab(J)/sum_H2
    ! with the exception of contributions from selective reactions 

    DY(bv_H2_lev:ev_H2_lev) = (YN_rovib_H2(1:NH2_lev) &
         + (YN(ind_H2) - Sel_tot_H2 - For_gr_H2) * v_variab(bv_H2_lev:ev_H2_lev) / sum_H2 &
         + H2_lev(1:NH2_lev)%density * Sel_ch_H2 &
         + For_gr_H2 * H2_lev(1:NH2_lev)%Form_gr &
         - v_variab(bv_H2_lev:ev_H2_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_H2_lev:ev_H2_lev) = 0._DP

    !------------------------------------------
    ! densities of rotational levels of CO
    !------------------------------------------
    !
    ! The contribution of chemical reactions YN(ind_CO) is included in proportion
    ! to the fractional population of the level : v_variab(J)/sum_CO
    sum_CO = SUM(v_variab(bv_CO_lev:ev_CO_lev))

    DY(bv_CO_lev:ev_CO_lev) = (YN_rovib_CO(1:NCO_lev) &
         + YN(ind_CO) * v_variab(bv_CO_lev:ev_CO_lev) / sum_CO &
         - v_variab(bv_CO_lev:ev_CO_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_CO_lev:ev_CO_lev) = 0._DP
    ! skip CO level populations if LVG = 0
      if(LVG.eq.0) DY(bv_CO_lev:ev_CO_lev) = 0._DP
         
!   DY(bv_CO_lev:ev_CO_lev) = (0.d0 &
!        + YN(ind_CO) * v_variab(bv_CO_lev:ev_CO_lev) / sum_CO &
!        - v_variab(bv_CO_lev:ev_CO_lev)*DY(iv_Vn) ) / Vn

    !print*, ICOUNT_CO, bv_CO_lev, ev_CO_lev, Y(bv_CO_lev:ev_CO_lev)
    !------------------------------------------
    ! densities of rotational levels of SiO
    !------------------------------------------
    !
    ! The contribution of chemical reactions YN(ind_SiO) is included in proportion
    ! to the fractional population of the level : v_variab(J)/sum_SiO
    sum_SiO = SUM(v_variab(bv_SiO_lev:ev_SiO_lev))

    DY(bv_SiO_lev:ev_SiO_lev) = (YN_rovib_SiO(1:NSiO_lev) &
         + YN(ind_SiO) * v_variab(bv_SiO_lev:ev_SiO_lev) / sum_SiO &
         - v_variab(bv_SiO_lev:ev_SiO_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_SiO_lev:ev_SiO_lev) = 0._DP
    ! skip SiO level populations if LVG = 0
      if(LVG.eq.0) DY(bv_SiO_lev:ev_SiO_lev) = 0._DP
         
!   DY(bv_SiO_lev:ev_SiO_lev) = (0.d0 &
!        + YN(ind_SiO) * v_variab(bv_SiO_lev:ev_SiO_lev) / sum_SiO &
!        - v_variab(bv_SiO_lev:ev_SiO_lev)*DY(iv_Vn) ) / Vn

    !------------------------------------------
    ! densities of rotational levels of o-H2O
    !------------------------------------------
    !
    ! The contribution of chemical reactions YN(ind_H2O) is included in proportion
    ! to the fractional population of the level : v_variab(J)/(sum_oH2O+sum_pH2O)
    sum_oH2O = SUM(v_variab(bv_oH2O_lev:ev_oH2O_lev))
    sum_pH2O = SUM(v_variab(bv_pH2O_lev:ev_pH2O_lev))

    DY(bv_oH2O_lev:ev_oH2O_lev) = (YN_rovib_oH2O(1:NoH2O_lev) &
         + YN(ind_H2O) * v_variab(bv_oH2O_lev:ev_oH2O_lev) / (sum_oH2O+sum_pH2O) &
         - v_variab(bv_oH2O_lev:ev_oH2O_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_oH2O_lev:ev_oH2O_lev) = 0._DP
    ! skip o-H2O level populations if LVG = 0
      if(LVG.eq.0) DY(bv_oH2O_lev:ev_oH2O_lev) = 0._DP
         
    !print*, ICOUNT_oH2O, bv_oH2O_lev, ev_oH2O_lev, Y(bv_oH2O_lev:ev_oH2O_lev)
    !------------------------------------------
    ! densities of rotational levels of p-H2O
    !------------------------------------------
    !
    ! The contribution of chemical reactions YN(ind_H2O) is included in proportion
    ! to the fractional population of the level : v_variab(J)/(sum_oH2O+sum_pH2O)

    DY(bv_pH2O_lev:ev_pH2O_lev) = (YN_rovib_pH2O(1:NpH2O_lev) &
         + YN(ind_H2O) * v_variab(bv_pH2O_lev:ev_pH2O_lev) / (sum_oH2O+sum_pH2O) &
         - v_variab(bv_pH2O_lev:ev_pH2O_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_pH2O_lev:ev_pH2O_lev) = 0._DP
    ! skip p-H2O level populations if LVG = 0
      if(LVG.eq.0) DY(bv_pH2O_lev:ev_pH2O_lev) = 0._DP
         
    !print*, ICOUNT_pH2O, bv_pH2O_lev, ev_pH2O_lev, Y(bv_pH2O_lev:ev_pH2O_lev)
    !------------------------------------------
    ! densities of rotational levels of o-NH3
    !------------------------------------------
    !
    ! The contribution of chemical reactions YN(ind_NH3) is included in proportion
    ! to the fractional population of the level : v_variab(J)/(sum_oNH3+sum_pNH3)
    sum_oNH3 = SUM(v_variab(bv_oNH3_lev:ev_oNH3_lev))
    sum_pNH3 = SUM(v_variab(bv_pNH3_lev:ev_pNH3_lev))

    DY(bv_oNH3_lev:ev_oNH3_lev) = (YN_rovib_oNH3(1:NoNH3_lev) &
         + YN(ind_NH3) * v_variab(bv_oNH3_lev:ev_oNH3_lev) / (sum_oNH3+sum_pNH3) &
         - v_variab(bv_oNH3_lev:ev_oNH3_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_oNH3_lev:ev_oNH3_lev) = 0._DP
    ! skip o-NH3 level populations if LVG = 0
      if(LVG.eq.0) DY(bv_oNH3_lev:ev_oNH3_lev) = 0._DP
         
    !print*, ICOUNT_oNH3, bv_oNH3_lev, ev_oNH3_lev, Y(bv_oNH3_lev:ev_oNH3_lev)
    !------------------------------------------
    ! densities of rotational levels of p-NH3
    !------------------------------------------
    !
    ! The contribution of chemical reactions YN(ind_NH3) is included in proportion
    ! to the fractional population of the level : v_variab(J)/(sum_oNH3+sum_pNH3)

    DY(bv_pNH3_lev:ev_pNH3_lev) = (YN_rovib_pNH3(1:NpNH3_lev) &
         + YN(ind_NH3) * v_variab(bv_pNH3_lev:ev_pNH3_lev) / (sum_oNH3+sum_pNH3) &
         - v_variab(bv_pNH3_lev:ev_pNH3_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_pNH3_lev:ev_pNH3_lev) = 0._DP
    ! skip p-NH3 level populations if LVG = 0
      if(LVG.eq.0) DY(bv_pNH3_lev:ev_pNH3_lev) = 0._DP
         
    !print*, ICOUNT_pNH3, bv_pNH3_lev, ev_pNH3_lev, Y(bv_pNH3_lev:ev_pNH3_lev)
    !------------------------------------------
    ! densities of rotational levels of OH
    !------------------------------------------
    !
    sum_OH = SUM(v_variab(bv_OH_lev:ev_OH_lev))

    DY(bv_OH_lev:ev_OH_lev) = (YN_rovib_OH(1:NOH_lev) &
         + YN(ind_OH) * v_variab(bv_OH_lev:ev_OH_lev) / sum_OH &
         - v_variab(bv_OH_lev:ev_OH_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_OH_lev:ev_OH_lev) = 0._DP
    ! skip OH level populations if LVG = 0
      if(LVG.eq.0) DY(bv_OH_lev:ev_OH_lev) = 0._DP
         
    !print*, ICOUNT_OH, bv_OH_lev, ev_OH_lev, Y(bv_OH_lev:ev_OH_lev)
    !print*, Z, Y

    !------------------------------------------
    ! densities of rotational levels of A-type
    !------------------------------------------
     
    ! The contribution of chemical reactions YN(ind_CH3OH) is included in proportion
    ! to the fractional population of the level : v_variab(J)/(sum_atype+sum_etype)
    sum_atype = SUM(v_variab(bv_atype_lev:ev_atype_lev))
    sum_etype = SUM(v_variab(bv_etype_lev:ev_etype_lev))

    DY(bv_atype_lev:ev_atype_lev) = (YN_rovib_atype(1:Natype_lev) &
         + YN(ind_CH3OH) * v_variab(bv_atype_lev:ev_atype_lev) / (sum_atype+sum_etype) &
         - v_variab(bv_atype_lev:ev_atype_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_atype_lev:ev_atype_lev) = 0._DP
    ! skip methanol level populations if imeth = 0
      if((imeth.eq.0).or.(LVG.eq.0)) DY(bv_atype_lev:ev_atype_lev) = 0._DP
         
    !print*, ICOUNT_atype, bv_atype_lev, ev_atype_lev, Y(bv_atype_lev:ev_atype_lev)

    !------------------------------------------
    ! densities of rotational levels of E-type
    !------------------------------------------
    !
    ! The contribution of chemical reactions YN(ind_CH3OH) is included in proportion
    ! to the fractional population of the level : v_variab(J)/(sum_atype+sum_etype)

    DY(bv_etype_lev:ev_etype_lev) = (YN_rovib_etype(1:Netype_lev) &
         + YN(ind_CH3OH) * v_variab(bv_etype_lev:ev_etype_lev) / (sum_atype+sum_etype) &
         - v_variab(bv_etype_lev:ev_etype_lev)*DY(iv_Vn) ) / Vn
      if (shock_type == "S") DY(bv_etype_lev:ev_etype_lev) = 0._DP
    ! skip methanol level populations if imeth = 0
      if((imeth.eq.0).or.(LVG.eq.0)) DY(bv_etype_lev:ev_etype_lev) = 0._DP
         
    !print*, ICOUNT_etype, bv_etype_lev, ev_etype_lev, Y(bv_etype_lev:ev_etype_lev)
    !print*, Z, Y

    !--------------------------------------------------------------
    ! save values of v_variab and its derivative DY
    ! to the arrays v_l_var and v_l_der.
    !--------------------------------------------------------------
    v_l_var(1:N)=v_variab(1:N)
    v_l_der(1:N)=DY(1:N)

    !---------------------------------
    ! divide DY by Y = v_variab
    !---------------------------------
    ! MHD variables : tiny avoids overflow if the denominator is very small
    DY(1:Nv_MHD) = &
         DY(1:Nv_MHD) / &
         v_variab(1:Nv_MHD)
    dyntime = DABS(1./(DY(iv_Tn) + tiny))

    ! species : tiny avoids overflow if the denominator is very small
    DY(bv_speci:ev_speci) = &
         DY(bv_speci:ev_speci) / &
         v_variab(bv_speci:ev_speci)
         
    ! H2 levels : tiny avoids overflow if the denominator is very small
    DY(bv_H2_lev:ev_H2_lev) = &
         DY(bv_H2_lev:ev_H2_lev) / &
         (v_variab(bv_H2_lev:ev_H2_lev) + tiny)

    ! CO levels : tiny avoids overflow if the denominator is very small
    DY(bv_CO_lev:ev_CO_lev) = &
         DY(bv_CO_lev:ev_CO_lev) / &
         (v_variab(bv_CO_lev:ev_CO_lev) + tiny)
    timescale_CO(1:NCO_lev) = DABS(1./(DY(bv_CO_lev:ev_CO_lev) + tiny))

    ! SiO levels : tiny avoids overflow if the denominator is very small
    DY(bv_SiO_lev:ev_SiO_lev) = &
         DY(bv_SiO_lev:ev_SiO_lev) / &
         (v_variab(bv_SiO_lev:ev_SiO_lev) + tiny)
    timescale_SiO(1:NSiO_lev) = DABS(1./(DY(bv_SiO_lev:ev_SiO_lev) + tiny))

    ! o-H2O levels : tiny avoids overflow if the denominator is very small
    DY(bv_oH2O_lev:ev_oH2O_lev) = &
         DY(bv_oH2O_lev:ev_oH2O_lev) / &
         (v_variab(bv_oH2O_lev:ev_oH2O_lev) + tiny)
    timescale_oH2O(1:NoH2O_lev) = DABS(1./(DY(bv_oH2O_lev:ev_oH2O_lev) + tiny))

    ! p-H2O levels : tiny avoids overflow if the denominator is very small
    DY(bv_pH2O_lev:ev_pH2O_lev) = &
         DY(bv_pH2O_lev:ev_pH2O_lev) / &
         (v_variab(bv_pH2O_lev:ev_pH2O_lev) + tiny)
    timescale_pH2O(1:NpH2O_lev) = DABS(1./(DY(bv_pH2O_lev:ev_pH2O_lev) + tiny))

    ! o-NH3 levels : tiny avoids overflow if the denominator is very small
    DY(bv_oNH3_lev:ev_oNH3_lev) = &
         DY(bv_oNH3_lev:ev_oNH3_lev) / &
         (v_variab(bv_oNH3_lev:ev_oNH3_lev) + tiny)
    timescale_oNH3(1:NoNH3_lev) = DABS(1./(DY(bv_oNH3_lev:ev_oNH3_lev) + tiny))

    ! p-NH3 levels : tiny avoids overflow if the denominator is very small
    DY(bv_pNH3_lev:ev_pNH3_lev) = &
         DY(bv_pNH3_lev:ev_pNH3_lev) / &
         (v_variab(bv_pNH3_lev:ev_pNH3_lev) + tiny)
    timescale_pNH3(1:NpNH3_lev) = DABS(1./(DY(bv_pNH3_lev:ev_pNH3_lev) + tiny))

    ! OH levels : tiny avoids overflow if the denominator is very small
    DY(bv_OH_lev:ev_OH_lev) = &
         DY(bv_OH_lev:ev_OH_lev) / &
         (v_variab(bv_OH_lev:ev_OH_lev) + tiny)
    timescale_OH(1:NOH_lev) = DABS(1./(DY(bv_OH_lev:ev_OH_lev) + tiny))

    ! A-type levels : tiny avoids overflow if the denominator is very small
    DY(bv_atype_lev:ev_atype_lev) = &
         DY(bv_atype_lev:ev_atype_lev) / &
         (v_variab(bv_atype_lev:ev_atype_lev) + tiny)
    timescale_atype(1:Natype_lev) = DABS(1./(DY(bv_atype_lev:ev_atype_lev) + tiny))

    ! E-type levels : tiny avoids overflow if the denominator is very small
    DY(bv_etype_lev:ev_etype_lev) = &
         DY(bv_etype_lev:ev_etype_lev) / &
         (v_variab(bv_etype_lev:ev_etype_lev) + tiny)
    timescale_etype(1:Netype_lev) = DABS(1./(DY(bv_etype_lev:ev_etype_lev) + tiny))

    ! Include molec_cool in Bn as before, so that Bn may be used without
    ! modifications in other parts of the code

    Bn = Bn - molec_cool

  END SUBROUTINE DIFFUN




END MODULE MODULE_EVOLUTION
