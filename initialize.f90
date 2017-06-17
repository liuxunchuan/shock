
MODULE MODULE_INITIALIZE
  !*****************************************************************************
  !** The module 'MODULE_INITIALIZE' calls subroutines for reading            **
  !** and initializing :                                                      **
  !**      * the shock parameters                                             **
  !**      * the chemical species                                             **
  !**      * the chemical reactions                                           **
  !**      * the H2 molecule (levels, collision rates, lines)                 **
  !**      * the Fe+ ion (emissivities and integrated intensities)            **
  !*****************************************************************************
  ! SC May 06: 
  ! - rg: Dens_GRAIN_init now calculated correctly from mass in cores only
  ! - ng: check that Dens(G,G+,G-) is consistent with Dens_GRAIN, else stop
  ! -     Automatic initialization of Dens(G,G+,G-) in the STEADY run case
  ! -     do not include mantles in grain/gas(H) mass ratio
  ! - lay: R_gr_scale12 = <rgrain>/<rsquare_grain> for current size distrib
  ! - G: C60 species replaced by G 

  IMPLICIT NONE
  INCLUDE "precision.f90"


  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  SUBROUTINE INITIALIZE
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    reads shock parameters, chemical species and reactions,
    !    initializes all physical variables
    ! subroutine/function needed :
    !    READ_PARAMETERS
    !    INITIALIZE_ELEMENTS
    !    READ_SPECIES
    !    CHECK_SPECIES
    !    READ_REACTIONS
    !    REACTION_TYPE
    !    CHECK_REACTIONS
    !    ADD_REVERSE_REACTIONS
    !    READ_H2_LEVELS
    !    INITIALIZE_ROVIB_H2
    !    READ_H2_RATES
    !    READ_H2_LINES
    !    READ_FE_DATA
    ! input variables :
    ! ouput variables :
    ! results :
    !   almost all physical variables are initialized here
    !---------------------------------------------------------------------------
    USE MODULE_CHEM_REACT
    USE MODULE_PHYS_VAR
    USE MODULE_VAR_VODE, ONLY : Tout_V
    USE MODULE_GRAINS
    USE MODULE_CONSTANTS, ONLY : pi, mP
    USE MODULE_H2,ONLY:READ_H2_LEVELS,INITIALIZE_ROVIB_H2,READ_H2_RATES,READ_H2_LINES
    USE MODULE_DEBUG_JLB
    USE MODULE_LINE_EXCIT
    USE MODULE_READ_FE_DATA
    IMPLICIT NONE

    INTEGER(KIND=LONG) :: i
    REAL(KIND=DP) :: Rgrain_cube ! useful for computing Dens_grain_init

    !------------------------
    !--- shock parameters ---
    !------------------------
    WRITE(*,'("parameters : ")',ADVANCE='NO')
    CALL READ_PARAMETERS
    WRITE(*,'(" read. ")')

    !---------------------------------------
    !--- chemical species                ---
    !---------------------------------------
    ! --- elemental abundances from Anders & Grevesse 1989 ---
    ! -> names, masses (g), abundances (cm-3)
    CALL INITIALIZE_ELEMENTS

    ! read in file_input + add 5 species + initialize indices ...
    WRITE(*,'("species : ")',ADVANCE='NO')
    CALL READ_SPECIES
    WRITE(*,'(" ",I3," +",I1," added ... ")',ADVANCE='NO') &
         Nspec, Nspec_plus-Nspec

    ! check the set of species
    WRITE(*,'("check.")')
    CALL CHECK_SPECIES

    !---------------------------------------
    !---       other variables           ---
    !---------------------------------------
    ! number of species in each fluid, indices at beginning and end
    DO i=1,Nspec
       IF (INDEX(speci(i)%name,'+') >0) THEN
          Nions=Nions+1
          e_ion=i
       ELSE IF (INDEX(speci(i)%name,'-') >0) THEN
          Nneg=Nneg+1
          e_neg=i
       ELSE IF (INDEX(speci(i)%name,'**') >0) THEN
          Noncores=Noncores+1
          e_cor=i
       ELSE IF (INDEX(speci(i)%name,'*') >0) THEN
          NonGrains=NonGrains+1
          e_gra=i
       ELSE
          Nneutrals=Nneutrals+1
          e_neu=i
       ENDIF
       ! verification : we must have speci(i)%index=i
       IF (speci(i)%index /= i) THEN
            print *, "  i =", i, speci(i)%index
            STOP "*** WARNING : problem in species indices"
       END IF
    END DO
    ! index at the beginning determined from the index at the end and the number of species
    b_neu=e_neu-Nneutrals+1
    b_ion=e_ion-Nions+1
    b_neg=e_neg-Nneg+1
    b_gra=e_gra-Nongrains+1
    b_cor=e_cor-Noncores+1

    ! nH=n(H)+2n(H2)+n(H+) proton density (cm-3)
    nH = speci(ind_H)%density + &
         2._DP*speci(ind_H2)%density + &
         speci(ind_Hplus)%density

    ! grains : mean square radius (cm2), density (cm-3)

    Rgrain_cube=5._DP*(Rmax_GRAIN**0.5_DP-Rmin_GRAIN**0.5_DP) / &
         (Rmin_GRAIN**(-2.5_DP)-Rmax_GRAIN**(-2.5_DP))

    Rsquare_GRAIN=5._DP*(Rmin_GRAIN**(-0.5_DP)-Rmax_GRAIN**(-0.5_DP)) / &
         (Rmin_GRAIN**(-2.5_DP)-Rmax_GRAIN**(-2.5_DP))

    R_gr_scale = Rsquare_GRAIN * Rgrain_cube**(-2.0_DP/3.0_DP)

! number of sites per grain
    Nsites_grain = 4._DP*pi*Rsquare_grain / (dsites_grain * dsites_grain)

! SC May 06: number-weighted grain radius and scale factor between <r> and sqrt<r^2> 
! needed in evolution to correct grain cross section for mantles
    r_grain = 2.5_DP*(Rmin_GRAIN**(-1.5_DP)-Rmax_GRAIN**(-1.5_DP)) / &
         (Rmin_GRAIN**(-2.5_DP)-Rmax_GRAIN**(-2.5_DP)) / 1.5_DP
    R_gr_scale12 = r_grain / sqrt(Rsquare_GRAIN)

    r_grain = sqrt(Rsquare_grain)

    MD_grain = DOT_PRODUCT(DBLE(speci(b_cor:e_cor)%density),&
         DBLE(speci(b_cor:e_cor)%mass))

! SC apr 06: mantles not included when computing dens_grain for MRN size distribution
    
    Dens_GRAIN_init = MD_grain * 3.0_DP / (4._DP*pi*Rho_GRAIN*Rgrain_cube)

    speci(ind_GRAIN)%density=Dens_GRAIN_init

    Dens_GRAIN = Dens_GRAIN_init

    if (shock_type == "S") then
! SC May 06: initialize the number of neutral grains to total grain number
       speci(ind_G)%density = Dens_GRAIN_init
       speci(ind_Gplus)%density = 1.D-20
       speci(ind_Gminus)%density = 1.D-20
     endif

    Dens_GRAIN = speci(ind_G)%density + speci(ind_Gplus)%density + &
                 speci(ind_Gminus)%density

    print*, ' n(G) + n(G+) + n(G-) =', Dens_GRAIN
    print*, ' Total ng from ** species   =', Dens_GRAIN_init

! check that input grain number densities match total grain density
    if (ABS(Dens_GRAIN/Dens_GRAIN_init - 1._DP) > 5e-2) then 
      print*, ' Grain number does not match grain mass and size distribution'
      print*, ' Run with shock type "S", Nfluids=1 to update species.in'
      stop
    endif

    ratio_GRAIN_gas = MD_grain / (nH * mP)

! SC: add to grain mass density the contribution of material in the grain mantles 

    MD_grain = MD_grain  + DOT_PRODUCT(DBLE(speci(b_gra:e_gra)%density),&
         DBLE(speci(b_gra:e_gra)%mass))

    ! number density (cm-3) of each fluid = sum(density)
    DensityN   = SUM(DBLE(speci(b_neu:e_neu)%density)) & ! neutrals
         + SUM(DBLE(speci(b_gra:e_gra)%density)) & ! neutrals on grains
         + SUM(DBLE(speci(b_cor:e_cor)%density)) ! neutrals on grain cores
    DensityI   = SUM(DBLE(speci(b_ion:e_ion)%density)) ! ions > 0
    DensityNEG = SUM(DBLE(speci(b_neg:e_neg)%density)) ! ions < 0
    speci(ind_e)%density = DensityI - DensityNEG ! electrons (charge neutrality)

    ! mass density (g.cm-3) of each fluid = sum(mass*density)
    ! include the contributions of the grains
    RhoN=DOT_PRODUCT(DBLE(speci(b_neu:e_neu)%density),&
         DBLE(speci(b_neu:e_neu)%mass)) &
         ! add neutrals on grains to neutrals
    + DOT_PRODUCT(DBLE(speci(b_gra:e_gra)%density),&
         DBLE(speci(b_gra:e_gra)%mass)) &
    + DOT_PRODUCT(DBLE(speci(b_cor:e_cor)%density),&
         DBLE(speci(b_cor:e_cor)%mass))
!   + MD_grain * speci(ind_G)%density/ &
!                      (speci(ind_Gminus)%density + speci(ind_Gplus)%density + speci(ind_G)%density)
    RhoI=DOT_PRODUCT(DBLE(speci(b_ion:e_ion)%density),&
         DBLE(speci(b_ion:e_ion)%mass))
!   + MD_grain * speci(ind_Gplus)%density/ &
!                      (speci(ind_Gminus)%density + speci(ind_Gplus)%density + speci(ind_G)%density)

    RhoNEG=DOT_PRODUCT(DBLE(speci(b_neg:e_neg)%density),&
         DBLE(speci(b_neg:e_neg)%mass))
!   + MD_grain * speci(ind_Gminus)%density/ &
!                      (speci(ind_Gminus)%density + speci(ind_Gplus)%density + speci(ind_G)%density)

    ! mean mass (g) of each fluid mu = rho/density
    muN=RhoN/DensityN
    muI=RhoI/DensityI
    muNEG=RhoNEG/DensityNEG

    ! fluid temperatures all the same at the beginning of the shock
    Ti=Tn
    Te=Tn

    ! distance into the cloud (cm) initialized to its value
    ! at the end of the first call to DRIVE
!   distance = Tout_V

    ! velocities (cm/s) of neutrals and ions
    Vn=Vs_cm
    if (Nfluids == 1) then
      DeltaV=0.0_DP
    else
      DeltaV=DeltaVmin ! = Vi-Vn, to induce the shock
    endif
    if (shock_type == 'S') then
      DeltaV = DeltaV * 10.0D3
    endif
    ABS_DeltaV=ABS(DeltaV)
    Vi=Vn-DeltaV

    ! velocity gradient (km.s-1.cm-1)
!   Vgrad = Vs_km/distance
    Vgrad = Vs_km / Tout_V

    !--------------------------------------
    !---      chemical reactions        ---
    !--------------------------------------

    WRITE(*,'("reactions : ")',ADVANCE='NO')
    ! read in file_chemistry
    CALL READ_REACTIONS
    WRITE(*,'(I4," ... ")',ADVANCE='NO')Nreact

    ! find reaction type and re-order the entire set
    ! calculate DE according to the type of the reaction
    CALL REACTION_TYPE

    ! check if the set of reactions is correct
    ! to do after REACTION_TYPE, because this subroutine needs react%type
    WRITE(*,'("check ... ")', ADVANCE='NO')
    CALL CHECK_REACTIONS

    ! add the missing endothermic reactions
    IF (do_we_add_reactions) CALL ADD_REVERSE_REACTIONS
    WRITE(*,'("+",I3," added.")')Nrever

    !---------------------------------------
    !---          H2 molecule            ---
    !---------------------------------------
    WRITE(*,'("H2 molecule : initialization.")')

    ! H2 levels
    CALL READ_H2_LEVELS

    ! populations of the levels, according to the choice of ortho:para
    CALL INITIALIZE_ROVIB_H2

    ! read collision rates for H-H2, He-H2, H2-H2
    CALL READ_H2_RATES

    ! read quadrupolar lines of H2
    CALL READ_H2_LINES

    ! read in atomic data for Fe+	
    CALL READ_FE_DATA

    ! Initialize line emissivities and integrated intensities to 0
    emicat   = 0.0_DP
    emicat_o = 0.0_DP
    intcat   = 0.0_DP
    eminat   = 0.0_DP
    eminat_o = 0.0_DP
    intnat   = 0.0_DP
    emioat   = 0.0_DP
    emioat_o = 0.0_DP
    intoat   = 0.0_DP
    emisat   = 0.0_DP
    emisat_o = 0.0_DP
    intsat   = 0.0_DP
    emisiat   = 0.0_DP
    emisiat_o = 0.0_DP
    intsiat   = 0.0_DP
    emicpl   = 0.0_DP
    emicpl_o = 0.0_DP
    intcpl   = 0.0_DP
    eminpl   = 0.0_DP
    eminpl_o = 0.0_DP
    intnpl   = 0.0_DP
    emiopl   = 0.0_DP
    emiopl_o = 0.0_DP
    intopl   = 0.0_DP
    emispl   = 0.0_DP
    emispl_o = 0.0_DP
    intspl   = 0.0_DP
    emisipl   = 0.0_DP
    emisipl_o = 0.0_DP
    intsipl   = 0.0_DP
    emifepl   = 0.0_DP
    emifepl_o = 0.0_DP
    intfepl   = 0.0_DP

  END SUBROUTINE INITIALIZE

END MODULE MODULE_INITIALIZE
