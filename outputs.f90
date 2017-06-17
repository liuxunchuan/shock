MODULE MODULE_OUTPUTS
  !*****************************************************************************
  !** The module 'MODULE_OUTPUTS' contains output subroutines (screen, files) **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  SUBROUTINE WRITE_INFO
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     to write information about shock parameters, H2 levels, grains,
    !     chemical species and reactions in the file file_info
    ! subroutine/function needed :
    !     INFO_SPECY
    !     INFO_REACTION
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CHEM_REACT
    USE MODULE_PHYS_VAR
    USE MODULE_H2
    USE MODULE_VAR_VODE, ONLY : duration_max, length_max, Eps_V, T0_V, H0_V, &
                                Tout_V, MF_V, Itask_V, Istate_V, Hdone, Hnext
    USE MODULE_GRAINS
    USE MODULE_ENERGETICS
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    IMPLICIT NONE

    INTEGER(KIND=LONG),DIMENSION(8) :: date
    INTEGER(KIND=LONG) :: i
    CHARACTER(LEN=*), PARAMETER :: name_file_info='output/info_mhd.out'
    INTEGER                     :: file_info

    ! file opening
    file_info = GET_FILE_NUMBER()
    OPEN(file_info,file=name_file_info,status='REPLACE',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    ! date, time
    CALL DATE_AND_TIME(values=date)
    WRITE(file_info,'("date : ",I2,"/",I2,"/",I4)',ADVANCE='NO')&
         date(3),date(2),date(1)
    WRITE(file_info,'(5X,"time : ",I2," H ",I2," min")')date(5),date(6)

    ! --- shock parameters ---
    WRITE(file_info,'(80("-"))')
    WRITE(file_info,'(T35,A1,"  shock")')shock_type
    WRITE(file_info,'(T35,"--------")')
    WRITE(file_info,*)
    WRITE(file_info,'("shock parameters")')
    WRITE(file_info,'(T11,"Nfluids          :",I9)')Nfluids
    WRITE(file_info,'(T11,"nH (cm-3)        :",ES9.2)')nH
    WRITE(file_info,'(T11,"Vs (km.s-1)      :",F9.1)')Vs_km
    WRITE(file_info,'(T11,"o/p H2           :",ES9.2)')op_H2
    WRITE(file_info,'(T11,"Tn (K)           :",F9.1)')Tn
    WRITE(file_info,'(T11,"Tgrain (K)       :",F9.1)')Tgrain
    WRITE(file_info,'(T11,"B (micro Gauss)  :",ES9.2)')Bfield*1.D6

    WRITE(file_info,*)
    WRITE(file_info,'("environmental parameters")')
    WRITE(file_info,'(T11,"RAD              :",F9.2)')RAD
    WRITE(file_info,'(T11,"Av (mag)         :",F9.2)')Av
    WRITE(file_info,'(T11,"Zeta (s-1)       :",ES9.2)')Zeta

    WRITE(file_info,*)
    WRITE(file_info,'("numerical parameters")')
    WRITE(file_info,'(T11,"Nstep_max        :",I9)')Nstep_max
    WRITE(file_info,'(T11,"Nstep_w          :",I9)')Nstep_w
    WRITE(file_info,'(T11,"duration_max (yr):",ES9.2)')duration_max
    WRITE(file_info,'(T11,"length_max (cm)  :",ES9.2)')length_max
    WRITE(file_info,'(T11,"XLL (cm)         :",ES9.2)',ADVANCE='NO')XLL
    IF (shock_type=='C') THEN
       WRITE(file_info,'(" (not used here)")')
    ELSE
       WRITE(file_info,*)
    ENDIF

    ! --- DVODE parameters ---
    WRITE(file_info,*)
    WRITE(file_info,'("DVODE parameters ")')
    WRITE(file_info,'(T11,"Eps_V         :",ES9.2)')Eps_V
    WRITE(file_info,'(T11,"T0_V          :",ES9.2)')T0_V
    WRITE(file_info,'(T11,"H0_V          :",ES9.2)')H0_V
    WRITE(file_info,'(T11,"Tout_V        :",ES9.2)')Tout_V
    WRITE(file_info,'(T11,"MF_V          :",I9)')MF_V
    WRITE(file_info,'(T11,"Itask_V       :",I9)')Itask_V
    WRITE(file_info,'(T11,"Istate_V      :",I9)')Istate_V
    WRITE(file_info,'(T11,"Hdone         :",ES9.2)')Hdone
    WRITE(file_info,'(T11,"Hnext         :",ES9.2)')Hnext

    !--- output choices ---
    WRITE(file_info,*)
    WRITE(file_info,'("output choices")')
    WRITE(file_info,'(T11,"speci_out      : ",A)')speci_out
    WRITE(file_info,'(T11,"H2_out         : ",A)')H2_out
    WRITE(file_info,'(T11,"line_out       : ",A)')line_out


    ! --- H2 levels ---
    WRITE(file_info,*)
    WRITE(file_info,'("H2 molecule ")')
    WRITE(file_info,'(T11,"number of levels : ",I7)')NH2_lev
    WRITE(file_info,'(T11,"E(v=",I2,",J=",I2,") (K) : ",F7.1,"  (last level)")') &
         H2_lev(NH2_lev)%V, &
         H2_lev(NH2_lev)%J, &
         H2_lev(NH2_lev)%energy
    WRITE(file_info,'(T11,"Vmax             : ", I7)')Vmax_H2
    WRITE(file_info,'(T11,"Jmax             : ", I7)')Jmax_H2
    WRITE(file_info,'(T11,"H-H2 collisions  : ", A4)')H_H2_flag
    WRITE(file_info,'(T11,"iforH2           : ", I7)')iforH2
    WRITE(file_info,'(T11,"ikinH2           : ", I7)')ikinH2
    WRITE(file_info,'(T11,"H2_int_E         : ", ES9.2," K")')H2_int_E

    ! --- grains ---
    WRITE(file_info,*)
    WRITE(file_info,'("grains")')
    WRITE(file_info,'(T11,"grain/gas        :", ES9.2," (computed from core abundances)")')ratio_GRAIN_gas
    WRITE(file_info,'(T11,"Rmin (cm)        :", ES9.2)')Rmin_grain
    WRITE(file_info,'(T11,"Rmax (cm)        :", ES9.2)')Rmax_grain
    WRITE(file_info,'(T11,"Rho (g.cm-3)     :", ES9.2)')Rho_grain

    ! --- elemental abundances ---
    WRITE(file_info,*)
    WRITE(file_info,'("elemental abundances (gas + mantles + PAH)")')
    DO i=1,Nelements
       WRITE(file_info,'(T11,A," : ",ES8.2, "   (ref : ",ES8.2,")")') &
            elements(i)%name, DBLE(elements(i)%ab_init), DBLE(elements(i)%ab_ref)
    END DO

    ! --- energetics ---
    WRITE(file_info,*)
    WRITE(file_info,'("energetics")')
    WRITE(file_info,'(T11,"mass flux (g/s/cm2)     :",ES9.2)')Mass_flux_init
    WRITE(file_info,'(T11,"momentum flux (erg/cm3) :",ES9.2)')Momentum_flux_init
    WRITE(file_info,'(T11,"energy flux (erg/s/cm2) :",ES9.2)')Energy_flux_init

    ! --- chemical species ---
    WRITE(file_info,*)
    WRITE(file_info,'(80("-"))')
    WRITE(file_info,'("-- ",I3," chemical species &
         &(+",I1," added)",T78," --")')Nspec, Nspec_plus-Nspec

    WRITE(file_info,'("-- incl. ",I2," neutrals, ",I2," on grains, &
         &", I2, " ions > 0, ", I2, " ions < 0",T78," --")') &
         Nneutrals, Nongrains, Nions, Nneg
    WRITE(file_info,'("-- name, composition, enthalpy(kCal/mol, -99.999=unknown)&
         &, Density (cm-3)",T78," --")')
    WRITE(file_info,'(80("-"))')
    DO i=1, Nspec
       CALL WRITE_SPECY(file_info,speci(i))
    ENDDO
    DO i=Nspec+1,Nspec_plus
       WRITE(file_info,'(I3,2X,A7)')speci(i)%index, speci(i)%name
    END DO

    ! --- chemical reactions ---
    WRITE(file_info,*)
    WRITE(file_info,'(80("-"))')
    WRITE(file_info,'("--- ",I4," chemical reactions",T77," ---")')Nreact
    WRITE(file_info,'("--- incl. ",I3," PHOTO, ",I3," CR_IO, ", &
         &I3," CR_DE, ",I3," H2_FO, ",I3," THREE, ",I3," SPUTT",T77," ---")') &
         &Nphoto, Ncr_io, Ncr_de, Nh2_fo, Nthree, Nsputt
    WRITE(file_info,'("---       ",I3," EROSI, ",I3," ADSOR, ",I3,&
         &" DISSO, ",I3," OTHER, and ",I3," REVER",T77," ---")') &
         &Nerosi, Nadsor, Ndisso, Nother, Nrever
    WRITE(file_info,'("--- R1 + R2 = P1 + P2 + P3 + P4, gamma (cm3.s-1), &
         &alpha, beta (K), DE (eV)",T77," ---")')
    WRITE(file_info,'(80("-"))')
    DO i=1, Nreact
       CALL WRITE_REACTION(file_info,react(i))
    END DO
    ! file closing
    CLOSE(file_info)

  END SUBROUTINE WRITE_INFO


  SUBROUTINE WRITE_SCREEN
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     to write information (counter, distance, time, temperatures) on screen
    !     during integration.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : screen
    USE MODULE_PHYS_VAR
    USE MODULE_VAR_VODE, ONLY : duration_max, length_max, Eps_V, T0_V, H0_V, Tout_V, MF_V, Itask_V, Istate_V, Hdone, Hnext
    IMPLICIT NONE
    WRITE(screen,'("step= ",I5, 2X)', ADVANCE='NO') counter / Nstep_w + 1

    WRITE(screen,'("step= ",I5, 2X)', ADVANCE='NO') counter / Nstep_w + 1
    WRITE(screen,'("z=", ES9.3, 2X)', ADVANCE='NO') distance
    WRITE(screen,'("Hdone=", ES9.3, 2X)', ADVANCE='NO') Hdone
    WRITE(screen,'("TimeN=", ES9.3, 2X)', ADVANCE='NO') timeN
    WRITE(screen,'("Tn=", ES9.3, 2X)', ADVANCE='NO') Tn
    WRITE(screen,'("Ti=", ES9.3, 2X)', ADVANCE='NO') Ti
    WRITE(screen,'("dV=", ES10.3, 2X)', ADVANCE='YES') DeltaV

  END SUBROUTINE WRITE_SCREEN



  SUBROUTINE WRITE_OUTPUT(again)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     to write several variables in corresponding output files. It opens the files, 
    !     when called for the first time, and closes the files, on the 
    !     last calculation step (again=.FALSE.).
    ! subroutine/function needed :
    !     GET_FILE_NUMBER
    ! input variables :
    !     again (logical)   -> tells if it is the last calculation step or not
    ! output variables :
    ! results :
    !
    ! 18 VII 2001 - JLB - File split as grown too large...
    !
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR
    USE MODULE_GRAINS
    USE MODULE_H2
    USE MODULE_CHEMICAL_SPECIES, ONLY : Nspec, speci, ind_H2, ind_H, ind_Hplus &
                                      , ind_Cplus, ind_C, ind_CO, ind_SiO, ind_O &
                                      ,ind_OH , ind_H2O, ind_NH3, ind_CH3OH
    USE MODULE_EVOLUTION, ONLY : B_inelastic_e_n
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_ENERGETICS
    USE MODULE_LINE_EXCIT
    IMPLICIT NONE
    INTEGER :: NCO_lev, NSiO_lev, NoH2O_lev, NpH2O_lev, NoNH3_lev, NpNH3_lev, NOH_lev &
               , Natype_lev, Netype_lev
! The values of the following PARAMETERS must be the same as those of the 
! PARAMETER (NLEVCO) in the corresponding LVG modules (CO.f, SiO.f, ... )
! See also evolution.f90 and mhd_vode.f90
    PARAMETER (NCO_lev=41)
    PARAMETER (NSiO_lev=41)
    PARAMETER (NoH2O_lev=45)
    PARAMETER (NpH2O_lev=45)
    PARAMETER (NoNH3_lev=17)
    PARAMETER (NpNH3_lev=24)
    PARAMETER (NOH_lev=20)
    PARAMETER (Natype_lev=256)
    PARAMETER (Netype_lev=256)

    LOGICAL, INTENT(in) :: again
    INTEGER(KIND=LONG) :: i, j
    CHARACTER(len=15) :: str,str2
    CHARACTER(LEN=*), PARAMETER :: format_header='(A13,1X)'
    CHARACTER(LEN=*), PARAMETER :: format_out='(ES12.3E3,1X)'
!   CHARACTER(LEN=*), PARAMETER :: format_lng='(ES17.8E3,1X)'
!   CHARACTER(LEN=*), PARAMETER :: format_lng='(ES18.9E3,1X)'
    CHARACTER(LEN=*), PARAMETER :: format_lng='(ES13.4E3,1X)'
    CHARACTER(LEN=*), PARAMETER :: name_file_phys='output/mhd_phys.out'
    INTEGER(KIND=LONG), SAVE    :: file_phys
    CHARACTER(LEN=*), PARAMETER :: name_file_speci='output/mhd_speci.out'
    INTEGER(KIND=LONG), SAVE    :: file_speci
    CHARACTER(LEN=*), PARAMETER :: name_file_H2_lev='output/H2_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_H2_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_H2_line='output/H2_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_H2_line
    CHARACTER(LEN=*), PARAMETER :: name_file_CO_line='output/CO_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_CO_line
    CHARACTER(LEN=*), PARAMETER :: name_file_CO_lev='output/CO_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_CO_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_CO_den='output/CO_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_CO_den
    CHARACTER(LEN=*), PARAMETER :: name_file_CO_cdens='output/CO_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_CO_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_CO_tau='output/CO_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_CO_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_SiO_line='output/SiO_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_SiO_line
    CHARACTER(LEN=*), PARAMETER :: name_file_SiO_lev='output/SiO_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_SiO_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_SiO_den='output/SiO_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_SiO_den
    CHARACTER(LEN=*), PARAMETER :: name_file_SiO_cdens='output/SiO_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_SiO_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_SiO_tau='output/SiO_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_SiO_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_oH2O_line='output/oH2O_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_oH2O_line
    CHARACTER(LEN=*), PARAMETER :: name_file_oH2O_lev='output/oH2O_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_oH2O_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_pH2O_line='output/pH2O_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_pH2O_line
    CHARACTER(LEN=*), PARAMETER :: name_file_pH2O_lev='output/pH2O_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_pH2O_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_oH2O_den='output/oH2O_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_oH2O_den
    CHARACTER(LEN=*), PARAMETER :: name_file_pH2O_den='output/pH2O_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_pH2O_den
    CHARACTER(LEN=*), PARAMETER :: name_file_oH2O_cdens='output/oH2O_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_oH2O_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_pH2O_cdens='output/pH2O_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_pH2O_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_oH2O_tau='output/oH2O_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_oH2O_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_pH2O_tau='output/pH2O_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_pH2O_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_oNH3_line='output/oNH3_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_oNH3_line
    CHARACTER(LEN=*), PARAMETER :: name_file_oNH3_lev='output/oNH3_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_oNH3_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_pNH3_line='output/pNH3_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_pNH3_line
    CHARACTER(LEN=*), PARAMETER :: name_file_pNH3_lev='output/pNH3_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_pNH3_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_OH_line='output/OH_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_OH_line
    CHARACTER(LEN=*), PARAMETER :: name_file_OH_lev='output/OH_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_OH_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_oNH3_den='output/oNH3_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_oNH3_den
    CHARACTER(LEN=*), PARAMETER :: name_file_pNH3_den='output/pNH3_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_pNH3_den
    CHARACTER(LEN=*), PARAMETER :: name_file_OH_den='output/OH_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_OH_den
    CHARACTER(LEN=*), PARAMETER :: name_file_oNH3_cdens='output/oNH3_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_oNH3_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_pNH3_cdens='output/pNH3_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_pNH3_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_OH_cdens='output/OH_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_OH_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_oNH3_tau='output/oNH3_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_oNH3_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_pNH3_tau='output/pNH3_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_pNH3_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_OH_tau='output/OH_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_OH_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_atype_line='output/atype_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_atype_line
    CHARACTER(LEN=*), PARAMETER :: name_file_etype_line='output/etype_line.out'
    INTEGER(KIND=LONG), SAVE    :: file_etype_line
    CHARACTER(LEN=*), PARAMETER :: name_file_atype_lev='output/atype_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_atype_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_etype_lev='output/etype_lev.out'
    INTEGER(KIND=LONG), SAVE    :: file_etype_lev
    CHARACTER(LEN=*), PARAMETER :: name_file_atype_den='output/atype_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_atype_den
    CHARACTER(LEN=*), PARAMETER :: name_file_etype_den='output/etype_den.out'
    INTEGER(KIND=LONG), SAVE    :: file_etype_den
    CHARACTER(LEN=*), PARAMETER :: name_file_atype_cdens='output/atype_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_atype_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_etype_cdens='output/etype_cdens.out'
    INTEGER(KIND=LONG), SAVE    :: file_etype_cdens
    CHARACTER(LEN=*), PARAMETER :: name_file_atype_tau='output/atype_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_atype_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_etype_tau='output/etype_tau.out'
    INTEGER(KIND=LONG), SAVE    :: file_etype_tau
    CHARACTER(LEN=*), PARAMETER :: name_file_cooling='output/cooling.out'
    INTEGER(KIND=LONG), SAVE    :: file_cooling
    CHARACTER(LEN=*), PARAMETER :: name_file_energetics='output/energetics.out'
    INTEGER(KIND=LONG), SAVE    :: file_energetics
    CHARACTER(LEN=*), PARAMETER :: name_file_intensity='output/intensity.out'
    INTEGER(KIND=LONG), SAVE    :: file_intensity
    CHARACTER(LEN=*), PARAMETER :: name_file_populations='output/populations.out'
    INTEGER(KIND=LONG), SAVE    :: file_populations
    CHARACTER(LEN=*), PARAMETER :: name_file_fe_pops='output/fe_pops.out'
    INTEGER(KIND=LONG), SAVE    :: file_fe_pops
    CHARACTER(LEN=*), PARAMETER :: name_file_fe_lines='output/fe_lines.out'
    INTEGER(KIND=LONG), SAVE    :: file_fe_lines
    CHARACTER(LEN=*), PARAMETER :: name_file_jlb='output/jlb.out'
    INTEGER(KIND=LONG), SAVE    :: file_jlb
    REAL(KIND=DP) :: dens_CO_lev,TAU_CO,TVGRAD_CO,CO_lines,CO_cdens
    REAL(KIND=DP) :: dens_SiO_lev,TAU_SiO,TVGRAD_SiO,SiO_lines,SiO_cdens
    REAL(KIND=DP) :: dens_oH2O_lev,TVGRAD_oH2O,oH2O_lines,oH2O_cdens
    REAL(KIND=DP) :: dens_pH2O_lev,TVGRAD_pH2O,pH2O_lines,pH2O_cdens
    REAL(KIND=DP) :: dens_oNH3_lev,TVGRAD_oNH3,oNH3_lines,oNH3_cdens
    REAL(KIND=DP) :: dens_pNH3_lev,TVGRAD_pNH3,pNH3_lines,pNH3_cdens
    REAL(KIND=DP) :: dens_OH_lev,TVGRAD_OH,OH_lines,OH_cdens
    REAL(KIND=DP) :: dens_atype_lev,TVGRAD_atype,atype_lines,atype_cdens
    REAL(KIND=DP) :: dens_etype_lev,TVGRAD_etype,etype_lines,etype_cdens
    REAL(KIND=DP) :: freq_oH2O,freq_pH2O &
                    ,tau_oH2O,tau_pH2O &
                    ,freq_oNH3,freq_pNH3 &
                    ,freq_OH &
                    ,tau_oNH3,tau_pNH3 &
                    ,tau_OH &
                    ,freq_atype,freq_etype &
                    ,tau_atype,tau_etype &
                    ,WWT_CO,WWT_SiO,WWT_oH2O,WWT_pH2O,WWT_oNH3,WWT_pNH3,WWT_OH,WWT_atype,WWT_etype &
                    ,G_CO,G_SiO,G_oH2O,G_pH2O,G_oNH3,G_pNH3,G_OH,G_atype,G_etype &
                    ,ELEV_CO(1:NCO_lev),ELEV_SiO(1:NSiO_lev),ELEV_oH2O(1:NoH2O_lev) &
                    ,ELEV_pH2O(1:NpH2O_lev),ELEV_oNH3(1:NoNH3_lev),ELEV_pNH3(1:NpNH3_lev) &
                    ,ELEV_OH(1:NOH_lev),NJLEV_OH(1:NOH_lev) &
                    ,ELEV_atype(1:Natype_lev),ELEV_etype(1:Netype_lev)
    REAL(KIND=DP) :: timescale_CO(1:NCO_lev),timescale_SiO(1:NSiO_lev),timescale_oH2O(1:NoH2O_lev) &
                    ,timescale_pH2O(1:NpH2O_lev),timescale_oNH3(1:NoNH3_lev) &
                    ,timescale_pNH3(1:NpNH3_lev),timescale_OH(1:NOH_lev),timescale_atype(1:Natype_lev) &
                    ,timescale_etype(1:Netype_lev),dyntime
    INTEGER :: NJLEV_CO(1:NCO_lev),NJLEV_SiO(1:NSiO_lev),NJLEV_oH2O(1:NoH2O_lev),NJLEV_pH2O(1:NpH2O_lev) &
              ,NJLEV_oNH3(1:NoNH3_lev),NJLEV_pNH3(1:NpNH3_lev) &
              ,NJLEV_atype(1:Natype_lev),NJLEV_etype(1:Netype_lev) &
              ,NKLEV_atype(1:Natype_lev),NKLEV_etype(1:Netype_lev)

      COMMON /FCTCO/ dens_CO_lev(1:NCO_lev),TAU_CO(NCO_lev,NCO_lev),WWT_CO, &
                     TVGRAD_CO(1:NCO_lev),CO_lines(1:NCO_lev),CO_cdens(1:NCO_lev)
      COMMON /FCTSiO/ dens_SiO_lev(1:NSiO_lev),TAU_SiO(NSiO_lev,NSiO_lev),WWT_SiO, &
                     TVGRAD_SiO(1:NSiO_lev),SiO_lines(1:NSiO_lev),SiO_cdens(1:NSiO_lev)
      COMMON /CO/ G_CO,ELEV_CO,NJLEV_CO
      COMMON /SiO/ G_SiO,ELEV_SiO,NJLEV_SiO
      COMMON /FCToH2O/ dens_oH2O_lev(1:NoH2O_lev),freq_oH2O(NoH2O_lev,NoH2O_lev), &
                     tau_oH2O(NoH2O_lev,NoH2O_lev),WWT_oH2O, &
                     TVGRAD_oH2O(NoH2O_lev,NoH2O_lev),oH2O_lines(NoH2O_lev,NoH2O_lev), &
                     oH2O_cdens(1:NoH2O_lev)
      COMMON /oH2O/ G_oH2O,ELEV_oH2O,NJLEV_oH2O
      COMMON /FCTpH2O/ dens_pH2O_lev(1:NpH2O_lev),freq_pH2O(NpH2O_lev,NpH2O_lev), &
                     tau_pH2O(NpH2O_lev,NpH2O_lev),WWT_pH2O, &
                     TVGRAD_pH2O(NpH2O_lev,NpH2O_lev),pH2O_lines(NpH2O_lev,NpH2O_lev), &
                     pH2O_cdens(1:NpH2O_lev)
      COMMON /pH2O/ G_pH2O,ELEV_pH2O,NJLEV_pH2O
      COMMON /FCToNH3/ dens_oNH3_lev(1:NoNH3_lev),freq_oNH3(NoNH3_lev,NoNH3_lev), &
                     tau_oNH3(NoNH3_lev,NoNH3_lev),WWT_oNH3, &
                     TVGRAD_oNH3(NoNH3_lev,NoNH3_lev),oNH3_lines(NoNH3_lev,NoNH3_lev), &
                     oNH3_cdens(1:NoNH3_lev)
      COMMON /oNH3/ G_oNH3,ELEV_oNH3,NJLEV_oNH3
      COMMON /FCTpNH3/ dens_pNH3_lev(1:NpNH3_lev),freq_pNH3(NpNH3_lev,NpNH3_lev), &
                     tau_pNH3(NpNH3_lev,NpNH3_lev),WWT_pNH3, &
                     TVGRAD_pNH3(NpNH3_lev,NpNH3_lev),pNH3_lines(NpNH3_lev,NpNH3_lev), &
                     pNH3_cdens(1:NpNH3_lev)
      COMMON /pNH3/ G_pNH3,ELEV_pNH3,NJLEV_pNH3
      COMMON /FCTOH/ dens_OH_lev(1:NOH_lev),freq_OH(NOH_lev,NOH_lev), &
                     tau_OH(NOH_lev,NOH_lev),WWT_OH, &
                     TVGRAD_OH(NOH_lev,NOH_lev),OH_lines(NOH_lev,NOH_lev), &
                     OH_cdens(1:NOH_lev)
      COMMON /OH/ G_OH,ELEV_OH,NJLEV_OH
      COMMON /FCTa_CH3OH/ dens_atype_lev(1:Natype_lev),freq_atype(Natype_lev,Natype_lev) &
                         ,tau_atype(Natype_lev,Natype_lev) &
                         ,WWT_atype,TVGRAD_atype(Natype_lev,Natype_lev) &
                         ,atype_lines(Natype_lev,Natype_lev),atype_cdens(Natype_lev)
      COMMON /Atype/ G_atype,ELEV_atype,NJLEV_atype,NKLEV_atype
      COMMON /FCTe_CH3OH/ dens_etype_lev(1:Netype_lev),freq_etype(Netype_lev,Netype_lev) &
                         ,tau_etype(Netype_lev,Netype_lev) &
                         ,WWT_etype,TVGRAD_etype(Netype_lev,Netype_lev) &
                         ,etype_lines(Netype_lev,Netype_lev),etype_cdens(Netype_lev)
      COMMON /Etype/ G_etype,ELEV_etype,NJLEV_etype,NKLEV_etype
      COMMON /tscales/ timescale_CO,timescale_SiO,timescale_oH2O,timescale_pH2O,timescale_oNH3 &
                      ,timescale_pNH3,timescale_OH,timescale_atype,timescale_etype,dyntime
    !--------------------------------------------------------------------------------
    ! When called for the first time, open the different files and write the headers
    !--------------------------------------------------------------------------------
    IF (counter == 1) THEN
       !--- MHD variables + grains ---
       !------------------------------
       file_phys = GET_FILE_NUMBER()
       OPEN(file_phys,file=name_file_phys,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       ! mhd variables
       WRITE(file_phys,format_header,ADVANCE='NO')'distance'
       WRITE(file_phys,format_header,ADVANCE='NO')'timeN'
       WRITE(file_phys,format_header,ADVANCE='NO')'timeI'
       WRITE(file_phys,format_header,ADVANCE='NO')'nH'
       WRITE(file_phys,format_header,ADVANCE='NO')'op_H2'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vn'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vi'
       WRITE(file_phys,format_header,ADVANCE='NO')'Tn'
       WRITE(file_phys,format_header,ADVANCE='NO')'Ti'
       WRITE(file_phys,format_header,ADVANCE='NO')'Te'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vsound'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vmagnet'
       WRITE(file_phys,format_header,ADVANCE='NO')'rhoN'
       WRITE(file_phys,format_header,ADVANCE='NO')'rhoI'
       WRITE(file_phys,format_header,ADVANCE='NO')'rhoNEG'
       WRITE(file_phys,format_header,ADVANCE='NO')'DensityN'
       WRITE(file_phys,format_header,ADVANCE='NO')'DensityI'
       WRITE(file_phys,format_header,ADVANCE='NO')'DensityNEG'
       WRITE(file_phys,format_header,ADVANCE='NO')'muN'
       WRITE(file_phys,format_header,ADVANCE='NO')'muI'
       WRITE(file_phys,format_header,ADVANCE='NO')'muNEG'
       WRITE(file_phys,format_header,ADVANCE='NO')'-dVn/dz'
       WRITE(file_phys,format_header,ADVANCE='NO')'-dVi/dz'
       ! grains
       WRITE(file_phys,format_header,ADVANCE='NO')'Nlayers_gr'
       WRITE(file_phys,format_header,ADVANCE='NO')'n(grain)'
       WRITE(file_phys,format_header,ADVANCE='NO')'MDens_gr'
       WRITE(file_phys,format_header,ADVANCE='NO')'MDens_gr_charge'
       WRITE(file_phys,format_header,ADVANCE='NO')'Rgrain'
       WRITE(file_phys,format_header,ADVANCE='NO')'-DissH2_t'
       WRITE(file_phys,format_header,ADVANCE='NO')'-DissH2_n'
       WRITE(file_phys,format_header,ADVANCE='NO')'-DissH2_i'
       WRITE(file_phys,format_header,ADVANCE='NO')'-DissH2_e'
       WRITE(file_phys,format_header,ADVANCE='NO')'F_gr_H2'
       WRITE(file_phys,format_header,ADVANCE='NO')'E_gr_H2_K'
       WRITE(file_phys,format_header,ADVANCE='NO')'x(H2)'
       WRITE(file_phys,format_header,ADVANCE='NO')'x(H)'
       WRITE(file_phys,format_header,ADVANCE='NO')'x(H+)'
       WRITE(file_phys,format_header,ADVANCE='NO')'sum(neut)'
       WRITE(file_phys,format_header,ADVANCE='NO')'sum(ions)'
       WRITE(file_phys,format_header,ADVANCE='NO')'sum(nega)'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vgrad'                !  Used in LVG calculations
       WRITE(file_phys,*)

       !--- chemical species ---
       !-------------------------------------------------
       file_speci = GET_FILE_NUMBER()
       OPEN(file_speci,file=name_file_speci,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       ! densities (cm-3) OR column densities (cm-2) OR fractional abundances
       WRITE(file_speci,format_header,ADVANCE='NO')'distance'
       WRITE(file_speci,format_header,ADVANCE='NO')'timeN'
       WRITE(file_speci,format_header,ADVANCE='NO')'Tn'
       WRITE(file_speci,format_header,ADVANCE='NO')'nH'
       DO i=1,Nspec
          SELECT CASE(TRIM(speci_out))
          CASE('AD')
             str = 'n(' // TRIM(speci(i)%name) // ')'
          CASE('CD')
             str = 'N(' // TRIM(speci(i)%name) // ')'
          CASE('FD')
             str = 'x(' // TRIM(speci(i)%name) // ')'
          CASE DEFAULT
             STOP "*** Error in output specification ***"
          END SELECT
          WRITE(file_speci,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_speci,*)


       !--- H2 levels and H2 lines ---
       !------------------------------
       ! populations
       file_H2_lev = GET_FILE_NUMBER()
       OPEN(file_H2_lev,file=name_file_H2_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_H2_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'Tn'
       ! densities (cm-3)
       DO i=1,NH2_lev
          WRITE(str,'(I2)') H2_lev(i)%V
          WRITE(str2,'(I2)') H2_lev(i)%J
          WRITE(file_H2_lev,format_header,ADVANCE='NO') '(' // TRIM(ADJUSTL(str)) // &
               ',' // TRIM(ADJUSTL(str2)) // ')'
       END DO
       WRITE(file_H2_lev,*)

       ! intensities (erg/s/cm2/sr)
       file_H2_line = GET_FILE_NUMBER()
       OPEN(file_H2_line,file=name_file_H2_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_H2_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_H2_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_H2_line,format_header,ADVANCE='NO')'Tn'
       DO i=1,NH2_lines_out
          WRITE(file_H2_line,format_header,ADVANCE='NO')TRIM(H2_lines(i)%name)
       END DO
       WRITE(file_H2_line,*)


       !--- CO levels and CO lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_CO_line = GET_FILE_NUMBER()
       OPEN(file_CO_line,file=name_file_CO_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_CO_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_CO_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_CO_line,format_header,ADVANCE='NO')'Tn'

       DO i=1,NCO_lev
          WRITE(str,'(I2)') NJLEV_CO(i)
          WRITE(file_CO_line,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_CO_line,*)

       ! integrated line temperatures (K km s-1)
       file_CO_lev = GET_FILE_NUMBER()
       OPEN(file_CO_lev,file=name_file_CO_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_CO_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_CO_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_CO_lev,format_header,ADVANCE='NO')'Tn'

       DO i=1,NCO_lev
          WRITE(str,'(I2)') NJLEV_CO(i)
          WRITE(file_CO_lev,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_CO_lev,*)

       ! population densities (cm-3)
       file_CO_den = GET_FILE_NUMBER()
       OPEN(file_CO_den,file=name_file_CO_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_CO_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_CO_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_CO_den,format_header,ADVANCE='NO')'Tn'

       DO i=1,NCO_lev
          WRITE(str,'(I2)') NJLEV_CO(i)
          WRITE(file_CO_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_CO_den,*)

       ! column densities (cm-2), divided by (2J+1) 
       file_CO_cdens = GET_FILE_NUMBER()
       OPEN(file_CO_cdens,file=name_file_CO_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_CO_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_CO_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_CO_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,NCO_lev
          WRITE(str,'(F10.3)') ELEV_CO(i)
          WRITE(file_CO_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_CO_cdens,*)
       
       ! line optical depths
       file_CO_tau = GET_FILE_NUMBER()
       OPEN(file_CO_tau,file=name_file_CO_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_CO_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_CO_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_CO_tau,format_header,ADVANCE='NO')'Tn'

       DO i=2,NCO_lev
          WRITE(str,'(I2)') NJLEV_CO(i)
          WRITE(file_CO_tau,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_CO_tau,*)

       !--- SiO levels and SiO lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_SiO_line = GET_FILE_NUMBER()
       OPEN(file_SiO_line,file=name_file_SiO_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_SiO_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_SiO_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_SiO_line,format_header,ADVANCE='NO')'Tn'

       DO i=1,NSiO_lev
          WRITE(str,'(I2)') NJLEV_SiO(i)
          WRITE(file_SiO_line,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_SiO_line,*)

       ! integrated line temperatures (K km s-1)
       file_SiO_lev = GET_FILE_NUMBER()
       OPEN(file_SiO_lev,file=name_file_SiO_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_SiO_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_SiO_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_SiO_lev,format_header,ADVANCE='NO')'Tn'

       DO i=1,NSiO_lev
          WRITE(str,'(I2)') NJLEV_SiO(i)
          WRITE(file_SiO_lev,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_SiO_lev,*)

       ! population densities (cm-3)
       file_SiO_den = GET_FILE_NUMBER()
       OPEN(file_SiO_den,file=name_file_SiO_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_SiO_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_SiO_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_SiO_den,format_header,ADVANCE='NO')'Tn'

       DO i=1,NSiO_lev
          WRITE(str,'(I2)') NJLEV_SiO(i)
          WRITE(file_SiO_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_SiO_den,*)

       ! column densities (cm-2), divided by (2J+1) 
       file_SiO_cdens = GET_FILE_NUMBER()
       OPEN(file_SiO_cdens,file=name_file_SiO_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_SiO_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_SiO_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_SiO_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,NSiO_lev
          WRITE(str,'(F10.3)') 1.4388*ELEV_SiO(i)
          WRITE(file_SiO_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_SiO_cdens,*)
       
       ! line optical depths
       file_SiO_tau = GET_FILE_NUMBER()
       OPEN(file_SiO_tau,file=name_file_SiO_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_SiO_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_SiO_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_SiO_tau,format_header,ADVANCE='NO')'Tn'

       DO i=2,NSiO_lev
          WRITE(str,'(I2)') NJLEV_SiO(i)
          WRITE(file_SiO_tau,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_SiO_tau,*)

       !--- o-H2O levels and o-H2O lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_oH2O_line = GET_FILE_NUMBER()
       OPEN(file_oH2O_line,file=name_file_oH2O_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oH2O_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_oH2O_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oH2O_line,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoH2O_lev
        DO j=1,NoH2O_lev
         if(freq_oH2O(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_oH2O(i,j)
          WRITE(file_oH2O_line,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_oH2O_line,*)

       ! integrated line temperatures (K km s-1)
       file_oH2O_lev = GET_FILE_NUMBER()
       OPEN(file_oH2O_lev,file=name_file_oH2O_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oH2O_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_oH2O_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oH2O_lev,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoH2O_lev
        DO j=1,NoH2O_lev
         if(freq_oH2O(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_oH2O(i,j)
          WRITE(file_oH2O_lev,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_oH2O_lev,*)

       ! population densities (cm-3) 
       file_oH2O_den = GET_FILE_NUMBER()
       OPEN(file_oH2O_den,file=name_file_oH2O_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oH2O_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_oH2O_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oH2O_den,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoH2O_lev
          WRITE(str,'(F10.3)') ELEV_oH2O(i)
          WRITE(file_oH2O_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_oH2O_den,*)

       ! column densities (cm-2), divided by (2J+1) 
       file_oH2O_cdens = GET_FILE_NUMBER()
       OPEN(file_oH2O_cdens,file=name_file_oH2O_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oH2O_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_oH2O_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oH2O_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoH2O_lev
          WRITE(str,'(F10.3)') ELEV_oH2O(i)
          WRITE(file_oH2O_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_oH2O_cdens,*)
       
       ! line optical depths 
       file_oH2O_tau = GET_FILE_NUMBER()
       OPEN(file_oH2O_tau,file=name_file_oH2O_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oH2O_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_oH2O_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oH2O_tau,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoH2O_lev
        DO j=1,NoH2O_lev
         if(freq_oH2O(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_oH2O(i,j)
          WRITE(file_oH2O_tau,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_oH2O_tau,*)
       
       !--- p-H2O levels and p-H2O lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_pH2O_line = GET_FILE_NUMBER()
       OPEN(file_pH2O_line,file=name_file_pH2O_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pH2O_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_pH2O_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pH2O_line,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpH2O_lev
        DO j=1,NpH2O_lev
         if(freq_pH2O(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_pH2O(i,j)
          WRITE(file_pH2O_line,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_pH2O_line,*)

       ! integrated line temperatures (K km s-1)
       file_pH2O_lev = GET_FILE_NUMBER()
       OPEN(file_pH2O_lev,file=name_file_pH2O_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pH2O_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_pH2O_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pH2O_lev,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpH2O_lev
        DO j=1,NpH2O_lev
         if(freq_pH2O(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_pH2O(i,j)
          WRITE(file_pH2O_lev,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_pH2O_lev,*)

       ! population densities (cm-3) 
       file_pH2O_den = GET_FILE_NUMBER()
       OPEN(file_pH2O_den,file=name_file_pH2O_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pH2O_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_pH2O_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pH2O_den,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpH2O_lev
          WRITE(str,'(F10.3)') ELEV_pH2O(i)
          WRITE(file_pH2O_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_pH2O_den,*)

       ! column densities (cm-2), divided by (2J+1) 
       file_pH2O_cdens = GET_FILE_NUMBER()
       OPEN(file_pH2O_cdens,file=name_file_pH2O_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pH2O_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_pH2O_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pH2O_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpH2O_lev
          WRITE(str,'(F10.3)') ELEV_pH2O(i)
          WRITE(file_pH2O_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_pH2O_cdens,*)
       
       ! line optical depths 
       file_pH2O_tau = GET_FILE_NUMBER()
       OPEN(file_pH2O_tau,file=name_file_pH2O_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pH2O_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_pH2O_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pH2O_tau,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpH2O_lev
        DO j=1,NpH2O_lev
         if(freq_pH2O(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_pH2O(i,j)
          WRITE(file_pH2O_tau,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_pH2O_tau,*)
       
       !--- o-NH3 levels and o-NH3 lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_oNH3_line = GET_FILE_NUMBER()
       OPEN(file_oNH3_line,file=name_file_oNH3_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oNH3_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_oNH3_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oNH3_line,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoNH3_lev
        DO j=1,NoNH3_lev
         if(freq_oNH3(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_oNH3(i,j)
          WRITE(file_oNH3_line,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_oNH3_line,*)

       ! integrated line temperatures (K km s-1)
       file_oNH3_lev = GET_FILE_NUMBER()
       OPEN(file_oNH3_lev,file=name_file_oNH3_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oNH3_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_oNH3_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oNH3_lev,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoNH3_lev
        DO j=1,NoNH3_lev
         if(freq_oNH3(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_oNH3(i,j)
          WRITE(file_oNH3_lev,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_oNH3_lev,*)

       ! population densities (cm-3) 
       file_oNH3_den = GET_FILE_NUMBER()
       OPEN(file_oNH3_den,file=name_file_oNH3_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oNH3_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_oNH3_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oNH3_den,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoNH3_lev
          WRITE(str,'(F10.3)') 1.4388*ELEV_oNH3(i)
          WRITE(file_oNH3_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_oNH3_den,*)

       ! column densities (cm-2), divided by (2J+1) 
       file_oNH3_cdens = GET_FILE_NUMBER()
       OPEN(file_oNH3_cdens,file=name_file_oNH3_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oNH3_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_oNH3_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oNH3_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoNH3_lev
          WRITE(str,'(F10.3)') 1.4388*ELEV_oNH3(i)
          WRITE(file_oNH3_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_oNH3_cdens,*)
       
       ! line optical depths 
       file_oNH3_tau = GET_FILE_NUMBER()
       OPEN(file_oNH3_tau,file=name_file_oNH3_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_oNH3_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_oNH3_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_oNH3_tau,format_header,ADVANCE='NO')'Tn'
       DO i=1,NoNH3_lev
        DO j=1,NoNH3_lev
         if(freq_oNH3(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_oNH3(i,j)
          WRITE(file_oNH3_tau,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_oNH3_tau,*)
       
       !--- p-NH3 levels and p-NH3 lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_pNH3_line = GET_FILE_NUMBER()
       OPEN(file_pNH3_line,file=name_file_pNH3_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pNH3_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_pNH3_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pNH3_line,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpNH3_lev
        DO j=1,NpNH3_lev
         if(freq_pNH3(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_pNH3(i,j)
          WRITE(file_pNH3_line,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_pNH3_line,*)

       ! integrated line temperatures (K km s-1)
       file_pNH3_lev = GET_FILE_NUMBER()
       OPEN(file_pNH3_lev,file=name_file_pNH3_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pNH3_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_pNH3_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pNH3_lev,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpNH3_lev
        DO j=1,NpNH3_lev
         if(freq_pNH3(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_pNH3(i,j)
          WRITE(file_pNH3_lev,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_pNH3_lev,*)

       ! population densities (cm-3) 
       file_pNH3_den = GET_FILE_NUMBER()
       OPEN(file_pNH3_den,file=name_file_pNH3_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pNH3_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_pNH3_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pNH3_den,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpNH3_lev
          WRITE(str,'(F10.3)') 1.4388*ELEV_pNH3(i)
          WRITE(file_pNH3_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_pNH3_den,*)

       ! column densities (cm-2), divided by (2J+1) 
       file_pNH3_cdens = GET_FILE_NUMBER()
       OPEN(file_pNH3_cdens,file=name_file_pNH3_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pNH3_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_pNH3_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pNH3_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpNH3_lev
          WRITE(str,'(F10.3)') 1.4388*ELEV_pNH3(i)
          WRITE(file_pNH3_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_pNH3_cdens,*)
       
       ! line optical depths 
       file_pNH3_tau = GET_FILE_NUMBER()
       OPEN(file_pNH3_tau,file=name_file_pNH3_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_pNH3_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_pNH3_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_pNH3_tau,format_header,ADVANCE='NO')'Tn'
       DO i=1,NpNH3_lev
        DO j=1,NpNH3_lev
         if(freq_pNH3(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_pNH3(i,j)
          WRITE(file_pNH3_tau,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_pNH3_tau,*)
       
       !--- OH levels and OH lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_OH_line = GET_FILE_NUMBER()
       OPEN(file_OH_line,file=name_file_OH_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_OH_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_OH_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_OH_line,format_header,ADVANCE='NO')'Tn'
       DO i=1,NOH_lev
        DO j=1,NOH_lev
         if(freq_OH(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_OH(i,j)
          WRITE(file_OH_line,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_OH_line,*)

       ! integrated line temperatures (K km s-1)
       file_OH_lev = GET_FILE_NUMBER()
       OPEN(file_OH_lev,file=name_file_OH_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_OH_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_OH_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_OH_lev,format_header,ADVANCE='NO')'Tn'
       DO i=1,NOH_lev
        DO j=1,NOH_lev
         if(freq_OH(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_OH(i,j)
          WRITE(file_OH_lev,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
  !    DO i=1,NOH_lev
  !       WRITE(str,'(F10.3)') 1.4388*ELEV_OH(i)
  !       WRITE(file_OH_lev,format_header,ADVANCE='NO')TRIM(str)
  !    END DO
       WRITE(file_OH_lev,*)

       ! population densities (cm-3) 
       file_OH_den = GET_FILE_NUMBER()
       OPEN(file_OH_den,file=name_file_OH_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_OH_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_OH_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_OH_den,format_header,ADVANCE='NO')'Tn'
       DO i=1,NOH_lev
          WRITE(str,'(F4.1)') 1.4388*ELEV_OH(i)
          WRITE(file_OH_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_OH_den,*)

       ! column densities (cm-2), divided by (2J+1) 
       file_OH_cdens = GET_FILE_NUMBER()
       OPEN(file_OH_cdens,file=name_file_OH_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_OH_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_OH_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_OH_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,NOH_lev
          WRITE(str,'(F10.3)') 1.4388*ELEV_OH(i)
          WRITE(file_OH_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_OH_cdens,*)
       
       ! line optical depths 
       file_OH_tau = GET_FILE_NUMBER()
       OPEN(file_OH_tau,file=name_file_OH_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_OH_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_OH_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_OH_tau,format_header,ADVANCE='NO')'Tn'
       DO i=1,NOH_lev
        DO j=1,NOH_lev
         if(freq_OH(i,j).ne.0.d0) then
          WRITE(str,'(f6.0)')freq_OH(i,j)
          WRITE(file_OH_tau,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_OH_tau,*)
       
       !--- A-type levels and A-type lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_atype_line = GET_FILE_NUMBER()
       OPEN(file_atype_line,file=name_file_atype_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_atype_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_atype_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_atype_line,format_header,ADVANCE='NO')'Tn'
       ! densities (cm-3) 
       DO i=1,Natype_lev
        DO j=1,Natype_lev
         if(freq_atype(i,j).ne.0.d0) then
          WRITE(str,'(f8.2)')freq_atype(i,j)
          WRITE(file_atype_line,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_atype_line,*)

       ! integrated line temperatures (K km s-1)
       file_atype_lev = GET_FILE_NUMBER()
       OPEN(file_atype_lev,file=name_file_atype_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_atype_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_atype_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_atype_lev,format_header,ADVANCE='NO')'Tn'
       ! densities (cm-3) 
       DO i=1,Natype_lev
        DO j=1,Natype_lev
         if(freq_atype(i,j).ne.0.d0) then
          WRITE(str,'(f8.2)')freq_atype(i,j)
          WRITE(file_atype_lev,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_atype_lev,*)

       ! densities (cm-3) 
       file_atype_den = GET_FILE_NUMBER()
       OPEN(file_atype_den,file=name_file_atype_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_atype_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_atype_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_atype_den,format_header,ADVANCE='NO')'Tn'
       DO i=1,Natype_lev
          WRITE(str,'(I3,I2)') NKLEV_atype(i), NJLEV_atype(i)
          WRITE(file_atype_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_atype_den,*)
       
       ! column densities (cm-2), divided by (2J+1) 
       file_atype_cdens = GET_FILE_NUMBER()
       OPEN(file_atype_cdens,file=name_file_atype_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_atype_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_atype_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_atype_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,Natype_lev
          WRITE(str,'(F10.3)') ELEV_atype(i)
          WRITE(file_atype_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_atype_cdens,*)
       
       ! line optical depths 
       file_atype_tau = GET_FILE_NUMBER()
       OPEN(file_atype_tau,file=name_file_atype_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_atype_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_atype_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_atype_tau,format_header,ADVANCE='NO')'Tn'
       DO i=1,Natype_lev
        DO j=1,Natype_lev
         if(freq_atype(i,j).ne.0.d0) then
          WRITE(str,'(f8.2)')freq_atype(i,j)
          WRITE(file_atype_tau,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_atype_tau,*)
       
       !--- E-type levels and E-type lines ---
       !------------------------------
       ! line temperatures (K)
       !------------------------------
       file_etype_line = GET_FILE_NUMBER()
       OPEN(file_etype_line,file=name_file_etype_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_etype_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_etype_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_etype_line,format_header,ADVANCE='NO')'Tn'
       ! densities (cm-3) 
       DO i=1,Netype_lev
        DO j=1,Netype_lev
         if(freq_etype(i,j).ne.0.d0) then
          WRITE(str,'(f8.2)')freq_etype(i,j)
          WRITE(file_etype_line,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_etype_line,*)

       ! integrated line temperatures (K km s-1)
       file_etype_lev = GET_FILE_NUMBER()
       OPEN(file_etype_lev,file=name_file_etype_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_etype_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_etype_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_etype_lev,format_header,ADVANCE='NO')'Tn'
       ! densities (cm-3) 
       DO i=1,Netype_lev
        DO j=1,Netype_lev
         if(freq_etype(i,j).ne.0.d0) then
          WRITE(str,'(f8.2)')freq_etype(i,j)
          WRITE(file_etype_lev,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_etype_lev,*)

       ! densities (cm-3) 
       file_etype_den = GET_FILE_NUMBER()
       OPEN(file_etype_den,file=name_file_etype_den,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_etype_den,format_header,ADVANCE='NO')'distance'
       WRITE(file_etype_den,format_header,ADVANCE='NO')'timeN'
       WRITE(file_etype_den,format_header,ADVANCE='NO')'Tn'
       DO i=1,Netype_lev
          WRITE(str,'(I3,I2)') NKLEV_etype(i), NJLEV_etype(i)
          WRITE(file_etype_den,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_etype_den,*)
       
       ! column densities (cm-2), divided by (2J+1) 
       file_etype_cdens = GET_FILE_NUMBER()
       OPEN(file_etype_cdens,file=name_file_etype_cdens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_etype_cdens,format_header,ADVANCE='NO')'distance'
       WRITE(file_etype_cdens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_etype_cdens,format_header,ADVANCE='NO')'Tn'
       DO i=1,Netype_lev
          WRITE(str,'(F10.3)') ELEV_etype(i)
          WRITE(file_etype_cdens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_etype_cdens,*)
       
       ! line optical depths 
       file_etype_tau = GET_FILE_NUMBER()
       OPEN(file_etype_tau,file=name_file_etype_tau,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_etype_tau,format_header,ADVANCE='NO')'distance'
       WRITE(file_etype_tau,format_header,ADVANCE='NO')'timeN'
       WRITE(file_etype_tau,format_header,ADVANCE='NO')'Tn'
       DO i=1,Netype_lev
        DO j=1,Netype_lev
         if(freq_etype(i,j).ne.0.d0) then
          WRITE(str,'(f8.2)')freq_etype(i,j)
          WRITE(file_etype_tau,format_header,ADVANCE='NO')TRIM(str)
         endif
        END DO
       END DO
       WRITE(file_etype_tau,*)
       

       !--- molecular and atomic cooling  ---
       !---          (erg/s/cm3)          ---
       !-------------------------------------
       file_cooling = GET_FILE_NUMBER()
       OPEN(file_cooling,file=name_file_cooling,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_cooling,format_header,ADVANCE='NO')'distance'
       WRITE(file_cooling,format_header,ADVANCE='NO')'timeN'
       WRITE(file_cooling,format_header,ADVANCE='NO')'Tn'
       ! H2
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(H2)'
       ! other molecules
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(13CO)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(OH)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(NH3)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(r,CO)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(v,CO)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(CO)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(CO)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(CO)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(r,o-H2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(r,p-H2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(v,H2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(H2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(oH2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(oH2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(pH2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(pH2O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(oNH3)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(oNH3)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(pNH3)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(pNH3)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(OH)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(OH)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(atype)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(atype)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'W(etype)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'G(etype)'

       ! fine structure
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(C+)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(Si+)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(C)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(Si)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(O)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(S+)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(N+)'
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(Fe+)'
    ! inelastic scattering of electrons on neutrals
    ! -- essentially on atomic and molecular hydrogen
       WRITE(file_cooling,format_header,ADVANCE='NO')'E(H+H2)'
       WRITE(file_cooling,*)

       !--- energetics ---
       !------------------
       file_energetics = GET_FILE_NUMBER()
       OPEN(file_energetics,file=name_file_energetics,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_energetics,format_header,ADVANCE='NO')'distance'
       WRITE(file_energetics,format_header,ADVANCE='NO')'timeN'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Tn'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Mass'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Mass_cons'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_kin'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_the'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_mag'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_vis'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Mom_cons'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_kin'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_the'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_int'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_mag'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_vis'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energ_cons'
       WRITE(file_energetics,*)

       !--- intensities---
       !------------------
       file_intensity = GET_FILE_NUMBER()
       OPEN(file_intensity,file=name_file_intensity,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')
 
       WRITE(file_intensity,format_header,ADVANCE='NO')'distance'
       WRITE(file_intensity,format_header,ADVANCE='NO')'timeN'
       WRITE(file_intensity,format_header,ADVANCE='NO')'Tn'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(158m)'                   ! 'C+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2324.7A)'                ! 'C+(3-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2323.5A)'                ! 'C+(4-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2328.1A)'                ! 'C+(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2326.9A)'                ! 'C+(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2325.4A)'                ! 'C+(5-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'Si+(34.8m)'                 ! 'Si+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(609.8m)'                  ! 'C(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(370.4m)'                  ! 'C(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(9850A)'                   ! 'C(4-3)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(9824A)'                   ! 'C(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'Si(129.7m)'                 ! 'Si(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'Si(68.5m)'                  ! 'Si(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(63.2m)'                   ! 'O(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(145.3m)'                  ! 'O(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(6300A)'                   ! 'O(4-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(6363A)'                   ! 'O(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'S+(6731A)'                  ! 'S+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'S+(6716A)'                  ! 'S+(3-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(205.3m)'                 ! 'N+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(121.8m)'                 ! 'N+(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(6527A)'                  ! 'N+(4-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(6548A)'                  ! 'N+(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(6583A)'                  ! 'N+(4-3)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N(5200A)'                   ! 'N(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N(5197A)'                   ! 'N(3-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'S(25.2m)'                   ! 'S(2-1)'
       WRITE(file_intensity,*)

       !--- populations---
       !------------------
       file_populations = GET_FILE_NUMBER()
       OPEN(file_populations,file=name_file_populations,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')
 
       WRITE(file_populations,format_header,ADVANCE='NO')'distance'
       WRITE(file_populations,format_header,ADVANCE='NO')'timeN'
       WRITE(file_populations,format_header,ADVANCE='NO')'Tn'
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(3P-J=0)'             ! C, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(3P-J=1)'             ! C, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(3P-J=2)'             ! C, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(1D-J=2)'             ! C, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(1S-J=0)'             ! C, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(4S-J=3/2)'           ! N, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2D-J=5/2)'           ! N, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2D-J=3/2)'           ! N, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2P-J=1/2)'           ! N, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2P-J=3/2)'           ! N, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(3P-J=2)'             ! O, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(3P-J=1)'             ! O, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(3P-J=0)'             ! O, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(1D-J=2)'             ! O, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(1S-J=0)'             ! O, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(3P-J=2)'             ! S, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(3P-J=1)'             ! S, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(3P-J=0)'             ! S, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(1D-J=2)'             ! S, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(1S-J=0)'             ! S, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(3P-J=0)'            ! Si, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(3P-J=1)'            ! Si, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(3P-J=2)'            ! Si, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(1D-J=2)'            ! Si, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(1S-J=0)'            ! Si, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(2P-J=1/2)'          ! C+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(2P-J=3/2)'          ! C+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(4P-J=1/2)'          ! C+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(4P-J=3/2)'          ! C+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(4P-J=5/2)'          ! C+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(3P-J=0)'            ! N+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(3P-J=1)'            ! N+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(3P-J=2)'            ! N+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(1D-J=2)'            ! N+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(1S-J=0)'            ! N+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(4S-J=3/2)'          ! O+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2D-J=5/2)'          ! O+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2D-J=3/2)'          ! O+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2P-J=3/2)'          ! O+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2P-J=1/2)'          ! O+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(4S-J=3/2)'          ! S+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2D-J=3/2)'          ! S+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2D-J=5/2)'          ! S+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2P-J=1/2)'          ! S+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2P-J=3/2)'          ! S+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(2P-J=1/2)'         ! Si+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(2P-J=3/2)'         ! Si+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(4P-J=1/2)'         ! Si+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(4P-J=3/2)'         ! Si+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(4P-J=5/2)'         ! Si+, lev = 5
       WRITE(file_populations,*)

       !--- Fe+ level populations ---
       !-----------------------------
       ! Relative to ground state
       file_fe_pops = GET_FILE_NUMBER()
       OPEN(file_fe_pops,file=name_file_fe_pops,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE (file_fe_pops, format_header, ADVANCE='NO')'Distance'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'timeN'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'Tn'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4P-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4P-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4P-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4P-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4P-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4P-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J13'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J11'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J11'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J1'
       WRITE (file_fe_pops,*)

       !--- [Fe II] lines ---
       !---------------------
       ! Intensities (erg//s/cm2/sr)
       file_fe_lines = GET_FILE_NUMBER()
       OPEN(file_fe_lines,file=name_file_fe_lines,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'Distance'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'timeN'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'Tn'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.248'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.275'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.271'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.279'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.295'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.298'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.321'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.328'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.534'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.600'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.644'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.664'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.677'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.711'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.745'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.798'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.800'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.810'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'17.936'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'25.988'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'35.777'
       WRITE (file_fe_lines,*)

!  JLB - 5 V 2003 specific output for HIFI cases

       file_jlb = GET_FILE_NUMBER()
       OPEN(file_jlb,file=name_file_jlb,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')
       WRITE (file_jlb, format_header, ADVANCE='NO') 'Distance'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'timeN'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'Tn'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'nH'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(H)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(H2)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(C+)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(C)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(CO)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(O)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(OH)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'n(H2O)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(H)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(H2)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(C+)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(C)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(CO)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(O)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(OH)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'N(H2O)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'C+(158m)'                   ! 'C+(2-1)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'C(609.8m)'                  ! 'C(2-1)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'C(370.4m)'                  ! 'C(3-2)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'O(63.2m)'                   ! 'O(2-1)'
       WRITE (file_jlb, format_header, ADVANCE='NO') 'O(145.3m)'                  ! 'O(3-2)'
       WRITE (file_jlb, *)

    END IF


    !-------------------------------------------
    ! MHD variables + grains
    !-------------------------------------------
    ! mhd variables
    WRITE(file_phys,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_phys,format_lng,ADVANCE='NO')timeN
    WRITE(file_phys,format_lng,ADVANCE='NO')timeI
    WRITE(file_phys,format_lng,ADVANCE='NO')nH
    WRITE(file_phys,format_lng,ADVANCE='NO')op_H2
    WRITE(file_phys,format_lng,ADVANCE='NO')Vn
    WRITE(file_phys,format_lng,ADVANCE='NO')Vi
    WRITE(file_phys,format_lng,ADVANCE='NO')Tn
    WRITE(file_phys,format_lng,ADVANCE='NO')Ti
    WRITE(file_phys,format_lng,ADVANCE='NO')Te
    WRITE(file_phys,format_lng,ADVANCE='NO')Vsound
    WRITE(file_phys,format_lng,ADVANCE='NO')Vmagnet
    WRITE(file_phys,format_lng,ADVANCE='NO')rhoN
    WRITE(file_phys,format_lng,ADVANCE='NO')rhoI
    WRITE(file_phys,format_lng,ADVANCE='NO')rhoNEG
    WRITE(file_phys,format_lng,ADVANCE='NO')DensityN
    WRITE(file_phys,format_lng,ADVANCE='NO')DensityI
    WRITE(file_phys,format_lng,ADVANCE='NO')DensityNEG
    WRITE(file_phys,format_lng,ADVANCE='NO')muN
    WRITE(file_phys,format_lng,ADVANCE='NO')muI
    WRITE(file_phys,format_lng,ADVANCE='NO')muNEG
    WRITE(file_phys,format_lng,ADVANCE='NO')-dVn
    WRITE(file_phys,format_lng,ADVANCE='NO')-dVi
    ! grains
    WRITE(file_phys,format_lng,ADVANCE='NO')Nlayers_grain
    WRITE(file_phys,format_lng,ADVANCE='NO')Dens_grain
    WRITE(file_phys,format_lng,ADVANCE='NO')MD_grain
    WRITE(file_phys,format_lng,ADVANCE='NO')Rho_GRAIN_charge
    WRITE(file_phys,format_lng,ADVANCE='NO')r_grain
    WRITE(file_phys,format_lng,ADVANCE='NO')-Sel_tot_H2
    WRITE(file_phys,format_lng,ADVANCE='NO')-Sel_tne_H2
    WRITE(file_phys,format_lng,ADVANCE='NO')-Sel_tio_H2
    WRITE(file_phys,format_lng,ADVANCE='NO')-Sel_tel_H2
    WRITE(file_phys,format_lng,ADVANCE='NO')For_gr_H2
    WRITE(file_phys,format_lng,ADVANCE='NO')H2_int_E
    WRITE(file_phys,format_lng,ADVANCE='NO')speci(ind_H2)%Density/nH
    WRITE(file_phys,format_lng,ADVANCE='NO')speci(ind_H)%Density/nH
    WRITE(file_phys,format_lng,ADVANCE='NO')speci(ind_Hplus)%Density/nH
    WRITE(file_phys,format_lng,ADVANCE='NO')SUM(v_variab(bv_neu:ev_neu))
    WRITE(file_phys,format_lng,ADVANCE='NO')SUM(v_variab(bv_ion:ev_ion))
    WRITE(file_phys,format_lng,ADVANCE='NO')SUM(v_variab(bv_neg:ev_neg))
    WRITE(file_phys,format_lng,ADVANCE='NO')Vgrad
    WRITE(file_phys,*)


    !-----------------
    ! chemical species
    !-----------------
    WRITE(file_speci,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_speci,format_lng,ADVANCE='NO')timeN
    WRITE(file_speci,format_lng,ADVANCE='NO')Tn
    WRITE(file_speci,format_lng,ADVANCE='NO')nH
    SELECT CASE(TRIM(speci_out))
    CASE('AD')
       ! densities (cm-3)
       DO i=1,Nspec
          WRITE(file_speci,format_lng,ADVANCE='NO')speci(i)%Density
       END DO
    CASE('CD')
       ! column densities (cm-2)
       DO i=1,Nspec
          WRITE(file_speci,format_lng,ADVANCE='NO')speci(i)%Col_dens
       END DO
    CASE('FD')
       ! fractional abundances
       DO i=1,Nspec
          WRITE(file_speci,format_lng,ADVANCE='NO')speci(i)%Density/nH
       END DO
    CASE DEFAULT
       STOP "*** Error in output specification ***"
    END SELECT
    WRITE(file_speci,*)


    !----------------------
    ! H2 levels + H2 lines
    !----------------------
    WRITE(file_H2_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_H2_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_H2_lev,format_lng,ADVANCE='NO') Tn

    SELECT CASE(TRIM(H2_out))
    CASE('AD')
       ! densities (cm-3)
       DO i=1,NH2_lev
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')H2_lev(i)%Density
       END DO
    CASE('CD')
       ! column densities (cm-2), divided by (2J+1)(2I+1) [I=0 (para), I=1 (ortho)]
       DO i=1,NH2_lev
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')H2_lev(i)%Col_dens
       END DO
    CASE('ln(N/g)')
       ! excitation diagram
       DO i=1,NH2_lev
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')log(H2_lev(i)%Col_dens/H2_lev(i)%weight)
       END DO
    CASE DEFAULT
       STOP "*** Error in output specification ***"
    END SELECT
    WRITE(file_H2_lev,*)

    WRITE(file_H2_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_H2_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_H2_line,format_lng,ADVANCE='NO') Tn
    SELECT CASE(TRIM(line_out))
    CASE('local')
       ! emissivities (erg/s/cm-3)
       DO i=1,NH2_lines_out
          WRITE(file_H2_line,format_lng,ADVANCE='NO')H2_lines(i)%emiss
       END DO
    CASE('integrated')
       ! intensities (erg/s/cm2/sr)
       DO i=1,NH2_lines_out
          WRITE(file_H2_line,format_lng,ADVANCE='NO')H2_lines(i)%intensity
       END DO
    CASE DEFAULT
       STOP "*** Error in output specification ***"
    END SELECT
    WRITE(file_H2_line,*)


    !----------------------
    ! CO levels 
    !----------------------
       ! line temperatures (K)
    WRITE(file_CO_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_CO_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_CO_line,format_lng,ADVANCE='NO') Tn

       DO i=1,NCO_lev
          WRITE(file_CO_line,format_lng,ADVANCE='NO')TVGRAD_CO(i)
       END DO
    WRITE(file_CO_line,*)

       ! integrated line temperatures (K km s-1) or
       ! timescales
    WRITE(file_CO_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_CO_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_CO_lev,format_lng,ADVANCE='NO') Tn

       DO i=1,NCO_lev
          WRITE(file_CO_lev,format_lng,ADVANCE='NO')CO_lines(i)
       END DO
    !  DO i=1,NCO_lev
    !     WRITE(file_CO_lev,format_lng,ADVANCE='NO')timescale_CO(i)
    !  END DO
    !     WRITE(file_CO_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_CO_lev,*)

       ! population densities (cm-3) or 
       ! fractional populations 
    WRITE(file_CO_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_CO_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_CO_den,format_lng,ADVANCE='NO') Tn

       DO i=1,NCO_lev
          WRITE(file_CO_den,format_lng,ADVANCE='NO')dens_CO_lev(i)
    !     WRITE(file_CO_den,format_lng,ADVANCE='NO')dens_CO_lev(i)/speci(ind_CO)%Density
       END DO
    WRITE(file_CO_den,*)

    ! 
    WRITE(file_CO_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_CO_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_CO_cdens,format_lng,ADVANCE='NO') Tn

       ! column densities (cm-2), divided by (2J+1)  
       DO i=1,NCO_lev
          WRITE(file_CO_cdens,format_lng,ADVANCE='NO')CO_cdens(i) &
                                                        /(2*NJLEV_CO(i) + 1)
       END DO
    WRITE(file_CO_cdens,*)

       ! line optical depths 
    WRITE(file_CO_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_CO_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_CO_tau,format_lng,ADVANCE='NO') Tn

       DO i=1,NCO_lev-1
          WRITE(file_CO_tau,format_lng,ADVANCE='NO')TAU_CO(i+1,i)
       END DO
    WRITE(file_CO_tau,*)

    !----------------------
    ! SiO levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_SiO_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_SiO_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_SiO_line,format_lng,ADVANCE='NO') Tn

       DO i=1,NSiO_lev
          WRITE(file_SiO_line,format_lng,ADVANCE='NO')TVGRAD_SiO(i)
       END DO
    WRITE(file_SiO_line,*)

       ! integrated line temperatures (K km s-1) or
       ! timescales
    WRITE(file_SiO_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_SiO_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_SiO_lev,format_lng,ADVANCE='NO') Tn

       DO i=1,NSiO_lev
          WRITE(file_SiO_lev,format_lng,ADVANCE='NO')SiO_lines(i)
       END DO
    !  DO i=1,NSiO_lev
    !     WRITE(file_SiO_lev,format_lng,ADVANCE='NO')timescale_SiO(i)
    !  END DO
    !     WRITE(file_SiO_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_SiO_lev,*)

       ! population densities (cm-3) or 
       ! fractional populations 
    WRITE(file_SiO_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_SiO_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_SiO_den,format_lng,ADVANCE='NO') Tn

       DO i=1,NSiO_lev
          WRITE(file_SiO_den,format_lng,ADVANCE='NO')dens_SiO_lev(i)
    !     WRITE(file_SiO_den,format_lng,ADVANCE='NO')dens_SiO_lev(i)/speci(ind_SiO)%Density
       END DO
    WRITE(file_SiO_den,*)

    ! 
    WRITE(file_SiO_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_SiO_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_SiO_cdens,format_lng,ADVANCE='NO') Tn

       ! column densities (cm-2), divided by (2J+1)  
       DO i=1,NSiO_lev
          WRITE(file_SiO_cdens,format_lng,ADVANCE='NO')SiO_cdens(i) &
                                                        /(2*NJLEV_SiO(i) + 1)
       END DO
    WRITE(file_SiO_cdens,*)

       ! line optical depths 
    WRITE(file_SiO_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_SiO_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_SiO_tau,format_lng,ADVANCE='NO') Tn

       DO i=1,NSiO_lev-1
          WRITE(file_SiO_tau,format_lng,ADVANCE='NO')TAU_SiO(i+1,i)
       END DO
    WRITE(file_SiO_tau,*)

    !----------------------
    ! o-H2O levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_oH2O_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oH2O_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_oH2O_line,format_lng,ADVANCE='NO') Tn

       DO i=1,NoH2O_lev
        DO j=1,NoH2O_lev
         if(freq_oH2O(i,j).ne.0.d0) then
          WRITE(file_oH2O_line,format_lng,ADVANCE='NO')TVGRAD_oH2O(i,j)
         endif
        END DO
       END DO
    WRITE(file_oH2O_line,*)

       ! integrated line temperatures (K km s-1) or
       ! product of line temperature (K) and velocity gradient (km s-1 cm-1) or
       ! timescales
    WRITE(file_oH2O_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oH2O_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_oH2O_lev,format_lng,ADVANCE='NO') Tn

       DO i=1,NoH2O_lev
        DO j=1,NoH2O_lev
         if(freq_oH2O(i,j).ne.0.d0) then
          WRITE(file_oH2O_lev,format_lng,ADVANCE='NO')oH2O_lines(i,j)
    !     WRITE(file_oH2O_lev,format_lng,ADVANCE='NO')TVGRAD_oH2O(i,j)
         endif
        END DO
       END DO
    !  DO i=1,NoH2O_lev
    !     WRITE(file_oH2O_lev,format_lng,ADVANCE='NO')timescale_oH2O(i)
    !  END DO
    !     WRITE(file_oH2O_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_oH2O_lev,*)

       ! population densities (cm-3)  
    WRITE(file_oH2O_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oH2O_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_oH2O_den,format_lng,ADVANCE='NO') Tn

       DO i=1,NoH2O_lev
          WRITE(file_oH2O_den,format_lng,ADVANCE='NO')dens_oH2O_lev(i)
       END DO
    WRITE(file_oH2O_den,*)

    ! 
    WRITE(file_oH2O_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oH2O_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_oH2O_cdens,format_lng,ADVANCE='NO') Tn

       ! column densities (cm-2), divided by (2J+1) and by  
       ! the statistical weight of the total proton spin, (2I+1) = 3
       DO i=1,NoH2O_lev
          WRITE(file_oH2O_cdens,format_lng,ADVANCE='NO')oH2O_cdens(i) &
                                                        /(2*NJLEV_oH2O(i) + 1)/3._DP
       END DO
    WRITE(file_oH2O_cdens,*)
    
       ! line optical depths
    WRITE(file_oH2O_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oH2O_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_oH2O_tau,format_lng,ADVANCE='NO') Tn

       DO i=1,NoH2O_lev
        DO j=1,NoH2O_lev
         if(freq_oH2O(i,j).ne.0.d0) then
          WRITE(file_oH2O_tau,format_lng,ADVANCE='NO')tau_oH2O(i,j)
         endif
        END DO
       END DO
    WRITE(file_oH2O_tau,*)

    !----------------------
    ! p-H2O levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_pH2O_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pH2O_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_pH2O_line,format_lng,ADVANCE='NO') Tn

       DO i=1,NpH2O_lev
        DO j=1,NpH2O_lev
         if(freq_pH2O(i,j).ne.0.d0) then
          WRITE(file_pH2O_line,format_lng,ADVANCE='NO')TVGRAD_pH2O(i,j)
         endif
        END DO
       END DO
    WRITE(file_pH2O_line,*)

       ! integrated line temperatures (K km s-1) or
       ! product of line temperature (K) and velocity gradient (km s-1 cm-1) or
       ! timescales
    WRITE(file_pH2O_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pH2O_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_pH2O_lev,format_lng,ADVANCE='NO') Tn

       DO i=1,NpH2O_lev
        DO j=1,NpH2O_lev
         if(freq_pH2O(i,j).ne.0.d0) then
          WRITE(file_pH2O_lev,format_lng,ADVANCE='NO')pH2O_lines(i,j)
    !     WRITE(file_pH2O_lev,format_lng,ADVANCE='NO')TVGRAD_pH2O(i,j)
         endif
        END DO
       END DO
    !  DO i=1,NpH2O_lev
    !     WRITE(file_pH2O_lev,format_lng,ADVANCE='NO')timescale_pH2O(i)
    !  END DO
    !     WRITE(file_pH2O_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_pH2O_lev,*)

       ! population densities (cm-3)  
    WRITE(file_pH2O_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pH2O_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_pH2O_den,format_lng,ADVANCE='NO') Tn

       DO i=1,NpH2O_lev
          WRITE(file_pH2O_den,format_lng,ADVANCE='NO')dens_pH2O_lev(i)
       END DO
    WRITE(file_pH2O_den,*)

       ! column densities (cm-2), divided by (2J+1)  
    ! 
    WRITE(file_pH2O_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pH2O_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_pH2O_cdens,format_lng,ADVANCE='NO') Tn

       DO i=1,NpH2O_lev
          WRITE(file_pH2O_cdens,format_lng,ADVANCE='NO')pH2O_cdens(i) &
                                                        /(2*NJLEV_pH2O(i) + 1)
       END DO
    WRITE(file_pH2O_cdens,*)

       ! line optical depths
    WRITE(file_pH2O_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pH2O_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_pH2O_tau,format_lng,ADVANCE='NO') Tn

       DO i=1,NpH2O_lev
        DO j=1,NpH2O_lev
         if(freq_pH2O(i,j).ne.0.d0) then
          WRITE(file_pH2O_tau,format_lng,ADVANCE='NO')tau_pH2O(i,j)
         endif
        END DO
       END DO
    WRITE(file_pH2O_tau,*)

    !----------------------
    ! o-NH3 levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_oNH3_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oNH3_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_oNH3_line,format_lng,ADVANCE='NO') Tn

       DO i=1,NoNH3_lev
        DO j=1,NoNH3_lev
         if(freq_oNH3(i,j).ne.0.d0) then
          WRITE(file_oNH3_line,format_lng,ADVANCE='NO')TVGRAD_oNH3(i,j)
         endif
        END DO
       END DO
    WRITE(file_oNH3_line,*)

       ! integrated line temperatures (K km s-1) or
       ! product of line temperature (K) and velocity gradient (km s-1 cm-1) or
       ! timescales
    WRITE(file_oNH3_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oNH3_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_oNH3_lev,format_lng,ADVANCE='NO') Tn

       DO i=1,NoNH3_lev
        DO j=1,NoNH3_lev
         if(freq_oNH3(i,j).ne.0.d0) then
          WRITE(file_oNH3_lev,format_lng,ADVANCE='NO')oNH3_lines(i,j)
    !     WRITE(file_oNH3_lev,format_lng,ADVANCE='NO')TVGRAD_oNH3(i,j)
         endif
        END DO
       END DO
    !  DO i=1,NoNH3_lev
    !     WRITE(file_oNH3_lev,format_lng,ADVANCE='NO')timescale_oNH3(i)
    !  END DO
    !     WRITE(file_oNH3_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_oNH3_lev,*)

       ! population densities (cm-3)  
    WRITE(file_oNH3_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oNH3_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_oNH3_den,format_lng,ADVANCE='NO') Tn

       DO i=1,NoNH3_lev
          WRITE(file_oNH3_den,format_lng,ADVANCE='NO')dens_oNH3_lev(i)
       END DO
    WRITE(file_oNH3_den,*)

    ! 
    WRITE(file_oNH3_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oNH3_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_oNH3_cdens,format_lng,ADVANCE='NO') Tn

       ! column densities (cm-2), divided by (2J+1)  
       DO i=1,NoNH3_lev
          WRITE(file_oNH3_cdens,format_lng,ADVANCE='NO')oNH3_cdens(i) &
                                                        /(2*NJLEV_oNH3(i) + 1)
       END DO
    WRITE(file_oNH3_cdens,*)
    
       ! line optical depths
    WRITE(file_oNH3_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_oNH3_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_oNH3_tau,format_lng,ADVANCE='NO') Tn

       DO i=1,NoNH3_lev
        DO j=1,NoNH3_lev
         if(freq_oNH3(i,j).ne.0.d0) then
          WRITE(file_oNH3_tau,format_lng,ADVANCE='NO')tau_oNH3(i,j)
         endif
        END DO
       END DO
    WRITE(file_oNH3_tau,*)

    !----------------------
    ! p-NH3 levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_pNH3_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pNH3_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_pNH3_line,format_lng,ADVANCE='NO') Tn

       DO i=1,NpNH3_lev
        DO j=1,NpNH3_lev
         if(freq_pNH3(i,j).ne.0.d0) then
          WRITE(file_pNH3_line,format_lng,ADVANCE='NO')TVGRAD_pNH3(i,j)
         endif
        END DO
       END DO
    WRITE(file_pNH3_line,*)

       ! integrated line temperatures (K km s-1) or
       ! product of line temperature (K) and velocity gradient (km s-1 cm-1) or
       ! timescales
    WRITE(file_pNH3_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pNH3_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_pNH3_lev,format_lng,ADVANCE='NO') Tn

       DO i=1,NpNH3_lev
        DO j=1,NpNH3_lev
         if(freq_pNH3(i,j).ne.0.d0) then
          WRITE(file_pNH3_lev,format_lng,ADVANCE='NO')pNH3_lines(i,j)
    !     WRITE(file_pNH3_lev,format_lng,ADVANCE='NO')TVGRAD_pNH3(i,j)
         endif
        END DO
       END DO
    !  DO i=1,NpNH3_lev
    !     WRITE(file_pNH3_lev,format_lng,ADVANCE='NO')timescale_pNH3(i)
    !  END DO
    !     WRITE(file_pNH3_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_pNH3_lev,*)

       ! population densities (cm-3)  
    WRITE(file_pNH3_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pNH3_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_pNH3_den,format_lng,ADVANCE='NO') Tn

       DO i=1,NpNH3_lev
          WRITE(file_pNH3_den,format_lng,ADVANCE='NO')dens_pNH3_lev(i)
       END DO
    WRITE(file_pNH3_den,*)

       ! column densities (cm-2), divided by (2J+1)  
    ! 
    WRITE(file_pNH3_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pNH3_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_pNH3_cdens,format_lng,ADVANCE='NO') Tn

       DO i=1,NpNH3_lev
          WRITE(file_pNH3_cdens,format_lng,ADVANCE='NO')pNH3_cdens(i) &
                                                        /(2*NJLEV_pNH3(i) + 1)
       END DO
    WRITE(file_pNH3_cdens,*)

       ! line optical depths
    WRITE(file_pNH3_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_pNH3_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_pNH3_tau,format_lng,ADVANCE='NO') Tn

       DO i=1,NpNH3_lev
        DO j=1,NpNH3_lev
         if(freq_pNH3(i,j).ne.0.d0) then
          WRITE(file_pNH3_tau,format_lng,ADVANCE='NO')tau_pNH3(i,j)
         endif
        END DO
       END DO
    WRITE(file_pNH3_tau,*)

    !----------------------
    ! OH levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_OH_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_OH_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_OH_line,format_lng,ADVANCE='NO') Tn

       DO i=1,NOH_lev
        DO j=1,NOH_lev
         if(freq_OH(i,j).ne.0.d0) then
          WRITE(file_OH_line,format_lng,ADVANCE='NO')TVGRAD_OH(i,j)
         endif
        END DO
       END DO
    WRITE(file_OH_line,*)

       ! integrated line temperatures (K km s-1) or
       ! product of line temperature (K) and velocity gradient (km s-1 cm-1) or
       ! timescales
    WRITE(file_OH_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_OH_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_OH_lev,format_lng,ADVANCE='NO') Tn

       DO i=1,NOH_lev
        DO j=1,NOH_lev
         if(freq_OH(i,j).ne.0.d0) then
          WRITE(file_OH_lev,format_lng,ADVANCE='NO')OH_lines(i,j)
    !     WRITE(file_OH_lev,format_lng,ADVANCE='NO')TVGRAD_OH(i,j)
         endif
        END DO
       END DO
    !  DO i=1,NOH_lev
    !     WRITE(file_OH_lev,format_lng,ADVANCE='NO')timescale_OH(i)
    !  END DO
    !     WRITE(file_OH_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_OH_lev,*)

       ! population densities (cm-3)  
    WRITE(file_OH_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_OH_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_OH_den,format_lng,ADVANCE='NO') Tn

       DO i=1,NOH_lev
          WRITE(file_OH_den,format_lng,ADVANCE='NO')dens_OH_lev(i)
       END DO
    WRITE(file_OH_den,*)

       ! column densities (cm-2), divided by (2J+1)  
    ! 
    WRITE(file_OH_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_OH_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_OH_cdens,format_lng,ADVANCE='NO') Tn

       DO i=1,NOH_lev
          WRITE(file_OH_cdens,format_lng,ADVANCE='NO')OH_cdens(i) &
                                                        /(2*NJLEV_OH(i) + 1)
       END DO
    WRITE(file_OH_cdens,*)

       ! line optical depths
    WRITE(file_OH_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_OH_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_OH_tau,format_lng,ADVANCE='NO') Tn

       DO i=1,NOH_lev
        DO j=1,NOH_lev
         if(freq_OH(i,j).ne.0.d0) then
          WRITE(file_OH_tau,format_lng,ADVANCE='NO')tau_OH(i,j)
         endif
        END DO
       END DO
    WRITE(file_OH_tau,*)

    !----------------------
    ! A-type levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_atype_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_atype_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_atype_line,format_lng,ADVANCE='NO') Tn

       DO i=1,Natype_lev
        DO j=1,Natype_lev
         if(freq_atype(i,j).ne.0.d0) then
          WRITE(file_atype_line,format_lng,ADVANCE='NO')TVGRAD_atype(i,j)
         endif
        END DO
       END DO
    WRITE(file_atype_line,*)

    WRITE(file_atype_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_atype_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_atype_lev,format_lng,ADVANCE='NO') Tn

       ! integrated line temperatures (K km s-1) or
       ! product of line temperature (K) and velocity gradient (km s-1 cm-1) or
       ! timescales
       DO i=1,Natype_lev
        DO j=1,Natype_lev
         if(freq_atype(i,j).ne.0.d0) then
          WRITE(file_atype_lev,format_lng,ADVANCE='NO')atype_lines(i,j)
    !     WRITE(file_atype_lev,format_lng,ADVANCE='NO')TVGRAD_atype(i,j)
         endif
        END DO
       END DO
    !  DO i=1,Natype_lev
    !     WRITE(file_atype_lev,format_lng,ADVANCE='NO')timescale_atype(i)
    !  END DO
    !     WRITE(file_atype_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_atype_lev,*)

    WRITE(file_atype_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_atype_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_atype_den,format_lng,ADVANCE='NO') Tn

       ! densities (cm-3)  
       DO i=1,Natype_lev
          WRITE(file_atype_den,format_lng,ADVANCE='NO')dens_atype_lev(i)
       END DO
    WRITE(file_atype_den,*)
    ! 
    WRITE(file_atype_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_atype_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_atype_cdens,format_lng,ADVANCE='NO') Tn

       ! column densities (cm-2), divided by (2J+1)  
       DO i=1,Natype_lev
          WRITE(file_atype_cdens,format_lng,ADVANCE='NO')atype_cdens(i) &
                                                        /(2*NJLEV_atype(i) + 1)
       END DO
    WRITE(file_atype_cdens,*)
    ! 
    WRITE(file_atype_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_atype_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_atype_tau,format_lng,ADVANCE='NO') Tn

       ! optical depths in the lines 
       DO i=1,Natype_lev
        DO j=1,Natype_lev
         if(freq_atype(i,j).ne.0.d0) then
          WRITE(file_atype_tau,format_lng,ADVANCE='NO')tau_atype(i,j)
         endif
        END DO
       END DO
    WRITE(file_atype_tau,*)

    !----------------------
    ! E-type levels 
    !----------------------
       ! line temperatures (K) 
    WRITE(file_etype_line,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_etype_line,format_lng,ADVANCE='NO') timeN
    WRITE(file_etype_line,format_lng,ADVANCE='NO') Tn

       DO i=1,Netype_lev
        DO j=1,Netype_lev
         if(freq_etype(i,j).ne.0.d0) then
          WRITE(file_etype_line,format_lng,ADVANCE='NO')TVGRAD_etype(i,j)
         endif
        END DO
       END DO
    WRITE(file_etype_line,*)

    WRITE(file_etype_lev,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_etype_lev,format_lng,ADVANCE='NO') timeN
    WRITE(file_etype_lev,format_lng,ADVANCE='NO') Tn

       ! integrated line temperatures (K km s-1) or
       ! product of line temperature (K) and velocity gradient (km s-1 cm-1) or
       ! timescales
       DO i=1,Netype_lev
        DO j=1,Netype_lev
         if(freq_etype(i,j).ne.0.d0) then
          WRITE(file_etype_lev,format_lng,ADVANCE='NO')etype_lines(i,j)
    !     WRITE(file_etype_lev,format_lng,ADVANCE='NO')TVGRAD_etype(i,j)
         endif
        END DO
       END DO
    !  DO i=1,Netype_lev
    !     WRITE(file_etype_lev,format_lng,ADVANCE='NO')timescale_etype(i)
    !  END DO
    !     WRITE(file_etype_lev,format_lng,ADVANCE='NO')dyntime
    WRITE(file_etype_lev,*)

    WRITE(file_etype_den,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_etype_den,format_lng,ADVANCE='NO') timeN
    WRITE(file_etype_den,format_lng,ADVANCE='NO') Tn

       ! densities (cm-3)  
       DO i=1,Netype_lev
          WRITE(file_etype_den,format_lng,ADVANCE='NO')dens_etype_lev(i)
       END DO
    WRITE(file_etype_den,*)
    ! 
    WRITE(file_etype_cdens,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_etype_cdens,format_lng,ADVANCE='NO') timeN
    WRITE(file_etype_cdens,format_lng,ADVANCE='NO') Tn

       ! column densities (cm-2), divided by (2J+1)  
       DO i=1,Netype_lev
          WRITE(file_etype_cdens,format_lng,ADVANCE='NO')etype_cdens(i) &
                                                        /(2*NJLEV_etype(i) + 1)
       END DO
    WRITE(file_etype_cdens,*)
    ! 
    WRITE(file_etype_tau,'(ES18.9E3,1X)',ADVANCE='NO') distance

    WRITE(file_etype_tau,format_lng,ADVANCE='NO') timeN
    WRITE(file_etype_tau,format_lng,ADVANCE='NO') Tn

       ! optical depths in the lines 
       DO i=1,Netype_lev
        DO j=1,Netype_lev
         if(freq_etype(i,j).ne.0.d0) then
          WRITE(file_etype_tau,format_lng,ADVANCE='NO')tau_etype(i,j)
         endif
        END DO
       END DO
    WRITE(file_etype_tau,*)


    !---------------------------------------------------
    ! molecular and atomic cooling rates (erg/s/cm3)
    !---------------------------------------------------
    WRITE(file_cooling,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_cooling,format_lng,ADVANCE='NO')timeN
    WRITE(file_cooling,format_lng,ADVANCE='NO')Tn
    ! H2
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_H2
    ! other molecules
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_13CO
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_OH
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_NH3
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_rot_CO
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_vib_CO
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_CO
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_CO
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_CO
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_rot_o_H2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_rot_p_H2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_vib_H2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')cooling_H2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_oH2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_oH2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_pH2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_pH2O
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_oNH3
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_oNH3
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_pNH3
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_pNH3
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_OH
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_OH
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_atype
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_atype
    WRITE(file_cooling,format_lng,ADVANCE='NO')WWT_etype
    WRITE(file_cooling,format_lng,ADVANCE='NO')G_etype

    ! fine structure (ions, atoms)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(emicpl)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(emisipl)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(emicat)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(emisiat)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(emioat)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(emispl)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(eminpl)
    WRITE(file_cooling,format_lng,ADVANCE='NO')SUM(emifepl)
    ! inelastic scattering of electrons on neutrals
    ! -- essentially on atomic and molecular hydrogen
    WRITE(file_cooling,format_lng,ADVANCE='NO')B_inelastic_e_n
    WRITE(file_cooling,*)


    !------------
    ! energetics
    !------------
    WRITE(file_energetics,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_energetics,format_lng,ADVANCE='NO')timeN
    WRITE(file_energetics,format_lng,ADVANCE='NO')Tn
    WRITE(file_energetics,format_lng,ADVANCE='NO')Mass_flux
    WRITE(file_energetics,format_lng,ADVANCE='NO')Mass_cons
    WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_kin
    WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_the
    WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_mag
    WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_vis
    WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux
    WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_cons
    WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_kin
    WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_the
    WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_int
    WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_mag
    WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_vis
    WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux
    WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_cons
    WRITE(file_energetics,*)

    !------------
    ! intensities
    !------------

    WRITE(file_intensity,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_intensity,format_lng,ADVANCE='NO')timeN
    WRITE(file_intensity,format_lng,ADVANCE='NO')Tn
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(1)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(5)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(6)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(2)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(3)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(4)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intsipl(1)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(1)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(2)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(4)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(5)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intsiat(1)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intsiat(2)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(2)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(1)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(6)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(5)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intspl(7)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intspl(8)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(1)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(2)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(6)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(5)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(4)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intnat(6)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intnat(7)
    WRITE(file_intensity,format_lng,ADVANCE='NO')intsat(2)
    WRITE(file_intensity,*)

    !------------
    ! populations
    !------------

    WRITE(file_populations,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_populations,format_lng,ADVANCE='NO')timeN
    WRITE(file_populations,format_lng,ADVANCE='NO')Tn
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(5)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(1)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(2)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(3)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(4)
    WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(5)
    WRITE(file_populations,*)

    !--- Fe+ level populations ---
    !-----------------------------
    WRITE(file_fe_pops,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_fe_pops,format_lng,ADVANCE='NO')timeN
    WRITE(file_fe_pops,format_lng,ADVANCE='NO')Tn
    DO i=1,nlvfepl
      WRITE(file_fe_pops,format_lng,ADVANCE='NO')pop_fepl(i)
    END DO
    WRITE (file_fe_pops,*)

    !--- [Fe II] lines ---
    !---------------------
    WRITE(file_fe_lines,'(ES18.9E3,1X)',ADVANCE='NO')distance
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')timeN
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')Tn
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(35)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(19)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(43)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(36)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(27)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(37)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(20)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(28)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(29)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(38)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(22)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(44)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(30)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(39)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(45)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(40)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(31)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(23)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(10)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(1)
    WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(18)
    WRITE (file_fe_lines,*)

! JLB - 5 V 2003

       WRITE (file_jlb, format_out, ADVANCE='NO') distance
       WRITE (file_jlb, format_out, ADVANCE='NO') timeN
       WRITE (file_jlb, format_out, ADVANCE='NO') Tn
       WRITE (file_jlb, format_out, ADVANCE='NO') nH
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_H)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_H2)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_Cplus)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_C)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_CO)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_O)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_OH)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_H2O)%Density
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_H)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_H2)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_Cplus)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_C)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_CO)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_O)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_OH)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') speci(ind_H2O)%Col_dens
       WRITE (file_jlb, format_out, ADVANCE='NO') intcpl(1)
       WRITE (file_jlb, format_out, ADVANCE='NO') intcat(1)
       WRITE (file_jlb, format_out, ADVANCE='NO') intcat(2)
       WRITE (file_jlb, format_out, ADVANCE='NO') intoat(2)
       WRITE (file_jlb, format_out, ADVANCE='NO') intoat(1)
       WRITE (file_jlb, *)

    !--------------------------------------------------------
    ! if this is the last calculation step, close the files
    !--------------------------------------------------------
    IF (.NOT.(again)) THEN
       CLOSE(file_phys)
       CLOSE(file_speci)
       CLOSE(file_H2_lev)
       CLOSE(file_H2_line)
       CLOSE(file_cooling)
       CLOSE(file_energetics)
       CLOSE(file_intensity)
       CLOSE(file_populations)
       CLOSE(file_fe_pops)
       CLOSE(file_fe_lines)
       CLOSE(file_jlb)
    END IF

  END SUBROUTINE WRITE_OUTPUT


  SUBROUTINE WRITE_EXCIT
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     to write results required for the H2 excitation diagram
    ! subroutine/function needed :
    !     GET_FILE_NUMBER
    ! input variables :
    !     none
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR, ONLY : NH2_lev
    USE MODULE_H2
    IMPLICIT NONE

    INTEGER(KIND=LONG) :: i
    CHARACTER(len=15) :: str,str2
    CHARACTER(LEN=*), PARAMETER :: format_head='(4X,"V   J   Energy(K)   log(N/g)")'
    CHARACTER(LEN=*), PARAMETER :: format_i='(i5,1X)'
    CHARACTER(LEN=*), PARAMETER :: format_r='(ES15.6E3,1X)'
    CHARACTER(LEN=*), PARAMETER :: name_file='output/excit.out'

       OPEN(30,file=name_file,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    WRITE(30,format_head)

    do i=1,NH2_lev

       WRITE(30,format_i,ADVANCE='NO') H2_lev(i)%V
       WRITE(30,format_i,ADVANCE='NO') H2_lev(i)%J
       WRITE(30,format_r,ADVANCE='NO') H2_lev(i)%Energy
       WRITE(30,format_r,ADVANCE='NO') log(H2_lev(i)%Col_dens/H2_lev(i)%weight)
       WRITE(30,*)

    end do

       CLOSE(30)

  END SUBROUTINE WRITE_EXCIT

  SUBROUTINE WRITE_SPECIES
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     to write final abundances (in same format as for input file)
    ! subroutine/function needed :
    !     GET_FILE_NUMBER
    ! input variables :
    !     none
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: name_file_in='input/species.in'
    CHARACTER(len=*), PARAMETER :: name_file_out='output/species.out'
    INTEGER(KIND=LONG) :: file_in, file_out
    CHARACTER(len=80) :: charact
    CHARACTER(len=28) :: beg_l
    CHARACTER(len=34) :: end_l
    INTEGER :: i
    REAL :: toto

    file_in = GET_FILE_NUMBER()
    OPEN(file_in,file=name_file_in,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    file_out = GET_FILE_NUMBER()
    OPEN(file_out,file=name_file_out,status='REPLACE',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')

     READ(file_in,'(A80)') charact
     WRITE(file_out,'(A55,F7.2,A18)') "!---- list of chemical species --- Steady state at T = ", Tn, " K ---------------"
    DO i=1,5
       READ(file_in,'(A80)') charact
       WRITE(file_out,'(A80)') charact
    END DO

    DO i=1,Nspec
       READ(file_in,'(A28,E10.3,A34)') beg_l, toto, end_l
       WRITE(file_out,'(A28,1P,D10.3,A34)') beg_l, speci(i)%Density/nH, end_l
    END DO

  END SUBROUTINE WRITE_SPECIES

END MODULE MODULE_OUTPUTS

