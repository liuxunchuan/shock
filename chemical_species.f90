MODULE MODULE_CHEMICAL_SPECIES
  !*****************************************************************************
  !** The module 'MODULE_CHEMICAL_SPECIES' contains variables and subroutines **
  !** related to the set of chemical species :                                **
  !**     * the data type TYPE_SPECY                                          **
  !**     * the vector containing the species                                 **
  !**     * the numbers of species in each fluid                              **
  !**     * the index useful to find one species in the vector                **
  !**     * subroutine to write information on one species                    **
  !**     * function to calculate the chemical formula of one species         **
  !**     * subroutine to calculate the elemental abundances                  **
  !*****************************************************************************
  ! SC May 06: 
  ! - C60 replaced by G everywhere. 
  ! - Element 'G' with mass = 1 amu added in INITIALIZE_ELEMENTS
  ! - Nelements = 12 instead of 11

  IMPLICIT NONE
  INCLUDE "precision.f90"


  INTEGER(KIND=LONG), PARAMETER :: name_length = 7     ! (7: old, 8: new) length of species or element name
  INTEGER(KIND=LONG), PARAMETER :: Nelements   = 12    ! number of 'basic' elements: 9 + 1(Mg) + 1(D) + 1(G)

  !---------------------------
  ! data type for one element
  !---------------------------
  TYPE TYPE_ELEMENT
     CHARACTER(len=name_length) :: name    ! name of the element
     REAL(KIND=DP)              :: mass    ! mass of the element (g)
     REAL(KIND=DP)              :: ab      ! abundance
     REAL(KIND=DP)              :: ab_init ! initial abundance
     REAL(KIND=DP)              :: ab_ref  ! abundance (from Anders & Grevesse 1989)
  END TYPE TYPE_ELEMENT
  ! elements are 'basic' components : H, C, N ..., initialized in INITIALIZE_ELEMENTS
  TYPE (TYPE_ELEMENT),DIMENSION(Nelements) :: elements
  INTEGER(KIND=LONG) :: ind_elem_H ! index of 'H' in the vector 'elements'

  !--------------------------
  ! data type for one specy
  !--------------------------
  TYPE TYPE_SPECY
     CHARACTER(LEN=name_length) :: name                  ! name (e.g. 'SiOH+')
     REAL(KIND=DP)              :: density               ! density (cm-3)
     REAL(KIND=DP)              :: Dens_old              ! density at last call to DRIVE
     REAL(KIND=DP)              :: Col_dens              ! column density (cm-2)
     REAL(KIND=DP)              :: mass                  ! mass (g)
     REAL(KIND=DP)              :: enthalpy              ! enthalpy of formation (kCal/mol)
     REAL(KIND=DP)              :: velocity              ! velocity (cm/s)
     REAL(KIND=DP)              :: temperature           ! temperature (K)
     INTEGER(KIND=LONG), DIMENSION(Nelements) :: formula ! chemical formula
     INTEGER(KIND=LONG), DIMENSION(12) :: useless        ! (old system for chemical formula)
     INTEGER(KIND=LONG)         :: index                 ! index of the species
  END TYPE TYPE_SPECY

  !--------------------------------------------------
  ! number of species in the model and in each fluid
  !--------------------------------------------------
  ! variables read in READ_SPECIES
  INTEGER(KIND=LONG) :: Nspec      ! number of chemical species
  INTEGER(KIND=LONG) :: Nspec_plus ! Nspec + 5 (e-, photon, CRP, grain, SECPHO)
  ! calculated in INITIALIZE
  INTEGER(KIND=LONG) :: Nneutrals=0   ! number of neutrals
  INTEGER(KIND=LONG) :: Nions=0       ! number of positive ions
  INTEGER(KIND=LONG) :: Nneg=0        ! number of negative ions
  INTEGER(KIND=LONG) :: Nongrains=0   ! number of species in grain mantle (name with a '*')
  INTEGER(KIND=LONG) :: Noncores=0    ! number of species in grain cores (name with a '**')


  !---------------------------------------------------------------------
  ! vector containing the species, read in READ_SPECIES
  ! species commences at index zero, as tab_specy(none=0)%name is ''
  ! the species are between the index 1 and Nspec
  ! between Nspec+1 and Nspec_plus, we find added species
  ! (electrons, grains,...)
  !---------------------------------------------------------------------
  TYPE(type_specy), SAVE, DIMENSION(:), ALLOCATABLE :: speci

  !--------------------------------------------------
  ! index of some species
  ! remark : CP -> C+, HP -> H+, SiP -> Si+
  !          SP -> S+, NP -> N+
  ! initialized in READ_SPECIES
  !--------------------------------------------------
  INTEGER(KIND=LONG) :: ind_H=0, ind_H2=0, ind_He=0, ind_O=0, ind_Oplus, ind_N=0, ind_C=0
  INTEGER(KIND=LONG) :: ind_S=0, ind_Si=0, ind_H2O=0, ind_OH=0, ind_CO=0, ind_NH3=0, ind_CH3OH=0
  INTEGER(KIND=LONG) :: ind_SiO=0, ind_Cplus=0, ind_Hplus=0, ind_Siplus=0, ind_Splus=0 &
                        , ind_Nplus=0, ind_Feplus=0
  INTEGER(KIND=LONG) :: ind_G=0, ind_Gplus=0, ind_Gminus=0
  INTEGER(KIND=LONG) :: ind_e=0, ind_PHOTON=0, ind_CRP=0, ind_GRAIN=0, ind_SECPHO=0
  INTEGER(KIND=LONG) :: ind_D=0

  !-------------------------------------------
  ! mass of some species
  ! remark : CP -> C+, HP -> H+, SiP -> Si+
  !          SP -> S+, NP -> N+
  ! initialized in READ_SPECIES
  !-------------------------------------------
  REAL(KIND=DP) :: mass_H=0.0_DP, mass_H2=0.0_DP, mass_He=0.0_DP, mass_O=0.0_DP, mass_Oplus=0.0_DP
  REAL(KIND=DP) :: mass_N=0.0_DP, mass_C=0.0_DP, mass_S=0.0_DP, mass_Si=0.0_DP, mass_H2O=0.0_DP
  REAL(KIND=DP) :: mass_OH=0.0_DP, mass_CO=0.0_DP, mass_NH3=0.0_DP, mass_CH3OH=0.0_DP, mass_SiO=0.0_DP
  REAL(KIND=DP) :: mass_Cplus=0.0_DP, mass_Hplus=0.0_DP, mass_Siplus=0.0_DP
  REAL(KIND=DP) :: mass_G=0.0_DP, mass_Gplus=0.0_DP, mass_Gminus=0.0_DP
  REAL(KIND=DP) :: mass_D=0.0_DP, mass_Splus=0.0_DP, mass_Nplus=0.0_DP, mass_Feplus=0.0_DP

  !-------------------------------------------
  ! density of some species
  ! remark : CP -> C+, HP -> H+, SiP -> Si+
  !          SP -> S+, NP -> N+
  ! calculated in DIFFUN
  !-------------------------------------------
  REAL(KIND=DP) :: Dens_H, Dens_H2, Dens_He, Dens_O, Dens_Oplus
  REAL(KIND=DP) :: Dens_N, Dens_C, Dens_S, Dens_Si, Dens_H2O
  REAL(KIND=DP) :: Dens_OH, Dens_CO, Dens_NH3, Dens_CH3OH, Dens_SiO
  REAL(KIND=DP) :: Dens_Cplus, Dens_Siplus, Dens_Hplus, Dens_Splus, Dens_Nplus, Dens_Feplus
  REAL(KIND=DP) :: Dens_G, Dens_Gplus, Dens_Gminus
  REAL(KIND=DP) :: Dens_ongrains, Dens_cor, Dens_e


  ! --------------------------------------------------------------------
  ! index of beginning and end of each type of species (calculated in INITIALIZE)
  ! --------------------------------------------------------------------
  INTEGER(KIND=LONG) :: b_neu=0, e_neu=0 ! neutrals
  INTEGER(KIND=LONG) :: b_ion=0, e_ion=0 ! positive ions
  INTEGER(KIND=LONG) :: b_neg=0, e_neg=0 ! negative ions
  INTEGER(KIND=LONG) :: b_gra=0, e_gra=0 ! species on grain mantles
  INTEGER(KIND=LONG) :: b_cor=0, e_cor=0 ! species in grain cores


  ! read/write format for one species
  CHARACTER(len=*), PRIVATE, PARAMETER :: format_specy= &
       '(I3,2X,A7,2X,2I2,10I1,ES10.3,F10.3)'

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS


  SUBROUTINE INITIALIZE_ELEMENTS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    initialize the elements and their abundances (from
    !    Anders & Grevesse 1989).
    !    elements are 'basic' components of molecules : H, C, N ...
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    elements
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : AMU, Zero

    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i

    ! initialization
    elements(:)%name    = ''
    elements(:)%mass    = Zero
    elements(:)%ab      = Zero
    elements(:)%ab_ref  = Zero
    elements(:)%ab_init = Zero

    ! first with mass in amu
    i = 1 ; elements(i)%name = 'H' ; elements(i)%mass = 1.00797_DP ;elements(i)%ab_ref = 1._DP
            ind_elem_H = i ! index of H element in this vector
    i = 2 ; elements(i)%name = 'O' ; elements(i)%mass = 15.9994_DP ;elements(i)%ab_ref = 8.53D-4
    i = 3 ; elements(i)%name = 'C' ; elements(i)%mass = 12.0111_DP ;elements(i)%ab_ref = 3.62D-4
    i = 4 ; elements(i)%name = 'N' ; elements(i)%mass = 14.0067_DP ;elements(i)%ab_ref = 1.12D-4
    i = 5 ; elements(i)%name = 'He'; elements(i)%mass = 4.00260_DP ;elements(i)%ab_ref = 1.0D-1
    i = 6 ; elements(i)%name = 'Na'; elements(i)%mass = 22.9898_DP ;elements(i)%ab_ref = 2.06D-6
    i = 7 ; elements(i)%name = 'Mg'; elements(i)%mass = 24.3_DP    ;elements(i)%ab_ref = 3.85D-5
    i = 8 ; elements(i)%name = 'S' ; elements(i)%mass = 32.0640_DP ;elements(i)%ab_ref = 1.85D-5
    i = 9 ; elements(i)%name = 'Si'; elements(i)%mass = 28.0860_DP ;elements(i)%ab_ref = 3.58D-5
    i = 10; elements(i)%name = 'Fe'; elements(i)%mass = 55.8470_DP ;elements(i)%ab_ref = 3.23D-5
    i = 11 ; elements(i)%name = 'D' ; elements(i)%mass = 2.00000_DP ;elements(i)%ab_ref = 2.00D-5
    ! SC May 06: replace C60 by element G (generic grain); small mass to keep muI unaffected
    i = 12 ; elements(i)%name = 'G' ; elements(i)%mass = 1._DP ;elements(i)%ab_ref = 1.00D-10

    ! conversion of mass : AMU -> g
    elements(:)%mass = AMU * elements(:)%mass

  END SUBROUTINE INITIALIZE_ELEMENTS

  SUBROUTINE READ_SPECIES
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     read chemical species from the input file and compute indices of
    !     current species; add e-, photon, CRP, grain, SECPHO
    !     method : the file is read twice (first to count the number of species)
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     speci, mass_X, ind_X (X is a species)
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    USE MODULE_PHYS_VAR, ONLY : nH_init
    IMPLICIT NONE
    CHARACTER(len=name_length) :: charact
    INTEGER(KIND=LONG) :: i, error, ii, Nlines_comment
    CHARACTER(len=*), PARAMETER :: name_file_speci='input/species.in'
    INTEGER(KIND=LONG) :: file_speci
    REAL(KIND=DP) :: density_limit = 1.D-20

    ! open file
    file_speci = GET_FILE_NUMBER()
    OPEN(file_speci,file=name_file_speci,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! comments
    nlines_comment = 6
    DO i=1,nlines_comment
       READ(file_speci,'(A1)') charact
    END DO

    ! count the number of chemical species
    Nspec = 0 ; error = 0
    DO WHILE (error == 0) ! stop at error or end of file
       charact = ''
!!$       READ(file_speci,'(A1)',iostat=error) charact        ! new
       READ(file_speci,format_specy,iostat=error) ii, charact ! old
       SELECT CASE (charact(1:1))
       CASE ('A':'Z') ! if the character is a capital letter, this is a species
          Nspec = Nspec + 1
       CASE DEFAULT ! if not a species, nor end of file -> error
          IF (error >=0) STOP "*** WARNING1, error in READ_SPECIES"
       END SELECT
    ENDDO

    ! allocation and initialization of the vector 'species'
    Nspec_plus = Nspec + 5          ! add e-, photon, grains, CRP, SECPHO
    ALLOCATE(speci(0:Nspec_plus)) ! start at zero : undetermined species
    speci(:)%name        = ''
    speci(:)%mass        = Zero
    speci(:)%enthalpy    = Zero
    speci(:)%density     = Zero
    speci(:)%Dens_old    = Zero
    speci(:)%Col_dens    = Zero
    speci(:)%velocity    = Zero
    speci(:)%temperature = Zero
! SC: initialize formula for non-existent species
    speci(0)%formula(:)  = Zero

    ! close and re-open file
    CLOSE(file_speci)
    OPEN(file_speci,file=name_file_speci,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! go directly to the chemical species
    DO i=1, Nlines_comment
       READ(file_speci,'(A1)') charact
    END DO
    i = 0 ; error = 0                     ! stop reading after Nspec or error
    DO WHILE (i < Nspec .AND. error == 0)
       i = i + 1
       READ(file_speci,format_specy, iostat=error) &
            ii, &                                ! useless
            speci(i)%name, &                     ! name
            speci(i)%useless, &                  ! useless, replaced by CHEMICAL_FORMULA
            speci(i)%density, &                  ! n(X) / nH
            speci(i)%enthalpy                    ! kCal/mol

    ! Convert initial densities to cm-3
       speci(i)%density = speci(i)%density * nH_init

       IF (error == 0) THEN
          speci(i)%index = i
          ! number of each element in the molecule
          speci(i)%formula = CHEMICAL_FORMULA(speci(i)%name)
          ! mass is determined from the formula and the mass of each element (g)
          speci(i)%mass = DOT_PRODUCT(DBLE(elements(:)%mass),DBLE(speci(i)%formula))
          ! indices and masses of most used chemical species
          SELECT CASE (TRIM(speci(i)%name))
          CASE ('H')
             ind_H  = i
             mass_H = speci(i)%mass
          CASE ('D')
             ind_D  = i
             mass_D = speci(i)%mass
          CASE ('H2')
             ind_H2  = i
             mass_H2 = speci(i)%mass
          CASE ('He')
             ind_He  = i
             mass_He = speci(i)%mass
          CASE ('O')
             ind_O  = i
             mass_O = speci(i)%mass
          CASE ('O+')
             ind_Oplus  = i
             mass_Oplus = speci(i)%mass
          CASE ('N')
             ind_N  = i
             mass_N = speci(i)%mass
          CASE ('C')
             ind_C  = i
             mass_C = speci(i)%mass
          CASE ('S')
             ind_S  = i
             mass_S = speci(i)%mass
          CASE ('Si')
             ind_Si  = i
             mass_Si = speci(i)%mass
          CASE ('H2O')
             ind_H2O  = i
             mass_H2O = speci(i)%mass
          CASE ('OH')
             ind_OH  = i
             mass_OH = speci(i)%mass
          CASE ('CO')
             ind_CO  = i
             mass_CO = speci(i)%mass
          CASE ('NH3')
             ind_NH3  = i
             mass_NH3 = speci(i)%mass
          CASE ('CH4O')
             ind_CH3OH  = i
             mass_CH3OH = speci(i)%mass
          CASE ('SiO')
             ind_SiO  = i
             mass_SiO = speci(i)%mass
          CASE ('G')
             ind_G  = i
             mass_G = speci(i)%mass
          CASE ('C+')
             ind_Cplus  = i
             mass_Cplus = speci(i)%mass
          CASE ('H+')
             ind_Hplus  = i
             mass_Hplus = speci(i)%mass
          CASE ('Si+')
             ind_Siplus  = i
             mass_Siplus = speci(i)%mass
          CASE ('S+')
             ind_Splus  = i
             mass_Splus = speci(i)%mass
          CASE ('N+')
             ind_Nplus  = i
             mass_Nplus = speci(i)%mass
          CASE ('Fe+')
             ind_Feplus  = i
             mass_Feplus = speci(i)%mass
          CASE ('G+')
             ind_Gplus  = i
             mass_Gplus = speci(i)%mass
          CASE ('G-')
             ind_Gminus  = i
             mass_Gminus = speci(i)%mass
          END SELECT
       ELSE
          IF (error > 0) STOP "*** WARNING2, error in READ_SPECIES"
       END IF
    END DO

    ! file closure
    CLOSE(file_speci)

    !-------------------------------------------------------
    ! avoid too small numbers : lower limit to the densities
    !-------------------------------------------------------
    WHERE (speci(1:Nspec)%density < density_limit)
       speci(1:Nspec)%density = density_limit
    END WHERE

    !--------------------------------------------------------------------------
    ! addition of e-, photon, CRP, GRAIN, SECPHO with : mass and enthalpy = 0.0
    !--------------------------------------------------------------------------
    i = Nspec
    ! e-, initial density is set in INITIALIZE
    i = i + 1
    speci(i)%name  = 'ELECTR'
    speci(i)%index = i
    ind_e = i
    ! GRAIN, initial density is set in INITIALIZE
    i = i + 1
    speci(i)%name  = 'GRAIN'
    speci(i)%index = i
    ind_GRAIN = i
    ! photon, density = constant = 1._DP (cf. CHEMISTRY)
    i = i + 1
    speci(i)%name  = 'PHOTON'
    speci(i)%index = i
    speci(i)%density = 1._DP
    ind_PHOTON = i
    ! CRP (cosmic ray proton), density = constant = 1._DP (cf. CHEMISTRY)
    i = i + 1
    speci(i)%name  = 'CRP'
    speci(i)%index = i
    speci(i)%density = 1._DP
    ind_CRP = i
    ! SECPHO (secondary photon), initial = constant = 1._DP (cf. CHEMISTRY)
    i = i + 1
    speci(i)%name  = 'SECPHO'
    speci(i)%index = i
    speci(i)%density = 1._DP
    ind_SECPHO = i

  END SUBROUTINE READ_SPECIES

  SUBROUTINE CHECK_SPECIES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    tests if the set of species contains identical species
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i, j

    DO i=1, Nspec
       DO j=i+1,Nspec
          IF (speci(i)%name==speci(j)%name) THEN
             WRITE(*,*) "*** WARNING, species ",i," and ",j," are identical"
             STOP
          END IF
       END DO
    END DO
  END SUBROUTINE CHECK_SPECIES

  SUBROUTINE WRITE_SPECY(num,one_specy)
    !--------------------------------------------------------------------------
    ! called by :
    !     WRITE_INFO
    ! purpose :
    !     write information about one chemical species
    ! subroutine/function needed :
    ! input variables :
    !     * num -> file number to which to write the information
    !     * one_specy -> (type TYPE_SPECY) chemical species
    ! output variables :
    ! results :
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(type_specy), INTENT(in) :: one_specy
    INTEGER(KIND=LONG),INTENT(in) :: num
    INTEGER(KIND=LONG) :: i

    WRITE(num,format_specy) &
         one_specy%index, &
         one_specy%name, &
         one_specy%useless, &
         one_specy%density, &
         one_specy%enthalpy
  END SUBROUTINE WRITE_SPECY


  SUBROUTINE ELEMENTAL_ABUNDANCES
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    computes elemental abundances in (H,C,N,...)
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    elements%ab
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i

    ! computation of abundances, using speci(i)%formula and speci(i)%density
    DO i=1, Nelements
       elements(i)%ab = SUM(DBLE(speci(:)%formula(i)) * speci(:)%density)
    END DO
    ! normalization to H abundance
    elements(:)%ab = elements(:)%ab / elements(ind_elem_H)%ab

!!$    ! check the calculated and the read masses
!!$    DO i=1,Nspec
!!$       speci(i)%mass2=DOT_PRODUCT(DBLE(elements(:)%mass),DBLE(speci(i)%formula))
!!$       IF (abs((speci(i)%mass-speci(i)%mass2)/speci(i)%mass) > 1.D-5) &
!!$            write(*,*)i,speci(i)%name,speci(i)%mass,speci(i)%mass2
!!$    END DO

  END SUBROUTINE ELEMENTAL_ABUNDANCES


  FUNCTION CHEMICAL_FORMULA (name) RESULT (formula)
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_SPECIES
    ! purpose :
    !    computes for one molecule its composition in elements H, C, ...
    ! subroutine/function needed :
    ! input variables :
    !    name -> molecule name
    ! ouput variables :
    ! results :
    !    formula -> vector containing the number of each element in the molecule
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(len=name_length), INTENT(in) :: name
    INTEGER(KIND=LONG), DIMENSION(Nelements) :: formula
    INTEGER(KIND=LONG) :: i,j,length
    INTEGER(KIND=LONG) :: digit1,digit2,number_element
    INTEGER(KIND=LONG) :: howmany,howmany_element,howmany_digit
    CHARACTER(len=8), DIMENSION(name_length) :: table_name
    CHARACTER(len=2) :: charact

    ! initialization
    table_name(:) = ''      ! describes what is in name(i)
    length = LEN_TRIM(name) ! length of name, without blanks
    formula = 0

!!$    NEW VERSION WITH "real" names : He, not HE ...
    ! filling the vector table_name
    DO i=1,length
       SELECT CASE (name(i:i))
       CASE ('0':'9')! digit
          table_name(i)="digit"
       CASE ('a':'z')! 'He', 'Si' ...
          table_name(i)="2letters"
       CASE ('+','-','*') ! ion or on grain
          table_name(i)="iongrain"
       CASE ('A':'Z') ! capital letters
          table_name(i)="element"
       CASE DEFAULT ! other
          STOP "*** WARNING, incorrect character in CHEMICAL_FORMULA"
       END SELECT
    END DO

!!$    OLD VERSION WITH 'bad' names : HE and not He
!!$    ! filling the vector table_name
!!$    DO i=1,length
!!$       SELECT CASE (name(i:i))
!!$       CASE ('0':'9')! digit
!!$          table_name(i) = "digit"
!!$       CASE ('E','A','I', 'G')! 'HE', 'FE', 'NA', 'SI', 'MG' ...
!!$          table_name(i) = "2letters"
!!$       CASE ('+','-','*') ! ion or on grain
!!$          table_name(i) = "iongrain"
!!$       CASE DEFAULT ! element
!!$          table_name(i) = "element"
!!$       END SELECT
!!$    END DO

    ! chemical formula
    i = 1
    DO WHILE (i <= length)
       digit1 = 1 ; digit2 = 0
       howmany_digit   = 1
       howmany_element = 1 ! 1 by default (if no digit)
       SELECT CASE (TRIM(table_name(i)))
       CASE ("element")
          ! determine how many characters to take into account -> howmany
          IF (i < length) THEN
             IF (table_name(i+1) == "2letters") howmany_element = 2
          END IF
          howmany = howmany_element
          ! extracts from name
          charact = name(i:i+howmany-1)
          ! research of the element
          DO j=1,Nelements
             IF (charact==TRIM(elements(j)%name)) EXIT
          END DO
          ! if not found ...
          IF ((j == Nelements) .AND. (charact /= TRIM(elements(j)%name))) THEN
             STOP "*** WARNING, no chemical element for this molecule"
          END IF
          number_element = j

       CASE ("digit")
          ! determine how many characters to take into acount -> howmany
          IF (i < length) THEN
             IF (table_name(i+1) == "digit") howmany_digit = 2
          END IF
          howmany = howmany_digit
          ! research of the number = digit1 or digit1*10+digit2
          digit1 = IACHAR(name(i:i)) - IACHAR('0')
          IF (howmany == 2) digit2 = IACHAR(name(i+1:i+1))-IACHAR('0')
       END SELECT

       ! determination of the number of elements in the molecule
       SELECT CASE (table_name(i))
       CASE('element','digit')
          formula(number_element) = digit1
          IF (howmany_digit == 2) formula(number_element) = digit1 * 10 + digit2
       CASE('iongrain') ! do nothing
       END SELECT

       ! look at the next characters in the name
       i = i + howmany
    END DO

  END FUNCTION CHEMICAL_FORMULA


END MODULE MODULE_CHEMICAL_SPECIES
