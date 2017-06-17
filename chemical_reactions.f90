MODULE MODULE_CHEM_REACT
  !*****************************************************************************
  !** The module 'MODULE_CHEM_REACT' contains variables and                   **
  !** subroutines related to chemical reactions but not the subroutine        **
  !** CHEMISTRY (see MODULE_EVOLUTION)                                        **
  !** REMARK :                                                                **
  !** --------                                                                **
  !**   the entire module 'MODULE_CHEMICAL_SPECIES' is included here.         **
  !*****************************************************************************
  USE MODULE_CHEMICAL_SPECIES
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !----------------------------
  ! data type of one reaction
  !----------------------------
  TYPE TYPE_REACTION
     INTEGER(KIND=LONG),DIMENSION(2) :: R                  ! index of 2 reactants
     INTEGER(KIND=LONG),DIMENSION(4) :: P                  ! index of 4 products
     REAL(KIND=DP)                   :: gamma, alpha, beta ! Arrhenius coefficients
     REAL(KIND=DP)                   :: DE                 ! exo(if >0)/endo(if <0)-thermicity
     REAL(KIND=DP)                   :: mass_prod          ! sum of the products' masses (g)
     INTEGER(KIND=LONG)              :: Nprod_m1           ! number of products minus one
     CHARACTER(LEN=5)                :: ref                ! reference of the reaction
     CHARACTER(LEN=5)                :: TYPE               ! type of the reaction
  END TYPE TYPE_REACTION

  !-----------------------------------------------------------------
  ! number of reactions in the model and in each type of reaction
  ! calculated in READ_REACTIONS
  ! Nreact is modified in ADD_INVERSE_REACTIONS
  ! as we add the missing endothermic reactions (if the boolean
  ! do_we_add_reactions is .TRUE.)
  !-----------------------------------------------------------------
  INTEGER(KIND=LONG),PARAMETER :: Nreact_max=2000 ! max. number of reactions
  INTEGER(KIND=LONG)           :: Nreact          ! number of reactions
  ! LOGICAL, PARAMETER :: do_we_add_reactions=.FALSE.   ! add the reverse reactions ?
  LOGICAL, PARAMETER :: do_we_add_reactions=.TRUE.   ! add the reverse reactions ?

  ! number of reactions in each type of reaction (PHOTO, SPUTTER, ...)
  ! calculated in INITIALIZE, except Nrever_new and Nrever : in ADD_INVERSE_REACTIONS
  INTEGER(KIND=LONG) :: Nphoto=0, Ncr_io=0, Ncr_de=0, Nh2_fo=0, Nthree=0
  INTEGER(KIND=LONG) :: Nsputt=0, Nerosi=0, Nadsor=0, Nother=0, Nrever=0
  INTEGER(KIND=LONG) :: Ndisso=0

  ! index of beginning and end of each type of reaction
  ! calculated in INITIALIZE, except b_rever and e_rever : in ADD_INVERSE_REACTIONS
  INTEGER(KIND=LONG) :: b_photo=0, e_photo=0
  INTEGER(KIND=LONG) :: b_cr_io=0, e_cr_io=0
  INTEGER(KIND=LONG) :: b_cr_de=0, e_cr_de=0
  INTEGER(KIND=LONG) :: b_h2_fo=0, e_h2_fo=0
  INTEGER(KIND=LONG) :: b_three=0, e_three=0
  INTEGER(KIND=LONG) :: b_sputt=0, e_sputt=0
  INTEGER(KIND=LONG) :: b_erosi=0, e_erosi=0
  INTEGER(KIND=LONG) :: b_adsor=0, e_adsor=0
  INTEGER(KIND=LONG) :: b_disso=0, e_disso=0
  INTEGER(KIND=LONG) :: b_other=0, e_other=0
  INTEGER(KIND=LONG) :: b_rever=0, e_rever=0


  !---------------------------------------------
  ! vector containing the chemical reactions
  ! stored between the index 1 and Nreact
  !---------------------------------------------
  TYPE(TYPE_REACTION),DIMENSION(Nreact_max) :: react

  ! useful for determining which reactions are correct -> READ_REACTIONS
  CHARACTER(LEN=name_length),PARAMETER   :: charact_blank=''
  INTEGER(KIND=LONG),PARAMETER :: RP_undefined = -1 ! if reactant or product undefined -> reaction is incorrect
  INTEGER(KIND=LONG),PARAMETER :: P_none = 0        ! if product is not present in chemical species

  ! read/write format for chemical reactions
  CHARACTER(len=*), PRIVATE, PARAMETER  :: format_reaction = &
       '(A5,1X,5(A7,1X),A7,ES8.2,1X,F5.2,1X,F8.1,1X,F5.2)'

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  SUBROUTINE WRITE_REACTION(num,one_reaction)
    !--------------------------------------------------------------------------
    ! called by :
    !     * WRITE_INFO
    !     * CHECK_REACTIONS
    ! purpose :
    !     write informations about one chemical reaction
    ! subroutine/function needed :
    ! input variables :
    !     * num -> file number to which to write the information
    !     * one_reaction -> (type TYPE_REACTION) chemical reaction
    ! ouput variables :
    ! results :
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG),INTENT(in)::num
    TYPE(TYPE_REACTION),INTENT(in)::one_reaction

    WRITE(num,format_reaction)&
         one_reaction%ref, &
         speci(one_reaction%R(1:2))%name, &
         speci(one_reaction%P(1:4))%name, &
         one_reaction%gamma, &
         one_reaction%alpha, &
         one_reaction%beta, &
         one_reaction%DE
  END SUBROUTINE WRITE_REACTION

  FUNCTION REACTIONS_IDENTICAL(reaction1, reaction2) RESULT(res)
    !---------------------------------------------------------------------------
    ! called by :
    !     * CHECK_REACTIONS
    !     * ADD_INVERSE_REACTIONS
    ! purpose :
    !     test if 2 reactions are identical (same reactants and products)
    !     whatever the order in which reactants and products are written
    ! subroutine/function needed :
    ! input variables :
    !     reaction1, reaction2 -> (type TYPE_REACTION) the 2 reactions
    ! output variables :
    ! results :
    !     (boolean) .TRUE. if the reactions are identical
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(TYPE_REACTION), INTENT(in) :: reaction1 , reaction2
    LOGICAL :: reactants_identical, produits_identical, res
    INTEGER(KIND=LONG) :: k,l
    LOGICAL,DIMENSION(4) :: prod_id1, prod_id2

    ! test if reactants are identical
    reactants_identical = &
         ((reaction1%R(1) == reaction2%R(1) .AND. &
         reaction1%R(2) == reaction2%R(2)) .OR. &
         (reaction1%R(1) == reaction2%R(2) .AND. &
         reaction1%R(2) == reaction2%R(1)))

    ! if reactants identical, test if products are identical
    produits_identical = .FALSE.
    IF (reactants_identical) THEN
       ! test if each product of reaction 1 is in reaction 2
       prod_id1(:) = .FALSE.
       DO k=1,4
          DO l=1,4
             IF (reaction1%P(k) == reaction2%P(l)) prod_id1(k) = .TRUE.
          END DO
       ENDDO
       ! test if each product of reaction 2 is in reaction 1
       prod_id2(:) = .FALSE.
       DO k=1,4
          DO l=1,4
             IF (reaction2%P(k)==reaction1%P(l)) prod_id2(k)=.TRUE.
          END DO
       ENDDO
       ! produits_identical is .TRUE. if all products are common to the 2 reactions
       produits_identical = ALL(prod_id1) .AND. ALL(prod_id2)
    END IF
    ! result is : (same reactants) and (sames products)
    res = reactants_identical .AND. produits_identical

  END FUNCTION REACTIONS_IDENTICAL

  SUBROUTINE READ_REACTIONS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read chemical reactions in the file file_chemistry
    !    and find index of each reactant and each product
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    react, Nreact
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE

    CHARACTER(len=5),PARAMETER :: e_comment='!end', e_file='END'
    CHARACTER(len=5) :: ref=''
    CHARACTER(len=8) :: R1, R2, P1, P2, P3, P4 ! names for reactants and products
    REAL(KIND=DP) :: alpha, beta, gamma, DE ! reaction coefficients
    INTEGER(KIND=LONG) :: i, error

    CHARACTER(LEN=*), PARAMETER :: name_file_chemistry='input/chemistry.in'
    INTEGER                     :: file_chemistry

    ! initialization
    react(:)%type = ''
    react(:)%ref = ''
    DO i=1,2
       react(:)%R(i)=RP_undefined
    END DO
    DO i=1,4
       react(:)%P(i) = RP_undefined
    END DO
    react(:)%gamma = Zero
    react(:)%alpha = Zero
    react(:)%beta = Zero
    react(:)%DE = Zero
    react(:)%Nprod_m1 = -1 ! don't forget the 'minus one' !
    react(:)%mass_prod = Zero

    ! file opening
    file_chemistry = GET_FILE_NUMBER()
    OPEN(file_chemistry,file=name_file_chemistry,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! advance until the end of the comments
    ref = '' ; error = 0
    DO WHILE (ref /= e_comment .AND. error==0)
       READ(file_chemistry,'(A5)',iostat=error) ref
    END DO
    IF (error /=0) STOP "*** WARNING, error in READ_REACTIONS"

    ! stop at end of file or at error
    Nreact = 0 ; ref = '' ; error = 0
    DO WHILE ((ref /= e_file) .AND. (Nreact < Nreact_max) .AND. (error==0))
       Nreact = Nreact + 1

       R1 = ''; R2 = ''; P1 = ''; P2 = ''; P3 = ''; P4 = ''
       READ(file_chemistry,format_reaction,iostat=error) &
            ref, &
            R1, R2, P1, P2, P3, P4, &      ! names and not indices
            gamma, alpha, beta             ! Arrhenius coefficients
!           DE                             ! DE is not read, but calculated after

       IF (ref /= e_file .AND. (error==0)) THEN
          ! search for indices of reactants and products
          DO i=0, Nspec_plus ! and not 1,Nspec, as we search also none, e-, photon, crp, grain
             IF (R1 == speci(i)%name) react(Nreact)%R(1) = i
             IF (R2 == speci(i)%name) react(Nreact)%R(2) = i
             IF (P1 == speci(i)%name) react(Nreact)%P(1) = i
             IF (P2 == speci(i)%name) react(Nreact)%P(2) = i
             IF (P3 == speci(i)%name) react(Nreact)%P(3) = i
             IF (P4 == speci(i)%name) react(Nreact)%P(4) = i
          END DO

          ! fill the fields of the reaction
          react(Nreact)%ref = ref
          react(Nreact)%gamma = gamma
          react(Nreact)%alpha = alpha
          react(Nreact)%beta = beta
          !react(Nreact)%DE = DE      ! DE is not read, but calculated after

          ! calculates the number of products and the sum of their masses (used in CHEMISTRY)
          IF (P1 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(1))%mass
          ENDIF
          IF (P2 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(2))%mass
          ENDIF
          IF (P3 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(3))%mass
          ENDIF
          IF (P4 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(4))%mass
          ENDIF

       END IF
    END DO
    ! if end of file has been reached, we counted 1 reaction too many
    IF (ref == e_file) Nreact = Nreact-1

    ! file closure
    CLOSE(file_chemistry)

  END SUBROUTINE READ_REACTIONS



  SUBROUTINE ENERGY_DEFECT_REACTION(i)
    !---------------------------------------------------------------------------
    ! called by :
    !    REACTION_TYPE
    ! purpose :
    !    Calculate the energy defect (DE) for the reaction i.
    !    If this reaction has one species with and unknown enthalpy, set DE = Zero
    !    except if this species appears (with the same occurrence) in BOTH the
    !    reactants and products. In this case, the uncertainty cancels and one
    !    can calculate DE = enthalpy(reactants)-enthalpy(products).
    !    REMARK :
    !       energy defect is set to zero for endothermic reactions
    !       (DE < Zero) => DE = Zero
    ! subroutine/function needed :
    ! input variables :
    !    i -> index of the chemical reaction
    ! ouput variables :
    ! results :
    !     react
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : kCaleV, Zero

    IMPLICIT NONE
    INTEGER(KIND=LONG), INTENT(in) :: i
    REAL(KIND=DP), PARAMETER       :: enthalpy_threshold=-99.9_DP
    LOGICAL                        :: enthalpy_unknown
    INTEGER(KIND=LONG), DIMENSION(0:Nspec_plus) :: occ_reactants
    INTEGER(KIND=LONG), DIMENSION(0:Nspec_plus) :: occ_products
    INTEGER :: j

    ! --- check if enthalpy of one species is unknown         ---
    ! --- and if this species appears with the same occurrence ---
    ! --- in reactants and in products                      ---
    ! initialization
    occ_reactants = 0 ! counts how many times one species appears in the reactants
    occ_products = 0  ! counts how many times one species appears in the products
    DO j=1,2          ! count the occurrence of the reactants
       occ_reactants(react(i)%R(j)) = occ_reactants(react(i)%R(j))+1
    END DO
    DO j=1,4          ! count the occurrence of the products
       occ_products(react(i)%P(j)) = occ_products(react(i)%P(j))+1
    END DO
    enthalpy_unknown = .FALSE.
    DO j=1,2          ! check for enthalpy of reactants
       IF (speci(react(i)%R(j))%enthalpy < enthalpy_threshold .AND. &
            occ_reactants(react(i)%R(j)) /= occ_products(react(i)%R(j))) &
            enthalpy_unknown = .TRUE.
    END DO
    DO j=1,4          ! check for enthalpy of products
       IF (speci(react(i)%P(j))%enthalpy < enthalpy_threshold .AND. &
            occ_reactants(react(i)%P(j)) /= occ_products(react(i)%P(j))) &
            enthalpy_unknown = .TRUE.
    END DO

    ! --- calculation of DE (eV, whereas speci%enthaply is in kCal/mol) ---
    IF (enthalpy_unknown) THEN
       ! DE = 0.0 in this case
       react(i)%DE = Zero
    ELSE
       ! DE = enthalpy (reactants) - enthalpy (products)
       react(i)%DE = kCaleV * ( &
            SUM(DBLE(speci(react(i)%R(:))%enthalpy)) - &
            SUM(DBLE(speci(react(i)%P(:))%enthalpy)))
    END IF

    ! set DE=Zero for endothermic reactions
    IF (react(i)%DE < Zero) react(i)%DE = Zero

  END SUBROUTINE ENERGY_DEFECT_REACTION


  SUBROUTINE REACTION_TYPE
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    to find the type of each reaction and re-order the entire set
    !    according to the different types. Compute the
    !    indices (beginning, end, number) of each reaction type.
    !    Also calculate DE for each reaction, according to its type.
    !    Order is :
    !       (1) : photo-reactions (type='PHOTO')
    !       (2) : cosmic ray ionization or dissociation (type='CR_IO')
    !       (3) : cosmic ray induced desorption from grains (type='CR_DE')
    !       (4) : H2 and HD formation (type='H2_FO')
    !       (5) : three-body reactions on grain surfaces (type='THREE')
    !       (6) : sputtering of grain mantles (type='SPUTT')
    !       (7) : erosion of grain cores (type='EROSI')
    !       (8) : adsorption on to grain surface (type='ADSOR')
    !       (9) : collisional dissociation of H2 (type='DISSO') - see below
    !      (10) : all other reactions (type='OTHER')
    !      (11) : reverse endothermic reactions (type='REVER')
    !             these reactions are defined and added in ADD_REVERSE_REACTIONS
    ! subroutine/function needed :
    !    ENERGY_DEFECT_REACTION
    ! input variables :
    ! ouput variables :
    ! results :
    !     reactions
    !
    !  22 juin 2001 - JLB - Add collisional dissociation of H2
    !       WARNING ! Required order in chemistry file : H2 + X -> X + H + H
    !                 identification is done on R1=H2, P2=H, P3=H, R2 = P1
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE

    INTEGER(KIND=LONG) :: i
    TYPE (TYPE_REACTION), DIMENSION(:), ALLOCATABLE :: react_aux

    !-------------------------------------
    ! Save the unsorted reaction set
    ! react_aux = old set (unsorted)
    !-------------------------------------
    ALLOCATE(react_aux(1:Nreact))
    react_aux(1:Nreact) = react(1:Nreact)

    !--- initialize the vector 'react' ---
    react(:)%type = ''
    react(:)%ref = ''
    DO i=1,2
       react(:)%R(i) = RP_undefined
    END DO
    DO i=1,4
       react(:)%P(i) = RP_undefined
    END DO
    react(:)%gamma = Zero
    react(:)%alpha = Zero
    react(:)%beta = Zero
    react(:)%DE = Zero
    react(:)%Nprod_m1 = -1 ! useless here ...
    react(:)%mass_prod = Zero

    !----------------------------------------------------------
    ! Find the type of each reaction.
    ! Re-order the entire set and calculate the energy defect
    ! according to this type
    ! NOW : 'react' = sorted reactions set
    !----------------------------------------------------------
    !--- photo-reactions (type='PHOTO')---
    b_photo = 1 ! index of beginning
    DO i=1, Nreact
       IF (react_aux(i)%R(2)==ind_PHOTON) THEN
          Nphoto = Nphoto + 1
          react(Nphoto+ b_photo-1) = react_aux(i)
          react(Nphoto+ b_photo-1)%DE = Zero
          react_aux(i)%type = 'PHOTO'
          react(Nphoto+ b_photo-1)%type = 'PHOTO'
       ENDIF
    END DO
    e_photo = b_photo + Nphoto - 1 ! index of end

    !--- cosmic ray ionization or dissociation (type='CR_IO') ---
    b_cr_io = e_photo + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            (react_aux(i)%R(2) == ind_CRP .OR. react_aux(i)%R(2) == ind_SECPHO) .AND. &
            react_aux(i)%P(2) /= ind_GRAIN) THEN
          Ncr_io = Ncr_io + 1
          react(Ncr_io+ b_cr_io-1) = react_aux(i)
          react(Ncr_io+ b_cr_io-1)%DE = Zero
          react_aux(i)%type = 'CR_IO'
          react(Ncr_io+ b_cr_io-1)%type = 'CR_IO'
          ! note that react(i)%beta is preset at a very large value
          ! for direct cosmic ray ionization : react(i)%beta = 1.0D8
          IF (react(Ncr_io+ b_cr_io-1)%beta < 1.D-9) &
               react(Ncr_io+ b_cr_io-1)%beta = 1.D8
       ENDIF
    END DO
    e_cr_io = b_cr_io + Ncr_io - 1 ! index of end

    !--- cosmic ray induced desorption from grains (type='CR_DE') ---
    b_cr_de = e_cr_io + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            react_aux(i)%R(2) == ind_CRP .AND. &
            react_aux(i)%P(2) == ind_GRAIN) THEN
          Ncr_de = Ncr_de + 1
          react(Ncr_de+ b_cr_de-1) = react_aux(i)
          react(Ncr_de+ b_cr_de-1)%DE = Zero
          react_aux(i)%type = 'CR_DE'
          react(Ncr_de+ b_cr_de-1)%type = 'CR_DE'
       ENDIF
    END DO
    e_cr_de = b_cr_de + Ncr_de - 1 ! index of end

    !--- H2 and HD formation (type='H2_FO') ---
    b_h2_fo = e_cr_de + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            react_aux(i)%R(1) == ind_H .AND. &
            react_aux(i)%R(2) == ind_H) THEN
          Nh2_fo = Nh2_fo + 1
          react(Nh2_fo+ b_h2_fo-1) = react_aux(i)
          CALL ENERGY_DEFECT_REACTION(Nh2_fo+ b_h2_fo-1) ! explicitly calculated
          react_aux(i)%type = 'H2_FO'
          react(Nh2_fo+ b_h2_fo-1)%type = 'H2_FO'
       ENDIF
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            react_aux(i)%R(1) == ind_H .AND. &
            react_aux(i)%R(2) == ind_D) THEN
          Nh2_fo = Nh2_fo + 1
          react(Nh2_fo+ b_h2_fo-1) = react_aux(i)
          CALL ENERGY_DEFECT_REACTION(Nh2_fo+ b_h2_fo-1) ! explicitly calculated
          react_aux(i)%type = 'H2_FO'
          react(Nh2_fo+ b_h2_fo-1)%type = 'H2_FO'
       ENDIF
    END DO
    e_h2_fo = b_h2_fo + Nh2_fo - 1 ! index of end

    !--- three-body reactions on grain surfaces (type='THREE') ---
    b_three = e_h2_fo + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            react_aux(i)%R(1) == ind_H .AND. &
            react_aux(i)%P(2) == ind_GRAIN) THEN
          Nthree = Nthree + 1
          react(Nthree+ b_three-1) = react_aux(i)
          react(Nthree+ b_three-1)%DE = Zero
          react_aux(i)%type = 'THREE'
          react(Nthree+ b_three-1)%type = 'THREE'
       ENDIF
    END DO
    e_three = b_three + Nthree - 1 ! index of end

    !--- sputtering of grain mantles (type='SPUTT') ---
    b_sputt = e_three + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            react_aux(i)%P(3) == ind_GRAIN .OR. &
            react_aux(i)%P(4) == ind_GRAIN) THEN
          Nsputt = Nsputt + 1
          react(Nsputt+ b_sputt-1) = react_aux(i)
          react(Nsputt+ b_sputt-1)%DE = Zero
          react_aux(i)%type = 'SPUTT'
          react(Nsputt+ b_sputt-1)%type = 'SPUTT'
       ENDIF
    END DO
    e_sputt = b_sputt + Nsputt - 1 ! index of end

    !--- erosion of grain cores (type='EROSI') ---
    b_erosi =  e_sputt + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            react_aux(i)%P(1) == ind_GRAIN) THEN
          Nerosi = Nerosi + 1
          react(Nerosi+ b_erosi-1) = react_aux(i)
          react(Nerosi+ b_erosi-1)%DE = Zero
          react_aux(i)%type = 'EROSI'
          react(Nerosi+ b_erosi-1)%type = 'EROSI'
       ENDIF
    END DO
    e_erosi = b_erosi + Nerosi - 1 ! index of end

    !--- adsorption on to a grain surface ---
    b_adsor = e_erosi + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
            react_aux(i)%R(2) == ind_GRAIN) THEN
          Nadsor = Nadsor + 1
          react(Nadsor+ b_adsor-1) = react_aux(i)
          react(Nadsor+ b_adsor-1)%DE = Zero
          react_aux(i)%type = 'ADSOR'
          react(Nadsor+ b_adsor-1)%type = 'ADSOR'
       ENDIF
    END DO
    e_adsor = b_adsor + Nadsor - 1 ! index of end

    !--- collisional dissociation of H2 ---
    b_disso = e_adsor + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '' .AND. &
             react_aux(i)%R(1) == ind_H2 .AND. &
             react_aux(i)%P(2) == ind_H .AND. &
             react_aux(i)%P(3) == ind_H .AND. &
             react_aux(i)%P(1) == react_aux(i)%R(2) ) THEN
          Ndisso = Ndisso + 1
          react(Ndisso+ b_disso-1) = react_aux(i)
          react(Ndisso+ b_disso-1)%DE = Zero
          react_aux(i)%type = 'DISSO'
          react(Ndisso+ b_disso-1)%type = 'DISSO'
       ENDIF
    END DO
    e_disso = b_disso + Ndisso - 1 ! index of end

    !--- other-reactions ---
    b_other = e_disso + 1 ! index of beginning
    DO i=1, Nreact
       IF (TRIM(react_aux(i)%type) == '') THEN
          Nother = Nother + 1
          react(Nother+ b_other-1) = react_aux(i)
          CALL ENERGY_DEFECT_REACTION(Nother+ b_other-1) ! explicitly calculated
          react(Nother+ b_other-1)%type = 'OTHER'
       ENDIF
    END DO
    e_other = b_other + Nother - 1 ! index of end

    !--- Nrever, b_rever, e_rever -> see ADD_REVERSE_REACTIONS ---

    !--- deallocation of the temporary variable ---
    DEALLOCATE(react_aux)

  END SUBROUTINE REACTION_TYPE


  SUBROUTINE CHECK_REACTIONS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    tests if the set of reactions is correct
    !        * each reaction must have at least 2 reactants and 1 product
    !        * reactants and products must be in the set of chemical species
    !        * species have to be conserved (except for reactions of the types
    !           THREE, EROSI or ADSOR)
    !        * charge has to be conserved
    ! subroutine/function needed :
    !    WRITE_REACTION
    !    REACTIONS_IDENTICAL
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : screen

    IMPLICIT NONE
    INTEGER(KIND=LONG)::i,j,ind
    INTEGER(KIND=LONG),DIMENSION(2,Nelements):: formula_reactants
    INTEGER(KIND=LONG),DIMENSION(4,Nelements):: formula_products
    INTEGER(KIND=LONG) :: Nions_reactants, Nneg_reactants, Ne_reactants, charge_reactants
    INTEGER(KIND=LONG) :: Nions_products, Nneg_products, Ne_products, charge_products
    CHARACTER(LEN=8) :: name
    LOGICAL :: incorrect

    DO i=1,Nreact
       ! --- checks if reaction is correct ---
       incorrect = &
            ! at least 2 reactants
       react(i)%R(1) == P_none .OR. &
            react(i)%R(2) == P_none .OR. &
            ! at least 1 product
       react(i)%Nprod_m1 < 0 .OR. &
            ! use only species in the current set
       react(i)%R(1) == RP_undefined .OR. &
            react(i)%R(2) == RP_undefined .OR. &
            react(i)%P(1) == RP_undefined .OR. &
            react(i)%P(2) == RP_undefined .OR. &
            react(i)%P(3) == RP_undefined .OR. &
            react(i)%P(4) == RP_undefined
       IF (incorrect) THEN
          WRITE(*,*) "*** WARNING, reaction ",i," is incorrect"
          CALL WRITE_REACTION(screen,react(i))
          STOP
       END IF

       ! --- species conservation ---
       formula_reactants(:,:) = 0
       formula_products(:,:)  = 0
       SELECT CASE (react(i)%type)
       CASE ('THREE', 'EROSI', 'ADSOR')
          ! for these reaction, one cannot check the conservation of species
       CASE DEFAULT
          DO j=1,2 ! chemical formula of reactants
             name=speci(react(i)%R(j))%name
             SELECT CASE (name)
             CASE ('PHOTON', 'CRP', 'ELECTR', 'GRAIN')
                ! do nothing
             CASE DEFAULT
                formula_reactants(j,:) = speci(react(i)%R(j))%formula
             END SELECT
          END DO
          DO j=1,4 ! chemical formula of products
             name=speci(react(i)%P(j))%name
             SELECT CASE (name)
             CASE ('PHOTON', 'CRP', 'ELECTR', 'GRAIN')
                ! do nothing
             CASE DEFAULT
                formula_products(j,:) = speci(react(i)%P(j))%formula
             END SELECT
          END DO
          ! checks the conservation of species between reactants and products
          DO j=1,Nelements
             IF (SUM(formula_reactants(:,j)) /= SUM(formula_products(:,j))) THEN
                WRITE(*,*)"*** WARNING, element ",elements(j)%name," is not conserved in the reaction ",i
                CALL WRITE_REACTION(screen,react(i))
                STOP
             END IF
          ENDDO
       END SELECT

       ! --- conservation of the charge ---
       Nions_reactants=0 ; Nneg_reactants=0 ; Ne_reactants=0 ; charge_reactants=0
       Nions_products=0 ; Nneg_products=0 ; Ne_products=0 ; charge_products=0
       ! counts the number of ions and electrons in the reactants
       DO j=1,2
          ind=react(i)%R(j)
          IF (ind == ind_e) THEN ! electron
             Ne_reactants = Ne_reactants + 1
          ELSE
             IF (ind >= b_ion .AND. ind <= e_ion) THEN ! ion
                Nions_reactants = Nions_reactants + 1
             ELSE
                IF (ind >= b_neg .AND. ind <= e_neg) THEN ! negative ion
                   Nneg_reactants = Nneg_reactants + 1
                END IF
             END IF
          END IF
       END DO
       ! counts the number of ions and electrons in the products
       DO j=1,4
          ind = react(i)%P(j)
          IF (ind == ind_e) THEN ! electron
             Ne_products = Ne_products + 1
          ELSE
             IF (ind >= b_ion .AND. ind <= e_ion) THEN !ion
                Nions_products = Nions_products + 1
             ELSE
                IF (ind >= b_neg .AND. ind <= e_neg) THEN ! negative ion
                   Nneg_products = Nneg_products + 1
                END IF
             END IF
          END IF
       END DO
       ! charge must be identical in reactants and products
       charge_reactants=Nions_reactants-Ne_reactants-Nneg_reactants
       charge_products=Nions_products-Ne_products-Nneg_products
       IF (charge_reactants /= charge_products) THEN
          WRITE(*,*)"*** WARNING, charge is not conserved in the reaction ",i
          CALL WRITE_REACTION(screen,react(i))
          STOP
       ENDIF

       ! --- look for the same reaction in the rest of the chemical set ---
       DO j=i+1,Nreact
          IF (REACTIONS_IDENTICAL(react(i),react(j))) THEN
             WRITE(*,*)"*** WARNING, reactions ",i," and ",j," are identical"
             CALL WRITE_REACTION(screen,react(i))
             CALL WRITE_REACTION(screen,react(j))
             STOP
          END IF
       END DO
    END DO
  END SUBROUTINE CHECK_REACTIONS



  SUBROUTINE ADD_REVERSE_REACTIONS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    add reverse reactions when :
    !        * ion-neutral reaction (no barrier)
    !        * 2 reactants <-> 2 products
    !        * 0 < DE <= DE_threshold
    !        * no radiative association
    !        * no recombination
    !        * reverse reaction is not already included
    !    the reverse reaction has the following :
    !        * same gamma, alpha as direct reaction
    !        * beta=DE(direct)*11600.
    !        * DE=0.0
    ! subroutine/function needed :
    !    REACTIONS_IDENTICAL
    ! input variables :
    ! ouput variables :
    ! results :
    !    react, Nreact, Nrever ...
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE

    REAL(KIND=DP), PARAMETER :: DE_threshold=2._DP ! eV
    INTEGER(KIND=LONG) :: i,j
    LOGICAL :: eliminate
    TYPE(TYPE_REACTION) :: reaction_aux

    ! the reactions are added at the end of the set
    b_rever=Nreact+1

    ! looks for ion-neutral reaction only -> type 'OTHER'
    DO i=b_other,e_other
       ! eliminates the following reactions
       eliminate= &
            (react(i)%Nprod_m1 /= 1) .OR. & ! more or less than 2 products
            (react(i)%DE <= Zero) .OR. &       ! DE <= 0
            (react(i)%DE > DE_threshold) .OR. &    ! DE > DE_threshold
            ((react(i)%R(1)>=b_neu) .AND. (react(i)%R(1)<=e_neu).AND.&
            (react(i)%R(2)>=b_neu) .AND. (react(i)%R(2)<=e_neu)).OR. & ! neutral-neutral (possibility of a barrier)
            (react(i)%R(1)==ind_e .OR. react(i)%R(2)==ind_e) .OR. & ! recombination
            (react(i)%P(1)==ind_PHOTON .OR. react(i)%P(2)==ind_PHOTON) ! radiative association

       ! tests if the reverse reaction is already in the set
       j=1
       DO WHILE ((j<=e_other) .AND. (.NOT. eliminate))
          ! build the reverse reaction
          reaction_aux%R(1)=react(i)%P(1)
          reaction_aux%R(2)=react(i)%P(2)
          reaction_aux%P(1)=react(i)%R(1)
          reaction_aux%P(2)=react(i)%R(2)
          reaction_aux%P(3:4)=P_none
          ! test if it exists already
          eliminate=REACTIONS_IDENTICAL(react(j),reaction_aux)
          j=j+1
       END DO

       ! if the reaction is not eliminated, we add the corresponding reverse reaction
       IF (.NOT. eliminate) THEN
          ! update of Nreact
          Nreact=Nreact+1
          ! Nreact_max must be large enough
          IF (Nreact > Nreact_max) STOP "*** WARNING, Nreact_max is too small"
          ! update of Nrever (e_rever is calculated at the end of the subroutine)
          Nrever=Nrever+1

          ! definition of the reverse reaction
          react(Nreact)=reaction_aux ! reactants and products
          react(Nreact)%gamma=react(i)%gamma
          react(Nreact)%alpha=react(i)%alpha
          react(Nreact)%beta=react(i)%DE*11600._DP
          react(Nreact)%DE=Zero
          react(Nreact)%Nprod_m1=1 ! only 2 products
          react(Nreact)%mass_prod= &
               speci(react(Nreact)%P(1))%mass + &
               speci(react(Nreact)%P(2))%mass
          react(Nreact)%type='REVER'
          react(Nreact)%ref='REVER'
       END IF
    END DO
    ! calculate e_rever
    e_rever=b_rever+Nrever-1

  END SUBROUTINE ADD_REVERSE_REACTIONS


END MODULE MODULE_CHEM_REACT
