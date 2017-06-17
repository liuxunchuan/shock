C---------------------------------------------------------------------------
      module m_atype
      
C   THIS MODULE CALCULATES THE POPULATION DENSITIES OF THE ROTATIONAL
C   LEVELS J.le.9 OF THE TORTIONAL GROUND STATE, nu = 0, OF A-type CH3OH
C   CALLED FROM SUBROUTINE DIFFUN

      public Atype_LVG
      
      contains

       SUBROUTINE Atype_LVG(ABH,ABH2,ABHE,ABCO,TN0,ROP,VGRAD,
     1                   NLEVCO,YTCO,ae_CH3OH,ICNT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C---------------------------------------------------------------------
C      XNH  IS NUMBER DENSITY OF HYDROGEN ATOMS (CM-3)
C      XNH2 IS NUMBER DENSITY OF HYDROGEN MOLECULES (CM-3)
C      XNHE IS NUMBER DENSITY OF HELIUM ATOMS (CM-3)
C      TN0 IS TEMPERATURE OF NEUTRAL FLUID
C      ROP IS THE ORTHO:PARA H_2 RATIO
C      TBB IS TEMPERATURE OF RADIATION FIELD

      REAL*8    EVLTOT(NLEVCO,NLEVCO)
      REAL*8    YTCO(NLEVCO)

      COMMON /CHKIT/ TOTO

C BACKGROUND RADIATION TEMPERATURE AND DILUTION FACTOR
      DATA TBB, WDIL / 2.73D0, 1.D0 /
      DATA ICOUNT       / 0 /

      ICOUNT=ICOUNT+1
      
C VELOCITY GRADIENT VGRAD (S-1)

C THE A:E CH3OH RATIO
       RAE_CH3OH = ae_CH3OH

C INITIALIZE THE A-type CH3OH ROTATIONAL LEVEL POPULATION DENSITIES TO THEIR
C VALUES IN LTE AT THE NEUTRAL TEMPERATURE (ON FIRST CALL ONLY)

      IF(ICOUNT.EQ.1) THEN
       CALL INIT_Atype(ABCO,RAE_CH3OH,TN0)
      ENDIF

      CALL MATRCO(ABH, ABH2, ABHE, ABCO, TN0, TBB, ROP, WDIL,
     1            VGRAD,RAE_CH3OH,EVLTOT)
      CALL EVOL_Atype (EVLTOT,YTCO)

C     WRITE(6, 2150) YTCO
C     WRITE(6, 2150) TOTO
C2150 FORMAT(6D12.4)
      ICNT = ICOUNT

      RETURN
      END SUBROUTINE Atype_LVG
C
C      %%%%%%%%%%%%%%%%%
       SUBROUTINE INIT_Atype(ABCO,RAE_CH3OH,TN0)
C      %%%%%%%%%%%%%%%%%

C  SOUS PROGRAMME D'INITIALISATION DE CO

C--------+---------+---------+---------+---------+---------+---------+-*-------+

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C   NLEVCO:  NUMBER OF ROTATIONAL LEVELS OF A-type CH3OH (nu = 0)

       PARAMETER (NLEVCO=256)
! The value of the PARAMETER NLEVCO must be the same as that of the corresponding 
! PARAMETER Natype_lev in evolution.f90, mhd_vode.f90 and outputs.f90
       PARAMETER (MLEVCO=1000)
       
C   NTEMP: NUMBER OF TEMPERATURES AT WHICH COLLISIONAL 
C   RATE COEFFICIENTS ARE TABULATED

       PARAMETER (NTEMP=20)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=20)

       CHARACTER*75       FICHIER
       REAL*8    ELHM(NLEVCO)
       REAL*8    ELK(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    freq(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    EVLDIP(MLEVCO,MLEVCO)
       REAL*8    frequ(MLEVCO,MLEVCO)
       INTEGER   NLEV(NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   JLEV(NLEVCO)
       INTEGER   KLEV(NLEVCO)
 
      COMMON /FCTa_CH3OH/ COVJ0(NLEVCO),freq
      COMMON /Atype/ GCO,ELK,JLEV,KLEV
      
      COMMON /TRANSF_a_CH3OH/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE,ELHM,NJLEV
      
      
       DATA MJCO / 50 /

       DATA IREAD2 / 51 /

C--------+---------+---------+---------+---------+---------+---------+-*-------+

C  INITIALIZE ARRAYS 

       DO 11 I=1,NLEVCO
         DO 10 J=1,NLEVCO
           DO 7 K=1,NTEMPH
             QHCO(K,J,I) = 0.0D0
    7    CONTINUE
           DO 8 K=1,NTMPHE
             QHECO(K,J,I) = 0.0D0
    8    CONTINUE
           DO 9 K=1,NTEMP
             QoH2CO(K,J,I) = 0.0D0
             QpH2CO(K,J,I) = 0.0D0
    9    CONTINUE
           EVLDI(J,I) = 0.0D0
           freq(J,I) = 0.0D0
   10    CONTINUE
   11  CONTINUE
   
       DO 13 I=1,MLEVCO
         DO 12 J=1,MLEVCO
           EVLDIP(J,I) = 0.0D0
           frequ(J,I) = 0.0D0
   12    CONTINUE
   13  CONTINUE

C  COLLISIONAL RATE COEFFICIENTS FOR HE + A-type CH3OH AND H2 + A-type CH3OH 

C  READ THE COLLISIONAL DE-EXCITATION RATE COEFFICIENTS (CM3 S-1)
c      FICHIER = 'input/QHCO.DAT              '
c      OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

c       READ (IREAD2,*) (TTABH(IT),IT=1,NTEMPH)
c      DO 17 ITRAN=1,136
c       READ (IREAD2,*)LEVU,LEVL,(QHCO(IT,LEVU,LEVL),
c    1        IT=1,NTEMPH) 
c  17  CONTINUE
c      CLOSE (IREAD2)
       
       FICHIER = 'input/Q_Atype_new.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

       DO 18 IT=1,NTEMP
        READ (IREAD2,*) TTABHE(IT)
       DO 17 LEVL=1,NLEVCO
        READ (IREAD2,*) NTRAN,(QHECO(IT,LEVU,LEVL),LEVU=1,NLEVCO)
c       print *, NTRAN,(QHECO(IT,LEVU,LEVL),LEVU=1,NLEVCO)
   17  CONTINUE
   18  CONTINUE
       CLOSE (IREAD2)
       
       FICHIER = 'input/QoH2_Atype_new.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

       DO 20 IT=1,NTEMP
         READ (IREAD2,*) TTAB(IT)
           DO 19 LEVL=1,NLEVCO
            READ (IREAD2,*) NTRAN,(QoH2CO(IT,LEVU,LEVL),LEVU=1,NLEVCO) 
c           print *, NTRAN,(QoH2CO(IT,LEVU,LEVL),LEVU=1,NLEVCO)
   19  CONTINUE
   20  CONTINUE
       CLOSE (IREAD2)

       FICHIER = 'input/QpH2_Atype_new.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

       DO 22 IT=1,NTEMP
         READ (IREAD2,*) TTAB(IT)
           DO 21 LEVL=1,NLEVCO
            READ (IREAD2,*) NTRAN,(QpH2CO(IT,LEVU,LEVL),LEVU=1,NLEVCO) 
c           print *, NTRAN,(QoH2CO(IT,LEVU,LEVL),LEVU=1,NLEVCO)
   21  CONTINUE
   22  CONTINUE
       CLOSE (IREAD2)


C*********************************************************************

C  LEVEL ENERGIES, IN CM-1, FOR A-TYPE, AND VALUES OF (J, K)
C  A+/- DENOTED BY THE SIGN OF K

       FICHIER = 'input/Atype_new.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

       DO 25 ILEV=1,NLEVCO
        READ (IREAD2,1000) NLEV(ILEV),JLEV(ILEV),KLEV(ILEV),ELHM(ILEV)
 1000 FORMAT(I4, I5, I4, F15.0)
c       print *, NLEV(ILEV),JLEV(ILEV),KLEV(ILEV),ELHM(ILEV)
        NJLEV(ILEV) = JLEV(ILEV)
C  convert level energies to K
        ELHM(ILEV)=ELHM(ILEV)*1.4388D0
        ELK(ILEV)=ELHM(ILEV)
   25  CONTINUE
      CLOSE (IREAD2)

C*********************************************************************

C  ELECTRIC DIPOLE TRANSITION PROBABILITIES FOR A-TYPE CH3OH
C  and line frequencies in GHz

       FICHIER = 'input/A_Atype_newer.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

        READ (IREAD2,*) NTRTOT
         IF(NTRTOT.GT.MLEVCO) THEN
          WRITE(6,2000)
 2000 FORMAT('CHECK DIMENSION MLEVCO OF THE A-TYPE TRANSITION 
     1        PROBABILITIES ARRAY EVLDIP')
          STOP
         ENDIF
       DO 30 ITRAN=1,NTRTOT
        READ (IREAD2,*) LEVU, LEVL, jup, kup, jlow, klow, 
     1                  EVLDIP(LEVU,LEVL), frequ(LEVU,LEVL), 
     2                  Eupcm, EupK
c       print *, LEVU, LEVL, jup, kup, jlow, klow, 
c    1                  EVLDIP(LEVU,LEVL), frequ(LEVU,LEVL), 
c    2                  Eupcm, EupK
   30  CONTINUE
      CLOSE (IREAD2)
      
C   ENSURE THAT THE INDICES OF THE TRANSITION PROBABILITIES EVLDI
C   AND OF THE LINE FREQUENCIES freq CORRESPOND TO 
C   THE INDICES OF THE COLLISIONAL RATE COEFFICIENTS
       DO 45 I=1,NLEVCO
         DO 40 J=1,NLEVCO
           EVLDI(J,I) = EVLDIP(NLEV(J),NLEV(I))
           freq(J,I) = frequ(NLEV(J),NLEV(I))
   40    CONTINUE
   45  CONTINUE

C*********************************************************************

C   INITIAL LEVEL POPULATIONS OF A-type CH3OH
C   HYPOTHESIS 1: LTE AT INITIAL TEMPERATURE

         ECH = 0.0D0
         DO 90 LEV=1,NLEVCO
           IJ = NJLEV(LEV)
           COVJ0(LEV) = (2.0D0*IJ+1.0D0) * DEXP(-ELHM(LEV)/TN0)
           ECH = ECH + COVJ0(LEV)
   90    CONTINUE

         DO 100 LEV=1,NLEVCO
           IJ = NJLEV(LEV)
           COVJ0(LEV) = ABCO * COVJ0(LEV) / ECH * 
     1                  RAE_CH3OH / (1. + RAE_CH3OH)
  100    CONTINUE

C   HYPOTHESIS 2: ENTIRE POPULATION IN GROUND LEVEL 


C        COVJ0(1) = ABCO * RAE_CH3OH / (1. + RAE_CH3OH) 

C        DO 110 LEV=2,NLEVCO
C            COVJ0(LEV) = 1.0D-25
C 110    CONTINUE

C     WRITE(6, 2150)(COVJ0(J),J=1,NLEVCO)
C2150 FORMAT(6D12.4)
C--------+---------+---------+---------+---------+---------+---------+-*-------+

       RETURN
       END SUBROUTINE INIT_Atype


C      %%%%%%%%%%%%%%%%%
       SUBROUTINE MATRCO(ABH, ABH2, ABHE, ABCO, TN0, TBB, ROP, WDIL,
     1                   VGRAD,RAE_CH3OH,EVLTOT)
C      %%%%%%%%%%%%%%%%%

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (NLEVCO=256)
! The value of the PARAMETER NLEVCO must be the same as 
! in SUBROUTINE INITCO 
       PARAMETER (NTEMP=20)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=20)

       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    freq(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    EVLCOL(NLEVCO,NLEVCO)
       REAL*8    EVLRAD(NLEVCO,NLEVCO)
       INTEGER   NJLEV(NLEVCO)

      COMMON /TRANSF_a_CH3OH/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE,ELHM,NJLEV

       DATA MJCO / 50 /

C------------------------------------------------------------------------------

C  CALCULATION OF THE ARRAY ELEMENTS NECESSARY TO DETERMINE THE EVOLUTION OF 
C  THE LEVEL POPULATIONS OF A-type CH3OH 

C     RADIATIVE TRANSITION PROBABILITIES (ELECTRIC DIPOLE): EVLDI
C     COLLISIONS: EVLCOL
C     TOTAL: EVLTOT

       DO 32 I1=1,NLEVCO
         DO 31 I2=1,NLEVCO
           EVLCOL(I2,I1) = 0.0D0
           EVLRAD(I2,I1) = 0.0D0
           EVLTOT(I2,I1) = 0.0D0
   31    CONTINUE
   32  CONTINUE

C      DETERMINE THE ELEMENTS OF THE ARRAYS EVLCOL, EVLRAD AND CEITQ  
C      IN SUBROUTINE COLCOV

       CALL COLCOV(ABH,ABH2,ABHE,ABCO,TN0,ROP,TBB,EVLCOL,
     1             EVLRAD,CEITQ,WDIL,VGRAD,RAE_CH3OH)

       DO 52 LEV1=1,NLEVCO
         IJ1 = NJLEV(LEV1)
         IF (IJ1 .GT. MJCO) GOTO 52
         DO 51 LEV2=1,NLEVCO
           IJ2 = NJLEV(LEV2)
           IF (IJ2 .GT. MJCO) GOTO 51
           EVLTOT(LEV2,LEV1) = EVLDI(LEV2,LEV1) + EVLCOL(LEV2,LEV1) 
     1                                           + EVLRAD(LEV2,LEV1)
   51    CONTINUE
           EVLTOT(LEV1,LEV1) = EVLTOT(LEV1,LEV1)
     1                       - AEITQ(LEV1) - CEITQ(LEV1)
   52  CONTINUE

       RETURN
       END SUBROUTINE MATRCO

C---------------------------------------------------------------------------

C      %%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE EVOL_Atype (EVLTOT,YTCO)
C      %%%%%%%%%%%%%%%%%%%%%%%%

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (NLEVCO=256)
! The value of the PARAMETER NLEVCO must be the same as 
! in SUBROUTINE INITCO 
       PARAMETER (NTEMP=20)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=20)

       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    freq(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    TTAB(NTEMP)
       REAL*8    YTCO(NLEVCO)
       INTEGER   NJLEV(NLEVCO)


      COMMON /FCTa_CH3OH/ FY(NLEVCO)

      COMMON /TRANSF_a_CH3OH/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE,ELHM,NJLEV

      COMMON /CHKIT/ TOTO

       DATA MJCO / 50 /

C------------------------------------------------------------------------------

C     WRITE(6, 2150)(FY(J),J=1,NLEVCO)
C2150 FORMAT(6D12.4)

C  CALCULATE THE CONTRIBUTIONS TO THE EVOLUTION OF THE LEVEL POPULATIONS OF 
C  A-type CH3OH (CM-3 S-1)

       DO 10 LEV=1,NLEVCO
         YTCO(LEV) = 0.0D0
   10  CONTINUE

       TOTO = 0.0D0
       LK1 = 0
       DO 42 LEV1=1,NLEVCO
         IJ1 = NJLEV(LEV1)
         IF (IJ1 .GT. MJCO) GOTO 42
         LK1 = LK1 + 1
         LK2 = 0
         DO 41 LEV2=1,NLEVCO
           IJ2 = NJLEV(LEV2)
           IF (IJ2 .GT. MJCO) GOTO 41
           LK2 = LK2 + 1
           YTCO(LK1) = YTCO(LK1) + FY(LK2) * EVLTOT(LEV2,LEV1)
   41    CONTINUE
         TOTO = TOTO + YTCO(LK1)
   42  CONTINUE

       RETURN
       END SUBROUTINE EVOL_Atype

C---------------------------------------------------------------------------

C      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE COLCOV (XNH,XNH2,XNHE,XNCO,TN0,ROP,TBB,EVLCOL,
     1                    EVLRAD,CEITQ,WDIL,VGRAD,RAE_CH3OH)
C      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C---------------------------------------------------------------------
C      XNH  IS NUMBER DENSITY OF HYDROGEN ATOMS (CM-3)
C      XNH2 IS NUMBER DENSITY OF HYDROGEN MOLECULES (CM-3)
C      XNHE IS NUMBER DENSITY OF HELIUM ATOMS (CM-3)
C      TN0 IS TEMPERATURE OF NEUTRAL FLUID
C      ROP IS THE ORTHO:PARA H_2 RATIO
C      TBB IS TEMPERATURE OF RADIATION FIELD
C      VGRAD IS ABSOLUTE VALUE OF LOCAL VELOCITY GRADIENT 

C      CALLED FROM SUBROUTINE MATRCO
C
C
C---------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       PARAMETER (NLEVCO=256)
! The value of the PARAMETER NLEVCO must be the same as 
! in SUBROUTINE INITCO 
       PARAMETER (NTEMP=20)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=20)

       REAL*8    ELHM(NLEVCO)
       REAL*8    ELK(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    freq(NLEVCO,NLEVCO)
       REAL*8    TAU(NLEVCO,NLEVCO)
       REAL*8    EP(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    HANDC(NLEVCO)
       REAL*8    ANDC(NLEVCO)
       REAL*8    TVGRAD(NLEVCO,NLEVCO)
       REAL*8    EVLCOL(NLEVCO,NLEVCO)
       REAL*8    EVLRAD(NLEVCO,NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   JLEV(NLEVCO)
       INTEGER   KLEV(NLEVCO)

      COMMON /FCTa_CH3OH/ COVJ(NLEVCO),freq,TAU,WWT,TVGRAD
      COMMON /Atype/ GCO,ELK,JLEV,KLEV
      
      COMMON /TRANSF_a_CH3OH/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE,ELHM,NJLEV
      
      DIMENSION check(NLEVCO)

      DATA XK,PI/1.38062D-16,3.1415926536D0/

       DATA MJCO / 50 /

C     PI=ACOS(-1.D0)

C     calculate the number densities (cm-3) of ortho- and para-H2
      XNoH2 = ROP * XNH2 / (1. + ROP)
      XNpH2 =       XNH2 / (1. + ROP)

C     NOTE: C = DE-EXCITATION; D = EXCITATION

c      DO 20 LEVU=2,NLEVCO
       DO 20 LEVU=1,NLEVCO
         JU = NJLEV(LEVU)
         IF (JU .GT. MJCO) GOTO 20

       check(LEVU) = 0.d0
c        DO 10 LEVL=1,LEVU-1
         DO 10 LEVL=1,NLEVCO
           IF (ELHM(LEVL) .GE. ELHM(LEVU)) GOTO 10
           JL = NJLEV(LEVL)
           IF (JL .GT. MJCO) GOTO 10

C  CTOD IS THE BOLTZMANN FACTOR WHICH RELATES THE
C  EXCITATION TO THE DE-EXCITATION RATE 

           NDEJ = ABS(JU - JL)
C  NEUTRALS: H, HE; H2
           FAC = (ELHM(LEVU) - ELHM(LEVL)) / TN0
            FAC=MIN(FAC,180.D0)
           CTOD = DEXP(-FAC) * (2.0D0*JU+1.0D0) / (2.0D0*JL+1.0D0)
C  INTRODUCE THE FACTOR 1/(EXP(HNU/KTBB)-1) NECESSARY TO CALCULATE
C  BUL FROM AUL
C  TBB IS THE TEMPERATURE OF THE BLACK-BODY RADIATION FIELD 
           FACBB = (ELHM(LEVU) - ELHM(LEVL)) / TBB
C  BEWARE OF OVERFLOW
           IF(FACBB.GT.200.D0) THEN
           FACBB = 0.D0
           ELSE
           FACBB = 1.D0/(DEXP(FACBB)-1.D0)
           ENDIF

C  CALCULATE AND SUM THE COLLISION RATES (S-1)

C  A-type - H

       IF(TN0.GE.TTABH(NTEMPH)) THEN
        GAMMA=QHCO(NTEMPH,LEVU,LEVL)
          GOTO 1
         ENDIF

       IF(TN0.LE.TTABH(1)) THEN
        GAMMA=QHCO(1,LEVU,LEVL)
          GOTO 1
         ENDIF

c  rate coefficient from linear interpolation
        GAMMA = spl(NTEMPH,TTABH(1),QHCO(1,LEVU,LEVL),TN0) 
      
    1  CCOH = GAMMA

C  A-type - He (Rabli & Flower 2010, MNRAS, 403, 2033; 2011, MNRAS, 411, 2093)
       IF(TN0.GE.TTABHE(NTMPHE)) THEN
        GAMMA=QHECO(NTMPHE,LEVU,LEVL)
          GOTO 2
         ENDIF

       IF(TN0.LE.TTABHE(1)) THEN
        GAMMA=QHECO(1,LEVU,LEVL)
          GOTO 2
         ENDIF

c  rate coefficient from linear interpolation
        GAMMA = spl(NTMPHE,TTABHE(1),QHECO(1,LEVU,LEVL),TN0) 

    2  CCOHE = GAMMA

C  A-type - o-H2 (Rabli & Flower 2010, MNRAS, 406, 95)

       IF(TN0.GE.TTAB(NTEMP)) THEN
        GAMMA=QoH2CO(NTEMP,LEVU,LEVL)
          GOTO 3
         ENDIF

       IF(TN0.LE.TTAB(1)) THEN
        GAMMA=QoH2CO(1,LEVU,LEVL)
          GOTO 3
         ENDIF

c  rate coefficient from linear interpolation
        GAMMA = spl(NTEMP,TTAB(1),QoH2CO(1,LEVU,LEVL),TN0) 
      
    3  CCOoH2 = GAMMA
    
C  A-type - p-H2 (Rabli & Flower 2010, MNRAS, 406, 95)

       IF(TN0.GE.TTAB(NTEMP)) THEN
        GAMMA=QpH2CO(NTEMP,LEVU,LEVL)
          GOTO 4
         ENDIF

       IF(TN0.LE.TTAB(1)) THEN
        GAMMA=QpH2CO(1,LEVU,LEVL)
          GOTO 4
         ENDIF

c  rate coefficient from linear interpolation
        GAMMA = spl(NTEMP,TTAB(1),QpH2CO(1,LEVU,LEVL),TN0) 
      
    4  CCOpH2 = GAMMA
    
C  SUM THE INDIVIDUAL CONTRIBUTIONS TO OBTAIN THE TOTAL COLLISION RATE (S-1):
           C = CCOoH2 * XNoH2
     1       + CCOpH2 * XNpH2
     2       + CCOHE  * XNHE
     3       + CCOoH2 * XNH
c     print *, TN0, levu, levl, CCOoH2, CCOpH2, CCOHE, EVLDI(LEVU,LEVL)
C     WRITE(6, 2150) TN0, cm
C2150 FORMAT(6D12.4)
C     EXCITATION RATES:
C     NEUTRALS
           D = C * CTOD

C     
           EVLCOL(LEVL,LEVU) = D 
C          
           EVLCOL(LEVU,LEVL) = C 
C          
C  OPTICAL DEPTH OF THE TRANSITION LEVU TO LEVL
C  VGRAD IS THE LOCAL VELOCITY GRADIENT (S-1)

C  EVALUATE THE CUBE OF THE WAVE NUMBER (CM-1) OF THE TRANSITION
C  DIVISION BY 1.4388 CONVERTS K TO CM-1
           WN3 = ((ELHM(LEVU) - ELHM(LEVL))/1.4388D0)**3.D0
C  EVALUATE THE OPTICAL DEPTH IN THE LINE
C  INTRODUCE THE FACTOR LAMBDA^3/(8*PI*VGRAD) 
           R_1 = COVJ(LEVL)*(2.0D0*JU+1.0D0)/
     2  (COVJ(LEVU)*(2.0D0*JL+1.0D0)) - 1.D0
           TAU(LEVU,LEVL) = 1.D0/(8.D0*PI*VGRAD*WN3) * 
     1     COVJ(LEVU) * R_1 * EVLDI(LEVU,LEVL)

         TS=TAU(LEVU,LEVL)
C           WRITE(6, *) WN3 , R_1, COVJ(LEVU), COVJ(LEVL), TS
C  EVALUATE THE ESCAPE PROBABILITY EP(LEVU,LEVL)
C  TAU < 1E-2:
         EP(LEVU,LEVL)=1.-TS/2.D0
C        EP(LEVU,LEVL)=1.D0
         TSA=ABS(TS)
         IF(TSA.LT.0.01D0) GO TO 11
         IF(TSA.LT.100.D0) GO TO 12
C  TAU > 100:
         EP(LEVU,LEVL)=1.D0/TS
         GO TO 11
C  1E-2 < TAU < 100:
 12      EP(LEVU,LEVL)=(1.D0-EXP(-TS))/TS
 11      E1=EP(LEVU,LEVL)
C  NEUFELD & KAUFMAN (1993) EXPRESSION FOR THE ESCAPE PROBABILITY
C        E1 = 1.D0/(1.D0 + 3.D0*TS)
C        EP(LEVU,LEVL) = E1
C           WRITE(6, *) VGRAD , TS, E1
         T1=(1.D0-E1)/R_1
C  T2 = ESCAPE_PROBABILITY/[EXP([E(J)-E(J-1)]/KT_BB)-1]
         T2=E1 * FACBB
C           WRITE(11, *) VGRAD, T1, T2
C  CALCULATE THE STIMULATED EMISSION AND ABSORPTION RATES (S-1)
C  [CF. THE SPONTANEOUS EMISSION RATE EVLDI (S-1)]
C  WHEN THE ESCAPE PROBABILITY EP(J) = 1, T1 = 0 AND
C  T2 = 1/[EXP([E(J)-E(J-1)]/KT_BB)-1]; THEN, THE SOURCE FUNCTION IS EQUAL 
C  TO THE BLACK-BODY RADIATION INTENSITY B_NU = [2HNU**3/C**2]/
C  [EXP([E(J)-E(J-1)]/KT_BB)-1]
C  WHEN THE ESCAPE PROBABILITY EP(J) = 0, T2 = 0 AND T1 = 1/(R-1); THEN 
C  THE SOURCE FUNCTION J_NU/KAPPA_NU = (2HC)(NU**3/C**3)/(R-1)
C  [2HC = 3.97297D-16]

C  INDUCED RADIATIVE TERM, BUL; WDIL IS THE RADIATION DILUTION FACTOR,
C  E1 IS THE ESCAPE PROBABILITY
           EVLRAD(LEVU,LEVL) = WDIL * T2
           EVLRAD(LEVU,LEVL) = (EVLRAD(LEVU,LEVL) + T1) * 
     1     EVLDI(LEVU,LEVL)

C  INTRODUCE STATISTICAL WEIGHT RATIO (FOR BLU)

           EVLRAD(LEVL,LEVU) = EVLRAD(LEVU,LEVL) * (2.0D0*JU+1.0D0)
     1                                           / (2.0D0*JL+1.0D0)

C  CALCULATE THE OBSERVED LINE TEMPERATURES (K): GUSDORF ET AL. (2008a)
       TLINE = EP(LEVU,LEVL)/(8.D0*PI*VGRAD*WN3)*
     1         EVLDI(LEVU,LEVL)*COVJ(LEVU)
       TLINE = TLINE - (1.D0-EXP(-TS))*FACBB
       TLINE = TLINE*(ELHM(LEVU) - ELHM(LEVL))
C  MULTIPLY BY THE VELOCITY GRADIENT, IN KM S-1 CM-1
C  THERE REMAINS THE INTEGRAL WRT Z (CM) IN ORDER TO OBTAIN
C  THE INTEGRAL LINE INTENSITIES, TdV (K KM S-1)
C      TVGRAD(LEVU,LEVL) = TLINE*VGRAD*1.D-5
C  OUTPUT THE LINE INTENSITIES (K)
       TVGRAD(LEVU,LEVL) = TLINE
       check(LEVU) = check(LEVU) + 8.d5*PI*WN3*TVGRAD(LEVU,LEVL)

   10    CONTINUE
   20  CONTINUE

       DO 32 LEVU=1,NLEVCO
         NJU = NJLEV(LEVU)
         IF (NJU .GT. MJCO) GOTO 32
         AEITQ(LEVU) = 0.0D0
         CEITQ(LEVU) = 0.0D0
         HANDC(LEVU) = 0.0D0
         ANDC(LEVU) = 0.0D0
         DO 31 LEVL=1,NLEVCO
           NJL = NJLEV(LEVL)
           IF (NJL .GT. MJCO) GOTO 31
C  COLLISIONAL AND RADIATIVE RATES (S-1)
             AEITQ(LEVU) = AEITQ(LEVU) + EVLDI(LEVU,LEVL) 
             CEITQ(LEVU) = CEITQ(LEVU) + EVLCOL(LEVU,LEVL) 
     1                                 + EVLRAD(LEVU,LEVL)
C  HEATING AND COOLING RATES (K S-1)
             HANDC(LEVU) = HANDC(LEVU) + EVLCOL(LEVU,LEVL) 
     1                                 * (ELHM(LEVU) - ELHM(LEVL))
             ANDC(LEVU) = ANDC(LEVU) + (EVLDI(LEVU,LEVL) 
     1            + EVLRAD(LEVU,LEVL)) * (ELHM(LEVU) - ELHM(LEVL))
   31  CONTINUE
   32  CONTINUE
c     print *, (i, check(i), ANDC(i), i=1,NLEVCO)
C  HEATING AND COOLING RATES (ERG CM-3 S-1)

      verif = 0.0d0
      heat = 0.0d0
      cool = 0.0d0
      GCO = 0.0
      GCOP = 0.0
      DO 53 LEV=1,NLEVCO
        IJ = NJLEV(LEV)
        IF (IJ .GT. MJCO) GOTO 53
c  COVJ are the level populations (cm-3); verify that their sum is equal
c  to the density of A-type CH3OH
          verif = verif + COVJ(LEV)
c  heating and cooling rates (erg cm-3 s-1)
          if(handc(LEV).gt.0.d0)then
           heat = heat + COVJ(LEV) * handc(LEV)
          else
           cool = cool + COVJ(LEV) * handc(LEV)
          endif
       GCO = GCO + COVJ(LEV)*ANDC(LEV)
       GCOP = GCOP + check(LEV)
   53 CONTINUE

          heat = heat * XK
          cool = cool * XK

C---RADIATIVE COOLING BY A-type CH3OH (ERG CM-3 S-1)
c  GCO is the rate of energy loss by A-type CH3OH to the radiation field
c  WWT is the net cooling rate, evaluated as (minus) the sum of
c  (positive) heating rate, HEAT, and the (negative) cooling 
c  rate, COOL 
c  GCO and WWT are both positive when cooling is taking place
c  and should be equal if the cooling is significant

          GCO = GCO * XK
          GCOP = GCOP * XK
          WWT = -(HEAT + COOL)

C  VERIFY THAT THE SUM OF THE A-type CH3OH POPULATION DENSITIES (CM-3) 
C  IS EQUAL TO THE DENSITY OF THIS SPECIES
      verif = verif / (XNCO * RAE_CH3OH / ( 1. + RAE_CH3OH))
      if (abs(verif-1.0d0) .gt. 1.0d-2) then
        print *, '  verif A-type CH3OH:', verif, verif-1.0d0
      endif

       RETURN
       END SUBROUTINE COLCOV
       
      function spl(n,x,y,t) 
      implicit real * 8 ( a-h , o-z) 
      dimension x(n), y(n) 
      k=2 
   10 if(t.le.x(k)) go to 20 
      k=k+1 
      go to 10 
   20 e=x(k)-x(k-1) 
      f=x(k)-t 
      g=t-x(k-1) 
      spl=(y(k)*g+y(k-1)*f)/e 
      return 
      end function spl

       end module m_atype

C-------------------------------------------------------------------------------
