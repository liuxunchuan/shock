C---------------------------------------------------------------------------
      module msio
      
C   THIS MODULE CALCULATES THE POPULATION DENSITIES OF THE 
C   ROTATIONAL LEVELS OF THE X VIBRATIONAL GROUND STATE OF SiO
C   CALLED FROM SUBROUTINE DIFFUN

      public SiO_LVG
C     public INITCO
C     public EVOLCO
      
      contains

       SUBROUTINE SiO_LVG(ABH,ABH2,ABHE,ABCO,TN0,ROP,VGRAD,
     1                   NLEVCO,YTCO,ICNT)
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

C INITIALIZE THE SiO ROTATIONAL LEVEL POPULATION DENSITIES TO THEIR
C VALUES IN LTE AT THE NEUTRAL TEMPERATURE (ON FIRST CALL ONLY)

      IF(ICOUNT.EQ.1) CALL INITCO(ABCO,TN0)
   
      CALL MATRCO(ABH, ABH2, ABHE, ABCO, TN0, TBB, ROP, WDIL,
     1            VGRAD,EVLTOT)
      CALL EVOLCO (EVLTOT,YTCO)

C     WRITE(6, 2150) YTCO
C     WRITE(6, 2150) TOTO
C2150 FORMAT(6D12.4)
      ICNT = ICOUNT
C      print *, ICOUNT, YTCO(1:NLEVCO)    

      RETURN
      END SUBROUTINE SiO_LVG
C
C      %%%%%%%%%%%%%%%%%
       SUBROUTINE INITCO(ABCO,TN0)
C      %%%%%%%%%%%%%%%%%

C  INITIALIZATION SUBROUTINE

C--------+---------+---------+---------+---------+---------+---------+-*-------+

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C   NLEVCO:  NUMBER OF ROVIBRATIONAL LEVELS OF SiO (X)
C   NVXCO :  HIGHEST VIBRATIONAL LEVEL OF SiO (X)

       PARAMETER (NLEVCO=41)
! The value of the PARAMETER NLEVCO must be the same as that of the corresponding 
! PARAMETER NSiO_lev in evolution.f90, mhd_vode.f90 and outputs.f90
       PARAMETER (NVXCO=0)
       
C   NTEMP: NUMBER OF TEMPERATURES AT WHICH COLLISIONAL 
C   RATE COEFFICENTS ARE TABULATED

       PARAMETER (NTEMP=30)
       PARAMETER (NTMPH2=8)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       CHARACTER*75       FICHIER
       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QH2CO(NTMPH2,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    AR(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    TTABH2(NTMPH2)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)

      COMMON /FCTSiO/ COVJ0(NLEVCO)
      COMMON /SiO/ GCO,ELHM,NJLEV
      
      COMMON /TRANSF_SiO/QHCO,QHECO,QoH2CO,QpH2CO,QH2CO,
     1               EVLDI,AEITQ,TTAB,
     2               TTABH2,TTABH,TTABHE
      
      
       DATA MJCO / 40 /


C  DATA FOR THE ENERGY LEVEL INDEX, LEVEL, THE ROTATIONAL QUANTUM NUMBER, 
C  NJLEV, THE LEVEL ENERGIES, ELHM, AND THE EINSTEIN A-COEFFICIENTS, AR, 
C  ARE SUPPLIED FOR J = 0, 1, ... 40 IN THE VIBRATIONAL MANIFOLD V = 0  
       DATA (LEVEL(0,J), J=0,NLEVCO-1) /   1, 2, 3, 4, 5, 6, 7, 8, 9,
     1          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
     2          22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
     3          34, 35, 36, 37, 38, 39, 40, 41 /

C  VALUES OF THE QUANTUM NUMBER, J   
       DATA (NJLEV(N), N=1,NLEVCO) /
     1       0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
     2      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
     3      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 /

C  ENERGIES, IN CM-1
       DATA (ELHM(N), N=1,NLEVCO) /
     1       0.0000, 1.4485, 4.3454, 8.6907, 14.4843, 21.7261,
     2      30.4161, 40.5540, 52.1397, 65.1730, 79.6537, 95.5815,
     3     112.9562, 131.7776, 152.0452, 173.7587, 196.9179,
     4     221.5222, 247.5714, 275.0649, 304.0024, 334.3833,
     5     366.2071, 399.4734, 434.1815, 470.3308, 507.9209,
     6     546.9511, 587.4206, 629.3289, 672.6753, 717.4590,
     7     763.6793, 811.3355, 860.4266, 910.9520, 962.9108,
     8    1016.3022, 1071.1252, 1127.3789, 1185.0624 /

       DATA(AR(N),N=1,NLEVCO)/ 0.0, 3.0490e-06,2.9273e-05,
     2      1.0584e-04,2.6019e-04,5.1965e-04,9.1161e-04,1.4638e-03,
     3      2.2030e-03,3.1570e-03,4.3529e-03,5.8168e-03,7.5787e-03,
     4      9.6622e-03,1.2097e-02,1.4912e-02,1.8129e-02,2.1779e-02,
     5      2.5882e-02,3.0479e-02,3.5583e-02,4.1223e-02,4.7434e-02,
     6      5.4225e-02,6.1643e-02,6.9699e-02,7.8441e-02,8.7868e-02,
     7      9.8014e-02,1.0890e-01,1.2057e-01,1.3302e-01,1.4631e-01,
     8      1.6048e-01,1.7548e-01,1.9138e-01,2.0823e-01,2.2603e-01,
     9      2.4477e-01,2.6450e-01,2.8527e-01 /

       DATA IREAD2 / 51 /

C--------+---------+---------+---------+---------+---------+---------+-*-------+


C  COLLISIONAL RATE COEFFICIENTS FOR H + SiO, HE + SiO AND H2 + SiO

       DO 12 I=1,NLEVCO
         DO 11 J=1,NLEVCO
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
           DO 10 K=1,NTMPH2
             QH2CO(K,J,I) = 0.0D0
   10    CONTINUE
           EVLDI(J,I) = 0.0D0
   11    CONTINUE
   12  CONTINUE

C  READ THE COLLISIONAL DE-EXCITATION RATE COEFFICIENTS (CM3 S-1)
       
       FICHIER = 'input/QpH2SiO.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

        READ (IREAD2,*) (TTAB(IT),IT=1,NTEMP)
       DO 22 ITRAN=1,820
        READ (IREAD2,*)NTRAN,LEVU,LEVL,(QpH2CO(IT,LEVU,LEVL),
     1        IT=1,NTEMP) 
   22  CONTINUE
       CLOSE (IREAD2)

       FICHIER = 'input/QH2SiO.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

        READ (IREAD2,*) (TTABH2(IT),IT=1,NTMPH2)
       DO 23 ITRAN=1,210
        READ (IREAD2,*)LEVU,LEVL,(QH2CO(IT,LEVU,LEVL),
     1        IT=1,NTMPH2) 
   23  CONTINUE
       CLOSE (IREAD2)

C*********************************************************************

C  RADIATIVE ELECTRIC DIPOLE TRANSITION PROBABILITIES OF SiO


         DO 31 I=1,NLEVCO
           DO 30 J=1,NLEVCO
           IF ((J-I) .NE. 1) GO TO 30
             EVLDI(J,I) = AR(J)
             AEITQ(J) = AR(J)
C            EVLDI(J,I) = 10. * AR(J)
C            AEITQ(J) = 10. * AR(J)
   30    CONTINUE
   31    CONTINUE

C*********************************************************************

C   INITIAL POPULATIONS OF SiO
C   HYPOTHESIS 1: LTE AT INITIAL TEMPERATURE

         ECH = 0.0D0
         DO 90 LEV=1,NLEVCO
           IJ = NJLEV(LEV)
           COVJ0(LEV) = (2.0D0*IJ+1.0D0) * DEXP(-(1.4388*ELHM(LEV))/TN0)
           ECH = ECH + COVJ0(LEV)
   90    CONTINUE

         DO 100 LEV=1,NLEVCO
           IJ = NJLEV(LEV)
           COVJ0(LEV) = ABCO * COVJ0(LEV) / ECH
  100    CONTINUE

C   HYPOTHESIS 2: ENTIRE POPULATION IN J = 0 


C        COVJ0(1) = ABCO 

C        DO 110 LEV=2,NLEVCO
C            COVJ0(LEV) = 1.0D-15
C 110    CONTINUE

C   HYPOTHESIS 3: STATISTICAL EQUILIBRIUM 


C      FICHIER = 'input/CO_levpop               '
C      OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

C       READ (IREAD2,*) (COVJ0(LEV),LEV=1,NLEVCO)
C      CLOSE (IREAD2)

C     WRITE(6, 2150)(COVJ0(J),J=1,NLEVCO)
C2150 FORMAT(6D12.4)
C--------+---------+---------+---------+---------+---------+---------+-*-------+

       RETURN
       END SUBROUTINE INITCO


C      %%%%%%%%%%%%%%%%%
       SUBROUTINE MATRCO(ABH, ABH2, ABHE, ABCO, TN0, TBB, ROP, WDIL,
     1                   VGRAD,EVLTOT)
C      %%%%%%%%%%%%%%%%%

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (NLEVCO=41)
       PARAMETER (NVXCO=0)
! The values of the PARAMETERS NLEVCO and NVXCO must be the same as 
! in SUBROUTINE INITCO 
       PARAMETER (NTEMP=30)
       PARAMETER (NTMPH2=8)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       REAL*8    ELHM(41)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QH2CO(NTMPH2,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    TTABH2(NTMPH2)
       REAL*8    EVLCOL(NLEVCO,NLEVCO)
       REAL*8    EVLRAD(NLEVCO,NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)

      COMMON /TRANSF_SiO/QHCO,QHECO,QoH2CO,QpH2CO,QH2CO,
     1               EVLDI,AEITQ,TTAB,
     2               TTABH2,TTABH,TTABHE

       DATA MJCO / 40 /


       DATA (LEVEL(0,J), J=0,NLEVCO-1) /   1, 2, 3, 4, 5, 6, 7, 8, 9,
     1          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
     2          22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
     3          34, 35, 36, 37, 38, 39, 40, 41 /

C  VALUES OF THE QUANTUM NUMBER, J   
       DATA (NJLEV(N), N=1,NLEVCO) /
     1       0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
     2      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
     3      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 /

C  ENERGIES, IN CM-1
       DATA (ELHM(N), N=1,NLEVCO) /
     1       0.0000, 1.4485, 4.3454, 8.6907, 14.4843, 21.7261,
     2      30.4161, 40.5540, 52.1397, 65.1730, 79.6537, 95.5815,
     3     112.9562, 131.7776, 152.0452, 173.7587, 196.9179,
     4     221.5222, 247.5714, 275.0649, 304.0024, 334.3833,
     5     366.2071, 399.4734, 434.1815, 470.3308, 507.9209,
     6     546.9511, 587.4206, 629.3289, 672.6753, 717.4590,
     7     763.6793, 811.3355, 860.4266, 910.9520, 962.9108,
     8    1016.3022, 1071.1252, 1127.3789, 1185.0624 /



C------------------------------------------------------------------------------

C  CALCULATION OF THE ARRAY ELEMENTS NECESSARY TO DETERMINE THE EVOLUTION OF 
C  THE LEVEL POPULATIONS OF SiO 

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
     1             EVLRAD,CEITQ,WDIL,VGRAD)

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

C     WRITE(6, 2150) EVLDI
C2150 FORMAT(6D12.4)
       RETURN
       END SUBROUTINE MATRCO

C---------------------------------------------------------------------------

C      %%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE EVOLCO (EVLTOT,YTCO)
C      %%%%%%%%%%%%%%%%%%%%%%%%

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (NLEVCO=41)
       PARAMETER (NVXCO=0)
! The values of the PARAMETERS NLEVCO and NVXCO must be the same as 
! in SUBROUTINE INITCO 
       PARAMETER (NTEMP=30)
       PARAMETER (NTMPH2=8)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       REAL*8    ELHM(41)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QH2CO(NTMPH2,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    TTABH2(NTMPH2)
       REAL*8    TTAB(NTEMP)
       REAL*8    YTCO(NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)


      COMMON /FCTSiO/ FY(NLEVCO)

      COMMON /TRANSF_SiO/QHCO,QHECO,QoH2CO,QpH2CO,QH2CO,
     1               EVLDI,AEITQ,TTAB,
     2               TTABH2,TTABH,TTABHE

      COMMON /CHKIT/ TOTO

       DATA MJCO / 40 /


       DATA (LEVEL(0,J), J=0,NLEVCO-1) /   1, 2, 3, 4, 5, 6, 7, 8, 9,
     1          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
     2          22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
     3          34, 35, 36, 37, 38, 39, 40, 41 /

C  VALUES OF THE QUANTUM NUMBER, J   
       DATA (NJLEV(N), N=1,NLEVCO) /
     1       0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
     2      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
     3      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 /

C  ENERGIES, IN CM-1
       DATA (ELHM(N), N=1,NLEVCO) /
     1       0.0000, 1.4485, 4.3454, 8.6907, 14.4843, 21.7261,
     2      30.4161, 40.5540, 52.1397, 65.1730, 79.6537, 95.5815,
     3     112.9562, 131.7776, 152.0452, 173.7587, 196.9179,
     4     221.5222, 247.5714, 275.0649, 304.0024, 334.3833,
     5     366.2071, 399.4734, 434.1815, 470.3308, 507.9209,
     6     546.9511, 587.4206, 629.3289, 672.6753, 717.4590,
     7     763.6793, 811.3355, 860.4266, 910.9520, 962.9108,
     8    1016.3022, 1071.1252, 1127.3789, 1185.0624 /


C------------------------------------------------------------------------------

C     WRITE(6, 2150)(FY(J),J=1,NLEVCO)
C2150 FORMAT(6D12.4)

C  CALCULATE THE CONTRIBUTIONS TO THE EVOLUTION OF THE LEVEL POPULATIONS OF SiO (CM-3 S-1)

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
       END SUBROUTINE EVOLCO

C---------------------------------------------------------------------------

C      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE COLCOV (XNH,XNH2,XNHE,XNCO,TN0,ROP,TBB,EVLCOL,
     1                    EVLRAD,CEITQ,WDIL,VGRAD)
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

       PARAMETER (NLEVCO=41)
       PARAMETER (NVXCO=0)
! The values of the PARAMETERS NLEVCO and NVXCO must be the same as 
! in SUBROUTINE INITCO 
       PARAMETER (NTEMP=30)
       PARAMETER (NTMPH2=8)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QH2CO(NTMPH2,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    TAU(NLEVCO,NLEVCO)
       REAL*8    EP(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    TTABH2(NTMPH2)
       REAL*8    HANDC(NLEVCO)
       REAL*8    ANDC(NLEVCO)
       REAL*8    TVGRAD(NLEVCO)
       REAL*8    EVLCOL(NLEVCO,NLEVCO)
       REAL*8    EVLRAD(NLEVCO,NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)

      COMMON /FCTSiO/ COVJ(NLEVCO),TAU,WWT,TVGRAD
      COMMON /SiO/ GCO,ELHM,NJLEV
      
      COMMON /TRANSF_SiO/QHCO,QHECO,QoH2CO,QpH2CO,QH2CO,
     1               EVLDI,AEITQ,TTAB,
     2               TTABH2,TTABH,TTABHE
      

      DATA XK,PI/1.38062D-16,3.1415926536D0/


       DATA MJCO / 40 /


       DATA (LEVEL(0,J), J=0,NLEVCO-1) /   1, 2, 3, 4, 5, 6, 7, 8, 9,
     1          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
     2          22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
     3          34, 35, 36, 37, 38, 39, 40, 41 /


C     PI=ACOS(-1.D0)

C     calculate the number densities (cm-3) of ortho- and para-H2
      XNoH2 = ROP * XNH2 / (1. + ROP)
      XNpH2 =       XNH2 / (1. + ROP)
      
C     NOTE: C = DE-EXCITATION; D = EXCITATION

       DO 20 LEVU=2,NLEVCO
         JU = NJLEV(LEVU)
         IF (JU .GT. MJCO) GOTO 20

         DO 10 LEVL=1,LEVU-1
           JL = NJLEV(LEVL)
           IF (JL .GT. MJCO) GOTO 10

C  CTOD IS THE BOLTZMANN FACTOR WHICH RELATES THE
C  EXCITATION TO THE DE-EXCITATION RATE 

           NDEJ = ABS(JU - JL)
C  NEUTRALS: H, HE; H2
           FAC = 1.4388 * (ELHM(LEVU) - ELHM(LEVL)) / TN0
            FAC=MIN(FAC,180.D0)
           CTOD = DEXP(-FAC) * (2.0D0*JU+1.0D0) / (2.0D0*JL+1.0D0)
C  INTRODUCE THE FACTOR 1/(EXP(HNU/KTBB)-1) NECESSARY TO CALCULATE
C  BUL FROM AUL
C  TBB IS THE TEMPERATURE OF THE BLACK-BODY RADIATION FIELD 
           FACBB = 1.4388 * (ELHM(LEVU) - ELHM(LEVL)) / TBB
C  BEWARE OF OVERFLOW
           IF(FACBB.GT.200.D0) THEN
           FACBB = 0.D0
           ELSE
           FACBB = 1.D0/(DEXP(FACBB)-1.D0)
           ENDIF

C  CALCULATE AND SUM THE COLLISION RATES (S-1)

C  SiO - H

       IF(TN0.GE.TTABH(NTEMPH)) THEN
        GAMMA=QHCO(NTEMPH,LEVU,LEVL)
          GOTO 1
         ENDIF

       IF(TN0.LE.TTABH(1)) THEN
        GAMMA=QHCO(1,LEVU,LEVL)
          GOTO 1
         ENDIF

c  rate coefficient from linear fit
        GAMMA = spl(NTEMPH,TTABH(1),QHCO(1,LEVU,LEVL),TN0) 
      
    1  CCOH = GAMMA

C  SiO - HE

       IF(TN0.GE.TTABHE(NTMPHE)) THEN
        GAMMA=QHECO(NTMPHE,LEVU,LEVL)
          GOTO 2
         ENDIF

       IF(TN0.LE.TTABHE(1)) THEN
        GAMMA=QHECO(1,LEVU,LEVL)
          GOTO 2
         ENDIF

c  rate coefficient from linear fit
        GAMMA = spl(NTMPHE,TTABHE(1),QHECO(1,LEVU,LEVL),TN0) 
      
    2  CCOHE = GAMMA

C  SiO - o-H2

       IF(TN0.GE.TTABH2(NTMPH2)) THEN
        GAMMA=QH2CO(NTMPH2,LEVU,LEVL)
          GOTO 3
         ENDIF

       IF(TN0.LE.TTABH2(1)) THEN
        GAMMA=QH2CO(1,LEVU,LEVL)
          GOTO 3
         ENDIF

c  rate coefficient from linear fit
        GAMMA = spl(NTMPH2,TTABH2(1),QH2CO(1,LEVU,LEVL),TN0) 
      
    3  CCOoH2 = GAMMA
c     print *,ELHM(LEVU),ELHM(LEVL),TN0,LEVU,LEVL,GAMMA,EVLDI(LEVU,LEVL)
    
C  SiO - p-H2 (Dayou & Balanca (2006, A&A, 459, 297), computed for collisions with He, 
C  multiplied by factor of 1.38 to simulate collisions with p-H2)

       IF(TN0.GE.TTAB(NTEMP)) THEN
        GAMMA=QpH2CO(NTEMP,LEVU,LEVL)
          GOTO 4
         ENDIF

       IF(TN0.LE.TTAB(1)) THEN
        GAMMA=QpH2CO(1,LEVU,LEVL)
          GOTO 4
         ENDIF

c  rate coefficient from linear fit
        GAMMA = spl(NTEMP,TTAB(1),QpH2CO(1,LEVU,LEVL),TN0) 
      
    4  CCOpH2 = GAMMA
    

C  SUM THE INDIVIDUAL CONTRIBUTIONS TO OBTAIN THE TOTAL COLLISION RATE (S-1):
C  QH2CO (para-)H2 rate coefficients from Turner et al. (1992, ApJ, 399, 114) are not used
C  the rate coefficients of Dayou & Balanca (2006, A&A, 459, 297) adopted for He:
C  their results were scaled by factor of 1.38 for para-H2 (QpH2CO)
C  para-H2 rate coefficients adopted for ortho-H2 and H
           C = CCOpH2 * XNoH2
     1       + CCOpH2 * XNpH2
C  remove scaling factor of 1.38 to recover rate coefficients for collisions with He
     2       + CCOpH2 * XNHE /1.38
     3       + CCOpH2 * XNH
     
c     print *, TN0,JU,JL,CCOoH2,CCOpH2,CCOHE,CCOH,C,EVLDI(LEVU,LEVL)

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
           WN3 = (ELHM(LEVU) - ELHM(LEVL))**3.D0
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
C  TLINE IS NON-ZERO ONLY FOR LEVL = LEVU - 1
       TLINE = EP(LEVU,LEVL)/(8.D0*PI*VGRAD*WN3)*
     1         EVLDI(LEVU,LEVL)*COVJ(LEVU)
       TLINE = TLINE - (1.D0-EXP(-TS))*FACBB
       TLINE = TLINE * 1.4388 * (ELHM(LEVU) - ELHM(LEVL))
C  MULTIPLY BY THE VELOCITY GRADIENT, IN KM S-1 CM-1
C  THERE REMAINS THE INTEGRAL WRT Z (CM) IN ORDER TO OBTAIN
C  THE INTEGRAL LINE INTENSITIES, TdV (K KM S-1)
C      TVGRAD(LEVU) = TLINE*VGRAD*1.D-5
C  OUTPUT THE LINE INTENSITIES (K)
       TVGRAD(LEVU) = TLINE

   10    CONTINUE
   20  CONTINUE

       DO 32 LEVU=1,NLEVCO
         NJU = NJLEV(LEVU)
         IF (NJU .GT. MJCO) GOTO 32
         CEITQ(LEVU) = 0.0D0
         HANDC(LEVU) = 0.0D0
         ANDC(LEVU) = 0.0D0
         DO 31 LEVL=1,NLEVCO
           NJL = NJLEV(LEVL)
           IF (NJL .GT. MJCO) GOTO 31
C  COLLISIONAL AND INDUCED RADIATIVE RATES (S-1)
             CEITQ(LEVU) = CEITQ(LEVU) + EVLCOL(LEVU,LEVL) 
     1                                 + EVLRAD(LEVU,LEVL)
C  HEATING AND COOLING RATES (K S-1)
             HANDC(LEVU) = HANDC(LEVU) + EVLCOL(LEVU,LEVL) 
     1                     * 1.4388 * (ELHM(LEVU) - ELHM(LEVL))
             ANDC(LEVU) = ANDC(LEVU) + (EVLDI(LEVU,LEVL) * EP(LEVU,LEVL)
     1       + EVLRAD(LEVU,LEVL)) * 1.4388 * (ELHM(LEVU) - ELHM(LEVL))
   31  CONTINUE
   32  CONTINUE
C  HEATING AND COOLING RATES (ERG CM-3 S-1)

      verif = 0.0d0
      heat = 0.0d0
      cool = 0.0d0
      GCO = 0.0
      DO 53 LEV=1,NLEVCO
        IJ = NJLEV(LEV)
        IF (IJ .GT. MJCO) GOTO 53
c  COVJ are the level populations (cm-3); verify that their sum is equal
c  to the density of SiO
          verif = verif + COVJ(LEV)
c  heating and cooling rates (erg cm-3 s-1)
          if(handc(LEV).gt.0.d0)then
           heat = heat + COVJ(LEV) * handc(LEV)
          else
           cool = cool + COVJ(LEV) * handc(LEV)
          endif
       GCO = GCO + COVJ(LEV)*ANDC(LEV)
   53 CONTINUE

          heat = heat * XK
          cool = cool * XK

C---RADIATIVE COOLING BY SiO (ERG CM-3 S-1)
c  GCO is the rate of energy loss by SiO to the radiation field
c  WWT is the net cooling rate, evaluated as (minus) the sum of
c  (positive) heating rate, HEAT, and the (negative) cooling 
c  rate, COOL 
c  GCO and WWT are both positive when cooling is taking place
c  and should be equal if the cooling is significant

          GCO = GCO * XK
          WWT = -(HEAT + COOL)

C  VERIFY THAT THE SUM OF THE SiO POPULATION DENSITIES (CM-3) 
C  IS EQUAL TO THE DENSITY OF THIS SPECIES
      verif = verif / XNCO
      if (abs(verif-1.0d0) .gt. 1.0d-2) then
        print *, '  verif SiO:', verif, verif-1.0d0
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

       end module msio

C---------------------------------------------------------------------------

