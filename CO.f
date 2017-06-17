C---------------------------------------------------------------------------
      module mco
      
C   THIS MODULE CALCULATES THE POPULATION DENSITIES OF THE 
C   ROTATIONAL LEVELS OF THE X VIBRATIONAL GROUND STATE OF CO
C   CALLED FROM SUBROUTINE DIFFUN

      public CO_LVG
C     public INITCO
C     public EVOLCO
      
      contains

       SUBROUTINE CO_LVG(ABH,ABH2,ABHE,ABCO,TN0,ROP,VGRAD,
     1                   NLEVCO,YTCO,ICNT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C---------------------------------------------------------------------
C      XNH  IS NUMBER DENSITY OF HYDROGEN ATOMS (CM-3)
C      XNH2 IS NUMBER DENSITY OF HYDROGEN MOLECULES (CM-3)
C      XNHE IS NUMBER DENSITY OF HELIUM ATOMS (CM-3)
C      TN0 IS TEMPERATURE OF NEUTRAL FLUID
C      ROP IS THE ORTHO:PARA H_2 RATIO
C      TBB IS TEMPERATURE OF RADIATION FIELD
C      WDIL IS DILUTION FACTOR OF RADIATION FIELD

      REAL*8    EVLTOT(NLEVCO,NLEVCO)
      REAL*8    YTCO(NLEVCO)

      COMMON /CHKIT/ TOTO

C BACKGROUND RADIATION TEMPERATURE AND DILUTION FACTOR
      DATA TBB, WDIL / 2.73D0, 1.D0 /
      
      DATA ICOUNT       / 0 /
      ICOUNT=ICOUNT+1
      
C VELOCITY GRADIENT VGRAD (S-1)

C INITIALIZE THE CO ROTATIONAL LEVEL POPULATION DENSITIES TO THEIR
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
      END SUBROUTINE CO_LVG
C
C      %%%%%%%%%%%%%%%%%
       SUBROUTINE INITCO(ABCO,TN0)
C      %%%%%%%%%%%%%%%%%

C  INITIALIZATION SUBROUTINE 

C--------+---------+---------+---------+---------+---------+---------+-*-------+

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C   NVXCO :  HIGHEST VIBRATIONAL LEVEL OF CO (X)
       PARAMETER (NVXCO=0)
C   NLEVCO :  NUMBER OF ROTATIONAL LEVELS OF CO 
       PARAMETER (NLEVCO=41)
! The value of the PARAMETER NLEVCO must be the same as that of the corresponding 
! PARAMETER NCO_lev in evolution.f90, mhd_vode.f90 and outputs.f90
       
C   NTEMP: NUMBER OF TEMPERATURES AT WHICH COLLISIONAL 
C   RATE COEFFICIENTS ARE TABULATED
       PARAMETER (NTEMP=25)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       CHARACTER*75       FICHIER
       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    AR(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)

      COMMON /FCTCO/ COVJ0(NLEVCO)
      COMMON /CO/ GCO,ELHM,NJLEV
      
      COMMON /TRANSF_CO/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE
      
      
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

C  ENERGIES, IN K
       DATA (ELHM(N), N=1,NLEVCO) /
     1             0.000, 5.53221, 16.59642, 33.19221,
     2       55.31894, 82.97576, 116.1616, 154.8753, 
     3       199.1152, 248.8797, 304.1670, 364.9748,  
     4       431.3009, 503.14272, 580.49752, 663.36183, 
     5      751.73353, 845.60870, 944.98378, 1049.8549, 
     6      1160.2182, 1276.05, 1397.38, 1524.19, 1656.47,
     7      1794.23, 1937.44, 2086.12, 2240.24, 2399.82, 
     8      2564.83, 2735.28, 2911.15, 3092.45, 3279.15,
     9      3471.27, 3668.78, 3871.69, 4079.98, 4293.64, 4512.67 /

       DATA(AR(N),N=1,NLEVCO)/ 0.0, 7.203e-08,6.910e-07,
     1      2.497e-06,6.126e-06,1.221e-05,2.137e-05,
     2      3.422e-05,5.134e-05,7.330e-05,1.006e-04,
     3      1.339e-04,1.735e-04,2.200e-04,2.739e-04,
     4      3.354e-04,4.050e-04,4.829e-04,5.695e-04,
     5      6.650e-04,7.695e-04,8.833e-04,1.006e-03,
     6      1.139e-03,1.281e-03,1.432e-03,1.592e-03,
     7      1.761e-03,1.940e-03,2.126e-03,2.321e-03,
     8      2.524e-03,2.735e-03,2.952e-03,3.175e-03,
     9      3.404e-03,3.638e-03,3.878e-03,4.120e-03,
     1      4.365e-03,4.613e-03 /


       DATA IREAD2 / 51 /

C--------+---------+---------+---------+---------+---------+---------+-*-------+


C  COLLISIONAL RATE COEFFICIENTS FOR H + CO, HE + CO AND H2 + CO

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
   10    CONTINUE
   11  CONTINUE

C  READ THE COLLISIONAL DE-EXCITATION RATE COEFFICIENTS (CM3 S-1)
       FICHIER = 'input/QHCO.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

        READ (IREAD2,*) (TTABH(IT),IT=1,NTEMPH)
       DO 17 ITRAN=1,136
        READ (IREAD2,*)LEVU,LEVL,(QHCO(IT,LEVU,LEVL),
     1        IT=1,NTEMPH) 
   17  CONTINUE
       CLOSE (IREAD2)
       
       FICHIER = 'input/QHECO.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

        READ (IREAD2,*) (TTABHE(IT),IT=1,NTMPHE)
       DO 18 ITRAN=1,105
        READ (IREAD2,*)LEVU,LEVL,(QHECO(IT,LEVU,LEVL),
     1        IT=1,NTMPHE) 
   18  CONTINUE
       CLOSE (IREAD2)
       
       FICHIER = 'input/QoH2CO.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

        READ (IREAD2,*) (TTAB(IT),IT=1,NTEMP)
       DO 19 ITRAN=1,820
        READ (IREAD2,*)INDEX,LEVU,LEVL,(QoH2CO(IT,LEVU,LEVL),
     1        IT=1,NTEMP) 
   19  CONTINUE
       CLOSE (IREAD2)

       FICHIER = 'input/QpH2CO.DAT              '
       OPEN(IREAD2, FILE = FICHIER, STATUS = 'OLD')

        READ (IREAD2,*) (TTAB(IT),IT=1,NTEMP)
       DO 20 ITRAN=1,820
        READ (IREAD2,*)INDEX,LEVU,LEVL,(QpH2CO(IT,LEVU,LEVL),
     1        IT=1,NTEMP) 
   20  CONTINUE
       CLOSE (IREAD2)

C*********************************************************************

C  RADIATIVE ELECTRIC DIPOLE TRANSITION PROBABILITIES OF CO


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

C   INITIAL POPULATIONS OF CO
C   HYPOTHESIS 1: LTE AT INITIAL TEMPERATURE

         ECH = 0.0D0
         DO 90 LEV=1,NLEVCO
           IJ = NJLEV(LEV)
           COVJ0(LEV) = (2.0D0*IJ+1.0D0) * DEXP(-ELHM(LEV)/TN0)
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

       PARAMETER (NTEMP=25)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    EVLCOL(NLEVCO,NLEVCO)
       REAL*8    EVLRAD(NLEVCO,NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)

      COMMON /TRANSF_CO/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE

       DATA MJCO / 40 /


       DATA (LEVEL(0,J), J=0,NLEVCO-1) /   1, 2, 3, 4, 5, 6, 7, 8, 9,
     1          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
     2          22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
     3          34, 35, 36, 37, 38, 39, 40, 41 /

       DATA (NJLEV(N), N=1,NLEVCO) /
     1       0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
     2      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
     3      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 /

       DATA (ELHM(N), N=1,NLEVCO) /
     1             0.000, 5.53221, 16.59642, 33.19221,
     2       55.31894, 82.97576, 116.1616, 154.8753, 
     3       199.1152, 248.8797, 304.1670, 364.9748,  
     4       431.3009, 503.14272, 580.49752, 663.36183, 
     5      751.73353, 845.60870, 944.98378, 1049.8549, 
     6      1160.2182, 1276.05, 1397.38, 1524.19, 1656.47,
     7      1794.23, 1937.44, 2086.12, 2240.24, 2399.82, 
     8      2564.83, 2735.28, 2911.15, 3092.45, 3279.15,
     9      3471.27, 3668.78, 3871.69, 4079.98, 4293.64, 4512.67 /



C------------------------------------------------------------------------------

C  CALCULATION OF THE ARRAY ELEMENTS NECESSARY TO DETERMINE THE EVOLUTION OF 
C  THE LEVEL POPULATIONS OF CO 

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

       PARAMETER (NTEMP=25)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    TTAB(NTEMP)
       REAL*8    YTCO(NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)


      COMMON /FCTCO/ FY(NLEVCO)

      COMMON /TRANSF_CO/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE

      COMMON /CHKIT/ TOTO

       DATA MJCO / 40 /


       DATA (LEVEL(0,J), J=0,NLEVCO-1) /   1, 2, 3, 4, 5, 6, 7, 8, 9,
     1          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
     2          22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
     3          34, 35, 36, 37, 38, 39, 40, 41 /

       DATA (NJLEV(N), N=1,NLEVCO) /
     1       0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
     2      15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
     3      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 /

       DATA (ELHM(N), N=1,NLEVCO) /
     1             0.000, 5.53221, 16.59642, 33.19221,
     2       55.31894, 82.97576, 116.1616, 154.8753, 
     3       199.1152, 248.8797, 304.1670, 364.9748,  
     4       431.3009, 503.14272, 580.49752, 663.36183, 
     5      751.73353, 845.60870, 944.98378, 1049.8549, 
     6      1160.2182, 1276.05, 1397.38, 1524.19, 1656.47,
     7      1794.23, 1937.44, 2086.12, 2240.24, 2399.82, 
     8      2564.83, 2735.28, 2911.15, 3092.45, 3279.15,
     9      3471.27, 3668.78, 3871.69, 4079.98, 4293.64, 4512.67 /


C------------------------------------------------------------------------------

C     WRITE(6, 2150)(FY(J),J=1,NLEVCO)
C2150 FORMAT(6D12.4)

C  CALCULATE THE CONTRIBUTIONS TO THE EVOLUTION OF THE LEVEL POPULATIONS OF CO (CM-3 S-1)

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
C      WDIL IS DILUTION FACTOR OF RADIATION FIELD
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

       PARAMETER (NTEMP=25)
       PARAMETER (NTEMPH=20)
       PARAMETER (NTMPHE=10)

       REAL*8    ELHM(NLEVCO)
       REAL*8    QHCO(NTEMPH,NLEVCO,NLEVCO)
       REAL*8    QHECO(NTMPHE,NLEVCO,NLEVCO)
       REAL*8    QoH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    QpH2CO(NTEMP,NLEVCO,NLEVCO)
       REAL*8    EVLDI(NLEVCO,NLEVCO)
       REAL*8    EVLTOT(NLEVCO,NLEVCO)
       REAL*8    TAU(NLEVCO,NLEVCO)
       REAL*8    EP(NLEVCO,NLEVCO)
       REAL*8    AEITQ(NLEVCO)
       REAL*8    CEITQ(NLEVCO)
       REAL*8    TTAB(NTEMP)
       REAL*8    TTABH(NTEMPH)
       REAL*8    TTABHE(NTMPHE)
       REAL*8    HANDC(NLEVCO)
       REAL*8    ANDC(NLEVCO)
       REAL*8    TVGRAD(NLEVCO)
       REAL*8    EVLCOL(NLEVCO,NLEVCO)
       REAL*8    EVLRAD(NLEVCO,NLEVCO)
       INTEGER   NJLEV(NLEVCO)
       INTEGER   LEVEL(0:NVXCO,0:NLEVCO-1)

      COMMON /FCTCO/ COVJ(NLEVCO),TAU,WWT,TVGRAD
      COMMON /CO/ GCO,ELHM,NJLEV
      
      COMMON /TRANSF_CO/QHCO,QHECO,QoH2CO,QpH2CO,EVLDI,AEITQ,TTAB,
     1               TTABH,TTABHE
      

      DATA XK,PI/1.38062D-16,3.1415926536D0/

      DATA MJCO / 40 /


       DATA (LEVEL(0,J), J=0,NLEVCO-1) /   1, 2, 3, 4, 5, 6, 7, 8, 9,
     1          10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
     2          22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
     3          34, 35, 36, 37, 38, 39, 40, 41 /



C     PI=ACOS(-1.D0)
C

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

C  CO - H (Balakrishnan et al.2002, ApJ, 568, 443)

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

C  CO - HE (Cecchi-Pestellini et al. 2002, ApJ, 571, 1015)

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

C  CO - o-H2 (Yang et al. 2010, ApJ, 718, 1062)

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
    
C  CO - p-H2 (Yang et al. 2010, ApJ, 718, 1062)

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
     2       + CCOHE * XNHE
C    3       + CCOH * XNH
C  ADOPT THE RATE COEFFICIENTS OF o-H2 FOR COLLISIONS WITH H
     3       + CCOoH2 * XNH
     
C     print *, TN0,LEVU,LEVL,EVLDI(LEVU,LEVL),CCOoH2,CCOpH2,CCOHE,CCOH

C     EXCITATION RATES: NEUTRALS
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
C  TLINE IS NON-ZERO ONLY FOR LEVL = LEVU - 1
       TLINE = EP(LEVU,LEVL)/(8.D0*PI*VGRAD*WN3)*
     1         EVLDI(LEVU,LEVL)*COVJ(LEVU)
       TLINE = TLINE - (1.D0-EXP(-TS))*FACBB
       TLINE = TLINE*(ELHM(LEVU) - ELHM(LEVL))
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
C  COLLISIONAL AND RADIATIVE RATES (S-1)
             CEITQ(LEVU) = CEITQ(LEVU) + EVLCOL(LEVU,LEVL) 
     1                                 + EVLRAD(LEVU,LEVL)
C  HEATING AND COOLING RATES (K S-1)
             HANDC(LEVU) = HANDC(LEVU) + EVLCOL(LEVU,LEVL) 
     1                                 * (ELHM(LEVU) - ELHM(LEVL))
             ANDC(LEVU) = ANDC(LEVU) + (EVLDI(LEVU,LEVL) * EP(LEVU,LEVL)
     1            + EVLRAD(LEVU,LEVL)) * (ELHM(LEVU) - ELHM(LEVL))
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
c  to the density of CO
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

C---RADIATIVE COOLING BY CO (ERG CM-3 S-1)
c  GCO is the net rate of energy loss by CO to the radiation field
c  WWT is the net cooling rate, evaluated as (minus) the sum of
c  (positive) heating rate, HEAT, and the (negative) cooling 
c  rate, COOL 
c  GCO and WWT are both positive when cooling is taking place
c  and should be equal if the cooling is significant

          GCO = GCO * XK
          WWT = -(HEAT + COOL)

C  VERIFY THAT THE SUM OF THE CO POPULATION DENSITIES (CM-3) 
C  IS EQUAL TO THE DENSITY OF THIS SPECIES
      verif = verif / XNCO
      if (abs(verif-1.0d0) .gt. 1.0d-2) then
        print *, '  verif CO:', verif, verif-1.0d0
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

       end module mco

C---------------------------------------------------------------------------

