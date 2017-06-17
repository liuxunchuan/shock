 MODULE MODULE_LINE_EXCIT
  !***************************************************************************
  ! Should replace part or all MODULE_FINE_STRUCTURE
  ! Sylvie Cabrit & Jacques Le Bourlot - Sept. 2002
  ! Shorter version SC & GPdF - Nov. 2002
  !***************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

! Atomic or molecular parameters:
!
!  All variable names use the same scheme
!  first 3 caracters : denote which parameter
!                nlv : number of levels
!                gst : statistical weight (usually, g = 2J + 1)
!                elk : energy in Kelvin

!                ntr : number of radiative transitions
!                aij : Aij Einstein coefficient (s-1)
!                wlk : transition energy (in K)
!                emi : emissivity of transition (erg cm-3 s-1)
!                iup : index of upper level
!                jlo : index of lower level

!  next 3 caracters   : denote species
!                cat  : atomic carbon    (C I)
!                nat  : atomic nitrogen  (N I)
!                oat  : atomic oxygen    (O I)
!                sat  : atomic sulphur   (S I)
!                siat : atomic silicon   (Si I)
!                cpl  : ionized carbon   (C II)
!                npl  : ionized nitrogen (N II)
!                opl  : ionized oxygen   (O II)
!                spl  : ionized sulphur   (S II)
!                sipl : ionized silicon  (Si II)
!                fepl : ionized iron  (Fe II)

!  real (kind=dp), public, parameter :: hcsurk = 1.43883442658354_dp    ! Converts from Wave number (cm-1) to temperature (K)

! C
! Data from NIST (12 Sept 2002)

! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

  INTEGER,        PUBLIC, PARAMETER          :: nlvcat = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvcat) :: gstcat
  DATA gstcat                                 / 1.0_DP, 3.0_DP, 5.0_DP, 5.0_DP, 1.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvcat) :: elkcat
  DATA elkcat                                 / 0.0_DP, 23.5969_DP, 62.4454_DP, 14665.507_DP, 31147.902_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvcat)        :: pop_cat   ! Atomic Carbon populations

! Transitions
!       1: 2-1 -> 609.75 micron
!       2: 3-2 -> 370.37 micron
!       3: 3-1 -> 230.41 micron
!       4: 4-3 -> 9850.26 Angstrom
!       5: 4-2 -> 9824.13 Angstrom
!       6: 4-1 -> 9808.32 Angstrom
!       7: 5-4 -> 8727.13 Angstrom
!       8: 5-3 -> 4627.346 Angstrom
!       9: 5-2 -> 4621.570 Angstrom
!      10: 5-1 

  INTEGER,        PUBLIC, PARAMETER          :: ntrcat = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntrcat) :: wlkcat
  DATA wlkcat                                 / 23.5969_DP, 38.8485_DP, 62.4454_DP, &
                                                14603.0616_DP, 14641.9101_DP, 14665.507_DP, &
                                                16482.395_DP, 31085.4566_DP, 31124.3051_DP, 31147.902_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrcat) :: aijcat
  DATA aijcat                                 / 7.93e-8_DP, 2.65e-7_DP, 1.71e-14_DP, &
                                                2.44e-4_DP, 8.21e-5_DP, 7.77e-8_DP, &
                                                5.28e-1_DP, 2.00e-5_DP, 2.71e-3_DP, 0.0e0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntrcat) :: emicat, emicat_o, intcat
  INTEGER,        PUBLIC, DIMENSION (ntrcat) :: iupcat
  DATA iupcat                                 / 2, 3, 3, 4, 4, 4, 5, 5, 5, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntrcat) :: jlocat
  DATA jlocat                                 / 1, 2, 1, 3, 2, 1, 4, 3, 2, 1 /


! N
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 4S J = 3/2
!        2 - 2D J = 5/2
!        3 - 2D J = 3/2
!        4 - 2P J = 1/2
!        5 - 2P J = 3/2

  INTEGER,        PUBLIC, PARAMETER          :: nlvnat = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvnat) :: gstnat
  DATA gstnat                                 / 4.0_DP, 6.0_DP, 4.0_DP, 2.0_DP, 4.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvnat) :: elknat
  DATA elknat                                 / 0.0_DP, 27660.82_DP, 27673.36_DP, 41494.431_DP, 41494.986_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvnat)        :: pop_nat   ! Atomic Nitrogen populations

! Transitions
!       1: 3-2 -> 1.147 mm (261.21 GHz)
!       2: 4-3 -> 1.04104 micron
!       3: 5-3 -> 1.04100 micron
!       4: 4-2 -> 1.04010 micron
!       5: 5-2 -> 1.04006 micron
!       6: 2-1 -> 5200.257 Angstrom
!       7: 3-1 -> 5197.902 Angstrom
!       8: 4-1 -> 3466.543 Angstrom
!       9: 5-1 -> 3466.497 Angstrom
!      10: 5-4

  INTEGER,        PUBLIC, PARAMETER          :: ntrnat = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntrnat) :: wlknat
  DATA wlknat                                 / 12.54_DP, 13821.071_DP, 13821.626_DP, &
                                                13833.611_DP, 13834.166_DP, 27660.82_DP, &
                                                27673.36_DP, 41494.431_DP, 41494.986_DP, 0.555_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrnat) :: aijnat
  DATA aijnat                                 / 1.27e-8_DP, 5.29e-2_DP, 2.76e-2_DP, &
                                                3.45e-2_DP, 6.14e-2_DP, 7.27e-6_DP, &
                                                2.02e-5_DP, 2.71e-3_DP, 6.58e-3_DP, 0.0e0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntrnat) :: eminat, eminat_o, intnat
  INTEGER,        PUBLIC, DIMENSION (ntrnat) :: iupnat
  DATA iupnat                                 / 3, 4, 5, 4, 5, 2, 3, 4, 5, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntrnat) :: jlonat
  DATA jlonat                                 / 2, 3, 3, 2, 2, 1, 1, 1, 1, 4 /

! O
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 3P J = 2
!        2 - 3P J = 1
!        3 - 3P J = 0
!        4 - 1D J = 2
!        5 - 1S J = 0

  INTEGER,        PUBLIC, PARAMETER          :: nlvoat = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvoat) :: gstoat
  DATA gstoat                                 / 5.0_DP, 3.0_DP, 1.0_DP, 5.0_DP, 1.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvoat) :: elkoat
  DATA elkoat                                 / 0.0_DP, 227.717_DP, 326.582_DP, 22831.226_DP, 48621.932_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvoat)        :: pop_oat   ! Atomic Oxygen populations

! Transitions
!       1: 3-2 -> 145.53 micron
!       2: 2-1 -> 63.19  micron
!       3: 3-1 -> 44.06  micron
!       4: 4-3 -> 6391.73 Angstrom
!       5: 4-2 -> 6363.78 Angstrom
!       6: 4-1 -> 6300.30 Angstrom
!       7: 5-4 -> 5577.34 Angstrom
!       8: 5-2 -> 2972.29 Angstrom
!       9: 5-1 -> 2958.37 Angstrom
!      10: 5-3

  INTEGER,        PUBLIC, PARAMETER          :: ntroat = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntroat) :: wlkoat
  DATA wlkoat                                 / 98.865_DP, 227.717_DP, 326.582_DP, &
                                                22504.644_DP, 22603.509_DP, 22831.226_DP, &
                                                25790.706_DP, 48394.215_DP, 48621.932_DP, 48295.35_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntroat) :: aijoat
  DATA aijoat                                 / 1.66e-5_dp, 8.46e-5_dp, 1.10e-10_dp, &
                                                1.20e-6_dp, 2.17e-3_dp, 6.71e-3_dp, &
                                                1.02e-0_dp, 7.83e-2_dp, 2.87e-4_dp, 0.0e0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntroat) :: emioat, emioat_o, intoat
  INTEGER,        PUBLIC, DIMENSION (ntroat) :: iupoat
  DATA iupoat                                 / 3, 2, 3, 4, 4, 4, 5, 5, 5, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntroat) :: jlooat
  DATA jlooat                                 / 2, 1, 1, 3, 2, 1, 4, 2, 1, 3 /

! S
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 3P J = 2
!        2 - 3P J = 1
!        3 - 3P J = 0
!        4 - 1D J = 2
!        5 - 1S J = 0

  INTEGER,        PUBLIC, PARAMETER          :: nlvsat = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsat) :: gstsat
  DATA gstsat                                 / 5.0_DP, 3.0_DP, 1.0_DP, 5.0_DP, 1.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsat) :: elksat
  DATA elksat                                 / 0.0_DP, 569.86_DP, 825.37_DP, 13292.83_DP, 31913.28_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsat)        :: pop_sat   ! Atomic Sulphur populations

! Transitions
!       1: 3-2 -> 56.31 micron
!       2: 2-1 -> 25.25 micron
!       3: 3-1 -> 17.43 micron
!       4: 4-3 -> 1.15407 micron
!       5: 4-2 -> 1.13089 micron
!       6: 4-1 -> 1.08241 micron
!       7: 5-4 -> 7725.05 Angstrom
!       8: 5-2 -> 4589.26 Angstrom
!       9: 5-1 -> 4507.31 Angstrom
!      10: 5-3 

  INTEGER,        PUBLIC, PARAMETER          :: ntrsat = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsat) :: wlksat
  DATA wlksat                                 / 255.51_DP, 569.86_DP, 825.37_DP, &
                                                12467.46_DP, 12722.97_DP, 13292.83_DP, &
                                                18620.45_DP, 31343.42_DP, 31913.28_DP, 31087.91_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsat) :: aijsat
  DATA aijsat                                 / 3.02e-4_DP, 1.39e-3_DP, 6.71e-8_DP, &
                                                3.84e-6_DP, 8.16e-3_DP, 2.78e-2_DP, &
                                                1.53e0_DP, 3.50e-1_DP, 8.23e-3_DP, 0.0e0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsat) :: emisat, emisat_o, intsat
  INTEGER,        PUBLIC, DIMENSION (ntrsat) :: iupsat
  DATA iupsat                                 / 3, 2, 3, 4, 4, 4, 5, 5, 5, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntrsat) :: jlosat
  DATA jlosat                                 / 2, 1, 1, 3, 2, 1, 4, 2, 1, 3 /

! Si
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

  INTEGER,        PUBLIC, PARAMETER          :: nlvsiat = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsiat) :: gstsiat
  DATA gstsiat                                 / 1.0_DP, 3.0_DP, 5.0_DP, 5.0_DP, 1.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsiat) :: elksiat
  DATA elksiat                                 / 0.0_DP, 110.951_DP, 321.086_DP, 9062.998_DP, 22149.939_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsiat)         :: pop_siat   ! Atomic Silicon populations

! Transitions
!       1: 2-1 -> 129.682 micron
!       2: 3-2 -> 68.472 micron
!       3: 3-1 -> 44.812 micron
!       4: 4-3 -> 1.6459 micron
!       5: 4-2 -> 1.6073 micron
!       6: 4-1 -> 1.5876 micron
!       7: 5-4 -> 1.0994 micron
!       8: 5-3 -> 6589.61 Angstrom
!       9: 5-2 -> 6526.78 Angstrom
!      10: 5-1

  INTEGER,        PUBLIC, PARAMETER          :: ntrsiat = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsiat) :: wlksiat
  DATA wlksiat                                 / 110.951_DP, 210.135_DP, 321.086_DP, &
                                                 8741.912_DP, 8952.047_DP, 9062.998_DP, &
                                                 13086.941_DP, 21828.853_DP, 22038.988_DP, 21828.853_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsiat) :: aijsiat
  DATA aijsiat                                 / 8.25e-6_DP, 4.21e-5_DP, 3.56e-13_DP, &
                                                 2.25e-3_DP, 7.93e-4_DP, 4.70e-7_DP, &
                                                 1.14e-0_DP, 9.02e-4_DP, 3.13e-2_DP, 0.0e0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsiat) :: emisiat, emisiat_o, intsiat
  INTEGER,        PUBLIC, DIMENSION (ntrsiat) :: iupsiat
  DATA iupsiat                                 / 2, 3, 3, 4, 4, 4, 5, 5, 5, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntrsiat) :: jlosiat
  DATA jlosiat                                 / 1, 2, 1, 3, 2, 1, 4, 3, 2, 1 /

! C+
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 2P J = 1/2
!        2 - 2P J = 3/2
!        3 - 4P J = 1/2
!        4 - 4P J = 3/2
!        5 - 4P J = 5/2

  INTEGER,        PUBLIC, PARAMETER          :: nlvcpl = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvcpl) :: gstcpl
  DATA gstcpl                                 / 2.0_DP, 4.0_DP, 2.0_DP, 4.0_DP, 6.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvcpl) :: elkcpl
  DATA elkcpl                                 / 0.0_DP, 91.25_DP, 61874.63_DP, 61906.28_DP, 61947.00_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvcpl)        :: pop_cpl   ! Ionized Carbon populations

! Transitions
!       1: 2-1 -> 157.68 micron
!       2: 3-2 -> 2328.12 Angstrom
!       3: 4-2 -> 2326.93 Angstrom
!       4: 5-2 -> 2325.40 Angstrom
!       5: 3-1 -> 2324.69 Angstrom
!       6: 4-1 -> 2323.50 Angstrom
!	    7: 5-1  

  INTEGER,        PUBLIC, PARAMETER          :: ntrcpl = 7
  REAL (kind=DP), PUBLIC, DIMENSION (ntrcpl) :: wlkcpl
  DATA wlkcpl                                 / 91.25_DP, 61783.38_DP, 61815.03_DP, &
                                                61855.75_DP, 61874.63_DP, 61906.28_DP, 61947.00_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrcpl) :: aijcpl
  DATA aijcpl                                 / 2.29e-6_DP, 6.55e+1_DP, 5.24e+0_DP, &
                                                4.32e+1_DP, 5.53e+1_DP, 1.71e-0_DP, 0.0e0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntrcpl) :: emicpl, emicpl_o, intcpl
  INTEGER,        PUBLIC, DIMENSION (ntrcpl) :: iupcpl
  DATA iupcpl                                 / 2, 3, 4, 5, 3, 4, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntrcpl) :: jlocpl
  DATA jlocpl                                 / 1, 2, 2, 2, 1, 1, 1 /

! N+
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

  INTEGER,        PUBLIC, PARAMETER          :: nlvnpl = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvnpl) :: gstnpl
  DATA gstnpl                                 / 1.0_DP, 3.0_DP, 5.0_DP, 5.0_DP, 1.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvnpl) :: elknpl
  DATA elknpl                                 / 0.0_DP, 70.07_DP, 188.20_DP, 22037.48_DP, 47033.77_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvnpl)        :: pop_npl   ! Ionized Nitrogen populations

! Transitions
!       1: 2-1 -> 205.34 micron
!       2: 3-2 -> 121.80 micron
!       3: 3-1 -> 76.45 micron
!       4: 4-3 -> 6583.45 Angstrom
!       5: 4-2 -> 6548.05 Angstrom
!       6: 4-1 -> 6527.23 Angstrom
!       7: 5-4 -> 5754.59 Angstrom
!       8: 5-3 -> 3070.55 Angstrom
!       9: 5-2 -> 3062.83 Angstrom
!      10: 5-1 

  INTEGER,        PUBLIC, PARAMETER          :: ntrnpl = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntrnpl) :: wlknpl
  DATA wlknpl                                 / 70.07_DP, 118.13_DP, 188.20_DP, &
                                                21849.28_DP, 21967.41_DP, 22037.48_DP, &
                                                24996.29_DP, 46845.57_DP, 46963.7_DP, 47033.77_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrnpl) :: aijnpl
  DATA aijnpl                                 / 2.08e-6_DP, 7.46e-6_DP, 1.16e-12_DP, &
                                                2.99e-3_DP, 1.01e-3_DP, 5.35e-7_DP, &
                                                1.12e-0_DP, 1.51e-4_DP, 3.38e-2_DP, 0.0e0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntrnpl) :: eminpl, eminpl_o, intnpl
  INTEGER,        PUBLIC, DIMENSION (ntrnpl) :: iupnpl
  DATA iupnpl                                 / 2, 3, 3, 4, 4, 4, 5, 5, 5, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntrnpl) :: jlonpl
  DATA jlonpl                                 / 1, 2, 1, 3, 2, 1, 4, 3, 2, 1 /

! O+
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 4S J = 3/2
!        2 - 2D J = 5/2
!        3 - 2D J = 3/2
!        4 - 2P J = 3/2
!        5 - 2P J = 1/2

  INTEGER,        PUBLIC, PARAMETER          :: nlvopl = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvopl) :: gstopl
  DATA gstopl                                 / 4.0_DP, 6.0_DP, 4.0_DP, 4.0_DP, 2.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvopl) :: elkopl
  DATA elkopl                                 / 0.0_DP, 38575.94_DP, 38604.75_DP, 58226.77_DP, 58229.63_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvopl)        :: pop_opl   ! Ionized Oxygen populations

! Transitions
!       1: 5-4 -> 5.03 mm (59.66 GHz)
!       2: 3-2 -> 0.4995 mm (600.18 GHz)
!       3: 4-3 -> 7330.73 Angstrom
!       4: 5-3 -> 7329.67 Angstrom
!       5: 4-2 -> 7319.99 Angstrom
!       6: 5-2 -> 7318.92 Angstrom
!       7: 2-1 -> 3728.815 Angstrom
!       8: 3-1 -> 3726.032 Angstrom
!       9: 4-1 -> 2470.341 Angstrom
!      10: 5-1 -> 2470.219 Angstrom

  INTEGER,        PUBLIC, PARAMETER          :: ntropl = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntropl) :: wlkopl
  DATA wlkopl                                 / 2.86_DP, 28.81_DP, 19622.02_DP, &
                                                19624.88_DP, 19650.83_DP, 19653.69_DP, &
                                                38575.94_DP, 38604.75_DP, 58226.77_DP, &
                                                58229.63_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntropl) :: aijopl
  DATA aijopl                                 / 2.08e-11_DP, 1.20e-7_DP, 6.14e-2_DP, &
                                                1.02e-1_DP, 1.17e-1_DP, 6.15e-2_DP, &
                                                3.82e-5_DP, 1.65e-4_DP, 5.64e-2_DP, &
                                                2.32e-2_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntropl) :: emiopl, emiopl_o, intopl
  INTEGER,        PUBLIC, DIMENSION (ntropl) :: iupopl
  DATA iupopl                                 / 5, 3, 4, 5, 4, 5, 2, 3, 4, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntropl) :: jloopl
  DATA jloopl                                 / 4, 2, 3, 3, 2, 2, 1, 1, 1, 1 /

! S+
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 4S J = 3/2
!        2 - 2D J = 3/2
!        3 - 2D J = 5/2
!        4 - 2P J = 1/2
!        5 - 2P J = 3/2

  INTEGER,        PUBLIC, PARAMETER          :: nlvspl = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvspl) :: gstspl
  DATA gstspl                                 / 4.0_DP, 4.0_DP, 6.0_DP, 2.0_DP, 4.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvspl) :: elkspl
  DATA elkspl                                 / 0.0_DP, 21370.92_DP, 21416.66_DP, 35287.17_DP, 35354.38_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvspl)        :: pop_spl   ! Ionized Sulphur populations

! Transitions
!       1: 3-2 -> 0.31456 mm (953.04 GHz)
!       2: 5-4 -> 0.21409 mm (1400.33 GHz)
!       3: 4-3 -> 1.0373 micron
!       4: 4-2 -> 1.0339 micron
!       5: 5-3 -> 1.0323 micron
!       6: 5-2 -> 1.0290 micron
!       7: 2-1 -> 6730.82 Angstrom
!       8: 3-1 -> 6716.44 Angstrom
!       9: 4-1 -> 4076.35 Angstrom
!      10: 5-1 -> 4068.60 Angstrom

  INTEGER,        PUBLIC, PARAMETER          :: ntrspl = 10
  REAL (kind=DP), PUBLIC, DIMENSION (ntrspl) :: wlkspl
  DATA wlkspl                                 / 45.74_DP, 67.21_DP, 13870.51_DP, &
                                                13916.25_DP, 13937.72_DP, 13983.46_DP, &
                                                21370.92_DP, 21416.66_DP, 35287.17_DP, &
                                                35354.38_DP /

! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrspl) :: aijspl
  DATA aijspl                                 / 3.35e-7_DP, 1.03e-6_DP, 7.79e-2, &
                                                1.63e-1_DP, 1.79e-1_DP, 1.33e-1_DP, &
                                                8.82e-4_DP, 2.60e-4_DP, 9.06e-2_DP, &
                                                2.25e-1_DP /

  REAL (kind=DP), PUBLIC, DIMENSION (ntrspl) :: emispl, emispl_o, intspl
  INTEGER,        PUBLIC, DIMENSION (ntrspl) :: iupspl
  DATA iupspl                                 / 3, 5, 4, 4, 5, 5, 2, 3, 4, 5 /
  INTEGER,        PUBLIC, DIMENSION (ntrspl) :: jlospl
  DATA jlospl                                 / 2, 4, 3, 2, 3, 2, 1, 1, 1, 1 /

! Si+
! Data from NIST (12 Sept. 2002)

! Energy levels.
!        1 - 2P J = 1/2
!        2 - 2P J = 3/2
!        3 - 4P J = 1/2
!        4 - 4P J = 3/2
!        5 - 4P J = 5/2

  INTEGER,        PUBLIC, PARAMETER          :: nlvsipl = 5
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsipl) :: gstsipl
  DATA gstsipl                                 / 2.0_DP, 4.0_DP, 2.0_DP, 4.0_DP, 6.0_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsipl) :: elksipl
  DATA elksipl                                 / 0.0_DP, 413.29_DP, 61617.06_DP, 61772.93_DP, 62025.14_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (nlvsipl)         :: pop_sipl   ! Ionized Silicon populations

! Transitions
!       1: 2-1 -> 34.814 micron
!       2: 3-2 -> 2350.17 Angstrom
!       3: 4-2 -> 2344.20 Angstrom
!       4: 5-2 -> 2334.60 Angstrom
!       5: 3-1 -> 2334.40 Angstrom
!       6: 4-1 -> 2328.52 Angstrom

  INTEGER,        PUBLIC, PARAMETER          :: ntrsipl = 6
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsipl) :: wlksipl
  DATA wlksipl                                 / 413.29_DP, 61203.77_DP, 61359.64_DP, &
                                                 61611.85_DP, 61617.06_DP, 61772.93_DP /
! Revision DRF - Nov. 2002
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsipl) :: aijsipl
  DATA aijsipl                                 / 2.2e-4_DP, 4.9e+3_DP, 1.7e+3_DP, &
                                                2.7e+3_DP, 6.3e+3_DP, 2.0e+1_DP /
  REAL (kind=DP), PUBLIC, DIMENSION (ntrsipl) :: emisipl, emisipl_o, intsipl
  INTEGER,        PUBLIC, DIMENSION (ntrsipl) :: iupsipl
  DATA iupsipl                                 / 2, 3, 4, 5, 3, 4 /
  INTEGER,        PUBLIC, DIMENSION (ntrsipl) :: jlosipl
  DATA jlosipl                                 / 1, 2, 2, 2, 1, 1 /

! Fe+
   INTEGER (KIND=LONG), PUBLIC, PARAMETER       :: nlvfepl=35
   REAL (kind=DP), PUBLIC, DIMENSION (nlvfepl)  :: gstfepl
   DATA gstfepl                  / 10.0_DP, 8.0_DP, 6.0_DP, 4.0_DP, 2.0_DP, & 
                                   10.0_DP, 8.0_DP, 6.0_DP, 4.0_DP, 8.0_DP, &
                                   6.0_DP, 4.0_DP, 2.0_DP, 6.0_DP, 4.0_DP,  & 
                                   2.0_DP, 6.0_DP, 4.0_DP, 2.0_DP, 14.0_DP, &
                                   12.0_DP, 10.0_DP, 8.0_DP, 10.0_DP, 8.0_DP, &
                                   6.0_DP, 4.0_DP, 12.0_DP, 10.0_DP, 8.0_DP, &
                                   6.0_DP, 8.0_DP, 6.0_DP, 4.0_DP, 2.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvfepl)  :: elkfepl
   DATA elkfepl                       / 0.0_DP, 553.95083_DP, 959.550155_DP, &
             1240.471089_DP, 1404.604668_DP, 2692.421982_DP, 3494.151388_DP, &
           4081.244575_DP, 4483.687486_DP, 11440.426112_DP, 12068.552694_DP, & 
          12483.62126_DP, 12723.508803_DP, 19378.809799_DP, 19664.465355_DP, &
          19997.46714_DP, 29957.534620_DP, 31370.030326_DP, 32228.575202_DP, & 
          30563.50_DP, 30820.75_DP, 31038.54_DP, 31224.77_DP, 32556.77_DP, &
          32804.55_DP, 32990.78_DP, 33123.35_DP, 36570.14_DP, 37113.05_DP, &
          37365.56_DP, 37471.30_DP, 45278.68_DP, 45141.37_DP, 45106.65_DP, &
          45112.96_DP /

   REAL (kind=DP), PUBLIC, DIMENSION (nlvfepl)  :: pop_fepl

! Observed transitions (microns).    
!  1: 12-3 -> 1.248	12: 13-8 -> 1.664
!  2: 10-1 -> 1.257	13: 11-7 -> 1.677
!  3: 13-5 -> 1.271	14: 12-8 -> 1.711
!  4: 12-4 -> 1.279	15: 13-9 -> 1.745
!  5: 11-3 -> 1.295	16: 12-9 -> 1.798
!  6: 12-5 -> 1.298	17: 11-8 -> 1.800
!  7: 10-2 -> 1.321	18: 10-7 -> 1.810
!  8: 11-4 -> 1.328	19: 7-6  -> 17.936
!  9: 11-6 -> 1.534	20: 2-1  -> 25.988
! 10: 12-7 -> 1.600	21: 9-8  -> 35.777
! 11: 10-6 -> 1.644

    INTEGER (KIND=LONG), PUBLIC, PARAMETER      :: ntrfepl=256
    REAL (KIND=DP), PUBLIC, DIMENSION(ntrfepl)  :: wlkfepl
    REAL(KIND=DP),PUBLIC,DIMENSION(ntrfepl)     :: emifepl,emifepl_o,intfepl
!
! emissivity (erg/cm3/s), emissivity at last call to DRIVE, integrated intensity (erg/cm2/s/sr)

! End of public database on levels and transitions
! ================================================

! Terms coupling fluids. Give global energy transfer by collisions
!          cool_X < 0 if the fluid X loses energy

  REAL (kind=DP), PUBLIC  :: Cool_n, Cool_i, Cool_e  ! erg cm-3 s-1

! This is the minimum collision rate (used when nothing else is known)

  REAL (kind=DP), PUBLIC, PARAMETER                   :: colr_min = 1.0e-30_DP

! changes made in Nov.2002:

! * Collision matrix defined only inside of the subroutine that calculates
! the level populations
! => Avoids adding a new set of variables for each new atom/ion.
! => routine can easily be cut and pasted for each new atom
! * Teff's are declared here instead of inside each atomic routine
! * alpha's and Teff's will be calculated automatically (condition on charge_cool)
!
! These are private variables, used in LINE_THERMAL_BALANCE and THERMAL_LOSS

  REAL (kind=DP) :: alp_H, alp_H2, alp_He, alp_Hp
  REAL (kind=DP) :: Teff_H, Teff_H2, Teff_He, Teff_Hp, Teff_e
  REAL (kind=DP) :: Dens_cool, mass_cool, charge_cool
!  INTEGER (KIND=LONG), PRIVATE :: i

CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLCAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4                           ! auxiliary factor
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

! ec matrices will be multiplied by collider density in THERMAL_LOSS
! units: cm3 s-1

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do


! Build collision matrix for atomic carbon
! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

! collisions with H
! Fine structure: reference unknown

  ec_H(1,2) = 1.01e-10_DP * Teff_H**0.117_DP
  ec_H(1,3) = 4.49e-11_DP * Teff_H**0.194_DP
  ec_H(2,3) = 1.06e-10_DP * Teff_H**0.234_DP

! 1D - 3P
! DRF : do not use O I + H rates

! collisions with H2
! fit from Schroder, Staemmler, Smith, Flower & Jaquet (1991)

  ec_paraH2(1,2)  = 0.80e-10_DP
  ec_orthoH2(1,2) = 0.75e-10_DP
  ec_paraH2(1,3)  = 0.90e-10_DP
  ec_orthoH2(1,3) = 3.54e-11_DP * Teff_H2**0.167_DP
  ec_paraH2(2,3)  = 2.00e-10_DP
  ec_orthoH2(2,3) = 5.25e-11_DP * Teff_H2**0.244_DP

! collisions with He
! Fine structure: Staemmler & Flower (1991)

  ec_He(1,2) = 8.38e-12_DP * Teff_He**0.159_DP
  ec_He(1,3) = 5.98e-11_DP * Teff_He**0.078_DP
  ec_He(2,3) = 3.68e-11_DP * Teff_He**0.041_DP

! collisions avec H+
! from A&A, 236, 515 (1990)  Roueff & Le Bourlot (fit: DRF)

  ec_Hp(1,2) = 10._DP**(-10.359_DP + 0.7959_DP*log10(Teff_Hp) - 0.08748_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(1,3) = 10._DP**(-13.232_DP + 2.4171_DP*log10(Teff_Hp) - 0.29151_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(2,3) = 10._DP**(-11.290_DP + 1.7915_DP*log10(Teff_Hp) - 0.23010_DP*(log10(Teff_Hp))**2._DP)

!--- effective temperature with electrons : Te

! collisions with electrons
! Should we use Pequignot 1990...
! From Pequignot & Aldrovandi, A&A 50, 141, 1976; as compiled by Mendoza (1983)

  auxt4 = 1.83e-4_DP * Te ** 0.444_DP
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstcat(4) * gstcat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstcat(4) * gstcat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstcat(4) * gstcat(3) / 9.0_DP
  auxt4 = 9.86e-5_DP * Te ** 0.343_DP
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstcat(5) * gstcat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstcat(5) * gstcat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstcat(5) * gstcat(3) / 9.0_DP
  auxt4 = 2.77e-3_DP
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstcat(5)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  end subroutine COLCAT


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLOAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2
  REAL (kind=DP)                      :: toto


  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic oxygen
! Energy levels.
!        1 - 3P J = 2
!        2 - 3P J = 1
!        3 - 3P J = 0
!        4 - 1D J = 2
!        5 - 1S J = 0


! All collision rates should be checked against Pequignot 1990

! collisions with H

! Fine structure

!  ec_H(1,2) = 4.37e-12_dp * Teff_H**0.66_dp
!  ec_H(1,3) = 1.06e-12_dp * Teff_H**0.80_dp
!  ec_H(2,3) = 1.35e-11_dp * Teff_H**0.45_dp

!  fit from Abrahamsson et al. 2007, ApJ, 654, 1171 (done by Evelyne - XII 08)
  ec_H(1,2) = 5.9601e-11_dp * Teff_H**0.39561_dp
  ec_H(1,3) = 6.1719e-11_dp * Teff_H**0.36291_dp
  IF (Teff_H <= 1000.0_dp) THEN
     ec_H(2,3) = (2.5544e-11_dp - 1.232e-14_dp * Teff_H) * Teff_H**0.62996_dp
  ELSE
     ec_H(2,3) =  1.063e-9_dp
  ENDIF

! 1D - 3P
! fit by JLB, from data of Federman & Shipsey, 1983, Ap J, 269, 791 (using potential V2)

  ec_H(1,4) = 9.85e-14_DP * Teff_H**0.245_DP
  ec_H(2,4) = 9.85e-14_DP * Teff_H**0.245_DP
  ec_H(3,4) = 9.85e-14_DP * Teff_H**0.245_DP

! collisions with H2
! fit from Jaquet, Staemmler, Smith & Flower, J.Phys.B (1991)

  ec_paraH2(1,2)  = 3.46e-11_DP * Teff_H2**0.316_DP
  ec_orthoH2(1,2) = 2.70e-11_DP * Teff_H2**0.362_DP
  ec_paraH2(1,3)  = 7.07e-11_DP * Teff_H2**0.268_DP
  ec_orthoH2(1,3) = 5.49e-11_DP * Teff_H2**0.317_DP
  ec_paraH2(2,3)  = 1.44e-14_DP * Teff_H2**1.109_DP
  ec_orthoH2(2,3) = 4.64e-14_DP * Teff_H2**0.976_DP

! collisions with He
! fit from Monteiro & Flower, MNRAS 228, 101 (1987)

  ec_He(1,2) = 1.55e-12_DP * Teff_He**0.694_DP
  ec_He(1,3) = 2.52e-12_DP * Teff_He**0.681_DP
  ec_He(2,3) = 2.00e-15_DP * Teff_He**1.428_DP

! collisions avec H+

! guess (Evelyne)
! ec_Hp(1,2) = 1.0e-8_DP
! ec_Hp(1,3) = 1.0e-8_DP
! ec_Hp(2,3) = 1.0e-8_DP

! Pequignot A&A 231, 499 (1990) + erratum A&A 313, 1026 (1996)

  if (Teff_Hp < 194.664_dp) then
    toto = 1.40e-3_dp * (Teff_Hp*1.0e-2_dp)**0.90_dp
  else if (Teff_Hp < 3686.2414_dp) then
    toto = 2.14e-2_dp * (Teff_Hp*1.0e-3_dp)**1.30_dp
  else
    toto = 2.78e-1_dp * (Teff_Hp*1.0e-4_dp)**0.87_dp
  endif
  ec_Hp(1,2) = toto * 8.629e-6_dp / (sqrt(Teff_Hp) * gstoat(2))

  if (Teff_Hp < 511.9563_dp) then
    toto = 1.12e-4_dp * (Teff_Hp*1.0e-2_dp)**1.60_dp
  else if (Teff_Hp < 7510.155_dp) then
    toto = 3.90e-3_dp * (Teff_Hp*1.0e-3_dp)**1.40_dp
  else
    toto = 8.25e-2_dp * (Teff_Hp*1.0e-4_dp)**0.80_dp
  endif
  ec_Hp(1,3) = toto * 8.629e-6_dp / sqrt(Teff_Hp)

  if (Teff_Hp < 2090.2558_dp) then
    toto = 3.10e-4_dp * (Teff_Hp*1.0e-2_dp)**1.06_dp
  else
    toto = 2.29e-2_dp * (Teff_Hp*1.0e-4_dp)**0.69_dp
  endif
  ec_Hp(2,3) = toto * 8.629e-6_dp / sqrt(Teff_Hp)

!--- effective temperature with electrons : Te

! collisions with electrons
! From Pequignot, 1990, A&A, 231, 499 and Berrington, 1998, MNRAS, 293, L83

  auxt4 = 4.28e-1_DP * (Te/1.e4_DP) ** 1.43_DP / (1.00_DP + 6.05e-1_DP * (Te/1.e4_DP) ** 1.105_DP)
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstoat(4) / sqrt(Te) * gstoat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstoat(4) / sqrt(Te) * gstoat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstoat(4) / sqrt(Te) * gstoat(3) / 9.0_DP
  auxt4 = 5.88e-2_DP * (Te/1.e4_DP) ** 1.50_DP / (1.00_DP + 8.00e-1_DP * (Te/1.e4_DP) ** 1.125_DP)
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstoat(5) / sqrt(Te) * gstoat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstoat(5) / sqrt(Te) * gstoat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstoat(5) / sqrt(Te) * gstoat(3) / 9.0_DP
  auxt4 = 1.16e-1_DP * (Te/1.e4_DP) ** 0.53_DP / (1.00_DP + 1.11e-1_DP * (Te/1.e4_DP) ** 0.160_DP)
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstoat(5) / sqrt(Te)

  auxt4 = 1.49e-4_DP * Te ** 0.4565_DP
  ec_e(1,2) = 8.629e-6_DP * auxt4 / gstoat(2) / sqrt(Te)
  auxt4 = 4.98e-5_DP * Te ** 0.4955_DP
  ec_e(1,3) = 8.629e-6_DP * auxt4 / gstoat(3) / sqrt(Te)
  auxt4 = 1.83e-9_DP * Te ** 1.347_DP
  ec_e(2,3) = 8.629e-6_DP * auxt4 / gstoat(3) / sqrt(Te)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLOAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLNAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic nitrogen
! Energy levels.
!        1 - 4S J = 3/2
!        2 - 2D J = 5/2
!        3 - 2D J = 3/2
!        4 - 2P J = 1/2
!        5 - 2P J = 3/2

! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_orthoH2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions with H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperature with electrons : Te

! collisions with electrons
! From compilation of Mendoza (1983)

  auxt4 = 7.92e-6_DP * Te**0.589_DP
  ec_e(1,3) = 8.629e-6_DP * auxt4 / gstnat(3)
  auxt4 = 1.5_DP * auxt4
  ec_e(1,2) = 8.629e-6_DP * auxt4 / gstnat(2)
  auxt4 = 1.74e-6_DP * Te**0.621_DP
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstnat(4)
  auxt4 = 2.0_DP * auxt4
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 4.78e-5_DP * Te**0.431_DP
  ec_e(2,3) = 8.629e-6_DP * auxt4 / gstnat(3)
  auxt4 = 2.61e-6_DP * Te**0.609_DP
  ec_e(4,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 3.59e-4_DP * Te**0.217_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 1.13e-4_DP * Te**0.279_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 6.82e-5_DP * Te**0.301_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstnat(4)
  auxt4 = 1.65e-4_DP * Te**0.193_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstnat(4)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLNAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLSAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic sulphur

! Energy levels.
!        1 - 3P J = 3
!        2 - 3P J = 1
!        3 - 3P J = 0
!        4 - 1D J = 2
!        5 - 1S J = 0

! All collision rates as for O, but CT with H is excluded, as it is non-resonant

! collisions with H

! Fine structure

  ec_H(1,2) = 4.37e-12_DP * Teff_H**0.66_DP
  ec_H(1,3) = 1.06e-12_DP * Teff_H**0.80_DP
  ec_H(2,3) = 1.35e-11_DP * Teff_H**0.45_DP

! collisions with H2
! fit from Jaquet, Staemmler, Smith & Flower, J.Phys.B (1991)

  ec_paraH2(1,2) = 3.46e-11_DP * Teff_H2**0.316_DP
  ec_orthoH2(1,2) = 2.70e-11_DP * Teff_H2**0.362_DP
  ec_paraH2(1,3) = 7.07e-11_DP * Teff_H2**0.268_DP
  ec_orthoH2(1,3) = 5.49e-11_DP * Teff_H2**0.317_DP
  ec_paraH2(2,3) = 1.44e-14_DP * Teff_H2**1.109_DP
  ec_orthoH2(2,3) = 4.64e-14_DP * Teff_H2**0.976_DP

! collisions with He
! fit from Monteiro & Flower, MNRAS 228, 101 (1987)

  ec_He(1,2) = 1.55e-12_DP * Teff_He**0.694_DP
  ec_He(1,3) = 2.52e-12_DP * Teff_He**0.681_DP
  ec_He(2,3) = 2.00e-15_DP * Teff_He**1.428_DP

! collisions avec H+
! guess (Evelyne)

  ec_Hp(1,2) = 1.0e-8_DP
  ec_Hp(1,3) = 1.0e-8_DP
  ec_Hp(2,3) = 1.0e-8_DP

!--- effective temperature with electrons : Te

! collisions with electrons
! From Pequignot, 1990, A&A, 231, 499 and Berrington, 1998, MNRAS, 293, L83

  auxt4 = 4.28e-1_DP * (Te/1.e4_DP) ** 1.43_DP / (1.00_DP + 6.05e-1_DP * (Te/1.e4_DP) ** 1.105_DP)
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstsat(4) / sqrt(Te) * gstsat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstsat(4) / sqrt(Te) * gstsat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstsat(4) / sqrt(Te) * gstsat(3) / 9.0_DP
  auxt4 = 5.88e-2_DP * (Te/1.e4_DP) ** 1.50_DP / (1.00_DP + 8.00e-1_DP * (Te/1.e4_DP) ** 1.125_DP)
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstsat(5) / sqrt(Te) * gstsat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstsat(5) / sqrt(Te) * gstsat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstsat(5) / sqrt(Te) * gstsat(3) / 9.0_DP
  auxt4 = 1.16e-1_DP * (Te/1.e4_DP) ** 0.53_DP / (1.00_DP + 1.11e-1_DP * (Te/1.e4_DP) ** 0.160_DP)
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstsat(5) / sqrt(Te)

  auxt4 = 1.49e-4_DP * Te ** 0.4565_DP
  ec_e(1,2) = 8.629e-6_DP * auxt4 / gstsat(2) / sqrt(Te)
  auxt4 = 4.98e-5_DP * Te ** 0.4955_DP
  ec_e(1,3) = 8.629e-6_DP * auxt4 / gstsat(3) / sqrt(Te)
  auxt4 = 1.83e-9_DP * Te ** 1.347_DP
  ec_e(2,3) = 8.629e-6_DP * auxt4 / gstsat(3) / sqrt(Te)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLSAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLSiAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic silicon
! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

! All collision rates as for C

! collisions with H
! Fine structure: Launay & Roueff, 1977, A&A, 56, 289

  ec_H(1,2) = 1.01e-10_DP * Tn**0.117_DP
  ec_H(1,3) = 4.49e-11_DP * Tn**0.194_DP
  ec_H(2,3) = 1.06e-10_DP * Tn**0.234_DP

! collisions with H2
! fit from Schroder, Staemmler, Smith, Flower & Jaquet (1991)

  ec_paraH2(1,2) = 0.80e-10_DP
  ec_orthoH2(1,2)= 0.75e-10_DP
  ec_paraH2(1,3) = 0.90e-10_DP
  ec_orthoH2(1,3)= 3.54e-11_DP * Teff_H2**0.167_DP
  ec_paraH2(2,3) = 2.00e-10_DP
  ec_orthoH2(2,3)= 5.25e-11_DP * Teff_H2**0.244_DP

! collisions with He
! Fine structure: Staemmler & Flower (1991)

  ec_He(1,2) = 8.38e-12_DP * Teff_He**0.159_DP
  ec_He(1,3) = 5.98e-11_DP * Teff_He**0.078_DP
  ec_He(2,3) = 3.68e-11_DP * Teff_He**0.041_DP

! collisions with H+

  ec_Hp(1,2) = 10._DP**(-10.359_DP + 0.7959_DP*log10(Teff_Hp) - 0.08748_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(1,3) = 10._DP**(-13.232_DP + 2.4171_DP*log10(Teff_Hp) - 0.29151_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(2,3) = 10._DP**(-11.290_DP + 1.7915_DP*log10(Teff_Hp) - 0.23010_DP*(log10(Teff_Hp))**2._DP)

!--- effective temperature with electrons : Te

! collisions with electrons
! Should we use Pequignot 1990...?
! From Pequignot & Aldrovandi, 1976, A&A, 50, 141, as compiled by Mendoza (1983)

  auxt4 = 1.83e-4_DP * Te ** 0.444_DP
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstsiat(4) * gstsiat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstsiat(4) * gstsiat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstsiat(4) * gstsiat(3) / 9.0_DP
  auxt4 = 9.86e-5_DP * Te ** 0.343_DP
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstsiat(5) * gstsiat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstsiat(5) * gstsiat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstsiat(5) * gstsiat(3) / 9.0_DP
  auxt4 = 2.77e-3_DP
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstsiat(5)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  end subroutine COLSiAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLCPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2
  REAL (kind=DP)                      :: toto

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv)  = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionized carbon
! Energy levels.
!        1 - 2P J = 1/2
!        2 - 2P J = 3/2
!        3 - 4P J = 1/2
!        4 - 4P J = 3/2
!        5 - 4P J = 5/2


! collisions with H

  ec_H(1,2) = 8.86e-10_dp

! collisions with H2

  ec_paraH2(1,2) = 3.94e-10_dp
  ec_orthoH2(1,2)= 4.92e-10_dp

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions with H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperature with electrons : Te

! collisions with electrons
! k21 : Hayes & Nussbaumer, 1984, A&A, 134, 193
! Dec. 02, CM & DRF
! k21 changed to power law fit: Wilson & Bell, 2001, MNRAS, 37, 1027
! other : Lennon et al., 1985, ApJ, 294, 200

  toto = 0.65582_dp * Te**(0.13244_DP)
  ec_e(1,2) = 8.629e-6_DP * min(toto,2.5_dp) / (4.0_DP*sqrt(Te))
  ec_e(1,3) = 8.629e-6_DP * (0.288_DP - 5.89e-7_DP*Te) / (2.0_DP * sqrt(Te))
  ec_e(1,4) = 8.629e-6_DP * (0.424_DP - 8.35e-7_DP*Te) / (4.0_DP * sqrt(Te))
  ec_e(1,5) = 8.629e-6_DP * (0.257_DP - 3.95e-7_DP*Te) / (6.0_DP * sqrt(Te))
  ec_e(2,3) = 8.629e-6_DP * (0.197_DP - 3.15e-7_DP*Te) / (2.0_DP * sqrt(Te))
  ec_e(2,4) = 8.629e-6_DP * (0.544_DP - 9.84e-7_DP*Te) / (4.0_DP * sqrt(Te))
  ec_e(2,5) = 8.629e-6_DP * (1.200_DP - 2.40e-6_DP*Te) / (6.0_DP * sqrt(Te))

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLCPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLNPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv)  = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionized nitrogen
! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_H2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions with H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperature with electrons : Te

! collisions with electrons
! Mendoza, 1983, International Astrophysical Union, Symposium 103

  ec_e(1,2) =  8.629e-6_DP * (0.138_DP * Te**0.119_DP)   / (3.0_DP * sqrt(Te))
  ec_e(1,3) =  8.629e-6_DP * (0.0834_DP * Te**0.130_DP)  / (5.0_DP * sqrt(Te))
  ec_e(2,3) =  8.629e-6_DP * (0.359_DP*Te**0.125_DP)     / (5.0_DP * sqrt(Te))
  ec_e(1,4) =  8.629e-6_DP * (0.203_DP * Te**0.0414_DP)  / (5.0_DP * sqrt(Te))
  ec_e(2,4) =  8.629e-6_DP * (0.609_DP * Te**0.0414_DP)  / (5.0_DP * sqrt(Te))
  ec_e(3,4) =  8.629e-6_DP * (1.02_DP * Te**0.0414_DP)   / (5.0_DP * sqrt(Te))
  ec_e(4,5) =  8.629e-6_DP * (3.14_DP*Te**(-0.142_DP))     / (2.0_DP * sqrt(Te))
  ec_e(3,5) =  8.629e-6_DP * (0.158_DP + Te*4.26e-7_DP)  / (2.0_DP * sqrt(Te))
  ec_e(2,5) =  8.629e-6_DP * (0.0949_DP + Te*2.56e-7_DP) / (2.0_DP * sqrt(Te))
  ec_e(1,5) =  8.629e-6_DP * (0.0316_DP + Te*8.52e-8_DP) / (2.0_DP * sqrt(Te))

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLNPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLOPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2 = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv)  = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionized oxygen
! Energy levels.
!        1 - 4S J = 3/2
!        2 - 2D J = 5/2
!        3 - 2D J = 3/2
!        4 - 2P J = 3/2
!        5 - 2P J = 1/2


! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_H2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions with H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperature with electrons : Te

! collisions with electrons

! ec_e(1,2) = 0.0_DP

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLOPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLSPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2 = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionized sulphur
! Energy levels.
!        1 - 4S J = 3/2
!        2 - 2D J = 3/2
!        3 - 2D J = 5/2
!        4 - 2P J = 1/2
!        5 - 2P J = 3/2

! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_H2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions with H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperature with electrons : Te

! collisions with electrons
! Keenan et al. 1996, MNRAS, 281, 1073

  ec_e(1,2) =  8.629e-6_DP * (12.195_DP * Te**(-0.16174_DP)) /(4.0_DP*sqrt(Te))
  ec_e(1,3) =  8.629e-6_DP * (18.311_DP * Te**(-0.16186_DP)) /(6.0_DP*sqrt(Te))
  ec_e(1,4) =  8.629e-6_DP * (1.2736_DP - (1.0000e-5_DP*Te)) /(2.0_DP*sqrt(Te))
  ec_e(1,5) =  8.629e-6_DP * (2.5629_DP - (2.1143e-5_DP*Te)) /(4.0_DP*sqrt(Te))
  ec_e(2,3) =  8.629e-6_DP * (8.3293_DP - (8.8206e-5_DP*Te)) /(6.0_DP*sqrt(Te))
  ec_e(2,4) =  8.629e-6_DP * (8.0121_DP * Te**(-0.16233_DP)) /(2.0_DP*sqrt(Te))
  ec_e(2,5) =  8.629e-6_DP * (9.3033_DP * Te**(-0.12316_DP)) /(4.0_DP*sqrt(Te))
  ec_e(3,4) =  8.629e-6_DP * (6.4954_DP * Te**(-0.11783_DP)) /(2.0_DP*sqrt(Te))
  ec_e(3,5) =  8.629e-6_DP * (19.378_DP * Te**(-0.14738_DP)) /(4.0_DP*sqrt(Te))
  ec_e(4,5) =  8.629e-6_DP * (11.584_DP * Te**(-0.15835_DP)) /(4.0_DP*sqrt(Te))


  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLSPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLSiPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2 = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionized silicon
! Energy levels.
!        1 - 2P J = 1/2
!        2 - 2P J = 3/2
!        3 - 4P J = 1/2
!        4 - 4P J = 3/2
!        5 - 4P J = 5/2

! collisions with H

  ec_H(1,2) = 6.00e-10_DP

! collisions with H2
! Same as C+

  ec_paraH2(1,2) = 3.94e-10_DP
  ec_orthoH2(1,2)= 4.92e-10_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions with H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperature with electrons : Te

! collisions with electrons
! Dufton & Kingston, 1991, MNRAS, 248, 827

  ec_e(1,2) =  8.629e-6_DP * 5.6_DP / (4.0_DP*SQRT(Te))

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLSiPL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  subroutine COLFEPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! The collision strengths of Zhang and Pradhan (1995) were used to calculate
! electron collision rates.  Radiative decay rates are those of
! Nussbaumer and Storey (1988) for the first four terms and Quinet, Le
! Dourneuf and Zeippen (1996) for the rest.
!
! Deexcitation rates for collisions with H, H2, He, H+ calculated using the
! orbiting approximation, assuming a transition probability of 1/2.

  USE MODULE_PHYS_VAR
  USE MODULE_READ_FE_DATA, ONLY : gamfepl
  USE MODULE_CONSTANTS
  USE MODULE_CHEMICAL_SPECIES
  implicit none

  INTEGER        :: i, j

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: aij
  REAL (kind=DP),              DIMENSION (ntr)   :: wlk
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2
  REAL (kind=DP) :: H2_red, H_red, He_red  ! reduced masses (au)


  do i=1, ntr 
    wlk(i) = elk(iup(i)) - elk(jlo(i))
  end do

! ec matrices will be multiplied by collider density in THERMAL_LOSS
! units: cm3 s-1

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

!   Reduced masses
     H2_red =(mP/me) * (mass_Feplus/amu*mass_H2/amu) / (mass_Feplus/amu+mass_H2/amu)  
     H_red  =(mP/me) * (mass_Feplus/amu*mass_H/amu) / (mass_Feplus/amu+mass_H/amu)
     He_red =(mP/me) * (mass_Feplus/amu*mass_He/amu) / (mass_Feplus/amu+mass_He/amu)

! H2, He

	do i=1,4
	  ec_orthoH2(i,i+1:5) = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red)
          ec_paraH2(i,i+1:5)  = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red) 
            ec_He(i,i+1:5)    = AUcgs*pi* sqrt((alpha_He/bohr**3)/He_red)
	end do

	do i=6,nlv-1
         ec_orthoH2(i,i+1:nlv) = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red) 
         ec_paraH2(i,i+1:nlv)  = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red) 
           ec_He(i,i+1:nlv)    = AUcgs*pi* sqrt((alpha_He/bohr**3)/He_red) 
	end do


        do i=1,nlv-1
!  H                                                                      
          ec_H(i,i+1:nlv) = AUcgs*pi* sqrt((alpha_H/bohr**3)/H_red)   

! Electron collisional deexcitation rate coefficients (cm3 s-1)

          ec_e(i,i+1:nlv) = 8.629D-6 *  gamfepl(i,i+1:nlv) /  &
                           (gst(i)*sqrt(Te)) 
	
          end do
	

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  end subroutine COLFEPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine LINE_THERMAL_BALANCE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!
! Nov. 2002: calls routines that compute the atomic/ionic
! line emissivities and the corresponding thermal loss terms
!
! Principle: - one routine COLxxx for each atom, where specific
!  collisional deexcitation rates are calculated; this routine
!  calls:  ALPHATEFF (calculates alph's and Teff's)
!          THERMAL_LOSS (matrix inversion and cooling terms)
! To add a new coolant xxx:
!      - set Dens_cool, mass_cool, and charge_cool for species xxx
!      - copy routine COLCAT into a new routine COLxxx; only
!  numerical expressions for the rates need to be changed.
!      - replace 'cat' by 'xxx' in arguments to COLxxx
!      - give atomic parameters at beginning of this module

  USE MODULE_CHEMICAL_SPECIES
  USE MODULE_READ_FE_DATA 
  implicit none
  integer  :: i

  Cool_n = 0.0_DP
  Cool_i = 0.0_DP
  Cool_e = 0.0_DP

! ---------------------------
! Excitation of atomic carbon
! ---------------------------

  Dens_cool = Dens_C
  mass_cool = mass_C
  charge_cool = 0.0_DP

  call COLCAT (ntrcat,nlvcat,gstcat,elkcat,wlkcat,aijcat, &
               iupcat,jlocat,emicat,pop_cat)

! ---------------------------
! Excitation of atomic oxygen
! ---------------------------
  Dens_cool = Dens_O
  mass_cool = mass_O
  charge_cool = 0.0_DP

  call COLOAT (ntroat,nlvoat,gstoat,elkoat,wlkoat,aijoat, &
               iupoat,jlooat,emioat,pop_oat)

! ---------------------------
! Excitation of atomic nitrogen
! ---------------------------

  Dens_cool = Dens_N
  mass_cool = mass_N
  charge_cool = 0.0_DP

  call COLNAT (ntrnat,nlvnat,gstnat,elknat,wlknat,aijnat, &
               iupnat,jlonat,eminat,pop_nat)

! ---------------------------
! Excitation of atomic sulphur
! ---------------------------
  Dens_cool = Dens_S
  mass_cool = mass_S
  charge_cool = 0.0_DP

  call COLSAT (ntrsat,nlvsat,gstsat,elksat,wlksat,aijsat, &
               iupsat,jlosat,emisat,pop_sat)

! ---------------------------
! Excitation of atomic silicon
! ---------------------------
  Dens_cool = Dens_Si
  mass_cool = mass_Si
  charge_cool = 0.0_DP

  call COLSiAT (ntrsiat,nlvsiat,gstsiat,elksiat,wlksiat,aijsiat, &
                iupsiat,jlosiat,emisiat,pop_siat)

! ---------------------------
! Excitation of ionized carbon
! ---------------------------
  Dens_cool = Dens_Cplus
  mass_cool = mass_Cplus
  charge_cool = 1.0_DP

  call COLCPL (ntrcpl,nlvcpl,gstcpl,elkcpl,wlkcpl,aijcpl, &
               iupcpl,jlocpl,emicpl,pop_cpl)

! ---------------------------
! Excitation of ionized nitrogen
! ---------------------------
  Dens_cool = Dens_Nplus
  mass_cool = mass_Nplus
  charge_cool = 1.0_DP

  call COLNPL (ntrnpl,nlvnpl,gstnpl,elknpl,wlknpl,aijnpl, &
               iupnpl,jlonpl,eminpl,pop_npl)

! ---------------------------
! Excitation of ionized oxygen
! ---------------------------
  Dens_cool = Dens_Oplus
  mass_cool = mass_Oplus
  charge_cool = 1.0_DP

  call COLOPL (ntropl,nlvopl,gstopl,elkopl,wlkopl,aijopl, &
               iupopl,jloopl,emiopl,pop_opl)

! ---------------------------
! Excitation of ionized sulphur
! ---------------------------
  Dens_cool = Dens_Splus
  mass_cool = mass_Splus
  charge_cool = 1.0_DP

  call COLSPL (ntrspl,nlvspl,gstspl,elkspl,wlkspl,aijspl, &
               iupspl,jlospl,emispl,pop_spl)


! ---------------------------
! Excitation of ionized silicon
! ---------------------------
  Dens_cool = Dens_Siplus
  mass_cool = mass_Siplus
  charge_cool = 1.0_DP

  call COLSiPL (ntrsipl,nlvsipl,gstsipl,elksipl,wlksipl,aijsipl, &
                iupsipl,jlosipl,emisipl,pop_sipl)

! ---------------------------
! Excitation of ionized iron
! ---------------------------

  Dens_cool = Dens_Feplus
  mass_cool = mass_Feplus
  charge_cool = 1.0_DP

  call COLFEPL (ntrfepl,nlvfepl,gstfepl,elkfepl,wlkfepl,aijfepl, &
                iupfepl,jlofepl,emifepl,pop_fepl)

  end subroutine LINE_THERMAL_BALANCE

! %%%%%%%%%%%%%%%%%%%%
  subroutine ALPHATEFF
! %%%%%%%%%%%%%%%%%%%%

! Nov. 2002: initializes automatically the alpha's and Teff's
!            according to the charge and mass of the coolant

  USE MODULE_CHEMICAL_SPECIES
  USE MODULE_PHYS_VAR
  USE MODULE_CONSTANTS

  implicit none

  REAL (kind=DP) :: T1, T2     ! partial effective temperatures

! neutral atom

  if (charge_cool == 0.0_DP) then

    alp_H  = 1.0_DP
    alp_H2 = 1.0_DP
    alp_He = 1.0_DP
    alp_Hp = mass_Hplus / (mass_Hplus + mass_cool)

!--- effective temperature with neutrals : Tn

    Teff_H  = Tn
    Teff_H2 = Tn
    Teff_He = Tn

!--- effective temperature with H+ :

    T1 = ABS_DeltaV * ABS_DeltaV * mass_Hplus * mass_cool &
     / (3.0_DP * kB * (mass_Hplus + mass_cool))
    T2 = (mass_cool * Tn + mass_Hplus * Ti) / (mass_Hplus + mass_cool)
    Teff_Hp = T1 + T2

  else

! positive ion

    alp_H  = mass_cool / (mass_H + mass_cool)
    alp_H2 = mass_cool / (mass_H2 + mass_cool)
    alp_He = mass_cool / (mass_He + mass_cool)
    alp_Hp = 0.0_DP

!--- effective temperature with neutrals :

    T1 = ABS_DeltaV * ABS_DeltaV * mass_H * mass_cool &
     / (3.0_DP * kB * (mass_H + mass_cool))
    T2 = (mass_cool * Ti + mass_H * Tn) / (mass_H + mass_cool)
    Teff_H  = T1 + T2
    T1 = ABS_DeltaV * ABS_DeltaV * mass_H2 * mass_cool &
     / (3.0_DP * kB * (mass_H2 + mass_cool))
    T2 = (mass_cool * Ti + mass_H2 * Tn) / (mass_H2 + mass_cool)
    Teff_H2 = T1 + T2
    T1 = ABS_DeltaV * ABS_DeltaV * mass_He * mass_cool &
     / (3.0_DP * kB * (mass_He + mass_cool))
    T2 = (mass_cool * Ti + mass_He * Tn) / (mass_He + mass_cool)
    Teff_He = T1 + T2

!--- effective temperature with H+ : Ti

    Teff_Hp = Ti

  end if

!--- effective temperature for electron collisions : Te

  Teff_e = Te

  end subroutine ALPHATEFF

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                           ec_paraH2,ec_He,ec_Hp,ec_e,emi,pop,ntr,nlv)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Nov. 2002:
! * finish calculating statistical equilibrium matrix
!   using de-excitation rates computed in subroutine COLxxx
!       - multiply by collider density
! 	    - derive excitation rates from detailed balance
!       - compute diagonal terms
!       - add radiative transitions
! * solve for level populations using LAPACK routine dgesv
! * compute radiative emissivities in each line => emi
! * compute thermal losses of all three fluids arising from
!   collisional excitation of these lines => Cool_n, Cool_i, Cool_e

  USE MODULE_CONSTANTS
  USE MODULE_CHEMICAL_SPECIES
  USE MODULE_H2

  implicit none

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi

  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He,ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2, ec_paraH2
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H2, aa ! total H2, total matrix
  REAL (kind=DP), DIMENSION (nlv)     :: pop       ! fractional level populations

!  Used by LAPACK

  INTEGER                            :: nrhs = 1
  INTEGER                            :: info
  REAL (kind=DP), dimension (nlv)    :: indx

  INTEGER        :: i, j, itr
  REAL (kind=DP) :: radi, boltz, lasum
  REAL (kind=DP) :: bth_H, bth_H2, bth_He, bth_Hp, bth_e

  bth_H  = 0.0_DP
  bth_H2 = 0.0_DP
  bth_He = 0.0_DP
  bth_Hp = 0.0_DP
  bth_e  = 0.0_DP

! finish calculating the collisional excitation matrix:
!       - multiply by collider density
! 	    - derive excitation rates from detailed balance
!       - compute diagonal terms
!       - add radiative transitions

  ec_H = ec_H * Dens_H
  ec_He = ec_He * Dens_He
  ec_Hp = ec_Hp * Dens_Hplus
  ec_e  = ec_e * Dens_e
  ec_paraH2  = ec_paraH2 * Dens_paraH2
  ec_orthoH2 = ec_orthoH2 * Dens_orthoH2

  ec_H2 = ec_paraH2 + ec_orthoH2

  do i=2,nlv
    do j=i-1,1,-1
      boltz = exp(-(elk(i)-elk(j)) / Teff_H) * gst(i) / gst(j)
      ec_H(i,j)  = ec_H(j,i)  * boltz
      boltz = exp(-(elk(i)-elk(j)) / Teff_H2) * gst(i) / gst(j)
      ec_H2(i,j) = ec_H2(j,i) * boltz
      boltz = exp(-(elk(i)-elk(j)) / Teff_He) * gst(i) / gst(j)
      ec_He(i,j) = ec_He(j,i) * boltz
      boltz = exp(-(elk(i)-elk(j)) / Teff_Hp) * gst(i) / gst(j)
      ec_Hp(i,j) = ec_Hp(j,i) * boltz
      boltz = exp(-(elk(i)-elk(j)) / Teff_e) * gst(i) / gst(j)
      ec_e(i,j) = ec_e(j,i) * boltz
    end do
  end do

  do i=1,nlv
    ec_H(i,i)  = - SUM(ec_H(:,i))
    ec_H2(i,i) = - SUM(ec_H2(:,i))
    ec_He(i,i) = - SUM(ec_He(:,i))
    ec_Hp(i,i) = - SUM(ec_Hp(:,i))
    ec_e(i,i)  = - SUM(ec_e(:,i))
  end do

! Total collision matrix

  aa = ec_H + ec_H2 + ec_He + ec_Hp + ec_e

! add radiative decay terms

  do itr=1,ntr
    i = iup(itr)
    j = jlo(itr)
    aa(i,i) = aa(i,i) - aij(itr)
    aa(j,i) = aa(j,i) + aij(itr)
  end do

! last row of transition matrix from closure relation
! fill the matrix on the right-hand side of the equation

  aa(nlv,:) = 1.0_DP
  pop = 0.0_DP
  pop(nlv) = 1.0_DP

! Solve with LAPACK routine dgesv
 
  call dgesv (nlv, nrhs, aa, nlv, indx, pop, nlv, info)

  lasum = 0.0_DP
  do i=1,nlv
    pop(i) = max(pop(i),0.0_DP)
    lasum = lasum + pop(i)
  end do
  pop = pop / lasum

! print *, '    level energies (K):', elk
! print *, '    fractional populations:', pop

! calculate kinetic energy lost by each collider in excitation process

  do j=1,nlv-1
    do i=j+1,nlv
      bth_H  = bth_H  + (ec_H(j,i)  * pop(i) - ec_H(i,j)  * pop(j)) * (elk(i) - elk(j))
      bth_H2 = bth_H2 + (ec_H2(j,i) * pop(i) - ec_H2(i,j) * pop(j)) * (elk(i) - elk(j))
      bth_He = bth_He + (ec_He(j,i) * pop(i) - ec_He(i,j) * pop(j)) * (elk(i) - elk(j))
      bth_Hp = bth_Hp + (ec_Hp(j,i) * pop(i) - ec_Hp(i,j) * pop(j)) * (elk(i) - elk(j))
      bth_e  = bth_e  + (ec_e(j,i)  * pop(i) - ec_e(i,j)  * pop(j)) * (elk(i) - elk(j))
    end do
  end do


! calculate total radiative losses

  radi = 0.0_DP
  do i=1,ntr
    emi(i) = wlk(i) * aij(i) * pop(iup(i))
    radi = radi + emi(i)
  end do

!  print *, '  total radiative losses =', radi

! Check that total thermal loss = total radiative emission to within 1 part in 1000

  if (abs(bth_H+bth_H2+bth_He+bth_Hp+bth_e+radi) > 1.0e-3*abs(radi)) then
    print *, bth_H, bth_H2, bth_He, bth_Hp, bth_e
    print *, bth_H+bth_H2+bth_He+bth_Hp+bth_e, radi, mass_cool/amu 
  endif

  Cool_n = Cool_n + kB * (bth_H * alp_H + bth_H2 * alp_H2 &
                  + bth_He * alp_He + bth_Hp * alp_Hp) * Dens_cool
  Cool_i = Cool_i + kB * (bth_H * (1.0_DP-alp_H) + bth_H2 * (1.0_DP-alp_H2) &
           + bth_He * (1.0_DP-alp_He) + bth_Hp * (1.0_DP-alp_Hp)) * Dens_cool
  Cool_e = Cool_e + kB * bth_e * Dens_cool

! Nov. 2002: multiply emi (K s-1) by Dens_cool (cm-3) and kB (erg K-1) 
!            to yield result in erg cm-3 s-1 

  emi = emi * kB * Dens_cool

  end subroutine THERMAL_LOSS

! %%%%%%%%%%%%%%%%%%%%%
  subroutine LINE_INTEG
! %%%%%%%%%%%%%%%%%%%%%

  USE MODULE_CONSTANTS
  USE MODULE_PHYS_VAR
  USE MODULE_CHEMICAL_SPECIES
  implicit none

  intcat = intcat + (emicat_o + emicat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intnat = intnat + (eminat_o + eminat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intoat = intoat + (emioat_o + emioat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intsat = intsat + (emisat_o + emisat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intsiat = intsiat + (emisiat_o + emisiat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intcpl = intcpl + (emicpl_o + emicpl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intnpl = intnpl + (eminpl_o + eminpl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intopl = intopl + (emiopl_o + emiopl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intspl = intspl + (emispl_o + emispl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intsipl = intsipl + (emisipl_o + emisipl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intfepl = intfepl + (emifepl_o + emifepl) * &
         0.5_DP * dist_step / (4._DP*pi) 

  emicat_o  = emicat
  eminat_o  = eminat
  emioat_o  = emioat
  emisat_o  = emisat
  emisiat_o = emisiat
  emicpl_o  = emicpl
  eminpl_o  = eminpl
  emiopl_o  = emiopl
  emispl_o  = emispl
  emisipl_o = emisipl
  emifepl_o = emifepl

  end subroutine LINE_INTEG

END MODULE MODULE_LINE_EXCIT
