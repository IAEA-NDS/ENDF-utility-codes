! **********************************************************************
* Copyright (C) 2004 Brookhaven National Laboratory
* Copyright (C) 2005 International Atomic Energy Agency
* -----------------------------------------------------------------------------
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is furnished
* to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*
*-----------------------------------------------------------------------------
! *
!+++MDC+++
!...VMS, ANS, WIN, UNX
!
!     Main program for non-windows implementation of FIZCON
!
      PROGRAM FIZCON
!
      IMPLICIT NONE
!...LWI, DVF, MOD
!/!
!/!     Module implementation of FIZCON for MODLIB and WINDOWS
!/!
!/      MODULE FIZCON
!/!
!/      IMPLICIT NONE
!/!
!/      PRIVATE
!/!
!/      PUBLIC :: RUN_FIZCON
!/      PUBLIC :: FIZCON_INPUT, FIZCON_DATA, FIZCON_SUCCESS
!...LWI, DVF
!/      PUBLIC :: Default_epsiln, epsiln3
!---MDC---
!-T Program FIZCON
!-P Check procedures and data in evaluated nuclear data files
!-P in ENDF-5 or ENDF-6 format
!-V
!-V         Version 8.19   December 2017 D. Brown
!-V                        - Add checks of P(nu) for fission
!-V                        - Add checks of fission energy release tables
!-V         Version 8.18   February 2015   A. Trkov
!-V                        Check that NLS for the unresolved resonance
!-V                        range in MF32 less or equal NLS in MF2
!-V         Version 8.17   October 2014   A. Trkov
!-V                        Check product ZA in MF10, MF40
!-V         Version 8.16   September 2014   A. Trkov
!-V                        Skip warning on AWI test for electrons
!-V                        Deactivate E-range test for MT=505,506
!-V         Version 8.15   March 2014   A. Trkov
!-V                        Increase the upper limit of nu-bar from 10 to 20
!-V                        (relevant for high-energy files extending to 200 MeV)
!-V         Version 8.14   October 2012   A. Trkov
!-V                        Allow E-dependent scattering radius in URR
!-V         Version 8.13   October 2012   A. Trkov
!-V                        Trivial fix of an IF statement in CK35
!-V         Version 8.12   October 2012   A. Trkov
!-V                        Check for negative resonance widths
!-V         Version 8.11   September 2012   A. Koning
!-V                        Cleanup unused variables.
!-V         Version 8.10   August 2011   M. White
!-V                        Upgrade for MF1/MT458 extension
!-V         Version 8.09   August 2011   A. Trkov
!-V                        - Improved error counting
!-V                        - Fix checking No. of gammas in MF14 when all are
!-V                          isotropic (found by A. Koning)
!-V         Version 8.08   February 2011   A. Trkov
!-V                        Improve testing of cross-reaction covariances
!-V         Version 8.07   January 2011    A. Trkov
!-V                        Fix option LRF=7, LCOMP=2
!-V         Version 8.06   January 2011    A. Trkov
!-V                        Implement resolved resonance option LRF=7
!-V         Version 8.05   December 2010    A. Trkov
!-V                        1. Disable testing for file completeness
!-V                           in derived files (LDRV>0)
!-V                        2. Fix bug reading/writing scratch files.
!-V         Version 8.04   October 2010    A. Trkov
!-V                        Fix test on ISR
!-V         Version 8.03   June 2009   A. Trkov
!-V                        1. Make error printout conditional for IZAPT>0 in MF10
!-V                        2. Skip test on QM in MF10 for MT5
!-V         Version 8.02   April 2009   A. Koning
!-V                        1. Fix ZA of the (z,3n+p) residual.
!-V         Version 8.01   December  2008     A. Trkov
!-V                        1. Increased precision of proton mass ratio.
!-V         Version 8.00   August  2008     A. Trkov
!-V                        1. Major updating of the code.
!-V                        2. Further reduction of non-essential output.
!-V                        3. Read real variables in double precision to
!-V                           avoid reporting errors reading small numbers.
!-V                        4. Implement extended features of the format
!-V                           endorsed by CSEWG.
!-V                        5. Include corrections to problems identified
!-V                           by A. Koning.
!-V         Version 7.04   April 2006     M. Herman
!-V                        1. 'LB=8 SECTION MISSING' NOT CONSIDERED AN ERROR
!-V         Version 7.03   February 2006     M. Herman
!-V                        1. Increased dimensions for workspace and number
!-V                           of gamma rays
!-V         Version 7.02   May 2005     C.L.Dunford
!-V                        1. ONLY ERRORS REPORTED IN OUTPUT
!-V         Version 7.01   April 2005     C.L.Dunford
!-V                        1. Set success flag after return from begin
!-V                        2. Fixed valid level check for an isomer
!-V                        3. Fix subsection energy range test in CKF9
!-V                        4. changed lower limit on potential
!-V                           scattering test
!-V                        5. Fixed error in J-value test WHEN L=0 and I=0
!-V                        6. Added one more significant figure to union
!-V                           grid check and sum mup output messages
!-V                        7. Partial fission cross sections MT=19,20,21
!-V                           and 38 do not reqire secondary energy
!-V                           distributions in file 5.
!-V                        8. Correct product test for elastic scattering
!-V                        9. Move potential scattering test to PSYCHE.
!-V         Version 7.00   October 2004     C.L.Dunford
!-V                        1. Modified to provide a module for the NEA
!-V                           MODLIB project
!-V                        2. Allow energy dependent delayed fission
!-V                           group parameters.
!-V                        4. Permit user to supply batch input file
!-V                           name
!-V                        5. Removed fortran line controls from output
!-V                        6. Added command line input to Unix and
!-V                           windows versions. note: only input and
!-V                           output file names can be given. default
!-V                           options are assumed unles third
!-V                           parameter is N.
!-V
!-V
!-V      Refer all comments and inquiries to
!-V
!-V         NATIONAL NUCLEAR DATA CENTER
!-V         BUILDING 197D
!-V         BROOKHAVEN NATIONAL LABORATORY
!-V         P.O. BOX 5000
!-V         UPTON, NY 11973-5000
!-V         USA
!-V
!-V      TELEPHONE           631-344-2902
!-V      E-MAIL              NNDC@BNL.GOV
!-V
!-M
!-M FIZCON - Execute the ENDF-file checking process Phase-II
!-M ========================================================
!-M
!-M FIZCON is a program for checking that an evaluated data file
!-M has valid data and conforms to recommended procedures. It can
!-M recognize the difference between ENDF-6 and ENDF-5 formats
!-M and performs its tests accordingly. Some of the tests performed
!-M include:
!-M   1. Data arrays are in increasing energy order,
!-M   2. Resonance parameter widths add up to the total,
!-M   3. Q-values are reasonable and consistent,
!-M   4. No required sections are missing and all cover the proper
!-M      energy range,
!-M   5. Secondary distributions are normalized to 1.0,
!-M   6. Energy conservation in decay spectra.
!-M
!-M Optional tests can be performed to check that redundant cross
!-M sections such as the inelastic cross section has an energy grid
!-M which is the union of all its components and the cross section
!-M values are the sum of the component values at each energy
!-M (SUMUP test). Also optionally, algorithms are used to check for
!-M possible incorrect entry of data values (Deviant Point test).
!-M It is assumed the file being checked has passed the CHECKR
!-M program without any errors being detected.
!-M
!-M Fortran Logical Units Used:
!-M     5  Default (keyboard) input
!-M     6  Default output (terminal screen)
!-M    20  Input data file, ENDF format
!-M    21  Message file for program checking results
!-M 22,23  Temporary paging files for large data arrays
!-M 24,25  Temporary files for the SUMUP tests
!-M
!-M Input Requirements:
!-M In batch mode operation, the user must supply the following
!-M control information repeated for each input file to be processed.
!-M
!-M Record  Description
!-M      1  Source ENDF filename
!-M      2  Output report filename
!-M         (if blank, messages go to standard output file on unit 6)
!-M      3  Flag to select standard or customised options
!-M             Y  Yes, adopt standard options (default)
!-M             N  No, proceed with the selection of special options
!-M      4  Options selection (5 fields)
!-M          - Material number where processing starts (integer)
!-M            (If zero or blank, then checking begins with the
!-M            first material)
!-M          - Material number where processing ends (integer)
!-M            (If zero or blank, then checking continues to end
!-M            of the file)
!-M          - Deviant point test control (character)
!-M             Y  Do the test
!-M             N  Do not do the test (default)
!-M          - SUMUP test control (character)
!-M             Y Do the test
!-M             N Do not do the test (default)
!-M          - Fractional acceptable difference (real)
!-M            The floating point number entered here represents the
!-M            maximum fractional difference tolerated in an equality
!-M            test such as a SUMUP test. If none is entered,
!-M            the default value is .001 (1/10 of a percent).
!-M
!-M If Record 4 is left entirely blank, then the 'default' options are
!-M executed. Those are to process the entire input file, to omit the
!-M SUMUP and Deviant Point tests and to assume an allowed fractional
!-M error of .001.
!-M
!-M Multiple input files can be processed to produce multiple output
!-M files by repeating the above input data sequence. The program
!-M execution is terminated with a record containing the word DONE.
!-M
!-M In interactive mode operation, the above data is supplied in
!-M response to the appropriate query; in graphical mode, via a
!-M dialog box.
!-M
!***********************************************************************
!
!     TO CUSTOMIZE THIS SOURCE RUN SETMDC
!        ANS  -  ANSI STANDARD BATCH MODE VERSION
!        VMS  -  COMMAND MODE FOR VMS OPERATING SYSTEM
!        WIN  -  COMMAND MODE FOR PC USING DIGITAL VISUAL FORTRAN
!        UNX  -  COMMAND MODE FOR UNIX USING LAHEY FORTRAN
!        DVF  -  GRAPHICAL MODE FOR PC USING DIGITAL VISUAL FORTRAN
!        LWI  -  GRAPHICAL MODE FOR UNIX USING LAHEY WINTERACTER
!        MOD  -  MODULE FOR THE MODLIB PROJECT OF NEA WPEC
!
!     THE "ANS" VERSION MEETS F95 STANDARDS FOR FIXED OR FREE FORMAT
!       SOURCE
!     THE "VMS" VERSION WILL COMPILE WITH EITHER THE FORTRAN-77 OR
!       FORTRAN-90 VMS COMPILER
!     THE "DVF" VERSION HAS A WINDOWS GRAPHICAL INTERFACE. IT WILL
!       COMPILE WITH THE DIGITAL VISUAL FORTRAN COMPILER RUNNING
!       UNDER WINDOWS
!     THE "LWI" VERSION HAS A X-WINDOWS GRAPHICAL INTERFACE. IT WILL
!       COMPILE WITH THE LAHEY FORTRAN COMPILER WITH WINTERACTER
!       RUNNING UNDER UNIX
!
!***********************************************************************
!
!
!     FIZCON VERSION NUMBER
!
      CHARACTER(LEN=*), PARAMETER :: VERSION = '8.19'
!
!     DEFINE VARIABLE PRECISION
!
      INTEGER(KIND=4), PARAMETER :: I4 = SELECTED_INT_KIND(8)
      INTEGER(KIND=4), PARAMETER :: R4 = SELECTED_REAL_KIND(6,37)
      INTEGER(KIND=4), PARAMETER :: R8 = SELECTED_REAL_KIND(15,307)
!
      REAL(KIND=R4), PARAMETER :: FACTOR=1.008665
      REAL(KIND=R4), PARAMETER :: OTHIRD=1./3.
!
!     STANDARD FORTRAN INPUT AND OUTPUT UNITS
!
      INTEGER(KIND=I4) :: NIN
      INTEGER(KIND=I4), PARAMETER :: INPUT0 = 5, IOUT=6
      INTEGER(KIND=I4), PARAMETER :: ISCRX = 22, ISCRY = 23, ISCRXY = 24
      INTEGER(KIND=I4), PARAMETER :: ISCRU1 = 25, ISCRU2 = 26
!
!     ENDF DISK FILE INPUT AND CHECKING OUTPUT FORTRAN UNITS
!
      INTEGER(KIND=I4), PARAMETER :: JIN=20,JOUT=21
!
!     FINAL FORTRAN OUTPUT UNIT
!
      INTEGER(KIND=I4) :: NOUT
!
!     IMDC  FLAG FOR COMPILER OPTION
!     TFMT  FORMAT FOR INTERACTIVE INPUT PROMPT
!     STATUS PARAMETER FOR OPENING NEW FILE
!
!+++MDC+++
!...ANS
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 0
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = ' '
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!...VMS
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 1
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = '(/A,$)'
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'NEW'
!...WIN
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 2
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = '(/A,$)'
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!...UNX
      INTEGER(KIND=I4), PARAMETER :: IMDC = 3
      CHARACTER(LEN=*), PARAMETER :: TFMT = '(/A,$)'
      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!...DVF
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 4
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = '(A)'
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!...LWI
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 5
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = '(A)'
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!...MOD
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 6
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = '(A)'
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!---MDC---
!
!     COMMAND LINE INPUT TEXT AND TEXT LENGTH
!
      CHARACTER(LEN=100) :: INPAR
      INTEGER(KIND=I4) :: ILENP
!
!     INPUT DATA STRUCTURE
!
      TYPE FIZCON_INPUT
         CHARACTER(LEN=100) :: INFIL
         CHARACTER(LEN=100) :: OUTFIL
         INTEGER(KIND=I4) :: MATMIN
         INTEGER(KIND=I4) :: MATMAX
         INTEGER(KIND=I4) :: ICKT
         INTEGER(KIND=I4) :: ISUM
         REAL(KIND=R4) :: EPSILN
      END TYPE FIZCON_INPUT
!
      TYPE(FIZCON_INPUT) :: FIZCON_DATA
!
      INTEGER (KIND=I4) :: NERROR         !  Counted number of errors
      INTEGER (KIND=I4) :: NWARNG         !  Counted number of warnings
!
!     FLAG TO INDICATE WHETHER MULTIPLE INPUT FILES CAN BE SELECTED
!
      INTEGER(KIND=I4) :: IONEPASS        !  0, YES;  1, NO
!
!     FLAG TO INDICATE SUCCESS OR FAILURE OF STANEF EXECUTION
!
      INTEGER(KIND=I4) :: FIZCON_SUCCESS
!
!     END OF FILE FLAG
!
      INTEGER(KIND=I4) :: IFIN
!
!     FILE (TAPE) LABEL FROM FIRST RECORD
!
      CHARACTER(LEN=66) :: TLABEL
      INTEGER(KIND=I4) :: LABEL
!
      CHARACTER(LEN= 5) :: ASEQ           ! Flag to indicate a sequenced file
!                                           unsequenced if blank
!
!     LIBRARY, VERSION, SUBLIBRARY, MOD NUMBER AND FORMAT OF
!         MATERIAL BEING PROCESSED
!
      INTEGER(KIND=I4) :: NLIB,NVER,NSUB,NMOD,NFOR
!
!     MATERIAL, FILE, AND SECTION NUMBER CURRENTLY BEING PROCESSED
!
      INTEGER(KIND=I4) :: MATO,MFO,MTO
!
!     MATERIAL, FILE, AND SECTION NUMBER OF LAST ERROR DETECTED
!
      INTEGER(KIND=I4) :: MATERR,MFERR,MTERR,MISFERR
!
!     1000*Z + A OF MATERIAL CURRENTLY BEING PROCESSED
!        AWR IS THE RATIO OF THE MATERIAL MASS TO THAT OF THE NEUTRON
!        AWI IS THE RATIO OF THE PROJECTILE MASS TO THE THAT OF NEUTRON
!        STA = 0.0, STABLE MATERIAL; STA = 1.0 RADIOACTIVE MATERIAL
!        ELIS IS THE EXCITATION ENERGY OF THE TARGET NUCLEUS
!
      CHARACTER(LEN=11) :: ZSA            ! Symbolic material designation
      REAL(KIND=R4)    :: ZA,AWR,AWI,STA,ELIS
!
!     ENERGY LIMITS FOR THE MATERIAL
!
      REAL(KIND=R4), PARAMETER :: ENMIN = 1.0E-05
      REAL(KIND=R4) :: ENMAX
!
!     STORES MTS(SECTIONS) AND THEIR ENERGY SPANS
!
      INTEGER(KIND=I4), PARAMETER :: NSECMAX=1000
      INTEGER(KIND=I4) :: NXC
      INTEGER(KIND=I4), DIMENSION(NSECMAX,2) :: INDX,ENGS
!
!     LIS   IS THE STATE NUMBER (0 FOR GROUND) OF THE MATERIAL
!     LISO  IS THE ISOMER STATE NUMBER OF THE MATERIAL
!
      INTEGER(KIND=I4) :: LIS,LISO
!
!     LDRV   IS THE DERIVED FILE FLAG
!     LRP    IS THE RESONANCE PARAMETER FLAG
!     LFI IS THE FISSION FLAG

      INTEGER(KIND=I4) :: LDRV,LRP,LFI
!
!     CONTENTS OF FIELDS ON A HEAD/CONT RECORD
!
      INTEGER(KIND=I4) :: L1H,L2H,N1H,N2H
      REAL(KIND=R4)    :: C1H,C2H
!
!     MAXIMUM SIZE OF AN INTERPOLATION TABLE
!
      INTEGER(KIND=I4), PARAMETER :: INTABMAX=20
!
!     CONTENTS OF FIRST RECORD AND INTERPOLATION TABLE FOR A TAB1 RECORD
!
      INTEGER(KIND=I4) :: L1,L2,NR,NP
      INTEGER(KIND=I4), DIMENSION(INTABMAX) :: NBT,JNT
      REAL(KIND=R4)    :: C1,C2
!
!     CONTENTS OF FIRST RECORD AND INTERPOLATION TABLE FOR A TAB2 RECORD
!
      INTEGER(KIND=I4) :: L12,L22,NR2,NP2
      INTEGER(KIND=I4), DIMENSION(INTABMAX) :: NBT2,JNT2
      REAL(KIND=R4)    :: C12,C22
!
!     CONTENTS OF FIRST RECORD OF A LIST RECORD
!
      INTEGER(KIND=I4) :: L1L,L2L,NPL,N2L
      REAL(KIND=R4)    :: C1L,C2L
!
!     POSSIBLE DATA REPETITION RATES ON A LIST RECORD
!
      INTEGER(KIND=I4), PARAMETER :: NREP6 = 6,NREP12 = 12
!
!     TEXT CONTENTS ON A TEXT RECORD
!
      CHARACTER(LEN=66) :: TEXT
!
!     MATERIAL, FILE, SECTION, AND SEQUENCE NUMBER OF CURRENT RECORD
!
      INTEGER(KIND=I4) :: MATP,MFP,MTP,NSEQP
!
!     SEQUENCE NUMBER OF THE CONT-LIKE PART OF A TAB OR LIST RECORD
!
      INTEGER(KIND=I4) :: NSEQP1
!
!     FLAG INDICATING WHETHER A SUMUP TEST HAS BEEN PERFORMED
!
      INTEGER(KIND=I4) :: ITEST
!
!     SUMUP TESTS
!
      INTEGER(KIND=I4) :: ITFLE,IPC,NMTO
      INTEGER(KIND=I4), DIMENSION(250) :: MTOO
!
!     TOTAL DATA STORAGE ARRAYS FOR SUMUP TESTS
!
      INTEGER(KIND=I4), PARAMETER :: SZDAT=500000
      INTEGER(KIND=I4) :: NTOT
      REAL(KIND=R4), DIMENSION(SZDAT) :: XT,YT,YINT
      REAL(KIND=R4), DIMENSION(4,3) :: COEFS
!
!     ENERGY LIMITS OF THE RESONANCE REGION
!
      REAL(KIND=R4) :: E1,E2
!
!     MAXIMUM NUMBR OF L-STATES IN THE UNRESOLVED RESONANCE REGION
!
      INTEGER(KIND=I4) :: MF2URL
!
!     ARRAY FOR UNRESOLVED ENERGY GRID
!
      INTEGER(KIND=I4), PARAMETER :: NEUR=250
      REAL(KIND=R4), DIMENSION(NEUR) :: EURGRID
!
!     SCATTERING RADIUS CHECKING DATA
!
      INTEGER(KIND=I4) :: NRO
      REAL(KIND=R4) :: AWRI1,AWRI2
!
!     FLAG INDICATING THE PRESENCE OF FILE 3
!
      INTEGER(KIND=I4) :: IFL3
!
!     FLAG IN ALL CHARGED PARTICLE ELASTICS ARE SET TO 1.
!
      INTEGER(KIND=I4) :: CPELAS
!
!     ARRAY STORING Q-VALUES FROM FILE 3 FOR LATER TESTS
!
      INTEGER(KIND=I4), PARAMETER ::  SZMT3=250
      INTEGER(KIND=I4) :: NMT3
      INTEGER(KIND=I4), DIMENSION(SZMT3) :: MT3
      REAL(KIND=R4), DIMENSION(SZMT3) :: QMVAL,QVAL
      REAL(KIND=R4), PARAMETER :: QUNK= 7.777E+07
      REAL(KIND=R4), PARAMETER :: SPIUNK= -77.777
!
!     LIGHT PARTICLE DEFINITIONS
!
      INTEGER(KIND=I4), PARAMETER :: NPARTS=7
      INTEGER(KIND=I4), DIMENSION(NPARTS), PARAMETER ::                  &
     &         IPARTS=(/0,1,1001,1002,1003,2003,2004/)
      REAL(KIND=R4), DIMENSION(NPARTS), PARAMETER ::                     &
     &     AWPART=(/0.,1.,0.998623,1.996256,2.989596,2.989033,3.967131/)
!
!     SIGNALS FOR PRESENCE OF FILES 5 AND 6
!
      INTEGER(KIND=I4) :: NCKF5,NCKF6
!
!     NUMBER OF PARTIALS THAT CAN BE CHECK WITH A TOTAL FOR
!       REPRESENTATION CONSISTENCY
!
      INTEGER(KIND=I4), PARAMETER :: SZPAR=10
!
!     DATA FOR TEST OF TOTAL FISSION AGAINST PARTIALS
!
      INTEGER(KIND=I4), DIMENSION(SZPAR) :: ILTFIS
      INTEGER(KIND=I4) :: IMTFIS,IKTFIS
!
!     DATA FOR TEST OF TOTAL N,P AGAINST PARTIALS
!
      INTEGER(KIND=I4), DIMENSION(SZPAR) :: ILTNP
      INTEGER(KIND=I4) :: IMTNP,IKTNP
!
!     DATA FOR TEST OF TOTAL N,A AGAINST PARTIALS
!
      INTEGER(KIND=I4), DIMENSION(SZPAR) :: ILTNA
      INTEGER(KIND=I4) :: IMTNA,IKTNA
!
!     DECAY DATA CHECKING VARIABLES
!
      REAL(KIND=R4) :: T12,DT12
      REAL(KIND=R4), PARAMETER :: EMASS=.511006E+6  ! ELECTRON MASS
      REAL(KIND=R4), PARAMETER :: ALPHA=1./137.04   ! FINE STRUCTURE
      INTEGER(KIND=I4), PARAMETER ::  NDYTP=7
      REAL(KIND=R4), DIMENSION(NDYTP) :: QO,DQ,BR,DBR
      REAL(KIND=R4) QMAX,QQ,DQQ,BE,DBE,GE,DGE,AE,DAE
!
!     STORES INFOMATION ABOUT RADIOACTIVE PRODUCTS FOUND IN FILE 8
!
      INTEGER(KIND=I4), PARAMETER :: SZLMF=200
      INTEGER(KIND=I4) :: NLMF
      INTEGER(KIND=I4), DIMENSION(4,SZLMF) :: LMFS
      INTEGER(KIND=I4) :: NISSEC
      INTEGER(KIND=I4), DIMENSION(SZLMF) :: MTISO
!
!     FLAG INDICATING THAT CURRENTLY PROCESSING S(ALPHA,BETA) DATA
!
      INTEGER(KIND=I4) :: INEGC
!
!     PARAMETERS FOR THE FISSION ENERGY RELEASE TEST
!
      INTEGER(KIND=I4) :: MT458
      REAL(KIND=R4) :: ERQ
!
!     DISCRETE GAMMA RAYS SEEN IN FILES 12 AND/OR 13 AND 14
!
      INTEGER(KIND=I4), PARAMETER :: SZGAM=500
      REAL(KIND=R4), DIMENSION(SZGAM) :: EGAM
      INTEGER(KIND=I4), DIMENSION(SZGAM) :: MTGAM,MMGAM,NNGAM
      INTEGER(KIND=I4) :: NMTGAM,MMTGAM
!
!     STORES FLAG FOR MT'S SEEN IN FILE 12 AND/OR 13
!
      INTEGER(KIND=I4), PARAMETER :: SZMTS = 150
      INTEGER(KIND=I4), DIMENSION(SZMTS,2) :: ICON
      INTEGER(KIND=I4) :: NPMT
!
!     COVARAINCE TESTS

      INTEGER(KIND=I4), PARAMETER :: NCXMAX=25
      INTEGER(KIND=I4) :: NCX, NCXLAS
      INTEGER(KIND=I4), DIMENSION(NCXMAX,3) :: MTLEFT
      REAL(KIND=R4), DIMENSION(NCXMAX,2) :: EC
      INTEGER(KIND=I4), PARAMETER :: MTRMAX=100
      INTEGER(KIND=I4) :: MTR
      INTEGER(KIND=I4), DIMENSION(MTRMAX) :: MTRITE
      INTEGER(KIND=I4), PARAMETER :: NIXMAX=10
      INTEGER(KIND=I4) :: NIX
      REAL(KIND=R4), DIMENSION(NIXMAX,2) :: EI
      INTEGER(KIND=I4), PARAMETER :: NMTMAX=100
      INTEGER(KIND=I4) :: NMT33
      INTEGER(KIND=I4), DIMENSION(NMTMAX,2) :: MTNI
      INTEGER(KIND=I4), PARAMETER :: NEGMAX=120
      INTEGER(KIND=I4) :: NEG
      REAL(KIND=R4), DIMENSION(NEGMAX,2) :: EGR33
!
!     TAGS ON CURRENT RECORD
!
      INTEGER(KIND=I4) :: MAT,MF,MT,NSEQ
!
!     DATA PAGING ARRAYS
!
      INTEGER(KIND=I4), PARAMETER :: PAGESZ = 996
      INTEGER(KIND=I4) :: IPAGE
      INTEGER(KIND=I4) :: IPAGEX,ILOWX,IHIGHX
      REAL(KIND=R4) :: XP(PAGESZ)
      INTEGER(KIND=I4) :: IPAGEY,ILOWY,IHIGHY
      REAL(KIND=R8) :: YP(PAGESZ)
!
!     SUMUP DATA PAGING ARRAY
!
      INTEGER(KIND=I4) :: IPAGEXY,ILOWXY,IHIGHXY
      REAL(KIND=R4) :: YTOT(PAGESZ)
!
!     ERROR FLAG
!
      INTEGER(KIND=I4) :: IERX
!
!     ERROR MESSAGE TEXT
!
      CHARACTER(LEN=80) :: EMESS
      INTEGER(KIND=I4) :: MESS
      INTEGER(KIND=I4), PARAMETER :: MAXMES=25
!
      REAL(KIND=R4), PARAMETER :: PI=3.1415927
      REAL(KIND=R4), PARAMETER :: BIGNO=1.0E+20
      REAL(KIND=R4), PARAMETER :: EPSILN3=.001, EPSILN4=.0001,           &
     &                           EPSILN5=.00001, EPSILN6 =0.000001
      REAL(KIND=R4), PARAMETER :: DEFAULT_EPSILN=EPSILN3
!
!     COGEND DATA
!
      INTEGER(KIND=I4) :: IDDONE,IBAV,IBREM,IUNC
      INTEGER(KIND=I4) :: NZ
      REAL(KIND=R4)    :: R0,V0,W0
!
      INTEGER(KIND=I4), PARAMETER :: NZMAX=100
!
      REAL(KIND=R4), DIMENSION(6,NZMAX) :: XLEV
      DATA XLEV/0.0,0.0,0.0,0.0,0.0,0.0,                                &
     &   0.0,0.0,0.0,0.0,0.0,0.0,                                       &
     &   54.75,0.0,0.0,0.0,0.0,0.0,                                     &
     &   111.0,0.0,0.0,0.0,0.0,0.0,                                     &
     &   188.0,0.0,4.7,4.7,0.0,0.0,                                     &
     &   283.8,0.0,6.4,6.4,0.0,0.0,                                     &
     &   401.6,0.0,9.2,9.2,0.0,0.0,                                     &
     &   532.0,23.7,7.1,7.1,0.0,0.0,                                    &
     &   685.4,31.0,8.6,8.6,0.0,0.0,                                    &
     &   866.9,45.0,18.3,18.3,0.0,0.0,                                  &
     &   1072.1,63.3,31.1,31.1,0.0,0.0,                                 &
     &   1305.0,89.4,51.4,51.4,0.0,0.0,                                 &
     &   1559.6,117.7,73.1,73.1,0.0,0.0,                                &
     &   1838.9,148.7,99.2,99.2,0.0,0.0,                                &
     &   2145.5,189.3,132.2,132.2,0.0,0.0,                              &
     &   2472.0,229.2,164.8,164.8,0.0,0.0,                              &
     &   2822.4,270.2,201.6,200.0,17.5,6.2,                             &
     &   3202.9,320.0,247.3,245.2,25.3,10.2,                            &
     &   3607.4,377.1,296.3,293.6,33.9,13.9,                            &
     &   4038.1,437.8,350.0,346.4,43.7,18.9,                            &
     &   4492.8,500.4,406.7,402.2,53.8,26.3,                            &
     &   4966.4,563.7,461.5,455.5,60.3,27.4,                            &
     &   5465.1,628.2,520.5,512.9,66.5,29.3,                            &
     &   5989.2,694.6,583.7,574.5,74.1,32.7,                            &
     &   6539.0,769.0,651.4,640.3,83.9,37.5,                            &
     &   7112.0,846.1,721.1,708.1,92.9,41.6,                            &
     &   7708.9,925.6,793.6,778.6,100.7,45.1,                           &
     &   8332.8,1008.1,871.9,854.7,111.8,51.0,                          &
     &   8978.9,1096.1,951.0,931.1,119.8,54.0,                          &
     &   9658.6,1193.6,1042.8,1019.7,135.9,65.1,                        &
     &   10367.1,1297.7,1142.3,1115.4,158.1,80.5,                       &
     &   11103.2,1414.3,1247.8,1216.7,180.0,97.2,                       &
     &   11866.7,1526.5,1358.6,1323.1,203.5,115.0,                      &
     &   12657.8,1653.9,1476.2,1435.8,231.5,135.0,                      &
     &   13473.7,1782.0,1596.0,1549.9,256.5,153.0,                      &
     &   14325.6,1921.0,1727.2,1674.9,292.1,181.0,                      &
     &   15199.7,2065.1,1863.9,1804.4,322.1,206.0,                      &
     &   16104.6,2216.3,2006.8,1939.6,357.5,235.0,                      &
     &   17038.4,2372.5,2155.5,2080.0,393.6,264.0,                      &
     &   17997.6,2531.6,2306.7,2222.3,430.3,293.0,                      &
     &   18985.6,2697.7,2464.7,2370.5,468.4,323.0,                      &
     &   19999.5,2865.5,2625.1,2520.2,504.6,353.0,                      &
     &   21044.0,3042.5,2793.2,2676.9,544.0,385.0,                      &
     &   22117.2,3224.0,2966.9,2837.9,585.0,418.0,                      &
     &   23219.9,3411.9,3146.1,3003.8,627.1,453.0,                      &
     &   24350.3,3604.3,3330.3,3173.3,669.9,487.0,                      &
     &   25514.0,3805.8,3523.7,3351.1,717.5,526.0,                      &
     &   26711.2,4018.0,3727.0,3537.5,770.2,570.0,                      &
     &   27939.9,4237.5,3938.0,3730.1,825.6,617.0,                      &
     &   29200.1,4464.7,4156.1,3928.8,883.8,667.0,                      &
     &   30491.2,4698.3,4380.4,4132.2,943.7,717.0,                      &
     &   31813.8,4939.2,4612.0,4341.4,1006.0,770.0,                     &
     &   33169.4,5188.1,4852.1,4557.1,1072.1,826.0,                     &
     &   34561.4,5452.8,5103.7,4782.2,1145.0,889.0,                     &
     &   35984.6,5714.3,5359.4,5011.9,1217.1,949.0,                     &
     &   37440.6,5988.8,5623.6,5247.0,1292.8,1014.0,                    &
     &   38924.6,6266.3,5890.6,5482.7,1361.3,1074.0,                    &
     &   40443.0,6548.8,6164.2,5723.4,1434.6,1135.0,                    &
     &   41990.6,6834.8,6440.4,5964.3,1511.0,1195.0,                    &
     &   43568.9,7126.0,6721.5,6207.9,1575.3,1251.0,                    &
     &   45184.0,7427.9,7012.8,6459.3,1650.0,1311.0,                    &
     &   46834.2,7736.8,7311.8,6716.2,1722.8,1374.0,                    &
     &   48519.0,8052.0,7617.1,6976.9,1800.0,1437.0,                    &
     &   50239.1,8375.6,7930.3,7242.8,1880.8,1503.0,                    &
     &   51995.7,8708.0,8251.6,7514.0,1967.5,1573.0,                    &
     &   53788.5,9045.8,8580.6,7790.1,2046.8,1638.0,                    &
     &   55617.7,9394.2,8917.8,8071.1,2128.3,1707.0,                    &
     &   57485.5,9751.3,9264.3,8357.9,2206.5,1777.0,                    &
     &   59389.6,10115.1,9616.9,8648.0,2306.8,1853.0,                   &
     &   61332.3,10486.4,9978.2,8943.6,2398.1,1925.0,                   &
     &   63313.8,10870.4,10348.6,9244.1,2491.2,2001.0,                  &
     &   65350.8,11270.7,10739.4,9560.7,2600.9,2090.0,                  &
     &   67416.4,11681.6,11136.1,9881.1,2708.0,2180.0,                  &
     &   69525.0,12099.8,11544.0,10206.8,2819.6,2271.0,                 &
     &   71676.4,12526.7,11958.7,10535.3,2931.7,2362.0,                 &
     &   73870.8,12968.0,12385.0,10870.9,3048.5,2458.0,                 &
     &   76111.0,13418.5,12824.1,11215.2,3173.7,2558.0,                 &
     &   78394.8,13879.9,13272.6,11563.7,3296.0,2658.0,                 &
     &   80724.9,14352.8,13733.6,11918.7,3424.9,2763.0,                 &
     &   83102.3,14839.3,14208.7,12283.9,3561.6,2873.0,                 &
     &   85530.4,15346.7,14697.9,12657.5,3704.1,2990.0,                 &
     &   88004.5,15860.8,15200.0,13035.2,3850.7,3108.0,                 &
     &   90525.9,16387.5,15711.1,13418.6,3999.1,3228.0,                 &
     &   93105.0,16939.3,16244.3,13813.8,4149.4,3357.0,                 &
     &   95729.9,17493.0,16784.7,14213.5,4317.0,3489.0,                 &
     &   98404.0,18049.0,17337.1,14619.4,4482.0,3619.0,                 &
     &   101137.0,18639.0,17906.5,15031.2,4652.0,3756.0,                &
     &   103921.9,19236.7,18484.3,15444.4,4822.0,3891.0,                &
     &   106755.3,19840.0,19083.2,15871.0,5002.0,4031.0,                &
     &   109650.9,20472.1,19693.2,16300.3,5182.3,4176.0,                &
     &   112601.4,21104.6,20313.7,16733.1,5366.9,4319.0,                &
     &   115606.1,21757.4,20947.6,17166.3,5548.0,4463.0,                &
     &   118678.0,22426.8,21600.5,17610.0,5723.2,4608.0,                &
     &   121818.0,23097.2,22266.2,18056.8,5932.9,4756.0,                &
     &   125027.0,23772.9,22944.0,18504.1,6120.5,4859.0,                &
     &   128220.0,24460.0,23779.0,18930.0,6288.0,5036.0,                &
     &   131590.0,25275.0,24385.0,19452.0,6556.0,5236.0,                &
     &   135960.0,26110.0,25250.0,19930.0,6754.0,5394.0,                &
     &   139490.0,26900.0,26020.0,20410.0,6977.0,5561.0,                &
     &   143090.0,27700.0,26810.0,20900.0,7205.0,5732.0/
      REAL(KIND=R4), DIMENSION(4,NZMAX) ::  RDENS
      DATA RDENS/                                                       &
     &   0.0,0.0,0.0,0.0,                                               &
     &   0.0,0.0,0.0,0.0,                                               &
     &   0.0,0.0,0.0,0.0,                                               &
     &   0.0,0.0,0.0,0.0,                                               &
     &   0.0405,0.00008,0.0,0.0,                                        &
     &   0.0493,0.00018,0.0,0.0,                                        &
     &   0.0541,0.00024,0.00024,0.0,                                    &
     &   0.0564,0.00032,0.00063,0.0,                                    &
     &   0.0577,0.00041,0.00122,0.0,                                    &
     &   0.0584,0.00053,0.00207,0.0,                                    &
     &   0.0627,0.00069,0.00268,0.018,                                  &
     &   0.0666,0.00088,0.00339,0.037,                                  &
     &   0.0699,0.00108,0.00417,0.048,                                  &
     &   0.0729,0.00131,0.00498,0.063,                                  &
     &   0.0756,0.00155,0.00586,0.079,                                  &
     &   0.0781,0.00181,0.00682,0.093,                                  &
     &   0.0804,0.00209,0.00784,0.103,                                  &
     &   0.0824,0.00240,0.00891,0.110,                                  &
     &   0.0844,0.00272,0.0100,0.128,                                   &
     &   0.0862,0.00306,0.0112,0.144,                                   &
     &   0.0879,0.00343,0.0125,0.148,                                   &
     &   0.0896,0.00382,0.0138,0.150,                                   &
     &   0.0910,0.00424,0.0151,0.152,                                   &
     &   0.0924,0.00467,0.0165,0.147,                                   &
     &   0.0938,0.00512,0.0180,0.154,                                   &
     &   0.0950,0.00560,0.0194,0.155,                                   &
     &   0.0962,0.00610,0.0210,0.156,                                   &
     &   0.0974,0.00663,0.0225,0.156,                                   &
     &   0.0985,0.00717,0.0241,0.152,                                   &
     &   0.0995,0.00774,0.0258,0.158,                                   &
     &   0.1006,0.00834,0.0274,0.162,                                   &
     &   0.1015,0.00895,0.0291,0.166,                                   &
     &   0.1026,0.00958,0.0308,0.172,                                   &
     &   0.1035,0.0102,0.0325,0.177,                                    &
     &   0.1043,0.0109,0.0343,0.182,                                    &
     &   0.1053,0.0116,0.0361,0.187,                                    &
     &   0.1063,0.0124,0.0379,0.193,                                    &
     &   0.1071,0.0131,0.0396,0.199,                                    &
     &   0.1080,0.0139,0.0414,0.206,                                    &
     &   0.1089,0.0147,0.0432,0.210,                                    &
     &   0.1098,0.0156,0.0450,0.213,                                    &
     &   0.1107,0.0164,0.0469,0.216,                                    &
     &   0.1115,0.0173,0.0487,0.219,                                    &
     &   0.1124,0.0183,0.0505,0.222,                                    &
     &   0.1133,0.0192,0.0523,0.225,                                    &
     &   0.1142,0.0202,0.0540,0.228,                                    &
     &   0.1150,0.0212,0.0558,0.223,                                    &
     &   0.1159,0.0222,0.0576,0.238,                                    &
     &   0.1168,0.0233,0.0593,0.241,                                    &
     &   0.1178,0.0244,0.0610,0.246,                                    &
     &   0.1187,0.0255,0.0627,0.249,                                    &
     &   0.1196,0.0267,0.0644,0.253,                                    &
     &   0.1205,0.0278,0.0661,0.257,                                    &
     &   0.1215,0.0291,0.0679,0.261,                                    &
     &   0.1224,0.0303,0.0693,0.266,                                    &
     &   0.1234,0.0316,0.0708,0.271,                                    &
     &   0.1244,0.0329,0.0722,0.274,                                    &
     &   0.1254,0.0343,0.0739,0.277,                                    &
     &   0.1264,0.0357,0.0752,0.279,                                    &
     &   0.1275,0.0371,0.0767,0.281,                                    &
     &   0.1285,0.0386,0.0780,0.284,                                    &
     &   0.1296,0.0401,0.0794,0.286,                                    &
     &   0.1306,0.0417,0.0807,0.288,                                    &
     &   0.1317,0.0433,0.0819,0.290,                                    &
     &   0.1328,0.0449,0.0831,0.292,                                    &
     &   0.1340,0.0466,0.0842,0.294,                                    &
     &   0.1351,0.0483,0.0853,0.296,                                    &
     &   0.1362,0.0501,0.0864,0.298,                                    &
     &   0.1374,0.0519,0.0871,0.300,                                    &
     &   0.1386,0.0538,0.0883,0.303,                                    &
     &   0.1398,0.0557,0.0892,0.305,                                    &
     &   0.1410,0.0577,0.0900,0.308,                                    &
     &   0.1423,0.0597,0.0908,0.310,                                    &
     &   0.1436,0.0618,0.0916,0.313,                                    &
     &   0.1448,0.0639,0.0922,0.315,                                    &
     &   0.1462,0.0661,0.0929,0.318,                                    &
     &   0.1475,0.0684,0.0934,0.320,                                    &
     &   0.1489,0.0707,0.0939,0.323,                                    &
     &   0.1502,0.0730,0.0943,0.326,                                    &
     &   0.1517,0.0755,0.0948,0.329,                                    &
     &   0.1531,0.0780,0.0954,0.331,                                    &
     &   0.1546,0.0806,0.0954,0.334,                                    &
     &   0.1561,0.0833,0.0956,0.337,                                    &
     &   0.1576,0.0860,0.0958,0.341,                                    &
     &   0.1591,0.0888,0.0959,0.344,                                    &
     &   0.1607,0.0917,0.0960,0.347,                                    &
     &   0.1623,0.0947,0.0959,0.350,                                    &
     &   0.1639,0.0978,0.0955,0.353,                                    &
     &   0.1656,0.1010,0.0954,0.356,                                    &
     &   0.1673,0.1042,0.0954,0.359,                                    &
     &   0.1690,0.1076,0.0953,0.361,                                    &
     &   0.1708,0.1111,0.0952,0.364,                                    &
     &   0.1726,0.1147,0.0946,0.367,                                    &
     &   0.1744,0.1184,0.0941,0.370,                                    &
     &   0.1763,0.1222,0.0936,0.374,                                    &
     &   0.1782,0.1262,0.0931,0.377,                                    &
     &   0.1802,0.1303,0.0924,0.380,                                    &
     &   0.1829,0.1345,0.0917,0.384,                                    &
     &   0.186,0.140,0.0910,0.388,0.189,0.145,0.0903,0.392/
      REAL(KIND=R4), DIMENSION(3,NZMAX) :: BX
      DATA BX/                                                          &
     &   0.0,0.0,0.0,                                                   &
     &   0.0,0.0,0.0,                                                   &
     &   0.0,0.0,0.0,                                                   &
     &   0.0,0.0,0.0,                                                   &
     &   0.0,0.0,0.0,                                                   &
     &   0.938,0.0,0.0,                                                 &
     &   0.948,1.475,0.0,                                               &
     &   0.958,1.405,0.0,                                               &
     &   0.964,1.360,0.0,                                               &
     &   0.969,1.309,0.0,                                               &
     &   0.973,1.283,0.0,                                               &
     &   0.974,1.248,0.0,                                               &
     &   0.975,1.212,0.0,                                               &
     &   0.976,1.186,0.921,                                             &
     &   0.977,1.169,0.929,                                             &
     &   0.978,1.154,0.935,                                             &
     &   0.979,1.143,0.940,                                             &
     &   0.980,1.132,0.944,                                             &
     &   0.981,1.120,0.946,                                             &
     &   0.982,1.113,0.948,                                             &
     &   0.982,1.101,0.947,                                             &
     &   0.982,1.096,0.950,                                             &
     &   0.983,1.091,0.953,                                             &
     &   0.984,1.088,0.956,                                             &
     &   0.985,1.085,0.958,                                             &
     &   0.985,1.080,0.960,                                             &
     &   0.986,1.078,0.962,                                             &
     &   0.986,1.076,0.964,                                             &
     &   0.986,1.072,0.965,                                             &
     &   0.987,1.070,0.967,                                             &
     &   0.987,1.069,0.968,                                             &
     &   0.988,1.067,0.969,                                             &
     &   0.988,1.064,0.970,                                             &
     &   0.988,1.062,0.971,                                             &
     &   0.981,1.060,0.971,                                             &
     &   0.989,1.059,0.972,                                             &
     &   0.989,1.057,0.973,                                             &
     &   0.990,1.053,0.973,                                             &
     &   0.990,1.051,0.974,                                             &
     &   0.990,1.050,0.974,                                             &
     &   0.990,1.048,0.975,                                             &
     &   0.990,1.046,0.975,                                             &
     &   0.990,1.045,0.976,                                             &
     &   0.990,1.043,0.976,                                             &
     &   0.991,1.042,0.976,                                             &
     &   0.991,1.041,0.977,                                             &
     &   0.991,1.040,0.977,                                             &
     &   0.991,1.039,0.977,                                             &
     &   0.990,1.038,0.978,                                             &
     &   0.991,1.037,0.978,                                             &
     &   0.991,1.036,0.978,                                             &
     &   0.991,1.035,0.979,                                             &
     &   0.991,1.034,0.979,                                             &
     &   0.991,1.033,0.979,                                             &
     &   0.992,1.032,0.979,                                             &
     &   0.992,1.032,0.979,                                             &
     &   0.992,1.031,0.980,                                             &
     &   0.992,1.030,0.980,                                             &
     &   0.992,1.030,0.980,                                             &
     &   0.992,1.029,0.980,                                             &
     &   0.992,1.028,0.980,                                             &
     &   0.992,1.028,0.980,                                             &
     &   0.992,1.028,0.980,                                             &
     &   0.992,1.027,0.980,                                             &
     &   0.992,1.027,0.981,                                             &
     &   0.992,1.027,0.981,                                             &
     &   0.992,1.026,0.981,                                             &
     &   0.992,1.026,0.981,                                             &
     &   0.992,1.026,0.981,                                             &
     &   0.992,1.025,0.981,                                             &
     &   0.992,1.025,0.981,                                             &
     &   0.992,1.024,0.981,                                             &
     &   0.992,1.024,0.981,                                             &
     &   0.992,1.024,0.982,                                             &
     &   0.992,1.024,0.982,                                             &
     &   0.992,1.024,0.982,                                             &
     &   0.992,1.023,0.982,                                             &
     &   0.992,1.023,0.982,                                             &
     &   0.992,1.022,0.982,                                             &
     &   0.992,1.022,0.982,                                             &
     &   0.992,1.022,0.982,                                             &
     &   0.992,1.022,0.982,                                             &
     &   0.992,1.022,0.982,                                             &
     &   0.992,1.022,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,                                             &
     &   0.992,1.021,0.982,0.992,1.021,0.982/
!
!***********************************************************************
!
!+++MDC+++
!...VMS, ANS, WIN, UNX
!
      CALL RUN_FIZCON
!
!     TERMINATE JOB
!
      IF(FIZCON_SUCCESS.EQ.0) THEN
         WRITE(IOUT,'(/A)') ' '
         STOP '    FIZCON - Tests completed successfully'
      ELSE
         WRITE(IOUT,'(/A)') ' '
         STOP '    FIZCON - Tests terminated abnormally!'
      END IF
!---MDC---
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_FIZCON
!
!     EXECUTES FIZCON PROCESS
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     CHARACTER(LEN=*), INTRINSIC :: TRIM
!
      CHARACTER(LEN=80) :: IFIELD
      INTEGER(KIND=I4) :: IQUIT   ! FLAG TO INDICATE WHETHER OR NOT TO EXIT
      INTEGER(KIND=I4) :: IFIND   ! FLAGS WHETHER DESIRED MATERIAL FOUND
      INTEGER(KIND=I4) :: MFN
      INTEGER(KIND=I4) :: MATT,MFT,MTT,NSEQT,NSEQB
      INTEGER(KIND=I4) :: I
!
      CHARACTER(LEN=*), PARAMETER :: DASHES = REPEAT('-',80)
!
!     OUTPUT PROGRAM IDENTIFICATION
!
      NERROR=0
      NWARNG=0
!
      FIZCON_SUCCESS = 0
      IF(IMDC.LT.4) THEN
         WRITE(IOUT,'(/2A)') ' PROGRAM FIZCON VERSION ',VERSION
      END IF
!
!     CHECK FOR COMMAND LINE INPUT (VMS ONLY)
!
      IONEPASS = 0
      CALL GET_FROM_CLINE
!
!     INITIALIZE RUN
!
   10 CALL BEGIN(IQUIT)
      IF(IQUIT.GT.0)    THEN
         IF(IONEPASS.EQ.1) FIZCON_SUCCESS = 1
         GO TO 100
      END IF
!
!     CHECK LABEL AND FIND STARTING MATERIAL
!
      CALL SEARCH(IFIND)
      IF(IFIND.EQ.0)   GO TO 50
!
!     UNEXPECTED END OF FILE ENCOUNTERED
!
   20 IF(IERX.EQ.2) THEN
         IF(IMDC.LT.4) THEN
            WRITE(IOUT,'(//5X,2A)')  'END OF FILE ENCOUNTERED BEFORE ', &
     &                      'TEND RECORD FOUND!'
         END IF
         IF(NOUT.NE.IOUT)   THEN
            WRITE(NOUT,'(//5X,2A)')  'END OF FILE ENCOUNTERED BEFORE ', &
     &                      'TEND RECORD FOUND!'
         END IF
         IF(NOUT.NE.IOUT) THEN
           WRITE(NOUT,'(A)') ' Done FIZCON'
           CLOSE(UNIT=NOUT)
         END IF
         CLOSE(UNIT=JIN)
         FIZCON_SUCCESS = 1
         GO TO 100
      END IF
!
!     PROCESS NEXT SECTION
!
      IF(MAT.NE.MATO)  THEN    !NEW MATERIAL
         IF(FIZCON_DATA%MATMAX.NE.0.AND.MAT.GT.FIZCON_DATA%MATMAX)      &
     &                  GO TO 70
         INEGC  = 0
         NSEQP1 = NSEQP
         MATO = MAT
         MFO = 0
         IFL3 = 0
         E1 = BIGNO
         E2 = 0.
         NPMT = 0
         NMT3 = 0
         NLMF = 0
         NMTGAM = 0
         MMTGAM = 0
         NISSEC = 0
         NCKF5 = 0
         NCKF6 = 0
         MF2URL= 0
         REWIND (UNIT=ISCRX)
         REWIND (UNIT=ISCRY)
         REWIND (UNIT=ISCRXY)
         REWIND (UNIT=ISCRU1)
         REWIND (UNIT=ISCRU2)
!        Unconditional printing of material header
!        WRITE(NOUT,'(A/1X,A,I5)')  CHAR(12),'Check material',MAT
!        IF(NOUT.NE.IOUT)  THEN
!           IF(IMDC.LT.4) WRITE(IOUT,'(/A)')  '   '
!        END IF
      END IF
      IF(MF.NE.MFO)   THEN
         MFO = MF
         IF(MF.GE.31) THEN
            NCX = 0
            NCXLAS = 0
            MTR = 0
            NEG = 0
            NMT33 = 0
         END IF
      END IF
!
!     NEW SECTION
!
      MTO = MT
!
!     IN INTERACTIVE MODE OUTPUT CURRENT SECTION ID TO TERMINAL
!
      IF(NOUT.NE.IOUT) CALL OUT_STATUS
!
!     CHECK THE NEW SECTION
!
      CALL CHKSEC
      IF(IERX.EQ.2)  GO TO 20
!
!     IF FATAL ERROR FOUND, SKIP REST OF SECTION
!
   35 IF(IERX.NE.0) THEN
         IERX = 2
         NSEQB= NSEQ
         DO WHILE (MT.NE.0)
            READ(JIN,'(A)',END=20) IFIELD
            READ(IFIELD,'(66X,I4,I2,I3,I5)',ERR=40)  MAT,MF,MT,NSEQ
   40       CONTINUE
         END DO
         WRITE(EMESS,'(A,I3,A,I4,2A,I6,A,I6)')                          &
     &      'MF=',MFO,' MT=',MTO,' CAN NOT BE CHECKED FROM SEQUENCE ',  &
     &      'NUMBER ',NSEQB,' TO',NSEQ
         CALL ERROR_MESSAGE(0)
         IERX = 0
      END IF
!
!     READ UNTIL HEAD OR TEND RECORD FOUND
!
   50 IF(MAT.GE.0)  THEN
   55    CALL RDHEAD(I)
         IF(IERX.GE.1) GO TO 35
         IF(I.GT.1.AND.I.LT.5)   THEN
            GO TO 55
         ELSE IF(I.EQ.5) THEN
            IFIN = 1
         END IF
      ELSE
         GO TO 100
      END IF
!
!     END OF FILE SUM UP TESTS
!
      IF(MF.NE.MFO.OR.IFIN.NE.0)   THEN
         IF(ITEST.GT.0)  THEN
            MATT = MAT
            MFT = MF
            MTT = MT
            NSEQT = NSEQ
            IF(MFO.EQ.1)   THEN
               CALL SUM452(0)
            ELSE IF(MFO.EQ.3) THEN
               CALL SUMF3(0)
            ELSE IF(MFO.EQ.23)  THEN
               CALL SUMGAM(0)
            END IF
            ITEST = 0
            MAT = MATT
            MF = MFT
            MT = MTT
            NSEQ = NSEQT
         END IF
!
!        CHECK FOR MISSING SECTIONS
!
         IF(IFIN.EQ.1.OR.MAT.NE.MATO)  THEN
            MFN = 100
         ELSE
            MFN = MF - 1
         END IF
         IF(LDRV.EQ.0) CALL EFTEST(MFO,MFN)
      END IF
!
!     CHECK END OF TAPE FLAG
!
      IF(IFIN.EQ.0) THEN
         IF(FIZCON_DATA%MATMAX.EQ.0.OR.MAT.LE.FIZCON_DATA%MATMAX)       &
     &                   GO TO 20
      END IF
!
!     CLOSE FILES
!
   70 CONTINUE
!
         IF(NERROR.EQ.0 .AND. NWARNG.EQ.0) THEN
           WRITE(IOUT,'(/A)') ' No problems to report'
           IF(IOUT.NE.NOUT)
     &     WRITE(NOUT,'(/A)') ' No problems to report'
         ELSE
           WRITE(IOUT,'(/A,2(I6,A))')
     &           '     Encountered',NERROR,' errors,  '                 &
     &                             ,NWARNG,' warnings'
           IF(IOUT.NE.NOUT)                                             &
     &     WRITE(NOUT,'(A/A,2(I6,A))') CHAR(12),
     &           '     Encountered',NERROR,' errors,  '                 &
     &                             ,NWARNG,' warnings'
           NERROR=0
           NWARNG=0
         END IF
!
      IF(NOUT.NE.IOUT) THEN
        WRITE(NOUT,'(A)') ' Done FIZCON'
        CLOSE(UNIT=NOUT)
      END IF
      CLOSE(UNIT=JIN)
      CLOSE(UNIT=ISCRX,STATUS='DELETE')
      CLOSE(UNIT=ISCRY,STATUS='DELETE')
      CLOSE(UNIT=ISCRXY,STATUS='DELETE')
      CLOSE(UNIT=ISCRU1,STATUS='DELETE')
      CLOSE(UNIT=ISCRU2,STATUS='DELETE')
!
!     SEE IF ONE PASS LIMIT SET
!
      IF(IONEPASS.EQ.0) GO TO 10

  100 RETURN
      END SUBROUTINE RUN_FIZCON
!
!***********************************************************************
!
      SUBROUTINE BEGIN(IQUIT)
!
!     ROUTINE TO SET UP JOB
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IQUIT
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     CHARACTER(LEN=*), INTRINSIC :: TRIM
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=I4), INTRINSIC :: ICHAR
!
      CHARACTER(LEN=4) :: BUF1
      CHARACTER(LEN=12) :: BUF2
      CHARACTER(LEN=1) :: IW
      CHARACTER(LEN=11) :: ADATE
      CHARACTER(LEN=50) :: MATSIN
      LOGICAL(KIND=I4) :: IEXIST
      INTEGER(KIND=I4) :: IC
      REAL(KIND=R4) :: EPS
!
!     INITIALIZE PROCESSING CONTROL VARIABLES
!
      IERX = 0
      MATO = 0
      MFO = 0
      MTO = 0
      MATERR = 0
      MFERR = 0
      MTERR = 0
      MISFERR = 0
      IFIN = 0
      IBAV = 2
      IBREM = 1
      IUNC = 1
      ITEST = 0
      IDDONE = 0
      NOUT = IOUT
   10 IQUIT = 0
!
!     INITIALIZE TO STANDARD OPTIONS
!
      IF(IMDC.LT.4) THEN
         FIZCON_DATA%INFIL = '*'
         FIZCON_DATA%OUTFIL = '*'
         FIZCON_DATA%MATMIN = 0
         FIZCON_DATA%MATMAX = 0
         FIZCON_DATA%ICKT = 1
         FIZCON_DATA%ISUM = 0
         FIZCON_DATA%EPSILN = DEFAULT_EPSILN
      END IF
      SELECT CASE (IMDC)
         CASE (0)
            IW = 'N'
            IONEPASS = 0
         CASE(1,2,3)
            IF(ILENP.NE.0)  THEN
               CALL TOKEN(INPAR,'%',1,FIZCON_DATA%INFIL)
               CALL TOKEN(INPAR,'%',2,FIZCON_DATA%OUTFIL)
               CALL TOKEN(INPAR,'%',3,IW)
               IC = ICHAR(IW)
               IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
               IF(IW.EQ.' ') THEN
                  IW = 'Y'
               ELSE IF(IW.NE.'Y'.AND.IW.NE.'N') THEN
                  IW = '*'
               END IF
               IONEPASS = 1
            ELSE
               IW = '*'
               IONEPASS = 0
            END IF
         CASE (4,5,6)
            IW = 'N'
            IONEPASS = 1
      END SELECT
!
!     GET INPUT FILE SPECIFICATION
!
      IF(IMDC.LT.4) THEN
         IF(FIZCON_DATA%INFIL.EQ.'*') THEN
            IF(IMDC.NE.0) THEN
               WRITE(IOUT,FMT=TFMT)                                     &
     &             ' Input File Specification             - '
            END IF
            READ(NIN,'(A)') FIZCON_DATA%INFIL
         ELSE
            WRITE(IOUT,'(/2A)') ' Input file - ',                       &
     &               TRIM(FIZCON_DATA%INFIL)
         END IF
      END IF
!
!     SEE IF INPUT INDICATES FILE TERMINATION
!
      IF(FIZCON_DATA%INFIL.EQ.' '.OR.FIZCON_DATA%INFIL.EQ.'DONE') THEN
         IQUIT = 1
         GO TO 100
      END IF
!
!     MAKE SURE INPUT FILE EXISTS
!
      INQUIRE(FILE=FIZCON_DATA%INFIL,EXIST=IEXIST)
      IF(.NOT.IEXIST)  THEN
         IF(IMDC.LT.4) THEN
            WRITE(IOUT,'(/A/)')  '       COULD NOT FIND INPUT FILE'
         END IF
         SELECT CASE (IMDC)
            CASE (1,2,3)
               IF(IONEPASS.EQ.0) GO TO 10
         END SELECT
         IQUIT = 1
         FIZCON_SUCCESS = 1
         GO TO 100
      END IF
!
!     GET OUTPUT FILE SPECIFICATION
!
      IF(IMDC.LT.4) THEN
         IF(FIZCON_DATA%OUTFIL.EQ.'*' ) THEN
            IF(IMDC.NE.0) THEN
               WRITE(IOUT,FMT=TFMT)                                     &
     &           ' Output Message File Specification    - '
            END IF
            READ(NIN,'(A)') FIZCON_DATA%OUTFIL
         ELSE
            WRITE(IOUT,'(/2A)') ' Output file - ',                      &
     &              TRIM(FIZCON_DATA%OUTFIL)
         END IF
      END IF
      IF(FIZCON_DATA%OUTFIL.NE.' ')  THEN
         NOUT = JOUT             ! SETS FORTRAN OUTPUT UNIT IF DISK FILE
      END IF
!
!     CHECK FOR STANDARD OPTIONS
!
      IF(IW.EQ.'*') THEN
         IF(IMDC.GE.1.AND.IMDC.LE.3) THEN
   15       WRITE(IOUT,FMT=TFMT)  ' Standard Options (Y(es),N(o),?)?  '
            READ(NIN,'(A)')  IW
            IC = ICHAR(IW)
            IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
            IF(IW.EQ.'?')  THEN
               IW = '*'
               WRITE(IOUT,20)
   20          FORMAT(10X,' STANDARD OPTIONS ARE'/                      &
     &             10X,'    CHECK ENTIRE TAPE'/                         &
     &             10X,'    OMIT DEVIANT POINT CHECK'/                  &
     &             10X,'    OMIT SUM UP TESTS    '/                     &
     &             10X,'    EPSILON = .001       ')
               GO TO 15
            END IF
         END IF
      END IF
!
!     GET USER OPTION CHOICE WHEN NOT STANDARD
!
      IF(IMDC.EQ.0.OR.(IW.EQ.'N'.AND.IMDC.LT.4)) THEN
!
!        MATERIAL NUMBER RANGE SELECTION
!
         CALL SELECT_MATS(MATSIN)
!
!        DEVIANT POINT TEST?
!
         IF(IMDC.EQ.0) THEN
            CALL TOKEN(MATSIN,',',3,BUF1)
            IW = BUF1(1:1)
         ELSE
            WRITE(IOUT,FMT=TFMT)                                        &
     &                '     Deviant Point Check (Y(es),N(o))?   -  '
         END IF
         READ(NIN,'(A)') IW
         IC = ICHAR(IW)
         IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
         IF(IW.EQ.'Y')   FIZCON_DATA%ICKT = 0
!
!        SUM UP TESTS?
!
         IF(IMDC.EQ.0) THEN
            CALL TOKEN(MATSIN,',',4,BUF1)
            IW = BUF1(1:1)
         ELSE
            WRITE(IOUT,FMT=TFMT)                                        &
     &                '     Sum Up Tests (Y(es),N(o))?   - '
         END IF
         READ(NIN,'(A)') IW
         IC = ICHAR(IW)
         IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
!
!        SUM UP TESTS SELECTED, GET THE EPSILON TOLERANCE
!
         IF(IW.EQ.'Y')   THEN
            FIZCON_DATA%ISUM = 1
            IF(IMDC.EQ.0) THEN
               CALL TOKEN(MATSIN,',',5,BUF2)
               READ(BUF2,'(BN,E12.5)',ERR=45) EPS
               GO TO 50
   45          EPS = 0.0
            ELSE
               WRITE(IOUT,FMT=TFMT)  '     Enter Epsilon    - '
               READ(NIN,'(E12.5)',ERR=50)  EPS
            END IF
   50       IF(EPS.EQ.0.)   EPS = DEFAULT_EPSILN
            FIZCON_DATA%EPSILN = EPS
         END IF
      END IF
!
!     OPEN INPUT AND OUTPUT FILES
!
      OPEN(UNIT=JIN,ACCESS='SEQUENTIAL',STATUS='OLD',                   &
     &                FILE=FIZCON_DATA%INFIL,ACTION='READ')
      IF(NOUT.NE.6) THEN
!+++MDC+++
!...VMS
!/         OPEN(UNIT=NOUT,ACCESS='SEQUENTIAL',STATUS=OSTATUS,           &
!/     &       FILE=FIZCON_DATA%OUTFIL,CARRIAGECONTROL='LIST')
!...WIN, DVF, UNX, LWI, ANS, MOD
         OPEN(UNIT=NOUT,ACCESS='SEQUENTIAL',STATUS=OSTATUS,             &
     &       FILE=FIZCON_DATA%OUTFIL)
!---MDC---
      END IF
!
!     OPEN SCRATCH FILES
!
      OPEN(UNIT=ISCRX,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(UNIT=ISCRY,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(UNIT=ISCRXY,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(UNIT=ISCRU1,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(UNIT=ISCRU2,FORM='UNFORMATTED',STATUS='SCRATCH')
!
!     OUTPUT SELECTED OPTIONS
!
      IF(IMDC.LT.4)  WRITE(IOUT,'(/A)') ' '
      CALL DATE(ADATE)
      IF(NOUT.NE.IOUT) THEN
         WRITE(NOUT,'(A///2A,30X,2A/)') CHAR(12),                       &
     &            'PROGRAM FIZCON VERSION ',VERSION,                    &
     &            'Run on ',ADATE
      END IF
      WRITE(NOUT,'(2A)')                                                &
     &   'Input File Specification------------------------',            &
     &   TRIM(FIZCON_DATA%INFIL)
      IF(FIZCON_DATA%MATMIN.EQ.0.AND.FIZCON_DATA%MATMAX.EQ.0)   THEN
         WRITE(NOUT,'(A)')  'Check the Entire File'
      ELSE
         WRITE(NOUT,'(A,I4,A,I4)')                                      &
     &        'Check Materials---------------------------------',       &
     &             FIZCON_DATA%MATMIN,' to ',FIZCON_DATA%MATMAX
      END IF
      IF(FIZCON_DATA%ISUM.EQ.1)   THEN
         WRITE(NOUT,'(A)')  'Sum Up Tests will be Performed'
         WRITE(NOUT,'(A,F8.5)') '  Fractional Difference Allowed '//    &
     &             'is ',FIZCON_DATA%EPSILN
      ELSE
         WRITE(NOUT,'(A)')  'Sum Up Tests will be Omitted'
      END IF
      IF(FIZCON_DATA%ICKT.EQ.0)   THEN
         WRITE(NOUT,'(A)')  'Deviant Point Check will be Performed'
         WRITE(NOUT,'(A)')                                              &
     &       'Consecutive Equal Value Check will be Performed'
      ELSE
         WRITE(NOUT,'(A)')  'Deviant Point Check will be Omitted'
         WRITE(NOUT,'(A)')                                              &
     &       'Consecutive Equal Value Check will be Omitted'
      END IF
!
  100 RETURN
      END SUBROUTINE BEGIN
!
!***********************************************************************
!
      SUBROUTINE SELECT_MATS(MATSIN)
!
!     SUBROUTINE GET MATERIALS TO BE EXTRACTED FROM INPUT
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: MATSIN
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     INTEGER(KIND=I4), INTRINSIC :: INDEX, LEN_TRIM
!
      CHARACTER(LEN=10) :: BUF
      CHARACTER(LEN=4) :: BUF1,BUF2
      INTEGER(KIND=I4) :: IDASH
      INTEGER(KIND=I4) :: LBUF
!
!     GET THE USER INPUT
!
      WRITE(IOUT,'(A)') ' '
      WRITE(IOUT,FMT=TFMT) '     Enter Range of MAT Numbers - '
      READ(NIN,'(A)')  MATSIN
!
!     BLANK RESPONSE IS THE SAME AS SELECTING ALL
!
      IF(MATSIN.EQ.' ')  THEN
         FIZCON_DATA%MATMIN = 0
         FIZCON_DATA%MATMAX = 0
         GO TO 100
      END IF
!
!     ANALYZE THE USER INPUT
!
      CALL TOKEN(MATSIN,',',1,BUF)
      IDASH = INDEX(BUF,'-')
      IF(IDASH.GT.0) THEN
         LBUF = LEN_TRIM(BUF)
         IF(IDASH.EQ.1) THEN
            BUF1 = '   1'
            BUF2 = BUF(2:)
         ELSE IF(IDASH.EQ.LBUF) THEN
            BUF2 = '9999'
            BUF1 = BUF(1:LBUF-1)
         ELSE
            BUF1 = BUF(1:IDASH-1)
            BUF2 = BUF(IDASH+1:)
         END IF
      ELSE
         BUF1 = BUF
         CALL TOKEN(MATSIN,',',2,BUF2)
      END IF
!
!     CONVERT FROM ASCII
!
      FIZCON_DATA%MATMIN = 1
      FIZCON_DATA%MATMAX = 9999
      READ(BUF1,'(BN,I4)',ERR=20) FIZCON_DATA%MATMIN
   20 READ(BUF2,'(BN,I4)',ERR=25) FIZCON_DATA%MATMAX
!
!     SET THE MATERIAL NUMBER LIMITS
!
   25 IF(FIZCON_DATA%MATMIN.LE.0) THEN
         FIZCON_DATA%MATMIN = 1
      END IF
      IF(FIZCON_DATA%MATMAX.LT.FIZCON_DATA%MATMIN)  THEN
         FIZCON_DATA%MATMAX = FIZCON_DATA%MATMIN
      END IF
      IF(FIZCON_DATA%MATMIN.EQ.1.AND.FIZCON_DATA%MATMAX.EQ.9999) THEN
         FIZCON_DATA%MATMIN = 0
         FIZCON_DATA%MATMAX = 0
      END IF
!
  100 RETURN
      END SUBROUTINE SELECT_MATS
!
!***********************************************************************
!
      SUBROUTINE SEARCH(IFIND)
!
!     ROUTINE TO CHECK TAPE LABEL AND SEARCH FOR STARTING RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IFIND   ! FLAG IF FIRST DESIRED MATERIAL IS FOUND
!
      CHARACTER(LEN=80) :: IFIELD
!
!     INITIALIZE TO NOT FOUND
!
      IFIND = 0
!
!     READ AND PARSE FIRST CARD TO SEE IF IT IS A LABEL
!
      READ(JIN,'(A)',END=90) IFIELD
!
!     PARSE FIRST CARD TO SEE IF IT IS SEQUENCED
!
      ASEQ = IFIELD(76:80)
!
      READ(IFIELD,'(A,I4,I2,I3,I5)',ERR=20)  TLABEL,MAT,MF,MT,NSEQ
!
!     A LABELED TAPE?
!
      IF(MF.NE.0.OR.MT.NE.0)  THEN
         TLABEL = 'TAPE IS NOT LABELED'
         LABEL = 0
         WRITE(NOUT,'(/A/)') 'TAPE BEING PROCESSED IS NOT LABELED'
         GO TO 60
      ELSE
         LABEL = MAT
         WRITE(NOUT,'(/2A,I5/3X,2A)') 'TAPE BEING PROCESSED IS ',       &
     &         'NUMBERED',LABEL,'LABEL IS  ',TLABEL
      END IF
      GO TO 40
!
!     IF READING ERROR ASSUME A PROPER LABEL AND GO ON
!
   20 WRITE (NOUT,'(6X,A//)')                                           &
     &        'FORMAT ERROR IN FIRST RECORD, PROPER LABEL ASSUMED'
      TLABEL = 'LABEL RECORD IS NOT READABLE'
      LABEL = 0
!
!     READ NEXT CARD
!
   40 READ(JIN,'(A)',END=90) IFIELD
      READ(IFIELD,'(66X,I4,I2,I3,I5)',ERR=50) MAT,MF,MT,NSEQ
!
!     MT=0, FOUND ANOTHER LABEL
!
   50 IF(MT.EQ.0.AND.MF.EQ.0)   THEN
         WRITE(NOUT,'(36X,A)')  'TAPE HAS TOO MANY LABELS'
         LABEL = MAT
         GO TO 40
      END IF
!
!     LOOK FOR BEGINNING OF FIRST MATERIAL REQUESTED
!
   60 IF(FIZCON_DATA%MATMIN.GT.0)   THEN
         DO WHILE(MAT.LT.FIZCON_DATA%MATMIN)
            READ(JIN,'(A)',END=90)  IFIELD
            READ(IFIELD,'(66X,I4,I2,I3,I5)',ERR=65) MAT,MF,MT,NSEQ
   65       IF(MAT.LT.0) GO TO 70
         END DO
         IF(MAT.GT.FIZCON_DATA%MATMAX) GO TO 70
      END IF
      GO TO 75
!
!     FAILED TO FIND A MATERIAL
!
   70 IF(FIZCON_DATA%MATMIN.EQ.FIZCON_DATA%MATMAX) THEN
         IF(FIZCON_DATA%MATMIN.EQ.0) THEN
            EMESS = 'INPUT FILE DOES NOT CONTAIN ANY ENDF EVALUATIONS'
         ELSE
            WRITE(EMESS,'(A,I5)')                                       &
     &           'INPUT FILE DOES NOT CONTAIN MATERIAL',                &
     &           FIZCON_DATA%MATMIN
         END IF
      ELSE
         WRITE(EMESS,'(A,I5,A,I5)')                                     &
     &        'INPUT FILE DOES NOT CONTAIN ANY MATERIALS',              &
     &         FIZCON_DATA%MATMIN,' TO',FIZCON_DATA%MATMAX
      END IF
      WRITE(NOUT,'(/A)')  EMESS
      IF(NOUT.NE.IOUT) THEN
         IF(IMDC.LT.4) WRITE(IOUT,'(10X,A)')  EMESS
      END IF
      GO TO 100
!
!     FOUND BEGINNING OF FIRST MATERIAL REQUESTED
!
   75 READ(IFIELD,'(2E11.4,4I11)')  C1H,C2H,L1H,L2H,N1H,N2H
      IFIND = 1
      NSEQP = NSEQ
      GO TO 100
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE SEARCH
!
!***********************************************************************
!
      SUBROUTINE CHKSEC
!
!     CONTROLS CHECK OF A SECTION BASED ON ITS FILE NUMBER (MF)
!
      IMPLICIT NONE
!
      SELECT CASE (MF)   ! BRANCH BASE ON FILE
         CASE (1)
            CALL CKF1
         CASE (2)
            CALL CKF2
         CASE (3)
            CALL CKF3
         CASE (4)
            CALL CKF4
         CASE (5)
            CALL CKF5
         CASE (6)
            CALL CKF6
         CASE (7)
            CALL CKF7
         CASE (8)
            CALL CKF8
         CASE (9,10)
            CALL CKF9
         CASE (12,13)
            CALL CKF12
         CASE (14)
            CALL CKF14
         CASE (15)
            CALL CKF15
         CASE (23)
            CALL CKF23
         CASE (26)
            CALL CKF26
         CASE (27)
            CALL CKF27
         CASE (28)
            CALL CKF28
         CASE (32)
            CALL CKF32
         CASE (31,33)
            CALL CKF33
         CASE (34)
            CALL CKF34
         CASE (35)
            CALL CKF35
         CASE (40)
            CALL CKF40
         CASE DEFAULT
            IERX = 1
            WRITE(EMESS,'(A,I3,A)') 'MF= ',MF,' IS NOT PERMITTED'
            CALL ERROR_MESSAGE(0)
            NERROR = NERROR + 1
      END SELECT
!
      RETURN
      END SUBROUTINE CHKSEC
!
!***********************************************************************
!
      SUBROUTINE CKF1
!
!     CHECK FILE 1 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IENT
!
!     TEST THAT SECTION IS IN THE INDEX
!
      IF(MT.NE.451)  CALL TESTD(1000*MF+MT)
!
!     BRANCH ON MT NUMBER
!
      SELECT CASE (MT)
!
         CASE (451)            ! COMMENTS AND DIRECTORY
            CALL CKS451
!
         CASE (452)            ! TOTAL NU BAR
            IENT = 1
            CALL CHKNUB(IENT)
!
         CASE (455)            ! DELAYED NUBAR
            IENT = 2
            CALL CHKNUB(IENT)
!
         CASE (456)            ! PROMPT NUBAR
            IENT = 3
            CALL CHKNUB(IENT)
!
         CASE (458)            !ENERGY RELEASE IN FISSION
            CALL CKS458
!
         CASE (460)            !Delayed photon data
            CALL CKS460
!
      END SELECT
!
      RETURN
      END SUBROUTINE CKF1
!
!***********************************************************************
!
      SUBROUTINE CKS451
!
!     CHECK SECTION 451 DATA
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: AMOD, FLOAT, AINT
!
      INTEGER(KIND=I4) :: IZT
      INTEGER(KIND=I4) :: JPART
      INTEGER(KIND=I4) :: NCD,NID
      INTEGER(KIND=I4) :: N1,NN,NC,N
      REAL(KIND=R4) :: ZTRY,ASAV
      REAL(KIND=R4) :: ZT,ELISM
!
!     INITIALIZE
!
      ZSA ='           '
      ZA = C1H
      AWR = C2H
      LRP = L1H
      LFI = L2H
      NLIB = N1H
      NMOD = N2H
      IZT = IFIX(ZA)/1000
      ZT = IZT
!
!     TEST CHARGE-MASS REASONABILITY
!
      IF(MAT.LT.100) THEN
         ZTRY = FLOAT(MAT) + 100.
         CALL TEST3F(ZA,ZTRY,'ZA')
      ELSE
         ASAV = AMOD(ZA,1000.)
         IF(ASAV.GT.0.) THEN
            ZTRY = AINT((ASAV+1.)/2.) + 1.
            CALL TEST6(ZT,1.,ZTRY,'Z')
         END IF
      END IF
!
!     READ THE NEXT CONTROL RECORD AND SET PARAMETERS
!
      CALL RDCONT
      ELIS = C1H
      STA = C2H
      LIS = L1H
      LISO = L2H
      NFOR = N2H
      IF(LIS.NE.0.AND.ELIS.EQ.0.0)  THEN
         EMESS = 'ELIS SHOULD NOT BE ZERO FOR A METASTABLE STATE'
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!
!     ENDF-V FORMAT FILE
!
      IF(NFOR.EQ.0)   THEN
         NFOR = 5
         IF((NLIB.GE.2.AND.NLIB.LE.4).OR.NLIB.EQ.35)  THEN
            NVER = 1
         ELSE IF(NLIB.EQ.5)   THEN
            NVER = 2
         ELSE IF(NLIB.EQ.6)   THEN
            NVER = 3
         ELSE
            NVER = 5
         END IF
         ENMAX = 2.0E+7
         NSUB = 10
         AWI = 1.
      ELSE
!
!     ENDF-VI OR LATER FORMAT, READ ANOTHER CONTROL RECORD
!
         CALL RDCONT
         AWI = C1H
         ENMAX = C2H
         NSUB = N1H
         NVER = N2H
         NFOR = MAX0(6,NFOR)
      END IF
!
!     IS TARGET EXCITATION ENERGY REASONABLE?
!
      IF(LIS.EQ.0)  THEN
         ELISM = 0.
      ELSE
         IF(NSUB.GE.10) THEN
            ELISM = ENMAX
         ELSE
            ELISM = 3.0E+6
         END IF
      END IF
      CALL TEST6(ELIS,0.0,ELISM,'ELIS')
!
!     CHECK FOR CORRECT AWI VALUE
!
      JPART = NSUB/10
      DO NN=1,NPARTS
         IF(JPART.EQ.IPARTS(NN))   THEN
            CALL TEST3F(AWI,AWPART(NN),'AWI')
            GO TO 10
         END IF
      END DO
      IF(JPART.EQ.11) THEN
         EMESS = 'AWI TEST BYPASSED FOR ELECTRONS.'
!        Skip warning message and continue
         GO TO 10
       ELSE
         EMESS = 'AWI TEST BYPASSED FOR PARTICLE MASS GREATER THAN 4.'
      END IF
      CALL WARNING_MESSAGE(1)
!
!     PROCESS LAST CONTROL RECORD
!
   10 CALL RDCONT
      LDRV= L1H
!
!     READ IN COMMENT RECORDS
!
      NCD = N1H
      IF(NFOR.GE.6)  THEN
         NID = 5
      ELSE
         NID = 2
      END IF
      DO NC=1,NCD
         CALL RDTEXT
         IF(NC.EQ.1) ZSA=TEXT(1:11)
         IF(NC.LE.NID)   THEN
            IF(IMDC.LT.4) WRITE(IOUT,'(1X,A66)')   TEXT
         END IF
      END DO
!
!     SET DECAY OPTIONS AT FIRST MATERIAL ON A DECAY DATA TAPE
!
!      IF(NSUB.EQ.4.AND.IDDONE.EQ.0) THEN
!         CALL SETDCHK
!         IDDONE = 1
!      END IF
!
!     PROCESS DIRECTORY
!
      N1 = 0
      NXC = N2H
      DO N=1,NXC
         CALL RDCONT
         IF(L1H.EQ.9.OR.L1H.EQ.10)  THEN
            NISSEC = NISSEC + 1
!
            IF(NISSEC.GT.SZLMF)
     &         STOP 'FIZCON ERROR - SZLMF limit exceeded'
!
            MTISO(NISSEC) = L2H
         END IF
         N1 = N1 + 1
         INDX(N1,1) = 1000*L1H + L2H
         INDX(N1,2) = 1
         ENGS(N1,1) = 0.
         ENGS(N1,2) = 0.
      END DO
!
!     MAKE SURE SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     INITIALIZE FOR FISSION ENERGY RELEASE TEST
!
      MT458 = 0
      ERQ = 0.0
!
!     INITIALIZE FOR NUBAR SUMUP TEST
!
      IF(LFI.EQ.1.AND.FIZCON_DATA%ISUM.EQ.1) CALL SUM452(-1)
!
      RETURN
      END SUBROUTINE CKS451
!
!***********************************************************************
!
      SUBROUTINE CHKNUB(IENT)
!
!     CHECK NUBAR SECTIONS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IENT
!
      INTEGER(KIND=I4) :: NN,N
      REAL(KIND=R4) :: E,VNU
      REAL(KIND=R4) :: ELO,EHI
!
      INTEGER(KIND=I4), PARAMETER :: NNUS=3
      CHARACTER(LEN=4), DIMENSION(NNUS), PARAMETER ::                   &
     &         KNU = (/'NU  ','NUD ','NUP '/)
! Upper limit for the values of nu-bar
      REAL(KIND=R4), DIMENSION(NNUS), PARAMETER ::                       &
     &         UPR = (/20.0,1.0,20.0/)
!
      IF(IENT.EQ.2) THEN
!********READ DECAY CONSTANTS
         IF(L1H.EQ.0) THEN
            CALL RDLIST
            CALL TEST5Y(1,NPL,1,1)
         ELSE
            CALL RDTAB2
            DO N=1,NP2
               CALL RDLIST
            END DO
         END IF
!********CHECK IF LAMBDA-S ARE IN INCREASING ORDER
      END IF
!
!     PROCESS NU BAR
!
      IF(L2H.EQ.1)   THEN
!*****POLYNOMIAL REPRESENTATION
         CALL RDLIST
         VNU = 0.0
         E = 1.0
         DO NN=1,NPL
            VNU = VNU + Y(NN)*E
            E = E*ENMAX
         END DO
         CALL TEST6(VNU,0.0,UPR(IENT),KNU(IENT))
         ELO = ENMIN
         EHI = ENMAX
      ELSE
!*****TABULAR REPRESENTATION
         CALL RDTAB1
         CALL TEST6Y(0.0,UPR(IENT),KNU(IENT))
         ELO = X(1)
         EHI = X(NP)
      END IF
!
!     DO SUMUP TEST
!
      IF(ITEST.EQ.1)   CALL SUM452(MT)
!
!     STORE ENERGY SPAN OF THE SECTION
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     CHECK ENERGY RANGE OF NU BAR SECTIONS
!
      IF(MT.EQ.452)  THEN
         CALL TESTER(ELO,EHI,QUNK)
      ELSE
         CALL ISFIL(MF,MF,MT,452)
      END IF
!
      RETURN
      END SUBROUTINE CHKNUB
!
!***********************************************************************
!
      SUBROUTINE CKS458
!
!     CHECK ENERGY RELEASE PER FISSION
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: K,KK
      REAL(KIND=R4) :: YK,YKK
      REAL(KIND=R4) :: SSUM
      REAL(KIND=R4) :: ET,DELTA,ERBAR
      INTEGER(KIND=I4) :: N, LFC, NFC
!
      MT458 = 1
      IF(LFI.NE.1)  THEN
         EMESS = 'SECTION SHOULD BE USED FOR FERTILE AND '//            &
     &           'FISSIONABLE ISOTOPES ONLY'
         CALL ERROR_MESSAGE(0)
      END IF
      NFC = N2H
      LFC = L2H
      CALL RDLIST
!*****SUM PARTIAL ENERGIES AND CHECK VALUES
      SSUM = 0.0
      DO K=1,NPL,2
         YK = Y(K)
         IF(K.LE.13) SSUM = SSUM + YK
         KK = MOD(K,18)
         YKK = Y(KK)
         IF(ABS(YK).LT.Y(K+1)) THEN
            WRITE(EMESS,'(A,I3)')                                       &
     &                'ERROR GREATER THAN VALUE AT COMPONENT #',K
            CALL ERROR_MESSAGE(NSEQP)
         END IF
         IF(Y(K+1).LT.0.0)   THEN
            WRITE(EMESS,'(A,I3)')                                       &
     &                'NEGATIVE FISSION ENERGY UNCERTAINTY #',K+1
            CALL ERROR_MESSAGE(NSEQP)
         END IF
         IF(K.LE.18.AND.YK.LT.0.0)   THEN
            WRITE(EMESS,'(A,I3)')                                       &
     &                'NEGATIVE FISSION ENERGY COMPONENT #',K
            CALL ERROR_MESSAGE(NSEQP)
         END IF
         IF(YKK.LT.ABS(YK)) THEN
            WRITE(EMESS,'(A,I3)')                                       &
     &                'ABS(COMPONENT) GREATER THAN C0 VALUE AT #',K
            CALL ERROR_MESSAGE(NSEQP)
         END IF
      END DO
!*****TEST SUMS
      ERQ = Y(15)
      ET = Y(17)
      DELTA = ABS(ET-SSUM)/ET
      IF(DELTA.GT.FIZCON_DATA%EPSILN)   THEN
         WRITE(EMESS,'(A,1PE12.5,A,1PE12.5)')                           &
     &             'TOTAL ENERGY RELEASE PER FISSION=',ET,              &
     &             '  SUM OF PARTIALS=',SSUM
         CALL ERROR_MESSAGE(0)
      END IF
      ERBAR = Y(13)
      DELTA = ABS(SSUM-ERBAR-ERQ)/ERQ
      IF(DELTA.GT.FIZCON_DATA%EPSILN)  THEN
         WRITE(EMESS,'(A,1PE12.5,A,1PE12.5,A)')                         &
     &            'TOTAL ENERGY (',SSUM,') LESS NEUTRINO ENERGY (',     &
     &            ERBAR,')'
         CALL ERROR_MESSAGE(0)
         WRITE(EMESS,'(A,1PE12.5,A)')                                   &
     &            '    DOES NOT EQUAL RELEASE (',ERQ,')'
         CALL ERROR_MESSAGE(0)
      END IF
!
      IF (LFC.EQ.1) THEN ! Fission energy release as tables, just skipping over, not checking anything
        DO N=1,NFC
          CALL RDTAB1
        END DO
      END IF

      RETURN
      END SUBROUTINE CKS458
!
!***********************************************************************
!
      SUBROUTINE CKS460
!
!     Ccheck delayed photon data
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LO,NG,I
!
      IF(LFI.NE.1)  THEN
         EMESS = 'SECTION SHOULD BE USED FOR FERTILE AND '//            &
     &           'FISSIONABLE ISOTOPES ONLY'
         CALL ERROR_MESSAGE(0)
      END IF
!
      LO = L1H
      NG = N1H
!
!     Discrete representation
!
      IF(LO.EQ.1) THEN
        DO I=1,NG
          CALL RDTAB1
        END DO
!
!     Continuous representation
!
      ELSE IF(LO.EQ.2) THEN
        CALL RDLIST
      END IF
!
      RETURN
      END SUBROUTINE CKS460
!
!***********************************************************************
!
      SUBROUTINE CKF2
!
!     CHECK FILE 2 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: IFIX
      REAL(KIND=R4), INTRINSIC :: ABS, AMOD
!
      INTEGER(KIND=I4) :: NIS
      INTEGER(KIND=I4) :: INAT,IZI,IZH
      INTEGER(KIND=I4) :: LRU,LRF,LFW
      INTEGER(KIND=I4) :: NUMSQ1
      INTEGER(KIND=I4) :: NE,NI,NER
      REAL(KIND=R4) :: ZAH,AWRH,ABNTOT,SPI
      REAL(KIND=R4) :: AWRIT
      REAL(KIND=R4) :: ABNM
      REAL(KIND=R4) :: EL,EH,EUBN1,ELB,EUB
      REAL(KIND=R4) :: DELTA,DELTAL,DELTAU
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     STORE VARIABLES FROM HEAD RECORD
!
      NIS = N1H
      ZAH = C1H
      AWRH = C2H
!
!     ONLY ONE ISOTOPE FOR IF NOT A NATURAL ELEMENT
!
      IF(AMOD(ZAH,1000.).EQ.0.0)    THEN
         INAT = 1
      ELSE
         INAT = 0
         CALL TEST3(NIS,1,'NIS')
      END IF
!
!     LOOP ON ALL ISOTOPES
!
      ABNTOT = 0.0
      IZH = IFIX(ZAH)/1000
      DO NI=1,NIS
         CALL RDCONT
         ABNTOT = ABNTOT + C2H
!
!        SET LIMITS ON AWRI
!
         IZI = IFIX(C1H)/1000
         AWRIT = (C1H-1000.*FLOAT(IZI))/FACTOR
         AWRI1 = AWRIT + 1.
         AWRI2 = AWRIT - 1.
!
!        CHECKS FOR NATURAL ELEMENT
!
         IF(INAT.EQ.1)  THEN
            ABNM = 0.0
            IF(IZI.NE.IZH)  THEN
               EMESS = 'ISOTOPE Z SHOULD EQUAL MATERIAL Z'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
            IF(C2H.EQ.0.0)    THEN
               EMESS = 'ISOTOPE ABUNDANCE CANNOT BE 0.0.'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
!
!        CHECKS FOR SINGLE ISOTOPE
!
         ELSE
            ABNM = 1.
            IF(C1H.NE.ZA)  THEN
               EMESS = 'ISOTOPE ZA SHOULD EQUAL MATERIAL ZA'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
         END IF
!
!        CHECK ABUNDANCE
!
         CALL TEST6(C2H,ABNM,1.0,'ABN')
!
!        INITIALIZE FOR ISOTOPE
!
         LFW = L2H
         NER = N1H
         E1 = BIGNO
         E2 = 0.
!
!        LOOP ON ENERGY RANGES
!
         DO NE=1,NER
            CALL RDCONT
            EL = C1H
            EH = C2H
            LRU = L1H
            LRF = L2H
            NRO = N1H
            E1 = AMIN1(E1,EL)
            E2 = AMAX1(E2,EH)
            IF(EL.GT.EH)   THEN
               EMESS = 'ENERGY RANGE LIMITS WRONG'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
!
!           CHECK FOR CONTINUITY OF REGION BOUNDARIES
!
            IF(NE.NE.1)  THEN
               DELTA = ABS(1.-EL/EUBN1)
               IF(DELTA.GT.EPSILN5)   THEN
                  EMESS = 'RESONANCE ENERGY RANGE NOT CONTINUOUS'
                  CALL ERROR_MESSAGE(NSEQP1)
                  WRITE(EMESS,'(4X,A,I6)')                              &
     &                'PREVIOUS RANGE DEFINED AT RECORD',NUMSQ1
                  CALL ERROR_MESSAGE(0)
               END IF
            END IF
            EUBN1 = EH
            NUMSQ1 = NSEQP1
!
!           PROCESS EACH DIFFERENT RESONANCE REGION REPRESENTATION
!
            IF(LRU.EQ.0)   THEN
               CALL RDCONT
               SPI = C1H
               CALL TESTSP(SPI)
            ELSE IF(LRU.EQ.1) THEN
               IF((LRF.GE.1.AND.LRF.LE.3).OR.LRF.EQ.5) THEN
                  CALL CHKBW(LRF)
               ELSE IF(LRF.EQ.4) THEN
                  CALL CHKAA
               ELSE IF(LRF.EQ.6) THEN
                 CALL CHKHR
               ELSE IF(LRF.EQ.7) THEN
                 CALL CHKRL
               END IF
            ELSE IF(LRU.EQ.2)   THEN
               CALL CHKUR(LRF,LFW)
            END IF
         END DO
!
!        CHECK THAT ISOTOPES SPAN SAME ENERGY RANGE
!
         IF(NI.EQ.1)   THEN
            ELB = E1
            EUB = E2
         ELSE
            DELTAL = ABS(1.-E1/ELB)
            DELTAU = ABS(1.0-E2/EUB)
            IF(DELTAL.GT.EPSILN5.OR.DELTAU.GT.EPSILN5)  THEN
               WRITE(EMESS,'(A,I2,A,1PE12.5,A,1PE12.5,A)')              &
     &              'ISOTOPE #',NI,' ENERGY RANGE(',E1,' TO ',E2,')'
               CALL ERROR_MESSAGE(0)
               WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5,A)')                &
     &               'DIFFERS FROM THE FIRST(',ELB,' TO ',EUB,')'
               CALL ERROR_MESSAGE(0)
            END IF
         END IF
      END DO
!
!     TEST THAT ABUNDANCES ADD UP TO ONE
!
      IF((ABS(ABNTOT-1.)).GT.EPSILN3)   THEN
         EMESS = 'ISOTOPIC ABUNDANCES DO NOT ADD UP TO UNITY'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
      RETURN
      END SUBROUTINE CKF2
!
!***********************************************************************
!
      SUBROUTINE CHKBW(LRF)
!
!     CHECK BREIT-WIGNER, REICH-MOORE, AND R-MATRIX REPRESENTATION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LRF
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: NREP
      INTEGER(KIND=I4) :: NLS
      INTEGER(KIND=I4) :: IGN,JGN
      INTEGER(KIND=I4) :: ISEQ
      INTEGER(KIND=I4) :: NL,I
      REAL(KIND=R4)  :: AWRIL
      REAL(KIND=R4) :: FL,AJLO,AJHI,AJ
      REAL(KIND=R4) :: SPI,GN
!
!     READ AND TEST ENERGY DEPENDENT SCATTERING LENGTH
!
      IF(NRO.NE.0)   THEN
         CALL RDTAB1
      END IF
!
!     CHECK SPIN AND ENERGY INDEPENDENT SCATTERING LENGTH
!
      CALL RDCONT
      NLS = N1H
      SPI = C1H
      CALL TESTSP(SPI)
!
!     PROCESS PARAMETERS FOR ALL L VALUES
!
      DO NL=1,NLS
         CALL RDLIST
         NREP = NPL/N2L
!******* TEST AWRI
         IF(NL.EQ.1)  THEN
            CALL TEST6(C1L,AWRI2,AWRI1,'AWR')
            AWRIL = C1L
         ELSE
            CALL TEST3F(C1L,AWRIL,'AWR')
         END IF
!********CHECK FOR CORRECT L-VALUE FOR THIS SUBSECTION
         CALL TEST3(L1L,NL-1,'L')
!********GET RANGE OF VALID J-VALUES
         FL = L1L
         AJLO = ABS(ABS(SPI-FL)-0.5)
         AJHI = SPI+FL+0.5
!********TEST THAT RESONANCE ENERGIES ARE IN INCREASING ORDER
         CALL TEST5Y(1,NPL,NREP,1)
!********TEST IF PARTIAL WIDTHS ADD UP TO TOTAL
         IF(LRF.NE.3)   CALL TESTW(NPL,3,NREP)
!
!        TEST ON INDIVIDUAL PARAMETERS
!
         DO I=4,NPL,NREP
!***********POSSIBLE J-VALUE?
            AJ = Y(I-2)
            IF(AJ.LT.0.0.AND.LRF.EQ.3) THEN
               AJ = - AJ
               IF(FL.EQ.0.0.AND.SPI.EQ.0.0) THEN
                  ISEQ = NSEQP1 + (I+2)/NREP
                  EMESS = 'AJ CANNOT BE NEGATIVE FOR L AND SPI '//      &
     &                    'EQUAL ZERO '
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(4X,A,1PE12.5)') 'FOR RESONANCE',Y(I-3)
                  CALL ERROR_MESSAGE(ISEQ)
               END IF
            END IF
            CALL TEST6(AJ,AJLO,AJHI,'AJ')
!***********TEST FOR ZERO NEUTRON WIDTH
            IF(LRF.EQ.3) THEN
               IGN = I - 1
            ELSE
               IGN = I
            END IF
            GN = Y(IGN)
            IF(GN.EQ.0.)   THEN
               ISEQ = NSEQP1 + (I+2)/NREP
               WRITE(EMESS,'(A,1PE12.5)')                               &
     &               'NEUTRON WIDTH IS ZERO FOR RESONANCE',Y(I-3)
               CALL ERROR_MESSAGE(ISEQ)
            END IF
!***********TEST FOR NEGATIVE WIDTHS
            DO JGN=IGN,NREP
              GN = Y(JGN)
              IF(GN.LT.0.)   THEN
                 ISEQ = NSEQP1 + (I+2)/NREP
                 WRITE(EMESS,'(A,1PE12.5)')                               &
     &                 'NEGATIVE WIDTHS FOR RESONANCE',Y(I-3)
                 CALL ERROR_MESSAGE(ISEQ)
              END IF
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE CHKBW
!
!***********************************************************************
!
      SUBROUTINE CHKAA
!
!     CHECK ADLER-ADLER RESONANCE REPRESENTATION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NLS,NJS
      INTEGER(KIND=I4) :: NL,NJ
      REAL(KIND=R4) :: SPI
      REAL(KIND=R4) :: AWRI
      REAL(KIND=R4) :: FL,AJLO,AJHI
!
!     READ AND TEST ENERGY DEPENDENT SCATTERING LENGTH
!
      IF(NRO.NE.0)   THEN
         CALL RDTAB1
      END IF
!
!     CHECK SPIN AND ENERGY INDEPENDENT SCATTERING LENGTH
!
      CALL RDCONT
      NLS = N1H
      SPI = C1H
      CALL TESTSP(SPI)
!
!     PROCESS PARAMETERS FOR ALL L VALUES
!
      CALL RDLIST
      AWRI = C1L
      CALL TEST6(AWRI,AWRI2,AWRI1,'AWR')
!
!     PROCESS ALL L VALUES
!
      DO NL=1,NLS
         CALL RDCONT
!********CHECK FOR CORRECT L-VALUE FOR THIS SUBSECTION
         CALL TEST3(L1H,NL-1,'L')
!********GET RANGE OF VALID J-VALUES
         FL = L1H
         AJLO = ABS(ABS(SPI-FL)-0.5)
         AJHI = SPI+FL+0.5
!
!        PROCESS ALL J VALUES
!
         NJS = N1H
         DO NJ=1,NJS
            CALL RDLIST
!**********POSSIBLE J-VALUE?
            CALL TEST6(C1L,AJLO,AJHI,'AJ')
!***********TEST INCREASING ORDER OF RESONANCE ENERGIES
            CALL TEST5Y(1,NPL,12,1)
         END DO
      END DO
!
      RETURN
      END SUBROUTINE CHKAA
!
!***********************************************************************
!
      SUBROUTINE CHKHR
!
!     CHECK HYBRID R-FUNCTION REPRESENTATION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: INT
      REAL(KIND=R4), INTRINSIC :: FLOAT, ABS
!
      INTEGER(KIND=I4) :: NLS,NJS,NSS
      INTEGER(KIND=I4) :: NL,NJ,NS
      INTEGER(KIND=I4), DIMENSION(4) :: MTRE
      INTEGER(KIND=I4) :: NCRE,MTREC
      INTEGER(KIND=I4) :: NAW,NCR
      INTEGER(KIND=I4) :: LBK,LPS,NLSJ,NREP
      INTEGER(KIND=I4) :: ISEQ
      INTEGER(KIND=I4) :: IGN
      INTEGER(KIND=I4) :: I,II,III,LIL
      REAL(KIND=R4) :: SPI
      REAL(KIND=R4) :: QTLOW,QTHIGH
      REAL(KIND=R4) :: AWRI,AWRC
      REAL(KIND=R4) :: FL,FLP,AJLO,AJHI
      REAL(KIND=R4) :: ETST0,ETST
      REAL(KIND=R4) :: GN
      REAL(KIND=R4) :: AS,FAS,AJ,FAJ
      REAL(KIND=R4) :: AL,ALTEST,EPTEST
      REAL(KIND=R4), DIMENSION(4) ::  APART

!
!     READ AND TEST ENERGY DEPENDENT SCATTERING LENGTH
!
      IF(NRO.NE.0)   THEN
         CALL RDTAB1
      END IF
!
!     TEST SPIN
!
      CALL RDCONT
      NLS = N1H
      SPI = C1H
      CALL TESTSP(SPI)
!
!     PROCESS EACH REACTION CHANNEL
!
      CALL RDCONT
      NCRE = N2H
!*****READ REACTION CHANNEL DEFINITIONS
      CALL RDCONT
      NAW = 0
      MTRE = (/L1H,L2H,N1H,N2H/)
      DO II=1,4
         MTREC = MTRE(II)
         IF(MTREC.GT.102)  THEN
            NAW = NAW + 1
            APART(NAW) = AWPART(MTREC-100)
         END IF
      END DO
!
!     CHECK THAT Q VALUES ARE REASONABLE FOR THE CORRESPONDING
!       REACTION
!
      CALL RDLIST
      ETST0 = 0.
      DO III=1,4
         MTREC = MTRE(III)
         IF(MTREC.EQ.102.OR.MTREC.EQ.18.OR.MTREC.EQ.0)  THEN
            IF(Y(III).NE.0.0) THEN
               WRITE(EMESS,'(A,I1,A,I1,A,I3,A)')                        &
     &                'QRE',III,' FOR MTRE',III,' = ',MTREC,            &
     &                ' MUST BE ZERO'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
         ELSE
            QTLOW = -2.0E+7
            QTHIGH = 2.0E+7
            ETST = - Y(III)
            IF(MTREC.GE.51.AND.MTREC.LE.54) THEN
               IF(ETST.LE.0.0)  THEN
                  WRITE(EMESS,'(A,I1,A,I1,A,I3,A)')                     &
     &                  'QRE',III,' FOR MTRE',III,' = ',MTREC,          &
     &                  ' MUST BE ','NEGATIVE'
                  CALL ERROR_MESSAGE(NSEQP)
               END IF
               IF(ETST0.NE.0.0)   THEN
                  IF(ETST.LE.ETST0)  THEN
                     EMESS = 'Q-VALUES FOR INELASTIC CHANNELS '//       &
     &                       'ARE OUT OF ORDER'
                     CALL ERROR_MESSAGE(NSEQP)
                  END IF
               END IF
               ETST0 = ETST
               QTLOW = 1000.
            END IF
            IF(ETST.LT.QTLOW.OR.ETST.GT.QTHIGH) THEN
               WRITE(EMESS,'(A,I1,A,I1,A,I3,A)')                        &
     &              'QRE',III,' FOR MTRE',III,' = ',MTREC,              &
     &              ' IS UNREASONABLE'
!              CALL ERROR_MESSAGE(NSEQP)
            END IF
         END IF
      END DO
!*****READ ANY CHARGED PARTICLE PENETRABILITIES
      IF(NCRE.GT.0)  THEN
         DO NCR=1,NCRE
             DO LIL=1,4
               CALL RDTAB1
               AWRC =  C1
               CALL TEST3F(AWRC,APART(NCR),'AWRC')
             END DO
         END DO
      END IF
!
!     PROCESS EACH L, S, AND J VALUE
!
      FLP = -1.
      DO NL=1,NLS
         CALL RDCONT
         AWRI = C1H
         FL = L1H
!********CHECK FOR CORRECT L-VALUE FOR THIS SUBSECTION
         CALL TEST3(L1L,NL-1,'L')
!********CHECK ORDER OF LISTS
         IF(FL.GT.FLP)  THEN
            FLP = FL
         ELSE
            EMESS = 'RESONANCE PARAMETER LISTS OUT OF ORDER IN L'
            CALL ERROR_MESSAGE(NSEQP)
         END IF
         CALL TEST6(AWRI,AWRI2,AWRI1,'AWR')
         NSS = N1H
         FAS = -1.
!********CHANNEL SPIN
         DO NS=1,NSS
            CALL RDCONT
            AS = C1H
!********CHECK ORDER OF LISTS
            IF(AS.GT.FAS)  THEN
               FAS = AS
            ELSE
               EMESS = 'RESONANCE PARAMETER LISTS OUT OF ORDER IN S'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
!***********POSSIBLE S-VALUE?
            CALL TEST6(AS,SPI-.5,SPI+0.5,'AS')
            NJS = N1H
!***********GET RANGE OF VALID J-VALUES
            AJLO = ABS(FL-AS)
            AJHI = FL+AS
            FAJ = -1.
!***********TOTAL SPIN
            DO NJ=1,NJS
               CALL RDLIST
               AJ = C1L
!**************CHECK ORDER OF LISTS
               IF(AJ.GT.FAJ)  THEN
                  FAJ = AJ
               ELSE
                  EMESS = 'RESONANCE PARAMETER LISTS OUT OF ORDER '//   &
     &                          'IN L'
                  CALL ERROR_MESSAGE(NSEQP)
               END IF
!**************POSSIBLE J-VALUE?
               CALL TEST6(AJ,AJLO,AJHI,'AJ')
               LBK = L1L
               LPS = L2L
               NLSJ = N2L
               NREP = NPL/NLSJ
!**************TEST THAT RESONANCE ENERGIES ARE IN INCREASING ORDER
               CALL TEST5Y(1,NPL,NREP,1)
!
!              TEST ON INDIVIDUAL PARAMETERS
!
               DO I=4,NPL,NREP
                  ISEQ = NSEQP1 + (I+2)/NREP
!*****************TEST FOR ZERO NEUTRON WIDTH
                  IGN = I-2
                  GN = Y(IGN)
                  IF(GN.EQ.0.)   THEN
                     WRITE(EMESS,'(A,1PE12.5)')                         &
     &                    'NEUTRON WIDTH IS ZERO FOR RESONANCE',Y(IGN-1)
                     CALL ERROR_MESSAGE(ISEQ)
                  END IF
!*****************TEST THAT OUTGOING ANGULAR MOMENTUM VALUES ARE
!*****************INTEGRAL  AND REASONABLE
                  DO II=5,8
                     AL = Y(I+II)
                     IF(AL.LT.0.0.OR.AL.GT.3.0) THEN
                        WRITE(EMESS,'(A,I1,A)')                         &
     &                         'ALRE',II-2,' IS NOT ACCEPTABLE'
                        CALL ERROR_MESSAGE(ISEQ)
                     END IF
                     ALTEST = FLOAT(INT(AL))
                     EPTEST = ABS(AL-ALTEST)
                     IF(AL.NE.0.0)   EPTEST = EPTEST/AL
                     IF(EPTEST.GT.EPSILN5)  THEN
                        WRITE(EMESS,'(A,I1,A)')                         &
     &                         'ALRE',II-2,' IS NOT AN INTEGRAL NUMBER'
                        CALL ERROR_MESSAGE(ISEQ)
                     END IF
                  END DO
               END DO
!**************READ BACKGROUND
               IF(LBK.NE.0) THEN
                  CALL RDTAB1
                  CALL RDTAB1
               END IF
!**************READ PHASE SHIFTS
               IF(LPS.NE.0)  THEN
                  CALL RDTAB1
                  CALL RDTAB1
               END IF
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE CHKHR
!
!***********************************************************************
!
      SUBROUTINE CHKRL
!
!     CHECK R-MATRIX-LIMITED REPRESENTATION
!
      IMPLICIT NONE
      INTEGER(KIND=I4) :: J,NJS,NPP,IZA,IA1,IB1
!
!     Read and test energy dependent scattering radius
!
      IF(NRO.NE.0)   THEN
         CALL RDTAB1
      END IF
      CALL RDCONT
      NJS=N1H
      CALL RDLIST
      NPP=L1L
      IZA=ZA
      IA1=0
      IB1=0
!     Checking of the contents of the LIST records ???
!      DO J=1,NPP
!!         Check charge conservation
!          IA1=Y((J-1)*12+3)
!          IB1=Y((J-1)*12+4)
!          CALL TEST3(IA1+IB1,IZA/1000,'Charge sum')
!      END DO
!
      DO J=1,NJS
         CALL RDLIST
         CALL RDLIST
      END DO
!
      RETURN
      END SUBROUTINE CHKRL
!
!***********************************************************************
!
      SUBROUTINE CHKUR(LRF,LFW)
!
!     CHECK UNRESOLVED REGION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LRF,LFW
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: MUF
      INTEGER(KIND=I4) :: NLS,NJS
      INTEGER(KIND=I4) :: N
      INTEGER(KIND=I4) :: NL,NJ,MJ
      REAL(KIND=R4) :: SPI
      REAL(KIND=R4) :: AWRIL
      REAL(KIND=R4) :: FL,AJLO,AJHI,AJ
!
!     READ AND TEST ENERGY DEPENDENT SCATTERING LENGTH
!
      IF(NRO.NE.0)   THEN
         CALL RDTAB1
      END IF
!
!     ALL PARAMETERS ENERGY DEPENDENT
!
      IF(LRF.EQ.2)  THEN
         CALL RDCONT
         SPI = C1H
         CALL TESTSP(SPI)
!
!        PROCESS ALL L VALUES
!
         NLS = N1H
         MF2URL=MAX(MF2URL,NLS)
         DO NL=1,NLS
            CALL RDCONT
!***********TEST AWRI
            IF(NL.EQ.1)   THEN
               CALL TEST6(C1H,AWRI2,AWRI1,'AWR')
               AWRIL = C1H
            ELSE
               CALL TEST3F(C1H,AWRIL,'AWR')
            END IF
!***********CHECK FOR CORRECT L-VALUE FOR THIS SUBSECTION
            CALL TEST3(L1H,NL-1,'L   ')
!***********GET RANGE OF VALID J VALUES
            N = L1H
            FL = L1H
            AJLO = ABS(ABS(SPI-FL)-0.5)
            AJHI = SPI+FL+0.5
!
!           PROCESS ALL J VALUES
!
            NJS = N1H
            DO NJ=1,NJS
               CALL RDLIST
!**************POSSIBLE J- VALUE?
               CALL TEST6(C1L,AJLO,AJHI,'AJ')
!**************TEST AMUX, AMUN, AMUG AND AMUF
               CALL TESTDF(2,N2L)
!**************TEST ENERGY GRID
               CALL TESTE(NPL,6,N,NJ)
            END DO
         END DO
      ELSE
!
!     ALL PARAMETERS ENERGY INDEPENDENT
!
         IF(LFW.EQ.0) THEN
            CALL RDCONT
            SPI = C1H
            CALL TESTSP(SPI)
!
!           PROCESS ALL L VALUES
!
            NLS = N1H
            DO NL=1,NLS
               CALL RDLIST
!**************TEST AWRI
               IF(NL.EQ.1)   THEN
                  CALL TEST6(C1L,AWRI2,AWRI1,'AWR')
                  AWRIL = C1L
               ELSE
                  CALL TEST3F(C1L,AWRIL,'AWR')
               END IF
!**************CHECK FOR CORRECT L-VALUE FOR THIS SUBSECTION
               CALL TEST3(L1L,NL-1,'L')
!**************GET RANGE OF VALID J VALUES
               FL = L1L
               AJLO = ABS(ABS(SPI-FL)-0.5)
               AJHI = SPI+FL+0.5
!**************TEST AMUN
               CALL TESTDF(1,N2L)
!**************TEST J VALUES
               DO MJ=2,NPL,6
                  AJ = Y(MJ)
                  CALL TEST6(AJ,AJLO,AJHI,'AJ')
               END DO
            END DO
         ELSE
!
!     ONLY FISSION WIDTHS ENERGY DEPENDENT
!
            CALL RDLIST
            SPI = C1L
            CALL TESTSP(SPI)
!
!           PROCESS ALL J VALUES
!
            NLS = N2L
            DO NL=1,NLS
               CALL RDCONT
!**************TEST AWRI
               IF(NL.EQ.1)   THEN
                  CALL TEST6(C1H,AWRI2,AWRI1,'AWR')
                  AWRIL = C1H
               ELSE
                  CALL TEST3F(C1H,AWRIL,'AWR')
               END IF
!**************CHECK FOR CORRECT L-VALUE FOR THIS SUBSECTION
               CALL TEST3(L1H,NL-1,'L')
!**************GET RANGE OF VALID J VALUES
               FL = L1H
               AJLO = ABS(ABS(SPI-FL)-0.5)
               AJHI = SPI+FL+0.5
!
!              PROCESS ALL J VALUES
!
               NJS = N1H
               DO NJ=1,NJS
                  CALL RDLIST
!*****************TEST MUF AND AMUN
                  MUF = L2L
                  IF(MUF.LT.1.OR.MUF.GT.4)   THEN
                     WRITE(EMESS,'(A,I2,A)')                            &
     &                  'MUF =',MUF,' NOT IN RANGE 1 TO 4'
                     CALL ERROR_MESSAGE(NSEQP1)
                  END IF
                  CALL TESTDF(1,N2L)
!*****************CHECK FOR CORRECT L-VALUE FOR THIS SUBSECTION
                  CALL TEST3(L1L,NL-1,'L')
!*****************CHECK J VALUE
                  AJ = Y(2)
                  CALL TEST6(AJ,AJLO,AJHI,'AJ')
               END DO
            END DO
         END IF
      END IF
!
      RETURN
      END SUBROUTINE CHKUR
!
!***********************************************************************
!
      SUBROUTINE TESTW(NBEG,NVALS,NSTEP)
!
!     TEST THAT PARTIAL WIDTHS ADD UP TO TOTAL
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NBEG,NVALS,NSTEP
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: K,KK,K0,K1,K2
      INTEGER(KIND=I4) :: KT
      REAL(KIND=R4) :: TOT,BK,BKK
      REAL(KIND=R4) :: SSUM,DELTA
!
!     INITIALIZE ERROR COUNT
!
      MESS = 0
!
!     PROCESS EACH SET OF PARAMETERS
!
      DO K=NBEG,NVALS,NSTEP
         K0 = K
         K1 = K + 1
         K2 = K + NSTEP - NVALS
         KK = K - 2
         TOT = Y(K0)
!********ERROR IF TOTAL IS ZERO
         IF(TOT.LE.0.0)  THEN
            MESS = 0
            BKK = Y(KK)
            WRITE(EMESS,'(A,1PE12.5)')                                  &
     &        'TOTAL WIDTH LESS THAN OR EQUAL TO ZERO AT ENERGY=',BKK
            CALL ERROR_MESSAGE(0)
            GO TO 50
         END IF
!
!        ADD UP PARTIALS
!
         SSUM = 0.0
         DO KT=K1,K2
            SSUM = SSUM + ABS(Y(KT))
         END DO
!
!        CHECK SUM AGAINST PARTIAL
!
         DELTA = ABS(1.-SSUM/TOT)
         IF(DELTA.GT.EPSILN3)   THEN
            IF(MESS.EQ.0)  THEN
               EMESS = 'SUM OF PARTIALS DOES NOT ADD UP TO TOTAL '//    &
     &                    'AT THE FOLLOWING POINTS'
               CALL ERROR_MESSAGE(0)
            END IF
            MESS = MESS + 1
            BKK = Y(KK)
            BK = Y(K0)
            WRITE(EMESS,'(A,1PE12.5,A,1PE12.5,A,1PE12.5)')              &
     &         'ENERGY=',BKK,' GAMMA-TOTAL=',BK,' SUM=',SSUM
            CALL ERROR_MESSAGE(0)
         END IF
 50   CONTINUE
      END DO
!
      RETURN
      END SUBROUTINE TESTW
!
!***********************************************************************
!
      SUBROUTINE TESTE(NPLT,L,NL,J)
!
!     ROUTINE TO COMPARE UNRESOLVED ENERGY REGION GRIDS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NPLT,L,NL,J
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: I2,NPLFST
      INTEGER(KIND=I4) :: I,I1
      REAL(KIND=R4) :: ENU,DEV
      SAVE NPLFST,I2
!
!     SAVE GRID ON FIRST PASS
!
      IF(J.EQ.1)   THEN
         I2 = 0
         NPLFST = NPLT
         DO I=7,NPLT,L
            I2 = I2 + 1
            EURGRID(I2) = Y(I)
         END DO
         GO TO 100
      END IF
!
!     COMPARE WITH STORED DATA ON SUCCEEDING PASSES
!
      IF(NPLFST.EQ.NPLT)   THEN
         I1 = 1
         DO I=1,I2
            I1 = I1 + L
            ENU = Y(I1)
            DEV = ABS(1.-ENU/EURGRID(I))
            IF(DEV.GT.EPSILN5)   THEN
               WRITE(EMESS,'(A,1PE12.5,A,I2,A,I2)')                     &
     &             'ENERGY POINT',ENU,' L STATE',NL,' J STATE',J
               CALL ERROR_MESSAGE(0)
               WRITE(EMESS,'(4X,A,1PE12.5)')                            &
     &             'DIFFERS FROM VALUE FOR FIRST L AND J STATE',        &
     &              EURGRID(I)
               CALL ERROR_MESSAGE(0)
            END IF
         END DO
      ELSE
!
!        NUMBER OF POINTS DIFFER
!
         WRITE(EMESS,'(A,I2,A,I2)')                                     &
     &           'ENERGY POINTS FOR L STATE',NL,' J STATE',J
         CALL ERROR_MESSAGE(0)
         EMESS = '    DOES NOT EQUAL THE NUMBER OF ENERGY POINTS '//    &
     &           'FOR THE FIRST L AND J STATE'
         CALL ERROR_MESSAGE(0)
      END IF
!
!     TEST THAT ENERGIES ARE IN INCREASING ORDER
!
  100 CALL TEST5Y(1,NPLT,6,1)
!
      RETURN
      END SUBROUTINE TESTE
!
!***********************************************************************
!
      SUBROUTINE TESTDF(N,NW)
!
!     TEST FOR LEGAL DEGREES OF FREEDOM
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N,NW
!
      INTEGER(KIND=I4) :: LINENO
      INTEGER(KIND=I4) :: I,J
      REAL(KIND=R4) :: AM
!
      CHARACTER(LEN=4), DIMENSION(4), PARAMETER ::                      &
     &        A = (/'AMUX','AMUN','AMUG','AMUF'/)
      REAL(KIND=R4), DIMENSION(2,4) :: OLIMS
      DATA OLIMS/1.,2.,1.,2.,0.,0.,1.,4./
!
!     SAVE LINE NUMBER
!
      LINENO = NSEQP1 + 1
!
!     TEST ONLY NEUTRON WIDTH DEGREES OF FREEDOM
!
      IF(N.NE.2)   THEN
         AM = Y(3)
         IF(AM.LT.OLIMS(1,2).OR.AM.GT.OLIMS(2,2)) THEN
            WRITE(EMESS,'(2A,F4.1,A)')                                  &
     &           A(2),' = ',AM,' NOT IN SPECIFIED RANGE'
            CALL ERROR_MESSAGE(LINENO)
         END IF
         GO TO 100
      END IF
!
!     TEST FOR ALL WIDTHS
!
      DO I=1,4
!********SEE IF COMPETITIVE OR FISSION WIDTHS ALL ZERO
         IF(I.LE.1.OR.I.GE.4)   THEN
            DO J=1,NW
               IF(Y(6*J+2+I).GT.0.0)   GO TO 70
            END DO
            GO TO 90
         END IF
!
!     TEST FOR LEGAL DEGREES OF FREEDOM
!
   70    AM = Y(I+2)
         IF(AM.LT.OLIMS(1,I).OR.AM.GT.OLIMS(2,I))  THEN
            WRITE(EMESS,'(2A,F4.1,A)')                                  &
     &           A(I),' = ',AM,' NOT IN SPECIFIED RANGE'
            CALL ERROR_MESSAGE(LINENO)
         END IF
 90   CONTINUE
      END DO
!
  100 RETURN
      END SUBROUTINE TESTDF
!
!***********************************************************************
!
      SUBROUTINE TESTSP(SPIN)
!
!     ROUTINE TO TEST LIMITS ON SPIN
!       ALSO INSURE INTEGRAL OR HALF INTEGRAL
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: SPIN
!
      REAL(KIND=R4), INTRINSIC :: ABS, FLOAT
!
      INTEGER(KIND=I4) :: ISPI
      REAL(KIND=R4) :: DIF
!
!     TEST SPIN
!
      IF(SPIN.LT.0.) THEN
!********TEST FOR NEGATIVE SPIN
         EMESS = 'NEGATIVE SPIN NOT ALLOWED'
         CALL ERROR_MESSAGE(NSEQP1)
      ELSE
!********TEST SPIN LIMITS
         CALL TEST6(SPIN,0.0,22.0,'SPI')
      END IF
!
!     TEST SPIN TO SEE IF INTEGRAL OR HALF-INTEGRAL
!
      ISPI = NINT(SPIN)
      DIF = ABS(SPIN-FLOAT(ISPI))
      IF(DIF.NE.0.0.AND.DIF.NE.0.5)   THEN
         EMESS = 'SPIN SHOULD BE INTEGRAL OR HALF INTEGRAL'
         CALL ERROR_MESSAGE(NSEQP1)
      END IF
!
      RETURN
      END SUBROUTINE TESTSP
!
!***********************************************************************
!
      SUBROUTINE CKF3
!
!     CHECK FILE 3 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: JPART
      INTEGER(KIND=I4) :: LR
      INTEGER(KIND=I4) :: MTT,MTL
      INTEGER(KIND=I4) :: NBEG,NLMOD,NCONT,IPART
      INTEGER(KIND=I4) :: N,NNN
      REAL(KIND=R4) :: Q,QM,QT
      REAL(KIND=R4) :: ELO,EHI
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
      IFL3 = 1
!
!     INITIALIZE FOR SUMUP TEST FIRST TIME
!
      IF(ITEST.EQ.0.AND.FIZCON_DATA%ISUM.GT.0)   CALL SUMF3(-1)
!
!     READ DATA TABLE
!
      CALL RDTAB1
!
!     SET A FLAG IF ALL VALUES OF CHARGED PARTICLE ELASTIC SIGMA
!       ARE SET TO 1.0
!
      JPART = NSUB/10
      IF(MT.EQ.2.AND.JPART.NE.1)   THEN
         DO N=1,NP
            IF(Y(N).NE.1.0)   THEN
               CPELAS = 0
               GO TO 10
            END IF
         END DO
         CPELAS = 1
      END IF
!
!     TEST LR
!
   10 LR = L2
      Q = C2
      QM = C1
      CALL TESTLR(LR,QM)
!
!     DO Q VALUE TESTS
!
      CALL TESTQ(QM,Q,LR,X(NP))
!
!     CHECK ENERGY SPAN OF SECTION
!
      ELO = X(1)
      EHI = X(NP)
      IF((MT.GE.18.AND.MT.LE.21).OR.MT.EQ.38)   THEN
         QT = QUNK
      ELSE
         QT = Q
      END IF
      CALL TESTER(ELO,EHI,QT)
!
!     SAVE ENERGY SPAN OF SECTION
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     TEST FOR MISSING LEVELS
!
      IF(MT.GE.50)   THEN
         IF(MT.LE.91)    THEN
            NBEG = 50
            NLMOD = 50
            NCONT = 41
            IPART = 1
         END IF
         IF(NFOR.GE.6) THEN
            IF(MT.LT.600.OR.MT.GT.849)    GO TO 20
            NBEG = 600
            NLMOD = 50
            NCONT = 49
         ELSE
            IF(MT.LT.699.OR.MT.GT.799)    GO TO 20
            NBEG = 700
            NLMOD = 20
            NCONT = 18
         END IF
         MTT = MT - NBEG
         MTL = MOD(MTT,NLMOD)
         IF(NBEG.NE.50)   IPART = (MTT/NLMOD) + 3
         IF(MTL.GE.1.AND.MTL.LT.NCONT)   THEN
            JPART = NSUB/10
            IF(MTL.NE.1.OR.JPART.NE.IPARTS(IPART))  THEN
               CALL TESTP(MF,MT-1)
            END IF
         END IF
      END IF
!
!     SAVE SECTION IF NEEDED FOR FILE 9 AND 10 TESTS
!
   20 IF(NISSEC.NE.0)   THEN
!
         IF(NISSEC.GT.SZLMF) STOP 'FIZCON ERROR - SZLMF limit exceeded'
!
         DO NNN=1,NISSEC
            IF(MTISO(NNN).EQ.MT)  THEN
              CALL RDWRIT(ISCRU2,2)
              GO TO 25
            END IF
         END DO
      END IF
!
!     IF SUMUP DESIRED, DO IT
!
   25 IF(FIZCON_DATA%ISUM.NE.0)   CALL SUMF3(MT)
!
      RETURN
      END SUBROUTINE CKF3
!
!***********************************************************************
!
      SUBROUTINE TESTLR(LR,S)
!
!     SUBROUTINE TESTS FOR A VALID LR FLAG
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LR
      REAL(KIND=R4) :: S
!
!     LR = 0 ALWAYS ALLOWED
!
      IF(LR.EQ.0)   GO TO 100
!
!     LR = 1  VERSION 6 FORMAT ONLY
!
      IF(LR.EQ.1)   THEN
         IF(NFOR.GE.6)    GO TO 100
      ELSE
!
!        LR GT 1  NEUTRON INCIDENT FILES ONLY
!
         IF(NSUB/10.NE.1)    GO TO 90
!
!        VALID ONLY FOR DISCRETE LEVELS
!
         IF(MT.GE.600.AND.MT.LE.849)   THEN
            IF(MOD(MT,50).NE.49)     GO TO 50
            GO TO 90
         END IF
         IF(MT.LT.50.OR.MT.GT.91)   GO TO 90
!
!        CHECK FOR VALID LR VALUE
!
   50    IF(LR.EQ.16.OR.LR.EQ.17)  THEN
            IF(NFOR.LT.6) GO TO 100
         ELSE IF(LR.GE.22.AND.LR.LE.25)  THEN
            GO TO 100
         ELSE IF(LR.GE.28.AND.LR.LE.36)  THEN
            GO TO 100
         ELSE IF(LR.EQ.39.OR.LR.EQ.40)  THEN
            IF(NFOR.LT.6)   CALL TEST3F(S,0.,'S')
            GO TO 100
         END IF
      END IF
!
!     BAD LR FLAG
!
   90 WRITE(EMESS,'(A,I3,A)')  'LR=',LR,' INVALID'
      CALL ERROR_MESSAGE(NSEQP1)
!
  100 RETURN
      END SUBROUTINE TESTLR
!
!***********************************************************************
!
      SUBROUTINE TESTQ(QM,QI,LR,EHI)
!
!     SUBROUTINE TESTS Q-VALUE TO SEE IF REASONABLE
!     Q MUST BE ASCENDING FOR MTS 50-90, 600-849
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LR
      REAL(KIND=R4) :: QM,QI,EHI
!
      INTEGER(KIND=I4) :: ILEVC,IEQU,IQTEST
      REAL(KIND=R4) :: ELEV,QTLOW,Q,EXL
      REAL(KIND=R4), PARAMETER :: QLOW =-2.0E+07,QHIGH =2.0E+07
      REAL(KIND=R4), PARAMETER :: QFLOW= 1.7E+08,QFHIGH=2.3E+08
      REAL(KIND=R4), PARAMETER :: QLLOW= 1.0E+03,QLHIGH=2.0E+07
!
!     SET UP FOR Q TEST
!
      CALL SETUP_Q(LR,QM,QI,Q,EXL,ILEVC,IEQU,IQTEST)
!
!     Q TESTS
!
      IF(IQTEST.EQ.1) THEN
         IF(Q.EQ.0.)   GO TO 50
      ELSE IF(IQTEST.EQ.2) THEN
         IF(Q.GE.0.)   GO TO 50
      ELSE IF(IQTEST.EQ.3) THEN
         IF(Q.EQ.ELIS)    GO TO 50
      ELSE IF(IQTEST.EQ.4) THEN
         IF(Q.NE.ELIS)   GO TO 50
         QTLOW = AMIN1(-EHI,QLOW)
         IF(Q.GE.QTLOW.AND.Q.LE.QHIGH)   GO TO 50
      ELSE IF(IQTEST.EQ.5) THEN
         QTLOW = AMIN1(-EHI,QLOW)
         IF(Q.GE.QTLOW.AND.Q.LE.QHIGH)   GO TO 50
      ELSE IF(IQTEST.EQ.6) THEN
         IF(Q.LE.0.)   GO TO 50
      ELSE IF(IQTEST.EQ.7) THEN
         IF(Q.GE.QFLOW.AND.Q.LE.QFHIGH)   GO TO 50
      ELSE
!        No explicit test to be performed
         GO TO 100
      END IF
      WRITE(EMESS,'(A,1PE12.5,A)')                                      &
     &      'Q=',Q,' MIGHT BE UNREASONABLE'
      CALL WARNING_MESSAGE(1)
!
!     CHECK IMPLIED INTERMEDIATE LEVEL ENERGY
!
   50 IF(ILEVC.EQ.0)   THEN
         IF(EXL-ELIS.NE.0.)   THEN
            EMESS = 'IMPLIED INTERMEDIATE LEVEL ENERGY SHOULD BE 0.0'
            CALL ERROR_MESSAGE(1)
            GO TO 100
         END IF
      END IF
      ELEV = EXL + ELIS
      IF(IEQU.EQ.0.OR.ELEV.NE.ELIS)  THEN
         IF(ELEV.EQ.0..OR.(ELEV.GE.QLLOW.AND.ELEV.LE.QLHIGH)) GO TO 100
      END IF
      IF(ILEVC.EQ.1.OR.LR.EQ.0)   THEN
         WRITE(EMESS,'(A,1PE12.5,A)')                                   &
     &          'ELEVEL=',ELEV,' MIGHT BE UNREASONABLE'
         CALL WARNING_MESSAGE(1)
      END IF
!
  100 RETURN
      END SUBROUTINE TESTQ
!
!***********************************************************************
!
      SUBROUTINE SETUP_Q(LR,QM,QI,Q,EXL,ILEVC,IEQU,IQTEST)
!
!     ROUTINE TO SETUP Q TEST
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LR,ILEVC,IEQU,IQTEST
      REAL(KIND=R4) :: QM,QI,Q,EXL
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: IPART,JPART
      INTEGER(KIND=I4) :: NBEG,NLMOD,NCONT,MTL,MTT
      REAL(KIND=R4) :: DELTA
!
      INTEGER(KIND=I4), PARAMETER :: NPARS=6
      INTEGER(KIND=I4), DIMENSION(NPARS) , PARAMETER ::                  &
     &          ISAME = (/0,1001,1002,1003,2003,2004/)
!
!     SETUP FOR ENDF-5 OR -6 FORMAT
!
      IQTEST=0
      IF(NFOR.GE.6) THEN
         Q = QM
         EXL = QM - QI + ELIS
      ELSE
         Q = QI
         EXL = ELIS
      END IF
!
!     SAVE Q-VALUE FOR EACH MT IN FILE 3
!
      NMT3 = NMT3 + 1
      MT3(NMT3) = MT
      QMVAL(NMT3) = Q
      QVAL(NMT3) = QI
!
!     SET UP TO CHECK Q VALUES AND LEVEL ORDER
!
      IPART = NSUB/10
      JPART = -1
      ILEVC = 0
      IEQU = 0
!
!     MT = 1 - 49
!
      IF(MT.LE.49)   THEN
!        IF(MT.EQ.1.OR.MT.EQ.2.OR.MT.EQ.5) THEN
         IF(MT.EQ.1.OR.MT.EQ.2) THEN
            IQTEST = 1
         ELSEIF(MT.EQ.3) THEN
            IQTEST = 2
         ELSE IF(MT.EQ.4) THEN
            IF(NFOR.GE.6)   THEN
               ILEVC = 1
               IF(IPART.EQ.1) THEN
                  IQTEST = 3
               ELSE
                  IQTEST = 4
               END IF
            ELSE
               IQTEST = 5
            END IF
         ELSE IF(MT.GE.6.AND.MT.LE.10) THEN
            IQTEST = 4
         ELSE IF((MT.GE.11.AND.MT.LE.17).OR.MT.EQ.37)   THEN
            IF(IPART.EQ.1)  THEN
               IQTEST = 6
            ELSE
               IQTEST = 4
            END IF
         ELSE IF((MT.GE.18.AND.MT.LE.21).OR.MT.EQ.38)   THEN
!***********CHECK IF FISSION Q MATCHES ENERGY RELEASE IN MT458
            IF(MT458.NE.0)  THEN
               IF(Q.NE.0.0)  THEN
                  DELTA = ABS(ERQ-Q)/Q
               ELSE
                  DELTA = ABS(ERQ)
               END IF
               IF(DELTA.GT.EPSILN3)   THEN
                  EMESS = 'Q VALUE NOT COMPATIBLE WITH MF=1, MT=458'
                  CALL ERROR_MESSAGE(NSEQP1)
                  WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')               &
     &                'Q=',Q,' ENERGY RELEASE=',ERQ
                  CALL ERROR_MESSAGE(0)
               END IF
            END IF
            IQTEST = 7
         ELSE
            IQTEST = 4
         END IF
!
!     MT = 50 - 91, SINGLE OUTGOING NEUTRONS
!
      ELSE IF(MT.LE.91)   THEN
         NBEG = 50
         NLMOD = 50
         NCONT = 41
         ILEVC = 2
         MTL = MT - NBEG
         JPART = 1
         IF(IPART.EQ.JPART) THEN
            IEQU = 1
         ELSE
            IEQU = 0
         END IF
         CALL CHK_LEVEL(Q,QI,QM,LR,MTL,NCONT,IEQU,QMVAL(NMT3))
         IF((LR.EQ.0.OR.LR.EQ.31).OR.(LR.EQ.39.OR.LR.EQ.40)) THEN
            IF(IEQU.EQ.1)   THEN
               IQTEST = 3
            ELSE
               IQTEST = 4
            END IF
         ELSE
            IQTEST = 5
         END IF
!
!     NO MT'S BETWEEN 92 AND 100
!
      ELSE IF(MT.LE.100) THEN
         GO TO 100
!
!     MTS FROM 101 TO 207
!
      ELSE IF(MT.LE.207)   THEN
         IF(NFOR.GE.6)   ILEVC = 1
!********MT = 101 - 107
         IF(MT.EQ.101) THEN
            IQTEST = 2
         ELSE IF(MT.GE.102.AND.MT.LE.107) THEN
            IF(IPART.EQ.ISAME(MT-101))   THEN
               IQTEST = 3
            ELSE
               IQTEST = 4
            END IF
!********MT = 108 - 120
         ELSE IF(MT.GE.108.AND.MT.LT.120) THEN
            IQTEST = 4
         ELSE IF(MT.EQ.120) THEN
            IQTEST = 2
         ELSE IF(MT.GE.121.AND.MT.LT.201) THEN
            GO TO 100
!********MT = 201 - 207
         ELSE IF(MT.EQ.201) THEN
            IF(IPART.EQ.1)   THEN
               IQTEST = 3
            ELSE IF(IPART.EQ.ISAME(MT-201)) THEN
               IQTEST = 3
            END IF
         ELSE
            IQTEST = 4
         END IF
!
!     MT > 600, SINGLE OUTGOING CHARGED PARTICLES
!
      ELSE IF(MT.GE.600.AND.MT.LE.849) THEN
         IF(NFOR.GE.6)   THEN
            NBEG = 600
            NLMOD = 50
            NCONT = 49
            ILEVC = 1
         ELSE
            IF(MT.LT.700.OR.MT.GT.799)    GO TO 100
            NBEG = 700
            NLMOD = 20
            NCONT = 18
            ILEVC = 2
         END IF
         MTT = MT - NBEG
         MTL = MOD(MTT,NLMOD)
         JPART = ISAME((MTT/NLMOD)+2)
         IF(IPART.EQ.JPART) THEN
            IEQU = 1
         ELSE
            IEQU = 0
         END IF
         CALL CHK_LEVEL(Q,QI,QM,LR,MTL,NCONT,IEQU,QMVAL(NMT3))
         IF((LR.EQ.0.OR.LR.EQ.31).OR.(LR.EQ.39.OR.LR.EQ.40)) THEN
            IF(IEQU.EQ.1) THEN
               IQTEST = 3
            ELSE
               IQTEST = 4
            END IF
         ELSE
            IQTEST = 5
         END IF
      ELSE
         GO TO 100
      END IF
!
  100 RETURN
      END SUBROUTINE SETUP_Q
!
!***********************************************************************
!
      SUBROUTINE CHK_LEVEL(Q,QI,QM,LR,MTL,NCONT,IEQU,QMVALT)
!
!     CHECK Q VALUE AND LEVEL ORDER FOR SINGLE PARTICLE EMISSION
!       CHANNELS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LR,MTL,NCONT,IEQU
      REAL(KIND=R4) :: Q,QI,QM,QMVALT
!
      INTEGER(KIND=I4), SAVE :: IFLEV,ISETQM
      INTEGER(KIND=I4), SAVE :: LRPR
      REAL(KIND=R4), SAVE :: QMSAV,EXL
      REAL(KIND=R4), SAVE :: EXLP
!
!     INITIALIZE ON FIRST MT OF AN OUTGOING PARTICLE TYPE
!
      IF(MTL.EQ.0) THEN
         QMSAV = Q
         IFLEV = 1
         ISETQM = 1
      ELSE IF(MTL.EQ.1)  THEN
         IF(IEQU.EQ.1)   THEN
            QMSAV = ELIS
            IFLEV = 1
            ISETQM = 1
         ELSE
            IFLEV = 0
            ISETQM = 0
         END IF
         LRPR = 0
         EXLP = 0.0
      ELSE IF(MTL.GT.1) THEN
         IFLEV = 0
         ISETQM = 0
      END IF
!
!     IN ENDF-6 FORMAT CHECK CONSISTANCY OF QI AND QM
!
      IF(NFOR.GE.6)  THEN
         EXL = QMSAV - QI + ELIS
         IF((LR.EQ.0.OR.LR.EQ.31).OR.(LR.EQ.39.OR.LR.EQ.40)) THEN
            IF(MTL.EQ.0) THEN
               CALL TEST3F(QI,QM,'QI')
            ELSE
               IF(ISETQM.EQ.1) CALL TEST3F(QM,QMSAV,'QM')
            END IF
         END IF
      ELSE
         QMVALT = QMSAV
         Q = QMSAV
         EXL = Q - QI + ELIS
      END IF
!
!     CHECK ORDER OF LEVELS
!
      IF(MTL.LT.NCONT)   THEN
         IF(IFLEV.EQ.0)  THEN
            IF(EXLP.GT.EXL) THEN
               GO TO 90
            ELSE IF(EXLP.EQ.EXL) THEN
!**************LEVEL ENERGIES EQUAL OK ONLY IF LR FLAGS DIFFER
               IF(MF.NE.3.OR.LRPR.EQ.LR)   GO TO 90
            END IF
         END IF
!********LEVEL ENERGY ORDER OK
         EXLP = EXL
         LRPR = LR
      END IF
      GO TO 100
!
!     ERROR MESSAGE
!
   90 WRITE(EMESS,'(A,I4)')                                             &
     &   'SECTIONS ARE NOT IN INCREASING LEVEL ENERGY ORDER AT MT =',MT
      CALL ERROR_MESSAGE(0)
      NERROR = NERROR + 1
!
  100 RETURN
      END SUBROUTINE CHK_LEVEL
!
!***********************************************************************
!
      SUBROUTINE CKF4
!
!     CHECK FILE 4 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: LTT,LVT,LI,LCT
      INTEGER(KIND=I4) :: NM
      INTEGER(KIND=I4) :: ICONT,MTT,MTL
      INTEGER(KIND=I4) :: NBEG,NLMOD,NCONT
      INTEGER(KIND=I4) :: NE
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: ELO,EHI,FNORM
      REAL(KIND=R4), DIMENSION(2) :: X2
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
      LTT = L2H
!
!     NO TRANSFORMATION MATRIX
!
      LVT = L1H
      IF(LVT.EQ.0)   THEN
         CALL RDCONT
         LI = L1H
         LCT = L2H
!
!     WITH TRANSFORMATION MATRIX
!
      ELSE
         CALL RDLIST
         LI = L1L
         LCT = L2L
         NM = N2L
         IF(NM.GE.2) THEN
            IF(MOD(NM,2).NE.0)   THEN
               WRITE(EMESS,'(A,I3,A)')  'NM=',NM,' SHOULD BE EVEN'
               CALL WARNING_MESSAGE(NSEQP1)
            END IF
         END IF
      END IF
!
!     DETERMINE IF A CONTINUUM OR DISCRETE CHANNEL REACTION
!
      IF((MT.GE.50.AND.MT.LT.91).OR.MT.EQ.2)  THEN
         ICONT = 0
      ELSE
         ICONT = 1
      END IF
      IF(NFOR.GE.6)   THEN
         IF(MT.LT.600.OR.MT.GT.849)    GO TO 30
         NBEG = 600
         NLMOD = 50
         NCONT = 49
      ELSE
         IF(MT.LT.699.OR.MT.GT.799)    GO TO 30
         NBEG = 700
         NLMOD = 20
         NCONT = 18
      END IF
      MTT = MT - NBEG
      MTL = MOD(MTT,NLMOD)
      IF(MTL.LT.NCONT)   ICONT = 0
!
!     CHECK IF FRAME OF REFERENCE APPROPRIATE TO CHANNEL TYPE
!
   30 IF(ICONT.EQ.1.AND.LCT.EQ.2)  THEN
         EMESS = 'CONTINUUM REACTION RECOMMENDS LCT=1'
         CALL WARNING_MESSAGE(NSEQP1)
      ELSE IF(ICONT.EQ.0.AND.LCT.EQ.1)  THEN
         EMESS = 'DISCRETE CHANNEL REACTION REQUIRES LCT=2'
         CALL ERROR_MESSAGE(NSEQP1)
      END IF
!
!     ISOTROPIC SO ONLY MAKE SURE FILE 3 EXISTS
!
      IF(LI.NE.0)   THEN
         CALL TESTP(3,MT)
         GO TO 100
      END IF
!
!     LEGENDRE EXPANSIONS
!
      X2(1) = -BIGNO
      IF(LTT.EQ.1.OR.LTT.EQ.3)  THEN
         CALL RDTAB2
         NE = NP2
         DO N=1,NE
            CALL RDLIST
!
!           NL is determined by the behavior of the Blatt-Biedenharn Zbar coefficients.
!           Inside each Zbar coefficient, there is a Racah coefficient
!           and a Clebsh-Gordon coefficient.  The CG coefficient looks like this:
!               ( l1 l2  L )
!               (  0  0  0 )
!           So, this means two things:
!               1. The CG coeff (and hence Zbar) will be zero if l1+l2+L=odd.
!               2. The maximum value of L will be l1max+l2max.  Hence, NL=2*lmax.
!
            IF(NPL.GE.2.AND.MOD(NPL,2).NE.0) THEN
               WRITE(EMESS,'(A,I3,A)')                                  &
     &            'NL=',NPL,' SHOULD BE EVEN'
               CALL WARNING_MESSAGE(NSEQP1)
            END IF
            NP = NPL
            CALL TEST6Y (-1.0,1.0,'FL')
!
!           SAVE MIN AND MAX INCIDENT ENERGY
!
            IF(N.EQ.1) THEN
               ELO = C2L
            ELSE IF(N.EQ.NE) THEN
               EHI = C2L
            END IF
!
!           CHECK ENERGIES ARE IN INCREASING ORDER
!
            X2(2) = C2L
            CALL TEST5(X2,2,1)
            X2(1) = X2(2)
         END DO
      END IF
!
!     TABULAR EXPANSIONS
!
      IF(LTT.EQ.2.OR.LTT.EQ.3)  THEN
         CALL RDTAB2
         NE = NP2
         DO N=1,NE
            CALL RDTAB1
            CALL TEST6X (-1.0,1.0,'MU')
            CALL TEST7(FNORM,1)
!
!           SAVE MIN AND MAX INCIDENT ENERGY
!
            IF(LTT.EQ.2.AND.N.EQ.1)   ELO = C2
            IF(N.EQ.NE)   EHI = C2
!
!           CHECK ENERGIES ARE IN INCREASING ORDER
!
            IF(LTT.EQ.3.AND.N.EQ.1) THEN
               CALL TEST3F(C2,EHI,'E1')
            ELSE
               X2(2) = C2
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
            END IF
         END DO
      END IF
!
!     SAVE ENERGY SPAN OF SECTION
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     CHECK THAT RANGE SPANNED IS THE SAME AS FILE 3
!
      CALL ISFIL(MF,3,MT,MT)
!
  100 RETURN
      END SUBROUTINE CKF4
!
!***********************************************************************
!
      SUBROUTINE CKF5
!
!     CHECK FILE 5 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: IFISFL,IFMT
      INTEGER(KIND=I4) :: LF
      INTEGER(KIND=I4) :: NK,NE
      INTEGER(KIND=I4) :: NSEQH
      INTEGER(KIND=I4) :: N,NM
      REAL(KIND=R4) :: ELO,ELOS,EHI,EHIS,U
      REAL(KIND=R4) :: FNORM,EONE,ENE
      REAL(KIND=R4), DIMENSION(2) :: X2
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     INITIALIZE
!
      IF((MT.GE.18.AND.MT.LE.21).OR.(MT.EQ.38)) THEN
         IFISFL = 1
      ELSE
         IFISFL = 0
      END IF
      IF((MT.GE.18.AND.MT.LE.21).OR.(MT.EQ.38.OR.MT.EQ.455))  THEN
         IFMT = 1
      ELSE
         IFMT = 0
      END IF
      IF (NCKF5.EQ.0) THEN
         IMTFIS = 0
         IKTFIS = 0
         ILTFIS = 0
         NCKF5 = 1
      END IF
      ELO = BIGNO
      EHI = 0.0
!
!     STORE # SUBSECTIONS FOR TOTAL AND PARTIAL FISSION CROSS SECTIONS
!
      NK = N1H
      IF(NLIB.EQ.2) THEN
         IF(MT.EQ.18) THEN
            IMTFIS = 1
            IKTFIS = NK
         ELSE IF((MT.GE.19.AND.MT.LE.21).OR.(MT.EQ.38)) THEN
            IF (NK.NE.IKTFIS) THEN
               WRITE(EMESS,'(A,I4,A,I4)')                               &
     &            'The number of subsections in MT',MT,' equals',NK
               CALL ERROR_MESSAGE(0)
               WRITE(EMESS,'(4X,A,I4,2A)')                              &
     &              'MUST be ',IKTFIS,', the number of subsections for',&
     &              ' total fission cross section.'
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
            END IF
         END IF
      END IF
!
!     PROCESS EACH PARTIAL SECONDARY DISTRIBUTION
!
      DO N=1,NK
         CALL RDTAB1
         U = C1
         IF(-U.GE.2.0E+07.AND.-U.LE.3.0E+07)    IFISFL = 1
         NSEQH = NSEQP + 1
         ELOS = X(1)
         IF(ELO.GT.ELOS)   ELO = ELOS
         EHIS = X(NP)
         IF(EHI.LT.EHIS)   EHI = EHIS
         CALL TEST6Y (0.0,1.0,'PKE')
!
!        CHECK BASED ON REPRESENTATION
!
         LF = L2
!
!        STORE LAWS FOR TOTAL FISSION CROSS SECTION AND COMPARE
!        LAWS OF PARTIAL FISSION CROSS SECTIONS WITH THE ONES
!        FOR THE TOTAL FISSION CROSS SECTION
!
         IF(NLIB.EQ.2) THEN
            IF(MT.EQ.18) THEN
               ILTFIS(N) = LF
            ELSE IF((MT.GE.19.AND.MT.LE.21).OR.(MT.EQ.38)) THEN
               IF(ILTFIS(N).EQ.0) THEN
                  WRITE(EMESS,'(A,I2,2A,I3,A)')                         &
     &               'Cannot check consistency of law ',LF,' for ',     &
     &                 'subsection ',N,' with  MT = 18'
                  CALL ERROR_MESSAGE(0)
                  EMESS = '    since MT=18 does not have a '//          &
     &                    'corresponding subsection'
                  CALL ERROR_MESSAGE(0)
               ELSE IF(LF.NE.ILTFIS(N)) THEN
                  WRITE(EMESS,'(A,I2,A,I3,A,I3)')                       &
     &               'Law ',LF,' for subsection ',N,' NOT equal to law',&
     &                  ILTFIS(N)
                  CALL ERROR_MESSAGE(0)
                  EMESS = '     for corresponding subsection in total'//&
     &                    ' fission cross section'
                  CALL ERROR_MESSAGE(0)
               END IF
            END IF
         END IF
!
!        LF=1
!
         IF(LF.EQ.1) THEN
            CALL RDTAB2
            NE = NP2
            X2(1) = -BIGNO
            DO NM=1,NE
               CALL RDTAB1
               CALL TEST7(FNORM,1)
               IF(NM.EQ.1)  THEN
                  EONE = C2
               ELSE IF(NM.EQ.NE)  THEN
                  ENE = C2
               END IF
!**************TEST FOR INCREASING ENERGY ORDER
               X2(2) = C2
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
!**************TEST UPPER LIMIT OF EMITTED PARTICLE
               U = C2
               CALL UTEST(U,1,NP,IFMT)
            END DO
!
!        LF=5
!
         ELSE IF(LF.EQ.5) THEN
            CALL UTEST(U,LF,NP,IFMT)
            CALL RDTAB1
            EONE = X(1)
            ENE = X(NP)
            IF(MT.EQ.455) THEN
               CALL TEST6Y(1.0,1.0,'THT')
            ELSE
               CALL TEST6Y(1.0E+4,1.0E+7,'THT')
            END IF
            CALL RDTAB1
            CALL TEST7(FNORM,1)
!
!        LF=7
!
         ELSE IF(LF.EQ.7) THEN
            CALL UTEST(U,LF,NP,IFMT)
            CALL RDTAB1
            EONE = X(1)
            ENE = X(NP)
            CALL TEST6Y(2.0E+5,5.0E+6,'THT')
!
!        LF=9
!
         ELSE IF(LF.EQ.9) THEN
            CALL UTEST(U,LF,NP,IFMT)
            CALL RDTAB1
            EONE = X(1)
            ENE = X(NP)
            CALL TEST6Y(1.0E+4,1.0E+7,'THT')
!
!        LF=11
!
         ELSE IF(LF.EQ.11) THEN
            CALL UTEST(U,LF,NP,IFMT)
            IFISFL = 1
            CALL RDTAB1
            EONE = X(1)
            ENE = X(NP)
            CALL RDTAB1
!
!        LF=12
!
         ELSE IF(LF.EQ.12) THEN
            CALL UTEST(U,LF,NP,IFMT)
            IFISFL = 1
            CALL RDTAB1
            EONE = X(1)
            ENE = X(NP)
         ELSE
            GO TO 100
         END IF
!
!        CHECK LAW DATA COVERS SAME RANGE AS PROBABILITY
!
         IF(EONE.NE.ELOS.OR.ENE.NE.EHIS)  THEN
            EMESS = 'PARAMETER TABLE ENERGY RANGE INCORRECT'
            CALL ERROR_MESSAGE(NSEQH)
         END IF
      END DO
!
!     SEE THAT A FISSION REACTION HAS A FISSION SPECTRUM
!
      IF(IFMT.EQ.1.AND.IFISFL.EQ.0)   THEN
         EMESS = 'NO FISSION SPECTRUM FOR A FISSION REACTION'
         CALL ERROR_MESSAGE(1)
      END IF
!
!     SAVE ENERGY RANGE SPANNED
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     FOR ALL BUT 455,  ENERGY RANGE SPANNED MUST BE SAME AS FILE 3
!
      IF(MT.NE.455)    CALL ISFIL(MF,3,MT,MT)
!
!     FOR MT=455, CHECK ONLY THAT SAME SECTION EXISTS IN FILE 1
!
      IF(MT.EQ.455)    CALL TESTP(1,MT)
!
  100 RETURN
      END SUBROUTINE CKF5
!
!***********************************************************************
!
      SUBROUTINE UTEST(U,LF,NPT,IFMT)
!
!     ROUTINE TO TEST UPPER LIMIT OF SECONDARY NEUTRON ENERGY
!
      IMPLICIT NONE
!
      INTEGER LF,NPT,IFMT
      REAL(KIND=R4) :: U
!
      INTEGER(KIND=I4) :: NPART
      INTEGER(KIND=I4) :: INMT,NN
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) Q,Q1,QNEG,ESMAX,EAVAIL,ETHRS
!
      INTEGER(KIND=I4), PARAMETER :: NSPARTS=7
      INTEGER(KIND=I4), DIMENSION(2,NSPARTS) :: IPART
      DATA IPART/0,102,1,4,1001,103,1002,104,1003,105,2003,106,2004,107/
!
!     GET Q FROM FILE 3
!
      DO I=1,NMT3
         IF(MT.EQ.MT3(I))  THEN
            Q = QVAL(I)
            GO TO 20
         END IF
      END DO
      GO TO 100
!
!     LF = 1
!
   20 IF(LF.EQ.1)   THEN
         IF(IFMT.EQ.0.AND.MT.NE.91)   THEN
            ESMAX = X(NPT)
            EAVAIL = U + Q
            IF(EAVAIL.LT.ESMAX)  THEN
               WRITE(EMESS,'(A,1PE12.5,A,1PE12.5)')                     &
     &          'FOR LF=1 EPMAX FOUND TO BE',ESMAX,' SHOULD BE',EAVAIL
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
            END IF
         END IF
         GO TO 100
      ELSE
!
!     LF NE 1
!
         IF(IFMT.EQ.1)   GO TO 100
         ETHRS = X(1)
         IF(Q.GE.0.0)   ETHRS = Q
         Q1 = -Q
         IF(MT.NE.91)   THEN
            IF(AWR.LT.40.0)   THEN
               IF(ABS(ABS(U-Q1)/Q1).GT.EPSILN3)   THEN
                  WRITE(EMESS,'(A,I2,A,1PE12.5,A,1PE12.5)')             &
     &               'FOR LF=',LF,' U OF',U,' OUT OF RANGE FOR Q OF ',Q
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
            ELSE
               IF(U.LT.Q1.OR.U.GT.ETHRS)   THEN
                  WRITE(EMESS,'(A,I2,A,1PE12.5,A,1PE12.5)')             &
     &               'FOR LF=',LF,' U OF',U,' OUT OF RANGE FOR Q OF',Q
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
            END IF
            GO TO 100
         ELSE
!
!     INELASTIC CONTINUUM
!
            NPART = NSUB/10
            DO I=1,NSPARTS
               IF(NPART.EQ.IPART(1,I))   GO TO 60
            END DO
            GO TO 100
   60       INMT = IPART(2,I)
            NN = MIN0(INMT,NMT3)
            DO I=1,NN
               IF(MT3(I).EQ.INMT)   GO TO 70
            END DO
            GO TO 100
   70       QNEG = -QVAL(I)
            IF(U.LT.QNEG)   THEN
               WRITE(EMESS,'(A,I2,A,1PE12.5,A,1PE12.5)')                &
     &            'FOR LF=',LF,' U FOUND TO BE',U,' SHOULD BE .GT.',Q
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
            END IF
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE UTEST
!
!***********************************************************************
!
      SUBROUTINE CKF6
!
!     CHECK FILE 6 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD, FLOAT
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: NK,LCT,LF
      INTEGER(KIND=I4) :: JP,JPN,JPP
      INTEGER(KIND=I4) :: NE,INTS
      INTEGER(KIND=I4) :: ND,NEP,NW,NREPT,NDISC,IUPD
      INTEGER(KIND=I4) :: L,LTP,NMU,MM,II,NL
      INTEGER(KIND=I4) :: NSEQH,NSEQC,ICHKER
      INTEGER(KIND=I4) :: I,J,N,NM
      REAL(KIND=R4) :: ELO,ELOS,EHI,EHIS,EONE,ENE,EIN
      REAL(KIND=R4) :: ZAP,ZAPT
      REAL(KIND=R4) :: E,XL,XU,YL,YU,ANS,ANS1,ANSP,XYINT,XYINTI
      REAL(KIND=R4), DIMENSION(2) ::  X2,X3
      REAL(KIND=R4), DIMENSION(201) :: XX,YY
!
      REAL(KIND=R4), PARAMETER :: PERR=5.0*EPSILN4
!
!     INITIALIZE
!
      ELO = BIGNO
      EHI = 0.0
      IF (NCKF6.EQ.0) THEN
         IMTNP = 0
         IKTNP = 0
         IMTNA = 0
         IKTNA = 0
         ILTNP = 0
         ILTNA = 0
         NCKF6 = 1
      END IF
      NK = N1H
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     CHECK FOR PRESENCE OF TOTAL (N,P) AND (N,A) CROSS SECTION
!
      IF(MT.EQ.103) THEN
         IMTNP = 1
         IKTNP = NK
      ELSE IF(MT.EQ.107) THEN
         IMTNA = 1
         IKTNA = NK
      END IF
!
!     STORE # SUBSECTIONS FOR TOTAL AND PARTIAL (N,P) AND (N,A)
!          CROSS SECTIONS
!
      IF(NLIB.EQ.2) THEN
         IF(MT.GE.600.AND.MT.LE.649) THEN
            IF(IMTNP.EQ.1) THEN
               IF (NK.NE.IKTNP) THEN
                  WRITE(EMESS,'(A,I4,A,I3)')                            &
     &               'The number of subsections in MT ',MT,' equals ',NK
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(4X,A,I3)')                              &
     &               'NOT equal to # subsections for total (n,p) '//    &
     &                 'cross section: ',IKTNP
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
            ELSE
               EMESS = 'NO distribution given for total (n,p) cross '// &
     &              'section distribution '
               CALL ERROR_MESSAGE(0)
               WRITE(EMESS,'(4X,A,I3,A,I3)')                            &
     &             'with ',NK,' subsections given for MT ',MT
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
            END IF
         ELSE IF(MT.GE.800.AND.MT.LE.849) THEN
            IF (IMTNA.EQ.1) THEN
               IF (NK.NE.IKTNA) THEN
                  WRITE(EMESS,'(A,I3,A,I3)')                            &
     &               'The number of subsections in MT ',MT,' equals ',NK
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(4X,A,I3)')                              &
     &               'NOT equal to # subsections for total'//           &
     &                 ' (n,alpha) cross section: ',IKTNA
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
            ELSE
               EMESS = 'NO distribution given for total (n,alpha)'//    &
     &              ' cross section distribution '
               CALL ERROR_MESSAGE(0)
               WRITE(EMESS,'(4X,A,I3,A,I3)')                            &
     &             'with ',NK,' subsections given for MT ',MT
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
            END IF
         END IF
      END IF
!
!     SET P(NU) JP, JPP, JPN FLAGS FOR FISSION
!     THESE SHOULD ONLY BE SET FOR MT=18
!
      JP = L1H
      JPP = JP/10
      JPN = JP-JPP
      IF((JP.NE.0).AND.(MT.NE.18)) THEN
         WRITE(EMESS,'(A)') 'JP.GT.0 ONLY ALLOWED FOR MT=18'
      END IF
!
!     LOOP OVER SUBSECTIONS
!
      LCT = L2H
      DO N=1,NK
         CALL RDTAB1
         NSEQH = NSEQP + 1
         ZAP = C1
         IF(MT.EQ.2) THEN
            ZAPT = FLOAT(NSUB/10)
            IF (ZAP.NE.ZA.AND.ZAP.NE.ZAPT) THEN
               CALL TEST3F(ZAP,ZAPT,'ZAP')
            END IF
         END IF
         ELOS = X(1)
         IF(ELO.GT.ELOS) ELO = ELOS
         EHIS = X(NP)
         IF(EHI.LT.EHIS) EHI = EHIS
!
!        STORE LAWS FOR TOTAL (N,P) AND (N,ALPHA) CROSS SECTIONS
!        AND COMPARE LAWS OF PARTIAL CROSS SECTIONS WITH THE ONES
!        FOR THE TOTAL CROSS SECTIONS
!
         LF = L2
         IF(NLIB.EQ.2) THEN
            IF(MT.EQ.103) THEN
               ILTNP(N) = LF
            ELSE IF(MT.EQ.107) THEN
               ILTNA(N) = LF
            END IF
            IF(MT.GE.600.AND.MT.LE.649) THEN
               IF ((IMTNP.EQ.1).AND.(LF.NE.ILTNP(N)))  THEN
                  WRITE(EMESS,'(A,I2,A,I3,A,I3)')                       &
     &               'Law ',LF,' for subsection ',N,' for MT ',MT
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(4X,A,I2,A)')                            &
     &               'NOT equal to law ',ILTNP(N),' for corresponding ' &
     &                 //'subsection in total (n,p) cross section'
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
            ELSE IF(MT.GE.800.AND.MT.LE.849) THEN
               IF ((IMTNA.EQ.1).AND.(LF.NE.ILTNA(N))) THEN
                  WRITE(EMESS,'(A,I2,A,I3,A,I3)')                       &
     &               'Law ',LF,' for subsection ',N,' for MT ',MT
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(4X,A,I2,A)')                            &
     &               'NOT equal to law ',ILTNA(N),' for corresponding ' &
     &                //'subsection in total (n,alpha) cross section'
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
            END IF
         END IF
!
!        TABULAR LAW
!
         IF(LF.EQ.1) THEN
            CALL RDTAB2
            NE = NP2
            INTS = L22
            X2(1) = -BIGNO
            DO I=1,NE
               CALL RDLIST
               E = C2L
               IF(I.EQ.1)  THEN
                  EONE = E
               ELSE IF(I.EQ.NE)  THEN
                  ENE = E
               END IF
!**************TEST FOR INCREASING ENERGY ORDER
               X2(2) = E
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
!**************TEST THAT E-PRIME IS IN INCREASING ORDER
               ND = L1L
               NEP = N2L
               NW = NPL
               NREPT = NW/NEP
               NDISC = ND*NREPT
               IF(ND.NE.0) THEN
                  IF(ZAP.EQ.0.)   THEN
                     IUPD = 0
                  ELSE
                     IUPD = 1
                  END IF
                  CALL TEST5Y(1,NDISC,NREPT,IUPD)
               END IF
               IF(ND.NE.NEP)   THEN
                  CALL TEST5Y(NDISC+1,NW,NREPT,1)
               END IF
!**************TEST NORMALIZATION INTEGRAL
               ANS = 0.0
               IF(ND.NE.0)  THEN
                  DO J=1,ND
                     L = NREPT*(J-1) + 2
                     ANS = ANS + Y(L)
                  END DO
               END IF
               IF(ND+2.LE.NEP) THEN
                  DO J=ND+2,NEP
                     L = NREPT*(J-2) + 1
                     XL = Y(L)
                     XU = Y(L+NREPT)
                     YL = Y(L+1)
                     YU = Y(L+NREPT+1)
                     CALL ECSI(XL,YL,XU,YU,XL,XU,INTS,ANS1)
                     ANS = ANS + ANS1
                  END DO
               END IF
               IF(ABS(ANS-1.0).GT.PERR) THEN
                  WRITE(EMESS,'(A,F11.6,A,1PE11.4)')                    &
     &                  'CHECK NORMALIZATION=',ANS,' AT E=',E
                  CALL ERROR_MESSAGE(NSEQP1)
               END IF
            END DO
            ICHKER = 1
!
!        DISCRETE 2-BODY LAW
!
         ELSE IF(LF.EQ.2) THEN
            IF((MT.GE.50.AND.MT.LE.90).OR.MT.EQ.2)    GO TO 40
            IF(MT.GE.600.AND.MT.LE.849)   THEN
               IF(MOD(MT,50).NE.49)     GO TO 40
            END IF
            WRITE(EMESS,'(A,I4)')                                       &
     &             'DISCRETE 2-BODY LAW NOT PERMITTED FOR MT=',MT
            CALL ERROR_MESSAGE(NSEQP1)
   40       IF(LCT.NE.2)  THEN
               WRITE(EMESS,'(A,I1)')                                    &
     &                   'ONLY LCT=2 ALLOWED FOR LAW ',LF
               CALL ERROR_MESSAGE(NSEQH)
            END IF
            CALL RDTAB2
            NE = NP2
            X2(1) = -BIGNO
            DO I=1,NE
               CALL RDLIST
               E = C2L
               IF(I.EQ.1)  THEN
                  EONE = E
               ELSE IF(I.EQ.NE)  THEN
                  ENE = E
               END IF
!**************TEST FOR INCREASING ENERGY ORDER
               X2(2) = C2L
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
            END DO
            ICHKER = 1
!
!        ISOTROPIC DISCRETE EMISSION
!
         ELSE IF(LF.EQ.3) THEN
            IF(LCT.NE.2) THEN
               WRITE(EMESS,'(A,I1)')                                    &
     &                 'ONLY LCT=2 ALLOWED FOR LAW ',LF
               CALL ERROR_MESSAGE(NSEQH)
            END IF
            ICHKER = 0
!
!        COULOMB ELASTIC LAW
!
         ELSE IF(LF.EQ.5) THEN
            IF(NSUB/10.EQ.1)  THEN
               EMESS ='COULOMB LAW NOT ALLOWED FOR INCIDENT NEUTRONS'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
            IF(MT.NE.2)  THEN
               EMESS = 'COULOMB LAW ONLY ALLOWED FOR MT=2'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
            CALL RDTAB2
            NE = NP2
            X2(1) = -BIGNO
            DO I=1,NE
               CALL RDLIST
               LTP = L1L
               IF(LTP.LE.10.AND.CPELAS.NE.1)  THEN
                  WRITE(EMESS,'(A,I2,A)')                               &
     &               'LTP = ',LTP,' REQUIRES THAT ALL ELASTIC CROSS '// &
     &                   'SECTIONS IN FILE 3 BE SET TO 1.0'
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
               E = C2L
               IF(I.EQ.1)  THEN
                  EONE = E
               ELSE IF(I.EQ.NE) THEN
                  ENE = E
               END IF
!**************TEST FOR INCREASING ENERGY ORDER
               X2(2) = C2L
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
!**************TEST NORMALIZATION INTEGRAL
               IF(LTP.GT.10)  THEN
                  ANS = 0.0
                  NL = N2L
                  INTS = LTP - 10
                  DO J=2,NL
                     L = 2*(J-1) - 1
                     XL = Y(L)
                     XU = Y(L+2)
                     YL = Y(L+1)
                     YU = Y(L+3)
                     CALL ECSI(XL,YL,XU,YU,XL,XU,INTS,ANS1)
                     ANS = ANS + ANS1
                  END DO
                  IF(ABS(ANS-1.0).GT.PERR) THEN
                     WRITE(EMESS,'(A,F11.6,A,1PE11.4)')                 &
     &                 'CHECK NORMALIZATION=',ANS,' AT E=',E
                     CALL ERROR_MESSAGE(NSEQP1)
                  END IF
               END IF
            END DO
            ICHKER = 1
!
!        N-BODY PHASE SPACE
!
         ELSE IF(LF.EQ.6) THEN
            CALL RDCONT
            ICHKER = 0
!
!        ANGLE-ENERGY TABULAR LAW
!
         ELSE IF(LF.EQ.7)  THEN
            IF(LCT.NE.1) THEN
               EMESS = 'ONLY LCT=1 ALLOWED FOR THIS LAW'
               CALL ERROR_MESSAGE(NSEQH)
            END IF
            CALL RDTAB2
            NE = NP2
            X2(1) = -BIGNO
            DO I=1,NE
               CALL RDTAB2
               EIN = C22
               NSEQC = NSEQP
               IF(I.EQ.1)  THEN
                  EONE = EIN
               ELSE IF(I.EQ.NE)  THEN
                  ENE = EIN
               END IF
!**************TEST FOR INCREASING ENERGY ORDER
               X2(2) = C22
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
               NMU = NP2
               X3(1) = -BIGNO
               DO NM=1,NMU
                  CALL RDTAB1
                  XX(NM) = C2
!*****************TEST FOR INCREASING ANGLE COSINE ORDER
                  X3(2) = C2
                  CALL TEST5(X3,2,1)
                  X3(1) = X3(2)
!*****************GET INTEGRAL OVER E-PRIME AND STORE IT
                  CALL TEST7(ANSP,2)
                  YY(NM) = ANSP
               END DO
!***********CHECK THAT INTEGRAL OVER ALL ANGLES IS NORMALIZED TO 1.
               XYINT = 0.
               MM = 1
               DO NM=2,NMU
                  IF(NM.GT.NBT2(MM))   THEN
                     MM = MM + 1
                     IF(MM.GT.NR2)   GO TO 60
                  END IF
                  II = JNT2(MM)
                  CALL ECSI(XX(NM-1),YY(NM-1),XX(NM),YY(NM),XX(NM-1),   &
     &                 XX(NM),II,XYINTI)
                  XYINT = XYINT + XYINTI
               END DO
               IF(ABS(XYINT-1.0).GE.PERR) THEN
                  WRITE(EMESS,'(A,F11.6,A,1PE11.4)')                    &
     &                 'CHECK NORMALIZATION=',XYINT,' AT E=',EIN
                  CALL ERROR_MESSAGE(NSEQC)
               END IF
   60          CONTINUE
            END DO
            ICHKER = 1
         ELSE
            ICHKER = 0
         END IF
!
!        CHECK LAW DATA COVERS SAME RANGE AS PROBABILITY
!
         IF(ICHKER.EQ.1) THEN
            IF(EONE.NE.ELOS.OR.ENE.NE.EHIS)  THEN
               EMESS = 'ENERGY RANGE FOR DISTRIBUTIONS IN LIST RECORDS'
               CALL ERROR_MESSAGE(0)
               EMESS = '    INCONSISTENT WITH TAB1 RECORD'
               CALL ERROR_MESSAGE(NSEQH)
            END IF
         END IF
      END DO
!
!     SAVE ENERGY RANGE SPANNED
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     ENERGY RANGE SPANNED MUST BE SAME AS FILE 3
!
      CALL ISFIL(MF,3,MT,MT)
!
      RETURN
      END SUBROUTINE CKF6
!
!***********************************************************************
!
      SUBROUTINE CKF7
!
!     CHECK FILE 7 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LTHR,LLN,LT
      INTEGER(KIND=I4) :: NS,NB
      INTEGER(KIND=I4) :: K,N,NN,NNN
      REAL(KIND=R4), DIMENSION(2) :: X2
      REAL(KIND=R4), DIMENSION(3) :: BFLAG
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     INCOHERENT INELASTIC SCATTERING
!
      LTHR = L1H
      IF(LTHR.EQ.0)   THEN
         CALL RDLIST
         NS = N2L
         DO NNN=1,NS
            BFLAG(NNN) = Y(6*NNN+1)
         END DO
         LLN = L1L
         IF(Y(1).GT.0)   THEN
            CALL RDTAB2
            NB = NP2
            X2(1) = -BIGNO
            IF(LLN.EQ.1) INEGC = 0
            DO N=1,NB
               CALL RDTAB1
!**************CHECK FOR INCREASING VALUES OF BETA
               X2(2) = C2
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
               LT = L1
               IF(LT.GT.0)   THEN
                  DO K=1,LT
                     CALL RDLIST
                 END DO
               END IF
            END DO
            INEGC = 1
         END IF
!
!        PROCESS EFFECTIVE TEMPERATURE RECORD
!
         IF(NFOR.GE.6) THEN
            CALL RDTAB1
            IF(NS.GT.0)   THEN
               DO NN=1,NS
                  IF(BFLAG(NN).EQ.0.)   CALL RDTAB1
               END DO
            END IF
         END IF
!
!     COHERENT ELASTIC SCATTERING
!
      ELSE IF(LTHR.EQ.1)   THEN
         CALL RDTAB1
         LT = L1
         IF(LT.GE.1) THEN
            DO K=1,LT
               CALL RDLIST
            END DO
         END IF
!
!     INCOHERENT ELASTIC SCATTERING
!
      ELSE
         CALL RDTAB1
      END IF
!
      RETURN
      END SUBROUTINE CKF7
!
!***********************************************************************
!
      SUBROUTINE CKF8
!
!     CHECK FILE 8 DATA
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: AMOD
!
      INTEGER(KIND=I4) :: IA
      REAL(KIND=R4) :: A
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     FISSION PRODUCT YIELDS  (MT=454,459)
!
      IF(MT.EQ.454.OR.MT.EQ.459) THEN
         CALL CHK_FPY
!
!     RADIOACTIVE DECAY DATA   (MT=457)
!
      ELSE IF(MT.EQ.457)   THEN
!********SECTION CANNOT EXIST FOR A NATURAL ELEMENT
         A = AMOD(C1H,1000.)
         IA = A
         IF(IA.EQ.0)  THEN
            WRITE(EMESS,'(A,I3,A)')                                     &
     &         'SECTION',MT,' SHOULD NOT EXIST FOR A NATURAL ELEMENT'
            CALL ERROR_MESSAGE(0)
            NERROR = NERROR + 1
         END IF
         CALL CHK457
!
!     PROCESS RADIOACTIVE PRODUCT DATA
!
      ELSE
         A = AMOD(C1H,1000.)
         IA = A
         CALL CHK_RPD(IA)
      END IF
!
      RETURN
      END SUBROUTINE CKF8
!
!***********************************************************************
!
      SUBROUTINE CHK_FPY
!
!     ROUTINE TO CHECK FISSION PRODUCT YIELDS
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: LEP
      INTEGER(KIND=I4) :: K,N
      REAL(KIND=R4) :: SSUM
      REAL(KIND=R4), DIMENSION(2) :: X2
!
      LEP = L1H
      X2(1) = -BIGNO
!
!     PROCESS EACH YIELD ENERGY
!
      DO N=1,LEP
         CALL RDLIST
!
!        CHECK INTERPOLATION SCHEME
!
         IF(N.NE.1)   THEN
            IF(L1L.LT.1.OR.L1L.GT.5)   THEN
               WRITE(EMESS,'(A,I3,A)')                                  &
     &                'INVALID INTERPOLATION',L1L,' USED'
               CALL ERROR_MESSAGE(NSEQP1)
            END IF
         END IF
!
!        ARE ENERGIES IN INCREASING ORDER?
!
         X2(2) = C1L
         CALL TEST5(X2,2,1)
         X2(1) = X2(2)
!
!        INDEPENDENT YIELDS SHOULD SUM TO 2.
!
         IF(MT.EQ.454) THEN
            SSUM = 0.0
            DO K=3,NPL,4
               SSUM = SSUM + Y(K)
            END DO
            IF(ABS(SSUM-2.0).GT.EPSILN3)   THEN
               WRITE(EMESS,'(A,F9.2,A,1PE15.5,A)')                      &
     &              'FISSION PRODUCT YIELDS SUM TO',SSUM,' AT',C1L,' EV'
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
            END IF
         END IF
      END DO
!
      RETURN
      END SUBROUTINE CHK_FPY
!
!***********************************************************************
!
      SUBROUTINE CHK457
!
!     THIS SUBROUTINE PERFORMS CHECKS ON RADIOACTIVE DECAY DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD, INT, IFIX
      REAL(KIND=R4), INTRINSIC :: ALOG10,FLOAT
!
      INTEGER(KIND=I4) :: IZ,IA
      INTEGER(KIND=I4) :: IC,IS,IRTIC,IS1,ITY
      INTEGER(KIND=I4) :: LCON,LCOV
      INTEGER(KIND=I4) :: NSP,NERX,NSEQ2,NSEQH,NSEQP2,NSEQPP,NEPS
      INTEGER(KIND=I4) :: II,KK,KKK,K1,K2,N
      INTEGER(KIND=I4) :: I,K,III,JJJ
      REAL(KIND=R4) :: F,FT,STYP,SUMX,XDUM,EEE
      REAL(KIND=R4) :: EBREMT,EB,EIC,EBAR,EBREM,EBIND,RECF,ERT,EMAX,    &
     &                 QTEST
      REAL(KIND=R4) :: FD,DFD,FC,DFC,ERX,DERX
      REAL(KIND=R4) :: TY,RIX,DRIX,RISX,DRISX
      REAL(KIND=R4) :: RICC,DRICC,RICK,RICL,RICT,RICCD
      REAL(KIND=R4) :: EPL,EPU
!
      INTEGER(KIND=I4), PARAMETER :: NSPMAX=1000
      REAL(KIND=R4), DIMENSION(NSPMAX) :: E,DE,RI,DRI,EE,XX,DXX,RTYPX
!
      INTEGER(KIND=I4),PARAMETER :: NSPTP=10
      CHARACTER(LEN=9), DIMENSION(NSPTP) , PARAMETER ::                 &
     &     CSTYPE= (/'GAMMA    ','BETA     ','E.C.     ','         ',   &
     &               'ALPHA    ','NEUTRON  ','SF       ','PROTON   ',   &
     &               'ELECTRON ','X-RAY    '/)
      INTEGER(KIND=I4), DIMENSION(NSPTP), PARAMETER ::                   &
     &         LSTYPE = (/6,5,5,0,6,8,3,7,9,6/)
!
      INTEGER(KIND=I4),PARAMETER :: NSPTP1=NSPTP+1
      INTEGER(KIND=I4), DIMENSION(NSPTP1) :: NER
      REAL(KIND=R4), DIMENSION(NSPTP1) ::  ER,DER
!
!     SAVE NUMBER OF SPECTRA AND SEQUENCE NUMBER OF FIRST CARD
!
      NSP = N2H
      NSEQH = NSEQP1
!
!     GET Z AND A
!
      IZ = IFIX(ZA)/1000
      IA = MOD(IFIX(ZA),1000)
!
!     PROCESS DECAY MODES AND AVERAGE ENERGIES
!
      CALL CHKDEC
!
!     INITIALIZE FOR SPECTRUM PROCESSING
!
      IC = 0
      NER = 0
      ER = 0.0
      DER = 0.0
      EBREMT = 0.0
!
!     PROCESS DECAY SPECTRA
!
      IF(NSP.EQ.0)  THEN
         EMESS = 'NO DECAY SPECTRA GIVEN'
         CALL ERROR_MESSAGE(NSEQP)
         GO TO 1000
      END IF
      DO 100 I=1,NSP
      CALL RDLIST
      NERX = N2L
      STYP = C2L
      IS = IFIX(STYP + 1.)
      LCON = L1L
!*****CHECK DISCRETE NORMALIZATION
      FD = Y(1)
      DFD = Y(2)
      CALL TEST6(DFD,0.,FD,'DFD')
      IF(FD.LE.0..AND.LCON.NE.1.AND.NERX.GT.0) THEN
         WRITE(EMESS,'(A,1PE12.5)')                                     &
     &                'DISCRETE NORMALIZATION .LE. 0 FD=',FD
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!*****CHECK CONTINUUM NORMALIZATION
      FC = Y(5)
      DFC = Y(6)
      CALL TEST6(DFC,0.,FC,'DFC')
      IF(FC.LE.0..AND.LCON.NE.0) THEN
         WRITE(EMESS,'(A1,1PE12.5)')                                    &
     &               'CONTINUUM NORMALIZATION .LE. 0 FC=',FC
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!*****CHECK MEAN ENERGY UNCERTAINTY
      ERX = Y(3)
      DERX = Y(4)
      CALL TEST6(DERX,0.,ERX,'DERB')
      ER(IS) = ERX
      DER(IS) = DERX
!
!     NO SPECTRA - DISCRETE OR CONTINUOUS
!
      IF(LCON.EQ.0.AND.NERX.EQ.0) THEN
         IC = IC + 1
         XX(IC) = ERX
         DXX(IC) = DERX
         EE(IC) = ERX
         DE(IC) = DERX
         RI(IC) = 1.0
         DRI(IC) = 0.0
!*****SAVE SOURCE OF RADIATION
         IF(IS.EQ.7) THEN
            RTYPX(IC) = 6.0
         ELSE
            RTYPX(IC) = 10.0
         END IF
         NER(IS) = IC
         GO TO 100
      END IF
!
      IF(LCON.EQ.1) GO TO 80
!
!     DISCRETE SPECTRA
!
!*****CHECKING ERROR IF MORE THAN NSPMAX ENERGIES
      IF(IC+NERX.GT.NSPMAX) THEN
         EMESS = 'TOO MANY DISCRETE SPECTRA FOR CODE TO PROCESS'
         CALL ERROR_MESSAGE(NSEQP1)
      END IF
!*****PROCESS EACH DISCRETE SPECTRUM
      ERT = 0.
      DO K=1,NERX
         IC = MIN0(IC+1,NSPMAX)
         CALL RDLIST
         NSEQP = NSEQP1
         E(IC) = C1L
!********CHECK FOR POSITIVE ENERGY
         IF(C1L.LE.0.)   THEN
            EMESS = 'RADIATION ENERGY MUST BE GREATER THAN 0.0'
            CALL ERROR_MESSAGE(NSEQP1)
!********CHECK INCREASING ENERGY ORDER OF RADIATION
         ELSE IF(ERT.GT.C1L)   THEN
            WRITE(EMESS,'(A,1PE12.5)')                                  &
     &          'RADIATION ENERGY OUT OF ORDER AT',C1L
            CALL ERROR_MESSAGE(NSEQP1)
         ELSE
            ERT = C1L
         END IF
!********CHECK DISCRETE ENERGY UNCERTAINTY
         DE(IC) = C2L
         CALL TEST6(C2L,0.,C1L,'DER')
!********DISCRETE ENERGY MUST BE LESS THAN THE Q OF THE SOURCE MODE
         IRTIC = IFIX(Y(1))
         IS1 = IS - 1
         IF(IS1.GE.0.AND.IS1.LE.7)   THEN
            IF(IS1.EQ.0) THEN
               QTEST = QMAX
            ELSE
               QTEST = AMIN1(QO(IS1),QO(IRTIC))
            END IF
            IF((1.+EPSILN3)*QTEST.LT.E(IC)) THEN
               WRITE(EMESS,'(A,1PE12.5,A,1PE12.5)')                     &
     &              'E(DISCRETE) > Q  E=',E(IC),'  Q=',QTEST
               CALL ERROR_MESSAGE(NSEQP1)
            END IF
         END IF
!********CHECK RADIATION INTENSITY AND UNCERTAINTY
         TY = Y(2)
         RIX = Y(3)
         DRIX = Y(4)
         NSEQP2 = NSEQP1 + 1
         IF(RIX.LE.0.)  THEN
            EMESS ='RADIATION INTENSITY MUST BE GREATER THAN 0.0'
            CALL ERROR_MESSAGE(NSEQP2)
         END IF
         CALL TEST6(DRIX,0.,RIX,'DRI')
!********RIS (PAIR FORMATION COEFFICIENT FOR STYPE=0.0 OR POSITRON
!*****     INTENSITY FOR STYPE=2.0) AND ITS UNCERTAINTY
         RISX = Y(5)
         DRISX = Y(6)
         CALL TEST6(DRISX,0.,RISX,'DRIS')
!********SAVE AVERAGE ENERGY, ABSOLUTE INTENSITY, AND TOTAL ENERGY PER
!*****     DECAY AND UNCERTAINTIES
         EE(IC) = E(IC)
         RI(IC) = FD*RIX
         DRI(IC) = FD*DRIX
         XX(IC) = E(IC)*RI(IC)
         DXX(IC) = DE(IC)*RI(IC) + E(IC)*DRI(IC)
         RTYPX(IC) = Y(1)
         NSEQP2 = NSEQP1 + 1
!********BRANCH ON RADIATION TYPE
!*****
!*****GAMMA RAYS
!*****
         IF(IS.EQ.1) THEN
            IF(NPL.GT.6)  THEN
!*****CHECK ANY INTERNAL CONVERSION COEFFICIENTS
               RICC = Y(7)
               DRICC = Y(8)
               RICK = Y(9)
               RICL = Y(11)
               RICT = RICK + RICL
               RICCD = RICC + DRICC
               NSEQPP = NSEQP2 + 1
               IF(RICCD.LT.RICT)  THEN
                  EMESS = 'RICC MUST BE GREATER THAN OR EQUAL TO '//    &
     &                  'RICK + RICL '
                  CALL ERROR_MESSAGE(NSEQPP)
               END IF
            END IF
!*****CHECK THAT TYPE IS ZERO
            CALL TEST3F(TY,0.0,'TYPE')
!*****RIS SHOULD BE ZERO IF COMPONENTS GIVEN
            IF(NPL.GT.6.AND.RISX.NE.0.) THEN
               EMESS = 'RIS USUALLY IS 0. FOR STYPE  0.0'
               CALL ERROR_MESSAGE(NSEQP2)
            END IF
!*****
!*****BETA MINUS DECAY
!*****
         ELSE IF (IS.EQ.2) THEN
            EB = 0.0
            IF(RISX.NE.0.)   THEN
               WRITE(EMESS,'(A,F4.1)') 'RIS MUST BE 0. FOR STYPE ',STYP
               CALL ERROR_MESSAGE(NSEQP2)
            END IF
            IF(TY.LT.1.0.OR.TY.GT.3.)    THEN
               CALL TEST6(TY,1.,3.,'TYPE')
               GO TO 60
            END IF
!*****CHECK THAT INTENSITY IS POSITIVE
            IF(RI(IC).LE.0.0)   THEN
               EMESS = 'RADIATION INTENSITY LE 0.0'
               CALL ERROR_MESSAGE(NSEQP1)
               GO TO 60
            END IF
!*****CHECK THAT RADIATION ENERGY IS POSITIVE
            IF(E(IC).LE.0.0)   THEN
               EMESS = 'RADIATION ENERGY LE 0.0'
               CALL ERROR_MESSAGE(NSEQP1)
               GO TO 60
            END IF
!***********CALCULATE AVERAGE BETA ENERGY
            ITY = IFIX(TY)
            CALL AVG(E(IC),FLOAT(IZ),FLOAT(IA),ITY-1,EB,F)
            IF (IBAV.EQ.2) THEN
               CALL BETA(E(IC),IZ,IA,1,ITY-1,EB)
            END IF
!********CALCULATE INTERNAL BREMSSTRAHLUNG ENERGY IF REQUIRED
            EBREM = 0.0
            IF (IBREM.EQ.1) THEN
               IF (E(IC).GT.1.0E+6) THEN
                  EBREM = 0.0015*E(IC)*ALOG10(4.0E-6*E(IC)-2.2)
               END IF
            END IF
            EBREMT = EBREMT + EBREM*RIX*FD
            EB = EB - EBREM
            GO TO 50
!*****
!*****ELECTRON CAPTURE OR POSITRON EMISSION
!*****
         ELSE IF(IS.EQ.3) THEN
            EB = 0
!***********FOR ELECTRON CAPTURE, CALCULATE MEAN BINDING ENERGY
!***********AND (IF REQUIRED) INTERNAL BREMSSTRAHLUNG ENERGY
            IF(RIX.GT.RISX) THEN
               ITY = IFIX(TY)
               CALL ECAP(E(IC),ITY-1,IZ,EBIND,EBREM)
               XX(IC) = XX(IC) - EBIND*(RIX-RISX)*FD
               EBREMT = EBREMT + EBREM*(RIX-RISX)*FD
            END IF
            IF(RISX.LE.0.0) GO TO 60
            IF(TY.LT.1.0.OR.TY.GT.3.)    THEN
               CALL TEST6(TY,1.,3.,'TYPE')
               GO TO 60
            END IF
!***********REMOVE PAIR FORMATION MASS FROM TOTAL ENERGY PER DECAY
            XX(IC) = XX(IC) - 2.*EMASS*RISX*FD
!***********CHECK THAT INTENSITY IS POSITIVE
            IF(RI(IC).LE.0.0)   THEN
               EMESS = 'RADIATION INTENSITY LE 0.0'
               CALL ERROR_MESSAGE(NSEQP1)
               GO TO 60
            END IF
!***********CHECK THAT RADIATION ENERGY IS POSITIVE
            IF(E(IC).LE.0.0)   THEN
               EMESS = 'RADIATION ENERGY LE 0.0'
               CALL ERROR_MESSAGE(NSEQP1)
               GO TO 60
            END IF
!***********CALCULATE AVERAGE POSITRON ENERGY
            EIC = E(IC) - 2.*EMASS
            ITY = IFIX(TY)
            IF(EIC.GT.0.0)  THEN
               CALL AVG(EIC,-(FLOAT(IZ)-1.),FLOAT(IA),ITY-1,EB,F)
               IF (IBAV.EQ.2) THEN
                  CALL BETA(E(IC),IZ,IA,2,ITY-1,EB)
               END IF
!**************CALCULATE INTERNAL BREMSSTRAHLUNG ENERGY IF REQUIRED
               EBREM = 0.0
               IF (IBREM.EQ.1) THEN
                  IF (EIC.GT.1.0E+6) THEN
                     EBREM = 0.0015*EIC*ALOG10(4.0E-6*EIC-2.2)
                  END IF
               END IF
               EBREMT = EBREMT + EBREM*RISX*FD
               EB = EB - EBREM
               EB = EB*RISX/RIX
               GO TO 50
            END IF
            GO TO 60
!*****
!*****ANY HEAVY PARTICLE OR DISCRETE ELECTRONS OR X-RAYS
!*****
         ELSE IF(IS.GE.5.AND.IS.LE.10) THEN
!***********ADD RECOIL ENERGY FOR ALPHA DECAY
            IF(IS.EQ.5) THEN
               RECF = AWR/(AWR-3.9682)
               E(IC) = E(IC)*RECF
               XX(IC) = XX(IC)*RECF
               DXX(IC) = DXX(IC)*RECF
               EE(IC) = EE(IC)*RECF
               DE(IC) = DE(IC)*RECF
            END IF
!***********ADD RECOIL ENERGY FOR NEUTRON EMISSION (BUT NOT FISSION)
            IF(IS.EQ.6.AND.IRTIC.NE.6) THEN
               RECF = AWR/(AWR-1.0)
               E(IC) = E(IC)*RECF
               XX(IC) = XX(IC)*RECF
               DXX(IC) = DXX(IC)*RECF
               EE(IC) = EE(IC)*RECF
               DE(IC) = DE(IC)*RECF
            END IF
            IF(IS.NE.9.AND.IS.NE.10) THEN
               IF(NER(1).EQ.0) THEN
                  IF(E(IC).LT.QO(IRTIC)) THEN
                     WRITE(EMESS,'(A,I5)')                              &
     &                   'GAMMA RAY NEEDED, SOURCE MODE= ',IRTIC
                     CALL ERROR_MESSAGE(NSEQP1)
                  END IF
               END IF
            END IF
            CALL TEST3F(TY,0.0,'TYPE')
            IF(RISX.NE.0.)   THEN
               WRITE(EMESS,'(A,F4.1)') 'RIS MUST BE 0. FOR STYPE ',STYP
               CALL ERROR_MESSAGE(NSEQP2)
            END IF
         END IF
         GO TO 75
!*****CHECK LOGFT
   50    IF(ITY.GT.1) THEN
            IF(IS.EQ.2) THEN
               CALL AVG(E(IC),FLOAT(IZ)-1.,FLOAT(IA),0,EEE,F)
            ELSE
               CALL AVG(E(IC),-(FLOAT(IZ)-1.),FLOAT(IA),0,EEE,F)
            END IF
            FT = F*T12/RI(IC)
            IF(FT.LT.1.0E+6) THEN
               EMESS = 'FT VALUE TOO SMALL'
               CALL ERROR_MESSAGE(NSEQP1)
               WRITE(EMESS,'(2X,A,1PE12.5,A,1PE12.5,A,I5)')             &
     &              'FT=',FT,'  E=',E(IC),'  I=',IC
               CALL ERROR_MESSAGE(NSEQP1)
            END IF
         ELSE
            FT = F*T12/RI(IC)
            IF(FT.LT.240.) THEN
               EMESS = 'FT VALUE TOO SMALL'
               CALL ERROR_MESSAGE(NSEQP1)
               WRITE(EMESS,'(2X,A,1PE12.5,A,1PE12.5,A,I5)')             &
     &                 'FT=',FT,'  E=',E(IC),'  I=',IC
               CALL ERROR_MESSAGE(NSEQP1)
            END IF
         END IF
!*****SET THE AVERAGE ENERGY AND ITS UNCERTAINTY
   60    EE(IC) = EB
         IF(EB.GT.0.0) THEN
            DE(IC) = EB*DE(IC)/E(IC)
         ELSE
            DE(IC) = 0.0
         END IF
   75    CONTINUE
      END DO
   80 NER(IS) = IC
!
!     CONTINUUM SPECTRUM
!
      IF(LCON.EQ.0) GO TO 100
      IC = IC + 1
      CALL RDTAB1
      CALL TEST7(XDUM,1)
!*****CHECK SPECTRUM MAX AGAINST APPROPRIATE Q VALUE
      EMAX = X(NP)
      IRTIC = IFIX(C1)
      IS1 = IS - 1
      IF(IS1.GE.0.AND.IS1.LE.7)   THEN
         IF(IS1.EQ.0) THEN
            QTEST = QMAX
         ELSE
            QTEST = AMIN1(QO(IS1),QO(IRTIC))
         END IF
         IF((1.+EPSILN3)*QTEST.LT.EMAX)  THEN
            WRITE(EMESS,'(A,1PE12.5,A,1PE12.5)')                        &
     &           'E(MAXIMUM) > Q  E=',EMAX,'  Q=',QTEST
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
      END IF
!*****CALCULATE SPECTRUM AVERAGE AND NORMALIZATION
      CALL EAVE(EBAR,SUMX)
!*****SAVE AVERAGE ENERGY, ABSOLUTE INTENSITY, AND TOTAL ENERGY PER
!*****  DECAY AND UNCERTAINTIES
      E(IC) = EBAR
      EE(IC) = EBAR
      DE(IC) = 0.0
      RI(IC) = FC*SUMX
      DRI(IC) = DFC*SUMX
      XX(IC) = E(IC)*RI(IC)
      DXX(IC) = E(IC)*DRI(IC)
      RTYPX(IC) = C1
!*****ADD RECOIL ENERGY FOR ALPHA DECAY
      IF(IS.EQ.5) THEN
         RECF = AWR/(AWR-3.9682)
         XX(IC) = XX(IC)*RECF
         DXX(IC) = DXX(IC)*RECF
         EE(IC) = EE(IC)*RECF
         DE(IC) = DE(IC)*RECF
      END IF
!*****ADD RECOIL ENERGY FOR NEUTRON EMISSION (BUT NOT FISSION)
      IF(IS.EQ.6.AND. IRTIC.NE.6) THEN
         RECF = AWR/(AWR-1.0)
         XX(IC) = XX(IC)*RECF
         DXX(IC) = DXX(IC)*RECF
         EE(IC) = EE(IC)*RECF
         DE(IC) = DE(IC)*RECF
      END IF
!*****SET POINTER TO END OF SPECTRA FOR THIS RADIATION TYPE
      NER(IS) = IC
!*****CHECK CONTINUUM SPECTRUM COVARIANCES
      LCOV = L2
      IF(LCOV.GT.0)   THEN
         EPL = X(1)
         EPU = X(NP)
         CALL RDLIST
         IF(Y(1).NE.EPL.OR.Y(NPL-1).NE.EPU)   THEN
            EMESS = 'SPECTRA COVARIANCES DO NOT COVER COMPLETE ENERGY '
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         CALL TEST5Y(1,NPL,2,1)
         IF(Y(NPL).NE.0.)  THEN
            EMESS = 'LAST VALUE OF F(K) MUST BE 0.0'
            CALL ERROR_MESSAGE(NSEQP1+((NPL+5)/6))
         END IF
      END IF
  100 CONTINUE
!
!     SET TOTAL NUMBER OF SPECTRA OF ALL TYPES
!
      NEPS = IC
      IF(NEPS.LE.0)    GO TO 1000
!
!     CHECK THAT SUM OF ALL RADIATIONS IS THE EFFECTIVE Q
!
      CALL SUMCK('TOTAL ENERGY RELEASE',NEPS,QQ,XX,DQQ,DXX,NSEQH)
!
!     CHECK BE,GE,AE
!
      NSEQ2 = NSEQH + 2
      DO III=1,3
         IF(III.EQ.3) THEN
            JJJ = III + 6
         ELSE
            JJJ = III + 1
         END IF
         XX(III) = ER(JJJ)
         DXX(III) = DER(JJJ)
      END DO
      CALL SUMCK('BETA ENERGY (BE)',3,BE,XX,DBE,DXX,NSEQ2)
      DO III=1,2
         IF(III.EQ.2) THEN
            JJJ = III + 8
         ELSE
            JJJ = III
         END IF
         XX(III) = ER(JJJ)
         DXX(III) = DER(JJJ)
      END DO
      CALL SUMCK('GAMMA ENERGY (GE)',2,GE,XX,DGE,DXX,NSEQ2)
      DO III=1,4
         JJJ = III + 4
         XX(III) = ER(JJJ)
         DXX(III) = DER(JJJ)
      END DO
      CALL SUMCK('PARTICLE ENERGY (AE)',4,AE,XX,DAE,DXX,NSEQ2)
!
!     CHECK SPECTRA PARAMETERS
!
      K1 = 1
      DO I=1,NSPTP
         II = I - 1
!*****SEE IF ANY SPECTRA FOR THIS TYPE
         K2 = NER(I)
         IF(K1.LE.K2) THEN
!*****SET NUMBER OF SPECTRA AND OFFSET FOR THIS TYPE
            N = K2 - K1 + 1
            KK = K1 - 1
            K1 = K2 + 1
!*****CHECK MULTIPLICITY FOR THIS RADIATION TYPE
            IF(II.GT.0.AND.II.LT.8.AND.II.NE.6) THEN
               DO K=1,N
                  KKK = K + KK
                  XX(K) = RI(KKK)
                  DXX(K) = DRI(KKK)
!*****************EXCLUDE NEUTRONS FROM SPONTANEOUS FISSION
                  IF(II.EQ.5.AND.INT(RTYPX(KKK)).EQ.6) THEN
                     XX(K) = 0.0
                     DXX(K) = 0.0
                  END IF
               END DO
               CALL SUMCK(CSTYPE(I)(1:LSTYPE(I))//'MULTIPLICITY',N,     &
     &               BR(II),XX,DBR(II),DXX,0)
            END IF
!*****CHECK AVERAGE ENERGY FOR THIS RADIATION TYPE
            DO K=1,N
               KKK = K + KK
               XX(K) = EE(KKK)*RI(KKK)
               DXX(K) = DE(KKK)*RI(KKK) + EE(KKK)*DRI(KKK)
            END DO
            CALL SUMCK(CSTYPE(I)(1:LSTYPE(I))//'AVERAGE ENERGY',N,      &
     &               ER(I),XX,DER(I),DXX,0)
         END IF
      END DO
!
 1000 RETURN
      END SUBROUTINE CHK457
!
!***********************************************************************
!
      SUBROUTINE CHKDEC
!
!     CHECK DECAY MODES AND AVERAGE ENERGIES FROM RADIOACTIVE DECAY
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT
!
      INTEGER(KIND=I4) :: NDK,IRFS,IR
      INTEGER(KIND=I4) :: NSEQ1,NSEQP2
      INTEGER(KIND=I4) :: I,II
      REAL(KIND=R4) :: SPI,PAR,RTYPI
      REAL(KIND=R4) :: QX,DQC,BRX,DBRX
      REAL(KIND=R4) :: DUM1,DUM2
      REAL(KIND=R4), DIMENSION(2) :: BRT,DBRT
      REAL(KIND=R4), DIMENSION(5) :: XX,DXX
!
!     READ AND TEST HALF LIFE AND AVERAGE ENERGIES
!
      CALL RDLIST
      T12 = C1L
      DT12 = C2L
      CALL TEST6(T12,0.,1.0E+24,'T12')
      CALL TEST6(DT12,0.,T12,'DT12')
      BE = Y(1)
      DBE = Y(2)
      CALL TEST6(DBE,0.,BE,'DBE')
      GE = Y(3)
      DGE = Y(4)
      CALL TEST6(DGE,0.,GE,'DGE')
      AE = Y(5)
      DAE = Y(6)
      CALL TEST6(DAE,0.,AE,'DAE')
      IF(NPL.GT.6) THEN
         DO II=1,4
            XX(II) = Y(2*II+5)
            DXX(II) = Y(2*II+6)
         END DO
         CALL TEST6(DXX(1),0.,XX(1),'DBM')
         CALL TEST6(DXX(2),0.,XX(2),'DBP')
         CALL TEST6(DXX(3),0.,XX(3),'DBAE')
         CALL TEST6(DXX(4),0.,XX(4),'DBCE')
         CALL SUMCK('LIGHT PARTICLE ENERGY',4,BE,XX,DBE,DXX,NSEQ1)
         DO II=1,4
            XX(II) = Y(2*II+13)
            DXX(II) = Y(2*II+14)
         END DO
         CALL TEST6(DXX(1),0.,XX(1),'DGM')
         CALL TEST6(DXX(2),0.,XX(2),'DXR')
         CALL TEST6(DXX(3),0.,XX(3),'DIB')
         CALL TEST6(DXX(4),0.,XX(4),'DAR')
         CALL SUMCK('ELECTROMAGNETIC ENERGY',4,GE,XX,DGE,DXX,NSEQ1)
         DO II=1,5
            XX(II) = Y(2*II+21)
            DXX(II) = Y(2*II+22)
         END DO
         CALL TEST6(DXX(1),0.,XX(1),'DA')
         CALL TEST6(DXX(2),0.,XX(2),'DRE')
         CALL TEST6(DXX(3),0.,XX(3),'DSF')
         CALL TEST6(DXX(4),0.,XX(4),'DN')
         CALL TEST6(DXX(5),0.,XX(5),'DP')
         CALL SUMCK('HEAVY PARTICLE ENERGY',5,AE,XX,DAE,DXX,NSEQ1)
         CALL TEST6(Y(33),0.,Y(34),'DNU')
      END IF
!
!     TEST SPIN AND PARITY
!
      CALL RDLIST
      SPI = C1L
      PAR = C2L
      IF(SPI.NE.SPIUNK)   THEN
         CALL TESTSP(SPI)
         IF(ABS(PAR).NE.1.0) THEN
            WRITE(EMESS,'(A,1PE12.5,A)')                                &
     &             'PARITY=',PAR,' NOT +1.0 OR -1.0'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
      ELSE
         IF(ABS(PAR).NE.0.0) THEN
            EMESS = 'PARITY FIELD MUST BE 0.0 FOR UNKNOWN SPIN'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
      END IF
!
!     INITIALIZE FOR DECAY MODE PROCESSING
!
      QO = 0.
      DQ = 0.
      BR = 0.
      DBR = 0.
      QMAX = 0.
      QQ = 0.0
      DQQ = 0.0
      BRT(1) = 0.0
      DBRT(1) = 0.0
!
!     PROCESS AND TEST DECAY MODES
!
      NDK = N2L
      DO I=1,NDK
         NSEQP2 = NSEQP1 + I
         II = 6*(I-1)
!********CHECK Q AND UNCERTAINTY
         QX = Y(II+3)
         DQC = Y(II+4)
         CALL TEST6(DQC,0.,QX,'DQ')
!********SAVE MAXIMUM Q AS THIS WILL BE TESTED AGAINST MAXIMUM OF GAMMA
!********   STECTRUM
         QMAX = AMAX1(QMAX,QX)
!********CHECK BR AND UNCERTAINTY
         BRX = Y(II+5)
         DBRX = Y(II+6)
         BRT(1) = BRT(1) + BRX
         DBRT(1) = SQRT(DBRT(1)**2+DBRX**2)
         CALL TEST6(DBRX,0.,BRX,'DBR')
!********CALCULATE MAXIMUM ENERGY RELEASE PER DECAY
         QQ = QQ + QX*BRX
         DQQ = DQQ + QX*DBRX + DQC*BRX
!********CHECK FOR VALID DECAY MODE
         DUM1 = 0.
         DUM2 = 0.
         RTYPI = Y(II+1)
         IRFS = IFIX(Y(II+2))
         CALL TESTRT(RTYPI,DUM1,DUM2,0,I)
         IF(MESS.EQ.0)   THEN
!***********GET EFFECTIVE Q AND BR FOR EACH OUTGOING PARTICLE TYPE
   10       IR = IFIX(RTYPI)
!***********************************************************************
!     FOR THIS TEST TO WORK PROPERLY, THE DECAY MODE CARD FOR THE GROUND
!        STATE MUST PRECEDE THE CARD FOR ANY ISOMER
!***********************************************************************
            IF(IRFS.GT.0.AND.QO(IR).GT.0.) THEN
               IF(QX.GT.QO(IR))   THEN
                  EMESS = 'Q(ISOMER).GT.Q(GS)'
                  CALL ERROR_MESSAGE(NSEQP2)
                  WRITE(EMESS,'(A,1PE12.5,A,1PE12.5)')                  &
     &                  'Q(IS)=',QX,'  Q(GS)=',QO(IR)
                  CALL ERROR_MESSAGE(0)
                  NERROR = NERROR + 1
               END IF
            END IF
            QO(IR) = AMAX1(QO(IR),QX)
            DQ(IR) = AMAX1(DQ(IR),DQC)
            IF(IR.EQ.6) THEN
               QO(5) = AMAX1(QO(5),QX)
               DQ(5) = AMAX1(DQ(5),DQC)
            END IF
            BR(IR) = BR(IR) + BRX
            DBR(IR) = SQRT(DBR(IR)**2+DBRX**2)
!***********SEE IF ANY MORE PARTICLES IN THIS DECAY MODE
            IF(ABS(RTYPI-IR).GT.EPSILN3) THEN
               RTYPI = 10.*ABS(RTYPI-IR)
               GO TO 10
            END IF
          END IF
      END DO
!
!     CHECK BRANCHING RATIO SUM
!
      CALL SUMCK('BRANCHING RATIO',1,1.0,BRT,0.0,DBRT,NSEQP1)
!
      RETURN
      END SUBROUTINE CHKDEC
!
!***********************************************************************
!
      SUBROUTINE SUMCK(SYM,N,S,XX,DS,DX,MSEQ)
!
!     PERFORMS SUMMATION CHECK WITHIN STATED UNCERTAINTIES
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: SYM
      INTEGER(KIND=I4) :: N,MSEQ
      REAL(KIND=R4) :: S,DS
      REAL(KIND=R4), DIMENSION(N) :: XX,DX
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: SUMX,SUMDX,DEL,XDEL,RAT,RAT1
!
      SUMX = 0.
      SUMDX = 0.0
      DO I=1,N
         SUMX = SUMX + XX(I)
         SUMDX = SUMDX + DX(I)
      END DO
      XDEL = S - SUMX
      DEL = ABS(XDEL)
      IF(SUMX.NE.0.0)   THEN
         RAT = S/SUMX
      ELSE
         RAT = 0.0
      END IF
      RAT1 = ABS(RAT-1.)
      IF(DEL.LT.EPSILN5)   GO TO 100
      IF(RAT1.LE.EPSILN4)   GO TO 100
      IF(IUNC.EQ.1)   THEN
         IF(DEL.LE.DS.OR.DEL.LE.SUMDX)   GO TO 100
      END IF
!
!     SUM TEST FAILED
!
      IF(MSEQ.GT.0)   THEN
         WRITE(EMESS,'(2A)')  SYM,' SUMUP FAILURE'
         CALL ERROR_MESSAGE(0)
      ELSE
         WRITE(EMESS,'(2A)')  SYM,' SUMUP FAILURE'
         CALL ERROR_MESSAGE(0)
      END IF
      WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)') 'WHOLE=',S,'  SUM=',SUMX
      CALL ERROR_MESSAGE(MSEQ)
!
  100 RETURN
      END SUBROUTINE SUMCK
!
!***********************************************************************
!
      SUBROUTINE AVG(EB,Z,A,NFD,EA,F)
!
!     ROUTINE TO CALCULATE BETA SPECTRUM INTEGRALS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NFD
      REAL(KIND=R4) :: EB,Z,A,EA,F
!
      REAL(KIND=R4), INTRINSIC :: ABS,AMOD,SQRT
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: WRE,ETA0,V,W,B,ETA,DETA,WBV,ARG
      REAL(KIND=R4) :: G,G0,DL,X,Y,ETAPRM,P,Q,S,TWOR,TWOS,GSQ,PHI
!
      REAL(KIND=R4), PARAMETER :: RADIUS0=1.3,SCREEN=1.3
!
      WRE = 1. + EB/EMASS
      ETA0 = SQRT(WRE*WRE-1.)
      V = 1.1*SCREEN*ALPHA*ALPHA*(ABS(Z)**(4.*OTHIRD))
      B = -Z/ABS(Z)
      DETA = ETA0/300.
      G = 0.
      F = 0.
      DL = 0.0
      DO I=2,300
         X = I - 1
         Y = 2. + 2.*AMOD(X,2.)
         ETA = X*DETA
         W = SQRT(ETA*ETA+1.)
         WBV = W + B*V
         ARG = WBV*WBV - 1.
         IF(ARG.GT.0.) THEN
            ETAPRM = SQRT(ARG)
            P = ETA*ETA
            Q = WRE - W
            IF(NFD.LT.1)   THEN
               S = 1.
            ELSE IF (NFD.EQ.1) THEN
               S = P+Q*Q
            ELSE
               S = 1.5*(P*P+Q**4)+ 5.*(P*Q*Q)
            END IF
            X = FNBS(Z,ETAPRM)*SQRT(P/ARG)*((W+B*V)/W)*Q*Y*S
            F = F + X*Q
            G = G + X*Q*W
            DL = DL + 2.*X
         END IF
      END DO
      IF(F.EQ.0.)   THEN
         EA = 0.
         EMESS = 'ERROR CALCULATING BETA SPECTRUM INTEGRAL'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      ELSE
         G = G/F
         DL = DL/F/(EMASS*1.E-6)
         G0 = SQRT(1.-(ALPHA*Z)*(ALPHA*Z))
         TWOS = 2.*(G0-1.)
         TWOR = 0.0051792*RADIUS0*A**OTHIRD
         GSQ = (GAMA(CMPLX(3.+TWOS,0.)))**2
         PHI = 2.*(1.+ G0)*TWOR**TWOS/GSQ
         F = F*DETA/3.*PHI
         EA = EMASS*(G-1.)
      END IF
!
      RETURN
      END SUBROUTINE AVG
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION FNBS(Z,ETA)
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: Z,ETA
!
      REAL(KIND=R4), INTRINSIC :: SQRT,EXP
      COMPLEX(KIND=R4), INTRINSIC :: CABS
!
      REAL(KIND=R4) :: AZ,G0,S,A,ETA2,DEL
!
      AZ = ALPHA*Z
      G0 = SQRT(1.-AZ*AZ)
      S = G0 - 1.
      ETA2 = ETA*ETA
      A = 88./PI/AZ
      IF(ETA2.GE.1./(A*A-1.)) THEN
         DEL= AZ/ETA*SQRT(ETA2+1.)
         A = CABS(GAMA(CMPLX(G0,DEL)))
         FNBS = ETA**(2.*G0)*EXP(PI*DEL)*A*A
      ELSE
         IF(Z.GT.0.) THEN
            FNBS= 2.* PI* AZ* ETA** (1.+ 2.* S)
         ELSE
            FNBS= 0.
         END IF
      END IF
!
      RETURN
      END FUNCTION FNBS
!
!***********************************************************************
!
      COMPLEX(KIND=R4) FUNCTION GAMA(Z)
!
      IMPLICIT NONE
!
      COMPLEX(KIND=R4) :: Z
!
      REAL(KIND=R4), INTRINSIC :: EXP,SIN,COS,REAL,AIMAG,ALOG,ATAN
      COMPLEX(KIND=R4), INTRINSIC :: CMPLX
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: A1,A2,B1,B2,F1,F2,C1,C2,C3,C4,C5,C6,C7
!
      REAL(KIND=R4), PARAMETER :: A=2.5066284
      REAL(KIND=R4), DIMENSION(4), PARAMETER ::                          &
     &    S = (/ 0.083333333, 0.003472222, -0.0026813327,-0.0002294721/)
!
      A1 = REAL(Z)
      A2 = AIMAG(Z)
      GAMA = CMPLX(0.,0.)
      F1 = 1.
      F2 = 0.
      DO WHILE (A1.LT.10.)
         C1 = A1**2 + A2**2
         IF(C1.EQ.0.)    GO TO 100
         C2 = (F1*A1+F2*A2)/C1
         F2 = (F2*A1-F1*A2)/C1
         F1 = C2
         A1 = A1 + 1.
      END DO
      C1 = A1**2 + A2**2
      B1 = A1/C1
      B2 = -A2/C1
      C1 = 1.
      C2 = 0.
      C3 = B1
      C4 = B2
      DO I=1,4
         C1 = C1 + S(I)*C3
         C2 = C2 + S(I)*C4
         C5 = C3*B1 - C4*B2
         C4 = C3*B2 + C4*B1
         C3 = C5
      END DO
      C3 = F1*C1 - F2*C2
      F2 = F1*C2 + F2*C1
      F1 = C3
      C1 = EXP(-A1)
      C2 = COS(-A2)*C1
      C3 = SIN(-A2)*C1
      C4 = 0.5*ALOG(A1*A1+A2*A2)
      C5 = ATAN(A2/A1)
      C6 = (A1-0.5)*C4 - A2*C5
      C7 = (A1-0.5)*C5 + A2*C4
      C1 = EXP(C6)
      C4 = COS(C7)*C1*A
      C5 = SIN(C7)*C1*A
      C1 = F1*C4 - F2*C5
      F2 = F1*C5 + F2*C4
      F1 = C1*C2 - F2*C3
      F2 = C1*C3 + F2*C2
      GAMA = CMPLX(F1, F2)
!
  100 RETURN
      END FUNCTION GAMA
!
!***********************************************************************
!
      SUBROUTINE BETA(ES,IZ,MA,N,NN,EBETA)
!
!     ROUTINE TO COMPUTE MEAN BETA ENERGIES  (+ AND -)
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IZ,MA,N,NN
      REAL(KIND=R4) :: MASS,ES,EBETA
!
      INTEGER(KIND=I4), INTRINSIC :: IABS
      REAL(KIND=R4), INTRINSIC :: SQRT, FLOAT
!
      INTEGER(KIND=I4) :: J,K,JJ
      REAL(KIND=R4) :: E,SN
      REAL(KIND=R4) :: W1,W2,W3,P,X,Y0
      REAL(KIND=R4), DIMENSION(5) :: FI
      REAL(KIND=R4), DIMENSION(10) :: F,W
!
      REAL(KIND=R4), DIMENSION(5) :: T
      DATA T/0.4869533,0.4325317,0.3397048,0.2166977,0.07443717/
      REAL(KIND=R4), DIMENSION(5) :: WF
      DATA WF/0.03333567,0.07472567,0.1095432,0.1346334,0.1477621/
!
      MASS = FLOAT(MA)
      IF(N.EQ.1) THEN
         NZ = IZ + 1
      ELSE IF(N.EQ.2) THEN
         NZ = -(IZ-1)
      END IF
      V0 = 1.13*ALPHA*ALPHA*NZ*IABS(NZ)**OTHIRD
      R0 = 0.5*ALPHA*MASS**OTHIRD
      E = ES/1.0E+06
      IF(N.EQ.2) E = E - EMASS*2.E-6
      W0 = E/(EMASS*1.E-6)+1.0 - V0
      W1 = 0.5*(W0+1.-V0)
      W2 = W0 - 1. + V0
      DO J=1,5
         JJ = 2*J - 1
         W3 = T(J)*W2
         W(JJ) = W1 + W3
         W(JJ+1) = W1 - W3
      END DO
      DO K=1,10
         IF (W(K).LE.1.0) THEN
            P = 1.0
         ELSE
            X = W(K) - 1.0
            P = SQRT(X*(X+2))
         END IF
         Y0 = ALPHA*NZ*W(K)/P
         CALL BSHAPE(NN,P,W(K),Y0,SN)
         F(K) = SN*(W0-W(K))*(W0-W(K))*W(K)/P
      END DO
      DO K=1,2
         FI(K) = 0.0
         IF(K.NE.1) THEN
            DO J=1,10
               F(J) = F(J)*(W(J)+V0)
            END DO
         END IF
         DO J=1,5
            JJ = 2*J - 1
            FI(K) =FI(K)+WF(J)*(F(JJ)+F(JJ+1))
         END DO
      END DO
      EBETA = EMASS*(FI(2)/FI(1)-1.)
!
      RETURN
      END SUBROUTINE BETA
!
!***********************************************************************
!
      SUBROUTINE BSHAPE(N,P,W,Y0,SN)
!
!     THIS SUBROUTINE EVALUATES THE BETA INTENSITY TAKING
!     INTO ACCOUNT THE NATURE OF THE TRANSITION.
!
!     SUBROUTINE ARGUMENTS :-
!                N - NATURE OF TRANSITION: 0=ALLOWED,1=1 UNIQUE
!                P - MOMENTUM OF BETA PARTICLE
!                W - ENERGY OF BETA PARTICLE
!               Y0 - USEFUL RATIO OF ENERGY/MOMENTUM
!               SN - VALUE OF SHAPE FACTOR
!
!     WRITTEN BY:-  A. TOBIAS  : CEGB, BNL
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4)  :: P,W,Y0,SN
!
      REAL(KIND=R4), INTRINSIC :: SQRT, EXP, FLOAT
!
      INTEGER(KIND=I4) :: I,K
      REAL(KIND=R4) :: GN,U,V,X,Q,FF,G,GG,H1,H2,S1,S2,FIV
!
      REAL(KIND=R4), DIMENSION(9) ::  FCTL
      DATA FCTL /1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0/
!
      SN = 0.0
      K = N + 1
      DO I=1,K
         GN = SQRT((I+ALPHA*NZ)*(I-ALPHA*NZ))
         FIV = FLOAT(I)
         V = 2*GN + 1
!
!        EVALUATE GAMMA FUNCTION
!
         CALL GMMMA(V,G)
!
!        APPROXIMATION FOR SMALL VALUES OF W
!
         IF(W.LE.1.0) THEN
            Y0 = ALPHA*NZ
            X = 1.0 - 4.0*Y0*R0/V
            Q = 2*PI*(2*R0*Y0)**(2*GN)/(G*G*Y0)
            FF = X*Y0*Y0/2.0
!
!        RIGOROUS EVALUATION
!
         ELSE
            CALL HYPG(GN,P,Y0,H1,H2)
            CALL LOGAM(GN,Y0,U,V)
            Q = (2*P*R0)**(2*GN)/SQRT(P*P+1)*EXP(PI*Y0+2*U)/(G*G)
            CALL SVAL(ALPHA,GN,FIV,P,Y0,S1,S2)
            FF = (H1*S2+H2*S1)**2*(W-1)
         END IF
         FIV = -FLOAT(I)
!
         IF(W.GT.1.0) THEN
            CALL SVAL(ALPHA,GN,FIV,P,Y0,S1,S2)
            GG = (H1*S1-H2*S2)**2*(W+1)
!
!        APPROXIMATION FOR SMALL VALUES OF W
!
         ELSE
            X = -X*(FIV*(GN-FIV)+Y0*Y0/2.0)
            GG = X - 2.0*R0*Y0*(GN-FIV)/V
         END IF
         FF = Q*(FF+GG)
         FF = FF*FCTL(2*I)/(4**(I-1)*FCTL(I)**2*FCTL(2*K-2*(I-1)))
         FF = FF*(W0-W+V0)**(2*(K-I))*R0**(-2*(I-1))
         SN = SN + FF
      END DO
!
      RETURN
      END SUBROUTINE BSHAPE
!
!***********************************************************************
!
      SUBROUTINE GMMMA(XX,GX)
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: XX,GX
!
      INTEGER(KIND=I4), INTRINSIC :: INT
      REAL(KIND=R4), INTRINSIC :: ABS, FLOAT
!
      REAL(KIND=R4) :: X,Y,GY,ERR
!
      IF(XX.GT.33.0) THEN
         GX =1.0E+33
         GO TO 100
      ELSE
         X = XX
      END IF
      ERR = 1.0E-06
      GX = 1.0
      DO WHILE (X.GT.2)
         X = X - 1.0
         GX = GX*X
      END DO
      IF(X.GT.1.0) THEN
         GO TO 30
      ELSE IF(X.EQ.1.0) THEN
         GO TO 100
      ELSE
         IF(X.LE.ERR) THEN
            Y = FLOAT(INT(X)) - X
            IF(ABS(Y).LE.ERR) THEN
               GO TO 100
            ELSE
               IF (1.0-Y.LE.ERR) THEN
                  GO TO 100
               ELSE
                  IF (X.GT.1.0) GO TO 30
               END IF
            END IF
         END IF
         DO WHILE (X.LE.1.)
            GX = GX/X
            X = X + 1.0
         END DO
      END IF
!
  30  Y = X - 1.0
      GY = 1.0 + Y*(-0.5771017+Y*(0.9858540+Y*(-0.8764218+Y*(0.8328212+ &
     &          Y*(-0.5684729+Y*(0.2548205+Y*(-0.05149930)))))))
      GX = GX*GY
!
  100 RETURN
      END SUBROUTINE GMMMA
!
!***********************************************************************
!
      SUBROUTINE HYPG(GN,P,Y0,H1,H2)
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: GN,P,Y0,H1,H2
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: A,B,X,D1,D2,T1,T2,FI
!
      A = GN
      B = 2*GN
      X = 2*P*R0
      T1 = 1.0
      T2 = 0.0
      H1 = 1.0
      H2 = 0.0
      DO I=1,100
         FI = FLOAT(I)
         A = A+1
         B = B+1
         D1 = T1
         D2 = T2
         T1 = (D1*(-Y0*X)-D2*A*X)/(B*FI)
         T2 = (-Y0*X*D2+A*X*D1)/(B*FI)
         H1 = H1 + T1
         H2 = H2 + T2
         IF(ABS(T1).LT.1.0E-8*ABS(H1)) THEN
            IF (ABS(T2).LE.1.0E-8*ABS(H2)) GO TO 100
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE HYPG
!
!***********************************************************************
!
      SUBROUTINE LOGAM(X,Y,U,V)
!
!  THIS SUBROUTINE COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
!  COMPLEX ARGUMENTS.
!
!          X IS THE REAL PART OF THE ARGUMENT
!          Y IS THE IMAGINARY PART OF THE ARGUMENT
!          U IS THE REAL PART OF THE RESULT
!          V IS THE IMAGINARY PART OF THE RESULT
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: X,Y,U,V
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
      REAL(KIND=R4), INTRINSIC :: ALOG, ATAN
!
      INTEGER(KIND=I4) :: J,IE6
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: B1,B2,B6,B7,X2,E4,E5,T,T1,T2,T3,T4,T5
!
      REAL(KIND=R4), DIMENSION(7) :: H
      DATA H/2.269488974,1.517473649,1.011523068,5.256064690E-1,        &
     &       2.523809524E-1,3.3333333E-2,8.3333333E-2/
!
      B1 = 0.0
      B2 = 0.0
      J = 2
      X2 = X
!
!   X IS NEGATIVE
!
   10 IF(X.LT.0.0) THEN
         E4 = 0.0
         E5 = 0.0
         IE6 = 0
   15    E4 = E4 + 0.5*(ALOG(X*X+Y*Y))
         E5 = E5 + ATAN(Y/X)
         IE6 = IE6 + 1
         X = X + 1.0
         IF(X.LT.0.) GO TO 15
         IF(MOD(IE6,2).NE.0) THEN
            E5 = E5 + PI
         END IF
         GO TO 10
!
!  X IS ZERO
!
      ELSE IF(X.EQ.0.0) THEN
         T = 0.0
         IF(Y.LT.0.0) THEN
            B6 = -PI/2.
         ELSE IF(Y.GT.0.0) THEN
            B6 = PI/2.
         ELSE
            GO TO 90
         END IF
         GO TO 20
      END IF
!
!     X IS POSITIVE
!
      B6 = ATAN(Y/X)
      T = X*X
   20 B7= Y*Y + T
!
!     REAL PART OF LOG
!
      T1 = 0.5*ALOG(B7)
      IF(X.LE.2.0) THEN
         B1 = B1 + B6
         B2 = B2 + T1
         X = X + 1.0
         J = 1
         GO TO 10
      END IF
      T3 = -Y*B6 + (T1*(X-0.5)-X+9.189385332E-1)
      T2 = B6*(X-0.5) + Y*T1 - Y
      T4 = X
      T5 = -Y
      T1 = B7
      DO I=1,7
         T = H(I)/T1
         T4 = T*T4 + X
         T5 = -(T*T5+Y)
         T1 = T4*T4 + T5*T5
      END DO
      T3 = T4 - X+ T3
      T2 = -T5 -Y + T2
      IF(J.EQ.1) THEN
         T3 = T3 - B2
         T2 = T2 - B1
      END IF
      IF(X2.GE.0.) THEN
         U = T3
         V = T2
         X = X2
      ELSE
         U = T3 - E4
         V = T2 - E5
         X = X2
      END IF
      GO TO 100
!
!     X=0.0 AND Y=0.0
!
   90 U = 0.0
      V = 0.0
!
  100 RETURN
      END SUBROUTINE LOGAM
!
!***********************************************************************
!
      SUBROUTINE SVAL(ALPHA,GN,FN,P,Y0,S1,S2)
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: ALPHA,GN,FN,P,Y0,S1,S2
!
      REAL(KIND=R4), INTRINSIC :: SQRT, ATAN2
!
      REAL(KIND=R4) :: D,X,FR,FI,R,A
!
      D = GN*GN + Y0*Y0
      X = ALPHA*NZ/P
      FR = (GN*(-FN)+X*Y0)/D
      FI = (Y0*FN+GN*X)/D
      R = SQRT(SQRT(FR*FR+FI*FI))*SQRT(D)
      A = ATAN2(FI,FR)/2+ATAN2(Y0,GN)-P*R0
      S1 = R*COS(A)
      S2 = R*SIN(A)
!
      RETURN
      END SUBROUTINE SVAL
!
!***********************************************************************
!
      SUBROUTINE ECAP(ES,N,IZ,EBIND,EBREM)
!
!     THIS SUBROUTINE COMPUTES K,L,M SHELL CAPTURE RATIOS
!     AND DERIVES MEAN BINDING ENERGY
!     AND (IF REQUIRED) INTERNAL BREMSSTRAHLUNG ENERGY
!
!     IT IS BASED ON SUBROUTINE ECAP FROM THE COGEND CODE
!
!     ENERGIES ARE IN EV
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N,IZ
      REAL(KIND=R4) :: ES,EBIND,EBREM
!
      REAL(KIND=R4) :: BL1K,BLL1,EBNDK,EBNDL1,EBNDL2,EBNDL3,EBNDM
      REAL(KIND=R4) :: QK,QL1,QL2,QL3,R,RL1,RL2,RL3,RM
!
!     CALCULATE INTERNAL BREMSSTRAHLUNG ENERGY (IF REQUIRED)
!
      IF (IBREM.EQ.1) THEN
         EBREM=0.15E-9*ES*ES
      ELSE
         EBREM = 0.0
      END IF
!
!     CALCULATE CAPTURE FRACTIONS FOR K,L,M SHELLS
!
      BL1K = 0.0
      BLL1 = 0.0
      IF(BX(1,IZ).NE.0.0) THEN
         BL1K = BX(2,IZ)/BX(1,IZ)
      END IF
      IF(BX(2,IZ).NE.0.0) THEN
         BLL1 = BX(3,IZ)/BX(2,IZ)
      END IF
      EBNDK = XLEV(1,IZ-1)
      EBNDL1 = XLEV(2,IZ-1)
      EBNDL2 = XLEV(3,IZ-1)
      EBNDL3 = XLEV(4,IZ-1)
      EBNDM = XLEV(5,IZ-1)
      QK = ES - EBNDK - EBREM
      QL1 = ES - EBNDL1 - EBREM
      QL2 = ES - EBNDL2 - EBREM
      QL3 = ES - EBNDL3 - EBREM
      IF (N .NE. 0) THEN
         RL3 = RL3*(QL3/QL1)**(2*N)/QL1/QL1
      ELSE
         RL3 = 0.333333*N*(2*N+1)*RDENS(3,IZ)*BLL1
      END IF
      RL2 = BLL1*RDENS(2,IZ)*(QL2/QL1)**(2*N+2)
      RL1 = BL1K*RDENS(1,IZ)*(QL1/QK)**(2*N+2)
      R = 1.0 + RL1*(1.0+RL2+RL3+RDENS(4,IZ))
      R = 1.0/R
      RL1 = RL1*R
      RL2 = RL2*RL1
      RL3 = RL3*RL1
      RM = 1.0 - R - RL1 - RL2 - RL3
!
!     CALCULATE MEAN BINDING ENERGY
!
      EBIND = EBNDK*R + EBNDL1*RL1 + EBNDL2*RL2 + EBNDL3*RL3 + EBNDM*RM
!
      RETURN
      END SUBROUTINE ECAP
!
!***********************************************************************
!
      SUBROUTINE TESTRT(ARTYP,ZAN,ZAP,IFL,INC)
!
!     CHECK THE RELATIONSHIP AND VALIDITY OF RTYPE AND OPTIONALLY
!      VERIFYS THE RELATIONSHIP BETWEEN RTYPE, ZAN AND ZAP
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IFL,INC
      REAL(KIND=R4) :: ARTYP,ZAN,ZAP
!
      CHARACTER(LEN=10) :: RTYPC
      CHARACTER(LEN=2) :: ITEG
      INTEGER(KIND=I4) :: IHI,IND,IMESS,NEQ
      INTEGER(KIND=I4) :: I,II
      REAL(KIND=R4) :: RTYP
!
!     CALCULATE SEQUENCE NUMBER
!
      MESS = 0
      NEQ = NSEQP1 + INC
!
!     SEPARATE INDIVIDUAL DIGITS
!
      RTYP = ARTYP + EPSILN6
      WRITE(RTYPC,'(F10.6)')  RTYP
      DO II=9,5,-1
         ITEG = ' '
         ITEG(2:2) = RTYPC(II:II)
         IF(ITEG.NE.' 0')    GO TO 10
      END DO
      II = 4
!
!     TEST FIRST ELEMENT IN RTYP
!
   10 ITEG = RTYPC(2:3)
      READ(ITEG,'(I2)')  IHI
      IND = IHI + 1
      IMESS = 1
      IF((IHI.GT.0.AND.IHI.LT.8).OR.IHI.EQ.10)   THEN
         IMESS = 0
      ELSE IF(IHI.EQ.0.AND.MT.NE.457)   THEN
         IMESS = 0
      END IF
      IF(IMESS.EQ.1) THEN
         WRITE(EMESS,'(I2,A,1PE12.5,A)')                                &
     &       IHI,' IN RTYPE = ',RTYP,' IS INVALID       NEAR '
         CALL ERROR_MESSAGE(NEQ)
         MESS = 1
      END IF
!
!     SCAN THROUGH REMAINING ELEMENTS IN RTYP
!
      IF(II.GE.5)  THEN
         DO I=5,II
            ITEG = '  '
            ITEG(2:2) = RTYPC(I:I)
            READ(ITEG,'(I2)') IHI
            IMESS = 1
            IF(IHI.GT.0.AND.IHI.LT.8)   THEN
               IMESS = 0
            ELSE IF(IHI.EQ.0.AND.MT.NE.457)   THEN
               IMESS = 0
            END IF
            IF(IMESS.EQ.1) THEN
               WRITE(EMESS,'(I2,A,1PE12.5,A)')                          &
     &             IHI,' IN RTYPE = ',RTYP,' IS INVALID       NEAR '
               CALL ERROR_MESSAGE(NEQ)
               MESS = 1
            END IF
         END DO
      END IF
!
!     CHECK CONSISTENCY OF DECAY TYPE, PARENT AND DECAY PRODUCT
!
      IF(IFL.NE.0)   THEN
         IMESS = 0
         IF(IND.EQ.1.OR.IND.EQ.4) THEN
            IF(ZAN.NE.ZAP) IMESS = 1
         ELSE IF(IND.EQ.2) THEN
            IF(ZAN.NE.(ZAP+1.0E+03))   IMESS = 1
         ELSE IF (IND.EQ.3) THEN
            IF(ZAN.NE.(ZAP-1.0E+03))   IMESS = 1
         ELSE IF(IND.EQ.5) THEN
            IF(ZAN.NE.(ZAP-2.004E+03))   IMESS = 1
         ELSE IF(IND.EQ.6) THEN
            IF(ZAN.NE.(ZAP-1.0))   IMESS = 1
         ELSE IF(IND.EQ.7.OR.IND.EQ.11) THEN
            IF(ZAN.NE.0.)   IMESS = 1
         ELSE IF(IND.EQ.8) THEN
            IF(ZAN.NE.(ZAP-1.001E+03))   IMESS = 1
         END IF
         IF(IMESS.EQ.1) THEN
            EMESS = 'ZAN-ZAP RELATIONAL TEST FAILED              NEAR'
            CALL ERROR_MESSAGE(NEQ)
         END IF
      END IF
!
      RETURN
      END SUBROUTINE TESTRT
!
!***********************************************************************
!
      SUBROUTINE CHK_RPD(IA)
!
!     CHECK PRODUCTION OF RATIOACTIVE PRODUCTS DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IA
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: NO,NS
      INTEGER(KIND=I4) :: IZAP,LMF,LFSO
      INTEGER(KIND=I4) :: INCR
      INTEGER(KIND=I4) :: N,J
      REAL(KIND=R4) :: ZAP,HL,RTYP,ZAN
      REAL(KIND=R4) :: SSUM,CT
!
      NO = N2H
!
!     LOOP THRU FINAL STATES
!
      NS = N1H
      DO N=1,NS
         IF(NO.EQ.0)  THEN
            CALL RDLIST
            ZAP = C1L
            LMF = L1L
            LFSO = L2L
         ELSE
            CALL RDCONT
            ZAP = C1H
            LMF = L1H
            LFSO = L2H
         END IF
!
!        SAVE POINTER TO FILE CONTAINING PRODUCTION DATA
!
         NLMF = NLMF + 1
         IF(NLMF.GE.SZLMF) THEN
            STOP 'FIZCON ERROR - SZLMF limit exceeded'
         END IF
         IZAP = IFIX(ZAP)
         LMFS(1,NLMF) = MT
         LMFS(2,NLMF) = LMF
         LMFS(3,NLMF) = IZAP
         LMFS(4,NLMF) = LFSO
!
!        CHECK DECAY INFORMATION
!
         IF(NO.NE.1)   THEN
            SSUM = 0.0
            INCR = 1
            DO J=1,NPL,6
!**************HALF-LIFE TEST
               HL = Y(J)
               CALL TEST6(HL,0.0,1.0E+24,'HL')
!**************CK DECAY MODE VS NEXT NUCLIDE IN THE DECAY CHAIN
               RTYP = Y(J+1)
               ZAN = Y(J+2)
               CALL TESTRT(RTYP,ZAN,ZAP,IA,INCR)
!**************CHAIN TERMINATOR
               CT = Y(J+5)
               CALL TEST6(CT,1.0,3.0,'CT')
               INCR = INCR + 1
               SSUM = SSUM + Y(J+3)
            END DO
!
!           CHECK SUM OF BRANCHING RATIOS
!
            IF(ABS(SSUM-1.0).GT.EPSILN3)    THEN
               WRITE(EMESS,'(A,F10.6,A)')                               &
     &            'SUM OF BRANCHING RATIOS =',SSUM,' IN ERROR'
               CALL ERROR_MESSAGE(NSEQP)
               EMESS = '    SHOULD BE WITHIN RANGE 1.000 +/- .001'
               CALL ERROR_MESSAGE(0)
            END IF
         END IF
      END DO
!
      RETURN
      END SUBROUTINE CHK_RPD
!
!***********************************************************************
!
      SUBROUTINE CKF9
!
!     CHECK FILE 9 AND FILE 10 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NS,LFSO,IZAP,IZAPT
      INTEGER(KIND=I4) :: IZ,IA,IZA
      INTEGER(KIND=I4) :: M,N,NMTX
      REAL(KIND=R4) :: Q,QM
      REAL(KIND=R4) :: ELO,EHI,ELOPR,EHIPR
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     INITIALIZE SUMUP TEST
!
      IF(FIZCON_DATA%ISUM.EQ.1)   CALL SUMPAR(-1)
!
!     GET MF=3 Q VALUE FOR THIS REACTION
!
      DO NMTX=1,NMT3
         IF(MT.EQ.MT3(NMTX))   GO TO 10
      END DO
      NMTX = 0
!
!     PROCESS ALL SUBSECTIONS
!
   10 NS = N1H
      DO N=1,NS
!
!        READ IN SUBSECTION
!
         CALL RDTAB1
         Q = C2
         IF(NFOR.GE.6) then
            QM = C1
         ELSE
            QM = Q
         END IF
         IZAP = L1
!
!        TEST PRODUCT SPECIFICATION
!
         IZA = IFIX(ZA+.001)
         IA = MOD(IZA,1000)
         IZ = IZA/1000
!
!        Check that the value of IZAP is within range
!        Photons IZAP=0) can be the last but not the only
!        If subactinide fission, IZAP=-1 is expected
!
         IF((IZAP.EQ.-1).AND.(MT.NE.18)) THEN
             IF(N.LT.NS .OR. NS.EQ.1) THEN
                CALL TEST6I(IZAP,3000,120000,'IZAP')
             END IF
         END IF
!
!        If IZAP can be uniquely defined from IZ, IA, MT and projectile
!        check for valid IZAP
!
         IZAPT = GET_IZAP(IZ,IA,NSUB/10,MT)
         IF(IZAPT.NE.0 .AND. IZAPT.NE.IZAP) THEN
            WRITE(EMESS,'(A,I6)')                                       &
     &                'IZAP SHOULD BE SET TO ',IZAPT
            CALL ERROR_MESSAGE(NSEQP)
         END IF
!
!        TEST ENERGY RANGE
!
         ELO = X(1)
         EHI = X(NP)
         CALL TESTER(ELO,EHI,Q)
         IF(N.EQ.1)   THEN
            ELOPR = ELO
            EHIPR = EHI
         END IF
!
!        COMPARE Q WITH Q VALUE FROM FILE 3 IF MT EXISTS IN FILE 3
!
         IF(NMTX.GT.0)   THEN
            IF(NFOR.GE.6 .AND. MT.NE.5)  THEN
               IF(QM.NE.QMVAL(NMTX)) THEN
                  WRITE(EMESS,'(A,1PE12.5,A,1PE12.5,A)')                &
     &                'QM=',QM,'  MUST BE ',QMVAL(NMTX),', THE QM '//   &
     &                'FOR THIS SECTION IN FILE 3'
                  CALL ERROR_MESSAGE(0)
               END IF
            ELSE
               IF(Q.GT.QVAL(NMTX))   THEN
                  WRITE(EMESS,'(A,1PE12.5,A,1PE12.5,A)')                &
     &                  'Q=',Q,'  CANNOT BE GREATER THAN ',QVAL(NMTX),  &
     &                  ', THE Q IN FILE 3'
                  CALL ERROR_MESSAGE(0)
               END IF
            END IF
         END IF
!
!        CHECK THAT THIS FINAL PRODUCT AND STATE IS DEFINED IN FILE 8
!
         LFSO = L2
         DO M=1,NLMF
            IF(LMFS(2,M).NE.0) THEN
               IF(LMFS(1,M).EQ.MT.AND.LMFS(2,M).EQ.MF)   THEN
                  IF(LMFS(3,M).EQ.IZAP.OR.IZAP.EQ.0) THEN
                     IF(LMFS(4,M).EQ.LFSO)  THEN
                       LMFS(2,M) = 0.0
                       GO TO 60
                     END IF
                  END IF
               END IF
            END IF
         END DO
         WRITE(EMESS,'(A,I3,A)')                                        &
     &              'SUBSECTION',N,' NOT DESCRIBED IN FILE 8'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
!
!        SUM MULTIPLICITIES AT ENERGIES
!
   60    IF(ITEST.NE.0)   CALL SUMPAR(N)
!
      END DO
!
!     SAVE ENERGY RANGE SPANNED
!
      CALL STORF(MF,MT,ELOPR,EHIPR)
!
!     DO SUMUP TEST IF POSSIBLE
!
      IF(ITEST.NE.0.AND.NS.GT.1)  CALL SUMPAR(0)
!
      RETURN
      END SUBROUTINE CKF9
!
!***********************************************************************
!
      INTEGER(KIND=I4) FUNCTION GET_IZAP(IZ,IA,IPZA,MT)
!
!     FUNCTION TO CALCULATE THE PRODUCT IZA FROM THE TARGET, PROJECTILE
!       AND REACTION (MT).
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MT,IZ,IA,IPZA
!
      INTEGER(KIND=I4) :: N
!
!     REACTION PRODUCTS
!
      INTEGER(KIND=I4), PARAMETER :: NRECS=36
      INTEGER(KIND=I4), PARAMETER, DIMENSION(NRECS) :: MTS =            &
     &       (/ 4,11,16,17,22, 23,24,25,28,29, 30,32,33,34, 35, 36,37,  &
     &         41,42,44,45,102,103,104,105,106,107,108,109,111,112,113, &
     &         114,115,116,117/)
      INTEGER(KIND=I4), PARAMETER, DIMENSION(NRECS) :: DZ =             &
     &       (/ 0,-1, 0, 0,-2, -6,-2,-2,-1,-4, -4,-1,-1,-2, -5, -5, 0,  &
     &         -1,-1,-2,-3,  0, -1, -1, -1, -2, -2, -4, -6, -2, -3, -5, &
     &         -5, -2, -2, -3/)
      INTEGER(KIND=I4), PARAMETER, DIMENSION(NRECS) :: DA =             &
     &       (/-1,-4,-2,-3,-5,-13,-6,-7,-2,-9,-10,-3,-4,-4,-11,-12, -4, &
     &         -3,-4,-3,-6,  0, -1, -2, -3, -3, -4, -8,-12, -2, -5,-11, &
     &         -10, -3, -4, -6/)
!
      GET_IZAP = 0
      DO N=1,NRECS
         IF(MT.EQ.MTS(N)) THEN
            GET_IZAP = 1000*(IZ+DZ(N)) + (IA+DA(N)) + IPZA
            GO TO 100
         END IF
      END DO
!
  100 RETURN
      END FUNCTION GET_IZAP
!
!***********************************************************************
!
      SUBROUTINE CKF12
!
!     CHECK FILE 12 AND FILE 13 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: LO,NK,LF,LG
      INTEGER(KIND=I4) :: NE1,NE2,LG1
      INTEGER(KIND=I4) :: IPLACE,MTT,NCON
      INTEGER(KIND=I4) :: JJ,J2,N
      REAL(KIND=R4) :: EG,ES
      REAL(KIND=R4) :: ELO,EHI,EMIN,EMAX,SSUM,EST,TEST
      REAL(KIND=R4) :: ELEVI,ETEST1,ETEST2,YJ1,YJ2
      REAL(KIND=R4), DIMENSION(2) :: X2
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     MAKE SURE THAT SECTION DOES NOT APPEAR IN MF = 13
!
      IPLACE = NPMT + 1
      IF(NPMT.NE.0.AND.MF.EQ.13) THEN
         DO N=1,NPMT
            MTT = MOD(ICON(N,1),1000)
            IF(MTT.EQ.MT) THEN
               IPLACE = N
               NCON = ICON(N,2)
               WRITE(EMESS,'(A,I4,A)')                                  &
     &              'MT=',MT,' CANNOT EXIST IN FILE 12 AND FILE 13'
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
               GO TO 20
            END IF
         END DO
      END IF
!
!     SAVE INDEX VALUE
!
      NPMT = IPLACE
      ICON(IPLACE,1) = MT + 1000*MF
      ICON(IPLACE,2) = 0
      NCON = 0
!
!     FIND REACTION IMPLIED LEVEL ENERGY
!
   20 ELEVI = 0.0
      DO N=1,NMT3
         IF(MT3(N).EQ.MT) THEN
            ELEVI = QMVAL(N) - QVAL(N)
            GO TO 25
         END IF
      END DO
!
!     BRANCH FOR TRANSITION PROBABILITIES
!
   25 LO = L1H
      IF(LO.EQ.2)   GO TO 50
!
!     INITIALIZE
!
      NK = N1H
!
!     Save the number of photons
!
      MMTGAM=MMTGAM+1
      MMGAM(MMTGAM)=MT
      NNGAM(MMTGAM)=NK
!
!     INITIALIZE SUMUP TEST
!
      IF(FIZCON_DATA%ISUM.NE.0.AND.NK.GT.1)  CALL SUMPAR(-1)
!
!     PROCESS ALL PARTIALS
!
      DO J2=1,NK
         CALL RDTAB1
         IF(J2.EQ.1)  THEN
            ELO = X(1)
            EHI = X(NP)
            X2(1) = BIGNO
            EMIN = BIGNO
            EMAX = 0.0
            IF(NK.GT.1) THEN
               IF(ITEST.NE.0)  CALL STOR(NP,2)
               CALL RDTAB1
            END IF
         END IF
         LF = L2
         EG = C1
         ES = C2
         IF(LF.EQ.1)   THEN
            NCON = NCON + 1
            NMTGAM = NMTGAM + 1
            EGAM(NMTGAM) = X(1)
            MTGAM(NMTGAM) = MT
         ELSE
            CALL TESTEG(EG,ES,ELEVI)
         END IF
!
!        CHECK FOR DECREASING GAMMA RAY ENERGY ORDER
!
         X2(2) = C1
         CALL TEST5B(X2,2,1)
         X2(1) = X2(2)
!
!        SET RANGE OF ALL PARTIALS
!
         EMIN = AMIN1(EMIN,X(1))
         EMAX = AMAX1(EMAX,X(NP))
!********DO SUMUP TEST
         IF(ITEST.NE.0)    CALL SUMPAR(J2)
      END DO
!
!     SAVE ENERGY RANGE SPAN
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     SEE THAT DATA SPAN THE SAME RANGE AS IN FILE 3
!
      IF(MT.NE.460) CALL ISFIL(MF,3,MT,MT)
!
!     MAKE SURE ONLY ONE CONTINUUM
!
      ICON(NPMT,2) = NCON
      IF(NCON.GT.1)  THEN
         WRITE(EMESS,'(A,I3)')                                          &
     &       'ONLY ONE CONTINUUM SUB-SECTION ALLOWED IN FILE'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
!     TESTS WHEN THERE IS A TOTAL YIELD RECORD
!
      IF(NK.GE.2)  THEN
!********CHECK THAT PARTIALS SPAN SAME ENERGY RANGE AS TOTAL
         ETEST1 = (ELO-EMIN)/ELO
         ETEST2 = (EHI-EMAX)/EHI
         IF(ETEST1.GT.EPSILN4.OR.ABS(ETEST2).GT.EPSILN4)  THEN
            EMESS = ' ENERGY RANGE OF TOTAL YIELD TABLE DOES NOT '
            CALL ERROR_MESSAGE(0)
            WRITE(EMESS,'(4X,A,I6)')                                    &
     &          'SPAN THE COMBINED RANGE OF THE PARTIAL TABLES, MT=',MT
            CALL ERROR_MESSAGE(0)
            NERROR = NERROR + 1
         END IF
!********DO FINAL SUMUP TEST
         IF(ITEST.NE.0)   CALL SUMPAR(0)
      END IF
      GO TO 100
!
!     TRANSITION PROBABILITIES
!
   50 LG = L2H
      SSUM = 0.0
      CALL RDLIST
      X2(1) = C1L
!
!     Save the number of gamma transitions
!
!
!     Save the number of photons
!
      MMTGAM=MMTGAM+1
      MMGAM(MMTGAM)=MT
      NNGAM(MMTGAM)=N2L
!
!     COMPARE HIGHEST LEVEL ENERGY WITH Q VALUE
!
      EST = ELEVI
      TEST = ABS((C1L-EST)/C1L)
      IF(TEST.GE.EPSILN4)   THEN
         WRITE(EMESS,'(A,1PE11.4)')                                     &
     &        'THE ENERGY OF HIGHEST LEVEL MUST BE',EST
         CALL ERROR_MESSAGE(NSEQP1)
      END IF
!
!     PROCESS ALL TRANSITION MATRIX ELEMENTS
!
      NE1 = 0
      NE2 = 0
      LG1 = LG + 1
      DO JJ=1,NPL,LG1
!
!        CHECK DECREASING ORDER OF LEVEL ENERGIES
!
         X2(2) = Y(JJ)
         CALL TEST5B(X2,2,1)
         X2(1) = X2(2)
!
!        CHECK VALUE OF MATRIX ELEMENTS
!
         YJ1 = Y(JJ+1)
         IF(YJ1.LE.0.0.OR.YJ1.GT.1.0)       NE1 = NE1 + 1
         IF(LG.NE.1) THEN
            YJ2 = Y(JJ+2)
            IF(YJ2.LE.0.0.OR.YJ2.GT.1.0)       NE2 = NE2 + 1
         END IF
         SSUM = SSUM + YJ1
      END DO
!
!     TRANSITION MATRIX ELEMENTS SHOULD SUM TO ONE
!
      IF(ABS(SSUM-1.0).GT.EPSILN4) THEN
         WRITE(EMESS,'(A,1PE11.4,A,I4)')                                &
     &      'SUM OF TRANSITION PROBABILITIES=',SSUM,' FOR MT=',MT
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
!     INDICATE IF ANY ELEMENTS FOUND OUT OF RANGE
!
      IF(NE1.NE.0)   THEN
         WRITE(EMESS,'(I4,A)')                                          &
     &      NE1,' TRANSITION PROBABILITIES NOT IN RANGE 0.0 TO 1.0'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
      IF(NE2.NE.0)   THEN
         WRITE(EMESS,'(I4,A)')                                          &
     &         NE2,' CONDITIONAL TRANSITION PROBABILITIES NOT IN '//    &
     &        'RANGE 0.0 TO 1.0'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
  100 RETURN
      END SUBROUTINE CKF12
!
!***********************************************************************
!
      SUBROUTINE CKF14
!
!     CHECK FILE 14 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: LI,NI,LTT,NK,NE
      INTEGER(KIND=I4) :: NCON,ICONT,NBEG,NLMOD,NCONT
      INTEGER(KIND=I4) :: MTT,MTL,JINT2,MF1
      INTEGER(KIND=I4) :: J,M,N
      REAL(KIND=R4) :: ELO,EHI,ELEVI,EG,ES,TESTQ,ELO1,XDUM
      REAL(KIND=R4), DIMENSION(2) :: X2,XG
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     BRANCH IF ALL DISTRIBUTIONS FOR ALL GAMMAS ARE ISOTROPIC
!
      LI = L1H
      IF(LI.EQ.1) GO TO 60
!
!     Check that the number of gammas matches the number in MF12
!
      NK = N1H
      DO J=1,MMTGAM
        IF(MMGAM(J).EQ.MT) THEN
          IF(NNGAM(J).NE.NK) THEN
            CALL TEST3(NK,NNGAM(J),'NK')
          END IF
          GO TO 20
        END IF
      END DO
!     Complementary MT not found in MF12/MF13
      EMESS = 'Complementary section not found in MF12,13'
      CALL ERROR_MESSAGE(NSEQP1)
!
!     INITIALIZE
!
   20 NI = N2H
      LTT = L2H
      ELO = BIGNO
      EHI = 0.
      XG(1) = BIGNO
      XG(2) = 0.0
      NCON = 0
!
!     GET IMPLIED LEVEL ENERGY
!
      ELEVI = 0.0
      DO N=1,NMT3
         IF(MT3(N).EQ.MT) THEN
            ELEVI = QMVAL(N) - QVAL(N)
            GO TO 10
         END IF
      END DO
!
!     DETERMINE IF A CONTINUUM OR DISCRETE CHANNEL REACTION
!
   10 IF(MT.GE.50.AND.MT.LT.91)  THEN
         ICONT = 0
      ELSE
         ICONT = 1
      END IF
      IF(NFOR.GE.6)   THEN
         IF(MT.LT.600.OR.MT.GT.849)    GO TO 25
         NBEG = 600
         NLMOD = 50
         NCONT = 49
      ELSE
         IF(MT.LT.699.OR.MT.GT.799)    GO TO 25
         NBEG = 700
         NLMOD = 20
         NCONT = 18
      END IF
      MTT = MT - NBEG
      MTL = MOD(MTT,NLMOD)
      IF(MTL.LT.NCONT)   ICONT = 0
!
!     PROCESS EACH SUBSECTION
!
   25 NK = N1H
      DO 50 N=1,NK
!
!     ISOTROPIC DISTRIBUTION FOR THIS GAMMA RAY
!
      IF(N.LE.NI) THEN
         CALL RDCONT
         EG = C1H
         ES = C2H
         XG(2) = EG
         IF(ICONT.NE.0)  THEN
!***********CONTINUUM
            CALL TESTEG(EG,ES,ELEVI)
         ELSE
!***********DISCRETE CHANNEL
            TESTQ = (ELEVI-ES)/ELEVI
            IF(ABS(TESTQ).GT.EPSILN5)   THEN
               EMESS = 'SOURCE LEVEL ENERGY MUST BE EQUAL TO -QI'
               CALL ERROR_MESSAGE(NSEQP1)
            END IF
         END IF
         GO TO 30
      END IF
!
!     GAMMA RAY NOT ISOTROPIC
!
      CALL RDTAB2
      EG = C12
      ES = C22
      NE = NP2
      NR = NR2
      XG(2) = EG
      X2(1) = -BIGNO
      IF(ICONT.NE.0)   THEN
!********CONTINUUM
         CALL TESTEG(EG,ES,ELEVI)
      ELSE
!********DISCRETE CHANNEL
         TESTQ = (ELEVI-EG)/ELEVI
         IF(ABS(TESTQ).GT.EPSILN5)   THEN
            EMESS = 'SOURCE LEVEL ENERGY MUST BE EQUAL TO -QI'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
      END IF
!
!     CHECK FOR ACCEPTABLE INTERPOLATION CODE
!
      DO J=1,NR
         JINT2 = JNT2(J)
         IF(JINT2.LT.1.OR.JINT2.GT.3)   THEN
            WRITE(EMESS,'(A,I4,A,I4,A)')                                &
     &          'SUBSECTION',N,' INTERPOLATION TYPE ',JINT2,' INVALID'
            CALL ERROR_MESSAGE(NSEQP)
         END IF
      END DO
!
!     PROCESS EACH INCIDENT ENERGY
!
      DO M=1,NE
!
!        LEGENDRE REPRESENTATION
!
         IF(LTT.NE.2)   THEN
            CALL RDLIST
            IF(NPL.GT.2.AND.MOD(NPL,2).NE.0)   THEN
               WRITE(EMESS,'(A,I3)') 'NL=',NPL,' SHOULD BE EVEN'
               CALL ERROR_MESSAGE(NSEQP1)
            END IF
            NP = NPL
            CALL TEST6Y(-1.0,1.0,'FL')
            X2(2) = C2L
         ELSE
!
!        TABULAR REPRESENTATION
!
            CALL RDTAB1
!
!           CHECK INTERPOLATION CODE
!
            DO J=1,NR
               JINT2 = JNT(J)
               IF(JINT2.LT.1.OR.JINT2.GT.2)   THEN
                  WRITE(EMESS,'(A,I4,A,I4,A)')                          &
     &                 'SUBSECTION',N,' INTERPOLATION TYPE ',           &
     &                  JINT2,' INVALID'
                  CALL ERROR_MESSAGE(NSEQP)
               END IF
            END DO
            CALL TEST6X (-1.0,1.0,'MU')
            CALL TEST7(XDUM,1)
            X2(2) = C2
         END IF
!
!        CHECK THAT ENERGIES ARE IN INCREASING ORDER
!
         CALL TEST5(X2,2,1)
         IF(M.EQ.1)   ELO1 = X2(2)
         X2(1) = X2(2)
      END DO
      ELO = AMIN1(ELO,ELO1)
      EHI = AMAX1(EHI,X2(2))
!
!     UP CONTINUUM COUNT IF A CONTINUUM
!
   30 IF(XG(2).EQ.0.0)   NCON = NCON + 1
!
!     TEST THAT GAMMA ENERGIES IN DECREASING ORDER
!
      CALL TEST5B(XG,2,1)
      XG(1) = XG(2)
!*****RESET TEST AT END OF ISOTROPIC GAMMAS
      IF(N.EQ.NI) XG(1) = BIGNO
   50 CONTINUE
!
!     CHECK THAT ONLY ONE CONTINUUM DISTRIBUTION EXISTS
!
      IF(NCON.GT.1) THEN
         EMESS = 'ONLY ONE CONTINUUM SUB-SECTION ALLOWED '
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
!     SAVE ENERGY RANGE SPANNED
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     CHECK THAT RANGE SPANNED IS SAME AS THAT IN 12 OR 13
!
   60 IF(NPMT.GT.0)   THEN
         DO N=1,NPMT
            MTT = MOD(ICON(N,1),1000)
            IF(MT.EQ.MTT) THEN
               MF1 = ICON(N,1)/1000
               ICON(N,2) = ICON(N,2) + 2000
               GO TO 70
            END IF
         END DO
         MF1 = MF - 2
   70    CALL ISFIL(MF,MF1,MT,MT)
      END IF
!
      RETURN
      END SUBROUTINE CKF14
!
!***********************************************************************
!
      SUBROUTINE CKF15
!
!     CHECK FILE 15 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: NC,NE,NM
      INTEGER(KIND=I4) :: ICKTT,JINT2,MTT,MF1
      INTEGER(KIND=I4) :: I,J,N
      REAL(KIND=R4) :: ELO,EHI
      REAL(KIND=R4) :: XDUM
      REAL(KIND=R4), DIMENSION(2) :: X2
!
!     AVOID DEVIANT POINT TEST ON THIS FILE
!
      ICKTT = FIZCON_DATA%ICKT
      FIZCON_DATA%ICKT = 1
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     PROCESS EACH DISTRIBUTION
!
      NC = N1H
      ELO = BIGNO
      EHI = 0.0
      DO I=1,NC
         CALL RDTAB1
!********ALL PROBABILITIES MUST BE ONE SINCE ONLY ONE DISTRIBUTION
!********IS ALLOWED
         CALL TEST6Y(1.0,1.0,'PKE')
         CALL RDTAB2
!
!       PROCESS EACH INCIDENT ENERGY
!
         NE = NP2
         X2(1) = -BIGNO
         DO NM=1,NE
            CALL RDTAB1
            DO J=1,NR
               JINT2 = JNT(J)
               IF(JINT2.LT.1.OR.JINT2.GT.3)   THEN
                  WRITE(EMESS,'(A,I4,A,I3,A)')                          &
     &                 'SUBSECTION',NM,' INTERPOLATION TYPE',JINT2,     &
     &                    ' INVALID'
                  CALL ERROR_MESSAGE(NSEQP1)
               END IF
            END DO
!***********CHECK NORMALIZATION AT THIS ENERGY
            CALL TEST7(XDUM,1)
            IF(ELO.GT.C2) ELO = C2
            IF(EHI.LT.C2) EHI = C2
!***********INCIDENT ENERGIES MUST BE IN INCREASING ORDER
            X2(2) = C2
            CALL TEST5(X2,2,1)
            X2(1) = X2(2)
         END DO
      END DO
!
!     SAVE ENERGY RANGE SPANNED
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     CHECK THAT RANGE SPANNED IS SAME AS IN FILE 12 OR 13
!
      IF(NPMT.GT.0) THEN
         DO N=1,NPMT
            MTT = MOD(ICON(N,1),1000)
            IF(MT.EQ.MTT) THEN
               MF1 = ICON(N,1)/1000
               ICON(N,2) = ICON(N,2) - MOD(ICON(N,2),1000)
               CALL ISFIL(MF,MF1,MT,MT)
               GO TO 70
            END IF
         END DO
       END IF
!
!     RESTORE DEVIANT POINT CHECK FLAG
!
  70  FIZCON_DATA%ICKT = ICKTT
!
      RETURN
      END SUBROUTINE CKF15
!
!***********************************************************************
!
      SUBROUTINE TESTEG(EG,ES,EST)
!
!     TESTS MADE ON DISCRETE PHOTON ENERGY(EG) AND SOURCE ENERGY(ES)
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: EG,ES,EST
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      REAL(KIND=R4) :: TEST,ENMAXM
!
      ENMAXM = ENMAX*1.0E-06
      IF(EG.LT.0.0)   THEN
         EMESS = 'PHOTON ENERGY CANNOT BE NEGATIVE'
         CALL ERROR_MESSAGE(NSEQP1)
      END IF
      IF(ES.GT.ENMAX)  THEN
         WRITE(EMESS,'(A,F4.1,A)')                                      &
     &      'SOURCE LEVEL ENERGY CANNOT EXCEED',ENMAXM,' MEV'
         CALL ERROR_MESSAGE(NSEQP1)
      END IF
      IF(MT.LT.50.OR.(MT.GT.90.AND.MT.LT.600).OR.MT.GT.849) THEN
         IF(EG.NE.0.0.AND.EG.LT.1.0E+03)  THEN
            WRITE(EMESS,'(A,1PE12.5,A)')                                &
     &         'PHOTON ENERGY(EG=',EG,') SHOULD BE >1 KEV'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         IF(ES.LT.0.)   THEN
            WRITE(EMESS,'(A,1PE12.5,A)')                                &
     &          'SOURCE LEVEL(ES=',ES,') MUST BE >=0.'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         IF(ES.NE.0.0.AND.EG.GT.ES)    THEN
            EMESS = 'GAMMA RAY ENERGY CANNOT EXCEED THE SOURCE LEVEL'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
      ELSE
         IF(EG.LT.1.0E+03)    THEN
            WRITE(EMESS,'(A,1PE12.5,A)')                                &
     &         'PHOTON ENERGY(EG=',EG,') SHOULD BE >1 KEV'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         TEST = ABS((ES-EST)/ES)
         IF(TEST.GT.EPSILN3)     THEN
            WRITE(EMESS,'(A,1PE12.5)')                                  &
     &        'SOURCE LEVEL ENERGY(ES) SHOULD BE',EST
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         IF(ES.LT.1.0E+04)     THEN
            WRITE(EMESS,'(A,1PE12.5,A)')                                &
     &            'SOURCE LEVEL(ES=',ES,') SHOULD BE >10 KEV'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         IF(EG.GT.ES)   THEN
            EMESS = 'GAMMA RAY ENERGY CANNOT EXCEED THE SOURCE LEVEL'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
      END IF
!
      RETURN
      END SUBROUTINE TESTEG
!
!***********************************************************************
!
      SUBROUTINE CKF23
!
!     CHECK FILE 23 DATA
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: ELO,EHI
!
      REAL(KIND=R4), PARAMETER :: EDMAX=150000.
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     INITIALIZE FOR SUMUP TEST FIRST TIME
!
      IF(ITEST.EQ.0.AND.FIZCON_DATA%ISUM.GT.0)   CALL SUMGAM(-1)
!
!     READ DATA
!
      CALL RDTAB1
!
!     DO SUMUP FOR THIS MT IF NEEDED
!
      IF(FIZCON_DATA%ISUM.GT.0)  CALL SUMGAM(MT)
!
!     SAVE ENERGY RANGE SPAN
!
      ELO = X(1)
      EHI = X(NP)
      CALL STORF(MF,MT,ELO,EHI)
!
!     CHECK THAT FILE GOES TO ENDF MAXIMUM
!
      CALL TESTER(ELO,EHI,QUNK)
!
!     CHECK PHOTOELECTRIC EDGE AND FLUORESCENCE YIELD
!
      IF(MT.GE.534)  THEN
         IF(C1.LE.0.) THEN
            EMESS = 'P.E. EDGE ENERGY MUST BE GREATER THAN ZERO'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         IF(C1.GT.EDMAX) THEN
            EMESS = 'P.E. EDGE ENERGY MUST BE LESS THAN 150 KEV'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         IF(C2.LT.0.) THEN
            EMESS = 'FLUORESCENCE YIELD MUST NOT BE NEGATIVE   '
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
         IF(C2.GT.C1) THEN
            EMESS = 'FLUORESCENCE YIELD CANNOT EXCEDE THE P.E. EDGE'
            CALL ERROR_MESSAGE(NSEQP1)
         END IF
      END IF
!
      RETURN
      END SUBROUTINE CKF23
!
!***********************************************************************
!
      SUBROUTINE CKF26
!
!     ROUTINE TO CHECK FILE 26 DATA
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: ABS,FLOAT
!
      INTEGER(KIND=I4) :: NK,LAW,LANG,NE,NEP
      INTEGER(KIND=I4) :: NSEQH,INTS,NW,NREPT,ICHKER
      INTEGER(KIND=I4) :: I,J,L,N
      REAL(KIND=R4) :: ZAP,ZAPT
      REAL(KIND=R4) :: E,EONE,ENE
      REAL(KIND=R4) :: ELO,EHI,ELOS,EHIS
      REAL(KIND=R4) :: ANS,ANS1,XL,XU,YL,YU
      REAL(KIND=R4), DIMENSION(2) :: X2
!
      REAL(KIND=R4), PARAMETER :: PERR=5.0*EPSILN4
!
      ELO = BIGNO
      EHI = 0.0
!
!     SEE IF SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     LOOP OVER LAWS
!
      NK = N1H
      DO I=1,NK
         CALL RDTAB1
         NSEQH = NSEQP + 1
         IF(MT.EQ.526) THEN
            ZAP = C1
            ZAPT = FLOAT(NSUB/10)
            CALL TEST3F(ZAP,ZAPT,'ZAP')
         END IF
         ELOS = X(1)
         IF(ELO.GT.ELOS) ELO = ELOS
         EHIS = X(NP)
         IF(EHI.LT.EHIS) EHI = EHIS
         LAW = L2
         IF(MT.EQ.526) THEN
            IF(LAW.NE.2) THEN
               EMESS = 'ONLY LAW = 2 ALLOWED IN MT = 526'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
         ELSE IF (MT.EQ.527) THEN
            IF(LAW.EQ.2) THEN
               EMESS = 'LAW = 2 NOT ALLOWED IN MT = 527'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
         ELSE IF(MT.EQ.528) THEN
            IF(LAW.NE.8) THEN
               EMESS = 'ONLY LAW = 8 ALLOWED IN MT = 528'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
         END IF
!
!        CONTINUUM DISTRIBUTION LAW
!
         IF(LAW.EQ.1) THEN
            CALL RDTAB2
            LANG = L12
            IF(MT.EQ.527.AND.LANG.NE.1) THEN
               EMESS = 'ONLY LANG = 1 FOR LAW = 1 ALLOWED IN MT = 527'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
            INTS = L22
            X2(1) = -BIGNO
            NE = NP2
            DO N=1,NE
               CALL RDLIST
               E = C2L
               IF(N.EQ.1)  THEN
                  EONE = E
               ELSE IF(N.EQ.NE)  THEN
                  ENE = E
               END IF
!**************TEST FOR INCREASING ENERGY ORDER
               X2(2) = E
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
               CALL TEST3(L2L,0,'NA')
!**************TEST THAT E-PRIME IS IN INCREASING ORDER
               NEP = N2L
               NW = NPL
               NREPT = NW/NEP
               CALL TEST5Y(1,NW,NREPT,1)
!**************TEST NORMALIZATION INTEGRAL
               ANS = 0.0
               IF(2.LE.NEP) THEN
                  DO J=2,NEP
                     L = NREPT*(J-2) + 1
                     XL = Y(L)
                     XU = Y(L+NREPT)
                     YL = Y(L+1)
                     YU = Y(L+NREPT+1)
                     CALL ECSI(XL,YL,XU,YU,XL,XU,INTS,ANS1)
                     ANS = ANS + ANS1
                  END DO
               END IF
               IF(ABS(ANS-1.0).GT.PERR) THEN
                  WRITE(EMESS,'(A,F11.6,A,1PE11.4)')                    &
     &                  'CHECK NORMALIZATION=',ANS,' AT E=',E
                  CALL ERROR_MESSAGE(NSEQP1)
               END IF
            END DO
            ICHKER = 1
!
!        TWO-BODY ANGULAR DISTRIBUTION LAW
!
         ELSE IF(LAW.EQ.2) THEN
            NE = NP2
            X2(1) = -BIGNO
            CALL RDTAB2
            LANG = L12
            IF(MT.EQ.527.AND.LANG.NE.1) THEN
               EMESS = 'ONLY LANG = 12 FOR LAW = 2 ALLOWED IN MT = 526'
               CALL ERROR_MESSAGE(NSEQP)
            END IF
            NE = NP2
            DO N=1,NE
               CALL RDLIST
               E = C2L
               IF(N.EQ.1)  THEN
                  EONE = E
               ELSE IF(N.EQ.NE)  THEN
                  ENE = E
               END IF
!**************TEST FOR INCREASING ENERGY ORDER
               X2(2) = C2L
               CALL TEST5(X2,2,1)
               X2(1) = X2(2)
            END DO
            ICHKER = 1
!
!        ENERGY TRANSFER FOR EXCITATION
!
         ELSE IF(LAW.EQ.8) THEN
            CALL RDTAB1
            ICHKER = 0
         END IF
!
!        CHECK LAW DATA COVERS SAME RANGE AS PROBABILITY
!
         IF(ICHKER.EQ.1) THEN
            IF(EONE.NE.ELOS.OR.ENE.NE.EHIS)  THEN
               EMESS = 'ENERGY RANGE FOR DISTRIBUTIONS IN LIST RECORDS'
               CALL ERROR_MESSAGE(0)
               EMESS = '    INCONSISTENT WITH TAB1 RECORD'
               CALL ERROR_MESSAGE(NSEQH)
            END IF
         END IF
      END DO
!
!     SAVE ENERGY RANGE SPANNED
!
      CALL STORF(MF,MT,ELO,EHI)
!
!     ENERGY RANGE SPANNED MUST BE SAME AS FILE 3
!
      CALL ISFIL(MF,23,MT,MT)
!
      RETURN
      END SUBROUTINE CKF26
!
!***********************************************************************
!
      SUBROUTINE CKF27
!
!     CHECK FILE 27 DATA
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: ELO,EHI
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     READ DATA
!
      CALL RDTAB1
!
!     SAVE ENERGY RANGE SPAN
!
      ELO = X(1)
      EHI = X(NP)
      CALL STORF(MF,MT,ELO,EHI)
!
!     SEE THAT DATA SPAN THE SAME RANGE AS IN FILE 23
!     -- This test is deemed unnecessary by the author of the EPDL library
!     IF(MT.EQ.505.OR.MT.EQ.506) THEN
!        CALL ISFIL(MF,23,MT,502)
!     END IF
!
      RETURN
      END SUBROUTINE CKF27
!
!***********************************************************************
!
      SUBROUTINE CKF28
!
!     ROUTINE TO CHECK FILE 28 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NSS
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: SUBI,SUBIP
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     LOOP OVER SUBSHELLS
!
      NSS = N1H
      DO I=1,NSS
         CALL RDLIST
         SUBI = C1L
         IF(I.EQ.1) THEN
            SUBIP = -1.
         ELSE
            IF(SUBI.LE.SUBIP) THEN
               EMESS = 'SUBSHELLS OUT OF ORDER'
               CALL ERROR_MESSAGE(NSEQP1)
            ELSE
               SUBIP = SUBI
            END IF
         END IF
         CALL SUB_CHECK
      END DO
!
      RETURN
      END SUBROUTINE CKF28
!
!***********************************************************************
!
      SUBROUTINE SUB_CHECK
!
!     CHECK ORDER OF SECONDARY AND TERTIARY SUBSHELLS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NTR,NRAD,NSS,NST
      INTEGER(KIND=I4) :: J
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4), DIMENSION(200) ::  SS,ST
!
      NTR = N2L
      IF(NTR.EQ.0) GO TO 100
      J = 7
      NRAD = NTR
      DO N=1,NTR
         IF(Y(J+1).NE.0.0) THEN
            NRAD = N - 1
            GO TO 10
         ELSE
            SS(N) = Y(J)
         END IF
      END DO
   10 IF(NRAD.GT.1) CALL TEST5(SS,NRAD,1)
      NSS = 0
      IF(NRAD.LT.NTR-1) THEN
         J = 7
         DO N=NRAD+1,NTR
            IF(N.EQ.NRAD+1) THEN
               NSS = 1
               SS(1) = Y(J)
               NST = 1
               ST(1) = Y(J+1)
            ELSE
               IF(Y(J).EQ.SS(NSS)) THEN
                  NST = NST + 1
                  ST(NST) = Y(J+1)
               ELSE
                  IF(NST.GT.1) CALL TEST5(ST,NST,1)
                  NSS = NSS + 1
                  SS(NSS) = Y(J)
                  NST = 1
                  ST(1) = Y(J+1)
               END IF
            END IF
            J = J + 6
         END DO
      END IF
!
      IF(NSS.GT.1) CALL TEST5(SS,NSS,1)
!
  100 RETURN
      END SUBROUTINE SUB_CHECK
!
!***********************************************************************
!
      SUBROUTINE CKF32
!
!     CHECK FILE 32 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NIS,NER,LRU,LRF,NROO,NIT,LCOMP,ISR,NLS
      INTEGER(KIND=I4) :: NSRS,NLRS,NLSA,NDIGIT
      INTEGER(KIND=I4) :: NI,N,NN,I1,NM,NNN,NJSX,I3
      INTEGER(KIND=I4) :: NSEQP2
!
!     CHECK THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     PROCESS ALL ISOTOPES
!
      NIS = N1H
      DO NI=1,NIS
         CALL RDCONT
!
!        PROCESS ALLL ENERGY RANGES
!
         NER = N1H
         DO N=1,NER
            CALL RDCONT
            LRU = L1H
            LRF = L2H
!
!           RESOLVED ENERGY REGION
!
            IF(LRU.NE.2)   THEN
!
!              PROCESS ANY ENERGY DEPENDENT SCATTERING RADIUS
!
               NROO = N1H
               IF(NROO.NE.0)  THEN
                  CALL RDCONT
                  NIT = N2H
                  DO NN=1,NIT
                     CALL RDLIST
                  END DO
               END IF
!
!              READ RECORD WITH NUMBER OF PARTIAL WAVES
!
               CALL RDCONT
               LCOMP = L2H
               NLS   = N1H
               ISR   = N2H
!
!              ENDF-5 FORMAT
!
               IF(LCOMP.EQ.0)  THEN
!                 Process scattering radius uncertainty for LCOMP=0
                  IF(ISR.EQ.1) CALL RDCONT
                  CALL CKF32_V5(NLS)
!
!              ENDF-6 FORMAT
!
               ELSE IF(LCOMP.EQ.1) THEN
!                 Process scattering radius uncertainty for LCOMP=1
                  IF(ISR.EQ.1) THEN
                    IF(LRF.EQ.1 .OR. LRF.EQ.2) THEN
                      CALL RDCONT
                    ELSE IF(LRF.EQ.3 .OR. LRF.EQ.7) THEN
                      CALL RDLIST
                    END IF
                  END IF
                  CALL RDCONT
                  NSRS = N1H
                  NLRS = N2H
!
!                 SHORT RANGE CORRELATIONS
!
                  IF(NSRS.NE.0)  THEN
                     DO I1=1,NSRS
                        IF(LRF.EQ.7) THEN
                            CALL RDCONT
                            NJSX=L1H
                            DO I3=1,NJSX
                                CALL RDLIST
                            END DO
                            CALL RDLIST
                        ELSE
                            CALL RDLIST
                        END IF
                     END DO
!
!                 LONG RANGE CORRELATIONS
!
                  ELSE
                     IF(NLRS.NE.0)  THEN
                        DO I1=1,NLRS
                           CALL RDLIST
                        END DO
                     END IF
                  END IF
!
!              Compact format representation
!
               ELSE IF(LCOMP.EQ.2) THEN
                 NLSA=N1H
!                Process scattering radius uncertainty for LCOMP=2
                 IF(ISR.EQ.1) THEN
                   IF(LRF.EQ.1 .OR. LRF.EQ.2) THEN
                     CALL RDCONT
                   ELSE IF(LRF.EQ.3 .OR. LRF.EQ.7) THEN
                     CALL RDLIST
                   END IF
                 END IF
                 CALL RDLIST
                 IF(LRF.EQ.7) THEN
                    DO I1=1,NLS
                       CALL RDLIST
                       CALL RDLIST
                    END DO
                 END IF
                 CALL RDCONT
                 NDIGIT=L1H
!                If undefined, assume NDIGIT=2 to conform
!                with older convention
                 IF(NDIGIT.LE.0) NDIGIT=2
                 NNN   =L2H
                 NM    =N1H
                 DO I1=1,NM
!                   CALL RDTEXT
                    CALL RDINTG(NDIGIT)
                 END DO
               ELSE
!                Invalid LCOMP
               END IF
!
!           UNRESOLVED RESONANCE REGION
!
            ELSE
               CALL RDCONT
               NLS = N1H
               IF(NLS.GT.MF2URL) THEN
               NSEQP2 = NSEQP1
                 WRITE(EMESS,'(4X,A,I4,A,I4,A)')                        &
     &             'Number of l-states',NLS,' in MF32 >'
     &             ,MF2URL,' in MF 2'                                   &
                 CALL ERROR_MESSAGE(NSEQP2)
               END IF
               DO NN=1,NLS
                  CALL RDLIST
               END DO
               CALL RDLIST
            END IF
         END DO
      END DO
!
      RETURN
      END SUBROUTINE CKF32
!
!***********************************************************************
!
      SUBROUTINE CKF32_V5(NLS)
!
!     CHECK FILE 32 DATA IN ENDF-5 FORMAT
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NLS
!
      REAL(KIND=R4), INTRINSIC :: ABS,SQRT
!
      INTEGER(KIND=I4) :: NRS
      INTEGER(KIND=I4) :: LOC,JMIN,JMAX,NSEQP2,IMESS
      INTEGER(KIND=I4) :: NL,NN,I,J
      REAL(KIND=R4) :: ER,YLOC,CORR
!
      REAL(KIND=R4), DIMENSION(5) :: VAR
      CHARACTER(LEN=2), DIMENSION(5), PARAMETER ::                      &
     &         NAME = (/'ER','GN','GG','GF','AJ'/)
!
!     PROCESS EACH L VALUE FOR ENDF-5 STYLE FORMAT
!
      DO NL=1,NLS
         CALL RDLIST
         NRS = N2L
         LOC = 1
!
!        LOOP OVER RESONANCE ENERGIES.
!
         DO NN=1,NRS
            ER = Y(LOC)
            LOC = LOC + 5
!
!           LOCATE AND TEST VARIANCES.
!
            DO I=1,5
               JMAX = I
               JMIN = MIN0(JMAX,2)
               DO J=JMIN,JMAX
                  LOC = LOC + 1
                  IF(J.GE.I) THEN
                     VAR(I) = Y(LOC)
                     IF(VAR(I).LT.0..OR.VAR(I).GT.BIGNO) THEN
                        NSEQP2 = NSEQP1 + (LOC+5)/6
                        EMESS = 'VARIANCE INCORRECT'
                        CALL WARNING_MESSAGE(NSEQP2)
                        WRITE(EMESS,'(4X,3A,1PE12.5,A,1PE12.5)')        &
     &                        'VAR(',NAME(I),')=',VAR(I),'  ER=',ER
                        CALL WARNING_MESSAGE(0)
                     END IF
                  END IF
               END DO
            END DO
!
!           LOCATE AND TEST CORRELATION COEFFICIENTS.
!
            LOC = LOC - 9
            DO I=3,5
               DO J=2,I
                  LOC = LOC + 1
                  IF(J.NE.I)   THEN
                     YLOC = Y(LOC)
                     IF(YLOC.NE.0.) THEN
                        IMESS = 0
                        IF(VAR(I).EQ.0..OR.VAR(J).EQ.0.) THEN
                           CORR = BIGNO
                           IMESS = 1
                        ELSE IF(VAR(I).GT.0..AND.VAR(J).GT.0.) THEN
                           CORR = YLOC/SQRT(VAR(I)*VAR(J))
                           IF(ABS(CORR).GT.(1.+EPSILN4)) IMESS = 1
                        END IF
                        IF(IMESS.EQ.1) THEN
                           NSEQP2 = NSEQP1 + (LOC+5)/6
                           EMESS = 'CORRELATION COEFFICIENT INCORRECT'
                           CALL ERROR_MESSAGE(NSEQP2)
                           WRITE(EMESS,'(4X,5A,1PE12.5,A,1PE12.5,A)')   &
     &                           'CORR(',NAME(J),',',NAME(I),')=',      &
     &                            CORR,'  ER=',ER
                           CALL ERROR_MESSAGE(0)
                        END IF
                     END IF
                  END IF
               END DO
            END DO
            LOC = LOC + 2
         END DO
      END DO
!
      RETURN
      END SUBROUTINE CKF32_V5
!
!***********************************************************************
!
      SUBROUTINE CKF33
!
!     CHECK FILE 31 AND FILE 33 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NL,MAT1,MT1,NC,NI,MTL
      INTEGER(KIND=I4) :: ILB8,IMIN
      INTEGER(KIND=I4) :: I,I1,I2,I3,J
      REAL(KIND=R4) :: EL,ER,QI
!
      INTEGER(KIND=I4) :: IFRST
      DATA IFRST/1/
!
!     SET SOME LIMITS IN COMMON ON FIRST PASS
!
      IF(IFRST.NE.0)   THEN
         IFRST = 0
         ILB8 = 0
      END IF
!
!     CHECK THAT SECTION IS IN INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     PROCESS NON-LUMPED MT
!
      MTL=L2H
      IF(MTL.EQ.0) THEN
!
!        LOOP OVER SUBSECTIONS
!
         NL = N2H
         DO I1=1,NL
            CALL RDCONT
            MAT1 = L1H
            MT1 = L2H
            NC = N1H
            NI = N2H
            NIX = 0
!
!           NC TYPE SUB-SUBSECTIONS
!
            IF(NC.NE.0) THEN
               DO I2=1,NC
                  CALL RDCONT
                  CALL RDLIST
                  CALL NCTEST(I2)
               END DO
            END IF
!
!           NI TYPE SUB-SUBSECTIONS
!
            IF(NI.NE.0) THEN
               DO I3=1,NI
                  CALL RDLIST
                  IF(L2L.EQ.8) ILB8 = 1
                  CALL NITEST(MAT1,MT1)
               END DO
            END IF
         END DO
!
!        COMPARE NC AND NI GRIDS TO LOCATE OCCURENCES OF ILLEGAL
!        OVERLAP BETWEEN THEM.  NC GRID CONTAINS CONTRIBUTIONS
!        ONLY FROM LTY=0 SUB-SUBSECTIONS.
!
         IF(NCX.NE.NCXLAS.AND.NIX.NE.0) THEN
            IMIN = NCXLAS + 1
            DO I=IMIN,NCX
               DO J=1,NIX
                  ER = AMIN1(EI(J,2),EC(I,2))
                  EL = AMAX1(EI(J,1),EC(I,1))
                  IF(ER.GE.(1.+EPSILN4)*EL) THEN
                     WRITE(EMESS,'(A,I3)')                              &
     &                  'CONFLICT BETWEEN NC AND NI GRIDS FOR MT=',MT
                     CALL ERROR_MESSAGE(0)
                     WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')            &
     &                  'EC(I,1)=',EC(I,1),' EC(I,2)=',EC(I,2)
                     CALL ERROR_MESSAGE(0)
                     WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')            &
     &                  'EI(J,1)=',EI(J,1),' EI(J,2)=',EI(J,2)
                     CALL ERROR_MESSAGE(0)
                     NERROR = NERROR + 1
                  END IF
               END DO
            END DO
         END IF
!
!        CHECK THAT AN LB=8 SUB-SUBSECTION EXISTS
!        (REMOVED APRIL 2006)
!        IF(NFOR.GE.6.AND.ILB8.EQ.0)  THEN
!              EMESS = 'REQUIRED SUB-SUBSECTION WITH LB=8 IS MISSING'
!              CALL ERROR_MESSAGE( 0)
!        END IF
      ELSE
!
!     PROCESS CONTRIBUTORS TO LUMPED MT
!
!       Check if MTL is in the reaction index list - add if needed
        DO I=1,NMT3
          IF(MT .EQ.MT3(I)) THEN
            QI=QVAL(I)
          END IF
          IF(MTL.EQ.MT3(I)) THEN
!           Set Q-value as the largest of the contributing reactions
            QVAL(I)=MAX(QI,QVAL(I))
            GOTO 100
          END IF
        END DO
!       Add MTL and its Q-value to the reaction list
        NMT3=NMT3+1
        MT3(NMT3)=MTL
        QVAL(NMT3)=QI
  100   CONTINUE
      END IF
!
      NCXLAS = NCX
!
      RETURN
      END SUBROUTINE CKF33
!
!***********************************************************************
!
      SUBROUTINE CKF34
!
!     CHECK FILE 34 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NMT1,NSS,NI
      INTEGER(KIND=I4) :: J,N,NN
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     LOOP OVER SUBSECTIONS
!
      NMT1 = N2H
      DO N=1,NMT1
         CALL RDCONT
         NSS = N1H*N2H
         IF(MT.EQ.L2H) THEN
            NSS = N1H*(N1H+1)/2
         ELSE
            NSS = N1H*N2H
         END IF
!
!        LOOP OVER SUB-SUBSECTIONS
!
         DO NN=1,NSS
            CALL RDCONT
            NI = N2H
!
!           LOOP OVER NI-TYPE SUB-SUB-SUBSECTIONS
!
            DO J=1,NI
               CALL RDLIST
            END DO
!
         END DO
      END DO
!
      RETURN
      END SUBROUTINE CKF34
!
!***********************************************************************
!
      SUBROUTINE CKF35
!
!     CHECK FILE 35 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NK
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: ELO,EHI
!
!     TEST THAT SECTION IS IN THE INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     LOOP OVER SUBSECTIONS
!
      NK = N1H
      DO N=1,NK
         CALL RDLIST
         IF(N.EQ. 1) ELO = C1L
         IF(N.EQ.NK) EHI = C2L
      END DO
      CALL STORF(MF,MT,ELO,EHI)
!
!     SEE THAT DATA SPAN THE SAME RANGE AS IN FILE 5
!
      CALL ISFIL(MF,5,MT,MT)
!
      RETURN
      END SUBROUTINE CKF35
!
!***********************************************************************
!
      SUBROUTINE CKF40
!
!     CHECK FILE 40 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NS,NL,MAT1,MT1,NC,NI,IZAP
      INTEGER(KIND=I4) :: IMIN
      INTEGER(KIND=I4) :: N,I,I1,I2,I3,J
      REAL(KIND=R4) :: ER,EL
!
      INTEGER(KIND=I4) :: IFRST
      DATA IFRST/1/
!
!     SET SOME LIMITS IN COMMON ON FIRST PASS
!
      IF(IFRST.NE.0)   THEN
         IFRST = 0
      END IF
!
!     CHECK THAT SECTION IS IN INDEX
!
      CALL TESTD(1000*MF+MT)
!
!     LOOP OVER FINAL STATES
!
      NS = N1H
      DO N=1,NS
         CALL RDCONT
!
!        Check that the value of IZAP is within range
!        Photons IZAP=0) can be the last but not the only
!
         IZAP=L1H
         IF(N.LT.NS .OR. NS.EQ.1) CALL TEST6I(IZAP,3000,120000,'IZAP')
!
!        LOOP OVER SUBSECTIONS
!
         NL = N2H
         DO I1=1,NL
            CALL RDCONT
            MAT1 = L1H
            MT1 = L2H
            NC = N1H
            NI = N2H
            NIX = 0
!
!           NC TYPE SUB-SUBSECTIONS
!
            IF(NC.NE.0) THEN
               DO I2=1,NC
                  CALL RDCONT
                  CALL RDLIST
                  CALL NCTEST(I2)
               END DO
            END IF
!
!           NI TYPE SUB-SUBSECTIONS
!
            IF(NI.NE.0) THEN
               DO I3=1,NI
                  CALL RDLIST
                  CALL NITEST(MAT1,MT1)
               END DO
            END IF
!
!           COMPARE NC AND NI GRIDS TO LOCATE OCCURENCES OF ILLEGAL
!           OVERLAP BETWEEN THEM.  NC GRID CONTAINS CONTRIBUTIONS ONLY
!           FROM LTY=0 SUB-SUBSECTIONS.
!
            IF(NCX.NE.NCXLAS.AND.NIX.NE.0) THEN
               IMIN = NCXLAS + 1
               DO I=IMIN,NCX
                  DO J=1,NIX
                     ER = AMIN1(EI(J,2),EC(I,2))
                     EL = AMAX1(EI(J,1),EC(I,1))
                     IF(ER.GE.((1.+EPSILN4)*EL)) THEN
                        WRITE(EMESS,'(A,I3)')                           &
     &                     'CONFLICT BETWEEN NC AND NI GRIDS FOR MT=',MT
                        CALL ERROR_MESSAGE(0)
                        WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')         &
     &                     'EC(I,1)=',EC(I,1),' EC(I,2)=',EC(I,2)
                        CALL ERROR_MESSAGE(0)
                        WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')         &
     &                     'EI(J,1)=',EI(J,1),' EI(J,2)=',EI(J,2)
                        CALL ERROR_MESSAGE(0)
                        NERROR = NERROR + 1
                     END IF
                  END DO
               END DO
            END IF
         END DO
         NCXLAS = NCX
      END DO
!
      RETURN
      END SUBROUTINE CKF40
!
!***********************************************************************
!
      SUBROUTINE NCTEST(II)
!
!     CHECK ENERGY GRIDS OF ALL NC-TYPE SUB-SUBSECTIONS IN ONE
!     SUBSECTION OF FILE 31 OR FILE 33.  IN ADDITION, IF LTY=0,
!     CHECK FOR PRESENCE OF COVARIANCES FOR EACH MTI AND THEN CHECK FOR
!     ILLEGAL NESTING OF DERIVATION FORMULAS.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: II,I
!
      INTEGER(KIND=I4) :: LTY,NCI,MTI
      INTEGER(KIND=I4) :: KMIN,KMAX,NMIN,NMAX,NCXM
      INTEGER(KIND=I4) :: ICX,K,N
      INTEGER(KIND=I4) :: LOC,ISET
      REAL(KIND=R4) :: E2LAST,ER,EL
      SAVE E2LAST
!
!     INITIALIZE ON FIRST PASS
!
      IF(II.EQ.1) E2LAST = 0.
!
!     CHECK THAT E1 IS LESS THAN E2
!
      IF(C1L.GE.C2L) THEN
         EMESS = 'ENERGIES INCORRECT'
         CALL ERROR_MESSAGE(NSEQP1)
         WRITE (EMESS,'(4X,1PE12.5,A,1PE12.5)')                         &
     &       'E1=',C1L,' MUST BE LESS THAN E2=',C2L
         CALL ERROR_MESSAGE(0)
         GO TO 100
      END IF
!
!     TEST FOR NON-OVERLAP AND ASCENDING ORDER.
!
      IF(C1L.LT.E2LAST) THEN
         EMESS = 'ENERGIES INCORRECT'
         CALL ERROR_MESSAGE(NSEQP1)
         WRITE (EMESS,'(4X,1PE12.5,A,1PE12.5)')                         &
     &          'E2LAST=',E2LAST,' AND E1=',C1L
         CALL ERROR_MESSAGE(0)
         GO TO 100
      END IF
!
!     END OF TESTS IF LTY NE 0
!
      LTY = L2H
      IF(LTY.NE.0) GO TO 100
!
!     NO PROBLEMS FOUND.  IF LTY=0, STORE ENERGIES AND MT NUMBERS.
!
      IF(NCX.GE.NCXMAX) THEN
         EMESS = ' '
         CALL ERROR_MESSAGE(0)
         EMESS = ' ***ERROR IN NCTEST.  NCX EXCEEDS NCXMAX.'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
         IERX = 1
         GO TO 100
      END IF
      NCX = NCX + 1
      EC(NCX,1) = C1L
      EC(NCX,2) = C2L
      E2LAST = C2L
      MTLEFT(NCX,1) = MT
      NCI = N2L
      MTLEFT(NCX,2) = NCI
      MTLEFT(NCX,3) = NSEQP1
!
!     SAVE MTI-VALUES AND CHECK FOR PRESENCE OF COVARIANCES FOR EACH MTI
!
      LOC = 0
      DO I=1,NCI
         IF(MTR.GE.MTRMAX) THEN
            EMESS = ' '
            CALL ERROR_MESSAGE(0)
            EMESS = ' ***ERROR IN NCTEST.  MTR EXCEEDS MTRMAX.'
            CALL ERROR_MESSAGE(0)
            NERROR = NERROR + 1
            IERX = 1
            GO TO 100
         END IF
         MTR = MTR + 1
         LOC = LOC + 2
         MTI = IFIX(Y(LOC) + .1)
         MTRITE(MTR) = MTI
         CALL TESTS(1000*MF+MTI,ISET)
         IF(ISET.GE.3)   THEN
            WRITE(EMESS,'(A,I3,A)')                                     &
     &            'MTI=',MTI,' INCORRECT'
            CALL ERROR_MESSAGE(NSEQP1)
            EMESS = 'SELF-COVARIANCES FOR MTI MISSING'
            CALL ERROR_MESSAGE(0)
         END IF
      END DO
!
!     CHECK FOR ILLEGAL NESTING OF DERIVATION FORMULAS.
!
      IF(NCX.EQ.1) GO TO 100
      KMIN = MTR - NCI + 1
      KMAX = MTR
      NMAX = 0
      NCXM = NCX - 1
      DO ICX=1,NCXM
         NMIN = NMAX + 1
         NMAX = NMAX + MTLEFT(ICX,2)
!
!        CHECK FOR FINITE OVERLAP WITH PREVIOUS LTY=0 SUB-SUBSECTION.
!
         ER = EC(ICX,2)
         IF(C2L.LT.ER)   ER = C2L
         EL = EC(ICX,1)
         IF(C1L.GT.EL)    EL = C1L
         IF(ER.GE.(1.+EPSILN4)*EL) THEN
!
!        COMPARE PREVIOUS MT WITH NEW RIGHT HAND SIDE
!           (I.E., MTI VALUES).
!
            DO K=KMIN,KMAX
!**************TYPE 1 CONFLICT DETECTED.
               IF(MTLEFT(ICX,1).EQ.MTRITE(K)) THEN
                  EMESS = 'BAD NC-TYPE SUB-SUBSECTION'
                  CALL ERROR_MESSAGE(NSEQP1)
                  WRITE(EMESS,'(4X,A,I6)')                              &
     &              'CONFLICTS WITH SUB-SUBSECTION AT',MTLEFT(ICX,3)
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(4X,A,I3,A,I3)')                         &
     &                 'DERIVED MT=',MTLEFT(ICX,1),                     &
     &                 ' USED AS AN MTI IN MT=',MT
                  CALL ERROR_MESSAGE(0)
               END IF
            END DO
!
!           COMPARE NEW MT WITH PREVIOUS RIGHT HAND SIDE.
!
            DO N=NMIN,NMAX
!**************TYPE 2 CONFLICT DETECTED.
               IF(MT.EQ.MTRITE(N))  THEN
                  EMESS = 'BAD NC-TYPE SUB-SUBSECTION'
                  CALL ERROR_MESSAGE(NSEQP1)
                  WRITE(EMESS,'(4X,A,I6)')                              &
     &                 'CONFLICTS WITH SUB-SUBSECTION AT',MTLEFT(ICX,3)
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(4X,A,I3,A,I3)')                         &
     &                 'DERIVED MT=',MTLEFT(ICX,1),                     &
     &                  ' USED AS AN MTI IN MT=',MT
                  CALL ERROR_MESSAGE(0)
               END IF
            END DO
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE NCTEST
!
!***********************************************************************
!
      SUBROUTINE NITEST(MAT1,MT1)
!
!     CHECK THE ENERGY GRID OR GRIDS OF ONE NI-TYPE SUB-SUBSECTION OF
!     FILE 31 OR 33 TO SEE IF IT BEGINS AT 1.E-5 EV, ENDS AT 2.E+7 EV,
!     AND IS IN STRICT ASCENDING ORDER.  IN ADDITION, IF MAT1=0 AND
!     MT1=MT, CHECK THAT THE DIAGONAL ELEMENTS (VARIANCES) LIE IN THE
!     RANGE (0.,1.E4) AND THAT THE OFF-DIAGONAL ELEMENTS CORRESPOND
!     TO CORRELATION COEFFICIENTS IN THE RANGE (-1.,1.).  LOCATE THE
!     UPPER AND LOWER BOUNDARIES OF THOSE ENERGY RANGES WHERE THE
!     VARIANCE IS NON-ZERO FOR ONE OR MORE NI-TYPE SUB-SUBSECTIONS
!     OF THE MT1=MT SUBSECTION.  STORE THE BOUNDARIES IN COMMON FOR
!     LATER TESTS.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MAT1,MT1
!
      REAL(KIND=R4), INTRINSIC :: ABS,SQRT
!
      INTEGER(KIND=I4) :: LB,LS
      INTEGER(KIND=I4) :: NBAD,IVR,NEN,LOC,ILOC,JLOC
      INTEGER(KIND=I4) :: NSTEP,N,N1,N2,JMIN
      INTEGER(KIND=I4) :: NSEQP2,NSEQP3
      INTEGER(KIND=I4) :: NMTX,J,I
      REAL(KIND=R4) :: ENIMIN,ENJMIN,Q,YN,YNP
      REAL(KIND=R4) :: ELAST,ENEXT,VRR,VRI,VRJ,YLOC,CORR
!
      REAL(KIND=R4), PARAMETER :: CRIT=1.0E+4,ZERO=0.
!
!     GET MF=3 Q VALUE FOR THIS REACTION
!
      ENIMIN =-1
      ENJMIN =-1
      DO NMTX=1,NMT3
         IF(MT.EQ.MT3(NMTX))   THEN
            Q = QVAL(NMTX)
            ENIMIN = MAX(ENMIN,-Q*(AWR+1.)/AWR)
         ELSE IF(MT1.EQ.MT3(NMTX))   THEN
            Q = QVAL(NMTX)
            ENJMIN = MAX(ENMIN,-Q*(AWR+1.)/AWR)
C...        GO TO 10
         END IF
      END DO
!
!     INITIALIZE
!
   10 NBAD = 0
      IF(MAT1.EQ.0.AND.MT1.EQ.MT) THEN
         IVR = 1
      ELSE
         IVR = 0
      END IF
      NEN = 0
      LB = L2L
!
!     SET LB-DEPENDENT PARAMETERS.
!
      N2=0
      IF(LB.LE.4) THEN
         NSTEP = 2
         N = -1
         N1 = N2L - L1L
         N2 = L1L
      ELSE
         NSTEP = 1
         N = 0
         IF(LB.EQ.5) THEN
            N1 = N2L
            N2 = 0
         ELSE
            N1 = N2L
            IF(LB.EQ.6)  THEN
               N2 = (NPL-1)/N2L
            ELSE
               NSTEP = 2
               N = -1
            END IF
         END IF
      END IF
!
!     BEGIN LOOP OVER NI ENERGY GRID.
!
   20 IF(NBAD.GE.100)   GO TO 100
      NEN = NEN + 1
      N = N + NSTEP
!
!     FINISHED WITH ENERGIES?
!
      IF(NEN.GT.(N1+N2)) GO TO 30
      YN = Y(N)
      YNP = Y(N+1)
      NSEQP2 = NSEQP1 + (N+5)/6
      NSEQP3 = NSEQP1 + N/6 + 1
!
!     CHECK FOR CORRECT FIRST ENERGY.
!
      IF(NEN.EQ.1) THEN
         ELAST = YN
         IF(ABS(YN-ENMIN).GT.(CRIT*ENMIN) .AND. ENMIN.GT.0) THEN
            IF(ABS(YN-ENIMIN).GT.(CRIT*ENIMIN) .AND. ENIMIN.GT.0) THEN
               EMESS = 'ENERGY INCORRECT'
               CALL ERROR_MESSAGE(NSEQP2)
               WRITE (EMESS,'(4X,A,1PE12.5,A,1PE12.5)')                 &
     &              'EXPECT ',ENIMIN,', FIND ',YN
               CALL ERROR_MESSAGE(0)
               NBAD = NBAD + 1
            END IF
         END IF
         GO TO 25
      END IF
      IF(NEN.EQ.(N1+1)) THEN
         ELAST = YN
         IF(ABS(YN-ENMIN).GT.(CRIT*ENMIN)
     &      .AND. ENMIN.GT.0 .AND. MAT1.EQ.0) THEN
            IF(ABS(YN-ENJMIN).GT.(CRIT*ENJMIN) .AND. ENJMIN.GT.0) THEN
               EMESS = 'ENERGY INCORRECT'
               CALL ERROR_MESSAGE(NSEQP2)
               WRITE (EMESS,'(4X,A,1PE12.5,A,1PE12.5)')                 &
     &              'EXPECT ',ENJMIN,', FIND ',YN
               CALL ERROR_MESSAGE(0)
               NBAD = NBAD + 1
            END IF
         END IF
         GO TO 25
      END IF
!
!     CHECK FOR ASCENDING ORDER.
!
      IF(YN.LE.ELAST) THEN
         EMESS = 'ENERGIES OUT OF ORDER'
         CALL ERROR_MESSAGE(NSEQP2)
         WRITE (EMESS,'(4X,A,1PE12.5,A,1PE12.5)')                       &
     &           'ELO=',ELAST,' EHI=',YN
         CALL ERROR_MESSAGE(0)
         NBAD = NBAD + 1
      END IF
      ELAST = YN
!
!     CHECK FOR CORRECT LAST ENERGY.
!
      IF(NEN.EQ.N1.OR.NEN.EQ.(N1+N2)) THEN
         IF(ABS(YN-ENMAX).GE.(CRIT*ENMAX)) THEN
            EMESS = 'ENERGY INCORRECT'
            CALL ERROR_MESSAGE(NSEQP2)
            WRITE (EMESS,'(4X,A,1PE12.5,A,1PE12.5)')                    &
     &           'EXPECT ',ENMAX,', FIND ',YN
            CALL ERROR_MESSAGE(0)
         END IF
!
!     CHECK FOR CORRECT LAST COVARIANCE VALUE.
!
         IF(LB.LE.4.AND.YNP.NE.0.) THEN
            EMESS = 'VARIANCE INCORRECT'
            CALL ERROR_MESSAGE(NSEQP3)
            WRITE (EMESS,'(4X,A,1PE12.5,A,1PE12.5)')                    &
     &         'EXPECT ',ZERO,', FIND ',YNP
            CALL ERROR_MESSAGE(0)
            NBAD = NBAD + 1
         END IF
         GO TO 20
      END IF
!
!     CHECK VARIANCES AND STORE ENERGIES NOW IF LB=0, 1, OR 2.
!     LB=5 WILL BE TREATED LATER.
!
   25 IF(LB.GE.5.OR.IVR.EQ.0) GO TO 20
      IF(LB.EQ.3.OR.LB.EQ.4) THEN
         IF(NEN.EQ.1) THEN
            WRITE(EMESS,'(A,I2,A)')                                     &
     &         'NOT EXPECTING LB =',LB,', MATRIX IS SYMMETRIC.'
            CALL WARNING_MESSAGE(0)
            NWARNG = NWARNG + 1
         END IF
         GO TO 20
      END IF
      VRR = Y(N+1)
      IF(LB.EQ.2) VRR = VRR*VRR
      ENEXT = Y(N+2)
      IF(VRR.LT.0..OR.VRR.GT.1.E+04) THEN
         NP = N + 1
         EMESS = 'VARIANCE INCORRECT'
         CALL ERROR_MESSAGE(NSEQP3)
         WRITE (EMESS,'(4X,A,1PE12.5,A,I4)')                            &
     &         'VAR=',VRR,' AT LIST LOCATION',NP
         CALL ERROR_MESSAGE(0)
         NBAD = NBAD + 1
         GO TO 20
      END IF
      IF(VRR.GT.0.) CALL NEWGRD(ELAST,ENEXT)
      GO TO 20
   30 IF(LB.LE.4.OR.LB.GT.6.OR.IVR.EQ.0)    GO TO 100
!
!     CHECK VARIANCES AND STORE ENERGIES FOR LB=5.
!
      IF(LB.GE.6) THEN
         WRITE(EMESS,'(A,I2,A)')                                        &
     &         'NOT EXPECTING LB =',LB,', MATRIX IS SYMMETRIC.'
         CALL WARNING_MESSAGE(0)
         NWARNG = NWARNG + 1
         GO TO 100
      END IF
      LS = L1L
      NSTEP = N2L
      LOC = N - 1
      NEN = 1
   40 VRR = Y(N)
      ELAST = Y(NEN)
      ENEXT = Y(NEN+1)
      IF(VRR.LT.0..OR.VRR.GT.1.0E+04) THEN
         NSEQP2 = NSEQP1 + (N+5)/6
         EMESS = 'VARIANCE INCORRECT'
         CALL WARNING_MESSAGE(NSEQP2)
         WRITE (EMESS,'(4X,A,1PE12.5,A,I4)')                            &
     &         'VAR=',VRR,' AT LIST LOCATION',N
         CALL WARNING_MESSAGE(0)
         NBAD = NBAD + 1
         IF(NBAD.GT.100)    GO TO 100
      ELSE
         IF(VRR.GT.0.) CALL NEWGRD(ELAST,ENEXT)
      END IF
      IF(N.LT.NPL) THEN
         IF(LS.EQ.1)    NSTEP = NSTEP - 1
         N = N + NSTEP
         NEN = NEN + 1
         GO TO 40
      END IF
!
!     CHECK CORRELATION COEFFICIENTS FOR LB=5.
!
      DO I=1,NEN
         ILOC = N2L*I + 1 - LS*I*(I-1)/2
         VRI = Y(ILOC)
         JMIN = 1
         IF(LS.EQ.1)    JMIN = I
         DO J=JMIN,NEN
            LOC = LOC + 1
            IF(J.GT.I) THEN
               YLOC = Y(LOC)
               IF(YLOC.NE.0.) THEN
                  JLOC = N2L*J + 1 - LS*J*(J-1)/2
                  VRJ = Y(JLOC)
                  CORR = BIGNO
                  IF(VRI.GE.0..AND.VRJ.GE.0.) THEN
                     IF(VRI.EQ.0..OR.VRJ.EQ.0.) THEN
                        NSEQP2 = NSEQP1 + (LOC+5)/6
                        EMESS = 'CORRELATION COEFFICIENT INCORRECT'
                        CALL ERROR_MESSAGE(NSEQP2)
                        WRITE (EMESS,'(4X,A,1PE12.5,A,I4)')             &
     &                       'CORR=',CORR,' AT LIST LOCATION',LOC
                        CALL ERROR_MESSAGE(0)
                        NBAD = NBAD + 1
                        IF(NBAD.GT.100) GO TO 100
                     END IF
                     CORR = YLOC/SQRT(VRI*VRJ)
                     IF(ABS(CORR).GT.(1.+EPSILN4)) THEN
                        NSEQP2 = NSEQP1 + (LOC+5)/6
                        EMESS = 'CORRELATION COEFFICIENT INCORRECT'
                        CALL ERROR_MESSAGE(NSEQP2)
                        WRITE (EMESS,'(4X,A,1PE12.5,A,I4)')             &
     &                        'CORR=',CORR,' AT LIST LOCATION',LOC
                        CALL ERROR_MESSAGE(0)
                        NBAD = NBAD + 1
                        IF(NBAD.GT.100) GO TO 100
                     END IF
                  END IF
               END IF
            END IF
         END DO
      END DO
!
  100 RETURN
      END SUBROUTINE NITEST
!
!***********************************************************************
!
      SUBROUTINE NEWGRD(E1T,E2T)
!
!     INSERT AN ENERGY RANGE INTO ARRAYS EI AND EG.  ON RETURN,
!     EI CONTAINS ALL RANGES ENCOUNTERED THUS FAR IN THIS MT,
!     WITH ALL UNNECESSARY ENERGY POINTS STRIPPED OUT.
!     EG CONTAINS THE SAME INFORMATION FOR MT, PLUS THE FINAL SET
!     OF RANGES FOR ALL PREVIOUS MT-S IN THIS FILE.  THE ARRAY MTNI
!     CONTAINS THE LIST OF MT-S PROCESSED SO FAR, ALONG WITH THE
!     NUMBER OF ENERGY PAIRS STORED IN EG FOR EACH MT.
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: E1T,E2T
!
      INTEGER(KIND=I4) :: IOLD,IEF,K,ILOC
      INTEGER(KIND=I4) :: I,J
      REAL(KIND=R4) :: ELO,EHI
!
!     K = CURRENT READ POINT IN EG.
!     M = CURRENT LOAD POINT IN EI.
!
      IF(E2T.LE.E1T) GO TO 100
!
!     NEW MT.  ADD MT AND (E1T,E2T) TO THE ARRAYS.
!
      IF(NMT33.EQ.0.OR.MT.NE.MTNI(NMT33,1))  THEN
         IF(NMT33.GE.NMTMAX) THEN
            EMESS = ' '
            CALL ERROR_MESSAGE(0)
            EMESS = ' ***ERROR IN NEWGRD.  NMT EXCEEDS NMTMAX.'
            CALL ERROR_MESSAGE(0)
            NERROR = NERROR + 1
            IERX = 1
            GO TO 100
         END IF
         NMT33 = NMT33 + 1
         MTNI(NMT33,1) = MT
         MTNI(NMT33,2) = 1
         IF(NEG.GE.NEGMAX) THEN
            EMESS = ' '
            CALL ERROR_MESSAGE(0)
            EMESS = ' ***ERROR IN NEWGRD.  NEG EXCEEDS NEGMAX.'
            CALL ERROR_MESSAGE(0)
            NERROR = NERROR + 1
            IERX = 1
            GO TO 100
         END IF
         NEG = NEG + 1
         EGR33(NEG,1) = E1T
         EGR33(NEG,2) = E2T
         NIX = 1
         EI(1,1) = E1T
         EI(1,2) = E2T
         GO TO 100
      END IF
!
!     OLD MT.  INSERT E1T AND E2T INTO LAST SET OF RANGES, REMOVING
!     ANY INTERIOR POINTS WHEN OVERLAP OCCURS.
!
      IOLD = MTNI(NMT33,2)
      NIX = 0
      IEF = 0
      K = NEG - IOLD
      ELO = E1T
      EHI = E2T
!
!     INCREMENT THE ENERGY RANGE.
!
   20 IF(K.GE.NEG) THEN
!
!        STORE MODIFIED GRID BACK INTO EG.
!
         IF(IEF.EQ.1) THEN
            MTNI(NMT33,2) = NIX
            DO I=1,NIX
               ILOC = NEG - IOLD + I
               DO J=1,2
                  EGR33(ILOC,J) = EI(I,J)
               END DO
            END DO
            NEG = ILOC
            GO TO 100
         END IF
         GO TO 50
      END IF
      K = K + 1
      IF(EGR33(K,1).GT.EHI.AND.IEF.EQ.0) THEN
         K = K - 1
         GO TO 50
      END IF
      IF(EGR33(K,2).LT.ELO.OR.EGR33(K,1).GT.EHI) THEN
!
!        SAVE EG(K,1) AND EG(K,2)
!
         IF(NIX.LT.NIXMAX) THEN
            NIX = NIX + 1
            EI(NIX,1) = EGR33(K,1)
            EI(NIX,2) = EGR33(K,2)
            GO TO 20
         END IF
         GO TO 70
      END IF
!
!     THIS RANGE ADJOINS (ELO,EHI).  MODIFY ELO AND EHI ACCORDINGLY.
!
      IF(EGR33(K,1).LT.ELO.OR.EGR33(K,2).GT.EHI) THEN
         IF(EGR33(K,1).LT.ELO) ELO=EGR33(K,1)
         IF(EGR33(K,2).GT.EHI) EHI=EGR33(K,2)
      END IF
      GO TO 20
!
!     SAVE (ELO,EHI)
!
   50 IF(NIX.LT.NIXMAX) THEN
         NIX = NIX + 1
         EI(NIX,1) = ELO
         EI(NIX,2) = EHI
         IEF = 1
         GO TO 20
      END IF
!
!     FATAL ERROR
!
   70 EMESS = ' '
      CALL ERROR_MESSAGE(0)
      EMESS = ' ***ERROR IN NEWGRD.  NIX EXCEEDS NIXMAX.'
      CALL ERROR_MESSAGE(0)
      NERROR = NERROR + 1
      IERX = 1
!
  100 RETURN
      END SUBROUTINE NEWGRD
!
!***********************************************************************
!
      SUBROUTINE SUM452(N0)
!
!     PERFORMS SUMUP TEST ON MT = 452
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N0
!
      INTEGER(KIND=I4), SAVE :: ICO452
      DATA ICO452/0/
!
!     TEST INITIALIZATION
!
      IF(N0.LT.0) THEN
         NMTO = 0
         ICO452 = 0
         IPC = 0
         ITFLE = 452
         ITEST = 1
         GO TO 100
      END IF
!
!     PERFORM SUMUP TEST
!
      IF(N0.EQ.0)   THEN
         IF(NMTO.NE.3) THEN
            EMESS = ' '
            CALL ERROR_MESSAGE(0)
            EMESS = 'CANNOT DO NU BAR SUMUP TEST DUE TO '//             &
     &                   'MISSING PARTIALS'
            CALL ERROR_MESSAGE(0)
!********ALL COMPONENTS IN COEFFICIENT FORM
         ELSE IF(ICO452.NE.0)  THEN
            CALL COEDIF
         ELSE
!********COMPONENTS TESTED IN POINTWISE FORM
            CALL TSTTOT
         END IF
         GO TO 100
      END IF
!
!     MT = 452
!
      IF(MT.EQ.452)   THEN
!********TABULATED RECORDS
         IF(L2H.EQ.2)   THEN
            CALL STOR(NP,2)
!********COEFFICIENT FORM
         ELSE
            ICO452 = 1
            CALL STOCO
         END IF
         GO TO 50
!
!     MT=455
!
      ELSE IF(MT.EQ.455)   THEN
         IF(NMTO.EQ.1)   GO TO 20
!
!     MT=456
!
      ELSE IF(MT.EQ.456)   THEN
         IF(NMTO.EQ.2)   THEN
            ITEST = 1
            GO TO 20
         END IF
      END IF
      GO TO 100
!
!     COEFFICIENTS USED FOR EITHER 455 OR  456
!     USE  SUBROUTINE STOCO TO PUT COEFFICIENTS INTO CO ARRAY
!     CONVERT TO TABULAR FORM IF ICO452.EQ.0  (IF 452 WAS TABULAR)
!
   20 IF(L2H.EQ.1)   THEN
         CALL STOCO
         IF(ICO452.EQ.0) CALL CONV(MT)
      ELSE
!
!     TABULATED DATA
!
!********452 WAS COEFF
         IF(ICO452.GT.0)  THEN
            CALL STOR(NP,1)
            ICO452 = 0
            IF(MT.EQ.456)  CALL CONV(455)
         END IF
         CALL INTPX(XT,YT,NTOT,0)
      END IF
!
!     STORE MT NUMBER JUST PROCESSED
!
   50 NMTO = NMTO + 1
      MTOO(NMTO) = MT
!
  100 RETURN
      END SUBROUTINE SUM452
!
!***********************************************************************
!
      SUBROUTINE STOCO
!
!     STORES COEFFICIENTS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NUPPER,IAR
      INTEGER(KIND=I4) :: I
!
!     CHECK NUMBER OF COEFFICIENTS.  ONLY 4 ALLOWED.
!
      NUPPER = NPL
      IF(NPL.GT.4)  THEN
         NUPPER = 4
         WRITE(EMESS,'(A,I4,A,I2,A)')                                   &
     &            'MT=',MT,' HAS TOO MANY COEFFICIENTS(',NPL,')'
         CALL ERROR_MESSAGE(0)
      END IF
!
!     STORE THE COEFFICIENTS
!
      IF(MT.EQ.452)   THEN
         IAR = 1
      ELSE
         IAR = MT - 453
      END IF
      DO I=1,NUPPER
         COEFS(I,IAR) = YP(I)
      END DO
!
      RETURN
      END SUBROUTINE STOCO
!
!***********************************************************************
!
      SUBROUTINE CONV(IT)
!
!     ROUTINE TO CONVERT COEFFICIENTS TO POINTWISE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IT
!
      INTEGER(KIND=I4) :: K
      INTEGER(KIND=I4) :: J
      REAL(KIND=R4) :: XM,YPVAL
!
!     SELECT THE SET
!
      K = IT - 453
!
!     LOOP THRU THE POINTS
!
      DO J=1,NTOT
         XM = XT(J)
!********CONVERT TO THIRD ORDER POLYNOMIAL
         YPVAL = COEFS(1,K) +                                           &
     &               XM*(COEFS(2,K) + XM*(COEFS(3,K)+COEFS(4,K)*XM))
!********SUBTRACT FROM TOTAL
         YT(J) = YT(J) - YPVAL
      END DO
!
      RETURN
      END SUBROUTINE CONV
!
!***********************************************************************
!
      SUBROUTINE COEDIF
!
!     CALCULATES THE DIFFERENCE BETWEEN TOTAL AND PARTIALS WHEN ALL ARE
!      IN COEFFICIENT FORM AND COMPARES AGAINST EPSILN
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: ITLE
      INTEGER(KIND=I4) :: I,J
      REAL(KIND=R4) :: COT,CPART,DELTA
!
!     CHECK EACH COEFFICIENT
!
      ITLE = 0
      DO I=1,4
         COT = COEFS(I,1)
         CPART = COEFS(I,2) + COEFS(I,3)
         IF(COT.NE.0.0) THEN
            DELTA = (COT-CPART)/COT
         ELSE
            DELTA = COT
         END IF
         IF(ABS(DELTA).GE.FIZCON_DATA%EPSILN)  THEN
!***********HEADING FOR FIRST TIME ERROR DETECTED
            IF(ITLE.EQ.0) THEN
               ITLE=1
               EMESS = ' '
               CALL ERROR_MESSAGE(0)
               EMESS = 'THE FOLLOWING MTS WERE USED IN THE '//          &
     &                      'SUMUP TEST - '
               CALL ERROR_MESSAGE(0)
               WRITE(NOUT,'(15(I3,1X))')  (MTOO(J),J=1,NMTO)
               WRITE(EMESS,'(3X,A,F10.6)')                              &
     &            'SUM OF PARTIAL COEFFICIENTS DIFFERED FROM THE'//     &
     &         ' TOTAL BY MORE THAN ',FIZCON_DATA%EPSILN
               CALL ERROR_MESSAGE(0)
               EMESS = '         COEF #            TOTAL'//             &
     &                 '           SUM PARTIALS      DELTA'
               CALL ERROR_MESSAGE(0)
               EMESS = ' '
               CALL ERROR_MESSAGE(0)
            END IF
!********ERROR MESSAGE FOR EACH BAD COEFFICIENT SUM
            WRITE(NOUT,'(14X,I4,2X,2(5X,1PE15.5),0PF12.6)')             &
     &              I,COT,CPART,DELTA
         END IF
      END DO
!
      RETURN
      END SUBROUTINE COEDIF
!
!***********************************************************************
!
      SUBROUTINE SUMF3(N0)
!
!     SUMS UP PARTIAL CROSS SECTIONS IN FILE THREE AND
!       COMPARES THEM TO THE CORRESPONDING SUM CROSS SECTION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N0
!
      CHARACTER(LEN=*), INTRINSIC :: TRIM
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: IACTION,IDX,IMTT,ITS,MTT,MTL,IADV,NREC,NTM,K1
      INTEGER(KIND=I4) :: NBEG,NLMOD,NCONT
      INTEGER(KIND=I4) :: I,J,JJ,K,L
      REAL(KIND=R4) :: E1T,E2T,EHI,XINT
!
      INTEGER(KIND=I4), PARAMETER :: NTESTS=8
      INTEGER(KIND=I4), SAVE, DIMENSION(NTESTS) :: MTFLGS,MTTOT,MSKIP
      INTEGER(KIND=I4), DIMENSION(NTESTS), PARAMETER ::                  &
     &             IWHO = (/4,18,101,103,104,105,106,107/)
!
      INTEGER(KIND=I4), SAVE  :: IMT1,IMT2,NWR
!
!     INITIALIZE
!
      IF(N0.LT.0)   THEN
         REWIND (UNIT=ISCRU1)
         NWR = 0
         ITEST = 1
         IMT1 = 0
         IMT2 = 0
         MSKIP = 0
         MTTOT = 0
         MTFLGS = 0
         ITFLE = 1
         IPC = 0
         GO TO 100
      ELSE IF(N0.EQ.0) THEN
         GO TO 40
      END IF
!
!     PROCESS MT=1
!
      IACTION = 0
      IF(MT.EQ.1)   THEN
         IF((NSUB/10).EQ.1)    THEN
            IMT1 = 1
            NMTO = 1
            MTOO(1) = MT
            CALL STOR(NP,2)
         END IF
!
!     PROCESS MT=2
!
      ELSE IF(MT.EQ.2)  THEN
         IMT2 = 1
         IF(IMT1.NE.0)   IACTION = 2
!
!     MT=3
!
      ELSE IF(MT.EQ.3)  THEN
         IF(IMT2.EQ.0) THEN
            EMESS = '    MT=2 IS MISSING'
            CALL ERROR_MESSAGE(0)
!
!        INTERPOLATE AND SUBTRACT FROM RUNNING DIFFERENCE
!
         ELSE
            IF(IMT1.NE.0)   THEN
               NMTO = NMTO + 1
               MTOO(NMTO) = 3
               CALL INTPX(XT,YT,NTOT,0)
!
!              DO SUMUP TEST AS ALL PARTIALS HAVE BEEN PROCESSED
!                  FOR TEST MF(1) = MF(2) + MF(3)
!
               CALL TSTTOT
            END IF
         END IF
         ITFLE = 3
         IMT1 = 1
         NMTO = 1
         MTOO(1) = MT
         CALL STOR(NP,2)
!
!     MT=4
!
      ELSE IF(MT.EQ.4)  THEN
         MTFLGS(1) = 1
         NWR = NWR + 1
         MTTOT(1) = NWR
         CALL RDWRIT(ISCRU1,2)
         IACTION = 2
!
!     MT=5 THRU 17 EXCEPT MT=10
!
      ELSE IF(MT.LE.17) THEN
         IF(MT.NE.10)  IACTION = 2
!
!     MT=18
!
      ELSE IF(MT.EQ.18)  THEN
         MTFLGS(2) = 1
         NWR = NWR + 1
         MTTOT(2) = NWR
         CALL RDWRIT(ISCRU1,2)
         IACTION = 2
!
!     MT=19 THRU 21 AND 38
!
      ELSE IF((MT.GE.19.AND.MT.LE.21).OR.MT.EQ.38)  THEN
         IDX = 2
         IACTION = 1
!
!     REST OF MTS BELOW 42
!
      ELSE IF(MT.LT.45)   THEN
         IF(MT.NE.26.AND.MT.NE.27)  IACTION = 2
!
!     SINGLE NEUTRON CHANNELS
!
      ELSE IF(MT.LE.91) THEN
         IF(MT.GE.50)   THEN
            IDX = 1
            IACTION = 1
         END IF
!
!     UNASSIGNED MT'S 92-100
!
      ELSE IF(MT.LE.100)   THEN
!
!     MT=101
!
      ELSE IF(MT.EQ.101)   THEN
         MTFLGS(3) = 1
         NWR = NWR + 1
         MTTOT(3) = NWR
         CALL RDWRIT(ISCRU1,2)
         IACTION = 2
!
!     MT=102 THRU 107
!
      ELSE IF(MT.LE.107) THEN
         IDX = 3
         IF(MT.EQ.102)   THEN
            IACTION = 1
         ELSE
            IACTION = 2
            ITS = MT - 99
            MTFLGS(ITS) = 1
            NWR = NWR + 1
            MTTOT(ITS) = NWR
            CALL RDWRIT(ISCRU1,2)
            IMTT = MTFLGS(IDX)
            IF(IMTT.NE.0)   THEN
               IF(IMTT.EQ.1)   MSKIP(IDX) = NWR - MTTOT(IDX) - 1
               MTFLGS(IDX) = IMTT + 1
               IACTION = 0
            END IF
         END IF
!
!     MT = 108 THRU 117
!
      ELSE IF(MT.LE.117)   THEN
         IDX = 3
         IACTION = 1
!
!     MT = 600 THRU 849
!
      ELSE
         IF(NFOR.GE.6)   THEN
            IF(MT.LT.600.OR.MT.GT.849)    GO TO 100
            NBEG = 600
            NLMOD = 50
            NCONT = 49
         ELSE
            IF(MT.LT.699.OR.MT.GT.799)    GO TO 100
            NBEG = 700
            NLMOD = 20
            NCONT = 18
         END IF
         MTT = MT - NBEG
         MTL = MOD(MTT,NLMOD)
         IF(MTL.LE.NCONT)    THEN
            IDX = MTT/NLMOD + 4
            IACTION = 1
         END IF
      END IF
!
!     SAVE SECTION FOR LATER TEST
!
      IF(IACTION.EQ.1) THEN
         IMTT = MTFLGS(IDX)
         IF(IMTT.NE.0)   THEN
            MTFLGS(IDX) = IMTT + 1
            IF(IMTT.EQ.1)   MSKIP(IDX) = NWR - MTTOT(IDX)
            NWR = NWR + 1
            CALL RDWRIT(ISCRU1,2)
            GO TO 100
         END IF
         IACTION = 2
      END IF
!
!     INTERPOLATE AND SUBTRACT FORM RUNNING DIFFERENCE
!
      IF(IACTION.EQ.2) THEN
         IF(IMT1.EQ.1)   THEN
            CALL INTPX(XT,YT,NTOT,0)
            NMTO = NMTO + 1
            MTOO(NMTO) = MT
         END IF
      END IF
      GO TO 100
!
!     DO SUMUP TEST ON NONELASTIC (OR TOTAL IF NONELASTIC NOT PRESENT)
!
   40 IF(IMT1.EQ.1.AND.NMTO.GT.1)   CALL TSTTOT
!
!     DO THE REST OF FILE 3 TESTS
!
      DO JJ=1,NTESTS
         ITFLE = IWHO(JJ)
!
!        CAN TEST BE DONE?
!
         IF(MTFLGS(JJ).GT.1) THEN
!
!           LOCATE SUM IN SCRATCH FILE
!
            REWIND (UNIT=ISCRU1)
            IADV = MTTOT(JJ) - 1
            CALL ADVDSK(ISCRU1,IADV)
!
!           READ AND STORE MT AS MASTER
!
            CALL RDWRIT(ISCRU1,1)
            CALL STOR(NP,2)
            NMTO = 1
            MTOO(1) = MT
!
!           POSITION TO FIRST PARTIAL
!
            IADV = MSKIP(JJ)
            CALL ADVDSK(ISCRU1,IADV)
            NREC = MTFLGS(JJ) - 1
!
!           INITIALIZE PANEL INTEGRALS
!
            YINT = 0.0
!
!           PROCESS EACH PARTIAL
!
            DO I=1,NREC
               CALL RDWRIT(ISCRU1,1)
               CALL INTPX(XT,YT,NTOT,0)
               NMTO = NMTO + 1
               MTOO(NMTO) = MT
!
!              CALCULATE INTEGRAL OF PARTIAL CROSS SECTION OVER TOTAL
!              CROSS SECTION ENERGY POINTS.
!
               NTM = NTOT - 1
               K = 1
               DO J=1,NTM
                  E1T = XT(J)
                  EHI = XT(J+1)
                  IF(E1T.NE.EHI) THEN
   50                IF(X(K).LE.E1T .AND. K.LT.NTOT) THEN
                        K = K + 1
                        GO TO 50
                     END IF
   55                IF(K.NE.1) THEN
                        E2T = X(K)
                        IF(E2T.GT.EHI) E2T = EHI
                        DO L=1,NR
                           IF(K.LE.NBT(L)) GO TO 60
                        END DO
   60                   K1 = K - 1
                        CALL ECSI(X(K1),Y(K1),X(K),Y(K),E1T,E2T,        &
     &                                                  JNT(L),XINT)
                        YINT(J) = YINT(J) + XINT
                        IF(E2T.LT.EHI .AND. K.LT.NTOT) THEN
                           E1T = E2T
                           K = K + 1
                           GO TO 55
                        END IF
                     END IF
                  END IF
               END DO
            END DO
!
!           PERFORM SUMUP TEST
!
            CALL TSTTOT
!
!           PERFORM SUMUP TEST ON THE INTEGRAL
!
            REWIND (UNIT=ISCRU1)
            IADV = MTTOT(JJ) - 1
            CALL ADVDSK(ISCRU1,IADV)
            CALL RDWRIT(ISCRU1,1)
            IPC = 1
            CALL TSTTOT
            IPC = 0
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE SUMF3
!
!***********************************************************************
!
      SUBROUTINE SUMPAR(N0)
!
!     SUMS UP PARTIAL CROSS SECTIONS IN FILES 10, 12 AND 13. THEN
!       COMPARES THEM TO THE CORRESPONDING TOTAL CROSS SECTIONS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N0
!
      INTEGER(KIND=I4), SAVE :: ISTO,MATS,MFS,MTS,NSEQS,N00,ISCRU
!
!     INITIALIZE FOR SUMUP TEST
!
      IF(N0.LT.0) THEN
         IF(MF.GE.11)  THEN
            ITEST = 1
            ITFLE = MF
            GO TO 100
         ELSE
            IF(MF.EQ.9)  THEN
               ISTO = 3
            ELSE
               ISTO = 2
            END IF
            ISCRU = ISCRU2
         END IF
!
!     FIND THE FILE 3 SECTION IN THE SCRATCH FILE
!
         MATS = MAT
         MFS = MF
         MTS = MT
         NSEQS = NSEQ
         REWIND ISCRU
   10    CALL RDWRIT(ISCRU,1)
         IF(MT.NE.0)  THEN
            IF(MT.NE.MTS)  GO TO 10
            ITFLE = MFS
            ITEST = 1
            CALL STOR(NP,ISTO)
         ELSE
            ITEST = 0
            EMESS = 'CANNOT DO SUMUP TEST BECAUSE FILE 3 SECTION '//    &
     &          'MISSING'
            CALL ERROR_MESSAGE(0)
         END IF
         MAT = MATS
         MF = MFS
         MT = MTS
         NSEQ = NSEQS
!
!     PERFORM SUMUP TEST
!
      ELSE IF(N0.EQ.0) THEN
         CALL TSTTOT
         ITEST = 0
!
!     INTERPOLATE AND SUBTRACT
!
      ELSE
         N00 = N0
         CALL INTPX(XT,YT,NTOT,N00)
      END IF
!
  100 RETURN
      END SUBROUTINE SUMPAR
!
!***********************************************************************
!
      SUBROUTINE SUMGAM(N0)
!
!     SUMS UP PARTIAL CROSS SECTIONS IN FILE 23 AND
!       COMPARES THEM TO THE CORRESPONDING TOTAL CROSS SECTION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N0
!
      INTEGER(KIND=I4), SAVE :: IMT501,IMT516,I515
      INTEGER(KIND=I4) :: I
!
!     INITIALIZE FOR TEST
!
      IF(N0.LT.0)   THEN
         REWIND (UNIT=ISCRU1)
         IMT501 = 0
         IMT516 = 0
         I515 = 0
         ITFLE = 0
         IPC = 0
         NMTO = 0
         GO TO 1000
!
!     PROCESS ALL SECTIONS
!
      ELSE IF(N0.GT.0)   THEN
!
!     INITIALIZE WHEN TOTAL IS ENCOUNTERED
!
         IF(MT.EQ.501)   THEN
            ITEST = 1
            ITFLE = 501
            IMT501 = 1
            CALL STOR(NP,2)
            NMTO = NMTO + 1
            MTOO(NMTO) = 501
            GO TO 1000
         END IF
!
!     PAIR PRODUCTION
!
         IF(MT.GE.515.AND.MT.LE.517) THEN
            IF(MT.EQ.515.OR.MT.EQ.517)  I515 = I515 + 1
            IMT516 = IMT516 + 1
            ITEST = 1
!*****STORE ON DISK FOR FUTURE
            CALL RDWRIT(ISCRU1,2)
            IF(MT.EQ.516.AND.I515.EQ.1)   GO TO 1000
            IF(MT.EQ.517.AND.I515.NE.2)   GO TO 1000
         ELSE
            GO TO 1000
         END IF
!
!     INTERPOLATE AND SUBTRACT
!
         IF(IMT501.NE.0)   THEN
            CALL INTPX(XT,YT,NTOT,0)
            NMTO = NMTO + 1
            MTOO(NMTO) = MT
         END IF
         GO TO 1000
      END IF
!
!     DO TEST ON CURRENT 23 CONTENTS
!
      IF(IMT501.NE.0)   CALL TSTTOT
!
!     DO OTHER TEST
!
      REWIND (UNIT=ISCRU1)
      ITFLE = 516
!
!     SAVE TOTAL
!
      IF(IMT516.EQ.3)   THEN
         CALL ADVDSK(ISCRU1,1)
         CALL RDWRIT(ISCRU1,1)
         CALL STOR(NP,2)
         NMTO = 1
         MTOO(1) = MT
         REWIND (UNIT=ISCRU1)
         DO I=1,2
            CALL RDWRIT(ISCRU1,1)
            CALL INTPX(XT,YT,NTOT,0)
            NMTO = NMTO + 1
            MTOO(NMTO) = MT
            IF(I.EQ.1)  CALL ADVDSK(ISCRU1,1)
         END DO
         CALL TSTTOT
      ELSE
         IF(I515.NE.2)  THEN
             WRITE(EMESS,'(A,I3,A)')                                    &
     &              'MT=',ITFLE,' PARTIALS ARE NOT PRESENT --NO TEST'
             CALL ERROR_MESSAGE(0)
         END IF
      END IF
!
 1000 RETURN
      END SUBROUTINE SUMGAM
!
!***********************************************************************
!
      SUBROUTINE STOR(NPP,IT)
!
!     STORE A TOTAL FOR SUMUP TEST
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NPP,IT
!
      INTEGER(KIND=I4) :: I,IPG,IL
      INTEGER(KIND=I4) :: J,K
      REAL(KIND=R4) :: XM,YM
!
!     INITIALIZE
!
      I = 1
      IPG = 1
      NTOT = NPP
      REWIND (UNIT=ISCRXY)
!
!     PROCESS EACH POINT
!
      DO J=1,NPP
         XT(J)=X(J)
!********CREATE THIRD DEGREE POLYNOMIAL FROM COEFFICIENTS
         IF(IT.EQ.1)      THEN
            XM = XT(J)
            YM = COEFS(1,1) +                                           &
     &                  XM*(COEFS(2,1)+XM*(COEFS(3,1)+COEFS(4,1)*XM))
!********STRAIGHT TRANSFER
         ELSE IF(IT.EQ.2)   THEN
            YM = Y(J)
!*******FOR FILE 9, FILL TOTAL WITH 1'S
         ELSE IF(IT.EQ.3)   THEN
            YM = 1.0
         ELSE
            GO TO 100
         END IF
!
!        FILL OUT TOTAL
!
         YT(J) = YM
         YTOT(I) = YM
         I = I + 1
!
!        PAGE Y TOTAL OUT
!
         IF(I.GT.PAGESZ)  THEN
            IL = PAGESZ
            WRITE(ISCRXY)   IL,(YTOT(K),K=1,IL)
            DO K=1,6
               YTOT(K) = YTOT(K+PAGESZ-6)
            END DO
            I = 7
            IPG = IPG + 1
         END IF
      END DO
!
!     ALL PTS TRANSFERRED
!
      ILOWXY = (PAGESZ-6)*(IPG-1)
      IHIGHXY = ILOWXY + PAGESZ
      IPAGEXY = IPG
      IF(IPG.GT.1)   WRITE(ISCRXY)  I,(YTOT(K),K=1,I)
!
  100 RETURN
      END SUBROUTINE STOR
!
!***********************************************************************
!
      SUBROUTINE INTPX(AX,BY,NTOTS,N0)
!
!     INTERPOLATES X/Y ARRAY TO AX GRID AND SUBTRACTS FROM BY
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N0
      INTEGER(KIND=I4) :: NTOTS
      REAL(KIND=R4), DIMENSION(NTOTS) :: AX,BY
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     CHARACTER(LEN=*), INTRINSIC :: TRIM
!
      INTEGER(KIND=I4) :: KP,NST,NUS,NPART
      INTEGER(KIND=I4) :: NBEG,NPR
      INTEGER(KIND=I4) :: I,N,IK
      REAL(KIND=R4) :: XONE,XN,XA,YA
!
      INTEGER(KIND=I4), PARAMETER :: NPLIM=500
      REAL(KIND=R4), DIMENSION(NPLIM) :: XPART
!
!     INITIALIZE
!
      KP = 1
      NST = 1
      NUS = 0
      NPART = 0
!
!     FIND FIRST POINT ON AX GRID AT OR ABOVE FIRST POINT ON X GRID
!
      XONE = X(1)
      DO N=1,NTOTS
         IF(AX(N).GE.XONE)   THEN
            NBEG = N
            GO TO 10
         END IF
      END DO
      GO TO 100
!
!     PROCESS ALL REMAINING POINTS ON AX GRID
!
   10 DO N=NBEG,NTOTS
         XA = AX(N)
!
!     MOVE TO NEXT POINT ON X GRID
!
   30    XN = X(NST)
         IF(XA.GT.XN)  THEN
            IF(NST.NE.NUS)    THEN
               NPART = NPART + 1
               IF(NPART.LE.NPLIM)   XPART(NPART) = XN
            END IF
            NST = NST + 1
            IF(NST.LE.NP)   GO TO 30
            GO TO 50
!
!     GRID MATCH FOUND
!
         ELSE IF(XA.EQ.XN) THEN
            NUS = NST
            YA = Y(NST)
!***********HANDLE DOUBLE VALUE PROBLEM
            IF(N.LT.NTOTS.AND.NST.LT.NP)    THEN
               IF(XA.EQ.AX(N+1))   THEN
                  IF(XN.EQ.X(NST+1))   NST = NST + 1
               ELSE IF(XN.EQ.X(NST+1))   THEN
                  YA = 0.5*(YA+Y(NST+1))
                  NUS = NST + 1
                  NST = NST + 2
               END IF
            END IF
!
!        POINT NOT IN PARTIAL SO INTERPOLATE
!
         ELSE
   40       IF(NST.GT.NBT(KP))    THEN
               KP = KP + 1
               IF(KP.LT.NR)   GO TO 40
            END IF
            CALL TERP1(X(NST-1),Y(NST-1),XN,Y(NST),XA,YA,JNT(KP))
         END IF
!
!        DO SUBTRACTION
!
         BY(N) = BY(N) - YA
         IF(NST.GT.NP)   GO TO 50
      END DO
!
!     STORE ANY UNUSED POINTS ON THE X GRID
!
      IF(NUS.LT.NP)    THEN
         DO IK=NUS+1,NP
            NPART = NPART + 1
            IF(NPART.LE.NPLIM)  XPART(NPART) = X(IK)
         END DO
      END IF
!
!     WRITE OUT POINTS MISSING IN SUM GRID
!
   50 IF(NPART.GT.0)   THEN
         IF(N0.GT.0)   THEN
            WRITE(EMESS,'(4X,A,I3,A,I3,A)')                             &
     &           'MT= ',MT,' SUBSECTION ',N0,' BEING PROCESSED'
            WRITE(NOUT,'(/5X,A)')  TRIM(EMESS)
         ELSE
            WRITE(EMESS,'(4X,A,I3,A)')                                  &
     &           'MT= ',MT,' BEING PROCESSED'
            WRITE(NOUT,'(/5X,A)')  TRIM(EMESS)
         END IF
         WRITE(NOUT,60)
   60    FORMAT(14X,'POINTS IN PARTIAL, NOT IN TOTAL'/                  &
     &          14X,2(7X,'X',13X,'X',4X),7X,'X')
         NPR = MIN0(NPART,NPLIM)
         WRITE(NOUT,'((8X,5(1PE14.6)))')  (XPART(I),I=1,NPR)
         IF(NPR.LT.NPART) THEN
            EMESS = '         AND MORE'
            WRITE(NOUT,'(5X,A)')  TRIM(EMESS)
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE INTPX
!
!***********************************************************************
!
      SUBROUTINE TERP1(XA,YA,XB,YB,XI,YI,I)
!
!     INTERPOLATE ONE POINT=============================================
!     (XA,YA) AND(XB,YB) ARE THE END POINTS OF THE LINE
!     I IS THE INTERPOLATION CODE
!     (XI,YI) IS THE INTERPOLATED POINT
!     NOTE- IF A NEGATIVE OR ZERO ARGUMENT OF A LOG IS
!           DETECTED, THE INTERPOLATION CODE IS AUTOMATICALLY
!           CHANGED FROM LOG TO LINEAR
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: EXP, ALOG
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: XA,YA,XB,YB,XI,YI
!
!*****HISTOGRAM
      IF(I.EQ.1) THEN
         YI = YA
!*****LINEAR-LINEAR
      ELSE IF(I.EQ.2) THEN
         YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
!*****LOG-LINEAR
      ELSE IF(I.EQ.3) THEN
         IF(XA.LE.0..OR.XB.LE.0.) THEN
            YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
         ELSE
            YI = YA + ALOG(XI/XA)*(YB-YA)/ALOG(XB/XA)
         END IF
!*****LINEAR-LOG
      ELSE IF(I.EQ.4) THEN
         IF(YA.LE.0..OR.YB.LE.0.)   THEN
            YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
         ELSE
            YI = YA*EXP((XI-XA)*ALOG(YB/YA)/(XB-XA))
         END IF
!*****LOG-LOG
      ELSE IF(I.EQ.5) THEN
         IF(YA.LE.0..OR.YB.LE.0.)   THEN
            IF(XA.LE.0..OR.XB.LE.0.) THEN
               YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
            ELSE
               YI = YA + ALOG(XI/XA)*(YB-YA)/ALOG(XB/XA)
            END IF
         END IF
         IF(XA.LE.0..OR.XB.LE.0.)   THEN
            IF(YA.LE.0..OR.YB.LE.0.)   THEN
               YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
            ELSE
               YI = YA*EXP((XI-XA)*ALOG(YB/YA)/(XB-XA))
            END IF
         END IF
         IF(XI.LE.0.) THEN
            YI = YA + ALOG(XI/XA)*(YB-YA)/ALOG(XB/XA)
         ELSE
            YI = YA*EXP(ALOG(XI/XA)*ALOG(YB/YA)/ALOG(XB/XA))
         END IF
      END IF
!
      RETURN
      END SUBROUTINE TERP1
!
!***********************************************************************
!
      SUBROUTINE TSTTOT
!
!     ROUTINE TO TEST TOTAL AGAINST SUM OF PARTS
!
      IMPLICIT NONE
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     CHARACTER(LEN=*), INTRINSIC :: TRIM
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: INIT,J0,J1,LL
      INTEGER(KIND=I4) :: J,L,II,JJ
      REAL(KIND=R4) :: E1T,E2T,ETEST,XINT,DIFF,DELTA,YTOJ,YSUM
!
      INTEGER(KIND=I4), PARAMETER :: NTESTS=17
      CHARACTER(LEN=28), DIMENSION(NTESTS) ::  ITLE
      DATA ITLE                                                         &
     &         /'NU BAR                      ',                         &
     &          'TOTAL                       ',                         &
     &          'NONELASTIC                  ',                         &
     &          'SINGLE NEUTRON PRODUCTION   ',                         &
     &          'TOTAL                       ',                         &
     &          'TOTAL                       ',                         &
     &          'PHOTON PROD MULTIPLICITY    ',                         &
     &          'PHOTON PROD CROSS SECTION   ',                         &
     &          'FISSION CROSS SECTION       ',                         &
     &          'NEUTRON DISAPPEARENCE       ',                         &
     &          'SINGLE PROTON PRODUCTION    ',                         &
     &          'SINGLE DEUTERON PRODUCTION  ',                         &
     &          'SINGLE TRITON PRODUCTION    ',                         &
     &          'SINGLE HE-3 PRODUCTION      ',                         &
     &          'SINGLE ALPHA PRODUCTION     ',                         &
     &          'TOTAL PHOTON INTERACTION    ',                         &
     &          'PAIR PROD CROSS SECTIONS    '/
      INTEGER(KIND=I4), DIMENSION(NTESTS), PARAMETER ::                 &
     &         MTALL = (/452,1,3,4,9,10,12,13,18,101,103,104,105,106,   &
     &                   107,501,516/)
      INTEGER(KIND=I4), DIMENSION(NTESTS), PARAMETER ::                 &
     &         MFALL = (/1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,23,23/)
!
!     INITIALIZE
!
      INIT = 0
!
!     LOOP THRU THE POINTS
!
      DO 100 J=1,NTOT
      J0 = J
!
!     CALCULATE PANEL INTEGRAL FOR SUM
!
      IF(IPC.GE.1) THEN
         IF(J.NE.NTOT) THEN
            J1 = J + 1
            E1T = XT(J)
            E2T = XT(J1)
            ETEST = (E2T-E1T)/E1T
            IF(ETEST.GT..00001) THEN
               DO L=1,NR
                  LL = L
                  IF(J1.LE.NBT(L)) GO TO 10
               END DO
   10          CALL ECSI(E1T,YTO(J0),E2T,YTO(J1),E1T,E2T,JNT(LL),XINT)
               DIFF = XINT - YINT(J)
               DELTA = ABS(DIFF)
               IF(XINT.NE.0.0) DELTA = DELTA/XINT
               GO TO 15
            ELSE
               DELTA = 0.
            END IF
         END IF
         GO TO 100
!
!     NORMAL DIFFERENTIAL TEST
!
      ELSE
         YTOJ = YTO(J0)
         IF(YTOJ.NE.0.0)   THEN
            DELTA = YT(J)/YTOJ
         ELSE
            DELTA = YT(J)
         END IF
      END IF
!
!     MULTIPLICITIES SUM TO UNITY OR LESS IN FILE 9
!     CROSS SECTIONS SUM TO FILE 3 OR LESS IN FILE 10
!
   15 IF(ITFLE.EQ.9.OR.ITFLE.EQ.10)   THEN
         IF(DELTA.LT.-FIZCON_DATA%EPSILN)  THEN
            IF(INIT.EQ.0)    THEN
               INIT = 1
               IF(MF.EQ.9)  THEN
                  WRITE(NOUT,20)  FIZCON_DATA%EPSILN
   20             FORMAT(/5X,'SUM OF MULTIPLICITIES EXCEEDED UNITY ',   &
     &               'BY MORE THAN ',F10.6,' AT THE FOLLOWING POINTS'// &
     &                7X,'ENERGY             SUMMATION')
               ELSE
                  WRITE(NOUT,25)  FIZCON_DATA%EPSILN
   25             FORMAT(/5X,'SUM OF CROSS SECTIONS EXCEEDED FILE 3 ',  &
     &               'BY MORE THAN ',F10.6,' AT THE FOLLOWING POINTS'// &
     &               7X,'ENERGY             SUMMATION')
               END IF
            END IF
!
!        WRITE OUT POINT IN QUESTION
!
            YSUM = YTOJ*(1. - DELTA)
            WRITE(NOUT,'(1PE15.5,5X,1PE15.5)')  XT(J),YSUM
         END IF
!
!     COMPARE AGAINST PERCENTAGE ERROR
!
      ELSE
         IF(ABS(DELTA).LE.FIZCON_DATA%EPSILN) GO TO 100
!
!        TEST FAILED
!
         IF(INIT.EQ.0) THEN
!
!           OUTPUT TITLE ON FIRST FAILURE OF TEST
!
            INIT = 1
            DO II=1,NTESTS
               IF(ITFLE.EQ.MTALL(II))   THEN
                  WRITE(EMESS,'(3A,I2,A,I3,A)')                         &
     &                 'SUMUP TEST FAILED FOR ',TRIM(ITLE(II)),
     &                 ' (MF=',MFALL(II),' AND MT=',MTALL(II),')'
                  WRITE(NOUT,'(//5X,A)')  TRIM(EMESS)
                  IF(ITFLE.NE.12.AND.ITFLE.NE.13)   THEN
                     WRITE(NOUT,40)  (MTOO(JJ),JJ=1,NMTO)
   40                FORMAT(/5X,'THE FOLLOWING MTS WERE USED IN ',      &
     &                   'THE SUMUP TEST - ',6(I3,1X)/9X,18(I3,1X))
                  END IF
                  GO TO 50
               END IF
            END DO
            EMESS = ' '
            CALL ERROR_MESSAGE(0)
            WRITE(EMESS,'(A,I4,A)')                                     &
     &        'TITLE FOR FILE (OR SECTION)',ITFLE,' NOT FOUND IN TABLE'
            CALL ERROR_MESSAGE(0)
            EMESS = '    TITLE WILL BE INCORRECT ON REPORT'
!
!           OUTPUT TEST TYPE DEPENDENT TABLE HEADINGS
!
   50       IF(IPC.EQ.1) THEN
               WRITE(NOUT,55) FIZCON_DATA%EPSILN
   55          FORMAT(/5X,'THE SUM OF PARTIAL INTEGRALS DIFFERED FROM ',&
     &                  'THE TOTAL INTEGRAL BY MORE '/                  &
     &           9X,'THAN ',F10.6,' IN THE FOLLOWING INTERVALS'/        &
     &           5X,'IF SUMUP POINT TEST DID NOT FAIL, ',               &
     &               'INTERPOLATION TYPES IN PARTIALS'/                 &
     &           9X,' ARE DIFFERENT FROM INTERPOLATION TYPE IN TOTAL'// &
     &           8X,'ENERGY RANGE          TOTAL INTEGRAL  ',           &
     &                  'TOTAL-PARTIALS      DELTA')
            ELSE
               WRITE(NOUT,60) FIZCON_DATA%EPSILN
   60          FORMAT(/5X,'THE SUM OF PARTIALS DIFFERED FROM THE ',     &
     &           'TOTAL BY MORE THAN ',F10.6/                           &
     &           10X,'AT THE FOLLOWING POINTS'//                        &
     &           '       ENERGY         TOTAL C-SECTION    ',           &
     &           '(TOTAL)-(PARTIALS)         DELTA')
            END IF
         END IF
!
!        WRITE OUT POINT IN QUESTION
!
         IF(IPC.EQ.1) THEN
            WRITE(NOUT,'(1X,1PE11.4,A,1PE11.4,2(1PE16.5),0PF14.6)')     &
     &                    E1T,' TO ',E2T,XINT,DIFF,DELTA
         ELSE
            IF(ABS(DELTA).GT.10.) THEN
               WRITE(NOUT,'(3(1PE15.6,5X),3X,A)')                       &
     &                   XT(J),YTOJ,YT(J),'GREATER THAN 10'
            ELSE
               WRITE(NOUT,'(3(1PE15.6,5X),5X,0PF10.6)')                 &
     &                   XT(J),YTOJ,YT(J),DELTA
            end if
         END IF
      END IF
!
  100 CONTINUE
!
      RETURN
      END SUBROUTINE TSTTOT
!
!***********************************************************************
!
      SUBROUTINE ECSI(X3,Y3,X4,Y4,X1,X2,I,ANS)
!
!     COMPUTE INTEGRAL OF Y(X)==========================================
!     Y(X) DEFINED BY THE END POINTS (X3,Y3), (X4,Y4), AND THE
!     INTERPOLATION CODE I.  X1 AND X2 ARE THE INTEGRATION LIMITS.
!     X1 AND X2 MAY LIE OUTSIDE X3 AND X4
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: X3,Y3,X4,Y4,X1,X2,ANS
!
      REAL(KIND=R4), INTRINSIC :: ABS,EXP,ALOG
!
      REAL(KIND=R4) :: A,B,Z
!
      ANS = 0.0
      IF(X4.LE.X3)   GO TO 100
!
!     Y CONSTANT
!
      IF(I.EQ.1) THEN
         ANS = (X2-X1)*Y3
!
!     Y LINEAR IN X
!
      ELSE IF(I.EQ.2) THEN
         B = (Y4-Y3)/(X4-X3)
         A = Y3-B*X3
         ANS=(X2-X1)*(A+0.5*B*(X2+X1))
!
!     Y LINEAR IN LN(X)
!
      ELSE IF(I.EQ.3) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)  THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3-B*X3
            ANS = (X2-X1)*(A+0.5*B*(X2+X1))
         ELSE
            B = (Y4-Y3)/ALOG(X4/X3)
            Z = (X2-X1)/X1
            IF(ABS(Z).LE.0.1)    THEN
               ANS = (X2-X1)*(Y3+B*ALOG(X1/X3))+(0.5*B*X1*Z*Z)*         &
     &                 (1.0+Z*(-0.33333333+Z*(0.16666667-0.1*Z)))
            ELSE
               ANS = (X2-X1)*(Y3+B*ALOG(X1/X3))+B*X1*                   &
     &                 (1.0+(X2/X1)*(ALOG(X2/X1)-1.0))
            END IF
         END IF
!
!     LN(Y) LINEAR IN X
!
      ELSE IF(I.EQ.4) THEN
         IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3-B*X3
            ANS = (X2-X1)*(A+0.5*B*(X2+X1))
         ELSE
            B = ALOG(Y4/Y3)/(X4-X3)
            A = ALOG(Y3)-B*X3
            Z = (X2-X1)*B
            IF(ABS(Z).LE.0.1)THEN
               ANS = EXP(A+B*X1)*(X2-X1)*(1.0+Z*(0.5+0.16666667*Z))
            ELSE
               ANS = EXP(A+B*X1)*(EXP(Z)-1.0)/B
            END IF
         END IF
!
!     LN(Y) LINEAR IN LN(X)
!
      ELSE IF(I.EQ.5) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
            IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
               B = (Y4-Y3)/(X4-X3)
               A = Y3-B*X3
               ANS = (X2-X1)*(A+0.5*B*(X2+X1))
            ELSE
               B = ALOG(Y4/Y3)/(X4-X3)
               A = ALOG(Y3)-B*X3
               Z = (X2-X1)*B
               IF(ABS(Z).LE.0.1)THEN
                  ANS = EXP(A+B*X1)*(X2-X1)*(1.0+Z*(0.5+0.16666667*Z))
               ELSE
                  ANS = EXP(A+B*X1)*(EXP(Z)-1.0)/B
               END IF
            END IF
         ELSE IF(Y3.LE.0.0.OR.Y4.LE.0.0)  THEN
            B = (Y4-Y3)/ALOG(X4/X3)
            Z = (X2-X1)/X1
            IF(ABS(Z).LE.0.1)    THEN
               ANS = (X2-X1)*(Y3+B*ALOG(X1/X3))+(0.5*B*X1*Z*Z)*         &
     &               (1.0+Z*(-0.33333333+Z*(0.16666667-0.1*Z)))
            ELSE
               ANS = (X2-X1)*(Y3+B*ALOG(X1/X3))+B*X1*                   &
     &               (1.0+(X2/X1)*(ALOG(X2/X1)-1.0))
            END IF
         ELSE
            B = ALOG(Y4/Y3)/ALOG(X4/X3)
            Z = (B+1.0)*ALOG(X2/X1)
            IF(ABS(Z).LE.0.1)   THEN
               ANS = Y3*X1*((X1/X3)**B)*ALOG(X2/X1)*                    &
     &                   (1.+Z*(0.5+0.16666667*Z))
            ELSE
               ANS = Y3*X1*((X1/X3)**B)*(((X2/X1)**(B+1.0))-1.0)/(B+1.0)
            END IF
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE ECSI
!
!***********************************************************************
!
      SUBROUTINE XECSI(X3,Y3,X4,Y4,X1,X2,I,ANS)
!
!     COMPUTE INTEGRAL OF X*Y(X)=================================
!     Y(X) DEFINED BY THE END POINTS (X3,Y3), (X4,Y4), AND THE
!     INTERPOLATION CODE I.  X1 AND X2 ARE THE INTEGRATION LIMITS.
!     X1 AND X2 MAY LIE OUTSIDE X3 AND X4
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: X3,Y3,X4,Y4,X1,X2,ANS
!
      REAL(KIND=R4), INTRINSIC :: ABS, ALOG, EXP
!
      REAL(KIND=R4) :: A,B,Z,B3,X12
!
      ANS = 0.0
      IF(X4.LE.X3)   GO TO 100
!
!     BRANCH ON INTERPOLATION SCHEME
!
!     Y CONSTANT
!
      IF(I.EQ.1) THEN
        ANS = (X2-X1)*Y3*(X2+X1)/2.
!
!     Y LINEAR IN X
!
      ELSE IF(I.EQ.2) THEN
         B= (Y4-Y3)/(X4-X3)
         A= Y3-B*X3
         B3 = B/3.
         X12 = X1 + X2
         ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
!
!     Y LINEAR IN LN(X)
!
      ELSE IF(I.EQ.3) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3-B*X3
            B3 = B/3.
            X12 = X1 + X2
            ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
         ELSE
            B = (Y4-Y3)/ALOG(X4/X3)
            A = Y3 - B*ALOG(X3)
            X12 = X2/X1
            Z = (X2-X1)/X1
            IF(ABS(Z).LE.0.1)    THEN
               ANS = (X2-X1)*(A+B*ALOG(X1))*(X2+X1)/2. +                &
     &              (0.5*B*(X2-X1)*(X2-X1)*(1.+Z*(.3333333-Z*           &
     &               (0.08333333-0.03333333*Z))))
            ELSE
               ANS = (X2-X1)*(A+B*ALOG(X1))*(X1+X2)/2.                  &
     &              +0.25*B*X1*X1*(1.+X12*X12*(2.0*ALOG(X12)-1.0))
            END IF
         END IF
!
!     LN(Y) LINEAR IN X
!
      ELSE IF(I.EQ.4) THEN
         IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3-B*X3
            B3 = B/3.
            X12 = X1 + X2
            ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
         ELSE
            B = ALOG(Y4/Y3)/(X4-X3)
            A = ALOG(Y3)-B*X3
            Z = (X2-X1)*B
            IF(ABS(Z).LE.0.1)   THEN
               ANS = EXP(A+B*X1)*(X2-X1)*(1.+(B*X2-1.)*(1.+Z*(.5+Z*     &
     &             (.3333333+0.25*Z))))/B
            ELSE
               ANS = EXP(A+B*X1)*((B*X2-1.)*EXP(Z)-(B*X1-1.))/(B*B)
            END IF
         END IF
!
!     LN(Y) LINEAR IN LN(X)
!
      ELSE IF(I.EQ.5) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
            IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
               B = (Y4-Y3)/(X4-X3)
               A = Y3-B*X3
               B3 = B/3.
               X12 = X1 + X2
               ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
            ELSE
               B = ALOG(Y4/Y3)/(X4-X3)
               A = ALOG(Y3)-B*X3
               Z = (X2-X1)*B
               IF(ABS(Z).LE.0.1)   THEN
                  ANS = EXP(A+B*X1)*(X2-X1)*(1.+(B*X2-1.)*(1.+Z*(.5+Z*  &
     &                (.3333333+0.25*Z))))/B
               ELSE
                  ANS = EXP(A+B*X1)*((B*X2-1.)*EXP(Z)-(B*X1-1.))/(B*B)
               END IF
            END IF
         ELSE IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
            B = (Y4-Y3)/ALOG(X4/X3)
            A = Y3 - B*ALOG(X3)
            X12 = X2/X1
            Z = (X2-X1)/X1
            IF(ABS(Z).LE.0.1)    THEN
               ANS = (X2-X1)*(A+B*ALOG(X1))*(X2+X1)/2. +                &
     &              (0.5*B*(X2-X1)*(X2-X1)*(1.+Z*(.3333333-Z*           &
     &               (0.08333333-0.03333333*Z))))
            ELSE
               ANS = (X2-X1)*(A+B*ALOG(X1))*(X1+X2)/2.                  &
     &              +0.25*B*X1*X1*(1.+X12*X12*(2.0*ALOG(X12)-1.0))
            END IF
         ELSE
            B = ALOG(Y4/Y3)/ALOG(X4/X3)
            Z = (B+2.)*ALOG(X2/X1)
            IF(ABS(Z).LE.0.1)   THEN
               ANS=Y3*X1*X1**((X1/X3)**B)*ALOG(X2/X1)*                  &
     &                        (1.+Z*(0.5+0.16666667*Z))
            ELSE
               ANS = Y3*X1*X1*((X1/X3)**B)*(((X2/X1)**(B+2.))-1.)/(B+2.)
            END IF
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE XECSI
!
!***********************************************************************
!
      SUBROUTINE EAVE(EBAR,FNORM)
!
!     CALCULATES THE AVERAGE ENERGY FOR A TABLUAR ENERGY DISTRIBUTION
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: EBAR,FNORM
!
      INTEGER(KIND=I4) :: NTERP,INTT
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: XY,YN,X3,X4,Y3,Y4,ANS
!
      EBAR = 0.0
      FNORM = 0.0
      NTERP = 1
      INTT = JNT(1)
      XY = 0.0
      YN = 0.0
      DO N=2,NP
         IF(N.GT.NBT(NTERP))   THEN
            IF(NTERP.LT.NR)    NTERP = NTERP + 1
            INTT = JNT(NTERP)
         END IF
         X3 = X(N-1)
         X4 = X(N)
         Y3 = Y(N-1)
         Y4 = Y(N)
         CALL XECSI(X3,Y3,X4,Y4,X3,X4,INTT,ANS)
         XY = XY + ANS
         CALL ECSI(X3,Y3,X4,Y4,X3,X4,INTT,ANS)
         YN = YN + ANS
      END DO
      IF(YN.GT.0.0)   THEN
         FNORM = YN
         EBAR = XY/FNORM
      END IF
!
      RETURN
      END SUBROUTINE EAVE
!
!***********************************************************************
!
      SUBROUTINE RDWRIT(ISCRU,IFL)
!
!     READ OR WRITE A TAB1 RECORD IN A SCRATCH FILE
!
!             IFL = 1: READ
!             IFL = 2: WRITE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: ISCRU,IFL
!
      INTEGER(KIND=I4) :: NI,NPP,INP,IFST,NF,IALL
      INTEGER(KIND=I4) :: I,J,K,II
      REAL(KIND=R4), DIMENSION(12) ::  X1,X2
      REAL(KIND=R8), DIMENSION(12) ::  Y1,Y2
!
      IF(IFL.EQ.2)   GO TO 50
!
!     READ A TAB1 SECTION
!
      MT = 0
      READ(ISCRU,END=100)   C1,C2,L1,L2,NR,NP,MAT,MF,MT,NSEQ
      READ(ISCRU)  (NBT(I),JNT(I),I=1,NR)
!
!     PAGE POINTS OUT
!
      NI = 1
      NPP = 0
      IPAGE = 1
      IFST = 0
   10 NF = NI + 11
!*****ADDED POINTS IN INPUT RECORD DO NOT EXCEED PAGING ARRAY
      IF(NF.LE.996)   THEN
         READ(ISCRU)  (XP(I),YP(I),I=NI,NF)
         NI = NI + 12
         NPP = NPP + 12
         IALL = 1
!*****ADDED POINTS IN INPUT RECORD DO EXCEED PAGING ARRAY
      ELSE
         READ(ISCRU)  (XP(I),YP(I),I=NI,996),(X2(II),Y2(II),II=1,12)
         NF = PAGESZ
         NPP = NPP + 12
         IALL = 0
      END IF
!*****NO MORE POINTS SO OUTPUT FINAL PAGE
      IF(NPP.GE.NP)   THEN
         INP = -1
         CALL OUTXY(INP,1,NF)
         GO TO 100
!*****OUTPUT A FULL PAGE
      ELSE IF(NF.EQ.PAGESZ)   THEN
         CALL OUTXY(IFST,1,NF)
         NI = 13
!********STORE ANY POINTS WHICH DID NOT FIT ON PRIOR PAGE
         IF(IALL.EQ.0) THEN
            DO II=1,12
               XP(II+12) = X2(II)
               YP(II+12) = Y2(II)
            END DO
            NI = NI + 12
            NPP = NPP + 12
         END IF
      END IF
      GO TO 10
!
!     WRITE A TAB1 SECTION
!
   50 WRITE(ISCRU)   C1,C2,L1,L2,NR,NP,MAT,MF,MT,NSEQ
      WRITE(ISCRU)     (NBT(I),JNT(I),I=1,NR)
!
!     OUTPUT POINTS 12 AT A TIME
!
      DO I=1,NP,12
         DO J=1,12
            K = I - 1 + J
            IF(K.LE.NP)  THEN
               X1(J) = X(K)
               Y1(J) = Y(K)
            ELSE
               X1(J) = 0.0
               Y1(J) = 0.0
            END IF
         END DO
         WRITE(ISCRU)(X1(J),Y1(J),J=1,12)
      END DO
!
  100 RETURN
      END SUBROUTINE RDWRIT
!
!***********************************************************************
!
      SUBROUTINE ADVDSK(ISCRU,N)
!
!     ADVANCES N TAB1 RECORDS ON ISCRU
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: ISCRU,N
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2,MAT,MF,MT,NS,NBT,JNT
      INTEGER(KIND=I4) :: K,II,NN
      REAL(KIND=R4) :: C1,C2,X1
      REAL(KIND=R8) :: Y1
!
      IF(N.GT.0) THEN
         DO K=1,N
            READ(ISCRU) C1,C2,L1,L2,N1,N2,MAT,MF,MT,NS
            READ(ISCRU) (NBT,JNT,II=1,N1)
            DO NN=1,N2,12
               READ(ISCRU) (X1,Y1,II=1,12)
            END DO
         END DO
      END IF
!
      RETURN
      END SUBROUTINE ADVDSK
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION YTO(INDEXX)
!
!     LOGICAL PAGING SYSTEM FOR Y TOTAL ARRAY USED IN SUMUP
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDEXX
!
!     POINT NOT ON CURRENT PAGE
!
      IF(INDEXX.LE.ILOWXY.OR.INDEXX.GT.IHIGHXY)  CALL LOADYT(INDEXX)
!
!     RETURN REQUIRED DATA VALUE
!
      YTO = YTOT(INDEXX-ILOWXY)
!
      RETURN
      END FUNCTION YTO
!
!***********************************************************************
!
      SUBROUTINE LOADYT(INDEXX)
!
!     PAGING LOADER FOR Y TOTAL ARRAY USED IN SUMUP
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDEXX
!
      INTEGER(KIND=I4) :: JPAGE,ITOP
      INTEGER(KIND=I4) :: K
!
!     DETERMINE WHICH PAGE TO LOAD
!
      JPAGE = (INDEXX-7)/(PAGESZ-6) + 1
!
!     IF CURRENT PAGE IS PAST REQUESTED PAGE REWIND SCRATCH
!
      IF(JPAGE.LE.IPAGEXY)   THEN
         REWIND (UNIT=ISCRXY)
         IPAGEXY = 0
      END IF
!
!     SKIP TO PAGE PRECEDING ONE REQUIRED
!
   10 IPAGEXY = IPAGEXY + 1
      IF(IPAGEXY.NE.JPAGE)   THEN
         READ(ISCRXY)  ITOP
         GO TO 10
      END IF
!
!     LOAD CORRECT PAGE
!
      READ(ISCRXY)   ITOP,(YTOT(K),K=1,ITOP)
!
!     LOWER AND UPPER INDICES
!
      ILOWXY = (PAGESZ-6)*(IPAGEXY-1)
      IHIGHXY = PAGESZ + ILOWXY
!
      RETURN
      END SUBROUTINE LOADYT
!
!***********************************************************************
!
      SUBROUTINE STORF(MF1,MT1,ELO,EHI)
!
!     ROUTINE TO SAVE ENERGY LIMITS OF A SECTION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MF1,MT1
      REAL(KIND=R4) :: ELO,EHI
!
      INTEGER(KIND=I4) :: MFT
      INTEGER(KIND=I4) :: N
!
!     FIND SECTION IN INDEX
!
      MFT = 1000*MF1 + MT1
      DO N=1,NXC
!********SECTION FOUND SO STORE LIMITS
         IF(INDX(N,1).EQ.MFT)   THEN
            ENGS(N,1) = ELO
            ENGS(N,2) = EHI
            GO TO 100
         ELSE IF(INDX(N,1).GT.MFT)   THEN
            GO TO 100
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE STORF
!
!***********************************************************************
!
      SUBROUTINE EFTEST(MF1,MF2)
!
!     ROUTINE TO PERFORM TESTS FOR MISSING SECTIONS IN A FILE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MF1,MF2
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     CHARACTER(LEN=*), INTRINSIC :: TRIM
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: MFL,MFU
      INTEGER(KIND=I4) :: NMISS,LMFSN,MFR,MTR
      INTEGER(KIND=I4) :: I455,I456
      INTEGER(KIND=I4) :: M,N
!
!     INITIALIZE RANGE OF FILES TO BE CHECKED
!
      MFL = MF1
      MFU = MIN0(MF2,15)
!
!     PROCESS EACH FILE IN THE RANGE REQUESTED
!
      DO 100 M=MFL,MFU
      NMISS = 0
!
!     IN FILE 1, BOTH 455 AND 456 MUST BE ABSENT OR PRESENT
!
      IF(M.EQ.1)  THEN
         IF(LFI.EQ.1)   THEN
            CALL TESTS(1000*M+455,I455)
            I455 = MOD(I455,2)
            CALL TESTS(1000*M+456,I456)
            I456 = MOD(I456,2)
            IF(I455.NE.I456) THEN
               EMESS = 'BOTH SECTIONS MF=1, MT=455 AND MT=456 '//       &
     &                        'MUST BE PRESENT'
               WRITE(NOUT,'(/5X,A)')  TRIM(EMESS)
               NMISS = 1
!              Increment count for warnings
               NWARNG=NWARNG+1
            END IF
         ELSE
            NMISS = -1
         END IF
!
!     ARE FILE 4 SECTIONS REQUIRED BY FILE 3 PRESENT?
!
      ELSE IF(M.EQ.4) THEN
         IF(IFL3.EQ.0.OR.NSUB.NE.10)   THEN
            NMISS = -1
         ELSE
            CALL MISFIL(M,2,2,NMISS)
            CALL MISFIL(M,5,17,NMISS)
            CALL MISFIL(M,22,26,NMISS)
            CALL MISFIL(M,28,37,NMISS)
            CALL MISFIL(M,39,42,NMISS)
            CALL MISFIL(M,50,91,NMISS)
          END IF
!
!     ARE FILE 5 SECTIONS REQUIRED BY FILE 3 PRESENT?
!
      ELSE IF(M.EQ.5) THEN
         IF(IFL3.EQ.0.OR.NSUB.NE.10)   THEN
            NMISS = -1
         ELSE
            CALL MISFIL(M,5,17,NMISS)
            CALL MISFIL(M,22,26,NMISS)
            CALL MISFIL(M,28,37,NMISS)
            CALL MISFIL(M,39,42,NMISS)
            CALL MISFIL(M,91,91,NMISS)
         END IF
!
!     ARE FILE 6 SECTIONS REQUIRED BY FILE 3 PRESENT?
!
      ELSE IF(M.EQ.6) THEN
         IF(IFL3.EQ.0.OR.NSUB.EQ.10)   THEN
            NMISS = -1
         ELSE
            CALL MISFIL(M,2,2,NMISS)
            CALL MISFIL(M,5,26,NMISS)
            CALL MISFIL(M,28,42,NMISS)
            CALL MISFIL(M,50,91,NMISS)
            CALL MISFIL(M,102,120,NMISS)
            CALL MISFIL(M,600,849,NMISS)
         END IF
!
!     ARE FILE 9 SECTIONS REQUIRED BY FILE 8 PRESENT?
!       AND ARE FILE 10 SECTIONS REQUIRED BY FILE 8 PRESENT?
!
      ELSE IF(M.EQ.9.OR.M.EQ.10)   THEN
         IF(NLMF.EQ.0)  THEN
            NMISS = -1
         ELSE
            DO N=1,NLMF
               LMFSN = LMFS(2,N)
               IF(LMFSN.NE.0.AND.LMFSN.EQ.M) THEN
                  NMISS = NMISS + 1
                  IF(MISFERR.NE.M) THEN
                     WRITE(NOUT,'(/2X,A,I5,A,I3)')                      &
     &               'ERROR(S) - MISSING SECTIONS IN MAT',MAT,' MF',M
                     MISFERR = M
!                    Increment count for warnings
                  END IF
                  WRITE(EMESS,'(A,I4,A,I3)')                            &
     &               'FILE 8, MT=',LMFS(1,N),                           &
     &               ' REQUIRES SAME MT IN FILE',LMFS(2,N)
                  WRITE(NOUT,'(5X,A)')  TRIM(EMESS)
!                 Increment count for warnings
                  NERROR=NERROR+1
               END IF
            END DO
            IF(M.EQ.10)  NLMF = 0
         END IF
!
!     ARE FILE 14 SECTIONS REQUIRED BY FILES 12 AND 13 PRESENT?
!
      ELSE IF(M.EQ.14) THEN
         IF(NPMT.EQ.0.OR.NSUB.NE.10)  THEN
            NMISS = -1
         ELSE
            DO N=1,NPMT
               IF(ICON(N,2).LT.2000) THEN
                  MFR = ICON(N,1)/1000
                  MTR = MOD(ICON(N,1),1000)
                  NMISS = NMISS + 1
                  IF(MISFERR.NE.M) THEN
                     WRITE(NOUT,'(/2X,A,I5,A,I3)')                      &
     &               'ERROR(S) - MISSING SECTIONS IN MAT',MAT,' MF',M
                     MISFERR = M
!                    Increment count for warnings
                  END IF
                  WRITE(EMESS,'(A,I3,A,I3,A,I4)')                       &
     &                'CONTENTS OF FILE',MFR,' REQUIRE A SECTION MF=',  &
     &                M,'  AND MT=',MTR
                  WRITE(NOUT,'(5X,A)')  TRIM(EMESS)
!                 Increment count for errors
                  NERROR=NERROR+1
               ELSE
                  ICON(N,2) = ICON(N,2) - 2000
               END IF
            END DO
         END IF
!
!     ARE FILE 15 SECTIONS REQUIRED BY FILES 12 AND 13 PRESENT?
!
      ELSE IF(M.EQ.15) THEN
         IF(NPMT.EQ.0.OR.NSUB.NE.10)  THEN
            NMISS = -1
         ELSE
            DO N=1,NPMT
               IF(ICON(N,2).NE.0) THEN
                  MFR = ICON(N,1)/1000
                  MTR = MOD(ICON(N,1),1000)
                  NMISS = NMISS + 1
                  IF(MISFERR.NE.M) THEN
                     WRITE(NOUT,'(/2X,A,I5,A,I3)')                      &
     &               'ERROR(S) - MISSING SECTIONS IN MAT',MAT,' MF',M
                     MISFERR = M
                  END IF
                  WRITE(EMESS,'(A,I3,A,I3,A,I4)')                       &
     &                'CONTENTS OF FILE',MFR,' REQUIRE A SECTION MF=',  &
     &                M,'  AND MT=',MTR
                  WRITE(NOUT,'(5X,A)')  TRIM(EMESS)
!                 Increment count for errors
                  NERROR=NERROR+1
               END IF
            END DO
            NPMT = 0
         END IF
      ELSE
         NMISS = -1
      END IF
!
  100 CONTINUE
!
      RETURN
      END SUBROUTINE EFTEST
!
!***********************************************************************
!
      SUBROUTINE MISFIL(MF0,MT1,MT2,NMISS)
!
!     ROUTINE TO CHECK THAT FILE MFO CONTAINS SECTIONS MT1 THRU MT2
!       IF FILE 3 DOES.  A SECTION IN MF=6 ALSO SATISFIES TEST
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MF0,MT1,MT2,NMISS
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     CHARACTER(LEN=*), INTRINSIC :: TRIM
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: MFTL,MFTH,NNL,NDUM,MTT,MFT1,ISET
      INTEGER(KIND=I4) :: N,NN
!
      MFTL = 3000 + MT1
      MFTH = 3000 + MT2
      NNL = 1
      NDUM = 1
      DO N=1,NXC
         IF(MFTH.LT.INDX(N,1)) GO TO 100
         IF(MFTL.LE.INDX(N,1)) THEN
            IF(INDX(N,2).NE.1) THEN
               MTT = MOD(INDX(N,1),1000)
               MFT1 = 1000*MF0 + MTT
               DO NN=NNL,NXC
                  IF(INDX(NN,1).EQ.MFT1)  THEN
                     IF(INDX(NN,2).NE.1) THEN
                        NDUM = NN
                        GO TO 50
                     END IF
                  END IF
               END DO
               IF(MF0.NE.6) THEN
                  IF(NFOR.GE.6)   THEN
                     CALL TESTS(1000*6+MTT,ISET)
                     IF(ISET.LT.3)   GO TO 50
                  END IF
               END IF
               IF(MISFERR.NE.MF0) THEN
                  WRITE(NOUT,'(/2X,A,I5,A,I3)')                         &
     &            'ERROR(S) - MISSING SECTIONS IN MAT',MAT,' MF',MF0
                  MISFERR = MF0
               END IF
               WRITE(EMESS,'(A,I4,2A,I2)')                              &
     &              'PRESENCE OF FILE 3, MT=',MTT,' REQUIRES AN ',      &
     &              'EQUIVALENT SECTION IN FILE',MF0
               WRITE(NOUT,'(5X,A)')  TRIM(EMESS)
               NMISS = NMISS + 1
!              Increment count for errors
               NERROR=NERROR+1
            END IF
         END IF
   50    NNL = NDUM
      END DO
!
  100 RETURN
      END SUBROUTINE MISFIL
!
!***********************************************************************
!
      SUBROUTINE ISFIL(MF1,MF2,MT1,MT2)
!
!     CHECKS THAT MF2,MT2 EXISTS AND SPANS THE SAME ENERGY RANGE
!       AS MF1,MT1
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MF1,MF2,MT1,MT2
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: MFT,MFT1,IERR
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: ELO,EHI,ELOT,EHIT
!
      IERR = 0
!
!     CALCULATE INDEX VALUES FOR SECTIONS
!
      MFT = 1000*MF1 + MT1
      MFT1 = 1000*MF2 + MT2
!
!     FIND ENERGY RANGE SPANNED BY MF1/MT1
!
      DO N=1,NXC
         IF(INDX(N,1).GE.MFT) THEN
            IF(INDX(N,1).GT.MFT.OR.INDX(N,2).EQ.1) GO TO 20
            ELO = ENGS(N,1)
            EHI = ENGS(N,2)
            IF(EHI.EQ.0.0) GO TO 100
            GO TO 25
         END IF
      END DO
!
!     SECTION MF1/MT1 NOT IN MATERIAL
!
   20 WRITE(EMESS,'(A,I3,A,I4,A)')                                      &
     &     'FILE',MF1,', MT=',MT1,' IS MISSING FROM INDEX'
      CALL ERROR_MESSAGE(1)
      GO TO 100
!
!     FIND MF2/MT2 IN THE INDEX
!
   25 DO N=1,NXC
         IF(INDX(N,1).GE.MFT1)THEN
            IF(INDX(N,1).GT.MFT1.OR.INDX(N,2).EQ.1) GO TO 30
            GO TO 35
         END IF
      END DO
!
!     MESSAGE WHEN REQUIRED SECTION IS MISSING
!
   30 IF(MF1.NE.13)   THEN
         WRITE(EMESS,'(A,I3,A,I4,A)')                                   &
     &         'THIS SECTION REQUIRES THAT MISSING FILE',MF2,           &
     &         ', MT=',MT2,'  BE PRESENT'
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
      GO TO 100
!
!     CALCULATE FRACTIONAL DIFFERENCE BETWEEN EACH ENERGY LIMIT
!
   35 ELOT=ABS(ELO-ENGS(N,1))/AMAX1(ELO,ENMIN)
      EHIT=ABS(EHI-ENGS(N,2))/EHI
!
!     FOR FILE 15 NEED LOWEST GAMMA RAY ENERGY
!
      IF(MF1.EQ.15) THEN
         DO N=1,NMTGAM
            IF(MT1.EQ.MTGAM(N))   THEN
               ELOT=(ELO-EGAM(N))/AMAX1(ELO,ENMIN)
               GO TO 50
            END IF
         END DO
         GO TO 60
      END IF
!
!     TEST THE LIMITS
!
   50 IF(ELOT.GT.EPSILN4) IERR = 1
   60 IF(EHIT.GT.EPSILN4) IERR = 1
      IF(IERR.EQ.1) THEN
         WRITE(EMESS,'(A,I3,A,I4)')                                     &
     &        'SECTION DOES NOT SPAN THE SAME ENERGY RANGE AS FILE ',   &
     &        MF2,', MT=',MT2
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
  100 RETURN
      END SUBROUTINE ISFIL
!
!***********************************************************************
!
      SUBROUTINE TESTD(INDXT)
!
!     TEST THAT MF,MT IS IN INDEX
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDXT
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: MFC,MTC
      INTEGER(KIND=I4) :: IU,K
      INTEGER(KIND=I4) :: I,J,N
!
!     SEARCH INDEX FOR VALUE
!
      IF(NXC.GT.0) THEN
         DO N=1,NXC
            IF(INDX(N,1).EQ.INDXT) THEN
               INDX(N,2) = 0
               GO TO 100
            ELSE IF(INDX(N,1).GT.INDXT) THEN
               GO TO 20
            END IF
         END DO
      END IF
      N = NXC + 1
!
!     NEW ONE SO ADD IT TO THE INDEX
!
   20 MFC = INDXT/1000
      MTC = MOD(INDXT,1000)
      WRITE(EMESS,'(A,I3,A,I3,A)')                                      &
     &        'SECTION',MFC,'/',MTC,' NOT IN INDEX'
      CALL ERROR_MESSAGE(0)
      NERROR = NERROR + 1
      IF(NXC.GE.NSECMAX) GO TO 100
      NXC = NXC + 1
      IF(N.NE.NXC)  THEN
         IU = NXC - N
         DO I=1,IU
            K = NXC - I
            DO J=1,2
               ENGS(K+1,J) = ENGS(K,J)
               INDX(K+1,J) = INDX(K,J)
            END DO
         END DO
      END IF
      INDX(N,1) = INDXT
      INDX(N,2) = 2
      ENGS(N,1) = 0.
      ENGS(N,2) = 0.
!
  100 RETURN
      END SUBROUTINE TESTD
!
!***********************************************************************
!
      SUBROUTINE TESTS(INDXT,ISET)
!
!     CHECK MATERIAL INDEX TO SEE WHETHER OR NOT REACTION MT IS
!     IS PRESENT IN FILE MF.  NO CHANGES ARE MADE IN THE INDEX ARRAY.
!     THE STATUS IS INDICATED BY THE VARIABLE ISET AS FOLLOWS--
!       ISET=0, IF THE REACTION IS IN THE MATERIAL INDEX, AND HAS
!               BEEN LOCATED IN THE DATA FILE,
!           =1, IF THE REACTION IS IN THE INDEX, BUT HAS NOT YET BEEN
!               FOUND IN THE DATA,
!           =2, IF THE REACTION IS IN THE DATA, BUT NOT THE INDEX,
!           =3, IF THE REACTION IS IN NEITHER THE DATA NOR THE INDEX.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDXT,ISET
!
      INTEGER(KIND=I4) :: N
!
!     LOOK UP INDEX VALUE
!
      IF(NXC.GT.0)   THEN
         DO N=1,NXC
            IF(INDX(N,1).GT.INDXT) THEN
               GO TO 20
            ELSE IF(INDX(N,1).EQ.INDXT) THEN
               ISET = INDX(N,2)
               GO TO 100
            END IF
         END DO
      END IF
!
!     VALUE NOT FOUND
!
   20 ISET = 3
      GO TO 100
!
  100 RETURN
      END SUBROUTINE TESTS
!
!***********************************************************************
!
      SUBROUTINE TESTP(MFC,MTC)
!
!     TEST FOR THE PRESENCE OF REQUIRED SECTION MT IN FILE MF
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MFC,MTC
!
      INTEGER(KIND=I4) :: MFTEM,MTTEM,ISET
!
      MFTEM = MFC
      MTTEM = MTC
      CALL TESTS(1000*MFTEM+MTTEM,ISET)
      IF(ISET.EQ.1.OR.ISET.EQ.3)   THEN
         WRITE(EMESS,'(A,I3,A,I3)')                                     &
     &       'THIS SECTION REQUIRES THE PRESENCE OF SECTION',MFC,'/',MTC
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
      RETURN
      END SUBROUTINE TESTP
!
!***********************************************************************
!
      SUBROUTINE TESTER(ELO,EHI,Q)
!
!     TEST TO SEE IF THE MINIMUM INCIDENT PARTICLE ENERGY
!      IS AT THE THRESHOLD OR 1.0E-05 WHICHEVER IS HIGHER
!      AND MAXIMUM IS AT THE FILE MAXIMUM
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: ELO,EHI,Q
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      REAL(KIND=R4) :: ELOP,ELOT
!
!     CALCULATE THRESHOLD FROM Q-VALUE
!
      IF(Q.NE.QUNK)   THEN
         ELOP = -Q*(AWR+1.0)/AWR
!
!        COMPARE REQUIRED MINIMUM ENERGY WITH ACTUAL
!
         IF(NSUB/10.EQ.1)   THEN
            ELOT = AMAX1(ELOP,ENMIN)
!
!           LOWER ENERGY IN ERROR
!
            IF(ABS((ELO-ELOT)/ELOT).GT.EPSILN3) THEN
               WRITE(EMESS,'(A,1PE12.5,A)')                             &
     &            'THE MINIMUM INCIDENT ENERGY OF ',ELO,'(EV)'
               CALL ERROR_MESSAGE(0)
               WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5,A)')                &
     &            'SHOULD BE ',ELOT,'(EV) FOR Q= ',Q,'(EV)'
               CALL ERROR_MESSAGE(0)
               NERROR = NERROR + 1
            END IF
         END IF
      END IF
!
!     CHECK UPPER LIMIT
!
      IF(EHI.LT.ENMAX)   THEN
         WRITE(EMESS,'(A,1PE12.5,A)')                                   &
     &         'THE MAXIMUM INCIDENT ENERGY OF ',EHI,' (EV)'
         CALL ERROR_MESSAGE(0)
         WRITE(EMESS,'(3X,A,1PE12.5)')                                  &
     &         'SHOULD BE GREATER THAN OR EQUAL TO ',ENMAX
         CALL ERROR_MESSAGE(0)
         NERROR = NERROR + 1
      END IF
!
      RETURN
      END SUBROUTINE TESTER
!
!***********************************************************************
!
      SUBROUTINE TEST1
!
!     CHECK (X,Y) TABLE
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: EXP, ALOG
!
      INTEGER(KIND=I4) :: N1,N2,NL,NLP,NH,NA,NB,INONEG
      INTEGER(KIND=I4) :: I2,I3,II4,I5,JCKT
      INTEGER(KIND=I4) :: N,M,INEG,NSEQ1S
      REAL(KIND=R4) :: XN,X1,X2,X3,X4,X5,Y1,Y2,Y3,Y4,Y5
      REAL(KIND=R4) :: Y2P,Y4P,G,Z1,Z2,Z3,ZL,ZH,AL3,AL4,ARG
!
      REAL(KIND=R4), PARAMETER :: C=1.1
!
      INTEGER(KIND=I4), PARAMETER :: NMTNEG=9
      INTEGER(KIND=I4), DIMENSION(NMTNEG) , PARAMETER ::                &
     &           MTNEG = (/1,2,3,18,19,20,21,38,102/)
!
!     INITIALIZE DEVIANT POINT CHECK FLAG
!
      JCKT = FIZCON_DATA%ICKT
      n1 = 0
      N2 = 0
!
!     CHECK THAT DATA NOT ZERO OR NEGATIVE WHEN LOG INTERPOLATION
!       IS SPECIFIED
!
      CALL TEST2
      IF(MESS.NE.0)   JCKT = 1
!
!     CHECK THAT X ARE IN ORDER
!
      CALL TEST5X(1)
      IF(MESS.NE.0)   JCKT = 1
!
!     SKIP FOR ANOMALOUS SCATTERING FACTORS AND LN(S(ALPHA,BETA))
!
      IF(MT.EQ.505.OR.MT.EQ.506.OR.INEGC.EQ.0)  GO TO 30
!
!     FIND RESONANCE REGION
!
      N1 = NP + 1
      N2 = 0
      DO N=1,NP
         XN = X(N)
         IF(XN.GE.E1.AND.XN.LE.E2) THEN
            IF(N2.EQ.0)   N1 = N
            N2 = N
         END IF
      END DO
!
!     CHARGED PARTICLE ELASTIC; TREAT AS ONE BIG RESONANCE REGION
!
      IF(MT.EQ.2.AND.NSUB.GT.10.AND.NSUB/10.GT.1) THEN
         N1 = 1
         N2 = NP
      END IF
!
!     IF MF=3 AND MT=251, CHECK THAT ALL VALUES ARE IN THE RANGE
!       (-1.0,1.0).
!
      IF(MF.EQ.3.AND.MT.EQ.251)  THEN
         CALL TEST6Y(-1.0,1.0,'MUBR')
         GO TO 30
      END IF
!
!     CHECK THAT THERE ARE NO NEGATIVE CROSS SECTIONS, EXCEPT FOR
!       CERTAIN MT'S IN THE RESONANCE REGION
!
      INONEG = 1
      DO INEG=1,NMTNEG
         IF(MT.EQ.MTNEG(INEG))  THEN
            INONEG = 0
            GO TO 20
         END IF
      END DO
!
!     TEST FOR NEGATIVE NUMBERS
!
   20 DO N=1,NP
         IF(Y(N).LT.0.0) THEN
            IF(INONEG.EQ.1.OR.(N.LT.N1.OR.N.GT.N2)) THEN
               CALL RESEQ(NSEQP1,N,NR,NSEQ1S)
               WRITE(EMESS,'(A,I5)')                                    &
     &            'NEGATIVE FUNCTION VALUE NEAR N2=',N
               CALL ERROR_MESSAGE(NSEQ1S)
            END IF
         END IF
      END DO
!
!     MESSAGE IF DEVIANT POINT TEST SUPPRESSED
!
   30 IF(JCKT.NE.FIZCON_DATA%ICKT)  THEN
         EMESS = 'DEVIANT POINT CHECK SUPPRESSED DUE TO DATA ERRORS'
         CALL ERROR_MESSAGE(0)
      END IF
!
!     DO DEVIANT POINT TEST IF THERE ARE MORE THAN 4 POINTS
!
      IF(JCKT.EQ.1.OR.NP.LT.5) GO TO 90
      MESS = 0
      NL = 1
!
!     CHECK FOR DISCONTINUITIES
!
   40 NLP = NL + 1
      DO N=NLP,NP
         NH = N - 1
         IF(X(N).LE.X(N-1))  GO TO 45
      END DO
      NH = NP
   45 IF((NH-NL).LT.4)   GO TO 85
!
!     SKIP RESONANCE REGION IF IT EXISTS
!
      IF(NH.GT.N1.AND.NL.LT.N2) THEN
         IF(NL.LE.N1) THEN
            NH = N1
            GO TO 50
         END IF
         IF(NH.LT.N2) THEN
            NH = N2 - 1
            GO TO 85
         END IF
         NL = N2
      END IF
   50 IF((NH-NL).LT.4) GO TO 85
!
!     TEST ALL VALUES IN A CONTINUOUS REGION
!
      NA = NL + 2
      NB = NH - 2
      DO 80 N=NA,NB
!
!     SAVE DATA POINTS
!
      Y1 = Y(N-2)
      Y2 = Y(N-1)
      Y3 = Y(N)
      Y4 = Y(N+1)
      Y5 = Y(N+2)
      X1 = X(N-2)
      X2 = X(N-1)
      X3 = X(N)
      X4 = X(N+1)
      X5 = X(N+2)
!
!     GET INTERPOLATION LAWS BETWEEN THE POINTS
!
      M = 1
!*****BETWEEN POINTS 1 AND 2
   55 IF((N-1).GT.NBT(M)) THEN
         M = M + 1
         IF(M.GT.NR) GO TO 90
         GO TO 55
      END IF
      I2 = JNT(M)
      IF(I2.LT.1.OR.I2.GT.5)   GO TO 90
!*****BETWEEN POINTS 2 AND 3
   60 IF(N.GT.NBT(M))  THEN
         M = M + 1
         IF(M.GT.NR)   GO TO 90
         GO TO 60
      END IF
      I3 = JNT(M)
      IF(I3.LT.1.OR.I3.GT.5)   GO TO 90
!*****BETWEEN POINTS 3 AND 4
   65 IF((N+1).GT.NBT(M)) THEN
         M = M + 1
         IF(M.GT.NR)  GO TO 90
         GO TO 65
      END IF
      II4 = JNT(M)
      IF(II4.LT.1.OR.II4.GT.5)   GO TO 90
!*****BETWEEN POINTS 4 AND 5
   70 IF((N+2).GT.NBT(M)) THEN
         M = M + 1
         IF(M.GT.NR)  GO TO 90
         GO TO 70
      END IF
      I5 = JNT(M)
      IF(I5.LT.1.OR.I5.GT.5)   GO TO 90
!
!     TEST A SINGLE VALUE
!
      IF(I2.EQ.1) THEN
        Y2P = 0.0
      ELSE IF(I2.EQ.2) THEN
        Y2P = (Y2-Y1)/(X2-X1)
      ELSE IF(I2.EQ.3) THEN
        Y2P = (Y2-Y1)/X2/ALOG(X2/X1)
      ELSE IF(I2.EQ.4) THEN
        Y2P = Y2*ALOG(Y2/Y1)/(X2-X1)
      ELSE IF(I2.EQ.5) THEN
        Y2P = Y2*ALOG(Y2/Y1)/X2/ALOG(X2/X1)
      END IF
!
      IF(I3.EQ.1) THEN
         Z1 = Y2
      ELSE IF(I3.EQ.2) THEN
         Z1 = Y2+Y2P*(X3-X2)
      ELSE IF(I3.EQ.3) THEN
         Z1 = Y2+Y2P*X2*ALOG(X3/X2)
      ELSE IF(I3.EQ.4) THEN
         ARG = Y2P*(X3-X2)/Y2
         IF(ARG.GT.50.) THEN
            ARG = 50.
         ELSE IF(ARG.LT.-50.) THEN
            ARG  =  - 50.
         END IF
         Z1 = Y2*EXP(ARG)
      ELSE IF(I3.EQ.5) THEN
         ARG = Y2P*X2*ALOG(X3/X2)/Y2
         IF(ARG.GT.50.) THEN
            ARG = 50.
         ELSE IF(ARG.LT.-50.) THEN
            ARG  =  - 50.
         END IF
         Z1 = Y2*EXP(ARG)
      END IF
!
      IF(I5.EQ.1) THEN
         Y4P = 0.0
      ELSE IF(I5.EQ.2) THEN
         Y4P = (Y5-Y4)/(X5-X4)
      ELSE IF(I5.EQ.3) THEN
         Y4P = (Y5-Y4)/X4/ALOG(X5/X4)
      ELSE IF(I5.EQ.4) THEN
         Y4P = Y4*ALOG(Y5/Y4)/(X5-X4)
      ELSE IF(I5.EQ.5) THEN
         Y4P = Y4*ALOG(Y5/Y4)/X4/ALOG(X5/X4)
      END IF
!
      IF(II4.EQ.1) THEN
         Z2 = Y4
      ELSE IF(II4.EQ.2) THEN
         Z2 = Y4+Y4P*(X3-X4)
      ELSE IF(II4.EQ.3) THEN
         Z2 = Y4+Y4P*X4*ALOG(X3/X4)
      ELSE IF(II4.EQ.4) THEN
         ARG = Y4P*(X3-X4)/Y4
         IF(ARG.GT.50.) THEN
            ARG = 50.
         ELSE IF(ARG.LT.-50.) THEN
            ARG  =  - 50.
         END IF
         Z2 = Y4*EXP(ARG)
      ELSE IF(II4.EQ.5) THEN
         ARG = Y4P*X4*ALOG(X3/X4)/Y4
         IF(ARG.GT.50.) THEN
            ARG = 50.
         ELSE IF(ARG.LT.-50.) THEN
            ARG  =  - 50.
         END IF
         Z2 = Y4*EXP(ARG)
      END IF
!
      IF(I3.EQ.3.OR.I3.EQ.5) THEN
         AL3 = X3*ALOG(X3/X2)
      ELSE
         AL3 = X3-X2
      END IF
!
      IF(II4.EQ.3.OR.II4.EQ.5) THEN
         AL4 = X3*ALOG(X4/X3)
      ELSE
         AL4 = X4-X3
      END IF
      G = AL3/AL4
!
      IF(I3.LE.3.AND.II4.LE.3)   THEN
         Z3 = (Y2+G*Y4)/(1.0+G)
      ELSE
         IF(Y2.GT.0..AND.Y4.GT.0.)   THEN
            Z3 = EXP((ALOG(Y2)+G*ALOG(Y4))/(1.0+G))
         ELSE
            Z3 = (Y2+G*Y4)/(1.0+G)
         END IF
      END IF
!
      ZH = Z1
      IF(Z2.GT.ZH)   THEN
         ZH = Z2
      END IF
      IF(Z3.GT.ZH)   THEN
         ZH = Z3
      END IF
!
      ZL = Z1
      IF(Z2.LT.ZL)   THEN
         ZL = Z2
      END IF
      IF(Z3.LT.ZL)   THEN
         ZL = Z3
      END IF
!
      IF(ZH.GE.0.0)    THEN
         IF(Y3.GT.C*ZH)   GO TO 75
      ELSE
         IF(Y3.GT.ZH/C.OR.Y3.GT.C*ZH)   GO TO 75
      END IF
!
      IF(ZL.GE.0.0)    THEN
         IF(Y3.LT.ZL/C)   GO TO 75
      ELSE
         IF(Y3.LT.C*ZL.OR.Y3.LT.ZL/C)   GO TO 75
      END IF
      GO TO 80
!
!     DATA POINT FAILS THE TEST
!
   75 CALL RESEQ(NSEQP1,N,NR,NSEQ1S)
      WRITE(EMESS,'(A,I5)')  'CHECK POINT',N
      CALL ERROR_MESSAGE(NSEQ1S)
      WRITE(EMESS,'(4X,A,1PE13.5,A,1PE13.5)') '  X =',X3,'  Y =',Y3
      CALL ERROR_MESSAGE(0)
      MESS = MESS + 1
      IF(MESS.GE.MAXMES)    THEN
         EMESS = 'AND MAYBE OTHERS FAIL THE DEVIANT POINT CHECK'
         CALL ERROR_MESSAGE(0)
         GO TO 90
      END IF
!
!     END POINT LOOP
!
   80 CONTINUE
!
!     CHECK FOR MORE REGIONS
!
   85 NL = NH + 1
      IF((NP-NL).GE.4)  GO TO 40
!
!     TEST FOR CONSECUTIVE EQUAL VALUES
!
   90 IF(FIZCON_DATA%ICKT.EQ.0)   CALL TEST4
!
      RETURN
      END SUBROUTINE TEST1
!
!***********************************************************************
!
      SUBROUTINE TEST2
!
!     TEST THAT NO ZEROS OR NEGATIVE VALUES OCCUR IN LOG INTERPOLATION
!        REGIONS OF A TAB1 RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: M,I,INLOG,NSEQ1S
      INTEGER(KIND=I4) :: N
!
!     CHECK THAT DATA NOT ZERO OR NEGATIVE WHEN LOG INTERPOLATION
!       IS SPECIFIED
!
      MESS = 0
      M = 1
      DO N=2,NP
!*****DETERMINE INTERPOLATION CODE
   10    IF(N.LE.NBT(M))   GO TO 20
         M = M + 1
         IF(M.LE.NR)   GO TO 10
         M = NR
   20    I = JNT(M)
!
!        CHECK FOR LE ZERO FOR A LOG INTERPOLATION SCHEME
!
         INLOG = 0
         IF(I.EQ.3) THEN
            IF(X(N-1).LE.0..OR.X(N).LE.0.)   INLOG = 1
         ELSE IF(I.EQ.4) THEN
            IF(Y(N-1).LE.0..OR.Y(N).LE.0.)   INLOG = 1
         ELSE IF(I.EQ.5) THEN
            IF(X(N-1).LE.0..OR.X(N).LE.0.)   INLOG = 1
            IF(Y(N-1).LE.0..OR.Y(N).LE.0.)   INLOG = 1
         END IF
!
!     ERROR MESSAGE
!
         IF(INLOG.EQ.1) THEN
            CALL RESEQ(NSEQP1,N,NR,NSEQ1S)
            WRITE(EMESS,'(A,I5)')                                       &
     &              'NEG OR ZERO ARG OF LOG BELOW POINT',N
            CALL ERROR_MESSAGE(NSEQ1S)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)    THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               CALL ERROR_MESSAGE(0)
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST2
!
!***********************************************************************
!
      SUBROUTINE TEST3(N1,N2,KXXX)
!
!     ROUTINE TO TEST THAT N1=N2
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: KXXX
      INTEGER(KIND=I4) :: N1,N2
!
      IF(N1.NE.N2)  THEN
         WRITE(EMESS,'(A,A,I5)')  KXXX,' SHOULD BE SET TO ',N2
         CALL ERROR_MESSAGE(NSEQP1)
      END IF
!
      RETURN
      END SUBROUTINE TEST3
!
!***********************************************************************
!
      SUBROUTINE TEST3F(A,B,KXXX)
!
!     ROUTINE TO TEST FOR EQUALITY OF FLOATING POINT NUMBERS, A = B
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: KXXX
      REAL(KIND=R4) :: A,B
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      REAL(KIND=R4) :: RTEST
!
      IF(B.NE.0.)  THEN
         RTEST = ABS((B-A)/B)
      ELSE
         RTEST = ABS(A)
      END IF
      IF(RTEST.GT.EPSILN6)  THEN
         WRITE(EMESS,'(A4,A,1PE12.5)')  KXXX,' SHOULD BE SET TO ',B
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!
      RETURN
      END SUBROUTINE TEST3F
!
!***********************************************************************
!
      SUBROUTINE TEST4
!
!     TESTS ARRAY Y(N) TO SEE IF THERE ARE 3 CONSECUTIVE EQUAL VALUES
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NCON,NSEQ1S
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: YN,YN1
!
      MESS = 0
      IF(NP.LT.3) GO TO 100
!
!     INITIALIZE
!
      NCON = 1
      YN = Y(1)
!
!     CHECK EACH POINT
!
      DO N=2,NP
         YN1 = YN
         YN = Y(N)
         IF(YN1.NE.YN) NCON = 0
         NCON = NCON + 1
         IF(NCON.EQ.3) THEN
            CALL RESEQ(NSEQP1,N,NR,NSEQ1S)
            WRITE(EMESS,'(A,I5)')                                       &
     &           'Y LIST HAS 3 OR MORE EQUAL VALUES N=',N
            CALL ERROR_MESSAGE(NSEQ1S)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)    THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST4
!
!***********************************************************************
!
      SUBROUTINE TEST5(Z,NPT,L)
!
!     TEST THAT THE ARRAY Z(N),N=1,NPT,L IS IN INCREASING ORDER
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L
      INTEGER(KIND=I4) :: NPT
      REAL(KIND=R4), DIMENSION(NPT) :: Z
!
      INTEGER(KIND=I4) :: K
      INTEGER(KIND=I4) :: N
!
      MESS = 0
      IF(NPT.LE.1.OR.L.GE.NPT.OR.L.LT.1)    GO TO 100
      K = L + 1
      DO N=K,NPT,L
         IF(Z(N).LT.Z(N-L))   THEN
!
!           ERROR MESSAGE
!
            WRITE(EMESS,'(A,I5)')  'LIST OUT OF ORDER NEAR N=',N
            CALL ERROR_MESSAGE(NSEQP)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)   THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               CALL ERROR_MESSAGE(0)
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST5
!
!***********************************************************************
!
      SUBROUTINE TEST5B(Z,NPT,L)
!
!     TEST THAT ARRAY Z(N),N=1,NPT,L IS IN DECREASING ORDER
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L
      INTEGER(KIND=I4) :: NPT
      REAL(KIND=R4), DIMENSION(NPT) :: Z
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: K
      INTEGER(KIND=I4) :: N
!
      MESS = 0
      IF(NPT.LE.1.OR.L.GE.NPT.OR.L.LT.1)   GO TO 100
      K = L + 1
      DO N=K,NPT,L
         IF(ABS(Z(N-L)).LT.ABS(Z(N))) THEN
!
!           ERROR MESSAGE
!
            WRITE(EMESS,'(A,I5)')  'LIST OUT OF ORDER NEAR N=',N
            CALL ERROR_MESSAGE(NSEQP)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)   THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               CALL ERROR_MESSAGE(0)
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST5B
!
!***********************************************************************
!
      SUBROUTINE TEST5X(L)
!
!     TEST THAT ARRAY X(N),N=1,NP,L IS IN INCREASING ORDER
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L
!
      INTEGER(KIND=I4) :: K
      INTEGER(KIND=I4) :: N
!
      MESS = 0
      IF(NP.LE.1.OR.L.GE.NP.OR.L.LT.1)   GO TO 100
      K = L + 1
      DO N=K,NP,L
         IF(X(N).LT.X(N-L))   THEN
!
!           ERROR MESSAGE
!
            WRITE(EMESS,'(A,I5)')  'LIST OUT OF ORDER NEAR N=',N
            CALL ERROR_MESSAGE(NSEQP)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)   THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               CALL ERROR_MESSAGE(0)
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST5X
!
!***********************************************************************
!
      SUBROUTINE TEST5Y(NBEG,NPT,L,IUP)
!
!     TEST THAT ARRAY Y(N),N=NBEG,NPT,L IS IN INCREASING (IUP=1) OR
!      DECREASING (IUP=0) ORDER
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NBEG,NPT,L,IUP
!
      INTEGER(KIND=I4) :: K,ILOO
      INTEGER(KIND=I4) :: N
!
      MESS = 0
      IF(NPT.LE.NBEG.OR.NBEG+L.GT.NPT.OR.L.LT.1)   GO TO 100
      K = NBEG + L
      DO N=K,NPT,L
         ILOO = 0
         IF(IUP.NE.0)  THEN
            IF(Y(N).LT.Y(N-L))   ILOO = 1
         ELSE
            IF(Y(N).GT.Y(N-L))   ILOO = 1
         END IF
!
!        ERROR MESSAGE
!
         IF(ILOO.EQ.1) THEN
            WRITE(EMESS,'(A,I5)')  'LIST OUT OF ORDER NEAR N=',N
            CALL ERROR_MESSAGE(NSEQP)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)   THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               CALL ERROR_MESSAGE(0)
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST5Y
!
!***********************************************************************
!
      SUBROUTINE TEST6(Z,A,B,KXXX)
!
!     CHECK THAT Z LIES BETWEEN A AND B
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: KXXX
      REAL(KIND=R4) :: Z,A,B
!
      IF(Z.LT.A.OR.Z.GT.B)   THEN
!
!        ERROR MESSAGE
!
         WRITE(EMESS,'(A4,A,1PE13.5,A,1PE13.5)')                        &
     &               KXXX,' NOT IN RANGE',A,' TO',B
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!
      RETURN
      END SUBROUTINE TEST6
!
!***********************************************************************
!
      SUBROUTINE TEST6I(IZ,IA,IB,KXXX)
!
!     CHECK THAT IZ LIES BETWEEN IA AND IB
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: KXXX
      INTEGER(KIND=I4) :: IZ,IA,IB
!
      IF(IZ.LT.IA.OR.IZ.GT.IB)   THEN
!
!        ERROR MESSAGE
!
         WRITE(EMESS,'(A4,I8,A,I8,A,I8)')                               &
     &               KXXX,IZ,' NOT IN RANGE',IA,' TO',IB
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!
      RETURN
      END SUBROUTINE TEST6I
!
!***********************************************************************
!
      SUBROUTINE TEST6X(A,B,KXXX)
!
!     CHECK THAT ALL X(N),N=1,NP LIE BETWEEN A AND B
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: KXXX
      REAL(KIND=R4) :: A,B
!
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: ZN
!
      MESS = 0
      DO N=1,NP
         ZN = X(N)
         IF(ZN.LT.A.OR.ZN.GT.B) THEN
!
!           ERROR MESSAGE
!
            WRITE(EMESS,'(A4,A,1PE13.5,A,1PE13.5)')                     &
     &               KXXX,' NOT IN RANGE',A,' TO',B
            CALL ERROR_MESSAGE(NSEQP)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)   THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               CALL ERROR_MESSAGE(0)
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST6X
!
!***********************************************************************
!
      SUBROUTINE TEST6Y(A,B,KXXX)
!
!     CHECK THAT ALL Y(N),N=1,NP LIE BETWEEN A AND B
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: KXXX
      REAL(KIND=R4) :: A,B
!
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: ZN
!
      MESS = 0
      DO N=1,NP
         ZN = Y(N)
         IF(ZN.LT.A.OR.ZN.GT.B)   THEN
!
!           ERROR MESSAGE
!
            WRITE(EMESS,'(A4,A,1PE13.5,A,1PE13.5)')                     &
     &               KXXX,' NOT IN RANGE',A,' TO',B
            CALL ERROR_MESSAGE(NSEQP)
            MESS = MESS + 1
            IF(MESS.GE.MAXMES)   THEN
               EMESS = 'AND MAYBE OTHERS FAIL THE TEST'
               CALL ERROR_MESSAGE(0)
               GO TO 100
            END IF
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TEST6Y
!
!***********************************************************************
!
      SUBROUTINE TEST7(ANS,IPATH)
!
!     CHECK THAT Y(X) IS NORMALIZED TO UNITY
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPATH
      REAL(KIND=R4) :: ANS
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: M,I
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: XT1,XT2,YT1,YT2,AN
!
!     INITIALIZE INTEGRAL
!
      ANS = 0.0
!
!     PROCESS EACH PANEL
!
      M = 1
      DO N=2,NP
!
!        DETERMINE INTERPOLATION IN THE PANEL
!
         IF(N.LE.NBT(M))   GO TO 40
         M = M + 1
         IF(M.GT.NR)   GO TO 100
   40    I = JNT(M)
!
!        SAVE FUNCTION VALUES
!
         XT1 = X(N-1)
         XT2 = X(N)
         YT1 = Y(N-1)
         YT2 = Y(N)
!
!        DO INTEGRAL
!
         CALL ECSI(XT1,YT1,XT2,YT2,XT1,XT2,I,AN)
         ANS = ANS + AN
      END DO
!
!     CHECK THAT INTEGRAL IS 1.0
!
      IF(IPATH.NE.2)   THEN
         IF(ABS(ANS-1.0).LE.EPSILN3)   GO TO 100
         WRITE(EMESS,'(A,1PE12.5,A)')                                   &
     &          'NORMALIZATION CHECK INTEGRAL=',ANS,' BEFORE '
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!
  100 RETURN
      END SUBROUTINE TEST7
!
!***********************************************************************
!
      SUBROUTINE RESEQ(NSEQP1,N,NN,NSEQ1)
!
!     CALCULATES CARD SEQUENCE NUMBER GIVEN THE NUMBER OF AN
!       ELEMENT IN A TAB1 RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NSEQP1,N,NN,NSEQ1
!
      NSEQ1 = NSEQP1 + (NN+2)/3 + (N+2)/3
!
      RETURN
      END SUBROUTINE RESEQ
!
!***********************************************************************
!
      SUBROUTINE RDTEXT
!
!     PROCESS A TEXT RECORD
!
      IMPLICIT NONE
!
!     READ IN RECORD
!
      READ(JIN,'(A66,I4,I2,I3,I5)',END=90)  TEXT,MATP,MFP,MTP,NSEQP
      GO TO 100
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE RDTEXT
!
!***********************************************************************
!
      SUBROUTINE RDCONT
!
!     PROCESS A CONT RECORD
!
      IMPLICIT NONE
!
!     READ IN RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90)                      &
     &                   C1H,C2H,L1H,L2H,N1H,N2H,MATP,MFP,MTP,NSEQP
      NSEQP1 = NSEQP
      GO TO 100
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE RDCONT
!
!***********************************************************************
!
      SUBROUTINE RDHEAD(I)
!
!     PROCESS A HEAD RECORD
!
!     SUBROUTINE RETURNS I WITH THE FOLLOWING MEANINGS
!          I=1, HEAD RECORD
!            2, SEND RECORD
!            3, FEND RECORD
!            4, MEND RECORD
!            5, TEND RECORD
!            0, AN ERROR
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: I
!
!     INITIALIZE RETURN FLAG
!
      I = 0
!
!     READ IN RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90)                      &
     &                 C1H,C2H,L1H,L2H,N1H,N2H,MAT,MF,MT,NSEQP
      NSEQP1 = NSEQP
!
!     LOOK FOR A TEND RECORD
!
      IF(MAT.LT.0)   THEN
         I = 5
         GO TO 100
      END IF
!
!     LOOK FOR A HEAD CARD
!
      IF(MT.NE.0)   THEN
         IF(MT.LT.0.OR.MF.LE.0.OR.MAT.EQ.0)   GO TO 100
         I = 1
         GO TO 100
      END IF
!
!     LOOK FOR SEND RECORD
!
      IF(MF.NE.0)   THEN
         IF(MF.LT.0.OR.MAT.EQ.0)   GO TO 100
         I = 2
         GO TO 100
      END IF
!
!     LOOK FOR A FEND RECORD
!
      IF(MAT.NE.0) THEN
         I = 3
         GO TO 100
      END IF
!
!     A MEND RECORD
!
      I = 4
      GO TO 100
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE RDHEAD
!
!***********************************************************************
!
      SUBROUTINE RDTAB2
!
!     PROCESS A TAB2 RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NI,NF
      INTEGER(KIND=I4) :: N,NN
!
!     READ IN CONT-LIKE RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90)                      &
     &         C12,C22,L12,L22,NR2,NP2,MATP,MFP,MTP,NSEQP
      NSEQP1 = NSEQP
!
!     READ IN INTERPOLATION TABLE
!
      NI = 1
      DO N=1,NR2,3
         NF = NI + 2
         READ(JIN,'(6I11,I4,I2,I3,I5)',END=90)                          &
     &                (NBT2(NN),JNT2(NN),NN=NI,NF),NSEQ
         NI = NI + 3
      END DO
      GO TO 100
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE RDTAB2
!
!***********************************************************************
!
      SUBROUTINE RDTAB1
!
!     PROCESS A TAB1 RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NI,NF,INP,NPP,IFST
      INTEGER(KIND=I4) :: N,NN
!
!     READ CONT-LIKE RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90)                      &
     &                 C1,C2,L1,L2,NR,NP,MATP,MFP,MTP,NSEQP
      NSEQP1 = NSEQP
!
!     READ IN INTERPOLATION TABLE
!
      NI = 1
      DO N=1,NR,3
         NF = NI + 2
         READ(JIN,'(6I11,I4,I2,I3,I5)',END=90)                          &
     &             (NBT(NN),JNT(NN),NN=NI,NF),MATP,MFP,MTP,NSEQP
         NI = NI + 3
      END DO
!
!     READ IN DATA TABLE
!
      NI = 1
      NPP = 0
      IPAGE = 1
      IFST = 0
!
!     READ AND STORE DATA A PAGE AT A TIME
!
   10 NF = NI + 2
      READ(JIN,'(6E11.4,I4,I2,I3,I5)',END=90)                           &
     &            (XP(N),YP(N),N=NI,NF),MATP,MFP,MTP,NSEQP
      NI = NI + 3
      NPP = NPP + 3
!*****OUTPUT LAST PAGE
      IF(NPP.GE.NP) THEN
         INP = -1
         CALL OUTXY(INP,1,NF)
!********TEST DATA TABLE
         CALL TEST1
         GO TO 100
!*****OUTPUT NEXT PAGE
      ELSE IF(NF.EQ.PAGESZ)  THEN
         CALL OUTXY(IFST,1,NF)
         NI = 7
      END IF
      GO TO 10
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE RDTAB1
!
!***********************************************************************
!
      SUBROUTINE RDINTG(NDIGIT)
!
!     Process a compact covariance format record
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NDIGIT,NROW,KIJ,II,JJ
      PARAMETER (NROW=18)
      DIMENSION KIJ(NROW)
      IF (NDIGIT.EQ.2)
     & READ (JIN,20,END=90,ERR=90) II,JJ,KIJ(1:18),MATP,MFP,MTP
      IF (NDIGIT.EQ.3)
     & READ (JIN,30,END=90,ERR=90) II,JJ,KIJ(1:13),MATP,MFP,MTP
      IF (NDIGIT.EQ.4)
     & READ (JIN,40,END=90,ERR=90) II,JJ,KIJ(1:11),MATP,MFP,MTP
      IF (NDIGIT.EQ.5)
     & READ (JIN,50,END=90,ERR=90) II,JJ,KIJ(1:9),MATP,MFP,MTP
      IF (NDIGIT.EQ.6)
     & READ (JIN,60,END=90,ERR=90) II,JJ,KIJ(1:1),MATP,MFP,MTP
      GO TO 100
   20 FORMAT (2I5, 1X, 18I3, 1X, I4, I2, I3)
   30 FORMAT (2I5, 1X, 13I4, 3X, I4, I2, I3)
   40 FORMAT (2I5, 1X, 11I5, I4, I2, I3)
   50 FORMAT (2I5, 1X, 9I6, 1X, I4, I2, I3)
   60 FORMAT (2I5, 8I7, I4, I2, I3)
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE RDINTG
!
!***********************************************************************
!
      SUBROUTINE RDLIST
!
!     PROCESS A LIST RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NI,NPP,IFST,NF,INP
      INTEGER(KIND=I4) :: N
!
!     READ CONT-LIKE RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90)                      &
     &              C1L,C2L,L1L,L2L,NPL,N2L,MATP,MFP,MTP,NSEQP
      NSEQP1 = NSEQP
!
!     READ IN DATA TABLE
!
      NI = 1
      NPP = 0
      IPAGE = 1
      IFST = 0
!
!     READ AND STORE DATA A PAGE AT A TIME
!
   10 NF = NI + 5
      READ(JIN,'(6E11.4,I4,I2,I3,I5)',END=90)                           &
     &                   (YP(N),N=NI,NF),MATP,MFP,MTP,NSEQP
      NI = NI + 6
      NPP = NPP + 6
!*****OUTPUT LAST PAGE
      IF(NPP.GE.NPL) THEN
         INP = -1
         CALL OUTXY(INP,2,NF)
         GO TO 100
!*****OUTPUT NEXT PAGE
      ELSE IF(NF.EQ.PAGESZ) THEN
         CALL OUTXY(IFST,2,NF)
         NI = 7
      END IF
      GO TO 10
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
!
  100 RETURN
      END SUBROUTINE RDLIST
!
!***********************************************************************
!
      SUBROUTINE OUTXY(N0,IW,NF)
!
!     ROUTINE TO PAGE OUT THE X AND/OR Y ARRAYS
!
!       N0 = 0: FIRST PAGE
!          = POSITIVE: CONTINUING PAGE(S)
!          = NEGATIVE: LAST PAGE
!       IW = 1: DO BOTH X AND Y
!            2: DO Y ONLY
!       NF = NO. OF POINTS ON THE PAGE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N0,IW,NF
!
      INTEGER(KIND=I4) :: K
!
!     FIRST PAGE INITIALIZATION
!
      IF(N0.EQ.0) THEN
         IF(IW.EQ.1)   REWIND (UNIT=ISCRX)
         REWIND (UNIT=ISCRY)
         IPAGE = 1
      END IF
!
!     LOAD PAGE IN SCRATCH FILE
!
      IF(N0.GE.0) THEN
         IF(IW.EQ.1)  THEN
            WRITE(ISCRX)    NF,(XP(K),K=1,NF)
            DO K=1,6
               XP(K) = XP(K+PAGESZ-6)
            END DO
         END IF
         WRITE(ISCRY)   NF,(YP(K),K=1,NF)
         DO K=1,6
            YP(K) = YP(K+PAGESZ-6)
         END DO
         IPAGE = IPAGE + 1
         N0 = N0 + 1
         GO TO 100
!
!     LAST BLOCK
!
      ELSE IF(N0.LT.0) THEN
         IF(IPAGE.LE.0)   IPAGE = 1
!*****X PAGE
         IF(IW.EQ.1)  THEN
            IPAGEX = IPAGE
            ILOWX = (PAGESZ-6)*(IPAGE-1)
            IHIGHX = ILOWX + PAGESZ
            IF(IPAGE.GT.1)THEN
               WRITE(ISCRX)   NF,(XP(K),K=1,NF)
               ENDFILE (UNIT=ISCRX)
            END IF
         END IF
!*****Y PAGE
         IPAGEY = IPAGE
         ILOWY = (PAGESZ-6)*(IPAGE-1)
         IHIGHY = ILOWY + PAGESZ
         IF(IPAGE.GT.1)  THEN
            WRITE(ISCRY)   NF,(YP(K),K=1,NF)
            ENDFILE (UNIT=ISCRY)
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE OUTXY
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION X(INDEXP)
!
!     LOGICAL PAGING SYSTEM FOR X ARRAY.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDEXP
!
      INTEGER(KIND=I4) :: INDEXX
!
!     POINT NOT ON CURRENT PAGE
!
      INDEXX = INDEXP
      IF(INDEXX.LE.ILOWX.OR.INDEXX.GT.IHIGHX)     CALL LOADX(INDEXX)
!
!     RETURN REQUIRED DATA VALUE
!
      X = XP(INDEXX-ILOWX)
!
      RETURN
      END FUNCTION X
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION Y(INDEXP)
!
!     LOGICAL PAGING SYSTEM FOR Y ARRAY.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDEXP
!
      INTEGER(KIND=I4) :: INDEXX
!
!     POINT NOT ON CURRENT PAGE
!
      INDEXX = INDEXP
      IF(INDEXX.LE.ILOWY.OR.INDEXX.GT.IHIGHY)    CALL LOADY(INDEXX)
!
!     RETURN REQUESTED VALUE
!
      Y = YP(INDEXX-ILOWY)
!
      RETURN
      END FUNCTION Y
!
!***********************************************************************
!
      SUBROUTINE LOADX(INDEXX)
!
!     LOGICAL PAGE LOADER FOR X ARRAY
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDEXX
!
      INTEGER(KIND=I4) :: JPAGE,ITOP
      INTEGER(KIND=I4) :: K
!
!     DETERMINE WHICH PAGE TO LOAD.
!
      JPAGE = (INDEXX-7)/(PAGESZ-6) + 1
!
!     IF CURRENT PAGE IS PAST REQUESTED PAGE REWIND SCRATCH
!
      IF(JPAGE.LE.IPAGEX) THEN
         REWIND (UNIT=ISCRX)
         IPAGEX = 0
      END IF
!
!     SKIP UP TO PAGE PRECEDING ONE REQUIRED
!
   10 IPAGEX = IPAGEX + 1
      IF(IPAGEX.NE.JPAGE)  THEN
         READ(ISCRX) ITOP
         GO TO 10
      END IF
!
!     LOAD CORRECT PAGE.
!
      READ(ISCRX) ITOP,(XP(K),K=1,ITOP)
!
!     SET LOWER AND UPPER INDICES.
!
      ILOWX = (PAGESZ-6)*(IPAGEX-1)
      IHIGHX = ILOWX + PAGESZ
!
      RETURN
      END SUBROUTINE LOADX
!
!***********************************************************************
!
      SUBROUTINE LOADY(INDEXX)
!
!     LOGICAL PAGE LOADER FOR Y ARRAY
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INDEXX
!
      INTEGER(KIND=I4) :: JPAGE,ITOP
      INTEGER(KIND=I4) :: K
!
!     DETERMINE WHICH PAGE TO LOAD.
!
      JPAGE = (INDEXX-7)/(PAGESZ-6) + 1
!
!     IF CURRENT PAGE IS PAST REQUESTED PAGE REWIND SCRATCH
!
      IF(JPAGE.LE.IPAGEY) THEN
         REWIND (UNIT=ISCRY)
         IPAGEY = 0
      END IF
!
!     SKIP UP TO PAGE PRECEDING ONE REQUIRED
!
   10 IPAGEY = IPAGEY + 1
      IF(IPAGEY.NE.JPAGE) THEN
         READ(ISCRY) ITOP
         GO TO 10
      END IF
!
!     LOAD CORRECT PAGE.
!
      READ(ISCRY) ITOP,(YP(K),K=1,ITOP)
!
!     SET LOWER AND UPPER INDICES
!
      ILOWY = (PAGESZ-6)*(IPAGEY-1)
      IHIGHY = ILOWY + PAGESZ
!
      RETURN
      END SUBROUTINE LOADY
!
!***********************************************************************
!
      SUBROUTINE WARNING_MESSAGE(JSEQ)
!
!     ROUTINE TO OUTPUT ERROR MESSAGE IN STANDARD FORM
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: JSEQ
!
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM
!
      INTEGER(KIND=I4) :: NEMES
!       Material printout delimiter
      CHARACTER(LEN=1) :: CHFF
!     CHFF=CHAR(12)                        ! Form-feed
      CHFF=CHAR(10)                        ! Double line-feed
!
!     PUT OUT ERROR MESSAGE
!
      IF(MATO.NE.MATERR) THEN
        WRITE(NOUT,'(A/1X,3A,I5)') CHFF
     &       ,'Check material ',ZSA,' MAT',MATO
        IF(NOUT.NE.IOUT)  THEN
           IF(IMDC.LT.4) WRITE(IOUT,'(/A)')  '   '
        END IF
      END IF
      IF(MATO.NE.MATERR.OR.MFO.NE.MFERR.OR.MTO.NE.MTERR) THEN
         WRITE(NOUT,'(//A,I4,A,I2,A,I3 )') '  WARNING(S)     IN MAT=',  &
     &          MATO,', MF=',MFO,', MT=',MTO
         MATERR = MATO
         MFERR = MFO
         MTERR = MTO
      END IF
      IF(JSEQ.NE.0) THEN
         WRITE(NOUT,'(5X,2A,I6)')  EMESS(1:49),'SEQUENCE NUMBER',JSEQ
!        Increment count for warnings
         NWARNG=NWARNG+1
      ELSE
         IF(EMESS.EQ.' ') THEN
            NEMES = 1
         ELSE
            NEMES = LEN_TRIM(EMESS)
         END IF
         WRITE(NOUT,'(5X,A)')  EMESS(1:NEMES)
      END IF
      RETURN
      END SUBROUTINE WARNING_MESSAGE
!
!***********************************************************************
!
      SUBROUTINE ERROR_MESSAGE(JSEQ)
!
!     ROUTINE TO OUTPUT ERROR MESSAGE IN STANDARD FORM
!     (If JSEQ record identifier is non-zero, the error counter
!      is incremented)
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: JSEQ
!
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM
!
      INTEGER(KIND=I4) :: NEMES
!       Material printout delimiter
      CHARACTER(LEN=1) :: CHFF
!     CHFF=CHAR(12)                        ! Form-feed
      CHFF=CHAR(10)                        ! Double line-feed
!
!     PUT OUT ERROR MESSAGE
!
      IF(MATO.NE.MATERR) THEN
        WRITE(NOUT,'(A/1X,3A,I5)') CHFF
     &       ,'Check material ',ZSA,' MAT',MATO
        IF(NOUT.NE.IOUT)  THEN
           IF(IMDC.LT.4) WRITE(IOUT,'(/A)')  '   '
        END IF
      END IF
!
      IF(MATO.NE.MATERR.OR.MFO.NE.MFERR.OR.MTO.NE.MTERR) THEN
         WRITE(NOUT,'(//A,I4,A,I2,A,I3 )') '  ERROR(S) FOUND IN MAT=',  &
     &          MATO,', MF=',MFO,', MT=',MTO
         MATERR = MATO
         MFERR = MFO
         MTERR = MTO
      END IF
      IF(JSEQ.NE.0) THEN
         WRITE(NOUT,'(5X,2A,I6)')  EMESS(1:49),'SEQUENCE NUMBER',JSEQ
!        Increment count for errors
         NERROR=NERROR+1
      ELSE
         IF(EMESS.EQ.' ') THEN
            NEMES = 1
         ELSE
            NEMES = LEN_TRIM(EMESS)
         END IF
         WRITE(NOUT,'(5X,A)')  EMESS(1:NEMES)
      END IF
      RETURN
      END SUBROUTINE ERROR_MESSAGE
!
!***********************************************************************
!
      SUBROUTINE TOKEN(INSTR,DELIM,ITOK,OUTSTR)
!
!     SUBROUTINE TO EXTRACT A STRING FROM A STRING WITH DELIMITERS
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: INSTR,OUTSTR,DELIM
      INTEGER(KIND=I4) :: ITOK
!
      INTEGER(KIND=I4), INTRINSIC :: INDEX, LEN_TRIM
!
      INTEGER(KIND=I4) :: ILEN,JLEN
      INTEGER(KIND=I4) :: ITOKP
      INTEGER(KIND=I4) :: IBEG
      INTEGER(KIND=I4) :: I
!
!     INITIALIZE
!
      OUTSTR = ' '
      ILEN = LEN_TRIM(INSTR)
      JLEN = LEN_TRIM(DELIM)
      IF(ITOK.EQ.0.OR.ILEN.EQ.0)   GO TO 100
      IF(JLEN.EQ.0)  THEN
         IF(ITOK.EQ.1)    OUTSTR = INSTR
         GO TO 100
      END IF
!
!     FIND ITOK-TH DELIMITER
!
      ITOKP = 1 - JLEN
      DO I=1,ITOK
         IBEG = ITOKP + JLEN
         IF(IBEG.LE.ILEN)  THEN
            ITOKP = INDEX(INSTR(IBEG:),DELIM) + IBEG - 1
            IF(ITOKP.LT.IBEG) ITOKP = ILEN + 1
         ELSE
            GO TO 100
         END IF
         IF(I.EQ.ITOK)  THEN
            IF(ITOKP.GT.IBEG) OUTSTR = INSTR(IBEG:ITOKP-1)
            GO TO 100
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE TOKEN
!
!***********************************************************************
!
      SUBROUTINE DATE(ADATE)
!
!     RETURNS DATE AS A CHARACTER STRING OF 11 CHARACTERS IN THE
!          FORM  DD-MMM-YYYY
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: ADATE
!
      CHARACTER(LEN=8) :: PDATE
      CHARACTER(LEN=4) :: YR
      CHARACTER(LEN=2) :: DD
      INTEGER(KIND=I4) :: MON
!
      CHARACTER(LEN=3), DIMENSION(12), PARAMETER ::                     &
     &         MONTHS = (/'Jan','Feb','Mar','Apr','May','Jun','Jul',    &
     &                    'Aug','Sep','Oct','Nov','Dec'/)
!
!     GET THE DATE AND TIME AS A CHARACTER STRING
!
      CALL DATE_AND_TIME(PDATE)
      READ(PDATE,'(A4,I2,A2)') YR,MON,DD
      ADATE = DD//'-'//MONTHS(MON)//'-'//YR
!
      RETURN
      END SUBROUTINE DATE
!
!***********************************************************************
!
      SUBROUTINE GET_FROM_CLINE
!
!     GET CONTENTS OF COMMAND LINE FOR VMS
!
      IMPLICIT NONE
!
!+++MDC+++
!...VMS
!/      INTEGER(KIND=2) :: ILENP2
!...ANS
!/      CHARACTER(LEN=100) :: CFILE
!---MDC---
!
      INPAR = ' '
      ILENP = 0
      NIN = INPUT0
!+++MDC+++
!...VMS
!/      CALL LIB$GET_FOREIGN(INPAR,,ILENP2)
!/      ILENP = ILENP2
!...UNX
!/      CALL GETCL(INPAR)
!/      ilenp = LEN_TRIM(INPAR)
!...DVF
      CALL GETARG(1,INPAR)
      ilenp = LEN_TRIM(INPAR)
!...ANS
!/      WRITE(IOUT,'(A)')                                               &
!/     &    ' Control File Specification        - '
!/      READ(NIN,'(A)') CFILE
!/      NIN = 19
!/      OPEN(UNIT=NIN,FILE=CFILE,STATUS='OLD')
!---MDC---
!
      RETURN
      END SUBROUTINE GET_FROM_CLINE
!
!***********************************************************************
!
      SUBROUTINE OUT_STATUS
!
!     DISPLAYS THE IDENTIFICATION OF THE SECTION BEING PROCESSED
!
      IMPLICIT NONE
!
      IF(MAT.GT.0.AND.MF.GT.0.AND.MT.GT.0) THEN
!+++MDC+++
!...VMS, ANS, WIN, UNX, MOD
         WRITE(IOUT,'(5X,A,I5,A,I3,A,I4)')                              &
     &         'PROCESSING MAT=',MAT,', MF=',MF,', MT=',MT
!---MDC---
      END IF
!
      RETURN
      END SUBROUTINE OUT_STATUS
!
!+++MDC+++
!...VMS, ANS, WIN, UNX
      END PROGRAM FIZCON
!...LWI, DVF, MOD
!/      END MODULE FIZCON
!---MDC---
