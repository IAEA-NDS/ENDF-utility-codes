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
!+++MDC+++
!...VMS, ANS, WIN, UNX
!
!     Main program for non-windows implementation of PSYCHE
!
      PROGRAM PSYCHE
!
      IMPLICIT NONE
!
!...LWI, DVF, MOD
!/!
!/!     MODULE IMPLEMENTATION OF PSYCHE FOR MODLIB AND WINDOWS
!/!
!/      MODULE PSYCHE
!/!
!/      IMPLICIT NONE
!/!
!/      PRIVATE
!/!
!/      PUBLIC :: RUN_PSYCHE
!/      PUBLIC :: PSYCHE_INPUT, PSYCHE_DATA, PSYCHE_SUCCESS
!---MDC---
!-T Program PSYCHE
!-P Perform physics tests on data in evaluated nuclear data files
!-P in ENDF-5 or ENDF-6 format
!-V
!-V         Version 8.11   December 2015     A. Trkov
!-V                        Increase arrays to process MF6/MT5
!-V         Version 8.10   January 2014      A. Trkov
!-V                        Correct recoil energy in energy balance for LCT=3
!-V         Version 8.09   September 2013    D. Brown
!-V                        Fix value test of LNU when no fission data given
!-V         Version 8.08   October2012       A. Trkov
!-V                        Allow E-dependent scattering radius in URR
!-V         Version 8.07   September 2012    A. Koning
!-V                        Cleanup unused variables.
!-V         Version 8.06   August 2012       A. Trkov
!-V                        Fix initialisation of LNU
!-V         Version 8.05   February 2011     A. Trkov
!-V                        Fix energy balance calculation depending on
!-V                        on the value of the LCT flag in MF6.
!-v                        LCT=3 is treated as CM, assuming that other 
!-V                        particles fly off in different directions,
!-V                        leaving the residual approximately in the CM.
!-V         Version 8.04   February 2011     A. Trkov
!-V                        Add sign to energy balance printout
!-V         Version 8.03   January  2011     A. Trkov
!-V                        Correct the format for printing energy balance
!-V                        of MT5.
!-V         Version 8.02   January  2011     A. Trkov
!-V                        Implement resolved resonance option LRF=7
!-V                        (no extensive checking of parameters is done).
!-V         Version 8.01   October  2010     A. Trkov
!-V                        Cosmetic changes.
!-V         Version 8.00   August  2008     A. Trkov
!-V                        1. Major updating of the code.
!-V                        2. Further reduction of non-essential output.
!-V                        3. Read real variables in double precision to
!-V                           avoid reporting errors reading small numbers.
!-V                        4. Implement extended features of the format
!-V                           endorsed by CSEWG.
!-V         VERSION 7.02   MAY 2005     C.L.DUNFORD
!-V                        1. Fixed bug in calculation of L=2 penetrability
!-V         VERSION 7.01   APRIL 2005     C.L.DUNFORD
!-V                        1. SET SUCCESS FLAG AFTER RETURN FROM BEGIN
!-V                        2. ADD POTENTIAL SCATTERING TEST FORMERLY IN
!-V                           FIZCON
!-V         Version 7.0    OCTOBER 2004     C.L.DUNFORD
!-V                        1. MODIFIED TO PROVIDE A MODULE FOR THE NEA
!-V                           MODLIB PROJECT
!-V                        2. ALLOW ENERGY DEPENDENT DELAYED FISSION
!-V                           GROUP PARAMETERS.
!-V                        4. PERMIT USER TO SUPPLY BATCH INPUT FILE
!-V                             NAME
!-V                        5. REMOVED FORTRAN LINE CONTROLS FROM OUTPUT
!-V                        6. ADDED COMMAND LINE INPUT TO UNIX AND
!-V                           WINDOWS VERSIONS. NOTE: ONLY INPUT AND
!-V                           OUTPUT FILE NAMES CAN BE GIVEN. DEFAULT
!-V                           OPTIONS ARE ASSUMED UNLESS THIRD
!-V                           PARAMETER IS N.
!-V
!-V      REFER ALL COMMENTS AND INQUIRIES TO
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
!-M PSYCHE - Perform Physics Tests on Evaluated Nuclear Data Files
!-M ==============================================================
!-M
!-M PSYCHE is a program for checking the physics content of an 
!-M evaluated data file. It can recognize the difference between 
!-M ENDF-6 and ENDF-5 formats and performs its tests accordingly. 
!-M The present version checks for energy conservation for emitted 
!-M particles, checks Wick's limit for elastic scattering, analyzes 
!-M resonance parameter statistics, calculates thermal cross sections 
!-M and resonance integrals, examines continuity across resonance 
!-M region boundaries and checks "Q" values against mass tables. 
!-M It is assumed the file being checked has passed the CHECKR 
!-M program without any errors being detected.
!-M 
!-M Fortran Logical Units Used:
!-M     5  Default (keyboard) input 
!-M     6  Default output (terminal screen)
!-M    20  Input data file, ENDF format 
!-M    21  Message file for program checking results  
!-M 22-25  Temporary files for energy conservation tests 
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
!-M      4  Options selection (2 fields)
!-M         - Material number where processing starts (integer)
!-M           (If zero or blank, then checking begins with the 
!-M           first material)
!-M         - Material number where processing ends (integer)
!-M           (If zero or blank, then checking continues to end 
!-M			      of the file)
!-M         If Record 4 is left entirely blank, then the entire input
!-M         file is processed.
!-M
!-M Multiple input files can be processed to produce multiple output 
!-M files by repeating the above input data sequence. The program 
!-M execution is terminated with a record containing the word DONE.
!-M 
!-M In interactive mode operation, the above data is supplied in response 
!-M to the appropriate query; in graphical mode, via a dialog box.
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
!     PSYCHE VERSION NUMBER
!
      CHARACTER(LEN=*), PARAMETER :: VERSION = '8.11'
!
!     DEFINE VARIABLE PRECISION
!
      INTEGER(KIND=4), PARAMETER :: I4 = SELECTED_INT_KIND(8)
      INTEGER(KIND=4), PARAMETER :: R4 = SELECTED_REAL_KIND(6,37)
      INTEGER(KIND=4), PARAMETER :: R8 = SELECTED_REAL_KIND(15,307)
!
!     Half of maximum number of Legendre coefficients
      INTEGER(KIND=I4), PARAMETER :: MXLG=32
!
      REAL(KIND=R4), PARAMETER :: FACTOR=1.008665
      REAL(KIND=R4), PARAMETER :: OTHIRD=1./3.
!
!     STANDARD FORTRAN INPUT AND OUTPUT UNITS
!
      INTEGER(KIND=I4), PARAMETER :: IOUT=6, INPUT0=5
      INTEGER(KIND=I4), PARAMETER :: ISCR=22,JSCR=23,KSCR=24,LSCR=25
      INTEGER(KIND=I4) :: NIN
!
!     ENDF DISK FILE INPUT AND CHECKING OUTPUT FORTRAN UNITS
!
      INTEGER(KIND=I4), PARAMETER :: JIN=20,JOUT=21
!
!     INCIDENT PARTICLE ZA DESIGNATION
!
      INTEGER(KIND=I4) :: IZAIN
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
!/      INTEGER(KIND=I4), :: PARAMETER :: IMDC = 0
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
      TYPE PSYCHE_INPUT
         CHARACTER(LEN=100) :: INFIL
         CHARACTER(LEN=100) :: OUTFIL
         INTEGER(KIND=I4) :: MATMIN
         INTEGER(KIND=I4) :: MATMAX
      END TYPE PSYCHE_INPUT
!
      TYPE(PSYCHE_INPUT) PSYCHE_DATA
!
!     FLAG TO INDICATE WHETHER MULTIPLE INPUT FILES CAN BE SELECTED
!
      INTEGER(KIND=I4) :: IONEPASS        !  0, YES;  1, NO
!
!     FLAG TO INDICATE SUCCESS OR FAILURE OF STANEF EXECUTION
!
      INTEGER(KIND=I4) :: PSYCHE_SUCCESS
!
!     FILE (TAPE) LABEL FROM FIRST RECORD
!
      CHARACTER(LEN=66) :: TLABEL
      INTEGER(KIND=I4) :: LABEL
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
!     1000*Z + A OF MATERIAL CURRENTLY BEING PROCESSED
!        AWR IS THE RATIO OF THE MATERIAL MASS TO THAT OF THE NEUTRON
!        AWI IS THE RATIO OF THE PROJECTILE MASS TO THE THAT OF NEUTRON
!        STA = 0.0, STABLE MATERIAL; STA = 1.0 RADIOACTIVE MATERIAL
!        ELIS IS THE EXCITATION ENERGY OF THE TARGET NUCLEUS
!
      REAL(KIND=R4) :: ZA,AWR,AWI,STA,ELIS
!
!     LIS   IS THE STATE NUMBER (0 FOR GROUND) OF THE MATERIAL
!     LISO  IS THE ISOMER STATE NUMBER OF THE MATERIAL
!
      INTEGER(KIND=I4) :: LIS,LISO
!
!     LRP    IS THE RESONANCE PARAMETER FLAG
!     LFI    IS THE FISSION FLAG
!
      INTEGER(KIND=I4) :: LRP,LFI
!
!     CONTENTS OF FIRST RECORD OF A LIST RECORD
!
      INTEGER(KIND=I4) :: L1L,L2L,NPL,N2L
      REAL(KIND=R4)    :: C1L,C2L
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
!     CONTENTS OF FIELDS ON A HEAD/CONT RECORD
!
      INTEGER(KIND=I4) :: L1H,L2H,N1H,N2H
      REAL(KIND=R4)    :: C1H,C2H
!
!     MAXIMUM SIZE OF AN INTERPOLATION TABLE
!
      INTEGER(KIND=I4), PARAMETER :: INTABMAX=20
      INTEGER(KIND=I4), DIMENSION(INTABMAX) :: NTERP,INTERP
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
!     DATA TABLE
!
      REAL(KIND=R4), DIMENSION(50000) :: X,Y
!
      INTEGER(KIND=I4) :: MAT,MF,MT,NSEQ
!
!     END OF FILE FLAG
!
      INTEGER(KIND=I4) :: IFIN
!
!     ENERGY LIMITS FOR THE MATERIAL
!
      REAL(KIND=R4), PARAMETER :: ENMIN = 1.0E-05
      REAL(KIND=R4) :: ENMAX
!
!     NUBAR MULTIPLICITIES EXPANSION ARRAYS
!
      INTEGER(KIND=I4) :: LNU,NRC,NPT
      INTEGER(KIND=I4), DIMENSION(INTABMAX) :: ITERP,JTERP
      REAL(KIND=R4), DIMENSION(4) :: POLI
      REAL(KIND=R4), DIMENSION(5000) :: ETERP,ANUTRP
!
!     ISOTOPES IN RESONANCE REGION
!
      INTEGER(KIND=I4) :: NIS,NISN,NISMAT
      INTEGER(KIND=I4), PARAMETER :: NISMAX=20
      INTEGER(KIND=I4), DIMENSION(NISMAX) :: NERS
      REAL(KIND=R4), DIMENSION(NISMAX) :: ZAIS,ABNS,SPINS
      INTEGER(KIND=I4), DIMENSION(NISMAX) :: IAIS
      REAL(KIND=R4), DIMENSION(NISMAX) :: ZAIMAT,ABNIS,SPINIS
!
!     ENERGY REGION BOUNDARIES
!
      INTEGER(KIND=I4) :: IUNR, NBOUND
      REAL(KIND=R4) :: E1,E2
      INTEGER(KIND=I4), PARAMETER :: NREGMAX=14
      INTEGER(KIND=I4), DIMENSION(NREGMAX) :: IRES
      REAL(KIND=R4), DIMENSION(NREGMAX) :: EMIDLE
      INTEGER(KIND=I4), PARAMETER :: NBOUNDMAX=NREGMAX+3
      REAL(KIND=R4), DIMENSION(NBOUNDMAX) :: RBOUND
!
!     ENERGY DEPENDENT SCATTERING RADIUS EXPANSION ARRAYS
!
      INTEGER(KIND=I4), DIMENSION(NREGMAX) :: NRED,NPED
      INTEGER(KIND=I4), DIMENSION(INTABMAX,NREGMAX) ::  NBTED=0,JNTED=0
      REAL(KIND=R4), DIMENSION(100,NREGMAX) :: EP=0.,APED=0.
!
!     POTENTIAL SCATTERING RADIUS LIMITS
!
      INTEGER(KIND=4) :: NRO
      REAL(KIND=I4) :: APLO,APHI
!
!     RESONANCE PARAMETER STORAGE
!
      INTEGER(KIND=I4) :: IPT
      REAL(KIND=R4), DIMENSION(40000) :: RESPAR
      REAL(KIND=R4), DIMENSION(500) :: CONSA,SFA
!
!     ARRAY STORING Q-VALUES FROM FILE 3 FOR LATER TESTS
!
      INTEGER(KIND=I4) :: NMT3, NMTBAR
      INTEGER(KIND=I4), PARAMETER ::  SZMT3=150
      INTEGER(KIND=I4), DIMENSION(SZMT3) :: MT3
      REAL(KIND=R4), DIMENSION(SZMT3) :: QIVAL,QVAL
      REAL(KIND=R4), PARAMETER :: QUNK= 7.777E+07
!
!     LEGENDRE COEFFICIENTS
!
      INTEGER(KIND=I4) :: NCOEF, NMOM,NPOWS
      REAL(KIND=R4), DIMENSION(MXLG+1) :: MEAN
      REAL(KIND=R4), DIMENSION(MXLG) :: VAR
      INTEGER(KIND=I4), PARAMETER :: NCOEFMAX=64
      REAL(KIND=R4), DIMENSION(NCOEFMAX) :: F
      REAL(KIND=R4), DIMENSION(NCOEFMAX) :: MOMENT=0.
      INTEGER(KIND=I4), PARAMETER :: NCOEFMAX1=NCOEFMAX+1
      REAL(KIND=R4), DIMENSION(NCOEFMAX1,NCOEFMAX1) :: W
!
!     SECTIONS INVOLVED IN ENERGY BALANCE CALCULATIONS
!
      INTEGER(KIND=I4) :: MT1213
      INTEGER(KIND=I4), DIMENSION(10) :: MTSCR
      REAL(KIND=R4) :: ENERGY
      INTEGER(KIND=I4), DIMENSION(50) :: MTBAR
      INTEGER(KIND=I4), PARAMETER :: NSECMAX=1000
      INTEGER(KIND=I4) :: MTTOT  ! SECTIONS WITH PHOTONS
      INTEGER(KIND=I4), DIMENSION(NSECMAX) :: MTDCT
      INTEGER(KIND=I4), PARAMETER :: NSMF5=NSECMAX
      INTEGER(KIND=I4) :: MF5   ! SECTIONS IN MF5 WITH PHOTONS
      INTEGER(KIND=I4), DIMENSION(NSMF5) :: MT5
      INTEGER(KIND=I4), PARAMETER :: NBSMAX=50
      INTEGER(KIND=I4) :: NBS,NBSN  ! SECTIONS IN MF13
      INTEGER(KIND=I4), DIMENSION(NBSMAX) :: MTS
      INTEGER(KIND=I4), PARAMETER :: NGRIDMAX=10000
      REAL(KIND=R4), DIMENSION(NGRIDMAX) :: EGRID,AEBAR,EINT,PKINT
      REAL(KIND=R4), DIMENSION(NGRIDMAX) :: EBAV,YLD
!
      REAL(KIND=R4), PARAMETER :: BIGNO=1.0E+20
      REAL(KIND=R4), PARAMETER :: PI=3.14159265
      REAL(KIND=R4), PARAMETER :: CROC=2.196807E-03
      REAL(KIND=R4), PARAMETER :: ANEUTR=1.008664904
      REAL(KIND=R4), PARAMETER :: THERMA=0.0253
      REAL(KIND=R4), PARAMETER :: EPI6=1.0E-06
      REAL(KIND=R4), PARAMETER :: EPI4=1.0E-04
      REAL(KIND=R4), PARAMETER :: EPI3=1.0E-03
      REAL(KIND=R4), PARAMETER :: EPI2=1.0E-02
      REAL(KIND=R4), PARAMETER :: EPI1=1.0E-01
!
      INTEGER(KIND=I4), PARAMETER :: NSMASS=23
      INTEGER(KIND=I4), DIMENSION(NSMASS), PARAMETER ::                  &      
     &       IZAB =(/1,1001,1002,1003,2003,2004,3006,  3007,  4009,     &       
     &               5010,  7014, 20048, 22046, 22049, 24054, 56136,    &       
     &               56137, 57139, 60142, 66156, 72174, 72180, 82208/)
      REAL(KIND=R4), DIMENSION(NSMASS) ::                                &      
     &       AMASSB=(/8071.43,7289.03,13135.84,14949.94,14931.32,       &       
     &               2424.93,14087.28,14908.22,11348.02,12051.77,       &       
     &                2863.44,-44216.52,-44122.72,-45329.01,-56931.33,  &       
     &                -87870.41,-87732.55,-87230.47,-85948.72,-70526.78,&       
     &                -55829.55,-49778.51,-21758.95/)
!
!     ERROR MESSAGE TEXT
!
      CHARACTER(LEN=80) :: EMESS
!
!     ERROR FLAG
!
      INTEGER(KIND=I4) :: IERX
!
!***********************************************************************
!
!+++MDC+++
!...VMS, ANS, WIN, UNX
!
      CALL RUN_PSYCHE
!
!     TERMINATE JOB
!
      IF(PSYCHE_SUCCESS.EQ.0) THEN
         WRITE(IOUT,'(/A)') ' '
         STOP '    PSYCHE - Tests completed successfully'
      ELSE
         WRITE(IOUT,'(/A)') ' '
         STOP '    PSYCHE - Tests terminated abnormally!'
      END IF
!---MDC---
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_PSYCHE
!
      IMPLICIT NONE
!
!
      CHARACTER(LEN=*), PARAMETER :: DASHES = REPEAT('-',80)
!
      CHARACTER(LEN=80) :: IFIELD
      INTEGER(KIND=I4) :: IQUIT   ! FLAG TO INDICATE WHETHER OR NOT TO EXIT     
      INTEGER(KIND=I4) :: IFIND   ! FLAGS WHETHER DESIRED MATERIAL FOUND
      INTEGER(KIND=I4) :: I
!
!     OUTPUT PROGRAM IDENTIFICATION
!
      PSYCHE_SUCCESS = 0
      IF(IMDC.LT.4) THEN
         WRITE(IOUT,'(/2A)')' PROGRAM PSYCHE VERSION ',VERSION
      END IF
!
!     CHECK FOR COMMAND LINE OPERATION
!
      IONEPASS = 0
      CALL GET_FROM_CLINE
!
!     INITIALIZE RUN
!
   10 CALL BEGIN(IQUIT)
      IF(IQUIT.GT.0) THEN
         IF(IONEPASS.EQ.1) PSYCHE_SUCCESS = 1
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
            WRITE(IOUT,'(//5X,2A)')  'PSYCHE ERROR - END OF FILE ',     &       
     &                      'ENCOUNTERED BEFORE TEND RECORD FOUND!'
         END IF
         IF(NOUT.NE.IOUT)   THEN
            WRITE(IOUT,'(//5X,2A)')  'PSYCHE ERROR - END OF FILE ',     &       
     &                      'ENCOUNTERED BEFORE TEND RECORD FOUND!'
         END IF
         IF(NOUT.NE.IOUT) THEN
           WRITE(NOUT,'(A)') ' Done PSYCHE'
           CLOSE(UNIT=NOUT)
         END IF
         CLOSE(UNIT=JIN)
         PSYCHE_SUCCESS = 1
         GO TO 100
      END IF
!
!     SET UP TO PROCESS NEXT SECTION
!
      IF(MAT.NE.MATO) THEN
         IF(PSYCHE_DATA%MATMAX.NE.0.AND.MAT.GT.PSYCHE_DATA%MATMAX)      &       
     &            GO TO 60
         NSEQP1 = NSEQP
         MATO = MAT
         NMT3 = 0
         MTS(1) = 1
         MTS(2) = 2
         NBS = 2
         NBSN = 0
         MT1213 = 0
         NMTBAR = 0
         MTTOT = 0
         MF5 = 0
         LNU  = 0
         REWIND (UNIT=ISCR)
         REWIND (UNIT=JSCR)
         REWIND (UNIT=KSCR)
         REWIND (UNIT=LSCR)
         WRITE(NOUT,'(A/3X,A,I5)')  CHAR(12),'CHECK MATERIAL',MAT
         WRITE(NOUT,'(19X,A)')                                          &       
     &          '(NO WARNINGS DETECTED IN SECTIONS WITHOUT COMMENTS)'
         IF(NOUT.NE.IOUT) THEN
            IF(IMDC.LT.4) WRITE(IOUT,'(/A)')   '   '
         END IF
      END IF
      IF(MF.NE.MFO)   THEN
         WRITE(NOUT,'(/A/A,I2)')  DASHES,'FILE ',MF
         MFO = MF
      END IF
!
!     NEW SECTION
!
   30 WRITE(NOUT,'(3X,A,I3)')  'SECTION ',MT
      MTO = MT
!
!     IN INTERACTIVE MODE OUTPUT CURRENT SECTION ID TO TERMINAL
!
      IF(NOUT.NE.IOUT)  THEN
         IF(IMDC.LT.4) THEN
            WRITE(IOUT,'(5X,A,I5,A,I3,A,I4)')                           &       
     &        'PROCESSING MAT=',MATO,', MF=',MFO,', MT=',MTO
         ELSE
            CALL OUT_STATUS
         END IF
      END IF
!
!     CHECK THE NEW SECTION
!
      CALL CHKSEC
      IF(IERX.EQ.2) GO TO 20
!
!     IF FATAL ERROR FOUND OR SECTION NOT PROCESSED, SKIP REST OF
!        SECTION
!
   35 IF(IERX.NE.0) THEN
         IERX = 2
         DO WHILE (MT.NE.0)
            READ(JIN,'(A)',END=20) IFIELD
            READ(IFIELD,'(66X,I4,I2,I3,I5)',ERR=40)  MAT,MF,MT,NSEQ
   40       CONTINUE
         END DO
         IERX = 0
      END IF
!
!     READ UNTIL HEAD OR TEND RECORD FOUND
!
   50 IF(MAT.GE.0)  THEN
   55    CALL RDHEAD(I)
         IF(IERX.EQ.1) GO TO 35
         IF(I.GT.1.AND.I.LT.5)   THEN
            GO TO 55
         ELSE IF(I.EQ.5) THEN
            IFIN = 1
         END IF
      ELSE
         GO TO 100
      END IF
!
!     END OF MATERIAL TESTS
!
      IF(MAT.NE.MATO.OR.IFIN.EQ.1)   THEN
         IF(MTTOT.GT.0)   CALL CKNG
         IF(LRP.EQ.1)  REWIND (UNIT=ISCR)
         REWIND (UNIT=JSCR)
         GO TO 57
      END IF
!
!     END OF FILE TESTS
!
      IF(MF.NE.MFO.OR.IFIN.EQ.1) REWIND (UNIT=JSCR)
!
!     CHECK END OF TAPE FLAG
!
   57 IF(IFIN.EQ.0) THEN
        IF(PSYCHE_DATA%MATMAX.EQ.0.OR.MAT.LE.PSYCHE_DATA%MATMAX)        &       
     &                 GO TO 20
      END IF
!
!     CLOSE FILES
!
   60 IF(NOUT.NE.IOUT) THEN
         WRITE(NOUT,'(A)') ' Done PSYCHE'
         CLOSE(UNIT=NOUT)
      END IF
      CLOSE(UNIT=JIN)
      CLOSE(UNIT=ISCR,STATUS='DELETE')
      CLOSE(UNIT=JSCR,STATUS='DELETE')
      CLOSE(UNIT=KSCR,STATUS='DELETE')
      CLOSE(UNIT=LSCR,STATUS='DELETE')
!
!     SEE IF ONE PASS LIMIT SET
!
      IF(IONEPASS.EQ.0) GO TO 10
!
  100 RETURN
      END SUBROUTINE RUN_PSYCHE
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
      CHARACTER(LEN=1) :: IW
      CHARACTER(LEN=15) :: MATSIN
      CHARACTER(LEN=11) :: ADATE
      LOGICAL(KIND=I4) :: IEXIST
      INTEGER(KIND=I4) :: IC
!
!     INITIALIZE PROCESSING CONTROL VARIABLES
!
      IERX = 0
      MATO = 0
      MFO = 0
      MTO = 0
      IFIN = 0
      NOUT = IOUT
   10 IQUIT = 0
!
!     INITIALIZE TO STANDARD OPTIONS
!
      IF(IMDC.LT.4) THEN
         PSYCHE_DATA%INFIL = '*'
         PSYCHE_DATA%OUTFIL = '*'
         PSYCHE_DATA%MATMIN = 0
         PSYCHE_DATA%MATMAX = 0
      END IF
      SELECT CASE (IMDC)
         CASE (0)
            IW = 'N'
            IONEPASS = 0
         CASE(1,2,3)
            IF(ILENP.NE.0)  THEN
               CALL TOKEN(INPAR,'%',1,PSYCHE_DATA%INFIL)
               CALL TOKEN(INPAR,'%',2,PSYCHE_DATA%OUTFIL)
               CALL TOKEN(INPAR,'%',3,IW)
               IC = ICHAR(IW)
               IF(IC.GT.96.AND.IC.LT.123)  IW = CHAR(IC-32)
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
         IF(PSYCHE_DATA%INFIL.EQ.'*') THEN
            IF(IMDC.NE.0) THEN
               WRITE(IOUT,FMT=TFMT)                                     &       
     &             ' Input File Specification             - '
            END IF
            READ(NIN,'(A)') PSYCHE_DATA%INFIL
         ELSE
            WRITE(IOUT,'(/2A)') ' Input file - ',                       &       
     &                 TRIM(PSYCHE_DATA%INFIL)
         END IF
      END IF
!
!     SEE IF INPUT INDICATES FILE TERMINATION
!
      IF(PSYCHE_DATA%INFIL.EQ.' '.OR.PSYCHE_DATA%INFIL.EQ.'DONE') THEN
         IQUIT = 1
         GO TO 100
      END IF
!
!     MAKE SURE INPUT FILE EXISTS
!
      INQUIRE(FILE=PSYCHE_DATA%INFIL,EXIST=IEXIST)
      IF(.NOT.IEXIST)  THEN
         IF(IMDC.LT.4) THEN
            WRITE(IOUT,'(/7X,A/)')  'COULD NOT FIND INPUT FILE'
         END IF
         SELECT CASE (IMDC)
            CASE (1,2,3)
               IF(IONEPASS.EQ.0) GO TO 10
         END SELECT
         IQUIT = 1
         PSYCHE_SUCCESS = 1
         GO TO 100
      END IF
!
!     GET OUTPUT FILE SPECIFICATION
!
      IF(IMDC.LT.4) THEN
         IF(PSYCHE_DATA%OUTFIL.EQ.'*' ) THEN
            IF(IMDC.NE.0) THEN
               WRITE(IOUT,FMT=TFMT)                                     &       
     &           ' Output Message File Specification    - '
            END IF
            READ(NIN,'(A)') PSYCHE_DATA%OUTFIL
         ELSE
            WRITE(IOUT,'(/2A)') ' Output file - ',                      &       
     &                   TRIM(PSYCHE_DATA%OUTFIL)
         END IF
      END IF
      IF(PSYCHE_DATA%OUTFIL.NE.' ')  THEN
         NOUT = JOUT             ! SETS FORTRAN OUTPUT UNIT IF DISK FILE
      END IF
!
!     CHECK IF ENTIRE TAPE TO BE PROCESSED (INTERACTIVE MODE ONLY)
!
      IF(IMDC.NE.0) THEN
         IF(IW.EQ.'*') THEN
            IF(IMDC.LT.4) THEN
               WRITE(IOUT,FMT=TFMT)                                     &       
     &                ' Check Entire File (Y(es),N(o))?  '
               READ(NIN,'(A)')  IW
               IC = ICHAR(IW)
               IF(IC.GT.96.AND.IC.LT.123)  IW = CHAR(IC-32)
            END IF
         END IF
      END IF
!
!     GET MATERIAL NUMBER RANGE (ALL) IF DEFAULT NOT SELECTED
!
      IF(imdc.eq.0.or.(IW.EQ.'N'.AND.IMDC.LT.4)) THEN
         CALL SELECT_MATS(MATSIN)
      END IF
!
!     OPEN INPUT AND OUTPUT FILES
!
      OPEN(UNIT=JIN,ACCESS='SEQUENTIAL',STATUS='OLD',                   &       
     &               FILE=PSYCHE_DATA%INFIL,ACTION='READ')
      IF(NOUT.NE.6) THEN
!+++MDC+++
!...VMS
!/         OPEN(UNIT=NOUT,ACCESS='SEQUENTIAL',STATUS=OSTATUS,           &       
!/     &       FILE=PSYCHE_DATA%OUTFIL,CARRIAGECONTROL='LIST')
!...WIN, DVF, UNX, LWI, ANS, MOD
         OPEN(UNIT=NOUT,ACCESS='SEQUENTIAL',STATUS=OSTATUS,             &       
     &       FILE=PSYCHE_DATA%OUTFIL)
!---MDC---
      END IF
!
!     OPEN SCRATCH FILES
!
      OPEN(UNIT=ISCR,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(UNIT=JSCR,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(UNIT=KSCR,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(UNIT=LSCR,FORM='UNFORMATTED',STATUS='SCRATCH')
!
!     OUTPUT SELECTED OPTIONS
!
      CALL DATE(ADATE)
      IF(IMDC.LT.4)  WRITE(IOUT,'(/A)') ' '
      IF(NOUT.NE.IOUT) THEN
         WRITE(NOUT,'(A///2A,30X,2A/)') CHAR(12),                       &       
     &            'PROGRAM PSYCHE VERSION ',VERSION,                    &       
     &            'Run on ',ADATE
      END IF
      WRITE(NOUT,'(2A)')                                                &       
     &   'Input File Specification------------------------',            &       
     &   TRIM(PSYCHE_DATA%INFIL)
      IF(PSYCHE_DATA%MATMIN.EQ.0.AND.PSYCHE_DATA%MATMAX.EQ.0)   THEN
         WRITE(NOUT,'(A)')  'Check the Entire File'
      ELSE
         WRITE(NOUT,'(A,I4,A,I4)')                                      &       
     &        'Check Materials---------------------------------',       &       
     &             PSYCHE_DATA%MATMIN,' to ',PSYCHE_DATA%MATMAX
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
      INTEGER(KIND=I4), INTRINSIC :: INDEX, LEN_TRIM
 
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
         PSYCHE_DATA%MATMIN = 0
         PSYCHE_DATA%MATMAX = 0
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
      PSYCHE_DATA%MATMIN = 1
      PSYCHE_DATA%MATMAX = 9999
      READ(BUF1,'(BN,I4)',ERR=20) PSYCHE_DATA%MATMIN
   20 READ(BUF2,'(BN,I4)',ERR=25) PSYCHE_DATA%MATMAX
!
!     SET THE MATERIAL NUMBER LIMITS
!
   25 IF(PSYCHE_DATA%MATMIN.LE.0) THEN
         PSYCHE_DATA%MATMIN = 1
      END IF
      IF(PSYCHE_DATA%MATMAX.LT.PSYCHE_DATA%MATMIN)  THEN
         PSYCHE_DATA%MATMAX = PSYCHE_DATA%MATMIN
      END IF
      IF(PSYCHE_DATA%MATMIN.EQ.1.AND.PSYCHE_DATA%MATMAX.EQ.9999) THEN
         PSYCHE_DATA%MATMIN = 0
         PSYCHE_DATA%MATMAX = 0
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
         WRITE(NOUT,'(/2A,I5/4X,2A)') 'TAPE BEING PROCESSED IS ',       &       
     &         'NUMBERED',LABEL,'LABEL IS  ',TLABEL
      END IF
      GO TO 40
!
!     IF READING ERROR ASSUME A PROPER LABEL AND GO ON
!
   20 WRITE (NOUT,'(5X,A//)')                                           &       
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
   60 IF(PSYCHE_DATA%MATMIN.GT.0)   THEN
         DO WHILE(MAT.LT.PSYCHE_DATA%MATMIN)
            READ(JIN,'(A)',END=90)  IFIELD
            READ(IFIELD,'(66X,I4,I2,I3,I5)',ERR=65) MAT,MF,MT,NSEQ
   65       IF(MAT.LT.0) GO TO 70
         END DO
         IF(MAT.GT.PSYCHE_DATA%MATMAX) GO TO 70
      END IF
      GO TO 75
!
!     FAILED TO FIND A MATERIAL
!
   70 IF(PSYCHE_DATA%MATMIN.EQ.PSYCHE_DATA%MATMAX) THEN
         IF(PSYCHE_DATA%MATMIN.EQ.0) THEN
            EMESS = 'INPUT FILE DOES NOT CONTAIN ANY ENDF EVALUATIONS'
         ELSE
            WRITE(EMESS,'(A,I5)')                                       &       
     &           'INPUT FILE DOES NOT CONTAIN MATERIAL',                &       
     &                   PSYCHE_DATA%MATMIN
         END IF
      ELSE
         WRITE(EMESS,'(A,I5,A,I5)')                                     &       
     &        'INPUT FILE DOES NOT CONTAIN ANY MATERIALS',              &       
     &         PSYCHE_DATA%MATMIN,' TO',PSYCHE_DATA%MATMAX
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
      IF(MF.EQ.1) THEN
         CALL CKF1
      ELSE IF(MF.EQ.2) THEN
         CALL CKF2
      ELSE IF(MF.EQ.3) THEN
         CALL CKF3
      ELSE IF(MF.EQ.4) THEN
         CALL CKF4
      ELSE IF(MF.EQ.5) THEN
         CALL CKF5
      ELSE IF(MF.EQ.6) THEN
         CALL CKF6
      ELSE IF(MF.EQ.12) THEN
         CALL CKF12
      ELSE IF(MF.EQ.13) THEN
         CALL CKF13
      ELSE IF(MF.EQ.14) THEN
         CALL CKF14
      ELSE IF(MF.EQ.15) THEN
         CALL CKF15
      ELSE
         IERX = 1
      END IF
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
      INTEGER(KIND=I4) :: NCD,NXC
      INTEGER(KIND=I4) :: NID
      INTEGER(KIND=I4) :: NC,J,N,NN
!
      IERX = 0
!
!     COMMENTS AND DIRECTORY
!
      IF(MT.EQ.451) THEN
         ZA = C1H
         AWR = C2H
         LRP = L1H
         LFI = L2H
         NLIB = N1H
         NMOD = N2H
!********READ THE NEXT CONTROL RECORD AND SET PARAMETERS
         CALL RDCONT
         ELIS  = C1H
         STA = C2H
         LIS = L1H
         LISO = L2H
         NFOR = N2H
!********ENDF-V FORMAT FILE
         IF(NFOR.EQ.0)  THEN
            NFOR = 5
            IF(NLIB.GE.2.AND.NLIB.LE.4)  THEN
               NVER = 1
            ELSE IF(NLIB.EQ.5)  THEN
               NVER = 2
            ELSE IF(NLIB.EQ.6)  THEN
               NVER = 3
            ELSE IF(NLIB.EQ.35)  THEN
               NVER = 1
            ELSE
               NVER = 5
            END IF
            ENMAX = 2.0E+7
            NSUB = 10
            AWI = 1.
!********ENDF-VI OR LATER FORMAT, READ ANOTHER CONTROL RECORD
         ELSE
            CALL RDCONT
            AWI = C1H
            ENMAX = C2H
            NSUB = N1H
            NVER = N2H
            IF(NFOR.GT.6)   NFOR = 6
         END IF
!********PROCESS LAST CONTROL RECORD
         CALL RDCONT
!********READ IN COMMENT RECORDS
         IF(NFOR.GE.5) THEN
            NID = 5
         ELSE
            NID = 2
         END IF
         NCD = N1H
         DO NC=1,NCD
            CALL RDTEXT
            IF(NC.LE.NID)   THEN
               IF(IMDC.LT.4) WRITE(IOUT,'(1X,A)')   TEXT
               IF(IOUT.NE.NOUT)   WRITE(NOUT,'(5X,A)')   TEXT
            END IF
         END DO
!********PROCESS DIRECTORY
         NXC = N2H
         DO N=1,NXC
            CALL RDCONT
!***********LOOK THROUGH DICTIONARY FOR PHOTON REACTIONS
            IF(L1H.EQ.14)   THEN
               IF(MTTOT.GE.NSECMAX)  THEN
                  WRITE(EMESS,'(A,I4,A)')                               &       
     &                'CANNOT CHECK AVERAGE PHOTON ENERGY FOR  MT = ',  &       
     &                   L2H,'. TOO MANY MTS.'
                  CALL ERROR_MESSAGE(0)
               ELSE
                  MTTOT = MTTOT + 1
                  MTDCT(MTTOT) = L2H
               END IF
            ELSE IF(L1H.EQ.13)   THEN
               IF(NBS.GE.NBSMAX)   THEN
                  WRITE(EMESS,'(A,I4,A)')                               &       
     &                'CANNOT CHECK AVERAGE PHOTON ENERGY FOR  MT = ',  &       
     &                   L2H,'. TOO MANY MTS.'
                  CALL ERROR_MESSAGE(0)
               ELSE
                  NBS = NBS + 1
                  MTS(NBS) = L2H
               END IF
            END IF
         END DO
         GO TO 100
      END IF
!
!     OTHER FILR ONE SECTIONS
!
      SELECT CASE (MT)
!
         CASE (452)        ! FISSION NEUTRONS
            LNU = L2H
!
         CASE (455)        ! DELAYED FISSION NEUTRONS
            CALL RDLIST
!
         CASE (456)        ! PROMPT FISSION NEUTRONS
!
         CASE DEFAULT
            IERX = 1
            GO TO 100
!
      END SELECT
!
!     POLYNOMIAL REPRESENTATION
!
      IF(L2H.EQ.1)   THEN
         CALL RDLIST
         NRC = NPL
         DO J=1,NRC
            POLI(J) = Y(J)
         END DO
!
!     TABULAR REPRESENTATION
!
      ELSE
         CALL RDTAB1
         NPT = NP
         DO NN=1,NPT
            ETERP(NN) = X(NN)
            ANUTRP(NN) = Y(NN)
         END DO
         NRC = NR
         DO NN=1,NRC
            ITERP(NN) = NBT(NN)
            JTERP(NN) = JNT(NN)
         END DO
      END IF
!
  100 RETURN
      END SUBROUTINE CKF1
!
!***********************************************************************
!
      SUBROUTINE CKF2
!
!     CHECK FILE 2 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD, IFIX
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: LFW,LRU,LRF
      INTEGER(KIND=I4) :: NER,NAPS
      INTEGER(KIND=I4) :: IZ,IA
      INTEGER(KIND=I4) :: NI,NRE,ILI
      REAL(KIND=R4) :: EL,EU
      REAL(KIND=R4) :: AWRH,AHI,AWREL,DIFF
!
!     TEST ABUNDANCE WEIGHTED MASS OF A NATURAL ELEMENT
!
      IZ=NINT(C1H)/1000
      AWREL = AWRN(IZ)
      IF(N1H.GT.1) THEN
         IF(NISN.EQ.N1H)   THEN
            DIFF = ABS(1.-AWR/AWREL)
            IF(DIFF.GT.EPI4)   THEN
               WRITE(EMESS,'(A,F9.4)')  'AWR SHOULD BE ',AWREL
               CALL ERROR_MESSAGE(0)
            END IF
         END IF
      END IF
!
!     PROCESS ALL ISOTOPES
!
      E1 = ENMAX
      E2 = ENMIN
      NBOUND = 0
      NIS = N1H
      AWRH = C2H
      DO NI=1,NIS
!
!        INITIALIZE
!
         IPT = 1
         IRES = 0
         EMIDLE = 0.0
         NRED = 0
         NPED = 0
!
!        CONTROL RECORD
!
         CALL RDCONT
         ZAIS(NI) = C1H
         IA = MOD(IFIX(ZAIS(NI)),1000)
         ABNS(NI) = C2H
         LFW = L2H
!
!        SET LIMITS ON SCATTERING LENGTH
!
         IF(LRP.EQ.0) THEN
            AHI = AWRH*FACTOR
         ELSE
            AHI = AMOD(C1H,1000.)
         END IF
         APLO = AHI**OTHIRD
         APHI = APLO*0.17
         APLO = APLO*0.07
!
!        PROCESS EACH ENERGY RANGE
!
         IUNR = 0
         NER = N1H
         NERS(NI) = NER
         DO NRE=1,NER
            CALL RDCONT
            EL = C1H
            EU = C2H
            E1 = AMIN1(E1,EL)
            E2 = AMAX1(E2,EU)
!
!           STORE REGION BOUNDARIES
!
            IF(NBOUND.NE.0) THEN
               DO ILI=1,NBOUND
                  IF(RBOUND(ILI).EQ.EU)   GO TO 10
               END DO
            END IF
            NBOUND = NBOUND + 1
            RBOUND(NBOUND) = EU
!
!           STORE DATA FOR THIS REGION
!
   10       IRES(NRE) = IPT
            EMIDLE(NRE) = EU
            LRU = L1H
            CALL STOFP(FLOAT(LRU))
            IF(LRU.EQ.2)  IUNR = 1
            LRF = L2H
            CALL STOFP(FLOAT(LRF))
            NRO = N1H
            NAPS = N2H
!
!           ISOTOPE WITH NO RESONANCE PARAMETERS
!
            IF(LRU.EQ.0)   THEN
               CALL RDCONT
               IF(NIS.EQ.1) THEN
                  E1 = EL
                  E2 = EU
               END IF
               SPINS(NI) = C1H
               EMIDLE(NI) = EU
               CALL TESTAP(C2H,APLO,APHI)
               IF(IA.GT.0)  THEN
                  EMESS = ' '
                  CALL ERROR_MESSAGE(0)
                  WRITE(EMESS,'(A,I4,A)')                               &       
     &                'ISOTOPE MASS =',IA,'. HAS NO RESONANCES GIVEN.'
                  CALL ERROR_MESSAGE(0)
                  EMESS = ' '
                  CALL ERROR_MESSAGE(0)
               END IF
!
!           RESOLVED RESONANCE REGION
!
            ELSE IF(LRU.EQ.1)   THEN
               CALL CHKRR(IA,NI,NAPS,LRF,NRE,EL,EU)
!
!           UNRESOLVED REGION
!
            ELSE IF(LRU.EQ.2)   THEN
               CALL CHKURR(IA,NI,LFW,LRF,NRE)
            END IF
         END DO
!
!        SAVE RESONANCE REPRESENTATION FOR CURRENT ISOTOPE ON A SCRATCH
!          FILE FOR LATER USE
!
         CALL RWRES(1)
!
      END DO
!
!     REWIND STORAGE FILE SO IT IS READY TO BE READ
!
      REWIND (UNIT=ISCR)
!
!     TEST ABUNDANCES AND SPINS FOR EACH ISOTOPE
!
      CALL CKNKLD
!
!     DO RESONANCE CALCULATION FOR CAREN TEST
!
      CALL CAREN(0)
!
      RETURN
      END SUBROUTINE CKF2
!
!***********************************************************************
!
      SUBROUTINE CHKRR(IA,NI,NAPS,LRF,NRE,EL,EU)
!
!     ROUTINE TO CHECK RESOLVED ENERGY REGIONS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IA,NI,NAPS,LRF,NRE
      REAL(KIND=R4) :: EL,EU
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT
!
      INTEGER(KIND=I4) :: LRX,LI,LBK
      INTEGER(KIND=I4) :: NX,NCRE
      INTEGER(KIND=I4) :: L,NLS,NJS,NSS,NRS,NLJ,NLSJ,LPS
      INTEGER(KIND=I4) :: IRP,IPTU,IBEG,IRBEG,IBO
      INTEGER(KIND=I4) :: NC,NOFF
      INTEGER(KIND=I4) :: N,NN,III,NL,NJ,NS,LIL,NCR
      REAL(KIND=R4) :: APL,AWRI,ROC,RHOC,SER,PER,ROCED,ER,AER,RHO
      REAL(KIND=R4) :: AP,AS,AC
      REAL(KIND=R4) :: SPI,R,AJ,QX
!
!     GET DEPENDENT SCATTERING RADIUS
!
      IF(NRO.GT.0)   THEN
         CALL RDTAB1
         NRED(NRE) = NR
         NPED(NRE) = NP
         DO III=1,NR
            NBTED(III,NRE) = NBT(III)
            JNTED(III,NRE) = JNT(III)
         END DO
         DO III=1,NP
            CALL TESTAP(Y(III),APLO,APHI)
            EP(III,NRE) = X(III)
            APED(III,NRE) = Y(III)
         END DO
      END IF
!
!     GET RESONANCE PARAMETERS
!
      CALL RDCONT
      AP = C2H
      SPI = C1H
      SPINS(NI) = SPI
      CALL STOFP(SPI)
      CALL STOFP(AP)
      IF(NRO.EQ.0)  CALL TESTAP(C2H,APLO,APHI)
      NLS = N1H
      CALL STOFP(FLOAT(NLS))
!
!     PROCESS BREIT-WIGNER, REICH-MOORE, AND GENERALIZER R-MATRIX
!
      IF((LRF.GE.1.AND.LRF.LE.3).OR.LRF.EQ.5) THEN
!
!        PROCESS EACH L-VALUE
!
         DO NL=1,NLS
            CALL RDCONT
            IF(NL.EQ.1)   THEN
               AWRI = C1H
               CALL STOFP(AWRI)
            END IF
            IF(LRF.EQ.3) THEN
               IF(C2H.EQ.0.0)  THEN
                  APL = AP
               ELSE
                  APL = C2H
               END IF
               CALL STOFP(APL)
            ELSE
               APL = AP
            END IF
            IF(NAPS.EQ.0)  THEN
               R = .123*ANEUTR*AWRI**(1./3.) + 0.08
            ELSE
               R = APL
            END IF
            ROC = CROC*(AWRI/(AWRI+1.0))*R
            IF(NL.EQ.1.OR.LRF.EQ.3) CALL STOFP(R)
            QX = C2H
            CALL STOFP(QX)
            L = L1H
            CALL STOFP(FLOAT(L))
            LRX = L2H
            CALL STOFP(FLOAT(LRX))
            NRS = N2H
            CALL STOFP(FLOAT(NRS))
            IRP = IPT
!
!           PROCESS EACH RESONANCE
!
            DO N=1,NRS
               IPTU = IPT + 5
               READ(JIN,'(6E11.5)')   (RESPAR(NN),NN=IPT,IPTU)
               ER = RESPAR(IPT)
               AER = ABS(ER)
               RHO = ROC*SQRT(AER)
               IF(NRO.EQ.1.AND.NAPS.EQ.1)   THEN
                  IBEG = 1
                  CALL TERPR(AER,EP(1,NI),APED(1,NI),NPED(NI),          &       
     &                   NBTED(1,NI),JNTED(1,NI),NRED(NI),IBEG,ROCED)
                  RHO = RHO*ROCED/ROC
               END IF
               CALL FACTS(L,RHO,RESPAR(IPTU+1),RESPAR(IPTU+2))
               RESPAR(IPTU+3) = 0.0
               IF(LRF.EQ.3)   THEN
                  RESPAR(IPTU+3) = RESPAR(IPT+5)
                  RESPAR(IPT+5) = RESPAR(IPT+4)
                  RESPAR(IPT+4) = RESPAR(IPT+3)
                  RESPAR(IPT+3) = RESPAR(IPT+2)
               END IF
               IF(LRX.NE.0)   THEN
                  RHOC = ROC*SQRT(ABS(ER+QX))
                  CALL FACTS(L,RHOC,SER,PER)
                  RESPAR(IPTU+3) = (RESPAR(IPT+2)-RESPAR(IPT+3)-        &       
     &               RESPAR(IPT+4)-RESPAR(IPT+5))/PER
                END IF
                IPT = IPTU + 4
            END DO
!
!           PERFORM RESOLVED RESONANCE PARAMETER TESTS
!
            CALL TESTRP(IA,AWRI,SPI,AJ,R,L,RESPAR(IRP),9*NRS,NRS,       &       
     &                 EL,EU,LRF)
         END DO
!
!     PROCESS ADLER-ADLER REPRESENTATION
!
      ELSE IF(LRF.EQ.4) THEN
         DO NL=1,NLS
            CALL RDCONT
            IF(NL.LE.1)   THEN
               AWRI = C1H
               CALL STOFP(AWRI)
               IF(NAPS.NE.0)  THEN
                  R = AP
               ELSE
                  R = .123*ANEUTR*AWRI**(1./3.) + 0.08
               END IF
            END IF
            LI = L1H
            CALL STOFP(FLOAT(LI))
            NX = N2H
            CALL STOFP(FLOAT(NX))
            CALL PAKLIS(JIN,NX)
            CALL RDCONT
            L = L1H
            CALL STOFP(FLOAT(L))
            NJS = N1H
            CALL STOFP(FLOAT(NJS))
            DO NJ=1,NJS
               CALL RDCONT
               AJ = C1H
               CALL STOFP(AJ)
               NLJ = N2H
               CALL STOFP(FLOAT(NLJ))
               NC = 2*NLJ
               CALL PAKLIS(JIN,NC)
            END DO
         END DO
         EMESS = ' '
         CALL ERROR_MESSAGE(0)
         EMESS = 'NO RESONANCE PARAMETER TESTS IMPLEMENTED FOR '//      &       
     &              'ADLER-ADLER'
         CALL ERROR_MESSAGE(0)
!
!     HYBRID R-FUNCTION
!
      ELSE IF(LRF.EQ.6) THEN
         CALL RDCONT
         CALL STOFP(FLOAT(L1H))
         CALL STOFP(FLOAT(L2H))
         CALL STOFP(FLOAT(N1H))
         CALL STOFP(FLOAT(N2H))
         NCRE = N2H
!********READ REACTION CHANNEL DEFINITIONS
         CALL RDCONT
         CALL STOFP(FLOAT(L1H))
         CALL STOFP(FLOAT(L2H))
         CALL STOFP(FLOAT(N1H))
         CALL STOFP(FLOAT(N2H))
!********READ REACION CHANNEL Q-VALUES
         CALL RDLIST
         CALL PAKLIS(JIN,1)
         IPT = IPT - 2
         IRBEG = IPT
         IPT = IPT + 4*NCRE + 1
!********READ ANY CHARGED PARTICLE PENETRABILITIES
         IF(NCRE.GT.0)  THEN
            NOFF = 1
            DO NCR=1,NCRE
               DO LIL=1,4
                  CALL RDTAB1
                  RESPAR(IRBEG+NOFF) = IPT
                  CALL PKTAB1
                  NOFF = NOFF + 1
               END DO
            END DO
         END IF
!
!        PROCESS EACH L, S, AND J VALUE
!
         RESPAR(IRBEG) = IPT
!********L-VALUE
         DO NL=1,NLS
            CALL RDCONT
            AWRI = C1H
            CALL STOFP(AWRI)
            L = L1H
            CALL STOFP(FLOAT(L1H))
            NSS = N1H
            CALL STOFP(FLOAT(NSS))
!***********CHANNEL SPIN
            DO NS=1,NSS
               CALL RDCONT
               AS = C1H
               CALL STOFP(AS)
               NJS = N1H
               CALL STOFP(FLOAT(NJS))
!**************TOTAL SPIN
               DO NJ=1,NJS
                  CALL RDLIST
                  AJ = C1L
                  CALL STOFP(AJ)
                  AC = C2L
                  CALL STOFP(AC)
                  CALL TESTAP(C2L,APLO,APHI)
                  R = AC
                  LBK = L1H
                  CALL STOFP(FLOAT(LBK))
                  LPS = L2H
                  CALL STOFP(FLOAT(LPS))
                  NLSJ = N2H
                  CALL STOFP(FLOAT(NLSJ))
                  NRS = NLSJ
                  NC = 2*NLSJ
                  IRP = IPT
                  CALL PAKLIS(JIN,NC)
                  IBO = IPT
                  IPT = IPT + 1
!*****************READ BACKGROUND
                  IF(LBK.NE.0) THEN
                     CALL PKTAB1
                     CALL PKTAB1
                  END IF
!*****************READ PHASE SHIFTS
                  IF(LPS.NE.0)  THEN
                     CALL PKTAB1
                     CALL PKTAB1
                  END IF
                  RESPAR(IBO) = IPT
!
!                 PERFORM RESOLVED RESONANCE PARAMETER TESTS
!
                  CALL TESTRP(IA,AWRI,SPI,AJ,R,L,RESPAR(IRP),           &       
     &                     12*NRS,NRS,EL,EU,LRF)
               END DO
            END DO
         END DO
!
!     R-matrix-limited
!
      ELSE IF(LRF.EQ.7) THEN
         CALL RDLIST
         DO NJ=1,NLS
            CALL RDLIST
            CALL RDLIST
         END DO
      END IF
!
      RETURN
      END SUBROUTINE CHKRR
!
!***********************************************************************
!
      SUBROUTINE CHKURR(IA,NI,LFW,LRF,NRE)
!
!     TEST UNRESOLVED RESONANCE REGION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IA,NI,LFW,LRF,NRE,III
!
      INTEGER(KIND=I4) :: LSSF
      INTEGER(KIND=I4) :: L,NRS,NLS,NJS,MUF
      INTEGER(KIND=I4) :: IRP,NNINT,NE,NC
      INTEGER(KIND=I4) :: NL,NJ,J
      REAL(KIND=R4) :: A,AWRI,AJ,SPI,R
!
!     GET DEPENDENT SCATTERING RADIUS
!
      IF(NRO.GT.0)   THEN
         CALL RDTAB1
         NRED(NRE) = NR
         NPED(NRE) = NP
         DO III=1,NR
            NBTED(III,NRE) = NBT(III)
            JNTED(III,NRE) = JNT(III)
         END DO
         DO III=1,NP
            CALL TESTAP(Y(III),APLO,APHI)
            EP(III,NRE) = X(III)
            APED(III,NRE) = Y(III)
         END DO
      END IF
!
!
      CALL STOFP(FLOAT(LFW))
!
!     ALL PARAMETERS ENERGY DEPENDENT
!
      IF(LRF.EQ.2)   THEN
         CALL RDCONT
         SPI = C1H
         SPINS(NI) = SPI
         CALL STOFP(SPI)
         LSSF = L1H
         CALL STOFP(FLOAT(LSSF))
         A = C2H
         CALL STOFP(A)
         CALL TESTAP(C2H,APLO,APHI)
         NLS = N1H
         CALL STOFP(FLOAT(NLS))
!
!        PROCESS EACH L VALUE
!
         DO NL=1,NLS
            CALL RDCONT
            IF(NL.EQ.1)   THEN
               AWRI = C1H
               CALL STOFP(AWRI)
               R = .123*ANEUTR*AWRI**(1./3.) + 0.08
               CALL STOFP(R)
            END IF
            L = L1H
            CALL STOFP(FLOAT(L))
            NJS = N1H
            CALL STOFP(FLOAT(NJS))
            NRS = NJS
            IRP = IPT
            DO J=1,NJS
               CALL RDCONT
               AJ = C1H
               CALL STOFP(AJ)
               NNINT = L1H
               CALL STOFP(FLOAT(NNINT))
               NE = N2H
               CALL STOFP(FLOAT(NE))
               CALL PAKLIS(JIN,NE+1)
            END DO
!
!           PERFORM UNRESOLVED RESONANCE PARAMETER TESTS
!
            CALL TESTUR(3,IA,L,NRS,AJ,RESPAR(IRP),NE,IPT-IRP+1)
         END DO
!
!     ALL PARAMETERS ENERGY INDEPENDENT
!
      ELSE
         IF(LFW.EQ.0)   THEN
            CALL RDCONT
            SPI = C1H
            SPINS(NI) = SPI
            CALL STOFP(SPI)
            LSSF = L1H
            CALL STOFP(FLOAT(LSSF))
            A = C2H
            CALL STOFP(A)
            CALL TESTAP(C2H,APLO,APHI)
            NE = 0
            CALL STOFP(FLOAT(NE))
            NLS = N1H
            CALL STOFP(FLOAT(NLS))
            DO NL=1,NLS
               CALL RDCONT
               IF(NL.EQ.1)   THEN
                  AWRI = C1H
                  CALL STOFP(AWRI)
                  R = .123*ANEUTR*AWRI**(1./3.) + 0.08
                  CALL STOFP(R)
               END IF
               L = L1H
               CALL STOFP(FLOAT(L))
               NJS = N2H
               CALL STOFP(FLOAT(NJS))
               NRS = NJS
               IRP = IPT
               CALL PAKLIS(JIN,NJS)
!
!              PERFORM UNRESOLVED RESONANCE PARAMETER TESTS
!
               CALL TESTUR(1,IA,L,NRS,AJ,RESPAR(IRP),0,IPT-IRP+1)
            END DO
!
!        FISSION WIDTHS ARE ENERGY DEPENDENT
!
         ELSE
            CALL RDCONT
            SPI = C1H
            SPINS(NI) = SPI
            CALL STOFP(SPI)
            LSSF = L1H
            CALL STOFP(FLOAT(LSSF))
            A = C2H
            CALL STOFP(A)
            NE = N1H
            CALL STOFP(FLOAT(NE))
            NLS = N2H
            CALL STOFP(FLOAT(NLS))
            NC = (NE+5)/6
            CALL PAKLIS(JIN,NC)
            DO NL=1,NLS
               CALL RDCONT
               IF(NL.EQ.1)   THEN
                  AWRI = C1H
                  CALL STOFP(AWRI)
                  R = .123*ANEUTR*AWRI**(1./3.) + 0.08
                  CALL STOFP(R)
               END IF
               L = L1H
               CALL STOFP(FLOAT(L))
               NJS = N1H
               CALL STOFP(FLOAT(NJS))
               NRS = NJS
               IRP = IPT
               DO NJ=1,NJS
                  CALL RDCONT
                  MUF = L2H
                  CALL STOFP(FLOAT(MUF))
                  CALL PAKLIS(JIN,NC+1)
               END DO
!
!              PERFORM UNRESOLVED RESONANCE PARAMETER TESTS
!
               CALL TESTUR(2,IA,L,NRS,AJ,RESPAR(IRP),NE,IPT-IRP+1)
            END DO
         END IF
      END IF
!
      RETURN
      END SUBROUTINE CHKURR
!
!***********************************************************************
!
      SUBROUTINE TESTAP(Z,A,B)
!
!     CHECK THAT POTENTIAL SCATTERING RADIUS LIES BETWEEN A AND B
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: Z,A,B
!
      IF(Z.LT.A.OR.Z.GT.B)   THEN
!
!        ERROR MESSAGE
!
         WRITE(EMESS,'(A,1PE13.5,A,1PE13.5)')                           &       
     &               '  AP NOT IN RANGE',A,' TO',B
         CALL ERROR_MESSAGE(NSEQP)
      END IF
!
      RETURN
      END SUBROUTINE TESTAP
!
!***********************************************************************
!
      SUBROUTINE TESTRP(IA,AWRI,SPI,AJ,R,L,A,NN,NRS,EL,EU,LRF)
!
!     SUBROUTINE PERFORMS PSYCHE TESTS ON SINGLE- AND MULTI-LEVEL
!      BREIT-WIGNER AND R-MATRIX RESOLVED RESONANCE PARAMETERS
!
!     A IS AN ARRAY CONTAINING THE RESONANCE PARAMETERS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IA,L,NN,NRS,LRF
      REAL(KIND=R4) :: AWRI,R,AJ,SPI,EL,EU
      REAL(KIND=R4), DIMENSION(NN) :: A
!
      REAL(KIND=R4), INTRINSIC :: ABS,SQRT
!
      INTEGER(KIND=I4) :: NREP,NRSA,INOW
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: DE,FNRS,DNRS
      REAL(KIND=R4) :: GDEN,ROC,GNOEST,GGH,GG2
      REAL(KIND=R4) :: GTBAR,GGBAR,GFBAR,GGLOG,GNOBAR,GNOLOG,GFLOG,SF,D
      REAL(KIND=R4) :: E,FL,RHO,SE,PE,GNO,AJC,GG,GF
      REAL(KIND=R4) :: ENUG,DELNUG,ENUF,DELNUF,ENUN,DELNUN,YIE
!
      REAL(KIND=R4), PARAMETER :: CGNO=1.66336E+06
!
      NREP = NN/NRS
!
!     CALCULATE AVERAGES OF RESONANCE PARAMETERS AND THEIR LOGS
!
      GTBAR = 0.0
      GGBAR = 0.0
      GGLOG = 0.0
      GFBAR = 0.0
      GFLOG = 0.0
      NRSA = 0
      DO N=1,NRS
         INOW = NREP*(N-1) + 1
         IF(A(INOW).GE.EL.AND.A(INOW).LE.EU)  THEN
            NRSA = NRSA + 1
            IF(LRF.EQ.3) THEN
               GG = A(INOW+4)
               GF = ABS(A(INOW+5))
            ELSE IF(LRF.EQ.6) THEN
               GG = A(INOW+2)
               GF = A(INOW+3)
            ELSE
               GTBAR = GTBAR + A(INOW+2)
               GG = A(INOW+4)
               GF = ABS(A(INOW+5))
            END IF
            GGBAR = GGBAR + GG
            IF(GG.GT.0.) THEN
               GGLOG = GGLOG + ALOG(GG)
            END IF
            IF(GF.NE.0.)   THEN
               GFBAR = GFBAR + GF
               GFLOG = GFLOG + ALOG(GF)
            END IF
         END IF
      END DO
      IF(NRSA.EQ.0) GO TO 100
      DE = EU - EL
      FNRS = NRSA
      DNRS = SQRT (1./FNRS)
      IF(LRF.LT.3)  GTBAR = GTBAR/FNRS
      GGBAR = GGBAR/FNRS
      GGLOG = GGLOG/FNRS
      GFBAR = GFBAR/FNRS
      GFLOG = GFLOG/FNRS
!
!     CALCULATE CONSTANTS FOR THE TESTS
!
      GDEN = 2.*(2.*SPI+1.)
      ROC = CROC*(AWRI/(AWRI+1.0))*R
      GNOEST = CGNO*ROC/R/R
      GGH = 0.33*GGBAR
      GG2 = 3.0*GGBAR
!
!     DO TESTS ON GAMMA-G, REDUCED GAMMA-N, AND THE STRENGTH
!       FUNCTION
!
      GNOBAR = 0.0
      GNOLOG = 0.0
      SF = 0.0
      WRITE(NOUT,10)  IA,L,NRS,NRSA
   10 FORMAT(//5X,'ISOTOPE MASS =',I4,'. L = ',I2//                     &       
     &  9X,'TOTAL NUMBER OF RESONANCES IS ',16X,I4/                     &       
     &  9X,'NUMBER OF RESONANCES IN THE ENERGY REGION ARE ',I4/)
      DO N=1,NRS
         INOW = NREP*(N-1) + 1
         IF(A(INOW).GE.EL.AND.A(INOW).LE.EU)  THEN
!***********DETERMINE THE REDUCED NEUTRON WIDTH
            E = ABS(A(INOW))
            RHO = ROC*SQRT(E)
            CALL FACTS(L,RHO,SE,PE)
            IF(LRF.NE.6)  THEN
               GNO = A(INOW+3)*ROC/PE
               GG = A(INOW+4)
               AJC = A(INOW+1)
            ELSE
               GNO = A(INOW+1)*ROC/PE
               GG = A(INOW+2)
               AJC = AJ
            END IF
            SF = SF + (2.*AJC+1.)*GNO/GDEN
            GNOBAR = GNOBAR + GNO
            GNOLOG = GNOLOG + ALOG(GNO)
!
!           TEST FOR EXTREME VALUES OF GAMMA-G AND REDUCED GAMMA-N
!
            IF(GNO.GT.GNOEST)  THEN
               WRITE(NOUT,20)  E,GNO,GNOEST
   20          FORMAT(9X,'AT RESONANCE ENERGY ',1PE12.5,                &       
     &                ' EV. THE REDUCED NEUTRON WIDTH ',1PE12.5,        &       
     &                 ' EXCEEDS THE WIGNER ESTIMATE ',1PE12.5)
            END IF
            IF(GG.LT.GGH.OR.GG.GT.GG2)   THEN
               WRITE(NOUT,30) E,GG,GGBAR
   30          FORMAT(9X,'AT RESONANCE ENERGY ',1PE12.5,                &       
     &                ' EV. THE GAMMA WIDTH ',1PE12.5,                  &       
     &                ' DEVIATES TOO MUCH FROM THE AVERAGE   ',1PE12.5)
            END IF
         END IF
      END DO
!
!     CALCULATE ESTIMATES OF THE DEGREES OF FREEDOM FOR THE
!      RESONANCE PARAMETER DISTRIBUTIONS
!
      GNOBAR = GNOBAR/FNRS
      GNOLOG = GNOLOG/FNRS
      FL = FLOAT (2*L+1)
      D = 0.0
      IF(NRSA.GT.1) THEN
         SF = SF/(FL*DE)
         D = DE/(FNRS-1.0)
      END IF
      DELNUN = DNRS
      YIE = GNOLOG - ALOG(GNOBAR)
      CALL CHISQR(YIE,ENUN,DELNUN)
      DELNUG = DNRS
      YIE = GGLOG - ALOG(GGBAR)
      CALL CHISQR(YIE,ENUG,DELNUG)
      IF(GFBAR.NE.0.)  THEN
         DELNUF = DNRS
         YIE = GFLOG - ALOG (GFBAR)
         CALL CHISQR(YIE,ENUF,DELNUF)
      END IF
!
!     SUMMARIZE RESULTS OF CALCULATIONS
!
      IF(LRF.LT.3) THEN
         WRITE(NOUT,40) GTBAR
   40    FORMAT(9X,'AVERAGE TOTAL WIDTH IS ',1PE12.5)
      END IF
      WRITE(NOUT,45) GNOBAR
   45 FORMAT(9X,'AVERAGE REDUCED NEUTRON WIDTH IS ',1PE12.5)
      WRITE(NOUT,50)  ENUN,DELNUN
   50 FORMAT(14X,'THE NUMBER OF DEGREES OF FREEDOM IS',F7.3,' +OR-',    &       
     &    F6.3)
      WRITE(NOUT,55) GGBAR
   55 FORMAT(9X,'AVERAGE GAMMA WIDTH IS ',1PE12.5)
      IF(ENUG.LE.10.)   THEN
         WRITE(NOUT,50) ENUG,DELNUG
      ELSE
         WRITE(NOUT,60)
   60    FORMAT(14X,'NU IS GREATER THAN 10.')
      END IF
      IF(GFBAR.NE.0.0.AND.LRF.NE.3)  THEN
         WRITE(NOUT,65) GFBAR
   65    FORMAT(9X,'AVERAGE FISSION WIDTH IS ',1PE12.5)
         IF(ENUF.LE.10.)    THEN
            WRITE(NOUT,50) ENUF,DELNUF
         ELSE
            WRITE(NOUT,60)
         END IF
      END IF
      WRITE(NOUT,70)  D,SF
   70 FORMAT(9X,'AVERAGE LEVEL SPACING IS ',1PE12.5/                    &       
     &       9X,'STRENGTH FUNCTION IS ',1PE12.5)
!
!     COMPARE THE STRENGTH FUNCTION AND THE AVERAGE GAMMA WIDTH
!       WITH SYSTEMATICS
!
      CALL TESTSF(L,SF,GGBAR)
!
  100 RETURN
      END SUBROUTINE TESTRP
!
!***********************************************************************
!
      SUBROUTINE TESTSF(L,SF,GGBAR)
!
!     TESTS THE STRENGTH FUNCTION AND THE AVERAGE GAMMA WIDTH
!      AGAINST SYSTEMATICS FROM BNL-325 VOL I
!
!     L IS THE NEUTRON ANGULAR MOMENTUM
!     SF IS THE STRENGTH FUNCTION
!     GGBAR IS THE AVERAGE GAMMA WIDTH
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L
      REAL(KIND=R4) :: SF,GGBAR
!
      INTEGER(KIND=I4) :: MS,LP1,ITST
      INTEGER(KIND=I4) :: I,K
      REAL(KIND=R4) :: AW
!
      REAL(KIND=R4), DIMENSION(4), PARAMETER ::                          &      
     &         AWT = (/80.,140.,210.,260./)
      REAL(KIND=R4), DIMENSION(2,2) :: GT
      DATA GT/.04,.01,9.0,.05/
      REAL(KIND=R4), DIMENSION(4,2,2) :: SFL
      DATA SFL/1.0E-04,1.0E-05,1.0E-04,4.0E-05,1.0E-05,1.0E-04,         &       
     &         1.0E-05,1.0E-04,9.0E-04,1.0E-04,6.0E-04,2.0E-04,         &       
     &         2.0E-04,8.0E-04,1.5E-04,3.0E-04/
!
!     TEST AVERAGE GAMMA WIDTH
!
      AW = ANEUTR*AWR
      IF(AW.LE.210.)   THEN
         MS = 1
      ELSE
         MS = 2
      END IF
      IF(GGBAR.LT.GT(MS,1).OR.GGBAR.GT.GT(MS,2)) THEN
         WRITE(NOUT,10)  GGBAR,(GT(MS,I),I=1,2)
   10    FORMAT(14X,'AVERAGE GAMMA WIDTH',1PE12.5/                      &       
     &    19X,' LIES OUTSIDE LIMITS',1PE12.5,' TO',1PE12.5,' EV.')
      END IF
!
!     TEST STRENGTH FUNCTION
!
      LP1 = L + 1
      IF(LP1.LE.2)   THEN
         DO I=1,4
            IF(AW.LE.AWT(I)) THEN
               ITST = I
               IF(SF.LT.SFL(ITST,LP1,1).OR.SF.GT.SFL(ITST,LP1,2)) THEN
                  WRITE(NOUT,20)  SF,(SFL(ITST,LP1,K),K=1,2)
   20             FORMAT(14X,'STRENGTH FUNCTION',1PE12.5/               &       
     &                19X,' LIES OUTSIDE LIMITS',1PE12.5,' TO',1PE12.5)
               END IF
               GO TO 100
            END IF
         END DO
      END IF
!
  100 RETURN
!
      END SUBROUTINE TESTSF
!
!***********************************************************************
!
      SUBROUTINE CHISQR(YIE,ENU,DELNU)
!
!     ESTIMATES THE NUMBER OF DEGREES OF FREEDOM FOR A SET OF
!      RESONANCE PARAMETERS WHOSE DISTRIBUTION LAW IS ASSUMED
!      TO BE A CHI-SQUARE
!
!     YIE IS AVERAGE OF LOGS - LOG OF THE AVERAGE
!     ENU IS THE NUMBER OF DEGREES OF FREEDOM CALCULATED
!     DELNU IS THE INVERSE OF THE SQRT OF THE NUMBER OF
!      PARAMETERS ON INPUT AND THE STANDARD DEVIATION
!      ON OUTPUT
!
!     THE ARRAY YIEFAC CONTAINS VALUES OF YIE FOR VALUES OF ENU
!      FROM 0.25 TO 10. IN STEPS OF 0.25, SO THAT ENU FOR
!      A VALUE OF YIE CAN BE OBTAINED BY INTERPOLATION
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: YIE,ENU,DELNU
!
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) AIPT,YIEL,YIEU,ENUM
!
      INTEGER(KIND=I4), PARAMETER :: IYIEMAX=40
      REAL(KIND=R4), DIMENSION(IYIEMAX), PARAMETER ::                    &      
     &   YIEFAC=(/                                                      &       
     &       -6.30961, -2.84121, -1.77349, -1.27035, -.982946, -.798151,&       
     &       -.670626, -.577216, -.506833, -.450644, -.406104, -.368965,&       
     &       -.338458, -.312116, -.289909, -.270363, -.253933, -.238430,&       
     &       -.225375, -.213124, -.202646, -.192672, -.184019, -.175828,&       
     &       -.169007, -.161711, -.155720, -.149596, -.144467, -.139191,&       
     &       -.134686, -.130177, -.126639, -.122282, -.118935, -.115196,&       
     &       -.112227, -.108913, -.106196, -.103320/)
!
!     RETURN IF NU OUTSIDE LIMITS 0. TO 10.
!
      ENU = 0.
      IF(YIE.LT.-6.)   GO TO 100
      ENU = 99.
      IF(YIE.GT.-.104)   GO TO 100
!
!     OBTAIN ENU BY INTERPOLATION
!
      DO N=1,IYIEMAX
         YIEU = YIEFAC(N)
         IF(YIE.LE.YIEU)   THEN
            AIPT = FLOAT(N-1)
            YIEL = YIEFAC(N-1)
            ENUM = 0.25*AIPT
            ENU = ENUM + 0.25*(YIE-YIEL)/(YIEU-YIEL)
!           DELNU OBTAINED FROM LINEAR RELATION ON PAGE 224 OF J.E.LYNN
!          'THE THEORY OF NEUTRON RESONANCE REACTIONS'.
!           NOTE THAT DELNU IS NEGATIVE FOR ENU LE 0.5
            DELNU = (1.4421*(ENU-0.5))*DELNU
            GO TO 100
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE CHISQR
!
!***********************************************************************
!
      SUBROUTINE TESTUR(IPATH,IA,L,NRS,AJ,A,NE,NN)
!
!     SUBROUTINE PERFORMS PSYCHE TESTS ON THE UNRESOLVED RESONANCE
!     REGION PARAMETERS
!
!       IPATH = 1  LRF=1, LFW=0
!               2  LRF=1, LFW=1
!               3  LRF=2, LFW=ANYTHING
!       A IS THE ARRAY OF RESONANCE PARAMETERS
!       NE IS THE NUMBER OF ENERGY GRID POINTS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPATH,IA,NE,NN,L,NRS
      REAL(KIND=R4) :: AJ
      REAL(KIND=R4), DIMENSION(NN) :: A
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: IOF,INOW,LENN,NERRS,ITF1,ITF2
      INTEGER(KIND=I4) :: N,K
      REAL(KIND=R4) :: DENSD,AJP1,CONSD,SFC
      REAL(KIND=R4) :: DENS,AJPONE,CONS,DENOM,SF,DIF
!
!
!     WRITE OUT HEADING
!
      WRITE(NOUT,5)  IA,L
    5 FORMAT(//5X,'ISOTOPE MASS = ',I4,'. L = ',I2,                     &       
     &       ' UNRESOLVED RESONANCE REGION'/)
!
!     ONLY FISSION WIDTH CAN BE ENERGY DEPENDENT
!
      IF(IPATH.LE.2)   THEN
         IOF = IPATH - 1
         DENSD = A(IOF+1)
         AJP1 = 2.*A(IOF+2)+1
         CONSD = DENSD*AJP1
         SFC = A(IOF+4)/DENSD
!
!        TEST STRENGTH FUNCTION AND LEVEL SPACING
!
         LENN = IOF +6*((NE+5)/6+1)
         DO N=1,NRS
            INOW = LENN*(N-1) + IOF + 1
            DENS = A(INOW)
            AJ = A(INOW+1)
            AJPONE = 2.*AJ+1
            CONS = DENS*AJPONE
            SF = A(INOW+3)/DENS
            WRITE(EMESS,'(6X,A,F4.1,A,1PE12.5)')                        &       
     &             'J = ',AJ,'. STRENGTH FUNCTION IS ',SF
            CALL ERROR_MESSAGE(0)
            DIF = ABS(CONS-CONSD)/CONSD
            IF(DIF.GT.0.01)  THEN
               DENOM = AJP1*DENSD/AJPONE
               WRITE(EMESS,'(11X,A,1PE12.5,A,1PE12.5)')                 &       
     &           'DENSITY ',DENS,' SHOULD BE ',DENOM
               CALL ERROR_MESSAGE(0)
            END IF
            DIF = ABS(SF-SFC)/SFC
            IF(DIF.GT.EPI2) THEN
               EMESS = '         STRENGTH FUNCTION IS NOT '//           &       
     &               'INDEPENDENT OF J'
               CALL ERROR_MESSAGE(0)
            END IF
!
!           COMPARE AVERAGE GAMMA WIDTH AND STRENGTH FUNCTION WITH
!            SYSTEMATICS
!
            CALL TESTSF(L,SF,A(INOW+4))
         END DO
         GO TO 100
       END IF
!
!     ALL PARAMETERS ARE ENERGY DEPENDENT
!
      IOF = 0
      AJP1 = 2.*A(1)+1.
      DO K=1,NE
         CONSA(K) = A(IOF+11)*AJP1
         SFA(K) = A(IOF+13)/A(IOF+11)
         IOF = IOF + 6
      END DO
!
!     TEST STRENGTH FUNCTION AND LEVEL SPACING
!
      LENN = 6*NE + 9
      DO N=1,NRS
         INOW = LENN*(N-1) + 1
         AJ = A(INOW)
         AJPONE = 2.*AJ+1.
         DENS = A(INOW+10)
         CONS = DENS*AJPONE
         SF = A(INOW+12)/DENS
         WRITE(EMESS,'(A,F4.1)') ' J = ',AJ
         CALL ERROR_MESSAGE(0)
         WRITE(EMESS,'(A,1PE12.5,A,1PE12.5)')                           &       
     &           ' ENERGY = ',A(INOW+9),'. STRENGTH FUNCTION IS ',SF
         CALL ERROR_MESSAGE(0)
         NERRS = 0
         DO K=1,NE
            DENS = A(INOW+10)
            CONS = DENS*AJPONE
            SF = A(INOW+12)/DENS
            CONSD = CONSA(K)
            DIF = ABS(CONS-CONSD)/CONSD
            DENOM = CONSD/AJPONE
            IF(DIF.GT.EPI2)  THEN
               ITF1 = 1
            ELSE
               ITF1 = 0
            END IF
            SFC = SFA(K)
            DIF = ABS(SF-SFC)/SFC
            IF(DIF.GT.EPI2)  THEN
               ITF2 = 1
            ELSE
               ITF2 = 0
            END IF
            IF(ITF1.NE.0.OR.ITF2.NE.0)   THEN
               NERRS = NERRS + 1
               IF(NERRS.LE.20)   THEN
                  WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')               &       
     &               'ENERGY = ',A(INOW+9),'. STRENGTH FUNCTION IS ',SF
                  CALL ERROR_MESSAGE(0)
                  IF(ITF1.EQ.1)  THEN
                     WRITE(EMESS,'(9X,A,1PE12.5,A,1PE12.5)')            &       
     &                  'DENSITY ',DENS,' SHOULD BE ',DENOM
                     CALL ERROR_MESSAGE(0)
                  END IF
                  IF(ITF2.EQ.1)  THEN
                     EMESS = '         STRENGTH FUNCTION IS NOT '//     &       
     &                  'INDEPENDENT OF J'
                      CALL ERROR_MESSAGE(0)
                  END IF
               END IF
            END IF
!
!           COMPARE AVERAGE GAMMA WIDTH AND STRENGTH FUNCTION WITH
!            SYSTEMATICS
!
            CALL TESTSF(L,SF,A(INOW+13))
            INOW = INOW + 6
         END DO
         IF(NERRS.GT.20)  THEN
            WRITE(EMESS,'(4X,A,I5,A)')                                  &       
     &           'TOTAL OF ',NERRS,' POINTS IN ERROR FOR THIS J-STATE'
            CALL ERROR_MESSAGE(0)
         END IF
         WRITE(EMESS,'(A,1PE12.5,A,1PE12.5)')                           &       
     &           ' ENERGY = ',A(INOW+3),'. STRENGTH FUNCTION IS ',SF
         CALL ERROR_MESSAGE(0)
         EMESS = '  '
         CALL ERROR_MESSAGE(0)
      END DO
!
  100 RETURN
      END SUBROUTINE TESTUR
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION AWRN(IZ)
!
!     EXTRACTS ISOTOPIC PARAMETERS OF AN ELEMENT FROM MASTER TABLE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IZ
!
      INTEGER(KIND=I4), INTRINSIC :: IABS, MOD
      REAL(KIND=R4), INTRINSIC :: ABS, FLOAT
!
      INTEGER(KIND=I4) :: NTHEL,IZN,IZN1
!
      INTEGER(KIND=I4), PARAMETER :: KISO=287
      INTEGER(KIND=I4) :: I
      INTEGER(KIND=I4), DIMENSION(KISO) :: IZASP
      DATA (IZASP(I),I=1,79)/                                           &       
     &  100105,  100210,  200305,  200400,                              &       
     &  300610, -300715, -400915,  501030, -501115,                     &       
     &  601200, -601305,  701410, -701505,  801600,                     &       
     &  801725,  801800,  901905, 1002000, 1002115,                     &       
     & 1002200, 1102315, 1202400, 1202525, 1202600,                     &       
     & 1302725, 1402800, 1402905, 1403000, 1503105,                     &       
     & 1603200, 1603315, 1603400, 1603600, 1703515,                     &       
     & 1703715, 1803600, 1803800, 1804000, 1903915,                     &       
     &-1904040, 1904115, 2004000, 2004200,-2004335,                     &       
     & 2004400, 2004600, 2004800,-2104535, 2204600,                     &       
     &-2204725, 2204800,-2204935, 2205000, 2305060,                     &       
     &-2305135, 2405000, 2405200,-2405315, 2405400,                     &       
     &-2505525, 2605400, 2605600,-2605705, 2605800,                     &       
     &-2705935, 2805800, 2806000,-2806115, 2806200,                     &       
     & 2806400,-2906315,-2906515, 3006400, 3006600,                     &       
     &-3006725, 3006800, 3007000,-3106915,-3107115/
      DATA (IZASP(I),I=80,159)/                                         &       
     & 3207000, 3207200, 3207345, 3207400, 3207600,                     &       
     &-3307515, 3407400, 3407600,-3407705, 3407800,                     &       
     & 3408000, 3408200,-3507915,-3508115, 3607800,                     &       
     & 3608000, 3608200, 3608345, 3608400, 3608600,                     &       
     &-3708525, 3708715, 3808400, 3808600, 3808745,                     &       
     & 3808800,-3908905, 4009000, 4009125, 4009200,                     &       
     & 4009400, 4009600, 4109345, 4209200, 4209400,                     &       
     & 4209525, 4209600, 4209725, 4209800, 4210000,                     &       
     & 4409600, 4409800, 4409925, 4410000, 4410125,                     &       
     & 4410200, 4410400,-4510305, 4610200, 4610400,                     &       
     & 4610525, 4610600, 4610800, 4611000,-4710705,                     &       
     &-4710905, 4810600, 4810800, 4811000, 4811105,                     &       
     & 4811200, 4811305, 4811400, 4811600, 4911345,                     &       
     & 4911545, 5011200, 5011400, 5011505, 5011600,                     &       
     & 5011705, 5011800, 5011905, 5012000, 5012200,                     &       
     & 5012400, 5112125, 5112335, 5212000, 5212200/
      DATA (IZASP(I),I=160,239)/                                        &       
     & 5212305, 5212400, 5212505, 5212600, 5212800,                     &       
     & 5213000, 5312725, 5412400, 5412600, 5412800,                     &       
     & 5412905, 5413000, 5413115, 5413200, 5413400,                     &       
     & 5413600, 5513335, 5613000, 5613200, 5613400,                     &       
     & 5613515, 5613600, 5613715, 5613800,-5713850,                     &       
     & 5713935, 5813600, 5813800, 5814000, 5814200,                     &       
     & 5914125, 6014200,-6014335, 6014400,-6014535,                     &       
     & 6014600, 6014800, 6015000, 6214400,-6214735,                     &       
     & 6214800,-6214935, 6215000, 6215200, 6215400,                     &       
     & 6315125, 6315325, 6415200, 6415400,-6415515,                     &       
     & 6415600,-6415715, 6415800, 6416000, 6515915,                     &       
     & 6615600, 6615800, 6616000, 6616125, 6616200,                     &       
     &-6616325, 6616400,-6716535, 6816200, 6816400,                     &       
     & 6816600, 6816735, 6816800, 6817000, 6916905,                     &       
     & 7016800, 7017000,-7017105, 7017200,-7017325,                     &       
     & 7017400, 7017600, 7117535,-7117670, 7217400/
      DATA (IZASP(I),I=240,287)/                                        &       
     & 7217600,-7217735, 7217800, 7217945, 7218000,                     &       
     & 7318080, 7318135, 7418000, 7418200,-7418305,                     &       
     & 7418400, 7418600, 7518525, 7518725, 7618400,                     &       
     & 7618600,-7618705, 7618800,-7618915, 7619000,                     &       
     & 7619200, 7719115, 7719315, 7819000, 7819200,                     &       
     & 7819400,-7819505, 7819600, 7819800, 7919715,                     &       
     & 8019600, 8019800,-8019905, 8020000,-8020115,                     &       
     & 8020200, 8020400, 8120305, 8120505, 8220400,                     &       
     & 8220600,-8220705, 8220800,-8320945, 9023200,                     &       
     & 9223400,-9223535, 9223800/
      REAL(KIND=R4), DIMENSION(KISO) :: ABN
      DATA (ABN(I),I=1,79)/                                             &       
     & 99.98500,  0.01500,  0.00014, 99.99986,                          &       
     &  7.50000, 92.50000,100.00000, 20.00000, 80.00000,                &       
     & 98.89000,  1.11000, 99.63400,  0.36600, 99.75800,                &       
     &  0.03800,  0.20400,100.00000, 90.51000,  0.27000,                &       
     &  9.22000,100.00000, 78.99000, 10.00000, 11.01000,                &       
     &100.00000, 92.23000,  4.67000,  3.10000,100.00000,                &       
     & 95.02000,  0.75000,  4.21300,  0.01700, 75.77000,                &       
     & 24.23000,  0.33700,  0.06300, 99.60000, 93.25800,                &       
     &  0.01200,  6.73000, 96.94000,  0.64650,  0.13000,                &       
     &  2.09000,  0.00350,  0.19000,100.00000,  8.25000,                &       
     &  7.45000, 73.70000,  5.40000,  5.20000,  0.25000,                &       
     & 99.75000,  4.35000, 83.79000,  9.50000,  2.36000,                &       
     &100.00000,  5.81000, 91.75000,  2.15000,  0.29000,                &       
     &100.00000, 68.27000, 26.10000,  1.13000,  3.59000,                &       
     &  0.91000, 69.20000, 30.80000, 48.60000, 27.90000,                &       
     &  4.10000, 18.78000,  0.62000, 60.10000, 39.90000/
      DATA (ABN(I),I=80,159)/                                           &       
     & 20.50000, 27.40000,  7.80000, 36.50000,  7.80000,                &       
     &100.00000,  0.87000,  9.03000,  7.60000, 23.50000,                &       
     & 49.80000,  9.20000, 50.69000, 49.31000,  0.35000,                &       
     &  2.25000, 11.60000, 11.50000, 57.00000, 17.30000,                &       
     & 72.17000, 27.83000,  0.56000,  9.84000,  7.00000,                &       
     & 82.60000,100.00000, 51.50000, 11.20000, 17.10000,                &       
     & 17.40000,  2.80000,100.00000, 14.80000,  9.30000,                &       
     & 15.90000, 16.70000,  9.60000, 24.10000,  9.60000,                &       
     &  5.50000,  1.90000, 12.70000, 12.60000, 17.00000,                &       
     & 31.60000, 18.70000,100.00000,  1.00000, 11.00000,                &       
     & 22.20000, 27.30000, 26.70000, 11.80000, 51.82000,                &       
     & 48.18000,  1.30000,  0.89000, 12.50000, 12.80000,                &       
     & 24.11000, 12.20000, 28.70000,  7.50000,  4.30000,                &       
     & 95.70000,  1.00000,  0.67000,  0.38000, 14.70000,                &       
     &  7.75000, 24.30000,  8.60000, 32.40000,  4.60000,                &       
     &  5.60000, 57.30000, 42.70000,  0.09100,  2.50000/
      DATA (ABN(I),I=160,239)/                                          &       
     &  0.88900,  4.62000,  7.00000, 18.70000, 31.70000,                &       
     & 34.50000,100.00000,  0.10000,  0.09000,  1.91000,                &       
     & 26.40000,  4.10000, 21.20000, 26.90000, 10.40000,                &       
     &  8.90000,100.00000,  0.11000,  0.10000,  0.42000,                &       
     &  6.59000,  7.90000, 11.20000, 71.70000,  0.08900,                &       
     & 99.91000,  0.19000,  0.25000, 88.48000, 11.08000,                &       
     &100.00000, 27.20000, 12.20000, 23.80000,  8.30000,                &       
     & 17.20000,  5.70000,  5.60000,  3.10000, 15.10000,                &       
     & 11.30000, 13.90000,  7.40000, 26.60000, 22.60000,                &       
     & 47.90000, 52.10000,  0.20000,  2.10000, 14.80000,                &       
     & 20.60000, 15.70000, 24.80000, 21.80000,100.00000,                &       
     &  0.05700,  0.10300,  2.34000, 19.00000, 25.50000,                &       
     & 24.90000, 28.10000,100.00000,  0.14000,  1.56000,                &       
     & 33.40000, 22.90000, 27.10000, 14.90000,100.00000,                &       
     &  0.14000,  3.16000, 14.40000, 21.90000, 16.20000,                &       
     & 31.60000, 12.60000, 97.40000,  2.60000,  0.16000/
      DATA (ABN(I),I=240,287)/                                          &       
     &  5.20000, 18.60000, 27.10000, 13.74000, 35.20000,                &       
     &  0.01200, 99.98800,  0.13000, 26.30000, 14.30000,                &       
     & 30.67000, 28.60000, 37.40000, 62.60000,  0.01800,                &       
     &  1.58200,  1.60000, 13.30000, 16.10000, 26.40000,                &       
     & 41.00000, 37.30000, 62.70000,  0.01300,  0.78700,                &       
     & 32.90000, 33.80000, 25.30000,  7.20000,100.00000,                &       
     &  0.15000, 10.00000, 16.85000, 23.10000, 13.20000,                &       
     & 29.80000,  6.90000, 29.50000, 70.50000,  1.40000,                &       
     & 24.10000, 22.10000, 52.40000,100.00000,100.00000,                &       
     &  0.00540,  0.72000, 99.27460/
      REAL(KIND=R4), DIMENSION(83) :: ABMN
      DATA ABMN/                                                        &       
     &   0.9993,   3.9682,   6.8814,   8.9348,  10.7171,                &       
     &  11.9079,  13.8864,  15.8616,  18.8352,  20.0060,                &       
     &  22.7923,  24.0963,  26.7498,  27.8442,  30.7077,                &       
     &  31.7889,  35.1484,  39.6045,  38.7624,  39.7339,                &       
     &  44.5697,  47.4558,  50.5039,  51.5492,  54.4661,                &       
     &  55.3663,  58.4269,  58.1838,  62.9991,  64.8352,                &       
     &  69.1240,  72.0061,  74.2780,  78.3114,  79.2176,                &       
     &  83.0801,  84.7335,  86.8643,  88.1421,  90.4393,                &       
     &  92.1083,  95.1059, 100.2018, 102.0215, 105.5157,                &       
     & 106.9417, 111.4443, 113.8317, 117.6704, 120.7120,                &       
     & 126.5237, 125.8143, 130.1621, 131.7637, 136.1503,                &       
     & 137.7122, 138.9113, 139.6972, 143.0009, 149.0584,                &       
     & 150.6576, 155.8991, 157.5601, 161.0941, 163.5135,                &       
     & 165.8231, 167.4830, 171.5337, 173.4639, 176.9567,                &       
     & 179.3935, 182.2696, 184.6071, 188.6057, 190.5648,                &       
     & 193.4042, 195.2746, 198.8767, 202.6282, 205.4369,                &       
     & 207.1851, 230.0448, 235.9841/
!
      NTHEL = 1
      NISN = 0
      IZN1 = 0
      AWRN = 0.0
      DO I=1,KISO
         IZN = IABS(IZASP(I)/100000)
         IF(IZN.LT.IZ) THEN
            IF(IZN.NE.IZN1)    NTHEL = NTHEL+1
            IZN1 = IZN
         ELSE IF(IZN.EQ.IZ) THEN
            NISN = NISN+1
            ABNIS(NISN) = ABN(I)/100.
            SPINIS(NISN) = ABS(FLOAT(MOD(IZASP(I),100))/10.)
            IAIS(NISN) = IABS(MOD(IZASP(I)/100,1000))
         ELSE
            GO TO 50
         END IF
      END DO
   50 IF(NISN.GT.0)   AWRN = ABMN(NTHEL)
!
      RETURN
      END FUNCTION AWRN
!
!***********************************************************************
!
      SUBROUTINE CKNKLD
!
!     COMPARE ABUNDANCE SPIN AND NUCLEAR WEIGHT WITH EACH ISOTOPE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
      REAL(KIND=R4), INTRINSIC :: ABS, AMOD, FLOAT
!
      INTEGER(KIND=I4) :: ICE
      INTEGER(KIND=I4) :: I,N,NN
      REAL(KIND=R4) :: ANUM,DIFF
!
!     WRITE OUT HEADER
!
      WRITE(NOUT,5)
    5 FORMAT(/5X,'CHECK ON ISOTOPE PROPERTIES'/)
      ANUM = AMOD(ZA,1000.)
!
!     NO RESONANCE PARAMETERS
!
      IF(LRP.EQ.0)   THEN
         IF(ABNS(1).NE.1.0)   THEN
            EMESS = '    THE ABUNDANCE SHOULD BE 1.0 FOR THIS MATERIAL'
            CALL ERROR_MESSAGE(0)
         END IF
         ANUM = AMOD(ZA,1000.)
         IF(ANUM.LE..9)   THEN
            DO N=1,NISN
               ZAIMAT(N) = ZA+FLOAT(IAIS(N))
            END DO
            NISMAT = NISN
         ELSE
            ZAIMAT(1) = ZA
            NISMAT = 1
         END IF
         GO TO 100
      END IF
!
!     HAS RESONANCE PARAMETERS
!
      IF(NISN.EQ.0)   THEN
         EMESS = '    ELEMENT NOT NATURALLY OCCURRING, NUCLIDE '//      &       
     &           'TESTS DISCONTINUED'
         CALL ERROR_MESSAGE(0)
         NISMAT = 0
         GO TO 100
      END IF
      IF(ANUM.GT.0..AND.ABNS(1).NE.1.0)   THEN
         EMESS = '    THE ABUNDANCE SHOULD BE 1.0 FOR THIS MATERIAL'
         CALL ERROR_MESSAGE(0)
      END IF
!
!     THE NUMBER CAN BE LESS THAN OR EQUAL TO THE NUMBER OF ISOTOPES
!       NATURALLY AVAILABLE
!
      IF(NIS.GT.NISN)   THEN
         EMESS = '    THE NUMBER OF ISOTOPES EXCEEDS THE NUMBER '//     &       
     &           'AVAILABLE'
         CALL ERROR_MESSAGE(0)
      END IF
      DO N=1,NIS
         ICE = INT(AMOD(ZAIS(N),1000.))
!********SCAN THROUGH EACH NATURALLY OCCURRING ISOTOPE
         DO I=1,NISN
            IF(ICE.EQ.IAIS(I))  THEN
!***********ISOTOPE MATCHES CHART TEST SPIN AND ABUNDANCE
               DIFF = ABS(SPINS(N)-SPINIS(I))
               IF(DIFF.GT.EPI1)   THEN
                  WRITE(EMESS,'(4X,A,F8.0,A,F4.1)')                     &       
     &                  'SPIN OF ISOTOPE ',ZAIS(N),' DISAGREES WITH ',  &       
     &                          SPINIS(I)
                  CALL ERROR_MESSAGE(0)
               END IF
               IF(NIS.NE.1)   THEN
                  IF(ABNIS(N).GT.0.0)   THEN
                     DIFF = ABS(1.-ABNS(I)/ABNIS(N))
                     IF(DIFF.GT.0.000001)    THEN
                        WRITE(EMESS,'(4X,A,F8.0,A,F10.7)')              &       
     &                       'ABUNDANCE OF ISOTOPE ',ZAIS(N),           &       
     &                       ' DISSAGREES WITH ',ABNIS(I)
                        CALL ERROR_MESSAGE(0)
                     END IF
                  END IF
               END IF
            END IF
         END DO
      END DO
!
!     SAVE MATERIAL COMPONENT DESCRIPTION
!
      NISMAT = NIS
      DO NN=1,NISMAT
         ZAIMAT(NN) = ZAIS(NN)
      END DO
!
  100 RETURN
      END SUBROUTINE CKNKLD
!
!***********************************************************************
!
      SUBROUTINE CKF3
!
!     CHECK FILE 3 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPART,JPART
      INTEGER(KIND=I4) :: MTT,MTL,NLMOD,NBEG
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: Q,QI,QMSAV
!
      INTEGER(KIND=I4), PARAMETER :: NPARTS=6
      INTEGER(KIND=I4), DIMENSION(NPARTS), PARAMETER ::                  &      
     &               ISAME=(/0,1001,1002,1003,2003,2004/)
!
!     READ THE SECTION
!
      CALL RDTAB1
!
!     SEE IF WE MUST SAVE THIS SECTION FOR LATER USE
!
      IF(NBSN.NE.NBS)   THEN
         DO N=1,NBS
            IF(MT.EQ.MTS(N))  THEN
               NBSN = NBSN + 1
               CALL RDWRSC(1,JSCR)
               GO TO 10
            END IF
         END DO
      END IF
!
!     SAVE Q VALUE FOR ENERGY CONSERVATION TEST
!
   10 Q = C1
      QI = C2
      IF(NFOR.LT.6)  THEN
         IPART = NSUB/10
         IF(MT.GE.50.AND.MT.LE.91) THEN
            JPART = 1
            MTL = MT - 50
         END IF
         IF(MT.LT.699.OR.MT.GT.799)  THEN
            NBEG = 700
            NLMOD = 20
            MTT = MT - NBEG
            MTL = MOD(MTT,NLMOD)
            JPART = ISAME((MTT/NLMOD)+2)
         END IF
         IF(MTL.LT.1) THEN
            QMSAV = Q
         ELSE IF(MTL.EQ.0) THEN
            IF(IPART.EQ.JPART)   QMSAV = ELIS
            Q = QMSAV
         ELSE
            Q = QMSAV
         END IF
      END IF
      IF(NISMAT.EQ.1) CALL QTEST(Q)
      NMT3 = NMT3 + 1
      IF(NMT3.GT.SZMT3)   THEN
         NMT3 = SZMT3
         WRITE(EMESS,'(A,I4,A)')                                        &       
     &        'CANNOT CHECK SECONDARY ENERGY CONSERVATION FOR MT = ',   &       
     &         MT,' . TOO MANY MTS.'
         CALL ERROR_MESSAGE(0)
      ELSE
         MT3(NMT3) = MT
         QVAL(NMT3) = Q
         QIVAL(NMT3) = QI
      END IF
!
!     PERFORM CAREN TESTS IF REQUIRED
!
      CALL CAREN(1)
!
      RETURN
      END SUBROUTINE CKF3
!
!***********************************************************************
!
      SUBROUTINE QTEST(Q)
!
!     CALCULATES AND COMPARES Q VALUES
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: Q
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: IPART,IZA,IZARES
      INTEGER(KIND=I4) :: MTLIB,MTALL,MTPROD
      INTEGER(KIND=I4) :: I,N,INIS
      INTEGER(KIND=I4) :: IDIV
      REAL(KIND=R4) :: QNEW,QLIB
      REAL(KIND=R4) :: DIFF
      REAL(KIND=8) :: AMASO,AMASN,AMAST,AMASI,AMASR
!
      INTEGER(KIND=I4), PARAMETER :: IDIV0=100000
      REAL(KIND=R4), PARAMETER :: EPID=.05,EPIQ=2.5E+5
!
      INTEGER(KIND=I4), PARAMETER :: NPARTS=7
      INTEGER(KIND=I4), DIMENSION(NPARTS) ::                             &      
     &                JPART=(/0,1,1001,1002,1003,2003,2004/)
      INTEGER(KIND=I4), DIMENSION(NPARTS) ::                             &      
     &                MTELAS=(/102,4,103,104,105,106,107/)
!
      INTEGER(KIND=I4), PARAMETER :: NMTS=38
      INTEGER(KIND=I4), DIMENSION(NMTS) ::                               &      
     & LIBRY = (/     4100000,  16200000,  17300000,                    &       
     &  22100001,  23100003,  24200001,  25300001,  28110000,           &       
     &  29100002,  30200002,  32101000,  33100100,  34100010,           &       
     &  35101002,  36100102,  37400000,  41210000,  42310000,           &       
     & 102000000, 103010000, 104001000, 105000100, 106000010,           &       
     & 107000001, 108000002, 109000003, 111020000, 112010001,           &       
     & 113000102, 114001002, 115011000, 116010100, 201100000,           &       
     & 203010000, 204001000, 205000100, 206000010, 207000001/)
!
!     SKIP TEST IF Q IS UNKNOWN
!
      IF(Q.EQ.QUNK)   GO TO 100
!
!     DETERMINE INCIDENT PARTICLE
!
      IPART = NSUB/10
      DO N=1,NPARTS
         IF(IPART.EQ.JPART(N))   THEN
            IF(MTELAS(N).EQ.MT)    GO TO 100
            IZAIN = JPART(N)
            GO TO 10
         END IF
      END DO
      GO TO 100
!
!     CALCULATE Q FOR EACH ISOTOPE IN THIS MATERIAL
!
   10 QNEW = -BIGNO
      IF (NISMAT.EQ.0)  GO TO 100
      DO INIS=1,NISMAT
         IZA = ZAIMAT(INIS)
         IZARES = IZA + IZAIN
!
!        LOOK UP REACTION PRODUCTS
!
         DO I=1,NMTS
            MTLIB = LIBRY(I)/1000000
            IF(MTLIB.EQ.MT)   GO TO 20
         END DO
         GO TO 100
!
!        CALCULATE AND SUM MASS OF ALL OUTGOING PARTICLES
!
   20    MTALL = MOD(LIBRY(I),1000000)
         AMASO = 0.0D+0
         IDIV = IDIV0
         DO I=1,6
            MTPROD = MOD(MTALL,10*IDIV)/IDIV
            IDIV = IDIV/10
            IF(MTPROD.NE.0)  THEN
               AMASN = AMASS(IZAB(I))*MTPROD
               AMASO = AMASO + AMASN
               IZARES = IZARES-MTPROD*IZAB(I)
            END IF
         END DO
!
!        GET MASS OF TARGET NUCLIDE
!
         AMAST = AMASS(IZA)
!
!        GET MASS OF INCIDENT PARTICLE
!
         AMASI = AMASS(IZAIN)
!
!        GET RESIDUAL MASS
!
         AMASR = AMASS(IZARES)
!
!        CALCULATE Q VALUE
!
         QLIB = (AMAST+AMASI-AMASO-AMASR)*1000.0
         QNEW = AMAX1(QNEW,QLIB)
      END DO
!
!     TEST Q VALUE
!
      IF(QNEW.EQ.0.0)  THEN
         IF(Q.NE.0.0)  THEN
            WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')                     &       
     &        'THE CALCULATED Q ',QNEW,' DISSAGREES WITH THE GIVEN Q ',Q
            CALL ERROR_MESSAGE(0)
         END IF
      ELSE
         DIFF = ABS(1.0-Q/QNEW)
         IF(DIFF.GT.EPID.AND.QNEW.GT.EPIQ)  THEN
            WRITE(EMESS,'(4X,A,1PE12.5,A,1PE12.5)')                     &       
     &        'THE CALCULATED Q ',QNEW,' DISSAGREES WITH THE GIVEN Q ',Q
            CALL ERROR_MESSAGE(0)
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE QTEST
!
!***********************************************************************
!
      SUBROUTINE CKF4
!
!     CHECK FILE 4 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), PARAMETER :: MXEN=3000
!
      INTEGER(KIND=I4) :: LVT,LTT,LI
      INTEGER(KIND=I4) :: NE,NOFF
      INTEGER(KIND=I4) :: INTT,IP,MTT
      INTEGER(KIND=I4) :: L,N,NN
      REAL(KIND=R4) :: ENN,PO,CON,CONS,YT,PERERR
      REAL(KIND=R4), DIMENSION(4) :: SIGP
      REAL(KIND=R4), DIMENSION(MXEN) :: EN,WLT,STOT,SEL
!
      REAL(KIND=R4), PARAMETER :: WLCONS=3.056E-08
!
!     SAVE FLAGS
!
      LVT = L1H
      LTT = L2H
!
!     NO TRANSFORMATION MATRIX
!
      IF(LVT.EQ.0)   THEN
         CALL RDCONT
         LI = L1H
!
!     TRANSFORMATION MATRIX GIVEN
!
      ELSE
         CALL RDLIST
         LI = L1L
      END IF
!
!     END OF DATA IF ALL DISTRIBUTIONS ISOTROPIC
!
      IF(LI.EQ.1)   GO TO 100
!
!     READ IN ENERGY INTERPOLATION SCHEME
!
      CALL RDTAB2
      NE = NP2
      NOFF = 0
      IF(NE.GT.MXEN) STOP 'PSYCHE ERROR - MXEN limit exceeded'
!
!     LEGENDRE COEFFICIENTS
!
      IF(LTT.EQ.1.OR.LTT.EQ.3)    THEN
         DO N=1,NE
            CALL RDLIST
            ENERGY = C2L
            NCOEF = NPL
            DO L=1,NCOEF
               F(L) = Y(L)
            END DO
!********CHECK THAT LEGENDRE COEFFICIENTS DO NOT PRODUCE NEGATIVE
!********ANGULAR DISTRIBUTIONS
            CALL LEGCK
!*********CALCULATE ZERO DEGREE CROSS SECTION FOR WICK'S LIMIT TEST
            IF(MT.EQ.2.AND.NSUB/10.EQ.1)   THEN
               EN(N) = ENERGY
               PO = 1.
               CON = 0.5
               DO L=1,NCOEF
                  PO = PO + CON*F(L)
                  CON = CON + 1.
               END DO
               WLT(N) = PO
            END IF
         END DO
!
!        POSSIBLE DUAL REPRESENTATION
!
         IF(LTT.EQ.1) GO TO 30
         NOFF = NE - 1
         CALL RDTAB2
         NE = NP2
      END IF
!
!     TABULAR DISTRIBUTION
!
      DO N=1,NE
         CALL RDTAB1
         IF(MT.EQ.2.AND.NSUB/10.EQ.1) THEN
            EN(NOFF+N) = C2
            WLT(NOFF+N) = Y(NP)
         END IF
      END DO
!
!     TEST WICK'S LIMIT FOR ELASTIC ANGULAR DISTRIBUTIONS
!
   30 IF(MT.EQ.2.AND.NSUB/10.EQ.1)  THEN
         DO NN=1,NE
            STOT(NN) = 0.
            SEL(NN) = 0.
         END DO
!********CALCULATE RESONANCE CONTRIBUTION IF THERE IS ONE
         IF(LRP.EQ.1)    THEN
            DO NN=1,NE
               ENN = EN(NN)
               IF(ENN.GE.E1.AND.ENN.LE.E2)   THEN
                  CALL SIGMA(ENN,SIGP,1)
                  STOT(NN) = SIGP(1)
                  SEL(NN) = SIGP(2)
               END IF
            END DO
         END IF
         CONS = WLCONS*(AWR/(AWR+1.0))**2
!
!        GET TOTAL FROM SCRATCH FILE
!
         MTT = MT
         CALL RDWRSC(2,JSCR)
         IF(MT.NE.1)   GO TO 90
!
!        CALCULATE RATIO FOR TESTING
!
         INTT = 1
         IP = 1
         DO N=1,NE
            CALL INTER (IP,INTT,EN(N),YT)
            YT = YT + STOT(N)
            WLT(N) = CONS*EN(N)*YT*YT/WLT(N)
         END DO
!
!        GET ELASTIC CROSS SECTION FROM SCRATCH
!
         CALL RDWRSC(2,JSCR)
         IF(MT.NE.2)   GO TO 90
!
!        DO TEST AT EACH ENERGY
!
         INTT = 1
         IP = 1
         DO N=1,NE
            CALL INTER(IP,INTT,EN(N),YT)
            YT = YT + SEL(N)
            IF(YT.LE.0.0) THEN
               WRITE(EMESS,'(A,1PE12.5,2A)')                            &       
     &              'WICKS LIMIT AT ',EN(N),' CANNOT BE TESTED DUE ',   &       
     &                 'ZERO ELASTIC CROSS SECTION'
               CALL ERROR_MESSAGE(0)
            ELSE
               WLT(N) = WLT(N)/YT
               IF(WLT(N).GT.1.0)   THEN
                  PERERR = 100.*(WLT(N)-1.0)
                  WRITE(EMESS,'(A,1PE12.5,A,1PE12.5,A)')                &       
     &                 'WICKS LIMIT AT ',EN(N),' EV IS EXCEEDED BY ',   &       
     &                    PERERR,' PERCENT.'
                  CALL ERROR_MESSAGE(0)
               END IF
            END IF
         END DO
         MT = MTT
         REWIND (UNIT=JSCR)
      END IF
      GO TO 100
!
!     MESSAGE IF A CROSS SECTION NEEDED IN THE TEST IS MISSING
!
   90 WRITE(NOUT,95)
   95 FORMAT(5X,'CANNOT PERFORM WICK LIMIT TEST BECAUSE TOTAL ',        &       
     &  'AND/OR ELASTIC CROSS SECTIONS ARE MISSING.')
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
      CHARACTER(LEN=1) :: IFLAG
      INTEGER(KIND=I4) :: LF
      INTEGER(KIND=I4) :: NE,NA
      INTEGER(KIND=I4) :: NK,NREGS,NPTS
      INTEGER(KIND=I4) :: IFP
      INTEGER(KIND=I4) :: INC,I,II,L,NN,NUL
      REAL(KIND=R4) :: E,U,THETA
      REAL(KIND=R4) :: GEBAR,EBAR,EBAR0,FNORM,PK,PERNK
      REAL(KIND=R4) :: EFL,EFH
      REAL(KIND=R4) :: EMULT,ECMPQ
      REAL(KIND=R4) :: AWRX
      REAL(KIND=R4) :: QI
!
!     PROCESS ALL SUBSECTIONS
!
      NK = N1H
      IF(NK.GT.1)   REWIND (UNIT=KSCR)
      DO INC=1,NK
         CALL RDTAB1
         NREGS = NR
         DO I=1,NREGS
            NTERP(I) = NBT(I)
            INTERP(I) = JNT(I)
         END DO
         NPTS = NP
            IF(NP.GT.NGRIDMAX) THEN
              WRITE(EMESS,'(A,I6,A)')                                        &       
     &        'TOO MANY ENERGIES (',NP,')--CHECK SUPRESSED'
              CALL ERROR_MESSAGE(0)
              GO TO 100
            END IF
         DO I=1,NPTS
            EINT(I) = X(I)
            PKINT(I) = Y(I)
         END DO
!
!        BRANCH ON REPRESENTATION
!
         LF = L2
!
!        LF=1
!
         IF(LF.EQ.1) THEN
            CALL RDTAB2
            NE = NP2
            IFP = 1
            IF(NE.GT.NGRIDMAX) THEN
              WRITE(EMESS,'(A,I6,A)')                                        &       
     &        'TOO MANY ENERGIES (',NE,')--CHECK SUPRESSED'
              CALL ERROR_MESSAGE(0)
              GO TO 100
            END IF
            DO I=1,NE
               CALL RDTAB1
               U = 0.
               THETA = 0.
               E = C2
               CALL EAVE(X,Y,NP,NBT,JNT,NR,E,U,THETA,LF,GEBAR,          &       
     &                FNORM,EBAR)
               CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
               AEBAR(I) = PK*EBAR
               EGRID(I) = E
            END DO
            NR = NR2
            DO I=1,NR
               NBT(I) = NBT2(I)
               JNT(I) = JNT2(I)
            END DO
!
!        LF=3
!
         ELSE IF(LF.EQ.3) THEN
            THETA = C2
            NE = NPTS
            IFP = 1
            IF(NE.GT.NGRIDMAX) THEN
              WRITE(EMESS,'(A,I6,A)')                                        &       
     &        'TOO MANY ENERGIES (',NE,')--CHECK SUPRESSED'
              CALL ERROR_MESSAGE(0)
              GO TO 100
            END IF
            DO I=1,NE
               E = EINT(I)
               CALL EAVE(X,Y,NP,NBT,JNT,NR,E,U,THETA,LF,GEBAR,          &       
     &                FNORM,EBAR)
               CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
               AEBAR(I) = PK*EBAR
               EGRID(I) = E
            END DO
!
!        LF=5
!
         ELSE IF(LF.EQ.5) THEN
            U = C1
            CALL RDTAB1
            NE = NP
            IF(NE.GT.NGRIDMAX) THEN
              WRITE(EMESS,'(A,I6,A)')                                        &       
     &        'TOO MANY ENERGIES (',NE,')--CHECK SUPRESSED'
              CALL ERROR_MESSAGE(0)
              GO TO 100
            END IF
            DO I=1,NE
               X(I+2500) = X(I)
               Y(I+2500) = Y(I)
            END DO
            CALL RDTAB1
            IFP = 1
            DO I=1,NE
               E = X(I+2500)
               THETA = Y(I+2500)
               CALL EAVE(X,Y,NP,NBT,JNT,NR,E,U,THETA,LF,GEBAR,          &       
     &                     FNORM,EBAR)
               CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
               AEBAR(I) = PK*EBAR
               EGRID(I) = E
            END DO
!
!     LF=7,9
!
         ELSE IF(LF.EQ.7.OR.LF.EQ.9) THEN
            U = C1
            CALL RDTAB1
            NE = NP
            IF(NE.GT.NGRIDMAX) THEN
              WRITE(EMESS,'(A,I6,A)')                                        &       
     &        'TOO MANY ENERGIES (',NE,')--CHECK SUPRESSED'
              CALL ERROR_MESSAGE(0)
              GO TO 100
            END IF
            IFP = 1
            DO I=1,NE
               E = X(I)
               THETA = Y(I)
               CALL EAVE(X,Y,NP,NBT,JNT,NR,E,U,THETA,LF,GEBAR,          &       
     &              FNORM,EBAR)
               CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
               AEBAR(I) = PK*EBAR
               EGRID(I) = E
            END DO
!
!        LF=11
!
         ELSE IF(LF.EQ.11) THEN
            U = C1
            CALL RDTAB1
            NE = NP
            DO I=1,NP
               X(I+2500) = X(I)
               Y(I+2500) = Y(I)
            END DO
!********STORE A(E) AND USE GRID FROM B(E)
            CALL RDTAB1
            DO I=1,NE
               X(I) = Y(I+2500)
            END DO
            IFP = 1
            DO I=1,NE
               E = X(I+2500)
               CALL EAVE(X(I),Y(I),NP,NBT,JNT,NR,E,U,THETA,LF,GEBAR,    &       
     &                FNORM,EBAR)
               CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
               AEBAR(I) = PK*EBAR
               EGRID(I) = E
            END DO
!
!        LF=12
!
         ELSE IF(LF.EQ.12) THEN
            U = C1
            CALL RDTAB1
            EFL = C1
            EFH = C2
            NE = NP
            IFP = 1
            EBAR0 = 0.5*(EFL+EFH)
            DO I=1,NE
               E = X(I)
               EBAR = EBAR0 + 1.33333*Y(I)
               CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
               AEBAR(I) = PK*EBAR
               EGRID(I) = E
            END DO
         ELSE
            GO TO 20
         END IF
!
!        SAVE ENERGY GRID AND EBAR
!
         IF(NK.GT.1) THEN
            NP = NE
            DO I=1,NP
               X(I) = EGRID(I)
               Y(I) = AEBAR(I)
            END DO
            CALL RDWRSC(1,KSCR)
         END IF
   20 END DO
!
!     CREATE A COMMON GRID FOR EBARS
!
      IF(NK.GT.1)     THEN
         REWIND (UNIT=KSCR)
!
!        GET FIRST PARTIAL AND SAVE GRID
!
         CALL RDWRSC(2,KSCR)
         NA = NP
         DO I=1,NA
            X(I+2500) = X(I)
         END DO
!
!        PROCESS REMAINING PARTIALS
!
         DO I=2,NK
            CALL RDWRSC(2,KSCR)
            NE = 15000
            CALL GMERGE(X(2501),NA,X,NP,Y(2501),NE)
            IF(I.NE.NK)   THEN
               NA = NE
               DO II=1,NA
                  X(II+2500) = Y(II+2500)
               END DO
            END IF
         END DO
         DO NUL=1,NE
            AEBAR(NUL) = 0.0
            EGRID(NUL) = Y(NUL+2500)
         END DO
         REWIND (UNIT=KSCR)
!
!        MERGE EBARS ON THE COMMON GRID
!
         DO I=1,NK
            CALL RDWRSC(2,KSCR)
            IFP = 1
            DO L=1,NE
               CALL TERPR(EGRID(L),X,Y,NP,NBT,JNT,NR,IFP,PERNK)
               AEBAR(L) = AEBAR(L) + PERNK
            END DO
         END DO
         REWIND (UNIT=KSCR)
      END IF
!
!     LOOK FOR Q VALUE FOR THIS MT
!
      QI = 0.
      DO I=1,NMT3
         IF(MT3(I).EQ.MT)   THEN
            QI = QIVAL(I)
            GO TO 30
         END IF
      END DO
!
!     CALCULATE TOTAL ENERGY AND TEST AGAINST AVAILABLE ENERGY
!
   30 IF(NSUB.EQ.4)  GO TO 100
      IF(LNU.LT.0) THEN
        WRITE(NOUT,34)
        GO TO 40
      END IF
      WRITE(NOUT,35)
   34 FORMAT(/4X,'WARNING - NU-BAR DATA UNDEFINED')
   35 FORMAT(/4X,'ENERGY GRID     TOTAL SECONDARY      ',               &       
     &  'AVERAGE SECONDARY        ENERGY'/                              &       
     &  25X,'ENERGY            NEUTRON ENERGY',8X,'AVAILABLE'/          &       
     &  4X,11('-'),6X,14('-'),6X,17('-'),6X,10('-'))
      AWRX = AWR/(AWR+1.)
      DO I=1,NE
         E = EGRID(I)
         X(I) = E
         ECMPQ = AMAX1(E*AWRX+QI,0.0)
         EBAR = AEBAR(I)
         EMULT = EBAR*ANU(E)
!********SAVE TOTAL SECONDARY ENERGY
         Y(I) = EMULT
!********FLAG ALL TOTAL SECONDARIES G.T. ECM+QI
         IF(EMULT.GT.ECMPQ)  THEN
            IFLAG = '*'
         ELSE
            IFLAG = ' '
         END IF
         WRITE(NOUT,'(4X,1PE11.4,6X,1PE11.4,A1,10X,1PE11.4,9X,1PE11.4)')&       
     &            E,EMULT,IFLAG,EBAR,ECMPQ
      END DO
      WRITE(NOUT,'(/4X,A/)')                                            &       
     &        '*  AVERAGE ENERGY EXCEEDS AVAILABLE ENERGY'
!
!     STORE FILE NAME AND REACTION IF THERE IS PHOTON PRODUCTION DATA
!
   40 DO NN=1,MTTOT
         IF(MT.EQ.MTDCT(NN))  THEN
            MF5 = MF5+1
            MT5(MF5) = MF*1000+MT
            NR = 1
            NP = NE
            NBT(1) = NP
            JNT(1) = 2
            CALL RDWRSC(1,ISCR)
            DO I=1,NP
               Y(I) = AEBAR(I)
            END DO
            CALL RDWRSC(1,ISCR)
            GO TO 100
         END IF
      END DO
!
  100 RETURN
      END SUBROUTINE CKF5
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION ANU(E)
!
!     FUNCTION TO RETURN NEUTRON MULTIPLICITY
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: E
!
      INTEGER(KIND=I4) :: IFP
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: EE
!
      INTEGER(KIND=I4),PARAMETER :: NMUMTS=8
      INTEGER(KIND=I4), DIMENSION(NMUMTS) ::                             &      
     &         MUMTS=(/16,17,24,25,30,37,41,42/)
      REAL(KIND=R4), DIMENSION(NMUMTS) ::                                &      
     &         AMULT=(/2.,3.,2.,3.,2.,4.,2.,3./)
!
!     INITIALIZE RETURN VARIABLE
!
      ANU = 1.
!
!     FISSION NEUTRONS
!
      IF((MT.GE.18.AND.MT.LE.21).OR.MT.EQ.38)  THEN
         IF(LNU.EQ.1)   THEN
            ANU = POLI(1)
            EE = E
            DO I=2,NRC
               ANU = ANU + POLI(I)*EE
               EE = E*EE
            END DO
         ELSE IF(LNU.EQ.2) THEN
            IFP = 1
            CALL TERPR(E,ETERP,ANUTRP,NPT,ITERP,JTERP,NRC,IFP,ANU)
         ELSE
!           -- This should not happen (protected in the calling program)
            STOP 'LNU undefined'
         END IF
!
!     ALL BUT FISSION
!
      ELSE
         DO I=1,NMUMTS
            IF(MT.EQ.MUMTS(I))   THEN
               ANU = AMULT(I)
               GO TO 100
            END IF
         END DO
      END IF
!
  100 RETURN
      END FUNCTION ANU
!
!***********************************************************************
!
      SUBROUTINE CKF6
!
!     CHECK FILE 6 DATA
!                        THIS ROUTINE PROVIDED BY R MACFARLANE LANL
!                            JUNE 1985
!                        REVISED TO HANDLE LAW=7 BY C.L. DUNFORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: LCT,LANG,LAW,LTP
      INTEGER(KIND=I4) :: NK,NREGS,NPTS
      INTEGER(KIND=I4) :: IZAP,IZ,IA
      INTEGER(KIND=I4) :: NW,ND
      INTEGER(KIND=I4) :: LL,NOG,NOC,MM,NU
      INTEGER(KIND=I4) :: INT1,INT2,NE,IFP,NCYC,L,IST,IDISCR,IBZ,NMU
      INTEGER(KIND=I4) :: NSEQPN
      INTEGER(KIND=I4) :: I,II,IK,IE,IU,N,NM
      REAL(KIND=R4) :: E,EMIN,EMAX,ENOW,ENEXT
      REAL(KIND=R4) :: Q,QI,QK
      REAL(KIND=R4) :: PK,XY,YN,X3,Y3,X4,Y4,EBAR,GEBAR,FNORM
      REAL(KIND=R4) :: TEST,YIBZ,SSUM,ANS,EAVAIL,PERR
      REAL(KIND=R4) :: AWRS,AWIN
!
      INTEGER(KIND=I4), PARAMETER :: NPOUTMAX=400
      INTEGER(KIND=I4), PARAMETER :: NMUMAX=201
      CHARACTER(LEN=10), DIMENSION(NPOUTMAX) :: HL
      INTEGER(KIND=I4), DIMENSION(NPOUTMAX) :: LIK,JIK,NIK
      REAL(KIND=R4),DIMENSION(NPOUTMAX) :: VAL
      REAL(KIND=R4), DIMENSION(NMUMAX) :: XX,YY,ZZ
!
      CHARACTER(LEN=10), PARAMETER :: HE='      E   '
      CHARACTER(LEN=10), PARAMETER :: HS='     SUM  '
      CHARACTER(LEN=10), PARAMETER :: HA='    AVAIL '
      CHARACTER(LEN=7), PARAMETER :: HP='  %DIFF'
!
!     FIND Q FOR THIS REACTION
!
      Q = 0.
      QI= 0.
      DO I=1,NMT3
         IF(MT3(I).EQ.MT) THEN
            Q = QVAL(I)
            QI = QIVAL(I)
         END IF
      END DO
!
!     PROCESS ALL SUBSECTIONS
!
      LL = 0
      LCT = L2H
      NK = N1H
      IF(NK.GT.NPOUTMAX) STOP 'PSYCHE ERROR - NPOUTMAX limit exceeded'
      EMIN = 1.E+10
      EMAX = 0.
      NOG = 1
      NOC = 0
      DO IK=1,NK
         CALL RDTAB1
         IZAP = IFIX(C1)
         IF(IK.LE.NPOUTMAX) WRITE(HL(IK),'(I8.5,A2)')  IZAP,'  '
         NIK(IK) = 0
         NREGS = NR
         DO I=1,NREGS
            NTERP(I) = NBT(I)
            INTERP(I) = JNT(I)
         END DO
         NPTS = NP
         DO I=1,NPTS
            EINT(I) = X(I)
            PKINT(I) = Y(I)
         END DO
!
!        BRANCH ON REPRESENTATION
!
         LAW = L2
!
!        LAW 1 -- TABULATED ENERGY-ANGLE DISTRIBUTION
!
         IF(LAW.EQ.1) THEN
            CALL RDTAB2
            IF(IZAP.EQ.0) NOG = 0
            LANG = L12
            INT2 = L22
            NE = NP2
            IFP = 1
            IF(IK.LE.NPOUTMAX) THEN
               LIK(IK) = LL
               NIK(IK) = NE
            END IF
            DO IE=1,NE
               CALL RDLIST
               E = C2L
               IF(E.LT.EMIN) EMIN = E
               IF(E.GT.EMAX) EMAX = E
               ND = L1L
            NCOEF = L2L
            NP = N2L
            NW = NPL
            CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
            XY = 0.
            YN = 0.
            NCYC = NW/NP
            IF(ND.GT.0) THEN
!
!              Discrete lines
!
               DO I=1,ND
                  L = NCYC*(I-1) + 1
                  XY = XY + Y(L)*Y(L+1)
                  YN = YN + Y(L+1)
               END DO
            END IF
            IF(ND.NE.NP) THEN
!
!              Spectrum
!
               IST = 2 + ND
               DO I=IST,NP
                  L = NCYC*(I-2) + 1
                  X3 = Y(L)
                  Y3 = Y(L+1)
                  X4 = Y(L+NCYC)
                  Y4 = Y(L+NCYC+1)
                  CALL XECSI(X3,Y3,X4,Y4,X3,X4,INT2,ANS)
                  XY = XY + ANS
                  CALL ECSI(X3,Y3,X4,Y4,X3,X4,INT2,ANS)
                  YN = YN + ANS
               END DO
            END IF
            IF(YN.NE.0.) THEN
               EBAR = XY/YN
            ELSE
               EBAR = 0.0
            END IF
            IF(MT.GE.51.AND.MT.LE.90) THEN
               IDISCR = 1
            ELSE IF(MT.GE.600.AND.MT.LE.849) THEN
               IF(MOD(MT,50).EQ.49)  THEN
                  IDISCR = 0
               ELSE
                  IDISCR = 1
               END IF
            ELSE
               IDISCR = 0
            END IF
            IF(IDISCR.EQ.1.AND.IZAP.EQ.0) THEN
               TEST = -Q + QI + EBAR
               IF (TEST.GE.EPI3) THEN
                  WRITE(EMESS,'(3(A,1PE11.4))')                         &       
     &                 'GAMMA ERROR  E=',E,'  EBAR=',EBAR,'  QI=',QI
                  CALL ERROR_MESSAGE(0)
               END IF
               NOC = 1
            ELSE
               LL = LL + 1
               IF (LL.LE.NGRIDMAX) THEN
                  EGRID(LL) = E
                  AEBAR(LL) = PK*EBAR
               END IF
            END IF
            IF(LANG.EQ.1.AND.NCOEF.NE.0) THEN
               ENERGY = E
               NSEQPN = NSEQP1
               DO II=1,NP
                  IBZ = NCYC*(II-1) + 2
                  YIBZ = Y(IBZ)
                  DO I=1,NCOEF
                     F(I) = Y(IBZ+I)
                     IF (Y(IBZ).NE.0.) F(I) = F(I)/YIBZ
                  END DO
                  NSEQP1 = NSEQPN + (IBZ-2)/6
                  CALL LEGCK
                END DO
            END IF
         END DO
!
!        LAW 2 -- TWO-BODY ANGULAR DISTRIBUTION
!
         ELSE IF(LAW.EQ.2) THEN
            CALL RDTAB2
            NE = NP2
            DO IE=1,NE
               CALL RDLIST
               ENERGY = C2L
               LTP = L1L
               NCOEF = N2L
               IF(LTP.EQ.0.AND.NCOEF.NE.0) THEN
                  DO I=1,NCOEF
                     F(I)=Y(I)
                  END DO
                  CALL LEGCK
               END IF
            END DO
            EMESS = 'NO ENERGY-BALANCE TEST FOR TWO-BODY LAW'
            CALL ERROR_MESSAGE(0)
!
!        LAW 5 -- CHARGED-PARTICLE ELASTIC SCATTERING
!
         ELSE IF(LAW.EQ.5) THEN
            CALL RDTAB2
            NE = NP2
            DO IE=1,NE
               CALL RDLIST
            END DO
            NOC = 1
            EMESS = 'NO ENERGY-BALANCE TEST FOR CHARGED-PARTICLE '//    &       
     &              'ELASTIC LAW'
            CALL ERROR_MESSAGE(0)
!
!        LAW 6 -- PHASE-SPACE DISTRIBUTION
!
         ELSE IF(LAW.EQ.6) THEN
            CALL RDCONT
            NOC = 1
            EMESS = 'NO ENERGY-BALANCE TEST FOR PHASE-SPACE LAW'
            CALL ERROR_MESSAGE(0)
!
!        LAW 7 -- LABORATORY ANGLE-ENERGY DISTRIBUTION
!
         ELSE IF(LAW.EQ.7) THEN
            CALL RDTAB2
            IF(IZAP.EQ.0) NOG = 0
            NE = NP2
            IFP = 1
            IF(IK.LE.NPOUTMAX) THEN
               LIK(IK) = LL
               NIK(IK) = NE
            END IF
            NE = NP2
            DO N=1,NE
               CALL RDTAB2
               NMU = NP2
               IF(NMU.GT.NMUMAX)
     &            STOP 'PSYCHE ERROR - NMUMAX limit exceeded'
               E = C22
               IF(E.LT.EMIN) EMIN=E
               IF(E.GT.EMAX) EMAX=E
               CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,PK)
               DO NM=1,NMU
                  CALL RDTAB1
                  XX(NM) = C2
                  CALL EAVE(X,Y,NP,NBT,JNT,NR,0.0,0.0,0.0,1,GEBAR,      &       
     &                 FNORM,EBAR)
                  YY(NM) = FNORM
                  ZZ(NM) = GEBAR
               END DO
               GEBAR = 0.
               FNORM = 0.
               MM = 1
               DO NM=2,NMU
                  IF(NM.GT.NBT2(MM))    THEN
                     MM = MM + 1
                     IF(MM.GT.NR2)   GO TO 30
                  END IF
                  II = JNT2(MM)
                  CALL ECSI(XX(NM-1),YY(NM-1),XX(NM),YY(NM),XX(NM-1),   &       
     &                      XX(NM),II,YN)
                  FNORM = FNORM + YN
                  CALL ECSI(XX(NM-1),ZZ(NM-1),XX(NM),ZZ(NM),XX(NM-1),   &       
     &                      XX(NM),II,YN)
                  GEBAR = GEBAR + YN
               END DO
               IF(FNORM.NE.0.) THEN
                  EBAR = GEBAR/FNORM
               ELSE
                  EBAR = 0.
               END IF
               IF ((MT.GE.51.AND.MT.LE.90).AND.IZAP.EQ.0) THEN
                  TEST=ABS(EBAR+Q)
                  IF (TEST.GT.EPI3*EBAR) THEN
                     WRITE(EMESS,'(3(A,1PE11.4))')                      &       
     &                    'GAMMA ERROR  E=',E,'  EBAR=',EBAR,'  QI=',QI
                     CALL ERROR_MESSAGE(0)
                  END IF
                  NOC = 1
               ELSE
                  LL = LL + 1
                  IF(LL.LE.NGRIDMAX) THEN
                     EGRID(LL) = E
                     AEBAR(LL) = PK*EBAR
                  END IF
               END IF
   30       END DO
         END IF
!
!     END OF LOOP OVER SUBSECTIONS
!
      END DO
      IF(MT.LE.4) GO TO 90
      IF(NOG.EQ.1)  THEN
         WRITE(EMESS,'(A)') ' NO GAMMAS'
         CALL ERROR_MESSAGE(0)
      END IF
      IF(NOC.EQ.1.OR.EMIN.EQ.1.0E+10)   GO TO 90
      IF (NK.GT.NPOUTMAX) THEN
         WRITE(EMESS,'(A,I4,A)')                                        &       
     &           'TOO MANY SUBSECTIONS, NK = ',NK,                      &       
     &           ' --UNION CHECK SUPRESSED'
         CALL ERROR_MESSAGE(0)
         GO TO 100
      END IF
      IF (LL.GT.NGRIDMAX) THEN
         WRITE(EMESS,'(A,I6,A)')                                        &       
     &        'TOO MANY ENERGIES (',LL,')--UNION CHECK SUPRESSED'
         CALL ERROR_MESSAGE(0)
         GO TO 100
      END IF
!
!     GENERATE UNION ENERGY GRID
!
      ENOW = EMIN
      NU = 1
      X(NU) = EMIN
   35 ENEXT = EMAX
      DO IK=1,NK
         L = LIK(IK)
         NE = NIK(IK)
         DO IE=1,NE
            E = EGRID(L+IE)
            IF(E.GT.ENOW.AND.E.LT.ENEXT) ENEXT = E
         END DO
      END DO
      NU = NU + 1
      X(NU) = ENEXT
      ENOW = ENEXT
      IF(ENEXT.LT.EMAX) GO TO 35
!
!     CHECK FOR ENERGY BALANCE ON UNION GRID
!
      INT1 = 2
      DO IK=1,NK
         JIK(IK) = 2
      END DO
      IF(LCT.EQ.1) THEN
        WRITE(NOUT,40) Q,'(Lab)'
      ELSE IF(LCT.EQ.2 .OR. LCT.EQ.3) THEN
        WRITE(NOUT,40) Q,'(CM) '
      ELSE
        WRITE(NOUT,'(9X,A,I3,A)') 'Unknown coordinate frame',LCT
     &                           ,' (treated as CM)'
        WRITE(NOUT,40) Q,'(CM) '
        LCT=3
      END IF
   40 FORMAT(/9X,'ENERGY BALANCE SUMMARY: Q = ',1PE13.5//               &       
     &     20X,'TOTAL SECONDARY ENERGY BY EMITTED PARTICLE ',A5)
      IF (NK.GT.1) THEN
         WRITE(NOUT,'(2A10,A7,10A10:/1(37X, 9A10))')                    &
     &         HE,HA,HP,HS,(HL(I),I=1,NK)
         IF(NK.GT.12) WRITE(NOUT,'(A)') ' '
      ELSE IF (NK.EQ.1) THEN
         WRITE(NOUT,'(2A10,A7,(10A10))')HE,HA,HP,HL(1)
      END IF
      DO IU=1,NU
         E = X(IU)
         SSUM = 0.
         DO IK=1,NK
            VAL(IK) = 0.
            IF(NIK(IK).GT.0)   THEN
   45          L = LIK(IK) + JIK(IK)
               IF(E.LT.EGRID(L).OR.JIK(IK).EQ.NIK(IK)) GO TO 50
               JIK(IK) = JIK(IK) + 1
               GO TO 45
   50          CALL TERP1(EGRID(L-1),AEBAR(L-1),EGRID(L),AEBAR(L),      &       
     &              E,ANS,INT1)
               READ(HL(IK),'(I8)') IZAP
               IF(LCT.EQ.3 .AND. IZAP.GT.2004) THEN
!
!                Remove CM energy for recoils if LCT=3 to convert to CM
!                Assume: - masses approximated by mass numbers
!                        - ejectiles are emitted in all directions so
!                          that the recoil travels approximately in CM
!
                 IZ  =IZAP/1000
                 IA  =IZAP-1000*IZ
                 IF(IA.LE.0) THEN
                   AWRS=IZ
                   IF(IZ.GT.1) AWRS=AWRS*2.1
                 ELSE
                   AWRS=IA
                 END IF
                 AWIN=IZAIN-1000*(IZAIN/1000)
                 ANS=ANS*(1-AWRS/(AWRS+AWIN))
               END IF
               VAL(IK) = ANS
               SSUM = SSUM + ANS
            END IF
         END DO
         IF(LCT.EQ.1) THEN
C           Work in Lab coordinate system for LCT=1
            QK = E
         ELSE
C           Work in CM coordinate system for LCT=2,3
            QK = E*AWR/(AWR+1.)
         END IF
         EAVAIL = Q + QK
         IF(EAVAIL.LE.0.) THEN
            EAVAIL= 0.
            IF(SSUM.EQ.0.)  THEN
               PERR = 0.
            ELSE
               PERR = 999.99
            END IF
         ELSE
            PERR = 100.*(SSUM-EAVAIL)/EAVAIL
            IF(PERR.GT.999.99) PERR=999.99
         END IF
         IF(MT.EQ.5)   THEN
            IF(NK.EQ.1) THEN
               WRITE(NOUT,'(1P,E10.2,17X,(10E10.2))')  E,VAL(1)
            ELSE
               WRITE(NOUT,'(1P,E10.2,17X,10E10.2:/1(37X, 9E10.2))')     &       
     &                E,SSUM,(VAL(I),I=1,NK)
               IF(NK.GT.12) WRITE(NOUT,'(A)') ' '
            END IF
         ELSE
            IF(ABS(PERR).LT.999.99) THEN
               IF(NK.EQ.1) THEN
                  WRITE(NOUT,'(1P,2E10.2,0P,F7.2,1P,10E10.2:/
     &                        (37X, 9E10.2))')                          &       
     &                     E,EAVAIL,PERR,VAL(1)
               ELSE
                  WRITE(NOUT,'(1P,2E10.2,0P,F7.2,1P,10E10.2:/
     &                        (37X, 9E10.2))')                          &       
     &                     E,EAVAIL,PERR,SSUM,(VAL(I),I=1,NK)
                  IF(NK.GT.12) WRITE(NOUT,'(A)') ' '
               END IF
            ELSE
               IF(NK.EQ.1) THEN
                  WRITE(NOUT,'(2(1PE10.2),A,1PE10.2)')                  &       
     &                   E,EAVAIL,' ******',VAL(1)
               ELSE
                  WRITE(NOUT,'(2(1PE10.2),A,(10(1PE10.2)))')            &       
     &                   E,EAVAIL,' ******',SSUM,(VAL(I),I=1,NK)
                  IF(NK.GT.12) WRITE(NOUT,'(A)') ' '
               END IF
            END IF
         END IF
      END DO
!
   90 WRITE(NOUT,'(A)')  ' '
!
  100 RETURN
      END SUBROUTINE CKF6
!
!***********************************************************************
!
      SUBROUTINE CKF12
!
!     PROCESS FILE 12 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LO,LF,LP
      INTEGER(KIND=I4) :: NE
      INTEGER(KIND=I4) :: NK
      INTEGER(KIND=I4) :: IFP
      INTEGER(KIND=I4) :: I,N,NTH
      REAL(KIND=R4) :: EGK,AWRX,E,ECM,EBAR,YIELD
!
      LO = L1H
!
!     DO NOT PROCESS TRANSITION PROBABILITIES
!
      IF(LO.NE.1)   THEN
         CALL RDLIST
         GO TO 100
      END IF
!
!     PHOTON MULTIPLICITIES
!
      NK = N1H
      CALL RDTAB1
      EGK = C1
      LF = L2
!*****STORE COMMON GRID IN EGRID
      NE = NP
      DO I=1,NE
         EGRID(I) = X(I)
      END DO
!
!     CHECK IF ONLY SUBSECTION IS A CONTINUUM
!
      IF(NK.EQ.1.AND.EGK.EQ.0.)   GO TO 40
!*****STORE TOTAL YIELDS ON LSCR
      CALL RDWRSC(1,LSCR)
!
!     PROCESS EACH SUBSECTION
!
      AWRX = AWR/(AWR+1.)
      DO I=1,NE
         EBAV(I) = 0.0
      END DO
      DO N=1,NK
         IF(NK.GT.1)   CALL RDTAB1
         LP = L1
         LF = L2
         EGK = C1
!********CONTINUUM MUST BE LAST SUBSECTION
         IF(EGK.EQ.0.)   GO TO 40
!
!        CALCULATE EBAR AT EACH POINT AND KEEP RUNNING SUM
!
         IFP = 1
         DO NTH=1,NE
            E = EGRID(NTH)
            CALL TERPR(E,X,Y,NP,NBT,JNT,NR,IFP,YIELD)
            ECM = E*AWRX
            IF(LP.GT.1)  THEN
               EBAR = EGK + ECM
            ELSE
               EBAR = EGK
            END IF
            EBAV(NTH) = EBAV(NTH)+EBAR*YIELD
         END DO
      END DO
      GO TO 60
!
!     STORE CONTINUUM YIELDS ON KSCR FOR PROCESSING IN FILE 15
!
   40 IF(LF.EQ.1)   THEN
         CALL RDWRSC(1,KSCR)
         MT1213 = MT1213 + 1
         MTSCR(MT1213) = MT + MF*1000
      END IF
      IF(NK.EQ.1)  GO TO 100
!
!     STORE NU*EBAR ON LSCR
!
   60 DO I=1,NE
         X(I) = EGRID(I)
         Y(I) = EBAV(I)
      END DO
      NP = NE
      NR = 1
      NBT(1) = NP
      JNT(1) = 2
      NMTBAR = NMTBAR + 1
      MTBAR(NMTBAR) = MT
      CALL RDWRSC(1,LSCR)
!
  100 RETURN
      END SUBROUTINE CKF12
!
!***********************************************************************
!
      SUBROUTINE CKF13
!
!     PROCESS FILE 13 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), PARAMETER :: MXEN=3000
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: MTNOW,NK,LP,LF,NE
      INTEGER(KIND=I4) :: NX,IFP,NSPT
      INTEGER(KIND=I4) :: I,N,NTH
      REAL(KIND=R4) :: EGK,AWRX,E,ECM,EBAR,YIELD,TERP,ESIG,SIG3N
      REAL(KIND=R4), DIMENSION(MXEN) :: SIG3
!
!     SAVE CURRENT MT VALUE
!
      MTNOW = MT
!
!     READ FIRST SUBSECTION
!
      NK = N1H
      CALL RDTAB1
      EGK = C1
      LF = L2
      NE = NP
!
      IF(NE.GT.MXEN) STOP 'PSYCHE ERROR - MXEN limit exceeded'
!
!     SAVE GRID AND TOTAL YIELD
!
      NX = JNT(1)
      DO I=1,NE
         EGRID(I) = X(I)
         YLD(I) = Y(I)
      END DO
!
!     RETRIEVE CROSS SECTION OF FILE 3 FROM SCRATCH
!
      REWIND (UNIT=JSCR)
      DO I=1,NBSN
         CALL RDWRSC(2,JSCR)
!********INTERPOLATE FILE 3 CROSS SECTIONS ONTO EGRID
         IF(MT.EQ.MTNOW)   THEN
            IFP = 1
            DO N=1,NE
               CALL TERPR(EGRID(N),X,Y,NP,NBT,JNT,NR,IFP,SIG3(N))
            END DO
            GO TO 20
         END IF
      END DO
!
!     NO FILE 3 CROSS SECTION
!
      MT = - MTNOW
      DO I=1,NE
         SIG3(I) = YLD(I)
      END DO
!
!     CALCULATE TOTAL YIELD AND SAVE
!
   20 DO I=1,NE
         X(I) = EGRID(I)
         IF(SIG3(I).NE.0.)   THEN
            Y(I) = YLD(I)/SIG3(I)
         ELSE
            Y(I) = 0.0
         END IF
      END DO
      NR = 1
      NP = NE
!
!     CHECK IF ONLY SUBSECTION IS A CONTINUUM
!
      IF(EGK.EQ.0.0.AND.NK.EQ.1)   GO TO 50
!*****STORE TOTAL YIELDS ON LSCR
      JNT(1) = NX
      NBT(1) = NE
      CALL RDWRSC(1,LSCR)
!
!     PUT YIELDS INTO Y ARRAY
!
      DO I=1,NE
         Y(I) = YLD(I)
      END DO
!
!     PROCESS EACH SUBSECTION
!
      AWRX = AWR/(AWR+1.)
      DO I=1,NE
        EBAV(I) = 0.
      END DO
      DO N=1,NK
         IF(NK.GT.1)  CALL RDTAB1
         LP = L1
         LF = L2
         EGK = C1
         IF(EGK.EQ.0.)   GO TO 30
!
!        CALCULATE EBAR AT EACH POINT AND KEEP RUNNING SUM
!
         IFP = 1
         DO NTH=1,NE
            E = EGRID(NTH)
            IF(SIG3(NTH).GT.0.0)   THEN
               CALL TERPR(E,X,Y,NP,NBT,JNT,NR,IFP,TERP)
               YIELD = TERP/SIG3(NTH)
               ECM = E*AWRX
               IF(LP.GT.1)   THEN
                  EBAR = EGK + ECM
               ELSE
                  EBAR = EGK
               END IF
               EBAV(NTH) = EBAV(NTH) + EBAR*YIELD
            END IF
         END DO
      END DO
      GO TO 70
!
!     STORE CONTINUUM YIELDS ON KSCR FOR PROCESSING IN FILE 15
!
   30 IF(LF.NE.1)   GO TO 60
      NSPT = 1
      DO N=1,NP
         E = X(N)
   40    ESIG = EGRID(NSPT)
         IF(ABS((E-ESIG)/ESIG).GE.EPI4) THEN
            IF(E.GE.ESIG)   THEN
               NSPT = NSPT + 1
               IF(NSPT.LE.NE)   GO TO 40
            ELSE
               Y(N) = 0.0
            END IF
         ELSE
            SIG3N = SIG3(NSPT)
            IF(SIG3N.EQ.0.0)   THEN
               Y(N) = 0.0
            ELSE
               Y(N) = Y(N)/SIG3(NSPT)
            END IF
         END IF
      END DO
   50 JNT(1) = 2
      NBT(1) = NP
      CALL RDWRSC(1,KSCR)
      MT1213 = MT1213 + 1
      MTSCR(MT1213) = MTNOW + MF*1000
   60 IF(NK.EQ.1)   GO TO 80
!
!     STORE NU*EBAR ON LSCR
!
   70 DO I=1,NE
         X(I) = EGRID(I)
         Y(I) = EBAV(I)
      END DO
      NP = NE
      NR = 1
      NBT(1) = NP
      JNT(1) = 2
      NMTBAR = NMTBAR + 1
      MTBAR(NMTBAR) = MTNOW
      CALL RDWRSC(1,LSCR)
!
!     RESTORE CURRENT VALUE OF MT
!
   80 MT = MTNOW
!
      RETURN
      END SUBROUTINE CKF13
!
!***********************************************************************
!
      SUBROUTINE CKF14
!
!     CHECK FILE 14
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NK,LI,NI,LTT,NE
      INTEGER(KIND=I4) :: L,M,N
!
      NK = N1H
      LI = L1H
!
!     ALL DONE IF ISOTROPIC
!
      IF(LI.NE.1) THEN
         NI = N2H
         LTT = L2H
         DO N=1,NK
            IF(N.LE.NI) THEN
               CALL RDCONT
            ELSE
               CALL RDTAB2
               NE = NP2
               NR = NR2
               DO M=1,NE
                  IF(LTT.LE.1)  THEN
                     CALL RDLIST
                     ENERGY = C2L
                     NCOEF = NPL
                     DO L=1,NCOEF
                        F(L) = Y(L)
                     END DO
                     CALL LEGCK
                  ELSE
                     CALL RDTAB1
                  END IF
               END DO
            END IF
         END DO
      END IF
!
      RETURN
      END SUBROUTINE CKF14
!
!***********************************************************************
!
      SUBROUTINE CKF15
!
!     PROCESS FILE 15 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: IABS, MOD
!
      INTEGER(KIND=I4) :: MTNOW,MTT
      INTEGER(KIND=I4) :: ICONTG,IFP
      INTEGER(KIND=I4) :: NE,NTG,NPTS,NREGS
      INTEGER(KIND=I4) :: I,J,NGAM,NGAM1
      REAL(KIND=R4) :: U,THETA
      REAL(KIND=R4) :: E,GEBAR,FNORM,TERP,ANUE
      REAL(KIND=R4), DIMENSION(3000) ::  XT
!
      IERX = 0
!
!     SAVE CURRENT MT
!
      MTNOW = MT
!
!     SEE IF WE HAVE CONTINUUM YIELDS IN STORAGE
!
      DO J=1,MT1213
         MTT = MOD(MTSCR(J),1000)
         IF(MTT.EQ.MTNOW)   GO TO 20
      END DO
!
!     NO CONTINUUM YIELDS SKIP FILE 15 PROCESSING
!
      IERX = 1
      GO TO 100
!
!     PROCESS CONTINUUM SECTION
!
   20 CALL RDTAB1
      CALL RDTAB2
      NE = NP2
!
!     PROCESS ALL ENERGY POINTS
!
      DO I=1,NE
         CALL RDTAB1
         U = 0.
         THETA = 0.
         EGRID(I) = C2
         CALL EAVE(X,Y,NP,NBT,JNT,NR,C2,U,THETA,1,GEBAR,FNORM,AEBAR(I))
      END DO
!
!     SEE IF A TOTAL YIELD GRID EXISTS FROM DISCRETE PROCESSING
!
      ICONTG = 0
      DO NGAM=1,NMTBAR
         IF(MTNOW.EQ.MTBAR(NGAM))   THEN
!********RETRIEVE TOTAL YIELD GRID
            REWIND (UNIT=LSCR)
            DO NGAM1=1,NMTBAR
               CALL RDWRSC(2,LSCR)
               IF(MTNOW.EQ.MTBAR(NGAM1))   GO TO 30
               CALL RDWRSC(2,LSCR)
            END DO
   30       NTG = NP
            DO I=1,NTG
               XT(I) = X(I)
            END DO
            GO TO 40
         END IF
      END DO
      ICONTG = 1
!
!     GET YIELD FROM 12/13 SCRATCH FILE
!
   40 REWIND (UNIT=KSCR)
      DO J=1,MT1213
         CALL RDWRSC(2,KSCR)
         IF(IABS(MT).EQ.MTNOW)    GO TO 50
      END DO
   50 NPTS = NP
      DO I=1,NPTS
         EINT(I) = X(I)
         PKINT(I) = Y(I)
      END DO
      NREGS = NR
      DO I=1,NREGS
         NTERP(I) = NBT(I)
         INTERP(I) = JNT(I)
      END DO
      IF(ICONTG.GT.0)   THEN
         NTG = NP
         DO I=1,NTG
            XT(I) = EINT(I)
         END DO
      END IF
!
!     MOVE NUBARS FOR INTERPOLATION
!
      NR = NR2
      DO I=1,NR
         NBT(I) = NBT2(I)
         JNT(I) = JNT2(I)
      END DO
      NP = NE
      DO I=1,NP
         X(I) = EGRID(I)
         Y(I) = AEBAR(I)
      END DO
!
!     INTERPOLATE AVERAGE ENERGIES TO TOTAL GRID
!
      IFP = 1
      DO NGAM=1,NTG
         E = XT(NGAM)
         EGRID(NGAM) = E
         CALL TERPR(E,X,Y,NP,NBT,JNT,NR,IFP,TERP)
         AEBAR(NGAM) = TERP
      END DO
!
!     CALCULATE TOTAL ENERGIES
!
      NE = NTG
      IFP = 1
      DO I=1,NE
         E = EGRID(I)
         X(I) = E
         CALL TERPR(E,EINT,PKINT,NPTS,NTERP,INTERP,NREGS,IFP,ANUE)
         Y(I) = AEBAR(I)*ANUE
      END DO
!
!     STORE GAMMA TOTAL ENERGIES
!
      NR = 1
      NP = NE
      NBT(1) = NE
      JNT(1) = 2
      MT = MTNOW
      CALL RDWRSC(1,ISCR)
!
!     STORE AVERAGE ENERGIES
!
      DO I=1,NE
        Y(I) = AEBAR(I)
      END DO
      CALL RDWRSC(1,ISCR)
      MF5 = MF5 + 1
      MT5(MF5) = MF*1000 + MT
!
  100 RETURN
      END SUBROUTINE CKF15
!
!***********************************************************************
!
      SUBROUTINE CKNG
!
!           SUMS DISCRETE GAMMA ENERGIES WITH GAMMA CONTINUUM
!            AND CONSTRUCTS A UNION GRID OF GAMMAS AND NEUTRONS
!            CONSERVATION IS TESTED FOR GAMMAS AND GAMMAS PLUS NEUTRONS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: IABS, MOD
!
      CHARACTER(LEN=2) :: FLAG1,FLAG2
      CHARACTER(LEN=11), DIMENSION(7) :: ECHAR
      INTEGER(KIND=I4) :: MTT,MTNOW,MFT,MFN,MTN,INP,ICMPQ
      INTEGER(KIND=I4) :: NE,NENU,NEGC,NGT
      INTEGER(KIND=I4) :: I,J,K,N,N1213,ITOTAL,NQ,NMPTS,NZERO
      INTEGER(KIND=I4), DIMENSION(1) ::  NPLIN,ALLLIN
      REAL(KIND=R4) :: Q,ECMPQ,AWRX
      REAL(KIND=R4), DIMENSION(7) :: EFINAL
      REAL(KIND=R4), DIMENSION(2000) :: XGT,YGT,XGC,YGC,YGD
      REAL(KIND=R4), DIMENSION(1500) :: XN,YN,YNU
!
      CHARACTER(LEN=2), PARAMETER :: BLANK='  '
      CHARACTER(LEN=2), DIMENSION(2), PARAMETER :: FLAG=(/'* ','**'/)
      CHARACTER(LEN=11), PARAMETER :: XXX='     X     '
!
      ALLLIN(1) = 1
!
!     WRITE OUTPUT HEADING
!
      WRITE(NOUT,'(A/39X,A)') CHAR(12),                                 &       
     &           'COMPOSITE TEST OF FILES 5 12 13 15'
!
!     PROCESS ALL MT'S FOR WHICH COMPOSIT TEST POSSIBLE
!
      MTT = MT
      DO ITOTAL=1,MTTOT
         MTNOW = MTDCT(ITOTAL)
         NENU = 0
         NEGC = 0
         NGT = 0
!********GET Q VALUE FOR THIS REACTION
         DO NQ=1,NMT3
            IF(MTNOW.EQ.MT3(NQ))  THEN
               Q = QVAL(NQ)
               GO TO 20
            END IF
         END DO
         NQ = 0
         Q = 0.0
!
!        READ NEUTRON AND GAMMA CONTINUUM EBARS
!
   20    IF(MF5.GT.0)   THEN
            REWIND (UNIT=ISCR)
            DO J=1,MF5
               MFT = MT5(J)
               MTN = MOD(MFT,1000)
               IF(MTNOW.NE.MTN)   THEN
!*****************SKIP THIS BLOCK
                  CALL RDWRSC(2,ISCR)
                  CALL RDWRSC(2,ISCR)
               ELSE
                  CALL RDWRSC(2,ISCR)
                  MFN = MFT/1000
!*****************GAMMA CONTINUUM AVERAGES
                  IF(MFN.EQ.15)  THEN
                     NEGC = NP
                     DO I=1,NEGC
                        XGC(I) = X(I)
                        YGC(I) = Y(I)
                     END DO
                     GO TO 30
!*****************NEUTRON AVERAGES
                  ELSE
                     NENU = NP
                     DO I=1,NENU
                        XN(I) = X(I)
                        YNU(I) = Y(I)
                     END DO
                     CALL RDWRSC(2,ISCR)
                     DO I=1,NENU
                        YN(I) = Y(I)
                     END DO
                  END IF
               END IF
            END DO
         END IF
!
!        READ AND STORE NU TOTAL AND EBAR(GD)
!
   30    IF(NMTBAR.GT.0)  THEN
            REWIND (UNIT=LSCR)
            DO J=1,NMTBAR
               CALL RDWRSC(2,LSCR)
               IF(MTNOW.EQ.MTBAR(J))  THEN
                  NGT = NP
                  DO I=1,NGT
                     XGT(I) = X(I)
                     YGT(I) = Y(I)
                  END DO
                  GO TO 35
               ELSE
                  CALL RDWRSC(2,LSCR)
               END IF
            END DO
            GO TO 40
!
!           READ DISCRETE EBAR GAMMAS
!
   35       CALL RDWRSC(2,LSCR)
            NGT = NP
            DO I=1,NGT
               YGD(I) = Y(I)
            END DO
!
!           MERGE CONTINUUM WITH DISCRETE IF NECESSARY
!
            IF(NEGC.GT.0)   THEN
               DO K=1,NGT
                  YGD(K) = YGD(K) + YGC(K)
               END DO
            END IF
            GO TO 50
         END IF
!
!        NO NU TOTALS USE CONTINUUM NU
!
   40    IF(NEGC.GT.0)   THEN
            NGT = NEGC
            DO I=1,NGT
               XGT(I) = XGC(I)
               YGD(I) = YGC(I)
            END DO
            REWIND (UNIT=KSCR)
            DO N1213=1,MT1213
               CALL RDWRSC(2,KSCR)
               IF(MTNOW.EQ.IABS(MT))   THEN
                  DO I=1,NP
                     YGT(I) = Y(I)
                  END DO
                  GO TO 50
               END IF
            END DO
         END IF
         GO TO 90
!
   50    DO K=1,NGT
            IF(YGT(K).NE.0.)  THEN
               YGC(K) = YGD(K)/YGT(K)
            ELSE
               YGC(K) = 0.0
            END IF
         END DO
!
!        GET MERGER GRID
!
         NE = NENU + NGT
         IF(NENU.NE.0)   THEN
            IF(NGT.EQ.0)    THEN
               NE = NENU
               DO I=1,NE
                  X(I) = XN(I)
               END DO
            ELSE
               NE = 4000
               CALL GMERGE(XN,NENU,XGT,NGT,X,NE)
            END IF
         ELSE
            NE = NGT
            DO I=1,NE
               X(I) = XGT(I)
            END DO
         END IF
!
!        INTERPOLATE ALL EBARS TO A COMMON GRID
!
         IF(NENU.GT.0)  THEN
            INP = 1
            NPLIN(1) = NENU
            DO N=1,NE
               CALL TERPR(X(N),XN,YNU,NENU,NPLIN,ALLLIN,1,INP,Y(N))
            END DO
            DO I=1,NE
               YNU(I) = Y(I)
            END DO
            INP = 1
            DO N=1,NE
            CALL TERPR(X(N),XN,YN,NENU,NPLIN,ALLLIN,1,INP,Y(N))
            END DO
            DO I=1,NE
               XN(I) = Y(I)
            END DO
         END IF
         INP = 1
         NPLIN(1) = NGT
         DO N=1,NE
            CALL TERPR(X(N),XGT,YGD,NGT,NPLIN,ALLLIN,1,INP,Y(N))
         END DO
         DO I=1,NE
            YGT(I) = Y(I)
         END DO
         INP = 1
         DO N=1,NE
            CALL TERPR(X(N),XGT,YGC,NGT,NPLIN,ALLLIN,1,INP,Y(N))
         END DO
         DO I=1,NE
            YGC(I) = Y(I)
         END DO
!
!        DO TEST
!
         ICMPQ = 1
         IF(MTNOW.NE.4.AND.(MTNOW.LT.51.OR.MTNOW.GT.91))   ICMPQ = 2
         AWRX = AWR/(AWR+1.0)
         WRITE(NOUT,75)    MTNOW
   75    FORMAT(//3X,'SECTION ',I3//                                    &       
     &      13X,'E(IN)           TOTAL          AVERAGE          TOTAL',&       
     &      10X,'AVERAGE          ENERGY          TOTAL'/               &       
     &      20X,4(7X,'SECONDARY'),'       AVAILABLE        NEUTRON'/    &       
     &      19X,2(9X,'NEUTRON'),'          GAMMA           GAMMA',      &       
     &      25X,'AND GAMMA'/                                            &       
     &      18X,4(10X,'ENERGY'),26X,'ENERGY'/                           &       
     &      13X,7('(EV)',12X)/11X,7(10('-'),6X))
         DO NMPTS=1,NE
            FLAG1 = BLANK
            FLAG2 = BLANK
            EFINAL(1) = X(NMPTS)
            WRITE(ECHAR(1),'(1PE11.4)')  EFINAL(1)
            DO NZERO=2,7
               ECHAR(NZERO) = '           '
               EFINAL(NZERO) = 0.0
            END DO
!***********GET NEUTRON EBARS & NUEBARS
            IF(NENU.GT.0)   THEN
               EFINAL(2) = YNU(NMPTS)
               WRITE(ECHAR(2),'(1PE11.4)')  EFINAL(2)
               EFINAL(3) = XN(NMPTS)
               WRITE(ECHAR(3),'(1PE11.4)')   EFINAL(3)
            END IF
!***********GET GAMMA EBARS & NUEBARS
            EFINAL(4) = YGT(NMPTS)
            WRITE(ECHAR(4),'(1PE11.4)')   EFINAL(4)
            EFINAL(5) = YGC(NMPTS)
            WRITE(ECHAR(5),'(1PE11.4)')    EFINAL(5)
            EFINAL(7) = EFINAL(2) + EFINAL(4)
            WRITE(ECHAR(7),'(1PE11.4)')   EFINAL(7)
!***********GET AVAILABLE ENERGY
            IF(ICMPQ.EQ.2)  THEN
               ECMPQ = EFINAL(1)*AWRX + Q
            ELSE
               ECMPQ = EFINAL(1)*AWRX
            END IF
            IF(ECMPQ.LT.0.0)   THEN
               EFINAL(6) = 0.0
            ELSE
               EFINAL(6) = ECMPQ
            END IF
            WRITE(ECHAR(6),'(1PE11.4)')   EFINAL(6)
!***********SET FLAGS IF TESTS FAIL
            IF(EFINAL(4).GT.(EFINAL(6)*(1.+EPI4)))   FLAG1 = FLAG(ICMPQ)
            IF(EFINAL(7).GT.(EFINAL(6)*(1.+EPI4)).AND.NENU.NE.0)        &       
     &                   FLAG2 = FLAG(2)
!***********LOOK FOR GAMMMAS
            IF(NENU.EQ.0)  THEN
               ECHAR(2) = XXX
               ECHAR(3) = XXX
               ECHAR(7) = XXX
               FLAG2 = BLANK
               IF(MT.LE.0)   THEN
                  ECHAR(4) = XXX
                  ECHAR(6) = XXX
                  FLAG1 = BLANK
               END IF
            END IF
!
!           OUTPUT TEST RESULTS
!
            WRITE(NOUT,'(4X,4(5X,A11),1X,A2,2X,A11,2(5X,A11),2X,A2)')   &       
     &          (ECHAR(N),N=1,4),FLAG1,(ECHAR(N),N=5,7),FLAG2
         END DO
!
!        OUTPUT FLAG DEFINITION
!
         WRITE(NOUT,80)
   80    FORMAT(/4X,'X  INDICATES VALUE CANNOT BE CALCULATED'/          &       
     &           4X,'*  INDICATES ETOTAL GREATER THAN ECM'/             &       
     &           4X,'** INDICATES ETOTAL GREATER THAN ECM+Q')
   90 ENDDO
!
!     READY SCRATCH UNITS FOR ANOTHER MATERIAL
!
      REWIND (UNIT=ISCR)
      REWIND (UNIT=KSCR)
      REWIND (UNIT=LSCR)
      MT = MTT
!
      RETURN
      END SUBROUTINE CKNG
!
!***********************************************************************
!
      SUBROUTINE CAREN(IPATH)
!
!     SETS UP (IPATH=0) OR PERFORMS THE CAREN TESTS (IPATH GT 0)
!       TESTS CONTINUITY AT RESONANCE REGION BOUNDARIES
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPATH
!
      INTEGER(KIND=I4), SAVE :: MBOUND
      INTEGER(KIND=I4) :: IPATH1,IP,INTT,IND,NCURR
      INTEGER(KIND=I4) :: I,II,J,N
      INTEGER(KIND=I4), PARAMETER :: NGRUPS=203
      REAL(KIND=R4), SAVE, DIMENSION(NGRUPS) :: EGR
      INTEGER(KIND=I4), PARAMETER :: NREACS= 4
      REAL(KIND=R4), SAVE, DIMENSION(NREACS) :: RESINT=0.,SIGO
      REAL(KIND=R4), SAVE, DIMENSION(NGRUPS,NREACS) :: SIGP=0.
      REAL(KIND=R4), SAVE :: EPIM,EPIP,ECURR,XSEC
      REAL(KIND=R4) :: XLP,XHP,ANS
!
      INTEGER(KIND=I4), PARAMETER :: NCARN=4
      INTEGER(KIND=I4), DIMENSION(NCARN), PARAMETER ::                   &      
     &                ICARN = (/1,2,18,102/)
!
!     PREPARE FOR CAREN TESTS
!
      IF(IPATH.EQ.0)  THEN
!
!        SET UP ENERGY GRID FOR BOUNDARY CONTINUITY TESTS AND THERMAL
!          VALUE CALCULATION
!
         MBOUND = 1
         EGR(1) = THERMA
         IF(LRP.NE.0)   THEN
            EPIM = 1. - EPI6
            EPIP = 1. + EPI6
            MBOUND = 3
            EGR(2) = E1*EPIM
            EGR(3) = E1*EPIP
   10       ECURR = 10.E+6
            NCURR = 0
            DO I=1,NBOUND
               IF(RBOUND(I).NE.0.) THEN
                  IF(ECURR.GT.RBOUND(I)) THEN
                     ECURR = RBOUND(I)
                     NCURR = I
                  END IF
               END IF
            END DO
            IF(NCURR.NE.0) THEN
               MBOUND = MBOUND + 2
               EGR(MBOUND) = ECURR*EPIP
               EGR(MBOUND-1) = ECURR*EPIM
               RBOUND(NCURR) = 0.
               GO TO 10
            END IF
!
!           CALCULATE RESONANCE INTEGRAL FROM FILE 2 DATA
!
            CALL INDER(SIGO,RESINT,NREACS)
            IF(LRP.LE.1)   THEN
!
!              CALCULATE SIGMA ON THE BREAKPOINT GRID
!
               DO N=1,MBOUND
                  IF(N.EQ.MBOUND-3.AND.IUNR.EQ.1) THEN
                     IPATH1 = 2
                  ELSE
                     IPATH1 = 1
                  END IF
                  CALL SIGMA(EGR(N),SIGO,IPATH1)
                  DO J=1,NREACS
                     SIGP(N,J) = SIGO(J)
                  END DO
               END DO
            END IF
            GO TO 100
         ELSE
            DO J=1,NREACS
               SIGP(1,J) = 0.0
               RESINT(J) = 0.0
            END DO
         END IF
      END IF
!
!     PERFORM CAREN TESTS
!
      DO I=1,NCARN
         IF(ICARN(I).EQ.MT)  GO TO 20
      END DO
      GO TO 100
!
!     ADD SMOOTH BACKGROUND
!
   20 IND = I
      DO N=1,MBOUND
         IF(N.LT.3)  THEN
            IP = 1
            INTT = 1
         END IF
         CALL INTER(IP,INTT,EGR(N),XSEC)
         SIGP(N,IND) = SIGP(N,IND) + XSEC
      END DO
!
!     CALCULATE SMOOTH CONTRIBUTION TO RESONANCE INTEGRAL
!
      IF(IND.GE.3)   THEN
         XLP = 0.5
         XHP = ENMAX
         CALL GRATE(XLP,XHP,ANS)
         RESINT(IND) = RESINT(IND) + ANS
      END IF
!
!     OUTPUT THERMAL CROSS SECTIONS AND RESONANCE INTEGRALS
!
      WRITE(NOUT,'(/5X,A/14X,A,35X,A,1PE12.5)')                         &       
     &      'THERMAL CROSS SECTIONS AND RESONANCE INTEGRALS',           &       
     &      'E = THERMAL','SIGMA   = ',SIGP(1,IND)
      IF(IND.GT.2)   THEN
         WRITE(NOUT,'(14X,A,19X,1PE12.5)')                              &       
     &       'RESONANCE INTEGRAL, 0.5 EV CUTOFF IS ',RESINT(IND)
      END IF
!
!     OUTPUT BOUNDARY VALUES
!
      IF(MBOUND.GT.1)   THEN
         WRITE(NOUT,'(/5X,A)') 'RESONANCE REGION BOUNDARY TESTS'
         DO II=2,MBOUND,2
            IF(II.EQ.2.AND.EGR(II).LT.ENMIN)  THEN
               WRITE(NOUT,'(/14X,A,1PE12.5,30X,A,1PE12.5)')             &       
     &             'E = ',EGR(II)/EPIM,'SIGMA+  = ',SIGP(II+1,IND)
            ELSE
               WRITE(NOUT,'(/A,1PE12.5,30X,A,1PE12.5/60X,A,1PE12.5)')   &       
     &            '              E = ',EGR(II)/EPIM,'SIGMA-  = ',       &       
     &            SIGP(II,IND),'SIGMA+  = ',SIGP(II+1,IND)
            END IF
         END DO
      END IF
!
  100 RETURN
      END SUBROUTINE CAREN
!
!***********************************************************************
!
      SUBROUTINE INDER(SIG,RESINT,NREACS)
!
!     INFINITE DILUTE RESONANCE CROSS SECTION CALCULATION
!
      IMPLICIT NONE
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT
!
      INTEGER(KIND=I4) :: LRU,LRF
      INTEGER(KIND=I4) :: NE,NG,NER,LL,NLS,NRS,NSS,NJS
      INTEGER(KIND=I4) :: I,K,K1,M1,N1,N,NN,NNN,NNNN
      REAL(KIND=R4) :: EL,EU,EAB
      REAL(KIND=R4) :: EMIN,EMID,EMAX
      REAL(KIND=R4) :: ABN,KFAC,RFAC,ARAT,AWRI,AJ,SPI,AC,DEN,U,DU
      REAL(KIND=R4) :: E0,ER,ERK,G,GN,GG,GF,GX,GFISS,GNK
      REAL(KIND=R4) :: RHO,SF,SHIFT,PF,RHOG,SHIFTK,PFK
      REAL(KIND=R4) :: GAM,GAM1,GAM2,GAM3,GAM4,GAM5,GAM6
      REAL(KIND=R4) :: GSQR,GGG,GAR
      REAL(KIND=R4) :: SR0,SR,Y1,Y2,SRY1,SRY2
      REAL(KIND=R4) :: GY1,GY2,GY3,GY4
      REAL(KIND=R4) :: YG1,YG2,YG3,YG4
      REAL(KIND=R4) :: ARCTAN,ARGLOG,ARCCTN
      REAL(KIND=R4) :: GU,CEE,PENDEP,GAMRCY,CEESQR,CEE1
!
      INTEGER(KIND=I4) :: NREACS
      REAL(KIND=R4), DIMENSION(NREACS) :: RESINT,SIG
      INTEGER(KIND=I4), PARAMETER :: NEMAX=151
      REAL(KIND=R4), DIMENSION(NEMAX) :: EGR
!
!     PROCESS ALL ISOTOPES
!
      RESINT(3) = 0.
      RESINT(4) = 0.
      DO 85 NN=1,NIS
      ABN = ABNS(NN)
!
!     READ RESONANCE PARAMETERS FROM SCRATCH FILE
!
      CALL RWRES(2)
!
!     PROCESS ALL RESONANCE REGIONS
!
      NER = NERS(NN)
      DO 75 NNN=1,NER
      IPT = IRES(NNN)
      LRU = RESPAR(IPT)
      LRF = RESPAR(IPT+1)
      IF(LRF.GT.3.AND.LRF.NE.6)  GO TO 75
      EMID = EMIDLE(NNN)
      EMAX = AMIN1(EMID,E2)
      IF(NNN.EQ.1)  THEN
         EMIN = AMAX1(E1,0.5)
      ELSE
         EMIN = AMAX1(EMIDLE(NNN-1),0.5)
      END IF
      IF(EMIN.GE.EMAX)   GO TO 75
      DU = 0.1
!
!     A RESOLVED REGION
!
      IF(LRU.EQ.1) THEN
         U = 0.
         DO N=1,NEMAX
            NE = N
            EGR(N) = EMAX*EXP(-U)
            IF(EGR(N).LE.EMIN) THEN
               EGR(NE) = EMIN
               GO TO 10
            END IF
            U = U + DU
         END DO
   10    NG = NE - 1
!
!        RETRIEVE NUCLIDE INFORMATION
!
         SPI = RESPAR(IPT+2)
         NLS = RESPAR(IPT+4)
         DEN = 2.*(2.*SPI+1.)
!
!        RETRIEVE COMPETING CHANNEL DESCRIPTIONS FOR LRF=6
!
         IF(LRF.EQ.6)  THEN
            IPT = RESPAR(IPT+17)
            AWRI = RESPAR(IPT)
!
!        LRF = 1,2,3
!
         ELSE
            AWRI = RESPAR(IPT+5)
            IF(LRF.EQ.3) THEN
               IPT = IPT + 6
            ELSE
               AC = RESPAR(IPT+6)
               IPT = IPT + 7
            END IF
            NLS = MIN0(3,NLS)
         END IF
         ARAT = AWRI/(AWRI+1.0)
         KFAC = CROC*ARAT
!
!        LOOP OVER ALL INCIDENT ANGULAR MOMENTUM
!
         DO 65 K1=1,NLS
         IF(LRF.EQ.6) THEN
            LL = RESPAR(IPT+1)
            NSS = RESPAR(IPT+2)
            IPT = IPT + 3
         ELSE
            IF(LRF.EQ.3)  THEN
               AC = RESPAR(IPT)
               LL = RESPAR(IPT+3)
               NRS = RESPAR(IPT+5)
               IPT = IPT + 6
            ELSE
               LL = RESPAR(IPT+1)
               NRS = RESPAR(IPT+3)
               IPT = IPT + 4
            END IF
            NSS = 1
         END IF
         DO 60 N1=1,NSS
         IF(LRF.EQ.6)  THEN
            NJS = RESPAR(IPT+1)
            IPT = IPT + 2
         ELSE
            NJS = 1
         END IF
         DO 50 M1=1,NJS
         IF(LRF.EQ.6) THEN
            AJ = RESPAR(IPT)
            AC = RESPAR(IPT+1)
            NRS = RESPAR(IPT+4)
            IPT = IPT + 5
         END IF
!
!        LOOP OVER ALL RESONANCES
!
         RFAC = KFAC*AC
         DO K=1,NRS
            ER = RESPAR(IPT)
            IF(LRF.EQ.6)  THEN
               GN = RESPAR(IPT+1)
               GG = RESPAR(IPT+2)
               GF = RESPAR(IPT+3)
               IPT = IPT + 12
            ELSE
               AJ = RESPAR(IPT+1)
               GN = RESPAR(IPT+3)
               GG = RESPAR(IPT+4)
               GF = RESPAR(IPT+5)
               GX = RESPAR(IPT+8)
               IF(LRF.EQ.3) GF = ABS(GF) + ABS(GX)
               IPT = IPT + 9
            END IF
            E0 = ABS(ER)
            RHO = RFAC*SQRT(E0)
            SF = (2.*AJ+1.)/DEN
            CALL FACTS(LL,RHO,SHIFT,PF)
!
!           LOOP OVER ALL GROUPS
!
            DO I=1,NG
               EL = EGR(I)
               EU = EGR(I+1)
               EAB = (EL+EU)/2.
               RHOG = RFAC*SQRT(EAB)
               CALL FACTS(LL,RHOG,SHIFTK,PFK)
               GNK = GN*PFK/PF
               ERK = ER + (SHIFT-SHIFTK)*GN/(2.*PF)
               E0 = ABS(ERK)
!**************ASSUME THAT INELASTIC AND CP CHANNELS HAVE A NEGLIGIBLE
!**************  EFFECT ON FISSION AND CAPTURE RESONANCE INTEGRALS
               G = GG + GF + GNK + GX
               GNK = GNK/SQRT(EAB)
               IF(GG.GT.0.) THEN
                  GFISS = GF/GG
               ELSE
                  GFISS = 0.
               END IF
               GSQR = SQRT(G/2.)
               GAM = 2.*ERK/G
               GAM1 = 1. + GAM*GAM
               GAM2 = SQRT(GAM1)
               GGG = GAM2 + GAM
               IF(GGG.GE.0.0) THEN
                  GAM3 = SQRT(GAM2+GAM)
               ELSE
                  GAM3 = 0.0
               END IF
               IF(GAM2.GT.GAM) THEN
                  GAM4 = GAM2 - GAM
               ELSE
                  GAM4 = .5/ABS(GAM)
               END IF
               GAM4 = SQRT(GAM4)
               GAM5 = 0.5*(GAM3+GAM*GAM4)
               GAM6 = GAM*GAM3-GAM4
               GAR = ARAT*G
               SR0 = 2.6E+06*SF*GG*GNK/(GAR*GAR)
               Y1 = 2.*EL/G
               Y2 = 2.*EU/G
               SRY1 = SQRT(Y1)
               SRY2 = SQRT(Y2)
               GY1 = (1.414214*SRY1-GAM3)/GAM4
               GY2 = (1.414214*SRY1+GAM3)/GAM4
               GY3 = (1.414214*SRY2-GAM3)/GAM4
               GY4 = (1.414214*SRY2+GAM3)/GAM4
               ARCTAN = ATAN(GY1) + ATAN(GY2) - ATAN(GY3) - ATAN(GY4)
               YG1 = Y1 + GAM2
               YG2 = Y2 + GAM2
               YG3 = 1.414214*SRY1*GAM3
               YG4 = 1.414214*SRY2*GAM3
               SR = 0.0
               IF((YG1-YG3.NE.0.0).AND.(YG2-YG4.NE.0))   THEN
                  ARGLOG=ABS((YG1+YG3)*(YG2-YG4)/((YG1-YG3)*(YG2+YG4)))
                  IF(LL.EQ.0)   THEN
!
!                    S WAVE
!
                    GU = GAM1
                    SR = (SR0/(GSQR*GU))*(2.0*(1.0/SRY2-1.0/SRY1)       &       
     &                    +(1.0/(1.414214*GAM2))*                       &       
     &                    (GAM5*ALOG(ARGLOG)+GAM6*ARCTAN))
                  ELSE
!
!                    APPROXIMATION FOR P AND D WAVES
!
                     CEE = 2.*RHOG*RHOG/(Y1+Y2)
                     PENDEP = (1.+RHOG*RHOG)/EAB
                     GAMRCY = (1./CEE) + GAM
                     CEESQR = SQRT(CEE)
                     CEE1 = GAMRCY  + GAM + CEE*GAM1
                     ARCCTN = ATAN( CEESQR*SRY1) - ATAN( CEESQR*SRY2)
                     SR = ((SR0*PENDEP*GSQR)/CEE1)*(2.*CEESQR*ARCCTN    &       
     &                  +(1.0/(1.414214*GAM2))*((GAM3+GAMRCY*GAM4)      &       
     &                  *ALOG(ARGLOG) + ((GAMRCY*GAM3-GAM4)*ARCTAN)))
                   END IF
               END IF
               RESINT(3) = RESINT(3) + ABN*SR*GFISS
               RESINT(4) = RESINT(4) + ABN*SR
            END DO
         END DO
         IF(LRF.EQ.6)  IPT = RESPAR(IPT)
   50    CONTINUE
   60    CONTINUE
   65    CONTINUE
!
!     UNRESOLVED REGION
!
      ELSE IF(LRU.EQ.2) THEN
         U = 0.
         DO N=1,NEMAX
         NE = N
         EGR(N) = EMAX*EXP(-U)
         IF(EGR(N).LE.EMIN) THEN
            EGR(NE) = EMIN
            GO TO 70
         END IF
         U = U + DU
         END DO
   70    NG = NE - 1
         IF(NG.GE.1)   THEN
!
!           GENERATE CROSS SECTIONS AT MID LETHARGY OF GROUP
!
            DO N=1,NG
               EAB = SQRT(EGR(N)*EGR(N+1))
               IF(EAB.GE.E1.AND.EAB.LE.E2)   THEN
                  NNNN = NN
                  CALL SIG1(EAB,SIG,NNNN,1)
                  DU = ALOG(EGR(N)/EGR(N+1))
                  RESINT(3) = RESINT(3) + SIG(3)*DU*ABN
                  RESINT(4) = RESINT(4) + SIG(4)*DU*ABN
               END IF
            END DO
         END IF
      END IF
   75 CONTINUE
   85 CONTINUE
!
      REWIND(ISCR)
!
      RETURN
      END SUBROUTINE INDER
!
!***********************************************************************
!
      SUBROUTINE SIGMA(E,SIG,IPATH)
!
!     CONTROLS THE CALCULATION OF THE RESONANCE CROSS SECTIONS, SIG,
!       AT ENERGY, E FOR A MIXTURE OF ISOTOPES
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPATH
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIG
!
      INTEGER(KIND=I4) :: J,N
      REAL(KIND=R4) :: ABN
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      SIG = 0.
!
!     IS E IN THE RESONANCE REGION?
!
      IF(E.GE.E1.AND.E.LE.E2)   THEN
!
!        LOOP OVER ALL ISOTOPES
!
         DO N=1,NIS
            ABN = ABNS(N)
!
!           READ RESONANCE PARAMETERS FROM SCRATCH FILE
!
            CALL RWRES(2)
!
!           GET SIGMA AT E
!
            CALL SIG1(E,SIGP,N,IPATH)
!
!           SUM ISOTOPE CONTRIBUTIONS
!
            DO J=1,4
               SIG(J) = SIG(J) + ABN*SIGP(J)
            END DO
         END DO
         REWIND ISCR
      END IF
!
      RETURN
      END SUBROUTINE SIGMA
!
!***********************************************************************
!
      SUBROUTINE SIG1(EC,SIG,NNIS,IPATH)
!
!     CALCULATES E RESONANCE CROSS SECTIONS, SIG, AT ENERGY, E, FOR A
!       SINGLE ISOTOPE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NNIS,IPATH
      REAL(KIND=R4) :: EC
      REAL(KIND=R4), DIMENSION(4) :: SIG
!
      INTEGER(KIND=I4) :: LRU,LRF,LSSF,LRUP,NER
      INTEGER(KIND=I4) :: INOWP,NLOOP,IBEG
      INTEGER(KIND=I4) :: I,KK
      REAL(KIND=R4) :: E,E0,DENOM,DE,ANS
      REAL(KIND=R4), DIMENSION(4) :: SIGM,SIGO
!
      SIG = 0.
!
!     LOOK FOR ENERGY REGION FOR THIS ISOTOPE
!
      NER = NERS(NNIS)
      DO I=1,NER
         IF(EC.LE.EMIDLE(I))   GO TO 10
      END DO
      GO TO 100
!
!     SET REPRESENTATION
!
   10 IPT = IRES(I)
      LRU = RESPAR(IPT)
      LRF = RESPAR(IPT+1)
      IF(LRU.EQ.2)  THEN
         LSSF = RESPAR(IPT+4)
      ELSE
         LSSF = 0
      END IF
      INOWP = IPT
!
!     INITIALIZE FOR CALCULATION AT A SINGLE ENERGY
!
      NLOOP = 1
      E = EC
      IF(IPATH.EQ.2)   THEN
!
!     IF REQUIRED AT THE UPPER BOUNDARY OF ALL THE RESOLVED REGIONS,
!       CALCULATE AN AVERAGE JUST BELOW THE BOUNDARY
!
         IF(I.GE.NER-1)   THEN
            E0 = 0.95*EMIDLE(I)
            IF(E.GE.E0)   THEN
               LRUP = RESPAR(IRES(NER))
               IF(I.EQ.NER)  THEN
                  IF(LRUP.NE.1)  GO TO 20
               ELSE
                  IF(LRUP.EQ.1)  GO TO 20
               END IF
               NLOOP = 1000
               E = E0
               DENOM = 0.05*EMIDLE(I)
               DE = DENOM/NLOOP
            END IF
         END IF
      END IF
!
!     CALCULATE CROSS SECTIONS OR AVERAGES
!
   20 DO KK=1,NLOOP
         IPT = INOWP
!
!        RESOLVED RESONANCE REGION
!
         IF(LRU.EQ.1)   THEN
!
!           GET ANY ENERGY DEPENDENT SCATTERING LENGTH
!
            IF(NPED(I).GT.0)   THEN
               IBEG = 1
               CALL TERPR(E,EP(1,I),APED(1,I),NPED(I),NBTED(1,I),       &       
     &             JNTED(1,I),NRED(I),IBEG,RESPAR(IPT+3))
            END IF
!
!           SINGLE LEVEL BREIT-WIGNER
!
            IF(LRF.EQ.1) THEN
                CALL CSSLBW(E,SIGO,IPT)
!
!           MULTI-LEVEL BREIT-WIGNER
!
            ELSE IF(LRF.EQ.2) THEN
               CALL CSMLBW(E,SIGO,IPT)
!
!           R-MATRIX (REICH-MOORE)
!
            ELSE IF(LRF.EQ.3) THEN
               CALL CSRMAT(E,SIGO,IPT)
!
!           ADLER-ADLER
!
            ELSE IF(LRF.EQ.4) THEN
               CALL CSAA(E,SIGO,IPT)
!
!           HYBRID R-FUNCTION
!
            ELSE IF(LRF.EQ.6) THEN
               CALL CSHYBR(E,SIGO,IPT)
!
!           UNRECOGNIZED OPTION
!
            ELSE
               GO TO 100
            END IF
!
!        CALCULATION IN THE UNRESOLVED REGION
!
         ELSE IF(LRU.EQ.2) THEN
            IF(LSSF.EQ.1)  GO TO 100
!
!           REPRESENTATION ONE
!
            IF(LRF.EQ.1) THEN
               CALL CSUNR1(E,SIGO,IPT)
!
!           REPRESENTATION TWO
!
            ELSE IF(LRF.EQ.2) THEN
               CALL CSUNR2(E,SIGO,IPT)
            END IF
!
!        POTENTIAL ONLY
!
         ELSE IF(LRU.EQ.0)  THEN
            SIGO = 0.
         END IF
!
         IF(NLOOP.EQ.1)  THEN
            SIG = SIGO
            GO TO 100
         END IF
         IF(KK.NE.1)   THEN
            DO I=1,4
               CALL ECSI(E-DE,SIGM(I),E,SIGO(I),E-DE,E,2,ANS)
               SIG(I) = SIG(I) + ANS/DENOM
            END DO
            E = E + DE
         END IF
         SIGM = SIGO
      END DO
!
  100 RETURN
      END SUBROUTINE SIG1
!
!***********************************************************************
!
      SUBROUTINE CSSLBW(E,SIGP,INOW)
!
! **********************************************************************
! *   CALCULATES SINGLE LEVEL BREIT-WIGNER CROSS SECTIONS AT ENERGY E  *
! *   FOR ONE SECTION   (ONE ISOTOPE-ONE ENERGY RANGE)                 *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT, COS, SIN
!
      INTEGER(KIND=I4) :: NLS
      INTEGER(KIND=I4) :: LL,LRX,NRS
      INTEGER(KIND=I4) :: I,L
      REAL(KIND=R4) :: AWRI,AP,A,SPI,QX
      REAL(KIND=R4) :: K,KFAC,ARAT,PIFAC,RHO,RHOP,PHI
      REAL(KIND=R4) :: ER,AJ,GJ,GN,GG,GF,GC,GNE,GX,GTT
      REAL(KIND=R4) :: SE,PE,SER,PER,SEC,PEC,RHOC,COS2P,SINSQ,SPOT
      REAL(KIND=R4) :: ERP,EDELT,COMFAC,DEN
!
!     ZERO OUT CROSS SECTION CONTRIBUTION FROM THIS SECTION (SIGP)
!
      SIGP = 0.0
!
!     RETRIEVE NUCLIDE INFORMATION
!
      AWRI = RESPAR(INOW+5)
      AP = RESPAR(INOW+3)
      A = RESPAR(INOW+6)
      SPI = RESPAR(INOW+2)
      NLS = RESPAR(INOW+4)
      INOW = INOW + 7
!
!     CALCULATE WAVE NUMBER(K),RHO AND RHOCAP AT ENERGY (E)
!
      ARAT = AWRI/(AWRI+1.0)
      KFAC = CROC*ARAT
      K = KFAC*SQRT(ABS(E))
      PIFAC = PI/(K*K)
      RHO = K*A
      RHOP = K*AP
!
!     LOOP OVER L STATES
!
      DO L=1,NLS
         QX = RESPAR(INOW)
         LL = RESPAR(INOW+1)
         LRX= RESPAR(INOW+2)
         NRS = RESPAR(INOW+3)
         INOW = INOW + 4
!********CALCULATE SHIFT AND PENETRATION FACTORS AT CROSS SECTION ENERGY
         CALL FACTS(LL,RHO,SE,PE)
         PEC = 0.0
         IF(LRX.NE.0)  THEN
            RHOC = KFAC*SQRT(ABS(E+QX))*A
            CALL FACTS(LL,RHOC,SEC,PEC)
         END IF
!
!        CALCULATE CONSTANTS FOR POTENTIAL SCATTERING
!
         CALL FACPHI(LL,RHOP,PHI)
         COS2P = COS(2.*PHI)
         SINSQ = (SIN(PHI))**2
         SPOT = 4.*(2.*FLOAT(LL)+1.)*PIFAC*SINSQ
!
!        LOOP OVER ALL RESONANCES
!
         DO I=1,NRS
            ER = RESPAR(INOW)
            AJ = RESPAR(INOW+1)
            GN = RESPAR(INOW+3)
            GG = RESPAR(INOW+4)
            GF = RESPAR(INOW+5)
            SER = RESPAR(INOW+6)
            PER = RESPAR(INOW+7)
            GC = RESPAR(INOW+8)
            INOW = INOW + 9
            ERP = ER + GN*(SER-SE)*(0.5/PER)
            EDELT = E - ERP
            GNE = GN*PE/PER
            GX = GG + GF
            IF(LRX.NE.0) THEN
               GTT = GNE + GX + GC*PEC
            ELSE
               GTT = GNE + GX
            END IF
            GJ = 0.5*(2.0*AJ+1.0)/(2.0*SPI+1.0)
!
!           CROSS SECTION CALCULATIONS
!
            DEN = EDELT**2 + 0.25*GTT*GTT
            COMFAC = PIFAC*GJ*GNE/DEN
!***********ELASTIC
            SIGP(2) = SIGP(2) +                                         &       
     &            COMFAC*(GNE*COS2P-2.0*GX*SINSQ+2.0*EDELT*SIN(2.0*PHI))
!***********FISSION
            IF(GF.GT.0.0) THEN
               SIGP(3) = SIGP(3) + COMFAC*GF
            END IF
!***********CAPTURE
            SIGP(4) = SIGP(4) + COMFAC*GG
         END DO
!
!        ADD POTENTIAL SCATTERING
!
         SIGP(2) = SIGP(2) + SPOT
      END DO
!
!     SUM TO GET TOTAL
!
      SIGP(1) = SIGP(2) + SIGP(3) + SIGP(4)
!
      RETURN
      END SUBROUTINE CSSLBW
!
!***********************************************************************
!
      SUBROUTINE CSMLBW(E,SIGP,INOW)
!
! **********************************************************************
! *   CALCULATES MULTI-LEVEL BREIT-WIGNER CROSS SECTIONS AT ENERGY E   *
! *   FOR ONE SECTION   (ONE ISOTOPE-ONE ENERGY RANGE)                 *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT, COS, SIN
!
      INTEGER(KIND=I4) :: NLS
      INTEGER(KIND=I4) :: LL,LRX,NRS
      INTEGER(KIND=I4) :: J,NUMJ
      INTEGER(KIND=I4) :: I,L
      REAL(KIND=R4) :: AWRI,AP,A,SPI,QX
      REAL(KIND=R4) :: K,KFAC,ARAT,PIFAC,RHO,RHOP,PHI
      REAL(KIND=R4) :: ER,AJ,GN,GG,GF,GC,GNE,GX,GTT
      REAL(KIND=R4) :: SE,PE,SER,PER,SEC,PEC,RHOC,COS2P,SIN2P
      REAL(KIND=R4) :: SSUM,FL,AJMIN,AJMAX
      REAL(KIND=R4) :: ERP,EDELT,COMFAC,DEN,DIFF,XF
      REAL(KIND=R4), DIMENSION(20) :: GJ
      REAL(KIND=R4), DIMENSION(20,2) :: SIGJ
!
!     ZERO OUT CROSS SECTION CONTRIBUTION FROM THIS SECTION (SIGP)
!
      SIGP = 0.0
!
!     RETRIEVE NUCLIDE INFORMATION
!
      AWRI = RESPAR(INOW+5)
      AP = RESPAR(INOW+3)
      A = RESPAR(INOW+6)
      SPI = RESPAR(INOW+2)
      NLS = RESPAR(INOW+4)
      INOW = INOW + 7
!
!     CALCULATE WAVE NUMBER(K),RHO AND RHOCAP AT ENERGY (E)
!
      ARAT = AWRI/(AWRI+1.0)
      KFAC = CROC*ARAT
      K = KFAC*SQRT(ABS(E))
      PIFAC = PI/(K*K)
      RHO = K*A
      RHOP = K*AP
      DEN = 2.*(2.0*SPI+1.)
!
!     LOOP OVER L STATES
!
      DO L=1,NLS
         QX = RESPAR(INOW)
         LL = RESPAR(INOW+1)
         LRX = RESPAR(INOW+2)
         NRS = RESPAR(INOW+3)
         INOW = INOW + 4
!********CALCULATE SHIFT AND PENETRATION FACTORS AT CROSS SECTION ENERGY
         CALL FACTS(LL,RHO,SE,PE)
         PEC = 0.0
         IF(LRX.NE.0)   THEN
            RHOC=KFAC*SQRT(ABS(E+QX))*A
            CALL FACTS(LL,RHOC,SEC,PEC)
         END IF
!
!        CALCULATE CONSTANTS FOR POTENTIAL SCATTERING
!
         CALL FACPHI(LL,RHOP,PHI)
         COS2P = 1.0 - COS(2.*PHI)
         SIN2P = SIN(2.0*PHI)
!
!        STATISTICAL FACTORS RECONSTRUCTED
!
         SSUM = 0.
         FL = FLOAT(LL)
         AJMIN = ABS(ABS(SPI-FL)-0.5)
         AJMAX = SPI +  FL + 0.5
         NUMJ = IFIX(AJMAX-AJMIN) + 1
         AJ = AJMIN
         DO I=1,NUMJ
            GJ(I) = (2.*AJ+1.)/DEN
            AJ = AJ + 1.
            SSUM = SSUM + GJ(I)
         END DO
         DIFF = 2.*FL + 1. - SSUM
!
!        INITIALIZE CROSS SECTION PARTIALS
!
         SIGJ = 0.0
!
!        LOOP OVER ALL RESONANCES
!
         DO I=1,NRS
            ER = RESPAR(INOW)
            AJ = RESPAR(INOW+1)
            J = IFIX(AJ-AJMIN) + 1
            GN = RESPAR(INOW+3)
            GG = RESPAR(INOW+4)
            GF = RESPAR(INOW+5)
            SER = RESPAR(INOW+6)
            PER = RESPAR(INOW+7)
            GC = RESPAR(INOW+8)
            INOW = INOW + 9
            ERP = ER + GN*(SER-SE)*(0.5/PER)
            EDELT = E - ERP
            GNE = GN*PE/PER
            GX = GG + GF
            IF(LRX.NE.0)   THEN
               GTT = GNE + GX + GC*PEC
            ELSE
               GTT = GNE + GX
            END IF
            XF = 2.0*EDELT/GTT
!
!           CROSS SECTION CALCULATIONS
!
            COMFAC = 2.0*GNE/GTT/(1.0+XF*XF)
!***********ELASTIC COMPONENTS
            SIGJ(J,1) = SIGJ(J,1) + COMFAC
            SIGJ(J,2) = SIGJ(J,2) + COMFAC*XF
!***********FISSION
            IF(GF.GT.0.0) THEN
               SIGP(3) = SIGP(3) + COMFAC*GJ(J)*GF/GTT
            END IF
!***********CAPTURE
            SIGP(4)=SIGP(4)+COMFAC*GJ(J)*GG/GTT
         END DO
!
!        ADD POTENTIAL SCATTERING
!
         DO J=1,NUMJ
            SIGP(2) = SIGP(2) + GJ(J)*((COS2P-SIGJ(J,1))**2             &       
     &                 + (SIN2P+SIGJ(J,2))**2)
         END DO
         SIGP(2) = SIGP(2) + 2.*DIFF*COS2P
      END DO
!
!     FINAL CROSS SECTION VALUES
!
      SIGP(2) = SIGP(2)*PIFAC
      SIGP(3) = SIGP(3)*2.0*PIFAC
      SIGP(4) = SIGP(4)*2.0*PIFAC
      SIGP(1) = SIGP(2) + SIGP(3) + SIGP(4)
!
      RETURN
      END SUBROUTINE CSMLBW
!
!***********************************************************************
!
      SUBROUTINE CSRMAT(E,SIGP,INOW)
!
! **********************************************************************
! *   CALCULATES R-MATRIX(REICH-MOORE) CROSS SECTIONS AT ENERGY E      *
! *   FOR ONE SECTION  (ONE ISOTOPE-ONE ENERGY RANGE)                  *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT
      REAL(KIND=8), INTRINSIC :: DABS, DSQRT, DCOS, DSIN
!
      INTEGER(KIND=I4) :: NLS,NRS,LL
      INTEGER(KIND=I4) :: INOWB,JJL,NUMJ,IFIS
      INTEGER(KIND=I4) :: KPSTV,KNGTV,ISKIP,KKKKKK
      INTEGER(KIND=I4) :: I,J,JJ,L,KCHANL
      REAL(KIND=R4) :: SPI,GJD
      REAL(KIND=R4) :: A,AP,RHO,RHOP,SE,PE,PHI
      REAL(KIND=R4) :: ER,GN,GG,GF,GC,PER,A1,A2,A3
      REAL(KIND=R4) :: FL,AJMIN,AJMAX,AJC,AJ,AJPM
      REAL(KIND=8) :: DEN,DIFF,DE2,P1,P2,GJ,GG4,T1,T2,T3,T4
      REAL(KIND=8) :: TERMT,TERMN,TERMF,TERMG,SIGNNI,SIGNGI,SIGNFI,     &       
     &                 SIGNTI
      REAL(KIND=8) :: ED,ARAT,AWRI,KFAC,K,PIFAC,PHID,U11R,U11I
      REAL(KIND=8), DIMENSION(3,3) ::  R,S,RI,SI
!
!     ZERO OUT CROSS SECTION CONTRIBUTION FROM THIS SECTION (SIGP)
!
      SIGP = 0.0
!
!     RETRIEVE NUCLIDE INFORMATION
!
      AWRI = RESPAR(INOW+5)
      SPI = RESPAR(INOW+2)
      GJD = 2.0*(2.0*SPI+1.0)
      NLS = RESPAR(INOW+4)
      INOW = INOW + 6
!
!     CALCULATE WAVE NUMBER(K), RHO AND RHOCAP AT ENERGY (E)
!
      ED = E
      ARAT = AWRI/(AWRI+1.0)
      KFAC = CROC*ARAT
      K = KFAC*DSQRT(DABS(ED))
      PIFAC = PI/(K*K)
!
!     INITIALIZE PARTIAL CROSS SECTIONS
!
      SIGNNI = 0.0
      SIGNGI = 0.0
      SIGNFI = 0.0
      SIGNTI = 0.0
!
!     LOOP OVER L STATES
!
      DO L=1,NLS
         AP = RESPAR(INOW)
         A = RESPAR(INOW+1)
         RHO = K*A
         RHOP = K*AP
         LL = RESPAR(INOW+3)
         NRS = RESPAR(INOW+5)
         INOW = INOW + 6
         INOWB = INOW
!********CALCULATE SHIFT AND PENETRATION FACTORS AT CROSS SECTION ENERGY
         CALL FACTS(LL,RHO,SE,PE)
!
!        CONSTANTS FOR POTENTIAL SCATTERING
!
         CALL FACPHI(LL,RHOP,PHI)
         PHID = PHI
         P1 = DCOS(2.0*PHID)
         P2 = DSIN(2.0*PHID)
!
!        LOOP OVER POSSIBLE J VALUES
!
         FL = FLOAT(LL)
         AJMIN = ABS(ABS(SPI-FL)-0.5)
         AJMAX = SPI +  FL + 0.5
         NUMJ = IFIX(AJMAX-AJMIN) + 1
         AJC = AJMIN - 1.
         IF(LL.NE.0.AND.(FL.GT.SPI-.5.AND.FL.LE.SPI)) THEN
            JJL = 0
         ELSE
            JJL = 1
         END IF
         DO JJ=1,NUMJ
            AJC = AJC + 1.
            GJ = (2.0*AJC+1.0)/GJD
!
!           LOOP OVER POSSIBLE CHANNEL SPINS
!
            DO KCHANL=1,2
               INOW = INOWB
               KPSTV = 0
               KNGTV = 0
!**************INITIALIZE MATRIX
               S = 0.0
               R = 0.0
               DO J=1,3
                  R(J,J) = 1.0
               END DO
!
!              LOOP OVER ALL RESONANCES
!
               DO I=1,NRS
                  AJPM = RESPAR(INOW+1)
                  AJ = ABS(AJPM)
!*****************SELECT ONLY RESONANCES WITH CURRENT J VALUE
                  IF(ABS(AJ-AJC).LE.EPI2)   THEN
                     IF(AJPM.LT.0.0) THEN
                        KNGTV = KNGTV + 1
                     ELSE IF(AJPM.GT.0.0) THEN
                        KPSTV = KPSTV + 1
                     END IF
                     IF(KCHANL.EQ.1.AND.AJPM.LT.0.0) THEN
                        ISKIP = 1
                     ELSE IF(KCHANL.EQ.2.AND.AJPM.GT.0.0) THEN
                        ISKIP = 1
                     ELSE
                        ISKIP = 0
                     END IF
!
!                    PROCESS RESONANCE
!
                     IF(ISKIP.EQ.0) THEN
                        ER = RESPAR(INOW)
                        GN = RESPAR(INOW+3)
                        GG = RESPAR(INOW+4)
                        GF = RESPAR(INOW+5)
                        PER = RESPAR(INOW+7)
                        GC = RESPAR(INOW+8)
                        A1 = SQRT(GN*PE/PER)
                        IF(GF.NE.0.0) THEN
                           A2 = SQRT(ABS(GF))
                           IF(GF.LT.0.)  A2 = -A2
                        ELSE
                           A2 = 0.0
                        END IF
                        IF(GC.NE.0.0) THEN
                           A3 = SQRT(ABS(GC))
                           IF(GC.LT.0.)   A3 = -A3
                        ELSE
                           A3 = 0.0
                        END IF
                        DIFF = ER - E
                        DEN = DIFF*DIFF + 0.25*GG*GG
                        DE2 = 0.5*DIFF/DEN
                        GG4 = 0.25*GG/DEN
                        IF(GF.EQ.0.0.AND.GC.EQ.0.0)   THEN
                           IFIS = 0
                        ELSE
                           IFIS = 1
                        END IF
!
!                       CALCULATE UPPER TRIANGULAR MATRIX TERMS
!
                        R(1,1) = R(1,1) + GG4*A1*A1
                        S(1,1) = S(1,1) - DE2*A1*A1
                        IF(IFIS.GT.0)   THEN
                           R(1,2) = R(1,2) + GG4*A1*A2
                           S(1,2) = S(1,2) - DE2*A1*A2
                           R(1,3) = R(1,3) + GG4*A1*A3
                           S(1,3) = S(1,3) - DE2*A1*A3
                           R(2,2) = R(2,2) + GG4*A2*A2
                           S(2,2) = S(2,2) - DE2*A2*A2
                           R(3,3) = R(3,3) + GG4*A3*A3
                           S(3,3) = S(3,3) - DE2*A3*A3
                           R(2,3) = R(2,3) + GG4*A2*A3
                           S(2,3) = S(2,3) - DE2*A2*A3
                        END IF
                      END IF
                   END IF
                   INOW = INOW + 9
               END DO
               KKKKKK = 0
               IF(KCHANL.EQ.1) THEN
                  IF(KPSTV.GT.0)   THEN
                     IF(KNGTV.EQ.0) THEN
                        IF(JJ.GT.JJL.AND.JJ.LT.NUMJ) THEN
                           KKKKKK = 2
                        ELSE
                           KKKKKK = 1
                        END IF
                     ELSE IF(KNGTV.GT.0) THEN
                        KKKKKK = 1
                     END IF
                  ELSE IF(KPSTV.EQ.0) THEN
                     IF(KNGTV.EQ.0) THEN
                        IF(JJ.GT.JJL.AND.JJ.LT.NUMJ) THEN
                           KKKKKK = 2
                        ELSE
                           KKKKKK = 1
                        END IF
                     ELSE IF(KNGTV.GT.0) THEN
                        KKKKKK = 0
                     END IF
                  END IF
               ELSE IF(KCHANL.EQ.2) THEN
                  IF(KPSTV.GT.0) THEN
                     IF(KNGTV.GT.0) THEN
                        KKKKKK = 1
                     END IF
                  ELSE IF(KPSTV.EQ.0) THEN
                     IF(KNGTV.GT.0) THEN
                        IF(JJ.GT.JJL.AND.JJ.LT.NUMJ) THEN
                           KKKKKK = 2
                        ELSE
                           KKKKKK = 1
                        END IF
                     END IF
                  END IF
               END IF
!
!              PREPARE TO CALCULATE THE CROSS SECTION
!
               IF(KKKKKK.NE.0) THEN
                  IF(IFIS.EQ.0)    THEN
!********************SPECIAL CASE OF NO FISSION
                     DEN = R(1,1)*R(1,1) + S(1,1)*S(1,1)
                     RI(1,1) = R(1,1)/DEN
                     SI(1,1) = -S(1,1)/DEN
                  ELSE
!********************MAKE SYMMETRIC MATRIX
                     R(2,1) = R(1,2)
                     S(2,1) = S(1,2)
                     R(3,1) = R(1,3)
                     S(3,1) = S(1,3)
                     R(3,2) = R(2,3)
                     S(3,2) = S(2,3)
!********************INVERT THE COMPLEX MATRIX
                     CALL FROBNS(R,S,RI,SI)
                  END IF
!
!                 CALCULATE THE CROSS SECTIONS
!
                  U11R = P1*(2.0*RI(1,1)-1.0) + 2.0*P2*SI(1,1)
                  U11I = P2*(1.0-2.0*RI(1,1)) + 2.0*P1*SI(1,1)
                  TERMT = 2.0*GJ*(1.0-U11R)
                  TERMN = GJ*((1.0-U11R)**2 + (U11I)**2)
                  IF(KKKKKK.EQ.2) THEN
                     TERMN = TERMN + 2.*GJ*(1.-P1)
                     TERMT = TERMT + 2.*GJ*(1.-P1)
                  END IF
                  IF(IFIS.EQ.0)  THEN
                     TERMF = 0.0
                  ELSE
                     T1 = RI(1,2)
                     T2 = SI(1,2)
                     T3 = RI(1,3)
                     T4 = SI(1,3)
                     TERMF = 4.0*GJ*(T1*T1+T2*T2+T3*T3+T4*T4)
                  END IF
                  TERMG = TERMT-TERMF-TERMN
                  SIGNNI = SIGNNI + TERMN
                  SIGNGI = SIGNGI + TERMG
                  SIGNFI = SIGNFI + TERMF
                  SIGNTI = SIGNTI + TERMT
               END IF
            END DO
         END DO
      END DO
!
!     CALCULATE CROSS SECTION AND STORE FOR RETURN
!
      SIGP(1) = PIFAC*SIGNTI
      SIGP(2) = PIFAC*SIGNNI
      SIGP(3) = PIFAC*SIGNFI
      SIGP(4) = PIFAC*SIGNGI
!
      RETURN
      END SUBROUTINE CSRMAT
!
!***********************************************************************
!
      SUBROUTINE FROBNS(A,B,C,D)
!
!     THIS SUBROUTINE INVERTS A COMPLEX MATRIX WITH REAL AND IMAGINARY
!     PARTS A AND B AND GIVES C AND D THE REAL AND IMAGINARY PARTS OF
!     THE INVERSE. FROBENIUS-SCHUR METHOD OF INVERSION.
!
      IMPLICIT NONE
!
      REAL(KIND=8), DIMENSION(3,3) :: A,B,C,D
!
      INTEGER(KIND=I4) :: IND
      REAL(KIND=8), DIMENSION(3,3) :: Q(3,3)
!
      C = A
      CALL THRINV(A,3,IND)
      IF(IND.NE.1) THEN
         CALL ABCMAT(A,B,Q)
         CALL ABCMAT(B,Q,D)
         C = C + D
         CALL THRINV(C,3,IND)
         IF(IND.NE.1) THEN
            CALL ABCMAT(Q,C,D)
            D = -D
         END IF
      END IF
!
      RETURN
      END SUBROUTINE FROBNS
!
!***********************************************************************
!
      SUBROUTINE THRINV(D,N,KIMERR)
!
!     INVERTS SYMMETRIC MATRIX (D(I,J),J=1,N,I=1,J)
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N,KIMERR
      REAL(KIND=8), DIMENSION(3,3) :: D
!
      INTEGER :: I,J,LR
      REAL(KIND=8) :: FOOEY
      REAL(KIND=8), DIMENSION(3) ::  S
!
      KIMERR = 0
      DO J=1,N
         DO I=1,J
            D(I,J) = -D(I,J)
            D(J,I) = D(I,J)
         END DO
         D(J,J) = 1.D+0 + D(J,J)
      END DO
      DO LR=1,N
         FOOEY = 1.D+0 - D(LR,LR)
         IF(FOOEY.EQ.0.0) THEN
            KIMERR = 1
            GO TO 100
         END IF
         D(LR,LR) = 1.D+0/FOOEY
         DO J=1,N
            S(J) = D(LR,J)
            IF(J.NE.LR)    THEN
               D(J,LR) = D(J,LR)*D(LR,LR)
               D(LR,J) = D(J,LR)
            END IF
         END DO
         DO J=1,N
            IF(J.NE.LR)   THEN
               DO I=1,J
                  IF(I.NE.LR)   THEN
                     D(I,J) = D(I,J) + D(I,LR)*S(J)
                     D(J,I) = D(I,J)
                  END IF
               END DO
            END IF
         END DO
      END DO
!
  100 RETURN
      END SUBROUTINE THRINV
!
!***********************************************************************
!
      SUBROUTINE ABCMAT(A,B,C)
!
!     ROUTINE TO DO A MATRIX MULTIPLICATION
!
      IMPLICIT NONE
!
      REAL(KIND=8), DIMENSION(3,3) :: A,B,C
!
      INTEGER(KIND=I4) :: I,J,K
!
      DO I=1,3
         DO J=1,3
            C(I,J)=0.0
            DO K=1,3
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE ABCMAT
!
!***********************************************************************
!
      SUBROUTINE CSAA(E,SIGP,INOW)
!
! **********************************************************************
! *   CALCULATES MULTI-LEVEL ADLER-ADLER  CROSS SECTIONS AT ENERGY E   *
! *   FOR ONE SECTION   (ONE ISOTOPE-ONE ENERGY RANGE)                 *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT, COS, SIN
!
      INTEGER(KIND=I4) :: NLS,NJS,NLJ,LI,IST
      INTEGER(KIND=I4) :: I,J,L,N
      REAL(KIND=R4) :: AWRI,AP,ARAT,EROOT,ESQ,ECUB
      REAL(KIND=R4) :: K,C,OMG,SNF,CSF
      REAL(KIND=R4) :: BAKT,BAKF,BAKC
      REAL(KIND=R4) :: DE,DW,GR,GI,T1,T2
!
!     ZERO OUT CROSS SECTION CONTRIBUTION FROM THIS SECTION (SIGP)
!
      SIGP = 0.0
!
!     RETRIEVE NUCLIDE INFORMATION
!
      AWRI = RESPAR(INOW+5)
      AP = RESPAR(INOW+3)
      NLS = RESPAR(INOW+4)
      LI = RESPAR(INOW+6)
      INOW = INOW + 8
!
!     CALCULATE CONSTANTS
!
      ARAT = AWRI/(AWRI+1.0)
      EROOT = SQRT(ABS(E))
      ESQ = E*E
      ECUB = ESQ*E
      K = CROC*ARAT*EROOT
      C = 6.5099897E+5/(ARAT*ARAT)
      OMG = 2.*K*AP
      SNF = SIN(OMG)
      CSF = COS(OMG)
!
!     CALCULATE BACKGROUND
!
      BAKT = 0.0
      BAKF = 0.0
      BAKC = 0.0
      IF(LI.LE.4) GO TO 100
      IF(LI.EQ.5.OR.LI.EQ.7) THEN
         BAKT = RESPAR(INOW) + RESPAR(INOW+1)/E + RESPAR(INOW+2)/ESQ
         BAKT = BAKT + RESPAR(INOW+3)/ECUB + RESPAR(INOW+4)*E +         &       
     &                 RESPAR(INOW+5)*ESQ
      END IF
      INOW = INOW + 6
      IF(LI.NE.5)   THEN
         BAKF = RESPAR(INOW) + RESPAR(INOW+1)/E + RESPAR(INOW+2)/ESQ
         BAKF = BAKF + RESPAR(INOW+3)/ECUB + RESPAR(INOW+4)*E +         &       
     &          RESPAR(INOW+5)*ESQ
      END IF
      INOW = INOW + 6
      BAKC = RESPAR(INOW) + RESPAR(INOW+1)/E + RESPAR(INOW+2)/ESQ
      BAKC = BAKC + RESPAR(INOW+3)/ECUB + RESPAR(INOW+4)*E +            &       
     &          RESPAR(INOW+5)*ESQ
      INOW = INOW + 6
!
!     CALCULATE RESONANCE CONTRIBUTION
!
      DO L=1,NLS
         NJS = RESPAR(INOW+1)
         INOW = INOW + 2
         DO J=1,NJS
            NLJ = RESPAR(INOW+1)
            INOW = INOW + 2
!
!           LOOP OVER INDIVIDUAL RESONANCES FOR THIS L-J STATE
!
            DO N=1,NLJ
               IST = INOW
               DO I=1,4
                  IF(I.NE.2)   THEN
                     DE = RESPAR(IST)
                     DW = RESPAR(IST+1)
                     GR = RESPAR(IST+2)
                     GI = RESPAR(IST+3)
                     IST = IST+4
                     T2 = (DE-E)**2+DW*DW
                     IF(I.LE.1)   THEN
                        T1 = DW*(GR*CSF+GI*SNF)+(DE-E)*(GI*CSF-GR*SNF)
                     ELSE
                        T1 = DW*GR+(DE-E)*GI
                     END IF
                     SIGP(I) = SIGP(I)+T1/T2
                  END IF
               END DO
               INOW = IST
            END DO
         END DO
      END DO
!
!     ADD BACKGROUND AND FACTORS
!
      SIGP(1) = (SIGP(1)+BAKT)*C/EROOT
      SIGP(1) = SIGP(1)+2.0*C*(1.0-CSF)/E
      SIGP(3) = (SIGP(3)+BAKF)*C/EROOT
      SIGP(4) = (SIGP(4)+BAKC)*C/EROOT
      IF(LI.NE.6)   THEN
         SIGP(2) = SIGP(1) - (SIGP(3)+SIGP(4))
      ELSE
         SIGP(1) = SIGP(3) + SIGP(4)
      END IF
!
  100 RETURN
      END SUBROUTINE CSAA
!
!***********************************************************************
!
      SUBROUTINE CSHYBR(E,SIGP,INOW)
!
! **********************************************************************
! *   CALCULATES HYBRID R-FUNCTION CROSS SECTIONS AT ENERGY E          *
! *   FOR ONE SECTION   (ONE ISOTOPE-ONE ENERGY RANGE)                 *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      REAL(KIND=R4), INTRINSIC :: ABS, SQRT
!
      INTEGER(KIND=I4) :: NLS,NJS,NSS,NRS,LL
      INTEGER(KIND=I4) :: LBK,LPS,ITPT
      INTEGER(KIND=I4) :: J,I2,INOWR,INOWN
      INTEGER(KIND=I4) :: I,II,III,L,N,NS,NJ
      INTEGER(KIND=I4) :: NGRE,NFRE,NIRE,NCRE,NRCHAN
      INTEGER(KIND=I4), DIMENSION(4) :: LC
      INTEGER(KIND=I4), DIMENSION(4,4) :: ICRPT
      REAL(KIND=R4) :: K,KFAC,AWRI,ARAT,PIFAC,RHRES,SABS
      REAL(KIND=R4) :: SE,PE,SD,PD,SN,PN,SER,PER
      REAL(KIND=R4) :: RHO,PHRO,PHRI,RR0,RI0,PHIR,PHRES,PHEN
      REAL(KIND=R4) :: EEFF,EREFF
      REAL(KIND=R4) :: ER,GN,GG,GF,GELIM
      REAL(KIND=R4) :: AC,AJ,GJ,SPI
      REAL(KIND=R4) :: GNE,GNR,GTT,EDELT,DEN,COMFAC
      REAL(KIND=R4), DIMENSION(4) :: QRE,GC,SIGC
      COMPLEX(KIND=R4) :: ONE,UNIT,RFCNN,SLSJNN,PHI0,PHI,R0,EXPHI
!
      UNIT = CMPLX(0.0,1.0)
      ONE = CMPLX(1.0,0.0)
!
!     ZERO OUT CROSS SECTION CONTRIBUTION FROM THIS SECTION (SIGP)
!
      SIGP = 0.0
!
!     RETRIEVE NUCLIDE INFORMATION
!
      SPI = RESPAR(INOW+2)
      NLS = RESPAR(INOW+4)
!
!     RETRIEVE COMPETING CHANNEL DESCRIPTIONS
!
      NGRE = RESPAR(INOW+5)
      NFRE = RESPAR(INOW+6)
      NIRE = RESPAR(INOW+7)
      NCRE = RESPAR(INOW+8)
      NRCHAN = NIRE + NCRE
      J = INOW + 8
      K = J + 4
      DO I=1,4
         I2 = I + K
         QRE(I) = RESPAR(I2)
      END DO
      IF(NRCHAN.GT.0)  THEN
         DO N=1,NRCHAN
            SIGC(N) = 0.
         END DO
      END IF
!*****SAVE POINTER TO BEGINNING OF RESONANCE PARAMETERS
      INOWR = RESPAR(INOW+17)
!*****MAKE POINTER INDEX FOR CHARGED PARTICLE PENETRABILITIES
      INOW = INOW + 18
      IF(NCRE.NE.0)    THEN
         DO N=1,NCRE
            DO L=1,4
               ICRPT(L,N) = RESPAR(INOW)
               INOW = INOW + 1
            END DO
         END DO
      END IF
!
!     CALCULATE WAVE NUMBER(K) AT ENERGY (E)
!
      INOW = INOWR
      AWRI = RESPAR(INOW)
      ARAT = AWRI/(AWRI+1.0)
      KFAC = CROC*ARAT
      K = KFAC*SQRT(ABS(E))
      PIFAC = PI/(K*K)
!
!     LOOP OVER L STATES
!
      DO L=1,NLS
         LL = RESPAR(INOW+1)
         NSS = RESPAR(INOW+2)
         INOW = INOW + 3
!
!        LOOP OVER CHANNEL SPIN STATES
!
         DO NS=1,NSS
            NJS = RESPAR(INOW+1)
            INOW = INOW + 2
!
!           LOOP OVER J STATES
!
            DO NJ=1,NJS
               AJ = RESPAR(INOW)
               GJ = 0.5*(2.0*AJ+1.0)/(2.0*SPI+1.0)
               AC = RESPAR(INOW+1)
               LBK = RESPAR(INOW+2)
               LPS = RESPAR(INOW+3)
               NRS = RESPAR(INOW+4)
               INOW = INOW + 5
               RHO = K*AC
!**************CALCULATE HARD SPHERE PHASE SHIFT
               CALL FACPHI(LL,RHO,PHIR)
               PHI0 = CMPLX(PHIR,0.0)
!**************CALCULATE SHIFT AND PENETRATION FACTORS AT CROSS SECTION
!**************ENERGY
               CALL FACTS(LL,RHO,SE,PE)
!
!              LOOP OVER ALL RESONANCES
!
               RFCNN = CMPLX(0.,0.)
               IF(NRS.NE.0)   THEN
                  DO I=1,NRS
                     ER = RESPAR(INOW)
                     GN = RESPAR(INOW+1)
                     GG = RESPAR(INOW+2)
                     GF = RESPAR(INOW+3)
                     INOW = INOW + 3
                     GELIM = GG + GF
                     IF(NRCHAN.GT.0)   THEN
!***********************DETERMINE VALUE OF COMPETITIVE WIDTHS
                        DO II=1,NRCHAN
                           GC(II) = 0.0
                           LC(II) = IFIX(RESPAR(INOW+II+4))
                           EEFF = E + (QRE(II)/ARAT)
                           IF(EEFF.GT.0.0)   THEN
                              EREFF = ABS(ER + (QRE(II)/ARAT))
!*****************************INELASTIC PENETRABILITY
                              IF(II.LE.NIRE) THEN
                                 PHRES = KFAC*SQRT(EREFF)*AC
                                 CALL FACTS(LC(II),PHRES,SD,PD)
                                 PHEN = KFAC*SQRT(EEFF)*AC
                                 CALL FACTS(LC(II),PHEN,SN,PN)
!*****************************CHARGED PARTICLE PENETRABILITY
                              ELSE
                                 IF(LC(II)+1.GT.4)  GO TO 25
                                 ITPT = ICRPT(LC(II)+1,II)
                                 CALL TERPRP(ITPT,EREFF,PD)
                                 CALL TERPRP(ITPT,EEFF,PN)
                              END IF
                              GC(II) = PN*RESPAR(INOW+II)/PD
                           END IF
   25                      GELIM = GELIM + GC(II)
                        END DO
                     END IF
                     INOW = INOW + 9
!
!                    CROSS SECTION CALCULATIONS
!
                     RHRES = KFAC*SQRT(ABS(ER))*AC
                     CALL FACTS(LL,RHRES,SER,PER)
                     GNE = GN*PE/PER
                     GNR = GN/(2.*PER)
                     GTT = GNE + GELIM
                     EDELT = ER - E
                     DEN = EDELT**2+0.25*GTT*GTT
                     COMFAC = PIFAC*GJ*GNE/DEN
                     RFCNN = RFCNN + GNR/(EDELT-0.5*UNIT*GELIM)
                     IF(NGRE.NE.0) SIGP(4) = SIGP(4)+COMFAC*GG
                     IF(NFRE.NE.0) SIGP(3) = SIGP(3)+COMFAC*GF
                     IF(NRCHAN.GT.0) THEN
                        DO III=1,NRCHAN
                           SIGC(III) = SIGC(III) + COMFAC*GC(III)
                        END DO
                     END IF
                  END DO
!
!                 SAVE POINTER TO BEGINNING OF NEXT RP SET
!
                  INOWN = RESPAR(INOW)
!*****************BACKGROUND R-FUNCTION
                  INOW = INOW + 1
                  R0 = CMPLX(0.,0.)
                  IF (LBK.NE.0)  THEN
                     CALL TERPRP(INOW,EEFF,RR0)
                     CALL TERPRP(INOW,EEFF,RI0)
                     R0 = CMPLX(RR0,RI0)
                  END IF
                  RFCNN = RFCNN + R0
!*****************OPTICAL MODEL PHASE SHIFT
                  PHI = PHI0
                  IF (LPS.NE.0)  THEN
                     CALL TERPRP(INOW,EEFF,PHRO)
                     CALL TERPRP(INOW,EEFF,PHRI)
                     PHI = CMPLX(PHRO,PHRI)
                  END IF
!*****************R-FUNCTION CALCULATIONS
                  EXPHI = CEXP(-2.0*UNIT*PHI)
                  SLSJNN = (ONE-EXPHI) - UNIT*PE*RFCNN*(ONE+EXPHI)
                  SLSJNN = SLSJNN/(ONE-UNIT*PE*RFCNN)
                  SABS = CABS(SLSJNN)
                  SIGP(1) = SIGP(1) + 2.*PIFAC*GJ*REAL(SLSJNN)
                  SIGP(2) = SIGP(2) + PIFAC*GJ*SABS*SABS
                  SIGP(1) = SIGP(2) + SIGP(3) + SIGP(4)
                  INOW = INOWN
               END IF
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE CSHYBR
!
!***********************************************************************
!
      SUBROUTINE CSUNR1(E,SIGP,INOW)
!
! **********************************************************************
! *   UNRESOLVED RESONANCE REGION (FORMAT 2) PROGRAM                   *
! *   SINGLE LEVEL BREIT WIGNER FORMALISM                              *
! *   ENERGY INDEPENDENT PARAMETERS                                    *
! *      LFW=0  NO FISSION WIDTHS GIVEN                                *
! *      LFW=1  FISSION WIDTHS TABULATED                               *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      REAL(KIND=R4), INTRINSIC :: SQRT
!
      INTEGER(KIND=I4) :: LFW,LSSF,MUF
      INTEGER(KIND=I4) :: NE,NLS,NJS,LL
      INTEGER(KIND=I4) :: NUME,INEXT,NFS,IND,IDG
      INTEGER(KIND=I4) :: MU,NU,LAMDA
      INTEGER(KIND=I4) :: J,L,N,IN
      REAL(KIND=R4) :: SPI,A
      REAL(KIND=R4) :: AWRI,AA,RAT,CONST
      REAL(KIND=R4) :: ENU,ER2,K,RHO,RHOC,PS,VL
      REAL(KIND=R4) :: GFX,GXX,DX,AJ,GJ,AMUN,GNOX,GGX,GNX
      REAL(KIND=R4) :: GC,GS,GFF
      REAL(KIND=R4) :: DIFF,DEN,TEMP,TERG,TERS,TERF,SPOT
      REAL(KIND=R4), DIMENSION(2) :: EG
      REAL(KIND=R4), DIMENSION(2,3) :: SIGINT
!
!     RETRIEVE NUCLIDE INFORMATION
!
      LFW = RESPAR(INOW+2)
      SPI = RESPAR(INOW+3)
      LSSF = RESPAR(INOW+4)
      A = RESPAR(INOW+5)
      NE = RESPAR(INOW+6)
      NLS = RESPAR(INOW+7)
      INOW = INOW + 8
      IF(LSSF.EQ.1)  GO TO 100
!
!     FIND ENERGIES AT WHICH UR CROSS SECTIONS MUST BE CALCULATED
!
      NUME = 1
      EG(1) = E
      IF(LFW.NE.0)   THEN
         CALL FINDE(E,INOW,NE,1,EG,NUME,IDG)
         NFS = 6*((NE+5)/6)
         INOW = INOW + NFS
      END IF
!
!     CALCULATE CHANNEL RADIUS
!
      AWRI = RESPAR(INOW)
      AA = RESPAR(INOW+1)
      RAT = AWRI/(AWRI+1.0)
      CONST = (2.0*PI*PI)/(CROC*RAT)**2
      INOW = INOW+2
      INEXT = INOW
!
!     CALCULATE AT ALL ENERGIES
!
      DO N=1,NUME
         SIGP = 0.
         INOW = INEXT
         ENU = EG(N)
         ER2 = SQRT(ENU)
         K = CROC*RAT*ER2
         RHO = K*AA
         RHOC = K*A
!
!        DO LOOP OVER ALL L STATES
!
         DO L=1,NLS
            NJS = RESPAR(INOW+1)
            LL = RESPAR(INOW)
            INOW = INOW+2
!
!           DO LOOP OVER ALL J STATES
!
            DO J=1,NJS
               GFX = 0.0
               GXX = 0.0
               IF(LFW.EQ.0)   THEN
                  MUF= 0
                  DX = RESPAR(INOW)
                  AJ = RESPAR(INOW+1)
                  AMUN = RESPAR(INOW+2)
                  GNOX = RESPAR(INOW+3)
                  GGX = RESPAR(INOW+4)
                  INOW = INOW + 6
               ELSE
                  MUF = RESPAR(INOW)
                  DX = RESPAR(INOW+1)
                  AJ = RESPAR(INOW+2)
                  AMUN = RESPAR(INOW+3)
                  GNOX = RESPAR(INOW+4)
                  GGX = RESPAR(INOW+5)
                  INOW = INOW + 7
                  IND = INOW + IDG + N - 2
                  GFX = RESPAR(IND)
                  INOW = INOW + NFS
               END IF
               GJ = (2.0*AJ+1.0)/(4.0*SPI+2.0)
               MU = AMUN
               NU = MUF
               LAMDA = 0
!**************CALCULATE PENETRABILITY (VL) AND PHASE SHIFT(PS)
               CALL UNFAC(LL,RHO,RHOC,AMUN,VL,PS)
               VL=VL*ER2
!
!              CALCULATE POTENTIAL SCATTERING
!
               IF(J.LE.1)  THEN
                  SPOT = 4.0*PI*(2.0*FLOAT(LL)+1.0)*(SIN(PS)/K)**2
               END IF
               GNX = GNOX*VL
               DIFF = GXX
               DEN = ENU*DX
               TEMP = CONST*GJ*GNX/DEN
               TERG = TEMP*GGX
               TERS = TEMP*GNX
               TERF = TEMP*GFX
!**************CALCULATE FLUCTUATION INTEGRALS
               CALL GNRL(GNX,GFX,GGX,MU,NU,LAMDA,GS ,DIFF,1)
               CALL GNRL(GNX,GFX,GGX,MU,NU,LAMDA,GC ,DIFF,2)
               CALL GNRL(GNX,GFX,GGX,MU,NU,LAMDA,GFF,DIFF,3)
               GC = GC*TERG
               GFF = GFF*TERF
               GS = GS*TERS - CONST*GJ*2.0*GNX*(SIN(PS)**2)/(ENU*DX)
!
!              CROSS SECTIONS
!
               SIGP(2) = SIGP(2) + GS
               SIGP(3) = SIGP(3) + GFF
               SIGP(4) = SIGP(4) + GC
            END DO
!
!           ADD POTENTIAL SCATTERING
!
            SIGP(2) = SIGP(2) + SPOT
         END DO
!
!        INTERPOLATE TO REQUIRED ENERGY
!
         IF(NUME.EQ.1)   GO TO 50
         DO IN=1,3
            SIGINT(N,IN) = SIGP(IN+1)
         END DO
      END DO
      DO IN=1,3
         CALL TERP1(EG(1),SIGINT(1,IN),EG(2),SIGINT(2,IN),E,            &       
     &              SIGP(IN+1),2)
      END DO
!
!     SUM TO GET TOTAL
!
   50 SIGP(1) = SIGP(2) + SIGP(3) + SIGP(4)
!
  100 RETURN
      END SUBROUTINE CSUNR1
!
!***********************************************************************
!
      SUBROUTINE CSUNR2(E,SIGP,INOW)
!
! **********************************************************************
! *   UNRESOLVED RESONANCE REGION (FORMAT 2) PROGRAM                   *
! *   SINGLE LEVEL BREIT WIGNER FORMALISM                              *
! *   ENERGY   DEPENDENT PARAMETERS                                    *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(4) :: SIGP
!
      INTEGER(KIND=I4) :: NLS,NJS,LSSF,NE,LL
      INTEGER(KIND=I4) :: INTT,NUME,IDG,IND,INEXT
      INTEGER(KIND=I4) :: MU,NU,LAMDA
      INTEGER(KIND=I4) :: J,L,N,IN
      REAL(KIND=R4) :: SPI,A
      REAL(KIND=R4) :: ENU,ER2,K,RHO,RHOC,PS,VL
      REAL(KIND=R4) :: AWRI,AA,RAT,CONST
      REAL(KIND=R4) :: AJ,AMUX,AMUN,AMUF
      REAL(KIND=R4) :: GJ,DX,GXX,GNOX,GNX,GGX,GFX,GC,GS,GFF,SPOT
      REAL(KIND=R4) :: GAMMA,GALPHA,GBETA
      REAL(KIND=R4) :: DIFF,DEN,TEMP,TERG,TERS,TERF
      REAL(KIND=R4), DIMENSION(2) :: EG
      REAL(KIND=R4), DIMENSION(2,3) :: SIGINT
!
!     RETRIEVE NUCLIDE INFORMATION
!
      SPI = RESPAR(INOW+3)
      LSSF = RESPAR(INOW+4)
      A = RESPAR(INOW+5)
      NLS = RESPAR(INOW+6)
      AWRI = RESPAR(INOW+7)
      AA = RESPAR(INOW+8)
      INOW = INOW+9
      IF(LSSF.EQ.1)  GO TO 100
!
!     CALCULATE CHANNEL RADIUS
!
      RAT=AWRI/(AWRI+1.0)
      CONST=(2.0*PI*PI)/(CROC*RAT)**2
      INTT = RESPAR(INOW+3)
!
!     FIND ENERGIES AT WHICH UR CROSS SECTION MUST BE CALCULATED
!
      NUME = 1
      EG(1) = E
      NE = RESPAR(INOW+4)
      CALL FINDE(E,INOW+11,NE,6,EG,NUME,IDG)
      INEXT = INOW
!
!     CALCULATE AT ALL ENERGIES
!
      DO N=1,NUME
         SIGP = 0.
         INOW = INEXT
         ENU = EG(N)
         ER2 = SQRT(ENU)
         K = CROC*RAT*ER2
         RHO = K*AA
         RHOC = K*A
!
!        DO LOOP OVER ALL L STATES
!
         DO L=1,NLS
            NJS = RESPAR(INOW+1)
            LL = RESPAR(INOW)
            INOW = INOW + 2
!
!           DO LOOP  OVER ALL J STATES
!
            DO J=1,NJS
               AJ = RESPAR(INOW)
               AMUX = RESPAR(INOW+5)
               AMUN = RESPAR(INOW+6)
               AMUF = RESPAR(INOW+8)
               INOW = INOW + 9
               GJ = (2.0*AJ+1.0)/(4.0*SPI+2.0)
               MU = AMUN
               NU = AMUF
               LAMDA = AMUX
               IND = INOW + IDG + 6*(N-1) - 1
               DX = RESPAR(IND+1)
               GXX = RESPAR(IND+2)
               GNOX = RESPAR(IND+3)
               GGX = RESPAR(IND+4)
               GFX = RESPAR(IND+5)
!**************CALCULATE PENETRABILITY (VL) AND PHASE SHIFT(PS)
               CALL UNFAC(LL,RHO,RHOC,AMUN,VL,PS)
               VL = VL*ER2
!
!              CALCULATE POTENTIAL SCATTERING
!
               IF(J.LE.1)   THEN
                  SPOT = 4.0*PI*(2.0*FLOAT(LL)+1.0)*(SIN(PS)/K)**2
               END IF
               GNX = GNOX*VL
               GAMMA = GGX
               GALPHA = GNX
               GBETA = GFX
               DIFF = GXX
               DEN = ENU*DX
               TEMP = CONST*GJ*GNX/DEN
               TERG = TEMP*GAMMA
               TERS = TEMP*GNX
               TERF = TEMP*GBETA
!**************CALCULATE FLUCTUATION INTEGRALS
               CALL GNRL(GALPHA,GBETA,GAMMA,MU,NU,LAMDA,GS ,DIFF,1)
               CALL GNRL(GALPHA,GBETA,GAMMA,MU,NU,LAMDA,GC ,DIFF,2)
               CALL GNRL(GALPHA,GBETA,GAMMA,MU,NU,LAMDA,GFF,DIFF,3)
               GC = GC*TERG
               GFF = GFF*TERF
               GS = GS*TERS - CONST*GJ*2.0*GNX*(SIN(PS)**2)/(ENU*DX)
!
!              CROSS SECTIONS
!
               SIGP(2) = SIGP(2) + GS
               SIGP(3) = SIGP(3) + GFF
               SIGP(4) = SIGP(4) + GC
               INOW = INOW + 6*NE
            END DO
!
!           ADD POTENTIAL SCATTERING
!
            SIGP(2) = SIGP(2) + SPOT
         END DO
!
!        INTERPOLATE TO REQUIRED ENERGY
!
         IF(NUME.EQ.1)    GO TO 50
         DO IN=1,3
            SIGINT(N,IN) = SIGP(IN+1)
         END DO
      END DO
      DO IN=1,3
        CALL TERP1(EG(1),SIGINT(1,IN),EG(2),SIGINT(2,IN),E,             &       
     &          SIGP(IN+1),INTT)
      END DO
!
!     SUM TO GET TOTAL
!
   50 SIGP(1) = SIGP(2) + SIGP(3) + SIGP(4)
!
  100 RETURN
      END SUBROUTINE CSUNR2
!
!***********************************************************************
!
      SUBROUTINE FACTS(L,RHO,SE,PE)
!
! **********************************************************************
! *   CALCULATES PENETRATION AND SHIFT FACTORS                         *
! *   LIBERATED FROM GREGSON ET. AL.                                   *
! **********************************************************************
!
!     IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L
      REAL(KIND=R4) :: RHO,SE,PE
!
      REAL(KIND=R4) :: R2,RR4,R6,RR8,DEN
!
      IF(L.EQ.0) THEN
!        S-WAVE
         SE = 0.0
         PE = RHO
      ELSE IF(L.EQ.1) THEN
!        P-WAVE
         R2 = RHO*RHO
         DEN = 1.0 + R2
         PE = R2*RHO/DEN
         SE = -1.0/DEN
      ELSE IF(L.EQ.2) THEN
!        D-WAVE (L=2)
         R2 = RHO*RHO
         RR4 = R2*R2
         DEN = 3.0*R2 + RR4 + 9.0
         PE = RR4*RHO/DEN
         SE = -(18.0+3.0*R2)/DEN
      ELSE IF(L.EQ.3) THEN
!        F-WAVE (L=3)
         R2 = RHO*RHO
         RR4 = R2*R2
         R6 = RR4*R2
         DEN = 225.0 + 45.0*R2 + 6.0*RR4 + R6
         PE = R6*RHO/DEN
         SE = -(675.0+90.0*R2+6.0*RR4)/DEN
      ELSE IF(L.EQ.4) THEN
!        G-WAVE (L=4)
         R2 = RHO*RHO
         RR4 = R2*R2
         R6 = RR4*R2
         RR8 = RR4*RR4
         DEN = 11025.0 + 1575.0*R2 + 135.0*RR4 + 10.0*R6 + RR8
         PE = RR8*RHO/DEN
         SE = -(44100.0+4725.0*R2+270.0*RR4+10.0*R6)/DEN
      END IF
!
      RETURN
      END SUBROUTINE FACTS
!
!***********************************************************************
!
      SUBROUTINE FACPHI(L,RHO,PHI)
!
! **********************************************************************
! *   CALCULATES PHASE SHIFT                                           *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L
      REAL(KIND=R4) :: RHO,PHI
!
      REAL(KIND=R4) :: R2,RR4,TOP,BOT
!
      IF(L.EQ.0) THEN
!        S-WAVE
         PHI = RHO
      ELSE IF(L.EQ.1) THEN
!        P-WAVE
         PHI = RHO - ATAN(RHO)
      ELSE IF(L.EQ.2) THEN
!        D-WAVE
         R2 = RHO*RHO
         PHI = RHO - ATAN(3.0*RHO/(3.0-R2))
         IF((PHI/RHO).LT.0.000001) PHI=0.0
      ELSE IF(L.EQ.3) THEN
!        F-WAVE
         R2 = RHO*RHO
         PHI = RHO - ATAN((15.0*RHO-RHO*R2)/(15.0-6.0*R2))
         IF((PHI/RHO).LT.0.000001) PHI=0.0
      ELSE IF(L.EQ.4) THEN
!        G-WAVE
         R2 = RHO*RHO
         RR4 = R2*R2
         TOP = 105.0*RHO-10.0*R2*RHO
         BOT = 105.0 - 45.0*R2 + RR4
         PHI = RHO - ATAN(TOP/BOT)
         IF((PHI/RHO).LT.0.000001) PHI=0.0
      END IF
!
      RETURN
      END SUBROUTINE FACPHI
!
!***********************************************************************
!
      SUBROUTINE UNFAC(L,RHO,RHOC,AMUN,VL,PS)
!
! **********************************************************************
! *   CALCULATES THE PENETRABILITY FACTOR (VL) AND PHASE SHIFT (PS)    *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L
      REAL(KIND=R4) :: RHO,RHOC,AMUN,VL,PS
!
      REAL(KIND=R4) :: R2,RR4
!
      IF(L.EQ.0) THEN
!        S-WAVE
         VL = AMUN
         PS = RHOC
      ELSE IF(L.EQ.1) THEN
!        P-WAVE
         R2 = RHO*RHO
         VL = AMUN*R2/(1.0+R2)
         PS = RHOC - ATAN(RHOC)
      ELSE IF(L.EQ.2) THEN
!        D-WAVE
         R2 = RHO*RHO
         RR4 = R2*R2
         VL = AMUN*RR4/(9.0+3.0*R2+RR4)
         PS = RHOC - ATAN(3.0*RHOC/(3.0-RHOC*RHOC))
      END IF
!
      RETURN
      END SUBROUTINE UNFAC
!
!***********************************************************************
!
      SUBROUTINE FINDE(E,INOW,NE,ISKIP,EG,NUME,IDG)
!
!     SUBROUTINE FINDS THE BRACKETING ENERGY VALUES ON THE
!     UNRESOLVED RESONANCE PARAMETER ENERGY GRID
!
!     E     ENERGY AT WHICH CROSS SECTION IS TO BE CALCULATED
!     INOW  FIRST GRID LOCATION IN PACKED STORAGE
!     NE    NUMBER OF GRID POINTS
!     ISKIP REPETITION RATE FOR VALUES OF THE
!           ENERGY GRID IN THE PACKED ARRAY, RESPAR
!     EG    ENERGY GRID VALUES RETURNED
!     NUME  NUMBER OF VALUES IN EG
!     IDG   POINTER TO EG(1) RELATIVE TO INOW
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW,NE,ISKIP,NUME,IDG
      REAL(KIND=R4) :: E
      REAL(KIND=R4), DIMENSION(2) :: EG
!
      INTEGER(KIND=I4) :: IND,NEXT
      INTEGER(KIND=I4) :: N
!
!     FIND RELATIVE LOCATION OF E ON THE GRID
!
      IDG=0
      DO N=1,NE
         IND = INOW + ISKIP*(N-1)
         IF(E.EQ.RESPAR(IND))  THEN
!
!           ON FALL THRU E IS GREATER THAN LAST POINT ON GRID. EG IS
!             IS RETURNED CONTAINING ONLY LAST GRID VALUE
!           ON TRANSFER E MATCHES A GRID POINT
!
            IDG = IND -INOW + 1
            GO TO 100
         ELSE IF(E.LT.RESPAR(IND))  THEN
!
!           BRACKETING VALUES FOUND
!
            IF(N.EQ.1)   THEN
!
!           ENERGY IS BELOW LOWEST GRID POINT. EG IS RETURNED WITH ONLY
!             THE LOWEST GRID POINT
!
               IDG = 1
               GO TO 100
            ELSE
               IDG  = IND - INOW - ISKIP + 1
               NEXT = IND - ISKIP
               EG(1)= RESPAR(NEXT)
               EG(2)= RESPAR(IND)
               NUME = 2
               GO TO 100
            END IF
         END IF
      END DO
!
!     BRACKETING INTERVAL NOT FOUND - RETURN LAST POINT
      IDG = IND -INOW + 1
!
  100 RETURN
      END SUBROUTINE FINDE
!
!***********************************************************************
!
      SUBROUTINE GNRL(GALPHA,GBETA,GAMMA,MU,NU,LAMBDA,S,DF,ID)
!
! **********************************************************************
! *   CALCULATES FLUCTUATION INTEGRALS FOR UNRESOLVED RESONANCES       *
! *   LIBERATED FROM AVERAGE4 BY M.BHAT                                *
! *   MODIFIED 9/1/76 TO USE THE WEIGHTING SCHEME FROM MCC
! *
! **********************************************************************
!
! **********************************************************************
! *   CONSTANTS FOR NUMERICAL INTEGRATION OF FLUCTUATION INTEGRALS     *
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MU,NU,LAMBDA,ID
      REAL(KIND=R4) :: GALPHA,GBETA,GAMMA,S,DF
!
      INTEGER(KIND=I4) :: J,K,L
      REAL(KIND=R4) :: XJ,WJ,XK,WK,XL,WL
!
      REAL(KIND=R4), DIMENSION(4,10) :: XX
      DATA XX/3.0013465E-03,1.3219203E-02,1.0004488E-03,1.3219203E-02,  &       
     &        7.8592886E-02,7.2349624E-02,2.6197629E-02,7.2349624E-02,  &       
     &        4.3282415E-01,1.9089473E-01,1.4427472E-01,1.9089473E-01,  &       
     &        1.3345267E+00,3.9528842E-01,4.4484223E-01,3.9528842E-01,  &       
     &        3.0481846E+00,7.4083443E-01,1.0160615E+00,7.4083443E-01,  &       
     &        5.8263198E+00,1.3498293E+00,1.9421066E+00,1.3498293E+00,  &       
     &        9.9452656E+00,2.5297983E+00,3.3150885E+00,2.5297983E+00,  &       
     &        1.5782128E+01,5.2384894E+00,5.2607092E+00,5.2384894E+00,  &       
     &        2.3996824E+01,1.3821772E+01,7.9989414E+00,1.3821772E+01,  &       
     &        3.6216208E+01,7.5647525E+01,1.2072069E+01,7.5647525E+01/
      REAL(KIND=R4), DIMENSION(4,10) :: WW
      DATA WW/1.1120413E-01,3.3773418E-02,3.3376214E-04,1.7623788E-03,  &       
     &        2.3546798E-01,7.9932171E-02,1.8506108E-02,2.1517749E-02,  &       
     &        2.8440987E-01,1.2835937E-01,1.2309946E-01,8.0979849E-02,  &       
     &        2.2419127E-01,1.7652616E-01,2.9918923E-01,1.8797998E-01,  &       
     &        1.0967668E-01,2.1347043E-01,3.3431475E-01,3.0156335E-01,  &       
     &        3.0493789E-02,2.1154965E-01,1.7766657E-01,2.9616091E-01,  &       
     &        4.2930874E-03,1.3365186E-01,4.2695894E-02,1.0775649E-01,  &       
     &        2.5827047E-04,2.2630659E-02,4.0760575E-03,2.5171914E-03,  &       
     &        4.9031965E-06,1.6313638E-05,1.1766115E-04,8.9630388E-10,  &       
     &        1.4079206E-08,0.0          ,5.0989546E-07,0.0/
!
      IF(NU.LT.1.OR.NU.GT.4) NU=4
      S = 0.0
      IF(GALPHA.GT.0.0) THEN
         IF(GAMMA.GT.0.0) THEN
            IF(GBETA.GT.0.0) THEN
               GO TO 30
            ELSE IF(GBETA.EQ.0.0) THEN
               IF(DF.EQ.0.0) THEN
                  GO TO 10
               ELSE IF(DF.GT.0.0) THEN
                  GO TO 20
               END IF
            END IF
         END IF
      END IF
      GO TO 100
!
   10 DO J=1,10
         XJ = XX(MU,J)
         WJ = WW(MU,J)
         IF(ID.EQ.1) THEN
            S = S + WJ*((XJ*XJ)/(GALPHA*XJ+GAMMA))
         ELSE IF(ID.EQ.2) THEN
            S = S + WJ*(XJ/(GALPHA*XJ+GAMMA))
         END IF
      END DO
      GO TO 100
!
   20 DO J=1,10
         XJ = XX(MU,J)
         WJ = WW(MU,J)
         DO K=1,10
            XK = XX(LAMBDA,K)
            WK = WW(LAMBDA,K)
            IF(ID.EQ.1) THEN
               S = S + WJ*WK*((XJ*XJ)/(GALPHA*XJ+GAMMA+DF*XK))
            ELSE IF(ID.EQ.2) THEN
               S = S+WJ*WK*(XJ/(GALPHA*XJ+GAMMA+DF*XK))
            END IF
         END DO
      END DO
      GO TO 100
!
   30 IF(DF.EQ.0.0) THEN
         DO J=1,10
            XJ = XX(MU,J)
            WJ = WW(MU,J)
            DO K=1,10
               XK = XX(NU,K)
               WK = WW(NU,K)
               IF(ID.EQ.1) THEN
                  S = S + WJ*WK*((XJ*XJ)/(GALPHA*XJ+GBETA*XK+GAMMA))
               ELSE IF(ID.EQ.2) THEN
                  S = S + WJ*WK*(XJ/(GALPHA*XJ+GBETA*XK+GAMMA))
               ELSE IF(ID.EQ.3) THEN
                  S = S + WJ*WK*((XJ*XK)/(GALPHA*XJ+GBETA*XK+GAMMA))
               END IF
            END DO
         END DO
      ELSE IF(DF.GT.0.0) THEN
         DO J=1,10
            XJ = XX(MU,J)
            WJ = WW(MU,J)
            DO K=1,10
               XK = XX(NU,K)
               WK = WW(NU,K)
               DO L=1,10
                  XL = XX(LAMBDA,L)
                  WL = WW(LAMBDA,L)
                  IF(ID.EQ.1) THEN
                     S = S + WJ*WK*WL*((XJ*XJ)/                         &       
     &                    (GALPHA*XJ+GBETA*XK+GAMMA+DF*XL))
                  ELSE IF(ID.EQ.2) THEN
                     S = S + WJ*WK*WL*(XJ/                              &       
     &                    (GALPHA*XJ+GBETA*XK+GAMMA+DF*XL))
                  ELSE IF(ID.EQ.3) THEN
                     S = S + WJ*WK*WL*((XJ*XK)/                         &       
     &                    (GALPHA*XJ+GBETA*XK+GAMMA+DF*XL))
                  END IF
               END DO
            END DO
         END DO
      END IF
!
  100 RETURN
      END SUBROUTINE GNRL
!
!***********************************************************************
!
      SUBROUTINE LEGCK
!
!     CHECKS LEGENDRE COEFFICIENTS FOR POSSIBILITY AND POSITIVITY
!
      IMPLICIT NONE
!
      NMOM = NCOEF
      CALL LEGEND
      CALL POSTIV
      CALL NEGTIV
!
      RETURN
      END SUBROUTINE LEGCK
!
!***********************************************************************
!
      SUBROUTINE LEGEND
!
!     CONVERTS LEGENDRE COEFFICIENTS TO MOMENTS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NFM1
      INTEGER(KIND=I4) :: L,N
      REAL(KIND=R4) :: P10,P20,FL,FNF
      REAL(KIND=R4), DIMENSION(NCOEFMAX) :: P1,P2
!
      MOMENT(1) = F(1)
      IF(NMOM.LE.1)   GO TO 100
!
      NFM1=NMOM-1
      P10 = 0.
      P1(1) = 1.
      DO L=2,NMOM
         P1(L)=0.
         MOMENT(L)=0
      END DO
      DO N=2,NMOM
         P20 = P1(1)/3.
         P2(1) = P10+.4*P1(2)
         IF(NFM1.GE.2)   THEN
            DO L=2,NFM1
               FL = L
               P2(L) = FL/(2.*FL-1.)*P1(L-1)+(FL+1.)/(2.*FL+3.)*P1(L+1)
            END DO
         END IF
         FNF = NMOM
         P2(NMOM) = FNF/(2.*FNF-1.)*P1(NFM1)
         MOMENT(N) = P20
         DO L=1,NMOM
            MOMENT(N) = MOMENT(N) + P2(L)*F(L)
            P1(L) = P2(L)
         END DO
         P10 = P20
      END DO
!
  100 RETURN
      END SUBROUTINE LEGEND
!
!***********************************************************************
!
      SUBROUTINE POSTIV
!
!     FOR LEGENDRE CHECKING ROUTINE. DETERMINES IF MOMENTS
!     ARE PHYSICALLY POSSIBLE.
!     N = NUMBER OF EVEN MOMENTS ACCEPTED
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N,ITEST
      INTEGER(KIND=I4) :: L
      INTEGER(KIND=I4) :: LL
      REAL(KIND=R4) :: QCK
!
      CALL GETMUS(N,ITEST)
      QCK = -1.
      IF(N.LE.0)   THEN
         IF(ITEST.EQ.0)  GO TO 100
      ELSE
         DO LL=1,N
            L = LL
            IF(QUE(L,1.).LT.0..or.QUE(L,-1.)*QCK.LT.0.)  THEN
               N = L - 1
               GO TO 90
            END IF
            QCK = -QCK
         END DO
         IF(ITEST.EQ.0)   GO TO 90
      END IF
      IF(QUE(N+1,1.).GE.0..AND.QUE(N+1,-1.)*QCK.GE.0.)   GO TO 100
   90 IF(N.GT.0)   THEN
         MEAN(N+1) = 1. - VAR(N)*QUE(N-1,1.)/QUE(N,1.)
         IF(QUE(N+1,-1.)*QCK.GE.0.)      GO TO 100
         N = N - 1
      END IF
!
  100 RETURN
      END SUBROUTINE POSTIV
!
!***********************************************************************
!
      SUBROUTINE GETMUS(N,ITEST)
!
!     CALCULATES MEANS AND VARIANCES FROM MOMENTS
!     FOR LEGENDRE CHECKING ROUTINE.
!
!     Half of maximum number of Legendre coefficients
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N,ITEST
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: IM1,IP1,NP1,IHF4
      INTEGER(KIND=I4) :: I,K
      REAL(KIND=R4) :: RMO1,RMO2,RMO3
      REAL(KIND=R4), DIMENSION(MXLG) :: MU,NORM,L,Q,SIG
      REAL(KIND=R4), DIMENSION(MXLG,MXLG+1) :: A
!
      RMO1 = MOMENT(1)
      RMO2 = MOMENT(2)
      RMO3 = MOMENT(3)
!
      DO I=1,MXLG
         MU(I) = 0.
         SIG(I) = 0.
         NORM(I) = 0.
         L(I) = 0.
         Q(I) = 0.
         MEAN(I) = 0.
         VAR(I) = 0.
      END DO
      MEAN(MXLG+1) = 0.
      A = 0.
      ITEST = MOD(NMOM,2)
!
!     ITEST = +1      ODD NUMBER OF MOMENTS ACCEPTED.
!           =  0     EVEN NUMBER OF MOMENTS ACCEPTED.
!           = -1     EVEN MOMENT REJECTED DUE TO NEGATIVE VARIANCE
!
      N = NMOM/2
      IF(N.GT.MXLG) STOP 'PSYCHE ERROR - MXLG limit exceeded'
      MU(1) = RMO1
      Q(1) = RMO1
      L(1) = RMO1
      A(1,2) = 1.
      A(1,1) = -Q(1)
      IF(N.LE.0)   GO TO 90
      SIG(1) = RMO2 - RMO1**2
      NORM(1) = SIG(1)
      IF(SIG(1).LE.0.)   THEN
         I = 1
         N = I - 1
         ITEST = -1
         GO TO 80
      END IF
      IF(N.GT.1)   THEN
         L(2) = RMO3-RMO1*RMO2
         Q(2) = L(2)/NORM(1)
         MU(2) = Q(2)-Q(1)
         A(2,3) = 1.
         A(2,2) = -Q(2)
         A(2,1) = (RMO1*RMO3-RMO2**2)/SIG(1)
         NORM(2) = MOMENT(4) + A(2,2)*RMO3 + A(2,1)*RMO2
         SIG(2) = NORM(2)/NORM(1)
         IF(SIG(2).LE.0.)  THEN
            I=2
            N = I - 1
            ITEST = -1
            GO TO 80
         END IF
         IF(N.GT.2)   THEN
            DO I=3,N
               IM1 = I - 1
               IP1 = I + 1
               DO K=1,I
                  IHF4 = IM1 + K
                  L(I) = L(I) + A(IM1,K)*MOMENT(IHF4)
               END DO
               Q(I) = L(I)/NORM(IM1)
               MU(I) = Q(I) - Q(IM1)
               A(I,IP1) = 1.
               A(I,I) = -Q(I)
               DO K=2,IM1
                  A(I,K) = A(IM1,K-1) - MU(I)*A(IM1,K)                  &       
     &                  - SIG(IM1)*A(I-2,K)
               END DO
               A(I,1) = -MU(I)*A(IM1,1) - SIG(IM1)*A(I-2,1)
               DO K=1,IP1
                  IHF4 = IM1 + K
                  NORM(I) = NORM(I) + A(I,K)*MOMENT(IHF4)
               END DO
               SIG(I) = NORM(I)/NORM(IM1)
               IF(SIG(I).LE.0.)  THEN
                  N = I - 1
                  ITEST = -1
                  GO TO 80
               END IF
            END DO
         END IF
      END IF
      IF(ITEST.GT.0) THEN
         NP1 = N+1
         DO K=1,NP1
            IHF4 = N+K
            L(NP1) = L(NP1) + A(N,K)*MOMENT(IHF4)
         END DO
         Q(NP1) = L(NP1)/NORM(N)
         MU(NP1) = Q(NP1)-Q(N)
      END IF
   80 DO I=1,N
         MEAN(I)=MU(I)
         VAR(I)=SIG(I)
      END DO
   90 IF(ITEST.NE.0)   THEN
         MEAN(N+1)=MU(N+1)
      END IF
!
      RETURN
      END SUBROUTINE GETMUS
!
!***********************************************************************
!
      REAL(KIND=R4)FUNCTION QUE(NDEGRE,XXX)
!
!     CALCULATES VALUES OF THE Q POLYNOMIALS BY RECURRENCE RELATION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NDEGRE
      REAL(KIND=R4) :: XXX
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: A,B,C
!
      A = 1.
      B = XXX - MEAN(1)
      IF(NDEGRE.GT.1) THEN
         DO I=2,NDEGRE
            C = ((XXX - MEAN(I)) * B) - VAR(I-1) * A
            A = B
            B = C
         END DO
      ELSE IF(NDEGRE.EQ.1) THEN
         QUE = B
      ELSE
         QUE = A
      END IF
!
      RETURN
      END FUNCTION QUE
!
!***********************************************************************
!
      SUBROUTINE NEGTIV
!
!     DETERMINES IF LEGENDRE EXPANSION IS NEGATIVE ON
!     THE INTERVAL -1 TO +1 BY STURM SEQUENCE METHODS.
!     IF EXPANSION IS NEGATIVE, ROUTINE FINDS NEGATIVE
!     REGIONS AND INTEGRATES NEGATIVE PROBABILITY.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NPOW,NM1,NNINT,I,NRT,N1,N2,N3,NC
      REAL(KIND=R4) :: PNEG,PROB,F1,F2,F3,FM1,R1,R2,X1,X2,X3,SGN
!
      W = 0.
      PNEG = 0.
      CALL GENSTR(NPOW)
      IF(NPOW.EQ.0)   GO TO 100
      CALL STERM(1.,NPOW,F1)
      IF(F1.EQ.0.)   THEN
         NPOW = NPOW + 1
         F1 = -LEGVAL(NCOEF,+1.)
      END IF
      CALL STERM(-1.,NM1,FM1)
      IF(FM1.EQ.0.)   THEN
         FM1 = LEGVAL(NCOEF,-1.)
      END IF
      IF(NM1.EQ.NPOW)   GO TO 100
      NNINT = NM1 - NPOW + 1
      X1 = 1.
      N1 = NPOW
      I = 1
      NRT = 0
      R1 = X1
      X2 = -1.
      N2 = NM1
      F2 = FM1
      SGN = F1
      NC = 50
   20 IF(N2.NE.N1+1)  THEN
         R2 = FINDRT(X1,X2,F1)
         GO TO 50
      END IF
   30 NC=NC-1
      IF(NC.LE.0)   THEN
         NRT = N2 - N1
         R2 = X1
         GO TO 50
      END IF
      X3 = (X1+X2)/2.
      CALL STERM(X3,N3,F3)
      IF(F3.NE.0.)   THEN
         IF(N3.EQ.N1)   THEN
            X1 = X3
            F1 = F3
            GO TO 30
         END IF
         X2 = X3
         N2 = N3
         F2 = F3
         GO TO 20
      END IF
      IF(N3.EQ.N1)  THEN
         R2 = X3
         X2 = X3
         N2 = N3 + 1
         F2 = -LEGVAL(NCOEF,X3)
         GO TO 50
      END IF
      F3 = LEGVAL(NCOEF,X3)
      X2 = X3
      N2 = N3
      F2 = F3
      GO TO 20
   50 IF(SGN.LT.0.)   THEN
         PROB = PRBINT(R1,R2)
         EMESS = 'DISTRIBUTION IS NEGATIVE'
         CALL ERROR_MESSAGE(NSEQP1)
         WRITE(EMESS,'(4X,3(A,1PE11.4))')                               &       
     &         'FROM ',R1,' TO ',R2,' NEGATIVE PROBABILITY IS ',PROB
         CALL ERROR_MESSAGE(0)
         PNEG = PNEG + PROB
      END IF
      IF(NRT.NE.0)   THEN
         NNINT = NNINT - NRT + 1
         NRT = 0
      END IF
      I = I + 1
      IF(I.LT.NNINT)   THEN
         X1 = X2
         F1 = F2
         N1 = N2
         R1 = R2
         X2 = -1.
         N2 = NM1
         F2 = FM1
         SGN = F1
         NC = 50
         GO TO 20
      ELSE IF(I.EQ.NNINT)  THEN
         R1 = R2
         R2 = -1.
         SGN = FM1
         GO TO 50
      ELSE
         IF(PNEG.NE.0.)   THEN
            WRITE(EMESS,'(A,1PE12.5)')                                  &       
     &           'TOTAL NEGATIVE PROBABILITY = ',PNEG
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE NEGTIV
!
!***********************************************************************
!
      SUBROUTINE GENSTR(NPOW)
!
!     GENERATES LEGENDRE EXPANSION OF STURM POLYNOMIALS
!     FROM ORIGINAL EXPANSION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NPOW
!
      INTEGER(KIND=I4) :: N,N1,N2,N3
      INTEGER(KIND=I4) :: I,L
      REAL(KIND=R4) :: FI,FL,FN2,FN3,A,B
!
      DO I=1,NCOEF
         N = NCOEF + 1 - I
         IF(F(N).NE.0.)   GO TO 20
      END DO
      NPOWS = 0
      GO TO 100
   20 NPOWS = N + 1
      DO I=1,N
         FI = I
         W(I+1,NPOWS) = (FI+.5)*F(I)
      END DO
      W(1,NPOWS) = .5
      CALL LEGDRV(NPOWS,W(1,NPOWS),W(1,N))
      IF(N.GT.1)   THEN
         DO I=2,N
            N3 = NPOWS - I
            N2 = N3 + 1
            N1 = N2 + 1
            FN2 = N2
            IF(W(N2,N2).NE.0.) THEN
               A=(2.*FN2-1.)*W(N1,N1)/FN2/W(N2,N2)
            ELSE
               A = 0.
            END IF
            FN3 = N3
            IF(W(N2,N2).NE.0.) THEN
               B=(W(N2,N1)-FN3*A*W(N3,N2)/(2.*FN2-3.))/W(N2,N2)
            ELSE
               B = 0.
            END IF
            IF(N3.GT.1)   THEN
               DO L=2,N3
                  FL = L
                  W(L,N3)=A*((FL-1.)/(2.*FL-3.)*W(L-1,N2)+FL/(2.*FL+1.)*&       
     &            W(L+1,N2))+B*W(L,N2)-W(L,N1)
               END DO
            END IF
            W(1,N3)=A*W(2,N2)/3.+B*W(1,N2)-W(1,N1)
         END DO
      END IF
!
  100 NPOW = NPOWS
      RETURN
      END SUBROUTINE GENSTR
!
!***********************************************************************
!
      SUBROUTINE LEGDRV(NP1,F,FP)
!
!     GENERATES DERIVATIVE OF LEGENDRE EXPANSION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NP1
      REAL(KIND=R4), DIMENSION(NP1) :: F,FP
!
      INTEGER(KIND=I4) :: N,IP1
      INTEGER(KIND=I4) :: I,J
      REAL(KIND=R4) :: FI
!
      N = NP1 - 1
      DO I=1,N
         FP(I) = 0.
         IP1 = I + 1
         DO J=IP1,NP1,2
            FP(I) = FP(I) + F(J)
         END DO
         FI = I
         FP(I) = (2.*FI-1.)*FP(I)
      END DO
!
      RETURN
      END SUBROUTINE LEGDRV
!
!***********************************************************************
!
      SUBROUTINE STERM(XXX,NCS,FUNC)
!
!     CALCULATES NUMBER OF ROOTS GREATER THAN ARGUMENT
!     FROM STURM POLYNOMIALS.
!
      INTEGER(KIND=I4) :: NCS
      REAL(KIND=R4) :: XXX,FUNC
!
      INTEGER(KIND=I4) :: I,J
      REAL(KIND=R4) :: FF,FI,SSIGN
      REAL(KIND=R4), DIMENSION(65) :: P
!
      P(1) = 1.
      P(2) = XXX
      DO I=3,NPOWS
         FI = I
         P(I) = ((2.*FI-3.)*XXX*P(I-1)-(FI-2.)*P(I-2))/(FI-1.)
      END DO
      NCS = 0
      SSIGN = W(1,1)
      DO I=2,NPOWS
         FF = 0.
         DO J=1,I
            FF = FF + W(J,I)*P(J)
         END DO
         IF(FF*SSIGN.LT.0.) THEN
            NCS = NCS + 1
            SSIGN = FF
         END IF
      END DO
      FUNC = FF
!
      RETURN
      END SUBROUTINE STERM
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION LEGVAL(NP1,XXX)
!
!     CALCULATES VALUE OF LEGENDRE EXPANSION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NP1
      REAL(KIND=R4) :: XXX
!
      INTEGER(KIND=I4) :: I
!
      REAL(KIND=R4) :: FI
      REAL(KIND=R4), DIMENSION(65) :: P
!
      P(1) = 1.
      P(2) = XXX
      LEGVAL = W(1,NP1) + W(2,NP1)*XXX
      DO I=3,NP1
         FI = I
         P(I) = ((2.*FI-3.)*XXX*P(I-1)-(FI-2.)*P(I-2))/(FI-1.)
         LEGVAL = LEGVAL + W(I,NP1)*P(I)
      END DO
!
      RETURN
      END FUNCTION LEGVAL
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION FINDRT(X1,X2,F1)
!
!     FINDS ROOT GIVEN BOUNDS ON ROOT.
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: X1,X2,F1
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: XA,XB,XC,FA,FC
!
      XA = X1
      XB = X2
      FA = F1
      DO I=1,37
         XC = (XA+XB)/2.
         FC = LEGVAL(NPOWS,XC)
         IF(FC*FA.GT.0) THEN
            XA = XC
         ELSE IF(FC*FA.GT.0) THEN
            XB = XC
         ELSE
            GO TO 40
         END IF
      END DO
   40 FINDRT = XC
!
      RETURN
      END FUNCTION FINDRT
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION PRBINT(R1,R2)
!
!     RETURNS -(INTEGRAL OF LEGENDRE EXPANSION BETWEEN ARGUMENTS)
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: R1,R2
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: FI,PROB
      REAL(KIND=R4), DIMENSION(65) :: P
!
      P(1) = 1.
      P(2) = R1
      DO I=3,NPOWS
         FI = I
         P(I) = ((2.*FI-3.)*R1*P(I-1)-(FI-2.)*P(I-2))/(FI-1.)
      END DO
      PROB = W(1,NPOWS)*R1
      DO I=2,NPOWS
         FI = I
         PROB = PROB + W(I,NPOWS)/FI*(R1*P(I)-P(I-1))
      END DO
      P(2) = R2
      DO I=3,NPOWS
         FI = I
         P(I) = ((2.*FI-3.)*R2*P(I-1)-(FI-2.)*P(I-2))/(FI-1.)
      END DO
      PROB = PROB - W(1,NPOWS)*R2
      DO I=2,NPOWS
         FI = I
         PROB = PROB - W(I,NPOWS)/FI*(R2*P(I)-P(I-1))
      END DO
      PRBINT = -PROB
!
      RETURN
      END FUNCTION PRBINT
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
     &        C1H,C2H,L1H,L2H,N1H,N2H,MATP,MFP,MTP,NSEQP
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
      IERX = 0
!
!     READ IN RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90,ERR=95)               &       
     &              C1H,C2H,L1H,L2H,N1H,N2H,MAT,MF,MT,NSEQ
      NSEQP1 = NSEQ
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
      IF(MT.GT.0)   THEN
         IF(MF.GT.0.AND.MAT.NE.0)  I = 1
         GO TO 100
      END IF
!
!     LOOK FOR SEND RECORD
!
      IF(MF.GT.0)   THEN
         IF(MF.GT.0.AND.MAT.NE.0)   I = 2
         GO TO 100
      END IF
!
!     LOOK FOR A FEND OR MEND RECORD
!
      IF(MAT.NE.0)   THEN
         I = 3
      ELSE
         I = 4
      END IF
      GO TO 100
!
!     UNEXPECTED END OF FILE
!
   90 IERX = 2
      GO TO 100
!
!     BAD DATA
!
   95 IERX = 1
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
     &      C12,C22,L12,L22,NR2,NP2,MATP,MFP,MTP,NSEQP
      NSEQP1 = NSEQP
!
!     READ IN INTERPOLATION TABLE
!
      NI = 1
      DO N=1,NR2,3
         NF = NI + 2
         READ(JIN,'(6I11,I4,I2,I3,I5)',END=90)                          &       
     &           (NBT2(NN),JNT2(NN),NN=NI,NF),MATP,MFP,MTP,NSEQ
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
      INTEGER(KIND=I4) :: NI,NF
      INTEGER(KIND=I4) :: N,NN
!
!     READ CONT-LIKE RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90)                      &       
     &     C1,C2,L1,L2,NR,NP,MATP,MFP,MTP,NSEQP
      NSEQP1 = NSEQP
!
!     READ IN INTERPOLATION TABLE
!
      NI = 1
      DO N=1,NR,3
         NF = NI + 2
         READ(JIN,'(6I11,I4,I2,I3,I5)',END=90)                          &       
     &       (NBT(NN),JNT(NN),NN=NI,NF),MATP,MFP,MTP,NSEQP
         NI = NI + 3
      END DO
!
!     READ IN DATA TABLE
!
      NI = 1
   50 NF = NI + 2
      READ(JIN,'(6E11.4,I4,I2,I3,I5)',END=90)                           &       
     &       (X(N),Y(N),N=NI,NF),MATP,MFP,MTP,NSEQP
      NI = NI + 3
      IF(NF.LT.NP)   GO TO 50
      GO TO 100
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
      SUBROUTINE RDLIST
!
!     PROCESS A LIST RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NI
      INTEGER(KIND=I4) :: N
      REAL(KIND=R8), DIMENSION(6) :: Y8

!
!     READ CONT-LIKE RECORD
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)',END=90)                      &       
     &      C1L,C2L,L1L,L2L,NPL,N2L,MATP,MFP,MTP,NSEQP
      NSEQP1 = NSEQP
!
!     READ IN DATA TABLE
!
      NI = 0
   50 READ(JIN,'(6E11.4,I4,I2,I3,I5)',ERR=90,END=90)                    &       
     &        (Y8(N),N=1,6),MATP,MFP,MTP,NSEQP
      IF(MATP.LE.0 .OR. MFP.LE.0 .OR. MTP.LE.0) GO TO 90
      DO N=1,6
        NI=NI+1
        Y(NI)=Y8(N)
        IF(NI.GE.NPL) GO TO 100
      END DO
      GO TO 50
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
      SUBROUTINE RDWRSC(IRW,IUNIT)
!
!     READS OR WRITES A PACKED BINARY TAB1 RECORD ON IUNIT
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IRW,IUNIT
!
      INTEGER(KIND=I4) :: NI,N2
      INTEGER(KIND=I4) :: I,N
!
      INTEGER(KIND=I4), PARAMETER :: ISIZE=250
!
      IF(IRW.EQ.1)   THEN
!
!     WRITE TO SCRATCH
!
         NI = MIN0(NP,ISIZE)
         WRITE(IUNIT) MT,NR,(NBT(I),JNT(I),I=1,NR),NP,                  &       
     &             NI,(X(I),Y(I),I=1,NI)
         IF(NP.GT.NI)   THEN
            DO N=ISIZE+1,NP,ISIZE
               N2 = MIN0(N+ISIZE-1,NP)
               WRITE(IUNIT)  (X(I),Y(I),I=N,N2)
            END DO
         END IF
      ELSE
!
!     READ FROM SCRATCH
!
        READ(IUNIT) MT,NR,(NBT(I),JNT(I),I=1,NR),NP,                    &       
     &          NI,(X(I),Y(I),I=1,NI)
         IF(NP.GT.NI)   THEN
            DO N=ISIZE+1,NP,ISIZE
               N2 = MIN0(N+ISIZE-1,NP)
               READ(IUNIT)  (X(I),Y(I),I=N,N2)
            END DO
         END IF
      END IF
!
      RETURN
      END SUBROUTINE RDWRSC
!
!***********************************************************************
!
      SUBROUTINE STOFP (VARA)
!
!     STORE FLOATING POINT VARIABLE IN THE PACKED RESONANCE
!      PARAMETER ARRAY
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: VARA
!
      RESPAR(IPT) = VARA
      IPT = IPT + 1
!
      RETURN
      END SUBROUTINE STOFP
!
!***********************************************************************
!
      SUBROUTINE RWRES(IPATH)
!
!     SAVE RESONANCE REPRESENTATION FOR CURRENT ISOTOPE ON A SCRATCH
!          FILE FOR LATER USE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPATH
!
      INTEGER(KIND=I4) :: IMAX,N2
      INTEGER(KIND=I4) :: I,LLL,LI,ILI,N
!
      IF(IPATH.EQ.1)  THEN
         WRITE(ISCR) IPT,IRES,EMIDLE,NPED,NRED
         DO LLL=1,10
            IF(NPED(LLL).NE.0) THEN
               WRITE(ISCR) (EP(ILI,LLL),ILI=1,100),                     &       
     &                     (APED(ILI,LLL),ILI=1,100),                   &       
     &                     (NBTED(ILI,LLL),ILI=1,INTABMAX),             &       
     &                     (JNTED(ILI,LLL),ILI=1,INTABMAX)
            END IF
         END DO
         IMAX = IPT - 1
         DO N=1,IMAX,500
            N2 = MIN0(500,IMAX-N+1)
            N2 = N + N2 -1
            WRITE(ISCR) (RESPAR(I),I=N,N2)
         END DO
      ELSE
         READ(ISCR) IPT,IRES,EMIDLE,NPED,NRED
         DO LLL=1,10
            IF(NPED(LLL).NE.0) THEN
               READ(ISCR) (EP(LI,LLL),LI=1,100),                        &       
     &             (APED(LI,LLL),LI=1,100),                             &       
     &             (NBTED(LI,LLL),LI=1,INTABMAX),                       &       
     &             (JNTED(LI,LLL),LI=1,INTABMAX)
            END IF
         END DO
         IMAX = IPT - 1
         DO N=1,IMAX,500
            N2 = MIN0(500,IMAX-N+1)
            N2 = N2 + N -1
            READ(ISCR) (RESPAR(I),I=N,N2)
         END DO
      END IF
!
  100 RETURN
!
      END SUBROUTINE RWRES
!
!***********************************************************************
!
      SUBROUTINE PKTAB1
!
!     PROCESS A TAB1 RECORD AND PUT IN RESONANCE DATA ARRAY
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPT0,NC
      INTEGER(KIND=I4) :: N,NN
      INTEGER(KIND=I4), DIMENSION(6) :: IRESPR
!
!     READ CONT-LIKE RECORD
!
      READ(JIN,'(2E11.4,4I11)')   C1,C2,L1,L2,NR,NP
!
!     READ IN INTERPOLATION TABLE
!
      RESPAR(IPT) = NR
      IPT = IPT + 1
      IPT0 = IPT
      DO N=1,NR,3
         READ(JIN,'(6I11)')  (IRESPR(NN),NN=1,6)
         DO NN=1,6
            IPT = IPT + 1
            RESPAR(IPT) = IRESPR(NN)
         END DO
      END DO
      IPT = IPT0 + 2*NR
!
!     READ IN DATA TABLE
!
      RESPAR(IPT+1) = NP
      IPT = IPT + 1
      NC = (NP+2)/3
      IPT = IPT0
      CALL PAKLIS(JIN,NC)
      IPT = IPT0 + 2*NP
!
      RETURN
      END SUBROUTINE PKTAB1
!
!***********************************************************************
!
      SUBROUTINE PAKLIS(NT,NC)
!
!     STORES A LIST RECORD IN THE PACKED RESONANCE PARAMETER
!      ARRAY
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NT,NC
!
      INTEGER(KIND=I4) :: IPTU
      INTEGER(KIND=I4) :: N,NN
!
      DO N=1,NC
         IPTU = IPT + 5
         READ(NT,'(6E11.5)')  (RESPAR(NN),NN=IPT,IPTU)
         IPT = IPTU + 1
      END DO
!
      RETURN
      END SUBROUTINE PAKLIS
!
!***********************************************************************
!
      SUBROUTINE GMERGE(X1,N1,X2,N2,XF,NF)
!
!     ROUTINE TO MERGE TWO ENERGY GRIDS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N1,N2,NF
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4),DIMENSION(N1) :: X1
      REAL(KIND=R4),DIMENSION(N2) :: X2
      REAL(KIND=R4),DIMENSION(*) :: XF
!
      INTEGER(KIND=I4) :: N1C,N2C
!
!     INITIALIZE
!
      NF = 0
      N1C = 1
      N2C = 1
!
!     COMPARE CURRENT X1 VALUE WITH CURRENT X2 VALUE
!
   20 IF(X1(N1C).LT.X2(N2C))  THEN
!
!     X1 LESS THAN X2 SO STORE IT IN MERGED GRID
!
         NF = NF + 1
         XF(NF) = X1(N1C)
         GO TO 30
!
!     X1 VALUE SAME AS X2 VALUE
!
      ELSE IF(X1(N1C).EQ.X2(N2C))  THEN
         IF(N2C.GE.N2)   GO TO 50
         IF(X2(N2C).EQ.X2(N2C+1))   GO TO 30
         N2C = N2C + 1
!
!      X2 VALUE LESS THAN X1 VALUE SO STORE IT ON THE MERGED GRID
!
      ELSE
         NF = NF + 1
         XF(NF) = X2(N2C)
         IF(N2C.GE.N2)    GO TO 50
         N2C = N2C + 1
      END IF
      GO TO 20
!
!     READY FOR NEXT POINT ON X1 GRID
!
   30 N1C = N1C + 1
      IF(N1C.LE.N1)    GO TO 20
!
!     NO MORE X1 VALUES SO STORE REST OF X2'S
!
      DO N=N2C,N2
         NF = NF + 1
         XF(NF) = X2(N)
      END DO
      GO TO 100
!
!    NO MORE X2 VALUES SO STORE REST OF X1'S
!
   50 DO N=N1C,N1
         NF = NF + 1
         XF(NF) = X1(N)
      END DO
!
  100 RETURN
      END SUBROUTINE GMERGE
!
!***********************************************************************
!
      SUBROUTINE INTER(IP,INTT,XI,YI)
!
!     INTERPOLATION ROUTINE
!
!     XI IS THE INCOMING INDEPENDENT VARIABLE
!     YI IS THE CALCULATED DEPENDENT VARIABLE
!     IP IS THE STARTING GRID POINT FOR XI SEARCH
!     INTT IS THE STARTING POINT IN THE INTERPOLATION
!              REGION DEFINING ARRAY
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IP,INTT
      REAL(KIND=R4) :: XI,YI
!
      INTEGER(KIND=I4) :: I,N
!
      DO N=IP,NP
         IF(XI.LT.X(N))  THEN
            GO TO 20
         ELSE IF(XI.EQ.X(N))  THEN
            IP = N
            YI = Y(IP)
            GO TO 100
         END IF
      END DO
   20 IP = N - 1
      IF(IP.LT.1)   IP = 1
      DO I=INTT,NR
         IF(IP.LT.NBT(I))   GO TO 30
      END DO
   30 IF(X(IP).NE.X(IP+1))  THEN
         INTT = I
         IF(INTT.GT.NR)   INTT = NR
         CALL TERP1(X(IP),Y(IP),X(IP+1),Y(IP+1),XI,YI,JNT(INTT))
      ELSE
         YI = (Y(IP)+Y(IP+1))/2.
      END IF
!
  100 RETURN
      END SUBROUTINE INTER
!
!***********************************************************************
!
      SUBROUTINE TERPR(E,EINT,PKINT,NP,NTERP,INTERP,NR,INP,PK)
!
!     FUNCTION TO RETURN VALUE INTERPOLATED FROM PKINT
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INP
      INTEGER(KIND=I4) :: NP,NR
      INTEGER(KIND=I4), DIMENSION(NR) :: NTERP,INTERP
      REAL(KIND=R4) :: E,PK
      REAL(KIND=R4), DIMENSION(NP) :: EINT,PKINT
!
      REAL(KIND=R4), INTRINSIC :: ABS
!
      INTEGER(KIND=I4) :: J
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: DIFF
      REAL(KIND=R4) :: EOLD
      DATA EOLD/0.0/
!
      PK = 0.0
      DO I=INP,NP
         J = I - 1
         IF(EINT(I).NE.0.) THEN
            DIFF = 1. - E/EINT(I)
            IF(ABS(DIFF).LE.EPI6)   GO TO 50
         ELSE
            IF(E.EQ.EINT(I))  GO TO 50
         END IF
         IF(E.LT.EINT(I))   GO TO 20
      END DO
      GO TO 80
   20 IF(J.EQ.0)   GO TO 80
      INP = J
      DO I=1,NR
         IF(J.LT.NTERP(I))   GO TO 30
      END DO
      I = NR
   30 CALL TERP1(EINT(J),PKINT(J),EINT(J+1),PKINT(J+1),E,PK,INTERP(I))
      GO TO 80
   50 IF(EOLD.NE.E)   THEN
         PK = PKINT(I)
         IF(I.EQ.NP)   GO TO 80
         IF(EINT(I).NE.EINT(I+1))   GO TO 80
         EOLD = E
         GO TO 100
      END IF
      PK = PKINT(I+1)
   80 EOLD = 0.0
!
  100 RETURN
      END SUBROUTINE TERPR
!
!***********************************************************************
!
      SUBROUTINE TERPRP(INOW,E,ANS)
!
!     FUNCTION TO RETURN VALUE INTERPOLATED FROM PACKED RESONANCE DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INOW
      REAL(KIND=R4) :: E,ANS
!
      INTEGER(KIND=I4) :: LR,LP,ITRP
      INTEGER(KIND=I4) :: IIMIN,IIMAX
      INTEGER(KIND=I4) :: J,INP
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: DIFF
      REAL(KIND=R4) :: EOLD
      DATA EOLD/0.0/
!
      ANS = 0.0
      NR = RESPAR(INOW)
      LR = 2*NR + 1
      NP = RESPAR(INOW+LR)
      LP = 2*NP + 1
      IIMIN = INOW + LR + 1
      IIMAX = INOW + LR + LP - 1
      DO I=IIMIN,IIMAX,2
         J = I - 2
         DIFF = 1. - E/RESPAR(I)
         IF(ABS(DIFF).LE.EPI6)   GO TO 50
         IF(E.LT.RESPAR(I))   GO TO 20
      END DO
      GO TO 75
   20 IF(J.EQ.INOW+LR-1)   GO TO 75
      INP = (J + 1 - INOW - LR)/2
      IIMIN = INOW + 1
      IIMAX = INOW + 2*NR
      DO I=IIMIN,IIMAX,2
         IF(INP.LT.RESPAR(I))   GO TO 30
      END DO
      I = INOW + 2*NR - 1
   30 ITRP = RESPAR(I+1)
      CALL TERP1(RESPAR(J),RESPAR(J+1),RESPAR(J+2),RESPAR(J+3),         &       
     &            E,ANS,ITRP)
      GO TO 75
   50 IF(EOLD.NE.E)   THEN
         ANS = RESPAR(I+1)
         IF(I.EQ.INOW-2+LP+LR)   GO TO 75
         IF(RESPAR(I).NE.RESPAR(I+2))   GO TO 75
         EOLD = E
         GO TO 100
      END IF
      ANS = RESPAR(I+3)
   75 EOLD=0.0
  100 INOW = INOW + LR + LP
!
      RETURN
      END SUBROUTINE TERPRP
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
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: XA,YA,XB,YB,XI,YI
!
      IF(I.EQ.1) THEN
         YI = YA
      ELSE IF(I.EQ.2) THEN
         YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
      ELSE IF(I.EQ.3) THEN
         IF(XA.LE.0..OR.XB.LE.0.) THEN
            YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
         ELSE
            YI = YA + ALOG(XI/XA)*(YB-YA)/ALOG(XB/XA)
         END IF
      ELSE IF(I.EQ.4) THEN
         IF(YA.LE.0..OR.YB.LE.0.)   THEN
            YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
         ELSE
            YI = YA*EXP((XI-XA)*ALOG(YB/YA)/(XB-XA))
         END IF
      ELSE IF(I.EQ.5) THEN
         IF(YA.LE.0..OR.YB.LE.0.)   THEN
            IF(XA.LE.0..OR.XB.LE.0.) THEN
               YI = YA + (XI-XA)*(YB-YA)/(XB-XA)
            ELSE
               YI = YA + ALOG(XI/XA)*(YB-YA)/ALOG(XB/XA)
            END IF
         ELSE IF(XA.LE.0..OR.XB.LE.0.)   THEN
               YI = YA*EXP((XI-XA)*ALOG(YB/YA)/(XB-XA))
         ELSE
            IF(XI.LE.0.)   THEN
               YI = YA + ALOG(XI/XA)*(YB-YA)/ALOG(XB/XA)
            ELSE
               YI = YA*EXP(ALOG(XI/XA)*ALOG(YB/YA)/ALOG(XB/XA))
            END IF
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE TERP1
!
!***********************************************************************
!
      SUBROUTINE GRATE(XLP,XHP,ANS)
!
!     INTEGRATE TAB1 FUNCTION===========================================
!     XLP AND XHP ARE THE INTEGRATION LIMITS
!     ANS IS THE ANSWER
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: XLP,XHP,ANS
!
      INTEGER(KIND=I4) :: I,NL,NH,M
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: XL,XH,AN
!
!     SWITCH LIMITS IF NECESSARY AND INITIALIZE
!
      ANS = 0.0
      IF(XLP.EQ.XHP) GO TO 100
      IF(XLP.LT.XHP) THEN
         XL = XLP
         XH = XHP
      ELSE
         XL = XHP
         XH = XLP
      END IF
!
!     LOCATE XL IN TABLE, XL GTHN OR EQUAL TO X(NL)
!
      IF(XL.LT.X(1)) THEN
         XL = X(1)
         IF(XH.LE.XL)   GO TO 100
      END IF
      DO N=1,NP
         NL = N - 1
         IF(XL.LT.X(N))   GO TO 20
      END DO
      GO TO 100
!
!     LOCATE XH IN TABLE, XH GTHN X(NH)
!
   20 IF(XH.GT.X(NP))  THEN
         XH = X(NP)
         NH = NP - 1
      ELSE
         DO N=NL,NP
            NH = N-1
            IF(XH.LE.X(N))  GO TO 30
         END DO
      END IF
!
!     FIND STARTING INTERPOLATION CODE
!
   30 M = 0
   35 M = M + 1
      IF(NL+1.LE.NBT(M))   GO TO 40
      IF(M.LE.NR-1)    GO TO 35
   40 I = JNT(M)
!
!     SUM OVER PANELS
!
      IF(NH.LE.NL)   THEN
!********ONLY ONE PANEL
         CALL RECSI(X(NL),Y(NL),X(NL+1),Y(NL+1),XL,XH,I,ANS)
         GO TO 90
      END IF
!
!*****DO FIRST PANEL
   45 CALL RECSI(X(NL),Y(NL),X(NL+1),Y(NL+1),XL,X(NL+1),I,ANS)
      N = NL
!*****DO INTERMEDIATE PANELS
   50 N = N + 1
   55 IF(N+1.LE.NBT(M))   GO TO 60
      IF(M.GT.NR-1)  THEN
         I = JNT(M)
         GO TO 45
      END IF
      M = M + 1
      GO TO 55
   60 I = JNT(M)
      IF(N.LT.NH)   THEN
         CALL RECSI(X(N),Y(N),X(N+1),Y(N+1),X(N),X(N+1),I,AN)
         ANS = ANS + AN
         GO TO 50
      END IF
!*****DO LAST PANEL
      CALL RECSI(X(N),Y(N),X(N+1),Y(N+1),X(N),XH,I,AN)
      ANS = ANS + AN
!
!     FINISHED
!
   90 IF(XLP.GT.XHP)  ANS = -ANS
!
  100 RETURN
      END SUBROUTINE GRATE
!
!***********************************************************************
!
      SUBROUTINE RECSI(X3,Y3,X4,Y4,X1,X2,I,ANS)
!
!     COMPUTE INTEGRAL OF Y(X)/X==============================
!     Y(X) DEFINED BY THE END POINTS (X3,Y3), (X4,Y4), AND THE
!     INTERPOLATION CODE I.  X1 AND X2 ARE THE INTEGRATION LIMITS.
!     X1 AND X2 MAY LIE OUTSIDE X3 AND X4
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: X3,Y3,X4,Y4,X1,X2,ANS
!
      REAL(KIND=R4), INTRINSIC :: ALOG, EXP, ABS
!
      INTEGER(KIND=I4) :: J,K
      REAL(KIND=R4) :: A,B,X12,AX12,BX,BX1,Z,ZZ,FJ,FK
      REAL(KIND=R4) :: FACK,ARK,ANSP,CON,EXPA
!
      ANS = 0.0
      IF(X4.LE.X3)   GO TO 100
!
!     Y CONSTANT
!
      IF(I.EQ.1) THEN
         ANS = Y3*ALOG(X2/X1)
!
!     Y LINEAR IN X
!
      ELSE IF(I.EQ.2) THEN
         B = (Y4-Y3)/(X4-X3)
         A = Y3 - B*X3
         ANS = A*ALOG(X2/X1) + B*(X2-X1)
!
!     Y LINEAR IN LN(X)
!
      ELSE IF(I.EQ.3) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3 - B*X3
            ANS = A*ALOG(X2/X1) + B*(X2-X1)
         ELSE
            B = (Y4-Y3)/ALOG(X4/X3)
            A = Y3 - B*ALOG(X3)
            X12 = X2/X1
            AX12 = ALOG(X12)
            ANS = A*AX12 + B*AX12*(ALOG(X1) + ALOG(X2))/2.
         END IF
!
!     LN(Y) LINEAR IN X
!
      ELSE IF(I.EQ.4) THEN
         IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3 - B*X3
            ANS = A*ALOG(X2/X1) + B*(X2-X1)
         ELSE
            B = ALOG(Y4/Y3)/(X4-X3)
            A = ALOG(Y3)-B*X3
            EXPA = EXP(A)
            IF(EXPA.LE.1.E-12)  GO TO 100
            X12 = X2/X1
            Z = (X2-X1)/X1
            BX1 = B*X1
            FACK = 1.
            BX = 1.
            BX1 = B*X1
            ARK = 1.
            ANS = EXPA*ALOG(X12)
            DO K=1,25
               FK = K
               BX = BX1*BX
               FACK = FK*FACK
               IF(ABS(Z).GE.0.1)THEN
                  ARK = ARK*X12
                  ANSP = BX*(ARK-1.)/(FK*FACK)
               ELSE
                  ANSP = 1.
                  DO J=1,4
                     FJ = J
                     CON = (FK-FJ)/(FJ+1.)
                     IF(CON.GT.0.)   ANSP = ANSP*(1.0+CON*Z)
                  END DO
                  ANSP = BX*Z*ANSP/FACK
               END IF
               ANSP = ANSP*EXPA
               ANS = ANS + ANSP
              IF(ABS(ANSP/ANS).LT.EPI3)   GO TO 100
            END DO
         END IF
!
!     LN(Y) LINEAR IN LN(X)
!
      ELSE IF(I.EQ.5) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
            IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
               B = (Y4-Y3)/(X4-X3)
               A = Y3 - B*X3
               ANS = A*ALOG(X2/X1) + B*(X2-X1)
            ELSE
               B = ALOG(Y4/Y3)/(X4-X3)
               A = ALOG(Y3)-B*X3
               EXPA = EXP(A)
               IF(EXPA.LE.1.E-12)  GO TO 100
               X12 = X2/X1
               Z = (X2-X1)/X1
               BX1 = B*X1
               FACK = 1.
               BX = 1.
               BX1 = B*X1
               ARK = 1.
               ANS = EXPA*ALOG(X12)
               DO K=1,25
                  FK = K
                  BX = BX1*BX
                  FACK = FK*FACK
                  IF(ABS(Z).GE.0.1)THEN
                     ARK = ARK*X12
                     ANSP = BX*(ARK-1.)/(FK*FACK)
                  ELSE
                     ANSP = 1.
                     DO J=1,4
                        FJ = J
                        CON = (FK-FJ)/(FJ+1.)
                        IF(CON.GT.0.)   ANSP = ANSP*(1.0+CON*Z)
                     END DO
                     ANSP = BX*Z*ANSP/FACK
                  END IF
                  ANSP = ANSP*EXPA
                  ANS = ANS + ANSP
                 IF(ABS(ANSP/ANS).LT.EPI3)   GO TO 100
               END DO
            END IF
         ELSE IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
            IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
               B = (Y4-Y3)/(X4-X3)
               A = Y3 - B*X3
               ANS = A*ALOG(X2/X1) + B*(X2-X1)
            ELSE
               B = (Y4-Y3)/ALOG(X4/X3)
               A = Y3 - B*ALOG(X3)
               X12 = X2/X1
               AX12 = ALOG(X12)
               ANS = A*AX12 + B*AX12*(ALOG(X1) + ALOG(X2))/2.
            END IF
         ELSE
            B = ALOG(Y4/Y3)/ALOG(X4/X3)
            Z = B*ALOG(X2/X1)
            IF(ABS(Z).LE.0.1)   THEN
               ZZ = (1.+Z*(.5+0.16666667*Z))
               ANS = Y3*((X1/X3)**B)*ALOG(X2/X1)*ZZ
            ELSE
               ANS = Y3*((X1/X3)**B)*((X2/X1)**B-1.)/B
            END IF
         END IF
      END IF
!
  100 RETURN
      END SUBROUTINE RECSI
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
      REAL(KIND=R4), INTRINSIC :: ALOG, EXP, ABS
!
      REAL(KIND=R4) :: A, B, Z, ZZ
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
         A = Y3 - B*X3
         ANS = (X2-X1)*(A+0.5*B*(X2+X1))
!
!     Y LINEAR IN LN(X)
!
      ELSE IF(I.EQ.3) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3 - B*X3
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
            A = Y3 - B*X3
            ANS = (X2-X1)*(A+0.5*B*(X2+X1))
         ELSE
            B = ALOG(Y4/Y3)/(X4-X3)
            A = ALOG(Y3) - B*X3
            Z = (X2-X1)*B
            IF(ABS(Z).LE.0.1)   THEN
               ANS = EXP(A+B*X1)*(X2-X1)*(1.0+Z*(0.5+Z*                 &       
     &              (0.1666667+0.04166667*Z)))
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
               A = Y3 - B*X3
               ANS = (X2-X1)*(A+0.5*B*(X2+X1))
            ELSE
               B = ALOG(Y4/Y3)/(X4-X3)
               A = ALOG(Y3) - B*X3
               Z = (X2-X1)*B
               IF(ABS(Z).LE.0.1)   THEN
                  ANS = EXP(A+B*X1)*(X2-X1)*(1.0+Z*(0.5+Z*              &       
     &                 (0.1666667+0.04166667*Z)))
               ELSE
                  ANS = EXP(A+B*X1)*(EXP(Z)-1.0)/B
               END IF
            END IF
         ELSE IF(Y3.LE.0.0.OR.Y4.LE.0.0)   THEN
            IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
               B = (Y4-Y3)/(X4-X3)
               A = Y3 - B*X3
               ANS = (X2-X1)*(A+0.5*B*(X2+X1))
            ELSE
               B = (Y4-Y3)/ALOG(X4/X3)
               Z = (X2-X1)/X1
               IF(ABS(Z).LE.0.1)    THEN
                  ANS = (X2-X1)*(Y3+B*ALOG(X1/X3))+(0.5*B*X1*Z*Z)*      &       
     &                    (1.0+Z*(-0.33333333+Z*(0.16666667-0.1*Z)))
               ELSE
                  ANS = (X2-X1)*(Y3+B*ALOG(X1/X3))+B*X1*                &       
     &                    (1.0+(X2/X1)*(ALOG(X2/X1)-1.0))
               END IF
            END IF
         ELSE
            B=ALOG(Y4/Y3)/ALOG(X4/X3)
            Z=(B+1.0)*ALOG(X2/X1)
            IF(ABS(Z).LE.0.1) THEN
               ZZ = (1.+Z*(0.5+0.16666667*Z))
               ANS=Y3*X1*((X1/X3)**B)*ALOG(X2/X1)*ZZ
            ELSE
               ANS=Y3*X1*((X1/X3)**B)*(((X2/X1)**(B+1.0))-1.0)/(B+1.0)
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
      REAL(KIND=R4), INTRINSIC :: ABS
!
      REAL(KIND=R4) :: A,B,B3,X12,Z
!
      ANS = 0.0
      IF(X4.LE.X3)   GO TO 100
!
!     Y CONSTANT
!
      IF(I.EQ.1) THEN
         ANS = (X2-X1)*Y3*(X2+X1)/2.
!
!     Y LINEAR IN X
!
      ELSE IF(I.EQ.2) THEN
         B = (Y4-Y3)/(X4-X3)
         A = Y3 - B*X3
         B3 = B/3.
         X12 = X1 + X2
         ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
!
!     Y LINEAR IN LN(X)
!
      ELSE IF(I.EQ.3) THEN
         IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3 - B*X3
            B3 = B/3.
            X12 = X1 + X2
            ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
         ELSE
            B = (Y4-Y3)/ALOG(X4/X3)
            A = Y3 - B*ALOG(X3)
            X12 = X2/X1
            Z = (X2-X1)/X1
            IF(ABS(Z).LE.0.10)    THEN
               ANS = (X2-X1)*(A+B*ALOG(X1))*(X2+X1)/2. +                &       
     &               (0.5*B*(X2-X1)*(X2-X1)*(1.+Z*(.3333333-Z*          &       
     &               (0.08333333-0.03333333*Z))))
            ELSE
               ANS = (X2-X1)*(A+B*ALOG(X1))*(X1+X2)/2. +                &       
     &               0.25*B*X1*X1*(1.+X12*X12*(2.0*ALOG(X12)-1.0))
            END IF
         END IF
!
!     LN(Y) LINEAR IN X
!
      ELSE IF(I.EQ.4) THEN
         IF(Y3.LE.0.0.OR.Y4.LE.0.0.OR.Y3.EQ.Y4)   THEN
            B = (Y4-Y3)/(X4-X3)
            A = Y3 - B*X3
            B3 = B/3.
            X12 = X1 + X2
            ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
         ELSE
            B = ALOG(Y4/Y3)/(X4-X3)
            A = ALOG(Y3) - B*X3
            Z = (X2-X1)*B
            IF(ABS(Z).LE.0.1)   THEN
               ANS = EXP(A+B*X1)*(X2-X1)*(1.+(B*X2-1.)*(1.+Z*(.5+Z*     &       
     &               (0.1666667+0.04166667*Z))))/B
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
               A = Y3 - B*X3
               B3 = B/3.
               X12 = X1 + X2
               ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
            ELSE
               B = ALOG(Y4/Y3)/(X4-X3)
               A = ALOG(Y3) - B*X3
               Z = (X2-X1)*B
               IF(ABS(Z).LE.0.1)   THEN
                  ANS = EXP(A+B*X1)*(X2-X1)*(1.+(B*X2-1.)*(1.+Z*(.5+Z*  &       
     &                  (0.1666667+0.04166667*Z))))/B
               ELSE
                  ANS = EXP(A+B*X1)*((B*X2-1.)*EXP(Z)-(B*X1-1.))/(B*B)
               END IF
            END IF
         ELSE IF(Y3.LE.0.0.OR.Y4.LE.0.0.OR.Y3.EQ.Y4)   THEN
            IF(X3.LE.0.0.OR.X4.LE.0.0)   THEN
               B = (Y4-Y3)/(X4-X3)
               A = Y3 - B*X3
               B3 = B/3.
               X12 = X1 + X2
               ANS = (X2-X1)*(X12*(0.5*A+B3*X12) -B3*X1*X2)
            ELSE
               B = (Y4-Y3)/ALOG(X4/X3)
               A = Y3 - B*ALOG(X3)
               X12 = X2/X1
               Z = (X2-X1)/X1
               IF(ABS(Z).LE.0.1)    THEN
                  ANS = (X2-X1)*(A+B*ALOG(X1))*(X2+X1)/2. +             &       
     &                  (0.5*B*(X2-X1)*(X2-X1)*(1.+Z*(.3333333-Z*       &       
     &                  (0.08333333-0.03333333*Z))))
               ELSE
                  ANS = (X2-X1)*(A+B*ALOG(X1))*(X1+X2)/2. +             &       
     &                  0.25*B*X1*X1*(1.+X12*X12*(2.0*ALOG(X12)-1.0))
               END IF
            END IF
         ELSE
            B = ALOG(Y4/Y3)/ALOG(X4/X3)
            Z = (B+2.)*ALOG(X2/X1)
            IF(ABS(Z).LE.0.1)   THEN
               ANS = Y3*X1*X1*((X1/X3)**B)*ALOG(X2/X1)
               ANS = ANS*(1.+Z*(0.5+0.16666667*Z))
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
      SUBROUTINE EAVE(X,Y,NP,NBT,JNT,NR,E,U,THETA,LF,XY,YN,EBAR)
!
!     CALCULATES THE AVERAGE SECONDARY NEUTRON ENERGY
!       FOR A FILE 5 DISTRIBUTION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LF
      REAL(KIND=R4) :: E,U,THETA,XY,YN,EBAR
      INTEGER(KIND=I4) :: NP,NR
      INTEGER(KIND=I4), DIMENSION(NR) :: NBT,JNT
      REAL(KIND=R4), DIMENSION(NP) :: X,Y
!
      REAL(KIND=R4), INTRINSIC :: SQRT
!
      INTEGER(KIND=I4) :: NTERP,INTT,IFIN
      INTEGER(KIND=I4) :: N
      REAL(KIND=R4) :: X3,X4,Y3,Y4,X2,ANS
      REAL(KIND=R4) :: EX,XMAX,SXMAX
      REAL(KIND=R4) :: C1,C2,C3,C4,C5,C6,C7,C8,C9,C6P,CERF,CON
!
      REAL(KIND=R4), PARAMETER :: PI2=0.8863369
!
      EBAR = 0.0
!
!     TABULATED FUNCTION (LF=1)
!
      IF(LF.EQ.1) THEN
         NTERP = 1
         INTT = JNT(1)
         XY = 0.0
         YN = 0.0
         DO N=2,NP
            IF(N.GT.NBT(NTERP))  THEN
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
         IF(YN.GT.0.0)  EBAR = XY/YN
!
!     DISCRETE CHANNEL  (LF=3)
!
      ELSE IF(LF.EQ.3) THEN
         EBAR = E + THETA
!
!     GENERAL EVAPORATION SPECTRUM (LF=5)
!
      ELSE IF(LF.EQ.5) THEN
         NTERP = 1
         INTT = JNT(1)
         XY = 0.0
         YN = 0.0
         IFIN = 0
         XMAX = (E-U)/THETA
         DO N=2,NP
            IF(N.GT.NBT(NTERP)) THEN
               IF(NTERP.LT.NR)   NTERP = NTERP + 1
               INTT = JNT(NTERP)
            END IF
            X3 = X(N-1)
            X4 = X(N)
            Y3 = Y(N-1)
            Y4 = Y(N)
            IF(XMAX.LE.X4)  THEN
               IFIN = 1
               X2 = XMAX
            ELSE
               X2 = X4
            END IF
            CALL XECSI(X3,Y3,X4,Y4,X3,X2,INTT,ANS)
            XY = XY + ANS
            CALL ECSI(X3,Y3,X4,Y4,X3,X2,INTT,ANS)
            YN = YN + ANS
            IF(IFIN.EQ.1)   GO TO 20
         END DO
   20    IF(YN.GT.0.0)   EBAR = THETA*XY/YN
!
!     MAXWELLIAN (LF=7)
!
      ELSE IF(LF.EQ.7) THEN
         XMAX = (E-U)/THETA
         SXMAX = SQRT(XMAX)
         IF(XMAX.LE.80.0)   THEN
            EX = EXP(-XMAX)
         ELSE
            EX = 0.0
         END IF
         YN = (THETA**1.5)*(PI2*ERF(SXMAX)-SXMAX*EX)
         IF(YN.GT.0.0) EBAR = 1.5*THETA - (THETA**2.5)*(SXMAX**3)*EX/YN
!
!     EVAPORATION SPECTRUM  (LF=9)
!
      ELSE IF(LF.EQ.9) THEN
         XMAX = (E-U)/THETA
         IF(XMAX.LE.80.0)  THEN
            EX = EXP(-XMAX)
         ELSE
            EX = 0.0
         END IF
         YN = THETA*THETA*(1.-EX*(1.+XMAX))
         IF(YN.GT.0.0) EBAR = 2.*THETA - (THETA**3)*(XMAX*XMAX)*EX/YN
!
!     WATT SPECTRUM  (LF=11)
!
      ELSE IF(LF.EQ.11) THEN
         C1 = (E-U)/X(1)
         C2 = SQRT(C1)
         IF(C1.LE.25.)  THEN
            C7 = EXP(-C1)
         ELSE
            C7 = 0.
         END IF
         C3 = X(1)*Y(1)/4.
         C4 = SQRT(C3)
         C5 = EXP(C3)
         C6P = (Y(1)*(E-U))
         C6 = SQRT(C6P)
         C8 = HSIN(C6)
         C9 = HCOS(C6)
         CERF = (ERF(C2-C4)+ERF(C2+C4))/2.
         CON = PI2*X(1)*SQRT(X(1)*Y(1))
         YN = CON*C5*CERF  -  X(1)*C7*C8
         IF(YN.GT.0.0)   THEN
            EBAR = CON*X(1)*C5*(3.+2.*C3)*CERF                          &       
     &             - 6.*X(1)*X(1)*C4*C7*(C2*C9-C4*C8)                   &       
     &             - 2.*X(1)*X(1)*C7*((C1+C3)*C8-C6*C9)                 &       
     &             - 2.*X(1)*X(1)*(1.+3.*C3)*C7*C8
            EBAR = EBAR/(2.*YN)
         END IF
      END IF
!
      RETURN
      END SUBROUTINE EAVE
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION ERF(X)
!
!     CALCULATES THE ERROR FUNCTION
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: X
!
      REAL(KIND=R4) :: T
!
      REAL(KIND=R4), PARAMETER :: P=0.3275911
      REAL(KIND=R4), PARAMETER :: A1=0.254829592
      REAL(KIND=R4), PARAMETER :: A2=-0.284496736
      REAL(KIND=R4), PARAMETER :: A3=1.421413741
      REAL(KIND=R4), PARAMETER :: A4=-1.453152027
      REAL(KIND=R4), PARAMETER :: A5=1.06140542
!
      IF(X.GT.5.)   THEN
         ERF = 1.
      ELSE
         T = 1./(1.+P*X)
         ERF = 1. - ((((A5*T+A4)*T+A3)*T+A2)*T+A1)*T*EXP(-X*X)
      END IF
!
      RETURN
      END FUNCTION ERF
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION HSIN(X)
!
!     CALCULATES THE HYPERBOLIC SINE
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: X
!
      IF(X.GE.0.1)   THEN
         HSIN = (EXP(X) - EXP(-X))/2.
      ELSE
         HSIN = X*(1.+X*X*(.3333333+.2*X*X))
      END IF
!
      RETURN
      END FUNCTION HSIN
!
!***********************************************************************
!
      REAL(KIND=R4) FUNCTION HCOS(X)
!
!     CALCULATES THE HYPERBOLIC COSINE
!
      IMPLICIT NONE
!
      REAL(KIND=R4) :: X
!
      HCOS = (EXP(X) + EXP(-X))/2.
!
      RETURN
      END FUNCTION HCOS
!
!***********************************************************************
!
      REAL(KIND=8) FUNCTION AMASS(IZA)
!
!      MASS ESTIMATION SUBROUTINE BASED ON MASS EXCESS OF:
!      NEUTRON  8664967.E-09+-34
!      PROTON   7825037.E-09+-10
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IZA
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      INTEGER(KIND=I4) :: IZ,IA,IN
      INTEGER(KIND=I4) :: I
      REAL(KIND=R4) :: ALPHA,BZSQ,ETANZ
      REAL(KIND=8) :: AMASN,AMASP,ACOR
!
      INTEGER(KIND=I4), PARAMETER :: NMAX=156,ZMAX=100,AMAX=256
      REAL(KIND=R4), PARAMETER :: ETA=195.856
      REAL(KIND=8), PARAMETER :: DMPKV=7289.034D+0
      REAL(KIND=8), PARAMETER :: DMNKV=8071.43D+0
      REAL(KIND=8), PARAMETER :: AMTOKV=931501.6D+0
!
      REAL(KIND=R4), DIMENSION(NMAX) :: G1
      DATA (G1(I),I=1,80)/                                              &       
     &       0.0,       0.0,       0.0,  320879.3,  297139.5,           &       
     &  271128.4,  250761.8,  229720.2,  215067.6,  197882.5,           &       
     &  185297.4,  170782.1,  160681.9,  147878.5,  139703.2,           &       
     &  130243.1,  125186.6,  118324.5,  115717.5,  110858.1,           &       
     &  111135.8,  108806.8,  109664.5,  108054.5,  109624.9,           &       
     &  108976.3,  111058.7,  110826.7,  114462.8,  116659.0,           &       
     &  121583.9,  124502.0,  130293.8,  133854.9,  140506.7,           &       
     &  144468.8,  151617.6,  155940.6,  163488.4,  168271.1,           &       
     &  176153.2,  181094.4,  189060.2,  194039.1,  201905.6,           &       
     &  206880.6,  214480.8,  219448.1,  226977.6,  231998.2,           &       
     &  241028.8,  248368.7,  257234.9,  264288.9,  272992.8,           &       
     &  279791.8,  288668.4,  295289.1,  303934.2,  309989.5,           &       
     &  318346.8,  324005.8,  332008.7,  337377.7,  344931.0,           &       
     &  349998.1,  357237.1,  361754.9,  368528.7,  372556.9,           &       
     &  378813.7,  382315.1,  387995.8,  390907.4,  395852.0,           &       
     &  398110.4,  402260.8,  403839.8,  407186.9,  408025.7/
      DATA (G1(I),I=81,NMAX)/                                           &       
     &  410491.5,  410881.4,  414543.0,  416095.0,  418972.7,           &       
     &  419887.4,  422239.3,  422410.1,  424096.5,  423007.9,           &       
     &  423604.1,  421808.4,  421562.7,  419308.1,  418452.2,           &       
     &  415538.2,  414003.1,  410528.9,  408349.8,  404427.4,           &       
     &  401540.5,  397096.5,  393668.9,  388683.6,  384737.3,           &       
     &  379222.7,  374623.5,  368420.1,  363463.0,  356880.0,           &       
     &  351465.4,  344211.9,  338270.2,  330353.7,  323836.8,           &       
     &  315137.8,  307992.2,  298798.5,  290930.5,  280992.9,           &       
     &  272120.9,  261408.8,  251768.8,  240293.8,  229648.1,           &       
     &  217866.2,  208367.6,  197188.4,  186917.4,  174749.6,           &       
     &  163735.8,  150740.1,  138839.1,  125183.0,  112507.7,           &       
     &   98189.1,   84852.2,   69712.0,   55684.8,   39835.9,           &       
     &   24908.4,    8438.0,   -7186.2,  -24385.6,  -40806.9,           &       
     &  -58689.1,  -75899.1,  -94662.2, -112761.7, -132312.2,           &       
     & -151130.6, -171548.4, -190931.6, -211889.2, -232052.1,           &       
     & -254011.1/
!
      REAL(KIND=R4), DIMENSION(ZMAX) :: G2
      DATA (G2(I),I=1,80)/                                              &       
     &       0.0,  320879.3,  287052.6,  248337.0,  217704.0,           &       
     &  184227.1,  160600.0,  135290.5,  118009.5,   97845.8,           &       
     &   83454.0,   66390.7,   54883.2,   40810.0,   33119.2,           &       
     &   24486.0,   20012.6,   14055.9,   11846.3,    7634.9,           &       
     &    8093.1,    5736.9,    7035.3,    5758.5,    7692.6,           &       
     &    7571.2,   10943.5,   11979.2,   17235.7,   20928.8,           &       
     &   27311.2,   31718.0,   39147.5,   44036.7,   51939.4,           &       
     &   57185.2,   65438.0,   71173.1,   79859.5,   86354.2,           &       
     &   95650.1,  102132.5,  111200.8,  117480.7,  126496.7,           &       
     &  132876.0,  141940.3,  148313.6,  157308.5,  163554.2,           &       
     &  173644.9,  181464.0,  191139.1,  198502.6,  207860.5,           &       
     &  214735.3,  223670.6,  230088.2,  238617.6,  244596.5,           &       
     &  252660.8,  258242.1,  265820.0,  271033.4,  278177.6,           &       
     &  283105.9,  289876.9,  294556.5,  300848.0,  305290.9,           &       
     &  311337.7,  315619.4,  321180.2,  324907.0,  330172.5,           &       
     &  333763.8,  339154.2,  342451.6,  347144.3,  349875.9/
      DATA (G2(I),I=81,ZMAX)/                                           &       
     &  354120.3,  356653.3,  361929.1,  365394.7,  370103.6,           &       
     &  372748.4,  376796.6,  378780.5,  382285.2,  383810.8,           &       
     &  386771.8,  387749.9,  390099.1,  390644.3,  392345.3,           &       
     &  392227.7,  393341.2,  392623.4,  393116.2,  391842.6/
!
      REAL(KIND=R4), DIMENSION(AMAX) :: G3
      DATA (G3(I),I=1,90)/ 5*0.0,                                       &       
     & -576288.8, -536756.0, -499332.7, -462060.9, -427145.4,           &       
     & -392247.7, -359824.4, -328627.7, -300336.3, -271885.0,           &       
     & -244315.0, -217996.4, -193445.9, -169703.9, -148282.3,           &       
     & -126773.5, -107169.8,  -87432.8,  -68791.2,  -50233.9,           &       
     &  -32982.6,  -16188.2,    -601.3,   14813.5,   28840.3,           &       
     &   42508.5,   54871.0,   66960.2,   77818.7,   88504.9,           &       
     &   98039.2,  107649.5,  115927.4,  124196.1,  131384.3,           &       
     &  138616.8,  144802.4,  151075.3,  156083.5,  161340.2,           &       
     &  165473.3,  169811.5,  173306.0,  177150.6,  180451.7,           &       
     &  183940.3,  186774.2,  189629.6,  191707.8,  193774.7,           &       
     &  195338.4,  197001.7,  197909.3,  198998.6,  199475.6,           &       
     &  199978.1,  199742.4,  199796.8,  199078.2,  198615.5,           &       
     &  197483.2,  196619.8,  195104.7,  193807.2,  191910.9,           &       
     &  190281.4,  187967.0,  185956.3,  183135.2,  180607.8,           &       
     &  177554.9,  174792.0,  171687.0,  168726.9,  165302.2,           &       
     &  162190.7,  158509.7,  155280.7,  151429.1,  147784.6,           &       
     &  143534.7,  139438.1,  134965.1,  130629.4,  126048.2/
      DATA (G3(I),I=91,180)/                                            &       
     &  121779.1,  117184.0,  112823.5,  108152.0,  103635.7,           &       
     &   98911.1,   94280.6,   89742.0,   85237.0,   80253.4,           &       
     &   75688.2,   70663.3,   65940.0,   60942.9,   56254.5,           &       
     &   51282.9,   46619.3,   41679.9,   37066.4,   32152.9,           &       
     &   27597.7,   22614.8,   18012.4,   13060.4,    8457.5,           &       
     &    3615.0,    -915.7,   -5707.7,  -10225.1,  -14990.6,           &       
     &  -19433.5,  -24199.3,  -28661.7,  -33220.2,  -37534.3,           &       
     &  -42024.6,  -46104.7,  -50346.7,  -54222.9,  -58176.7,           &       
     &  -61783.7,  -65508.8,  -68939.1,  -72348.5,  -75576.2,           &       
     &  -78774.5,  -81639.1,  -84675.5,  -87357.7,  -90181.6,           &       
     &  -92832.9,  -95914.1,  -97978.3, -100561.5, -102814.2,           &       
     & -105240.1, -107337.4, -109574.0, -111417.2, -113457.0,           &       
     & -115193.1, -117180.2, -118873.7, -120615.7, -122118.6,           &       
     & -123645.3, -124879.3, -126174.5, -127067.1, -128130.4,           &       
     & -128913.8, -129790.7, -130391.4, -130929.2, -131315.8,           &       
     & -131733.1, -131974.3, -132191.1, -132090.4, -132087.6,           &       
     & -131855.2, -131674.0, -131325.4, -130988.8, -130336.0,           &       
     & -129686.8, -128875.7, -128153.3, -127183.9, -126331.1/
      DATA (G3(I),I=181,AMAX)/                                          &       
     & -125246.7, -124173.3, -122986.4, -121711.7, -120257.3,           &       
     & -118962.4, -117442.6, -116021.7, -114250.2, -112874.2,           &       
     & -111148.0, -109555.9, -107670.2, -106041.0, -104013.6,           &       
     & -102207.8, -100173.1,  -98028.1,  -95721.2,  -93344.9,           &       
     &  -90399.3,  -87809.3,  -84917.7,  -82009.3,  -78759.0,           &       
     &  -75597.0,  -72209.2,  -68615.0,  -64705.4,  -60724.2,           &       
     &  -56462.9,  -52182.9,  -47662.8,  -43221.6,  -38510.9,           &       
     &  -33974.6,  -29106.7,  -24404.5,  -19307.3,  -14512.2,           &       
     &   -9344.8,   -4343.3,     993.5,    6171.0,   11626.5,           &       
     &   16942.3,   22594.8,   28073.6,   33894.5,   39552.7,           &       
     &   45567.5,   51484.7,   57691.1,   63881.7,   70409.5,           &       
     &   76926.6,   83667.5,   90444.6,   97487.7,  104633.9,           &       
     &  112019.1,  119530.2,  127267.7,  135060.1,  143090.3,           &       
     &  151158.2,  159506.8,  167931.0,  176515.4,  185252.4,           &       
     &  194212.4,  203235.4,  212516.4,  221924.8,  231567.3,           &       
     &  241401.8/
!
      AMASS = 0.0
      IF(IZA.GT.0)   THEN
         IZ = IZA/1000
         IA = MOD(IZA,1000)
         IN = IA - IZ
         AMASN = IN*AMTOKV
         AMASP = IZ*AMTOKV
!
!        SEARCH TABLE FOR MASSES NOT CALCULATED
!
         DO I=1,NSMASS
            IF(IZAB(I).EQ.IZA)   THEN
               AMASS = AMASP + AMASN + AMASSB(I)
               GO TO 100
            END IF
         END DO
!
!        ADD NEUTRON AND PROTON MASS DEFECT
!
         AMASN = AMASN + IN*DMNKV
         AMASP = AMASP + IZ*DMPKV
         ALPHA = 16000.*IA
         BZSQ = 120.*IZ**2
         ETANZ = ETA*(IN-IZ)**2
         IF(IN.GT.NMAX.OR.IZ.GT.ZMAX.OR.IA.GT.AMAX) THEN
            ACOR = 0.0
         ELSE
            ACOR = G1(IN) + G2(IZ) + G3(IA)
         END IF
         AMASS = AMASP + AMASN - ALPHA + BZSQ + ETANZ + ACOR
      END IF
!
  100 RETURN
      END FUNCTION AMASS
!
!***********************************************************************
!
      SUBROUTINE ERROR_MESSAGE(ISEQ)
!
!     ROUTINE TO OUTPUT ERROR MESSAGE IN STANDARD FORM
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: ISEQ
!
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM
!
      INTEGER(KIND=I4) :: NEMES
!
!     PUT OUT ERROR MESSAGE
!
      IF(ISEQ.NE.0) THEN
         WRITE(NOUT,'(5X,2A,I6)')  EMESS(1:49),'SEQUENCE NUMBER',ISEQ
      ELSE
         IF(EMESS.EQ.' ') THEN
            NEMES = 1
         ELSE
            NEMES = LEN_TRIM(EMESS)
         END IF
         WRITE(NOUT,'(5X,A)')  EMESS(1:NEMES)
      END IF
!
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
      CHARACTER(LEN=3), DIMENSION(12) ::                                &       
     &        MONTHS = (/'Jan','Feb','Mar','Apr','May','Jun','Jul',     &       
     &                   'Aug','Sep','Oct','Nov','Dec'/)
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
!+++MDC+++
!...VMS
!/      INTEGER(KIND=2) ILENP2
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
      END PROGRAM PSYCHE
!...LWI, DVF, MOD
!/      END MODULE PSYCHE
!---MDC---
