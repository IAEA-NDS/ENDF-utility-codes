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
!     Main program for Non-Windows implementation of INTER
!
      PROGRAM INTER
!
      IMPLICIT NONE
!...LWI, DVF, MOD
!/!
!/!     Module implementation of INTER for MODLIB and Windows
!/!
!/      MODULE INTER
!/!
!/      IMPLICIT NONE
!/!
!/      PRIVATE
!/!
!/      PUBLIC :: RUN_INTER
!/      PUBLIC :: INTER_INPUT, INTER_DATA, INTER_SUCCESS
!...LWI, DVF
!/      PUBLIC :: DEFAULT_ERRX,DEFAULT_EZERO,DEFAULT_ELT
!/      PUBLIC :: DEFAULT_EHT,DEFAULT_ELRI,DEFAULT_EHRI
!/      PUBLIC :: DEFAULT_FTEMP,DEFAULT_ELFI,DEFAULT_EHFI
!/      PUBLIC :: DEFAULT_ERRXT,DEFAULT_EZEROT,DEFAULT_ELTT
!/      PUBLIC :: DEFAULT_EHTT,DEFAULT_ELRIT,DEFAULT_EHRIT
!/      PUBLIC :: DEFAULT_FTEMPT,DEFAULT_ELFIT,DEFAULT_EHFIT
!/      PUBLIC :: DEFAULT_E14,DEFAULT_E14T
!---MDC---
!-T Program INTER
!-P Calculate integral constants from cross sections
!-V
!-V         Version 8.09   August 2016      A. Trkov
!-V                        Correct the text defining the g-factor
!-V         Version 8.08   March 2015       A. Trkov
!-V                        Fix printing multiple isomeric states.
!-V         Version 8.07   October 2013     A. Trkov
!-V                        Improve numerical stability of Maxwellian
!-V                        average cross sections at high spectrum
!-V                        temperatures (e.g. Gd-157 at kT=0.1 MeV)
!-V         Version 8.06   March 2013       A. Trkov
!-V                        Make reading MAT range more robust.
!-V         Version 8.05   September 2012   A. Koning
!-V                        Cleanup unused variables.
!-V         Version 8.04   May 2012       A. Trkov
!-V                        Cases were identified where the precision of
!-V                        the energy variable was insufficient and led
!-V                        to differences in the resonance integral,
!-V                        particularly with zero K cross sections.
!-V                        All variables and floating point operations
!-V                        are now declared double-precision.
!-V         Version 8.03   May 2012       A. Trkov
!-V                        The ANS option to read the filename of the
!-V                        input file is no longer supported. The
!-V                        instructions are read from the default input
!-V                        like in all other cases. If needed, an
!-V                        input file can be piped-in as the default
!-V                        input (works on Windows and Unix).
!-V         Version 8.02   March 2012     R. Capote
!-V                        1. Suppress redefinition of thermal temperature
!-V                        2. Definitions of constants added to header
!-V                        3. Chaged printout to get all isotopes consecutively
!-V                           WARNING - print layout changed!
!-V         Version 8.01   February 2012     A. Trkov
!-V                        1. Process MF 10 data
!-V                        2. Minor improvements for more precision
!-V         Version 8.0    October 2009     A. Trkov
!-V                        Reorganize in style with other codes
!-V         Version 7.0    October 2004     C.L. Dunford
!-V                        1. Modified to provide a module for the NEA
!-V                           MODLIB project
!-V                        2. Permit user to supply batch input file
!-V                           name
!-V                        3. Removed fortran line controls in output
!-V         Version 7.01   January 2005     C.L. Dunford
!-V                        1. Set success flag after return from begin
!-V
!-V      Refer all comments and inquiries to
!-V
!-V         National Nuclear Data Center
!-V         Building 197D
!-V         Brookhaven National Laboratory
!-V         P.O. Box 5000
!-V         Upton, NY 11973-5000
!-V         USA
!-V
!-V      Telephone           631-344-2902
!-V      E-mail              NNDC@BNL.GOV
!-V
!-M
!-M INTER - Calculate integral constants from cross sections
!-M ========================================================
!-M 
!-M INTER is a program that reads a cross section file in ENDF format
!-M and calculates the following integral constants:
!-M      1)Cross sections averaged over a Maxwellian spectrum
!-M      2)The thermal  cross section
!-M      3)The g-factor
!-M      4)The resonance integrals
!-M      5)The 14 MeV cross section
!-M      6)Cross section averaged over a fission spectrum
!-M
!-M NOTE: - Materials with resonance parameter data in File 2 must
!-M         have cross sections reconstructed into tabular form
!-M       - All cross section tables must be linearly interpolable.
!-M
!-M Fortran Logical Units Used:
!-M     5  Default (keyboard) input 
!-M     6  Default output (terminal screen)
!-M    20  Input data file, ENDF format 
!-M    21  Message file for program checking results  
!-M 
!-M Input instructions:
!-M In batch mode operation, the user must supply the following control 
!-M information repeated for each input file to be processed.
!-M
!-M Record  Description
!-M      1  Source ENDF filename 
!-M 
!-M In interactive mode operation, the above data are supplied in 
!-M response to the appropriate query; in graphical mode, via a 
!-M dialog box.
!-M 
!***********************************************************************
!
!     To customize this source run SETMDC
!        ANS  -  ANSI standard batch mode version
!        VMS  -  Command mode for VMS operating system
!        WIN  -  Command mode for PC using Digital Visual Fortran
!        UNX  -  Command mode for Unix using Lahey Fortran
!        DVF  -  Graphical mode for PC using Digital Visual Fortran
!        LWI  -  Graphical mode for Unix using Lahey Winteracter
!        MOD  -  Module for the MODLIB project of NEA WPEC
!
!     The "ANS" VERSION MEETS F95 STANDARDS FOR FIXED OR FREE FORMAT
!       SOURCE
!     The "VMS" version will compile with either the Fortran-77 or
!       Fortran-90 VMS compiler
!     The "DVF" version has a Windows graphical interface. It will
!       compile with the Digital Visual Fortran compiler running
!       under Windows
!     The "LWI" version has an X-Windows graphical interface. It will
!       compile with the Lahey Fortran compiler with Winteracter
!       running under Unix
!
!***********************************************************************
!
!
!     INTER Version Number
!
!+++MDC+++
!...VMS, UNX, ANSI, WIN, LWI, DVF
      CHARACTER(LEN=*), PARAMETER :: VERSION = '8.08'
!...MOD
!/      CHARACTER(LEN=*), PARAMETER :: VERSION = '8.08'
!---MDC---
!
!     Define variable precision
!
      INTEGER(KIND=4), PARAMETER :: I4 = SELECTED_INT_KIND(8)
      INTEGER(KIND=4), PARAMETER :: R4 = SELECTED_REAL_KIND(6,37)
      INTEGER(KIND=4), PARAMETER :: R8 = SELECTED_REAL_KIND(15,307)
!
!     Standard Fortran input and output units
!
      INTEGER(KIND=I4) :: NIN
      INTEGER(KIND=I4), PARAMETER :: INPUT0 = 5, IOUT=6
!
!     ENDF file input and checking output Fortran units
!
      INTEGER(KIND=I4), PARAMETER :: JIN=20,JOUT=21
!
!     Final Fortran output unit
!
      INTEGER(KIND=I4) :: NOUT
!
!     IMDC    flag for compiler option
!     TFMT    format for interactive input prompt
!     OSTATUS status parameter for opening new file
!
!+++MDC+++
!...ANS
      INTEGER(KIND=I4), PARAMETER :: IMDC = 0
      CHARACTER(LEN=*), PARAMETER :: TFMT = '(A)'
      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!...VMS
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 1
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = '(/A,$)'
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'NEW'
!...WIN
!/      INTEGER(KIND=I4), PARAMETER :: IMDC = 2
!/      CHARACTER(LEN=*), PARAMETER :: TFMT = '(/A,$)'
!/      CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
!...UNX
!/    INTEGER(KIND=I4), PARAMETER :: IMDC = 3
!/    CHARACTER(LEN=*), PARAMETER :: TFMT = '(/A,$)'
!/    CHARACTER(LEN=*), PARAMETER :: OSTATUS = 'REPLACE'
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
      CHARACTER(LEN=100)    :: INPAR    ! Command line input text
      INTEGER(KIND=I4)      :: ILENP    ! Command line input text and text length
!
      TYPE INTER_INPUT
         CHARACTER(LEN=100) :: INFIL
         CHARACTER(LEN=100) :: OUTFIL
         INTEGER(KIND=I4)   :: MATMIN
         INTEGER(KIND=I4)   :: MATMAX
         REAL(KIND=R8)      :: E14
         REAL(KIND=R8)      :: ERRX
         INTEGER(KIND=I4)   :: ITHER
         REAL(KIND=R8)      :: EZERO
         REAL(KIND=R8)      :: ELT
         REAL(KIND=R8)      :: EHT
         INTEGER(KIND=I4)   :: IRESI
         REAL(KIND=R8)      :: ELRI
         REAL(KIND=R8)      :: EHRI
         INTEGER(KIND=I4)   :: IFISSI
         REAL(KIND=R8)      :: FTEMP
         REAL(KIND=R8)      :: ELFI
         REAL(KIND=R8)      :: EHFI
      END TYPE INTER_INPUT
!
      TYPE(INTER_INPUT) INTER_DATA
!
      INTEGER (KIND=I4) :: IONEPASS       !  Flag to indicate whether multiple input files can be selected
!                                         !  0, YES;  1, NO
!
!     FLAG TO INDICATE SUCCESS OR FAILURE OF STANEF EXECUTION
!
      INTEGER(KIND=I4) :: INTER_SUCCESS
!
!     FLAG TO INDICATE SUCCESS OR FAILURE OF STANEF EXECUTION
!
      REAL(KIND=R8), PARAMETER :: THER=293.16D0 ! DEFINE THERMAL TEMP.
      REAL(KIND=R8), PARAMETER :: ETH=0.0253D0 ! EV EQUIVALENT
      REAL(KIND=R8) :: TZERO
!
      REAL(KIND=R8), PARAMETER :: DEFAULT_EZERO=0.0D0
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_EZEROT='0.0'
      REAL(KIND=R8), PARAMETER :: DEFAULT_ELT=1.0D-05
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_ELTT='1.0E-05'
      REAL(KIND=R8), PARAMETER :: DEFAULT_EHT=10.0D0
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_EHTT='10.0'
      REAL(KIND=R8), PARAMETER :: DEFAULT_ELRI=0.5D0
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_ELRIT='0.5'
      REAL(KIND=R8), PARAMETER :: DEFAULT_EHRI=1.0D+05
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_EHRIT='1.0E+05'
      REAL(KIND=R8), PARAMETER :: DEFAULT_FTEMP=1.35D+06   ! MODIFIED FROM 1.02E
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_FTEMPT='1.35E+06'
      REAL(KIND=R8), PARAMETER :: DEFAULT_ELFI=1.0D+03
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_ELFIT='1.0E+03'
      REAL(KIND=R8), PARAMETER :: DEFAULT_EHFI=20.D+06
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_EHFIT='20.E+06'
      REAL(KIND=R8), PARAMETER :: DEFAULT_ERRX=0.001D0
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_ERRXT='0.001'
      REAL(KIND=R8), PARAMETER :: DEFAULT_E14=14.0D+6
      CHARACTER(LEN=*), PARAMETER :: DEFAULT_E14T='14.0E+6'
!
      INTEGER(KIND=I4), PARAMETER :: NIRMAX=20   ! INTERPOLATION REGIONS
      INTEGER(KIND=I4), DIMENSION(NIRMAX) :: NBT,JNT
      INTEGER(KIND=I4), PARAMETER :: NITMAX=997  ! POINTS IN PAGED TABLE
      REAL(KIND=R8), DIMENSION(NITMAX) :: X,Y
      INTEGER(KIND=I4) :: IBOT,ITOP ! POINT RANGE OF CURRENT PAGED TABLE
      INTEGER(KIND=I4) :: N1,N2
!
      INTEGER(KIND=I4) :: MATO,MFO,MTO ! CURRENT SECTION
!
!     CURRENT HEAD RECORD CONTENTS
!
      REAL(KIND=R8) :: C1H,C2H
      INTEGER(KIND=I4) :: L1H,L2H,N1H,N2H
      INTEGER(KIND=I4) :: MATH,MFH,MTH,NSP
!     FINAL STATE DESIGNATION
      CHARACTER(LEN=2)  TS,FS
!
!     CONTENTS OF CURRENT TEXT RECORD
!
      CHARACTER(LEN=66) :: TEXT
!
      INTEGER(KIND=I4), PARAMETER :: JMAX=12 ! MAX ITERATIONS TO CONVERGE       
      REAL(KIND=R8), DIMENSION(3,2) :: ELO,EHI ! INTEG PANEL LIMITS
!
      INTEGER(KIND=I4) :: LINES,MATPR  ! CURRENT LINE COUNT, MATERIAL NUMBER    
      INTEGER(KIND=I4), PARAMETER :: MAXLIN=60 !MAXIMUM LINES PER PAGE
!
!     CROSS SECTIONS AT THERMAL, 2200 METERS AND 14MEV
!
      REAL(KIND=R8) :: CZERO,C025,C14
!
!     UNNORMALIZED THERMAL CROSS SECTION, RESONANCE INTEGRAL AND
!       FISSION AVERAGED CROSS SECTION AND NORMALIZATION FACTORS
      REAL(KIND=R8) :: TIN1,TIN2,TIN3
      REAL(KIND=R8) :: PNORM1,PNORMX,PNORMY
!
      REAL(KIND=R8), PARAMETER :: RPI2=0.886226925D0   ! SQRT(PI)/2
!
      INTEGER(KIND=I4) :: NSECT, JSECT
!
!***********************************************************************
!
!+++MDC+++
!...VMS, ANS, WIN, UNX
!
      TEXT=' '
      CALL RUN_INTER
!
!     TERMINATE JOB
!
      IF(INTER_SUCCESS.EQ.0) THEN
         WRITE(IOUT,'(/A)') '   '
         STOP '    INTER - Tests completed successfully'
      ELSE
         WRITE(IOUT,'(/A)') '   '
         STOP '    INTER - Tests terminated abnormally!'
      END IF
!---MDC---
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_INTER
!
!     EXECUTES INTER PROCESS
!
      IMPLICIT NONE
!
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=I4), INTRINSIC :: MOD, IFIX, NINT
      REAL(KIND=R8), INTRINSIC :: SQRT
!
      INTEGER(KIND=I4) :: IQUIT,IFOUND
      INTEGER(KIND=I4) :: IZ,IA,LIS0
      INTEGER(KIND=I4) :: IDESC
      INTEGER(KIND=I4) :: MATP,MTP
      INTEGER(KIND=I4) :: MAT,MF,MT,NS
      INTEGER(KIND=I4) :: NNS,JNS
      INTEGER(KIND=I4) :: I,I000
      REAL(KIND=R8) :: TZ
      REAL(KIND=R8) :: Z1,Z2
      REAL(KIND=R8) :: SIGAV,SIGFAV,GFACT,RESINT
!
      INTEGER(KIND=I4), PARAMETER :: MTMAX=12
      CHARACTER(LEN=8), DIMENSION(MTMAX) :: MTDESC
      DATA MTDESC/'Total   ','Elastic ','Inelas  ','n,2n    ',          &       
     &            'n,3n    ','Fission ','n,gamma ','n,p     ',          &       
     &            'n,d     ','n,t     ','n,He3   ','n,alpha '/
      INTEGER(KIND=I4), DIMENSION(MTMAX) :: MTWANT

!
      DATA MTWANT/1,2,4,16,17,18,102,103,104,105,106,107/
      DATA I000 / 1000 /
!
!     OUTPUT PROGRAM IDENTIFICATION
!
      TZ = 0
      INTER_SUCCESS = 0
      IF(IMDC.LT.4) THEN
         WRITE(IOUT,'(/A/)')  ' PROGRAM INTER VERSION '//VERSION
      END IF
!
!     CHECK FOR COMMAND LINE INPUT (VMS ONLY)
!
      LINES = 0
      IONEPASS = 0
      CALL GET_FROM_CLINE
!
!     INITIALIZE FOR RUN
!
   10 CALL BEGIN(IQUIT)
      IF(IQUIT.GT.0)    THEN
         IF(IONEPASS.EQ.1) INTER_SUCCESS = 1
         GO TO 100
      END IF
      IFOUND = 0
!
!     Write header label
!
      WRITE(NOUT,'(/1(3A))')                                            &       
     &          '   Z   A LISO  LFS  MT  Reaction    Sig(2200)   ',     &       
     &          'Sig(Ezero)  Avg-Sigma  G-fact   Res Integ   ',         &       
     &          'Sig(Fiss)    Sig(E14)   MAT',                          &
     &          ' --- ------------- ---  ---------- ----------- -',     &
     &          '---------- ---------- -------  ----------- -',         &
     &          '---------- ----------- ----'
      WRITE(IOUT,'(/1(3A))')                                            &       
     &          '   Z   A LISO  LFS  MT  Reaction    Sig(2200)   ',     &       
     &          'Sig(Ezero)  Avg-Sigma  G-fact   Res Integ   ',         &       
     &          'Sig(Fiss)    Sig(E14)   MAT',                          &
     &          ' --- ------------- ---  ---------- ----------- -',     &
     &          '---------- ---------- -------  ----------- -',         &
     &          '---------- ----------- ----'
      LINES=10
!
!     READ TAPE LABEL
!
      READ(JIN,'(66X,I4,I2,I3,I5)')   MAT,MF,MT,NS
      IF(MF.NE.0.OR.MT.NE.0)  REWIND(UNIT=JIN)
!
!     READ CONTROL RECORD FOR NEXT MATERIAL
!
   20 CALL CONTIN
C...
C...  PRINT *,' '
C...  print *,'Next material',math,mfh,mth,c1h
C... &       ,INTER_DATA%MATMIN,INTER_DATA%MATMAX
C...
      IF(MATH.GE.INTER_DATA%MATMIN) THEN
         IF(MATH.LE.INTER_DATA%MATMAX) GO TO 30
         IF(INTER_DATA%MATMAX.GT.0) GO TO 90
      ELSE
         IF(MATH.LE.-1) GO TO 90
         CALL MEND
         GO TO 20
      END IF
!
!     MATERIAL SHOULD BE PROCESSED
!
   30 IFOUND = 1
      IA = NINT(C1H)
      IZ = IA/I000
      IA = MOD(IA,I000)
      MATP = MATH
      CALL CONTIN
      LIS0 = L2H
      FS = '  '
      IF     (LIS0.EQ.0) THEN
        TS = '  '
      ELSE IF(LIS0.EQ.1) THEN
        TS = ' m'
      ELSE IF(LIS0.EQ.2) THEN
        TS = ' n'
      ELSE
        TS = ' *'
      END IF
C...
C     print *,'      New Z,A',IZ,IA
C...
C...      IF(INTER_DATA%ITHER.NE.0.AND.TZERO.GE.0.) THEN
C...         IF(MFH.EQ.1.AND.MTH.EQ.451) THEN
C...!
C...!           REDEFINE EZERO FROM TEMPERATURE ON THE ENDF FILE
C...!
C...            CALL CONTIN
C...            CALL CONTIN
C...            IF(C1H.EQ.0) THEN
C...               TZ = THER
C...            ELSE
C...               TZ = C1H
C...            END IF
C...!
C...!           IF TEMPERATURE ON THE FILE DIFFERS BY MORE THAN 0.2 K
C...!
C...            IF(ABS(TZ-TZERO).GT.0.2) THEN
C...               TZERO = TZ
C...               INTER_DATA%EZERO = TZERO*ETH/THER
C...               INTER_DATA%EHT = 20*INTER_DATA%EZERO
C...               Z1 = INTER_DATA%ELT/INTER_DATA%EZERO
C...               Z2 = INTER_DATA%EHT/INTER_DATA%EZERO
C...               PNORM1 = (Z1+1.0)*EXP(-Z1)-(Z2+1.0)*EXP(-Z2)
C...               WRITE(NOUT,'(/A,1PE13.5,A)')                             &       
C...     &           'Temperature from ENDF file redefined to',TZERO,' K'
C...               WRITE(NOUT,'(A,1PE12.5,A)')                              &       
C...     &           'Maxwellian Spectrum at Temperature = ',               &       
C...     &           INTER_DATA%EZERO,' (eV)'
C...               WRITE(NOUT,'(A,1PE12.5,A,1PE12.5,A)')                    &       
C...     &           'Integration Limits from     ',INTER_DATA%ELT,         &       
C...     &           ' (eV) to',INTER_DATA%EHT,' (eV)'
C...               WRITE(NOUT,'(A,1PE12.5//)')                              &       
C...     &           'Integral of Spectrum =        ',PNORM1
C...            END IF
C...         END IF
C...      END IF
!
!     POSITION MATERIAL AT FILE 3
!
  35  DO WHILE (MFH.NE.3 .AND. MFH.NE.10)
         IF(MFH.EQ.0) GO TO 20
         CALL FEND
         CALL CONTIN
C...
C        print *,'Reading',math,mfh,mth
C...
      END DO
!
!     TEST IF THIS SECTION SHOULD BE PROCESSED (SEE MTWANT LIST)
!
  40  MTP = MTH
      DO I=1,MTMAX
         IF(MTP.EQ.MTWANT(I)) THEN
            IDESC = I
            MATO = MATH
            MFO = MFH
            MTO = MTH
            NSECT = N2H
            JSECT = 0
!           CALL OUT_STATUS
            GO TO 45
         END IF
      END DO
      CALL SEND
      GO TO 80
!
   45 JSECT = JSECT + 1
      IF(MFH.EQ.10) THEN
        NNS=N1H
      ELSE
        NNS=1
      END IF
      JNS=0
   46 JNS=JNS+1
      CALL PRSEC
C...
C     print *,'read one xs set',JSECT,NSECT
C...
!
!     READING DONE, READ SEND RECORD IF LAST SECTION
!
      IF(JNS.GE.NNS .AND. JSECT.GE.NSECT) CALL SEND
!
!     CHECK IF REQUESTED LIMITS WERE ACTUALLY USED IN THERMAL
!       INTEGRATION
!
      SIGAV = 0
      GFACT = 0
      IF(INTER_DATA%ITHER.NE.0 .AND. TIN1.NE.0) THEN
         IF(INTER_DATA%ELT.EQ.ELO(1,2).AND.                             &       
     &              INTER_DATA%EHT.EQ.EHI(1,2)) THEN
            PNORMX = PNORM1
         ELSE
            Z1 = ELO(1,2)/INTER_DATA%EZERO
            Z2 = EHI(1,2)/INTER_DATA%EZERO
!           Use Taylor series expansion for small arguments
            IF(Z1.LT.1.D-6 .AND. Z2.LT.1.D-6) THEN
              PNORMX = (Z2-Z1)*(Z2+Z1)/2
            ELSE
              PNORMX = (Z1+1.D0)*EXP(-Z1) - (Z2+1.D0)*EXP(-Z2)
           END IF
         END IF
!
!        CALCULATE AVERAGE CROSS SECTION FOR EZERO
!
         SIGAV = TIN1/(PNORMX*RPI2)
!
!        CALCULATE G-FACTOR
!
         IF(ABS(C025).LT.1.D-20) THEN
           GFACT=0
         ELSE
           GFACT=SQRT(INTER_DATA%EZERO/.0253)*SIGAV/C025
         END IF
         IF(GFACT.GT.99.9999D0) GFACT=99.9999D0
      END IF
!
!     CALCULATE RESONANCE INTEGRAL
!
      IF(INTER_DATA%IRESI.NE.0) THEN
         RESINT = TIN2
      ELSE
         RESINT = 0
      END IF
!
!     CALCULATE FISSION AVERAGE CROSS SECTION
!
      IF(INTER_DATA%IFISSI.NE.0 .AND. TIN3.NE.0) THEN
         SIGFAV=TIN3/PNORMY
      ELSE
         SIGFAV = 0
      END IF
!
!     OUTPUT RESULTS FOR SECTION
!
      IF(MATP.NE.MATPR)  THEN
!        -- Message to default output (screen)
         MATPR = MATP
         IF(LINES+7.GT.MAXLIN) THEN
            WRITE(IOUT,'(1X,A)')  CHAR(12)
            LINES = 1
         ELSE
            WRITE(IOUT,'(1X,A)') ' '
            LINES = LINES + 1
         END IF
C...
C        WRITE(IOUT,'(1x,A,I5)') 'Material number (MAT) = ',MATP
C        WRITE(IOUT,'(1x,3A)')                                          &
C    &          ' Z    A  LISO  LFS  MT  Reaction    Sig(2200)   ',     &       
C    &          'Sig(Ezero)  Avg-Sigma  G-fact   Res Integ   ',         &       
C    &          'Sig(Fiss)   Sig(E14)   '
C        LINES = LINES + 2
C...
      END IF
      IF(LINES.LE.1) THEN
        LINES = LINES + 3
        WRITE(IOUT,'(/1(3A))')                                          &       
     &          '   Z   A LISO  LFS  MT  Reaction    Sig(2200)   ',     &       
     &          'Sig(Ezero)  Avg-Sigma  G-fact   Res Integ   ',         &       
     &          'Sig(Fiss)    Sig(E14)   MAT',                          &
     &          ' --- ------------- ---  ---------- ----------- -',     &
     &          '---------- ---------- -------  ----------- -',         &
     &          '---------- ----------- ----'
      ENDIF
C...
   50 WRITE(IOUT,55) IZ,IA,TS,FS,MTP,MTDESC(IDESC),C025,CZERO,          &       
     &               SIGAV,GFACT,RESINT,SIGFAV,C14,MATP
      WRITE(NOUT,55) IZ,IA,TS,FS,MTP,MTDESC(IDESC),C025,CZERO,          &       
     &               SIGAV,GFACT,RESINT,SIGFAV,C14,MATP
   55 FORMAT(I4,I4,A2,6X,A2,I4,2X,A8,2X,2(1PE12.5),1PE11.4,             &       
     &        0PF8.5,1PE13.5,2(1PE12.5),I5)
      LINES = LINES+1
!
!     CHECK IF ALL SUB-SECTIONS WERE PROCESSED
!
      IF(JNS.LT.NNS) GO TO 46
!
!     CHECK IF ALL SECTIONS WERE PROCESSED
!
      IF(JSECT.LT.NSECT) GO TO 45
!
!     READ HEAD RECORD OF NEXT SECTION OR FEND CARD
!
   80 CALL CONTIN
C...
C        print *,'Next',math,mfh,mth
C...
      IF(MATH.EQ.0) GO TO 20
      IF(MFH.EQ. 0) GO TO 80
      IF(MFH.GT.10) THEN
         CALL MEND
         GO TO 20
      END IF
      GO TO 35
!
!     MESSAGE IF NO MATERIALS FOUND
!
   90 IF(IFOUND.EQ.0) THEN
         IF(IMDC.LT.4) THEN
            WRITE(IOUT,'(//5X,A//)')                                    &       
     &         'NO MATERIALS FOUND IN SELECTED RANGE'
         END IF
         WRITE(NOUT,'(//4X,A//)')                                       &       
     &         'NO MATERIALS FOUND IN SELECTED RANGE'
      END IF
!
!     CLOSE FILES
!
      CLOSE(UNIT=JIN)
      IF(NOUT.NE.6)  CLOSE(UNIT=NOUT)
!
!     SEE IF ONE PASS LIMIT SET
!
      IF(IONEPASS.EQ.0) GO TO 10
!
  100 RETURN
!
      END SUBROUTINE RUN_INTER
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
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM,ICHAR
!
      LOGICAL(KIND=I4) :: IEXIST
      CHARACTER(LEN=50) :: MATSIN,TDIN,RIIN,FIIN
      CHARACTER(LEN=80) :: TITL
      CHARACTER(LEN=1) :: IW
      CHARACTER(LEN=4) :: BUF
      CHARACTER(LEN=12) :: BUF1,BUF2
      INTEGER(KIND=I4) :: IC, ITLEN
      REAL(KIND=R8) :: Z1,Z2
      REAL(KIND=R8) :: EX1,EX2,ER1,ER2,W1,W2
      REAL(KIND=R8) :: PNORM2,PNORM3,RNORM3
      REAL(KIND=R8) :: SQRZ1,SQRZ2,SQRW1,SQRW2
!
      TITL='----------------------------------------'
     &   //'                                        '
!
      NOUT = IOUT
   10 IQUIT = 0
!
! **********************************************************************
! *
! *   INPUT CONSISTS OF THE FOLLOWING=
! *   INFIL        INPUT ENDF FILE SPECIFICATION
! *   OUTFIL       OUTPUT LIST FILE SPECIFICATION
! *   MATMIN/MATMAX RANGE OF MATERIALS TO BE PROCESSED
! *   E14          ANY ENERGY AT WHICH A CROSS SECTION IS TO BE
! *                EXTRACTED
! *   ERRX         FRACTIONAL ERROR ACCEPTABLE IN NUMERICAL INTEGRATION
! *                OF AREA BETWEEN DATA POINTS
! *   ITHER  =Y/N  DO/DO NOT CALCULATE MAXWELLIAN AVERAGES
! *   EZERO        TEMPERATURE OF MAXWELLIAN (EV)
! *   ELT/EHT      LOWER AND UPPER ENERGY LIMITS FOR MAXWELLIAN INT (EV)
! *   IRESI  =Y/N  DO/DO NOT CALCULATE RESONANCE INTEGRALS
! *   ELRI/EHRI    LOWER AND UPPER ENERGY LIMITS FOR RESONANCE INT (EV)
! *   IFISSI =Y/N  DO/DO NOT CALCULATE FISSION AVERAGES
! *   FTEMP        TEMPERATURE OF FISSON SPECTRUM (EV)
! *   ELFI/EHFI    LOWER AND UPPER ENERGY LIMITS FOR FISSION INT (EV)
! **********************************************************************
!
!     INITIALIZE FOR STANDARD OPTIONS
!
      IF(IMDC.LT.4) THEN
         INTER_DATA%INFIL = '*'
         INTER_DATA%OUTFIL = '*'
         INTER_DATA%MATMIN = 0
         INTER_DATA%MATMAX = 0
         INTER_DATA%ERRX = DEFAULT_ERRX
         INTER_DATA%ITHER = 1
         INTER_DATA%EZERO = DEFAULT_EZERO
         INTER_DATA%ELT = DEFAULT_ELT
         INTER_DATA%EHT = DEFAULT_EHT
         INTER_DATA%IRESI = 1
         INTER_DATA%ELRI = DEFAULT_ELRI
         INTER_DATA%EHRI = DEFAULT_EHRI
         INTER_DATA%IFISSI = 1
         INTER_DATA%FTEMP = DEFAULT_FTEMP
         INTER_DATA%ELFI = DEFAULT_ELFI
         INTER_DATA%EHFI = DEFAULT_EHFI
         INTER_DATA%E14 = DEFAULT_E14
      END IF
      SELECT CASE (IMDC)
!... Reading just the filename of the input file is no longer supported
!        CASE (0)
!           IW = 'N'
!           IONEPASS = 0
!        CASE(1,2,3)
         CASE(0,1,2,3)
            IF(ILENP.NE.0)  THEN
               CALL TOKEN(INPAR,'%',1,INTER_DATA%INFIL)
               CALL TOKEN(INPAR,'%',2,INTER_DATA%OUTFIL)
               CALL TOKEN(INPAR,'%',3,IW)
               IC = ICHAR(IW)
               IF(IC.GT.96 .AND. IC.LT.123)   IW = CHAR(IC-32)
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
         IF(INTER_DATA%INFIL.EQ.'*') THEN
!...        IF(IMDC.NE.0) THEN
               WRITE(IOUT,FMT=TFMT)                                     &       
     &             ' Input ENDF File Specification        - '
!...        END IF
            READ(NIN,'(A)') INTER_DATA%INFIL
         ELSE
            WRITE(IOUT,'(/2A)') ' Input file - ', TRIM(INTER_DATA%INFIL)
         END IF
      END IF
!
!     SEE IF INPUT INDICATES FILE TERMINATION
!
      IF(INTER_DATA%INFIL.EQ.' '.OR.INTER_DATA%INFIL.EQ.'DONE') THEN
         IQUIT = 1
         GO TO 100
      END IF
!
!     MAKE SURE INPUT FILE EXISTS
!
      INQUIRE(FILE=INTER_DATA%INFIL,EXIST=IEXIST)
      IF(.NOT.IEXIST)  THEN
         IF(IMDC.LT.4) THEN
            WRITE(IOUT,'(/7X,A/)')  'COULD NOT FIND INPUT FILE'
         END IF
         SELECT CASE (IMDC)
!...        CASE (1,2,3)
            CASE (0,1,2,3)
               IF(IONEPASS.EQ.0) GO TO 10
         END SELECT
         IQUIT = 1
         INTER_SUCCESS = 1
         GO TO 100
      END IF
!
!     GET OUTPUT FILE SPECIFICATION
!
      IF(IMDC.LT.4) THEN
         IF(INTER_DATA%OUTFIL.EQ.'*' ) THEN
!...        IF(IMDC.NE.0) THEN
               WRITE(IOUT,FMT=TFMT)                                     &       
     &           ' Output Message File Specification    - '
!...        END IF
            READ(NIN,'(A)') INTER_DATA%OUTFIL
         ELSE
            WRITE(IOUT,'(/2A)') ' Output file - ',                      &       
     &                    TRIM(INTER_DATA%OUTFIL)
         END IF
      END IF
      IF(INTER_DATA%OUTFIL.NE.' ')  THEN
         NOUT = JOUT             ! SETS FORTRAN OUTPUT UNIT IF DISK FILE
      END IF
!
!     CHECK FOR STANDARD OPTIONS
!
      IF(IW.EQ.'*') THEN
!...     IF(IMDC.GE.1.AND.IMDC.LE.3) THEN
         IF(IMDC.LT.4) THEN
   15       WRITE(IOUT,FMT=TFMT)  ' Standard Options (Y(es),N(o),?)?  '
            READ(NIN,'(A)')  IW
            IC = ICHAR(IW)
            IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
            IF(IW.EQ.'?')  THEN
               IW = '*'
               WRITE(IOUT,20)
   20          FORMAT(10X,' STANDARD OPTIONS ARE'/                      &       
     &             10X,'    PROCESS ENTIRE TAPE'/                       &       
     &             10X,'    CALCULATE THERMAL CROSS SECTION'/           &       
     &             10X,'    CALCULATE 14 MEV CROSS SECTION'/            &       
     &             10X,'    CALCULATE MAXWELL AVERAGE CROSS SECTION'/   &       
     &             10X,'    CALCULATE FISSION SPECTRUM AVERAGE'/        &       
     &             10X,'    CALCULATE RESONANCE INTEGRAL'/              &       
     &             10X,'    FRACTIONAL INTEGRATION ERROR = .001       ')
               GO TO 15
            END IF
         END IF
      END IF
!
!     GET USER OPTION CHOICE WHEN NOT STANDARD
!
      IF(IW.EQ.'N'.AND.IMDC.LT.4) THEN
!
!        MATERIAL NUMBER RANGE SELECTION
!
         CALL SELECT_MATS(MATSIN)
!
!        DEFINE SINGLE ENERGY
!
!...     IF(IMDC.EQ.0) THEN
!...        CALL TOKEN(MATSIN,',',6,BUF1)
!...     ELSE
            WRITE(IOUT,FMT=TFMT)  '  Single Energy (eV)---'
            READ(NIN,'(A)') BUF1
!...     END IF
         IF(BUF1.NE.' ')   THEN
            READ(BUF1,'(BN,E12.5)',ERR=25)  INTER_DATA%E14
         END IF
!
!        DEFINE ERROR
!
   25    CONTINUE
!...     IF(IMDC.EQ.0) THEN
!...        CALL TOKEN(MATSIN,',',7,BUF1)
!...     ELSE
            WRITE(IOUT,TFMT) '  Fractional Error--------'
            READ(NIN,'(A)') BUF1
!...     END IF
         IF(BUF1.NE.' ')   THEN
            READ(BUF1,'(BN,E12.5)',ERR=30)   INTER_DATA%ERRX
         END IF
!
!        THERMAL CROSS SECTIONS
!
   30    CONTINUE
!...     IF(IMDC.EQ.0) THEN
!...        CALL TOKEN(MATSIN,',',3,BUF)
!...        IW = BUF(1:1)
!...     ELSE
            WRITE(IOUT,TFMT) '  Maxwellian Average (Y(es),N(o)) ---'
            READ(NIN,'(A)')  IW
!...     END IF
         IC = ICHAR(IW)
         IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
         IF(IW.NE.'Y')    THEN
            INTER_DATA%ITHER = 0
         ELSE
!...        IF(IMDC.EQ.0) THEN
!...           READ(NIN,'(A)') TDIN
!...           CALL TOKEN(TDIN,',',1,BUF1)
!...        ELSE
               WRITE(IOUT,TFMT) '     Spectrum Temperature (eV)--'
               READ(NIN,'(A)') BUF1
!...        END IF
            IF(BUF1.NE.' ') THEN
               READ(BUF1,'(BN,E12.5)',ERR=35)   INTER_DATA%EZERO
            END IF
   35       CONTINUE
!...        IF(IMDC.EQ.0) THEN
!...           CALL TOKEN(TDIN,',',2,BUF1)
!...           CALL TOKEN(TDIN,',',3,BUF2)
!...        ELSE
               WRITE(IOUT,TFMT)                                         &       
     &            '     Integration Limits(ELOW,EHIGH)(eV)---'
               READ(NIN,'(A)') TDIN
               CALL TOKEN(TDIN,',',1,BUF1)
               CALL TOKEN(TDIN,',',2,BUF2)
!...        END IF
            IF(BUF1.NE.' ')  THEN
               READ(BUF1,'(BN,E12.5)',ERR=40) INTER_DATA%ELT
            END IF
   40       IF(BUF2.NE.' ') THEN
               READ(BUF2,'(BN,E12.5)',ERR=45) INTER_DATA%EHT
            END IF
         END IF
!
!        RESONANCE INTEGRALS
!
   45    CONTINUE
!...     IF(IMDC.EQ.0) THEN
!...        CALL TOKEN(MATSIN,',',4,BUF)
!...        IW = BUF(1:1)
!...     ELSE
            WRITE(IOUT,TFMT)                                            &       
     &           '  Resonance Integral (Y(es),N(o))---'
            READ(NIN,'(A)')  IW
!...     END IF
         IC = ICHAR(IW)
         IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
         IF(IW.NE.'Y')    THEN
            INTER_DATA%IRESI = 0
         ELSE
c...        IF(IMDC.NE.0) THEN
               WRITE(IOUT,TFMT)                                         &       
     &               '     Integration Limits(ELOW,EHIGH)(eV)---'
c...        END IF
            READ(NIN,'(A)') RIIN
            CALL TOKEN(RIIN,',',1,BUF1)
            CALL TOKEN(RIIN,',',2,BUF2)
            IF(BUF1.NE.' ')  THEN
               READ(BUF1,'(BN,E12.5)',ERR=50) INTER_DATA%ELRI
            END IF
   50       IF(BUF2.NE.' ') THEN
               READ(BUF2,'(BN,E12.5)',ERR=55)  INTER_DATA%EHRI
            END IF
         END IF
!
!        FISSION SPECTRUM AVERAGE
!
   55    CONTINUE
!...     IF(IMDC.EQ.0) THEN
!...        CALL TOKEN(MATSIN,',',5,BUF)
!...        IW = BUF(1:1)
!...     ELSE
            WRITE(IOUT,TFMT)                                            &       
     &           '  Fission Spectrum Average (Y(es),N(o)) ---'
            READ(NIN,'(A)')  IW
!...     END IF
         IC = ICHAR(IW)
         IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
         IF(IW.NE.'Y')    THEN
            INTER_DATA%IFISSI = 0
         ELSE
!...        IF(IMDC.EQ.0) THEN
!...           READ(NIN,'(A)') FIIN
!...           CALL TOKEN(FIIN,',',1,BUF1)
!...        ELSE
               WRITE(IOUT,TFMT) '     Spectrum Temperature (eV)--'
               READ(NIN,'(A)') BUF1
!...        END IF
            IF(BUF1.NE.' ') THEN
               READ(BUF1,'(BN,E12.5)',ERR=60)   INTER_DATA%FTEMP
            END IF
   60       CONTINUE
!...        IF(IMDC.EQ.0) THEN
!...           CALL TOKEN(FIIN,',',2,BUF1)
!...           CALL TOKEN(FIIN,',',3,BUF2)
!...        ELSE
               WRITE(IOUT,TFMT)                                         &       
     &              '     Integration Limits(ELOW,EHIGH)(eV)---'
               READ(NIN,'(A)') FIIN
               CALL TOKEN(FIIN,',',1,BUF1)
               CALL TOKEN(FIIN,',',2,BUF2)
!...        END IF
            IF(BUF1.NE.' ')  THEN
               READ(BUF1,'(BN,E12.5)',ERR=65) INTER_DATA%ELFI
            END IF
   65       IF(BUF2.NE.' ') THEN
               READ(BUF2,'(BN,E12.5)',ERR=70) INTER_DATA%EHFI
            END IF
         END IF
      END IF
   70 IF(INTER_DATA%EZERO.LE.0.0) THEN
         INTER_DATA%EZERO = ETH
         TZERO = THER
      ELSE
         TZERO = -1
      END IF
!
!     OPEN INPUT AND OUTPUT FILES
!
      OPEN(UNIT=JIN,ACCESS='SEQUENTIAL',STATUS='OLD',                   &       
     &                   FILE=INTER_DATA%INFIL,ACTION='READ')
      IF(NOUT.NE.6) THEN
!+++MDC+++
!...VMS
!/         OPEN(UNIT=NOUT,ACCESS='SEQUENTIAL',STATUS=OSTATUS,           &       
!/     &       FILE=INTER_DATA%OUTFIL,CARRIAGECONTROL='LIST')
!...WIN, DVF, UNX, LWI, ANS, MOD
         OPEN(UNIT=NOUT,ACCESS='SEQUENTIAL',STATUS=OSTATUS,             &       
     &       FILE=INTER_DATA%OUTFIL)
!---MDC---
      END IF
!
!     CALCULATE INTEGRAL OF THERMAL MAXWELLIAN
!
      IF(INTER_DATA%ITHER.NE.0) THEN
         Z1 = INTER_DATA%ELT/INTER_DATA%EZERO
         Z2 = INTER_DATA%EHT/INTER_DATA%EZERO
!        Use Taylor series expansion for small arguments
         IF(Z1.LT.1.D-6 .AND. Z2.LT.1.D-6) THEN
           PNORM1 = (Z2-Z1)*(Z2+Z1)/2
         ELSE
           PNORM1 = (Z1+1.D0)*EXP(-Z1)-(Z2+1.D0)*EXP(-Z2)
         END IF
      END IF
!
!     CALCULATE INTEGRAL OF SECOND WT FUNCTION (1/E)
!
      IF(INTER_DATA%IRESI.NE.0)   THEN
         PNORM2 = LOG(INTER_DATA%EHRI/INTER_DATA%ELRI)
      END IF
!
!     CALCULATE INTEGRAL OF FISSION MAXWELLIAN
!     Ref: R.J. Perry, C.J. Dean, JEFDOC-487, (April 1994)
!
      IF(INTER_DATA%IFISSI.NE.0) THEN
         Z1 = INTER_DATA%ELFI/INTER_DATA%FTEMP
         Z2 = INTER_DATA%EHFI/INTER_DATA%FTEMP
         SQRZ1 = SQRT(Z1)
         SQRZ2 = SQRT(Z2)
         EX1 = EXP(-Z1)
         EX2 = EXP(-Z2)
         ER1 = DERF(SQRZ1,1)
         ER2 = DERF(SQRZ2,1)
         PNORMY = SQRZ1*EX1 - SQRZ2*EX2 + RPI2*(ER2-ER1)
!***** Renormalisation constant
         W1 = 1.0D3/INTER_DATA%FTEMP
         W2 = 2.0D7/INTER_DATA%FTEMP
         SQRW1 = SQRT(W1)
         SQRW2 = SQRT(W2)
         EX1  =EXP(-W1)
         EX2  =EXP(-W2)
         ER1  = DERF(SQRW1,1)
         ER2  = DERF(SQRW2,1)
         RNORM3 = SQRW1*EX1-SQRW2*EX2 + RPI2*(ER2-ER1)
         PNORM3 = PNORMY/RNORM3
      END IF
!
      IF(IMDC.LT.4) WRITE(IOUT,'(A)') ' '
!
!     OUTPUT SELECTED OPTIONS
!
      WRITE(NOUT,'(A)') CHAR(12)
      ITLEN = LEN_TRIM('PROGRAM INTER VERSION '//VERSION)
      WRITE(NOUT,'(29X,2A)') 'PROGRAM INTER VERSION ',VERSION
      WRITE(NOUT,'(29X,A//A)')  TITL(1:ITLEN),                          &       
     &' Selected Integrations of ENDF File 3 and File 10 Cross Sections'
!
!     INITIALIZE FOR OUTPUT HEADING
!
      LINES = MAXLIN + 1
      MATPR = -1
!
!     OUTPUT INTEGRATION SUMMARY
!
         WRITE(NOUT,'(/A/A,1P,E12.5,A)')                                &       
     &     ' Thermal cross section : Sig(2200) = Sig(Eth)'              &
     &    ,' Thermal energy        : Eth= ',0.0253,' (eV)'
!
      IF(INTER_DATA%ITHER.NE.0)    THEN
         WRITE(NOUT,'(/A/A,1P,E12.5,A)')                                &       
     &     ' Ezero cross section   : Sig(Ezero)'                        &
     &    ,' Ezero energy (input)  : E0 = ',INTER_DATA%EZERO,' (eV)'
         WRITE(NOUT,'(/A,A/A/A,1P,E12.5,A/A,E12.5,A,E12.5,A/A,E12.5)')                      &       
     &     ' Maxwellian average    : Avg-Sigma = (2/sqrt(Pi))'          &
     &    ,' Intg[E1:E2] Sig(E) Phi_m(E) dE / Intg[E1:E2] Phi_m(E) dE'  &
     &    ,' Maxwellian spectrum   : Phi_m(E)  = (E/E0^2) exp(-E/E0)'   &
     &    ,' Spectrum Temperature  : E0 = ',INTER_DATA%EZERO,' (eV)'    &
     &    ,' Integration Limits    : E1 = ',INTER_DATA%ELT              &
     &    ,' (eV)  E2 = ',INTER_DATA%EHT,' (eV)'                        &
     &    ,' Integral of Spectrum       = ',PNORM1
         WRITE(NOUT,'(/A)')                                             &       
     &     ' Westcott g-factor     : G-fact = Avg-Sigma / Sig(2200)'
      END IF
      IF(INTER_DATA%IRESI.NE.0)  THEN
         WRITE(NOUT,'(/A/A,1P,E12.5,A,E12.5,A/A,E12.5,A)')              &
     &    ' Resonance Integral    : Res.Integ = Intg[E3:E4] Sig(E)/E dE'&
     &   ,' Integration Limits    : E3 = ',INTER_DATA%ELRI              &
     &   ,' (eV)  E4 = ',INTER_DATA%EHRI,' (eV)'                        &
     &   ,' Integral of Spectrum       = ',PNORM2
      END IF
      IF(INTER_DATA%IFISSI.NE.0) THEN
         WRITE(NOUT,'(/A,A/A/A,1P,E12.5,A/A,E12.5,A,E12.5,A/A,E12.5)')  &
     &   ' Fiss.spect. average   : Sig(Fiss) = '                        &
     &  ,'Intg[E1:E2] Sig(E) Phi_f(E) dE / Intg[E1:E2] Phi_f(E) dE'     &
     &  ,' Fission spectrum      : Phi_f(E)  = sqrt(E/Ef)/Ef exp(-E/E0)'&
     &  ,' Spectrum Temperature  : Ef = ',INTER_DATA%FTEMP,' (eV)'      &
     &  ,' Integration Limits    : E1 = ',INTER_DATA%ELFI               &
     &  ,' (eV)  E2 = ',INTER_DATA%EHFI,' (eV)'                         &
     &  ,' Integral of Spectrum       = ',PNORM3
      END IF
      WRITE(NOUT,'(/A/A,1P,E12.5,A/)')                                  &       
     &   ' E14 cross-section     : Sig(E14)'                            &
     &  ,' Selected Energy       : E14 = ',INTER_DATA%E14,' eV'
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
      CHARACTER(LEN=10) :: BUF
      CHARACTER(LEN=10) :: BUF1,BUF2
!
      INTEGER(KIND=I4) :: IDASH
      INTEGER(KIND=I4) :: LBUF
!
!     GET THE USER INPUT
!
      BUF='          '
      WRITE(IOUT,'(A)') ' '
      WRITE(IOUT,FMT=TFMT) '     Enter Range of MAT Numbers - '
      READ(NIN,'(A)')  MATSIN
!
!     BLANK RESPONSE IS THE SAME AS SELECTING ALL
!
      IF(MATSIN.EQ.' ')  THEN
         INTER_DATA%MATMIN = 0
         INTER_DATA%MATMAX = 0
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
      INTER_DATA%MATMIN = 1
      INTER_DATA%MATMAX = 9999
      READ(BUF1,'(BN,I10)',ERR=20) INTER_DATA%MATMIN
   20 READ(BUF2,'(BN,I10)',ERR=25) INTER_DATA%MATMAX
!
!     SET THE MATERIAL NUMBER LIMITS
!
   25 IF(INTER_DATA%MATMIN.LE.0) THEN
         INTER_DATA%MATMIN = 1
      END IF
      IF(INTER_DATA%MATMAX.LT.INTER_DATA%MATMIN)  THEN
         INTER_DATA%MATMAX = INTER_DATA%MATMIN
      END IF
      IF(INTER_DATA%MATMIN.EQ.1.AND.INTER_DATA%MATMAX.EQ.9999) THEN
         INTER_DATA%MATMIN = 0
         INTER_DATA%MATMAX = 0
      END IF
!
  100 RETURN
      END SUBROUTINE SELECT_MATS
!
!***********************************************************************
!
      SUBROUTINE CONTIN
!
!     READ IN ONE CONTROL RECORD FROM NT
!
      IMPLICIT NONE
!
      READ(JIN,'(2E11.4,4I11,I4,I2,I3,I5)')                             &       
     &         C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NSP
!
      RETURN
      END SUBROUTINE CONTIN
!
!***********************************************************************
!
      SUBROUTINE SEND
!
!     SKIP TO NEXT SEND CARD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MX
!
      MX = 1
      DO WHILE (MX.GT.0)
         READ(JIN,'(72X,I3,5X)') MX
      END DO
!
      RETURN
      END SUBROUTINE SEND
!
!***********************************************************************
!
      SUBROUTINE FEND
!
!     SKIP TO NEXT FEND CARD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MX
!
      MX = 1
      DO WHILE (MX.GT.0)
         READ(JIN,'(70X,I2,8X)') MX
      END DO
!
      RETURN
      END SUBROUTINE FEND
!
!***********************************************************************
!
      SUBROUTINE MEND
!
!     SKIP TO NEXT MEND CARD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MX
!
      MX = 1
      DO WHILE (MX.GT.0)
         READ(JIN,'(66X,I4,10X)') MX
      END DO
!
      RETURN
      END SUBROUTINE MEND
!
!***********************************************************************
!
      SUBROUTINE PRSEC
!
!     PROCESS ONE SECTION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: MIN0
!
      INTEGER(KIND=I4) :: IPMAX,IST,INTX,N
      INTEGER(KIND=I4) :: I
      REAL(KIND=R8) :: AAINT
!
!     NOTE SECTION HEAD RECORD ALREADY READ
!
      IPMAX = NITMAX
!
!     ZERO OUT INTEGRALS
!
      TIN1 = 0
      TIN2 = 0
      TIN3 = 0
!
!     ZERO OUT CROSS SECTIONS
!
      CZERO = 0
      C14 = 0
      C025 = 0
!
!     READ IN SECTION
!
      CALL CONTIN
      IF(MFH.EQ.10) THEN
        IF     (L2H.EQ.0) THEN
          FS = 'g'
        ELSE IF(L2H.EQ.1) THEN
          FS = 'm'
        ELSE IF(L2H.EQ.2) THEN
          FS = 'n'
!       ELSE IF(L2H.EQ.3) THEN
!         FS = ' o'
        ELSE
          FS = '*'
        END IF
      ELSE
        FS = ' '
      END IF
      N1 = N1H
      N2 = N2H
!*****READ INTERPOLATION TABLE
      READ(JIN,'(6I11,14X)')(NBT(I),JNT(I),I=1,N1)
!*****SET PAGE COUNTERS FOR FIRST PAGE
      IBOT = 0
      IST = 1
      ITOP = MIN0(IPMAX-1,N2)
      N2 = N2 - ITOP
!     
!        READ IN NEXT PAGE
!     
   10    READ(JIN,'(6E11.4,14X)')(X(I),Y(I),I=IST,ITOP)
!     
!        FIND 2200 M/S VALUE
!     
         IF(ETH.GE.X(1).AND.ETH.LE.X(ITOP))  THEN
            CALL FIND(ETH,C025,INTX)
         END IF
!     
!        EXTRACT THERMAL POINT
!     
         IF(INTER_DATA%EZERO.GE.X(1).AND.                               &       
     &                   INTER_DATA%EZERO.LE.X(ITOP)) THEN
            CALL FIND(INTER_DATA%EZERO,CZERO,INTX)
         END IF
!     
!        EXTRACT ANY OTHER POINT
!     
         IF(INTER_DATA%E14.GE.X(1).AND.INTER_DATA%E14.LE.X(ITOP)) THEN
            CALL FIND(INTER_DATA%E14,C14,INTX)
         END IF
!     
!        THERMAL MAXWELLIAN INTEGRATION
!     
         IF(INTER_DATA%ITHER.NE.0)   THEN
            IF(X(1).LT.INTER_DATA%EHT.AND.                              &       
     &                     X(ITOP).GT.INTER_DATA%ELT) THEN
               N = 1
               CALL INTEG(1,N,AAINT,INTER_DATA%ELT,INTER_DATA%EHT)
               TIN1 = TIN1 + AAINT
            END IF
         END IF
!     
!        RESONANCE INTEGRAL
!     
         IF(INTER_DATA%IRESI.NE.0) THEN
            IF(X(1).LT.INTER_DATA%EHRI.AND.                             &       
     &                          X(ITOP).GT.INTER_DATA%ELRI) THEN
               N = 2
               CALL INTEG(3,N,AAINT,INTER_DATA%ELRI,INTER_DATA%EHRI)
               TIN2 = TIN2 + AAINT
            END IF
         END IF
!     
!        FISSION MAXWELLIAN INTEGRATION
!     
         IF(INTER_DATA%IFISSI.NE.0)   THEN
            IF(X(1).LT.INTER_DATA%EHFI.AND.                             &       
     &                         X(ITOP).GT.INTER_DATA%ELFI) THEN
               N = 3
               CALL INTEG(2,N,AAINT,INTER_DATA%ELFI,INTER_DATA%EHFI)
               TIN3 = TIN3 + AAINT
            END IF
         END IF
!     
!        SEE IF MORE PAGES TO BE PROCESSED
!     
      IF(N2.GT.0) THEN
!*****SAVE LAST POINT OF PREVIOUS PAGE
         X(1) = X(ITOP)
         Y(1) = Y(ITOP)
!*****RESET PAGE COUNTERS
         IBOT = IBOT + ITOP - 1
         IST = 2
         ITOP = MIN0(IPMAX,N2+1)
         N2 = N2 + 1 - ITOP
         GO TO 10
      END IF
!
      RETURN
      END SUBROUTINE PRSEC
!
!***********************************************************************
!
      SUBROUTINE INTEG(F,N,TINT,EL,EH)
!
!     SET UP N-TH INTEGRATION
!     TEST INTEGRATION LIMITS - DISCONTINUITIES, POINTS TOO FAR APART
!     OR ZERO RANGES
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: N,F
      REAL(KIND=R8) :: TINT,EL,EH
!
      INTEGER(KIND=I4), INTRINSIC :: INT
      REAL(KIND=R8), INTRINSIC ::  LOG10
!
      INTEGER(KIND=I4) :: IFLAG,INTX,INTR,MLAST
      INTEGER(KIND=I4) :: I,J,K
      REAL(KIND=R8) :: AAINT,XP2,STEMP,XBASE,TEST
      INTEGER(KIND=I4), DIMENSION(2) :: INTER
      REAL(KIND=R8), DIMENSION(2) :: XP,YP,PARTS,GOOF
!
!     TEST LOWER LIMIT
!
      IF(EL.GE.X(1)) THEN
         ELO(N,1) = EL
         ELO(N,2) = EL
      ELSE
         ELO(N,1) = X(1)
!        TEST IF FIRST PANEL
         IF(IBOT.LE.0) ELO(N,2)=X(1)
      END IF
!
!     TEST IF UPPER LIMIT INSIDE PRESENT PANEL
!
      IF(EH.LE.X(ITOP))  THEN
         EHI(N,1) = EH
         EHI(N,2) = EH
      ELSE
         EHI(N,1) = X(ITOP)
         IF(N2.LE.0)   EHI(N,2) = X(ITOP)
      END IF
!
!     START INTEGRATION
!
      IFLAG = 0
      TINT = 0.0
      XP(1) = ELO(N,1)
      CALL FIND(XP(1),YP(1),INTX)
!*****SKIP OVER LOWER END OF PANEL
      DO I=1,ITOP
         IF(X(I).GT.ELO(N,1)) GO TO 10
      END DO
      GO TO 100
!
!     LOOP OVER POINTS INSIDE INTEGRATION RANGE
!
   10 DO J=I,ITOP
         IF(X(J).LT.EHI(N,1)) THEN
            XP(2) = X(J)
         ELSE
            XP(2) = EHI(N,1)
            IFLAG = 1
         END IF
   20    IF(XP(2).EQ.XP(1))THEN
            XP2 = XP(2)*1.000001
         ELSE
            XP2 = XP(2)
         END IF
         CALL FIND(XP2,YP(2),INTR)
!
!        TEST IF INTEGRATION SHOULD BE DONE ANALYTICALLY
!
         IF(INTR.LE.2) THEN
            IF(N.EQ.1)  THEN
               STEMP = INTER_DATA%EZERO
            ELSE
               STEMP = INTER_DATA%FTEMP
            END IF
            CALL AANINT(XP(1),YP(1),XP(2),YP(2),INTR,N,STEMP,AAINT)
            TINT = TINT + AAINT
            IF(IFLAG.GT.0) GO TO 100
            XP(1) = XP(2)
            YP(1) = YP(2)
            GO TO 90
         END IF
!
!        TEST IF FUNCTION IS ZERO AT BOTH PTS BEING INTEGRATED
!
         IF(YP(1).EQ.0 .AND. YP(2).EQ.0) THEN
            XP(1)=XP(2)
            YP(1)=YP(2)
            GO TO 90
         END IF
!
!        TEST IF POINTS ARE TOO FAR APART
!
         TEST = XP(2)/XP(1)
         IF(TEST.GT.1.D3) THEN
            XBASE = XP(1)
            MLAST = INT(LOG10(TEST)) - 1
            DO K=1,MLAST
               XP(2) = XBASE*10.D0**K
               CALL FIND(XP(2),YP(2),INTR)
               IF(YP(1).NE.0 .OR. YP(2).NE.0) THEN
                  CALL GREAT1(F,N,AAINT,2,XP,PARTS,GOOF,INTER,          &       
     &                   INTER_DATA%ERRX, NOUT)
                  TINT = TINT + AAINT
               END IF
               XP(1) = XP(2)
               YP(1) = YP(2)
            END DO
            IF(IFLAG.LE.0)   THEN
               XP(2) = X(J)
            ELSE
               XP(2) = EHI(N,1)
               IFLAG = 1
            END IF
            GO TO 20
!
!        NORMAL INTEGRATION-TEST FOR DISCONTINUITIES
!
         ELSE
            IF(TEST.NE.1) THEN
               CALL GREAT1(F,N,AAINT,2,XP,PARTS,GOOF,INTER,             &       
     &               INTER_DATA%ERRX,NOUT)
            END IF
            TINT = TINT + AAINT
            IF(IFLAG.GT.0) GO TO 100
            XP(1) = XP(2)
            YP(1) = YP(2)
         END IF
   90    CONTINUE
      END DO
!
  100 RETURN
      END SUBROUTINE INTEG
!
!***********************************************************************
!
      SUBROUTINE FIND(XP,YP,INTT)
!
!     GIVEN AN XP FIND RETRIEVES YP FROM PAGED TABLE
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INTT
      REAL(KIND=R8) :: XP,YP
!
      INTEGER(KIND=I4) :: IFLAG,IREAL
      INTEGER(KIND=I4) :: I,J
!
      IFLAG = 0
      YP = 0
      INTT = 0
!
!     LOCATE XP
!
      DO I=1,ITOP
         IF(X(I).EQ.XP) THEN
            YP = Y(I)
            IFLAG = 1
            GO TO 10
         ELSE IF(X(I).GT.XP) THEN
            IF(I.EQ.1) THEN
               YP = 0
               INTT = 0
               GO TO 100
            END IF
            GO TO 10
         END IF
      END DO
      GO TO 100
!*****FIND INTERPOLATION LAW
   10 IREAL = I + IBOT
      DO J=1,N1
         IF(NBT(J).GE.IREAL) GO TO 20
      END DO
      J = N1
   20 INTT = JNT(J)
!*****INTERPOLATE TO GET
      IF(IFLAG.EQ.0) CALL TERP1(X(I-1),Y(I-1),X(I),Y(I),XP,YP,INTT)
!
  100 RETURN
      END SUBROUTINE FIND
!
!***********************************************************************
!
      SUBROUTINE AANINT(X1,Y1,X2,Y2,INLAW,ITYPE,CONS1,AAINT)
!
!     PERFORM ANALYTIC INTEGRATIONS
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: INLAW,ITYPE
      REAL(KIND=R8) :: X1,Y1,X2,Y2,CONS1,AAINT
!
      REAL(KIND=R8), INTRINSIC :: SQRT
!
      INTEGER(KIND=I4) :: JM8
      INTEGER(KIND=I4) :: J
      REAL(KIND=R8) :: Z,FACLOG
      REAL(KIND=R8) :: Z1,Z2,STEMP,EX1,EX2,A,B,C,D,T1,T2,T3,T12
      REAL(KIND=R8) :: RPI2,SQRZ1,SQRZ2,DSQZ,ER21
!
      REAL(KIND=R8), DIMENSION(8) :: AM
      DATA AM/.9999964239D0,-.4998741238D0,.3317990258D0,-.2407338084D0,&
     &        .1676540711D0,-.0953293897D0,.0360884937D0,-.0064535442D0/
!
!     INITIALIZE
!
      AAINT = 0
      IF(ABS(X1-X2).LT.1.0D-30)  GO TO 100
!
!     DETERMINE TYPE OF INTEGRATION
!
      SELECT CASE (ITYPE)
!
!     MAXWELLIAN INTEGRATION -SET PARAMETERS
!
         CASE(1)
            STEMP = CONS1
            Z1 = X1/STEMP
            IF(Z1.LT.88.028D0) THEN
               EX1 = DEXP(-Z1)
            ELSE
               EX1 = 0
            END IF
            Z2 = X2/STEMP
            IF(Z2.LT.88.028D0) THEN
               EX2 = DEXP(-Z2)
            ELSE
               EX2 = 0
            END IF
            SELECT CASE (INLAW)
!
!             (Y) IS CONSTANT IN (X)
!
               CASE (1)
                  A = Y1
!                 Use taylor series expansion for small arguments
                  IF(Z1.LT.1.D-6 .AND. Z2.LT.1.D-6) THEN
                    AAINT = A*(Z2-Z1)*(Z2+Z1)/2
                  ELSE
                    AAINT = A*(EX1*(Z1+1)-EX2*(Z2+1))
                  END IF
!
!             (Y) IS LINEAR IN (X)
!
               CASE (2)
                   B = (DBLE(Y2)-DBLE(Y1))/(DBLE(X2)-DBLE(X1))
                   A = DBLE(Y2)-B*DBLE(X2)
                   C = B*STEMP
                   D = 2*C+A
!                  Use Taylor series expansion for small arguments
                   IF(Z1.LT.1.D-6 .AND. Z2.LT.1.D-6) THEN
                     T1 = (Z2-Z1)*(1-(Z2+Z1)/2)
                     T2 =-(Z2-Z1)*(1-(Z2+Z1))
                     T12= (Z2-Z1)*(Z2+Z1)/2
                     T3 =-(Z2-Z1)*(Z2+Z1)
                   ELSE
!                    T1 = EX1-EX2
!                    T2 = Z1*EX1 - Z2*EX2
                     T12= (1+Z1)*EX1 - (1+Z2)*EX2
                     T3 =  Z1*Z1*EX1 -  Z2*Z2*EX2
                   END IF
!                  AAINT = C*T3 + D*T2 + D*T1
                   AAINT = C*T3 + D*T12
          END SELECT
!
!         1/E INTEGRATION
!
          CASE (2)
c...         Z = (X2-X1)/X1
C...         Z =  dble(X2)/dble(X1)
             Z =  X2/X1
             IF(ABS(Z-1).GT.0.1D0 ) THEN
                FACLOG = LOG(Z)
                Z = Z - 1
             ELSE
                FACLOG = AM(8)
                Z = Z - 1
                DO J=1,7
                  JM8 = 8-J
                  FACLOG=FACLOG*Z+AM(JM8)
                END DO
                FACLOG=FACLOG*Z
             END IF
             SELECT CASE (INLAW)
!
!               (Y) IS CONSTANT IN (X)
!
                CASE (1)
                   AAINT = DBLE(Y1)*FACLOG
!
!               (Y) IS LINEAR IN (X)
!
                CASE (2)
                   AAINT = DBLE(Y2) - DBLE(Y1)
     &                   +((DBLE(Y1*X2)-DBLE(Y2*X1))/DBLE(X2-X1))*FACLOG
            END SELECT
!
!        MAXWELLIAN FISSION INTEGRATION -SET PARAMETERS
!
         CASE (3)
            STEMP = CONS1
            Z1 = DBLE(X1)/STEMP
            IF(Z1.LT.88.028D0) THEN
               EX1 = DEXP(-Z1)
            ELSE
               EX1 = 0
            END IF
            Z2 = DBLE(X2)/STEMP
            IF(Z2.LT.88.028D0)THEN
               EX2=DEXP(-Z2)
            ELSE
               EX2=0
            END IF
            SQRZ1=SQRT(Z1)
            SQRZ2=SQRT(Z2)
            DSQZ =SQRZ2-SQRZ1
            RPI2 =0.886226925D0
!
!           CALCULATE THE DIFFERENCE ERF(SQRZ2)-ERF(SQRZ1)
!
            IF(ABS(DSQZ).GT.1.D-2) THEN
               IF(Z1.LT.1) THEN
!
!              USE DIRECTLY THE ERROR FUNCTION FOR SMALL ARGUMENTS
!
                  ER21 = DERF(SQRZ2,1) - DERF(SQRZ1,1)
               ELSE
!
!              USE THE COMPLEMENTARY ERROR FUNCTION FOR LARGE ARGUMENTS
!
                  ER21 = DERF(SQRZ1,-1) - DERF(SQRZ2,-1)
               END IF
            ELSE
!
!              USE TAYLOR SERIES EXPANSION FOR SMALL DIFFERENCES
!
               ER21 = EXP(-Z1)*DSQZ                                     &       
     &            *( 1 - DSQZ*SQRZ1 + DSQZ*DSQZ*(2*Z1-1)/3 )/RPI2
            END IF
            SELECT CASE (INLAW)
!
!              (Y) IS CONSTANT IN (X)
!
               CASE (1)
                  AAINT =DBLE(Y1)*(SQRZ1*EX1-SQRZ2*EX2 + RPI2*ER21)
!
!              (Y) IS LINEAR IN (X)
!
               CASE(2)
                  A = (DBLE(Y2)-DBLE(Y1))/(DBLE(X2)-DBLE(X1))
                  B = DBLE(Y1)-A*DBLE(X1)
                  C = SQRZ1*EX1-SQRZ2*EX2 + RPI2*ER21
                  C = B*C
                  D = (Z1*SQRZ1*EX1-Z2*SQRZ2*EX2)*STEMP
                  T1 = SQRZ1*EX1-SQRZ2*EX2 + RPI2*ER21
                  T1 = 1.5D0*STEMP*T1
                  T2 = A*(D+T1)
                  AAINT = C+T2
              END SELECT
      END SELECT
!
  100 RETURN
      END SUBROUTINE AANINT
!
!***********************************************************************
!
      SUBROUTINE GREAT1(F,N,FINT,NTAB,XTAB,PARTS,GOOF,INTER,ERROR,NOUT)
!
!     CARRY OUT CONVERGENCE INTEGRATION SCHEME USING TRAPAZOIDAL RULE
!     AND DOUBLING THE NUMBER OF REGIONS PER SUBINTERVAL FOR EACH
!     ITERATION. ONLY DOUBLE UP IN THOSE INTERVALS THAT HAVE NOT ALREADY
!     CONVERGED.
!     F     =SINGLE PRECISION FUNCTION TO BE INTEGRATED
!     FINT  =THE RESULTING INTEGRAL
!     NTAB  =NUMBER OF ORDINATES SUPPLIED (THERE ARE N-1 INTERVALS)
!     XTAB  =TABLE OF THE ORDINATE VALUES. RANGE OF INTEGRATION IS
!            FROM XTAB(1) TO XTAB(NTAB)
!     PARTS =ARRAY OF DIMENSION NTAB, EQUAL TO THE PARTIAL INTEGRALS
!            OVER EACH OF THE NTAB-1 INTERVALS
!     INTER =ARRAY OF DIMENSION NTAB, SEPECIFYING THE NUMBER OF
!            OF THE NTAB-1 INTERVAL.
!            SUBINTERVALS IN EACH INTERVAL
!     GOOF  =ARRAY OF DIMENSION NTAB, EQUAL TO THE NORMAL ERROR IN EACH
!     ERROR =ACCEPTABLE NORMAL ERROR
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NOUT,N,F
      REAL(KIND=R8) :: FINT,ERROR
      INTEGER(KIND=I4) :: NTAB
      INTEGER(KIND=I4), DIMENSION(NTAB) :: INTER
      REAL(KIND=R8), DIMENSION(NTAB) :: XTAB,PARTS,GOOF
!
      REAL(KIND=R8), INTRINSIC :: FLOAT
!
      INTEGER(KIND=I4) :: NM1,II
      INTEGER(KIND=I4) :: I,J,K
      REAL(KIND=R8) :: ERRN,TOTAL,TOTAL1,REST,XNOW,DX
!
!     INITIALIZE VALUE OF THE INTEGRAL
!
      FINT = 0
!
!     CALCULATE THE NUMBER OF INTERVALS
!
      NM1 = NTAB - 1
!
!     CALCULATE ALLOWABLE ERROR PER INTERVAL
!
      ERRN = ERROR/FLOAT(NM1)
!
!     CALCULATE INITIAL APPROXIMATION
!
      TOTAL = 0
      DO I=1,NM1
         INTER(I) = 1
         SELECT CASE (F)
            CASE (1)
               PARTS(I) = 0.5D0*(XTAB(I+1)-XTAB(I))*                    &       
     &                    (FMAXW(XTAB(I+1))+FMAXW(XTAB(I)))
            CASE (2)
               PARTS(I) = 0.5D0*(XTAB(I+1)-XTAB(I))*                    &       
     &                    (FMAXW1(XTAB(I+1))+FMAXW1(XTAB(I)))
            CASE (3)
               PARTS(I) = 0.5D0*(XTAB(I+1)-XTAB(I))*                    &       
     &                    (FEM1(XTAB(I+1))+FEM1(XTAB(I)))
         END SELECT
         TOTAL = TOTAL + PARTS(I)
      END DO
      IF(TOTAL.EQ.0)   GO TO 100
!
!     CALCULATE INITIAL ERRORS
!
      DO I=1,NM1
         GOOF(I) = PARTS(I)
      END DO
!
!     LOOP OVER ITERATIONS
!
      DO J=1,JMAX
         TOTAL1 = TOTAL
!
!        LOOP OVER INTERVALS
!
         DO I=1,NM1
!
!           CHECK FOR CONVERGENCE IN THIS INTERVAL
!
            IF(ABS(GOOF(I)/TOTAL).GT.ERRN) THEN
!              CALCULATE DOUBLE INTERVAL
               DX = (XTAB(I+1)-XTAB(I))/FLOAT(INTER(I))
!              DOUBLE NUMBER OF STEPS
               INTER(I) = 2*INTER(I)
!              INITIALIZE CONTRIBUTION TO INTEGRAL
               REST=0
               II = INTER(I)
!              INITIALIZE ORDINATE
               XNOW = XTAB(I) + 0.5D0*DX
!
!              LOOP OVER ORDINATES
!
               DO K=1,II,2
                  SELECT CASE (F)
                     CASE (1)
                        REST = REST + FMAXW(XNOW)
                     CASE (2)
                        REST = REST + FMAXW1(XNOW)
                     CASE (3)
                        REST = REST + FEM1(XNOW)
                  END SELECT
                  XNOW = XNOW + DX
               END DO
!              CALCULATE NEXT PARTIAL INTEGRAL
               REST = 0.5D0*(PARTS(I) + DX*REST)
!              ADD NEW PARTIAL INTEGRAL AND SUBTRACT OLD
               TOTAL = TOTAL + REST - PARTS(I)
!              CALCULATE NEW ERROR AND SET PARTIAL INTEGRAL TO NEW
               GOOF(I) = REST - PARTS(I)
               PARTS(I) = REST
            END IF
         END DO
!
!        CHECK FOR CONVERGENCE
!
         IF(ABS(1-TOTAL1/TOTAL).LE.ERROR) THEN
            FINT = TOTAL
            GO TO 100
         END IF
      END DO
!
!     THE METHOD HAS NOT CONVERGED
!
      FINT = TOTAL
      WRITE(NOUT,'(9X,A,I2,A,1PE12.5,A,1PE12.5)')                       &       
     &     ' INTEGRATION',N,' NOT CONVERGED FROM E=',XTAB(1),           &       
     &     ' TO E=',XTAB(2)
!
  100 RETURN
      END SUBROUTINE GREAT1
!
!***********************************************************************
!
      SUBROUTINE TERP1(X1,Y1,X2,Y2,X,Y,I)
!
!     INTERPOLATE ONE POINT=============================================
!     (X1,Y1) AND (X2,Y2) ARE THE END POINTS OF THE LINE
!     (X,Y) IS THE INTERPOLATED POINT
!     I IS THE INTERPOLATION CODE
!     NOTE- IF A NEGATIVE OR ZERO ARGUMENT OF A LOG IS
!           DETECTED, THE INTERPOLATION CODE IS AUTOMATICALLY
!           CHANGED FROM LOG TO LINEAR
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: I
      REAL(KIND=R8) :: X1,Y1,X2,Y2,X,Y
!
      INTEGER(KIND=I4) :: II
      REAL(KIND=R8) :: XA,YA,XB,YB,XP,YP
!
      XA = X1
      YA = Y1
      XB = X2
      YB = Y2
      XP = X
      YP = 0
!
      IF(XA.NE.XB)  THEN
         IF(I.LE.0.OR.I.GT.5)  THEN
            II = 2
         ELSE
            II = I
         END IF
         SELECT CASE (II)
            CASE (1)
               YP = YA
            CASE (2)
               YP = YA + (XP-XA)*(YB-YA)/(XB-XA)
            CASE (3)
               IF(XA.LE.0 .OR.XB.LE.0)   THEN
                  YP = YA + (XP-XA)*(YB-YA)/(XB-XA)
               ELSE
                  IF(XP.GT.0) THEN
                     YP = YA + LOG(XP/XA)*(YB-YA)/LOG(XB/XA)
                  END IF
               END IF
            CASE (4)
               IF(YA.LE.0 .OR. YB.LE.0)   THEN
                  YP = YA + (XP-XA)*(YB-YA)/(XB-XA)
               ELSE
                  YP = YA*EXP((XP-XA)*LOG(YB/YA)/(XB-XA))
               END IF
            CASE (5)
               IF(YA.LE.0 .OR. YB.LE.0)  THEN
                  IF(XA.LE.0 .OR. XB.LE.0)   THEN
                     YP = YA + (XP-XA)*(YB-YA)/(XB-XA)
                  ELSE
                     IF(XP.GT.0.) THEN
                        YP = YA + LOG(XP/XA)*(YB-YA)/LOG(XB/XA)
                     END IF
                  END IF
               ELSE IF(XA.LE.0 .OR. XB.LE.0)  THEN
                  YP = YA*EXP((XP-XA)*LOG(YB/YA)/(XB-XA))
               ELSE
                  IF(XP.GT.0.)   THEN
                     YP = YA*EXP(LOG(XP/XA)*LOG(YB/YA)/LOG(XB/XA))
                  END IF
               END IF
 
         END SELECT
      END IF
      Y = YP
!
      RETURN
      END SUBROUTINE TERP1
!
!***********************************************************************
!
      REAL(KIND=R8) FUNCTION FMAXW(X)
!
!     COMPUTES CROSS SECTION TIMES WEIGHT (MAXWELLIAN) AT ENERGY = X
!
      IMPLICIT NONE
!
      REAL(KIND=R8) :: X
!
      REAL(KIND=R8) :: Z,WT,Y
      INTEGER(KIND=I4) :: IDUM
!
!     TEST IF X IS WITHIN REQUESTED RANGE
!
      FMAXW = 0
      IF(X.GE.INTER_DATA%ELT.AND.X.LE.INTER_DATA%EHT) THEN
!
!        CALCULATE WEIGHT
!
         Z = X/INTER_DATA%EZERO
         WT = (Z/INTER_DATA%EZERO)*EXP(-Z)
!
!        GET THE CROSS SECTION
!
         CALL FIND(X,Y,IDUM)
         FMAXW=Y*WT
      END IF
!
      END FUNCTION FMAXW
!
!***********************************************************************
!
      REAL(KIND=R8) FUNCTION FMAXW1(X)
!
!     COMPUTES CROSS SECTION TIMES WEIGHT (FISSION) AT ENERGY = X
!
      IMPLICIT NONE
!
      REAL(KIND=R8) :: X
!
      REAL(KIND=R8), INTRINSIC :: SQRT
!
      INTEGER(KIND=I4) :: IDUM
      REAL(KIND=R8) :: Z,WT,Y
!
!     TEST IF XX IS WITHIN REQUESTED RANGE
!
      FMAXW1 = 0
      IF(X.GE.INTER_DATA%ELFI.AND.X.LE.INTER_DATA%EHFI) THEN
!
!        CALCULATE WEIGHT
!
         Z = X/INTER_DATA%FTEMP
         WT = (SQRT(Z)/INTER_DATA%FTEMP)*EXP(-Z)
!
!        RETRIEVE CROSS SECTION
!
         CALL FIND(X,Y,IDUM)
         FMAXW1 = Y*WT
      END IF
!
      END FUNCTION FMAXW1
!
!***********************************************************************
!
      REAL(KIND=R8) FUNCTION FEM1(X)
!
!     COMPUTES CROSS SECTION TIMES WEIGHT (1/E) AT ENERGY = X
!
      IMPLICIT NONE
!
      REAL(KIND=R8) :: X
!
      INTEGER(KIND=I4) :: IDUM
      REAL(KIND=R8) :: WT,Y
!
!     TEST IF X IS WITHIN REQUESTED RANGE
!
      FEM1 = 0
      IF(X.GE.INTER_DATA%ELRI.AND.X.LE.INTER_DATA%EHRI) THEN
!
!        CALCULATE WEIGHT
!
         WT = 1/X
!
!        RETRIEVE CROSS SECTION
!
         CALL FIND(X,Y,IDUM)
         FEM1 = Y*WT
      END IF
!
      END FUNCTION FEM1
!
!***********************************************************************
!
      REAL(KIND=R8) FUNCTION DERF(XX,IPATH)
!
!     CALCULATES ERROR FUNCTION (IPATH=1) OR THE COMPLEMENTARY ERROR
!     FUNCTION (IPATH=-1) WITH ACCURACY Eps<1.5E-7
!     Ref.: M.Abramowitz,I.Stegun: Handbook of Mathematical Functions,
!     Sec.: 7.1.26, Dover Publications, New York, 9th Ed., Dec.1972.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPATH
      REAL(KIND=R8) :: XX
!
      REAL(KIND=R8) :: X,T,E
      REAL(KIND=R8), PARAMETER :: P=0.3275911D0
      REAL(KIND=R8), PARAMETER :: A1=0.254829592D0
      REAL(KIND=R8), PARAMETER :: A2=-0.284496736D0
      REAL(KIND=R8), PARAMETER :: A3=1.421413741D0
      REAL(KIND=R8), PARAMETER :: A4=-1.453152027D0
      REAL(KIND=R8), PARAMETER :: A5=1.061405429D0
!
      X = ABS(XX)
      T = 1/(1+P*X)
      IF(IPATH.GT.0) THEN
         E = 1 - ((((A5*T+A4)*T+A3)*T+A2)*T+A1)*T*EXP(-X*X)
      ELSE
         E = ((((A5*T+A4)*T+A3)*T+A2)*T+A1)*T*EXP(-X*X)
      END IF
!
      IF(XX.LT.0) THEN
         DERF = -E
      ELSE
         DERF = E
      END IF
!
      RETURN
      END FUNCTION DERF
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
      SUBROUTINE GET_FROM_CLINE
!
!     GET CONTENTS OF COMMAND LINE FOR VMS OR NAME OF BATCH INPUT FILE
!
      IMPLICIT NONE
!
!+++MDC+++
!...VMS
!/      INTEGER(KIND=2) ILENP2
!...ANS
      CHARACTER(LEN=100) :: CFILE
!---MDC---
!
      CFILE = ' '
      INPAR = ' '
      ILENP = 0
      NIN = INPUT0
!+++MDC+++
!...VMS
!/      CALL LIB$GET_FOREIGN(INPAR,,ILENP2)
!/      ILENP = ILENP2
!...UNX
!/     CALL GETCL(INPAR)
!/      ilenp = LEN_TRIM(INPAR)
!...WIN
!/      CALL GETARG(1,INPAR)
!/      ilenp = LEN_TRIM(INPAR)
!---MDC---
!...ANS coding no longer supported
!/      WRITE(IOUT,'(A)')                                                 &       
!/     &    ' Control File Specification        - '
!/      READ(NIN,'(A)') CFILE
!/      NIN = 19
!/      OPEN(UNIT=NIN,FILE=CFILE,STATUS='OLD')
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
      IF(MATO.GT.0.AND.MFO.GT.0.AND.MTO.GT.0) THEN
!+++MDC+++
!...VMS, ANS, WIN, UNX, MOD
         WRITE(IOUT,'(5X,A,I5,A,I3,A,I4)')                              &       
     &         'PROCESSING MAT=',MATO,', MF=',MFO,', MT=',MTO
!---MDC---
      END IF
!
      RETURN
      END SUBROUTINE OUT_STATUS
!
!+++MDC+++
!...VMS, ANS, WIN, UNX
      END PROGRAM INTER
!...LWI, DVF, MOD
!/      END MODULE INTER
!---MDC---
