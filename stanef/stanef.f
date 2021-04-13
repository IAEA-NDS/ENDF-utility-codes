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
! **********************************************************************
! *
!+++MDC+++
!...VMS, ANS, WIN, UNX
!
!     MAIN PROGRAM FOR NON-WINDOWS IMPLEMENTATION OF STANEF
!
      PROGRAM STANEF
!
      IMPLICIT NONE
!...LWI, DVF, MOD
!/!
!/!     MODULE IMPLEMENTATION OF STANEF FOR MODLIB AND WINDOWS
!/!
!/      MODULE STANEF
!/!
!/      IMPLICIT NONE
!/!
!/      PRIVATE
!/!
!/      PUBLIC :: RUN_STANEF
!/      PUBLIC :: STANEF_INPUT, STANEF_DATA, STANEF_SUCCESS
!/!
!---MDC---
!-T Program STANEF
!-P Convert an ENDF file into standard form
!-V         Version 8.10   October 2017 D. Brown
!-V                        - Add checks of P(nu) for fission
!-V                        - Add checks of fission energy release tables
!-V         Version 8.09   October 2012   A. Trkov
!-V                        Allow E-dependent scattering radius in URR
!-V         Version 8.08   September 2012   A. Koning
!-V                        Cleanup unused variables.
!-V         Version 8.07   May 2012, A. Trkov
!-V                        1. Increase POINTSMAX from 500000 to 800000.
!-V         Version 8.06   November 2011, A. Trkov
!-V                        1. Prevent stripping of the header record if
!-V                           tape ID equals zero.
!-V                        2. Special case: suppress redefining ZASYM 
!-V                           if ZA=MAT (special purpose files where
!-V                           the ZA is not define from Z and A)
!-V         Version 8.05   July 2011     A. Trkov
!-V                        Increase No. of points in MF3 from 250K to 500K.
!-V         Version 8.04   January 2011     A. Trkov
!-V                        1. Read data in double-precision to avoid underflow
!-V                        2. Fix processing of MF32, LCOMP=2, LRF=7
!-V         Version 8.03   January 2011     A. Trkov
!-V                        Implement resolved resonance option LRF=7
!-V         Version 8.02   January 2010     A. Trkov
!-V                        Allow long LIST records (ASCII only).
!-V         Version 8.01   November 2008     A. Trkov
!-V                        1. If NSEQ>99999, reset NSEQ=0
!-V         Version 8.00   August  2008     A. Trkov
!-V                        1. Fix minor bugs.
!-V                        2. Implement extended features of the format
!-V                           endorsed by CSEWG.
!-V         Version 7.02   FEBRUARY 2006     M.HERMAN & A.TRKOV  
!-V                        1. SEVERAL DIMENSIONS INCREASED               
!-V                        2. EXTENDED MF=32 PROCESSING
!-V         Version 7.01   February 2005     C.L.Dunford
!-V                        1. Set success flag when done
!-V                        2. Add new element Roentgenium (Rg)
!-V                        3. Corrected symbol generation for second
!-V                            third metastable state (Kellett, NEA)
!-V                        4. added symbol nn for neutron in symbol
!-V                            generation and xx for unnamed elements
!-V                        5. Corrections for fact that some compilers
!-V                            do not recognize the intrinsic function
!-V                            JFIX
!-V         Version 7.0    October 2004     C.L.Dunford
!-V                        1. Modified to provide a module for the NEA
!-V                           MODLIB module
!-V                        2. allow energy dependent delayed fission
!-V                           group parameters.
!-V                        3. Construct proper ZSA symbol for first
!-V                           field on record 5 of 1-451
!-V                        4. Permit user to supply batch input file
!-V                           name
!-V                        5. Added command line input to Unix and
!-V                           windows versions. NOTE: only input and
!-V                           output file names can be given. Default
!-V                           options are assumed unless third
!-V                           parameter is N.
!-V
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
!-M STANEF - Convert ENDF-files to standard form
!-M ============================================
!-M 
!-M STANEF is designed to perform bookkeeping operations on a data 
!-M file containing one or more material evaluations in ENDF format. 
!-M These operations include:
!-M  1.  Creation or modification of a ``tape ID'' record,
!-M  2.  Creation or update of the directory in MT=451,
!-M  3.  Create or modify special hollerith ID records in MT=451 
!-M      (ENDF-6 only),
!-M  4.  Resequencing,
!-M  5.  Conversion of integer and floating point fields to standard 
!-M      format,
!-M  6.  Creation of a binary (ENDF alternate format) file.
!-M 
!-M Fortran Logical Units Used:
!-M     5  Default (keyboard) input 
!-M     6  Default output (terminal screen)
!-M    20  Input data file, ENDF format 
!-M    21  Output ENDF file  
!-M    22  Temporary paging files for large data arrays 
!-M 
!-M Input Requirements:
!-M In batch mode operation, the user must supply the following 
!-M control information repeated for each input file to be processed.
!-M
!-M Record  Description
!-M      1  Input file specification 
!-M      2  Output file specification 
!-M      3  Tape label 
!-M         - Tape label number (integer up to 4 digits)
!-M            <0  means tape has no label
!-M            =0  means copy old label (default)
!-M            >0  means the value is used as the tape label number 
!-M         - Up to 66 character text for ENDF tape label 
!-M      4  Flag to select standard or customised options
!-M             Y  Yes, adopt standard options (default)
!-M             N  No, proceed with the selection of special options
!-M      5  Options selection (4 fields) 
!-M         - Mode of output file (integer) 
!-M             0  Character format output  (default)
!-M             1  Binary format (fields 2 and 3 ignored) 
!-M         - Update directory and special hollerith (character) 
!-M             Y execute this option (default)
!-M             N do not execute this option
!-M         - Numeric field standardization (character)
!-M             Y execute this option  (default)
!-M             N do not execute this option
!-M         - Numeric field standardization for 2/151(character)
!-M             Y  execute this option 
!-M             N do not execute this option	
!-M         If record 5 is left entirely blank, then the default
!-M         options are executed. Those are character format output 
!-M         file, update MT=451 and standardize numeric fields.
!-M 
!-M Multiple input files can be processed to produce multiple output 
!-M files by repeating the above input data sequence. The program 
!-M execution is terminated with a record containing the word DONE.
!-M 
!-M In interactive mode operation, the above data is supplied in 
!-M response to the appropriate query; in graphical mode, via a 
!-M dialog box.
!***********************************************************************
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
!       SOURCE AND USES A FILE FOR NON-INTERACTIVE INPUT (BATCH MODE)
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
!     STANEF VERSION NUMBER
!
!+++MDC+++
!...VMS, UNX, ANSI, WIN, LWI, DVF
      CHARACTER(LEN=*), PARAMETER :: VERSION = '8.09'
!...MOD
!/      CHARACTER(LEN=*), PARAMETER :: VERSION = '8.09'
!---MDC---
!
!     DEFINE VARIABLE PRECISION
!
      INTEGER (KIND= 4), PARAMETER :: I4 = SELECTED_INT_KIND(8)
      INTEGER (KIND= 4), PARAMETER :: R4 = SELECTED_REAL_KIND(6,37)
      INTEGER (KIND= 4), PARAMETER :: R8 = SELECTED_REAL_KIND(15,307)
!
!     STANDARD FORTRAN INPUT AND OUTPUT UNITS
!
      INTEGER(KIND=I4) ::  INPUT
      INTEGER(KIND=I4), PARAMETER :: INPUT0 = 5
      INTEGER(KIND=I4), PARAMETER :: OUTPUT=6
!
!     ENDF DISK FILE INPUT AND CHECKING OUTPUT FORTRAN UNITS
!
      INTEGER(KIND=I4), PARAMETER :: ITAPE=20,OTAPE=21
!
!     SCRATCH UNIT FOR CONSTRUCTING THE DIRECTORY
!
      INTEGER(KIND=I4), PARAMETER :: IDIR=22
!
!     IMDC  FLAG FOR IMPLEMENTATION OPTION
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
!     INPUT DEFINED DATA TYPE
!
      TYPE STANEF_INPUT
         CHARACTER(LEN=100) :: INFIL
         CHARACTER(LEN=100) :: OUTFIL
         INTEGER(KIND=I4) :: MODE
         INTEGER(KIND=I4) :: LABEL
         CHARACTER(LEN=66) :: LTEXT
         CHARACTER(LEN=1) :: INDX
         CHARACTER(LEN=1) :: IDLTES
         CHARACTER(LEN=1) :: I151
      END TYPE STANEF_INPUT
!
      TYPE(STANEF_INPUT) STANEF_DATA
!
!     FLAG TO INDICATE WHETHER MULTIPLE INPUT FILES CAN BE SELECTED
!
      INTEGER(KIND=I4) :: IONEPASS        !  0, YES;  1, NO
!
!     FLAG TO INDICATE SUCCESS OR FAILURE OF STANEF EXECUTION
!
      INTEGER(KIND=I4) :: STANEF_SUCCESS
!
!     FORMAT OF MATERIAL BEING PROCESSED
!
      INTEGER(KIND=I4) :: NFOR
!
!     CONTENTS OF FIELDS ON A HEAD/CONT RECORD
!
      INTEGER(KIND=I4) :: L1H,L2H,N1H,N2H
      REAL(KIND=R8)    :: C1H,C2H
!
!     TEXT PORTION OF A MT=451 RECORD
!
      CHARACTER(LEN=66) :: TEXT
!
      INTEGER(KIND=I4), PARAMETER :: INTERPMAX=400
      INTEGER(KIND=I4), DIMENSION(INTERPMAX) ::  NBT,INTRP
!original start
!     INTEGER(KIND=I4), PARAMETER :: POINTSMAX=90000
!original start
!trkov    start
      INTEGER(KIND=I4), PARAMETER :: POINTSMAX=800000
!trkov    end
      REAL(KIND=R8), DIMENSION(POINTSMAX) :: XP,YP
!
!     TAGS ON CURRENT RECORD
!
      INTEGER(KIND=I4) :: MAT,MF,MT,NSEQ
!
      REAL(KIND=R4), PARAMETER :: ZERO=0.0
!
!+++MDC+++
!...VMS, ANS, WIN, UNX
!
!     EXECUTE THE STANEF PROGRAM WHEN A STAND ALONE PROGRAM
!
      CALL RUN_STANEF
!
!     TERMINATE JOB
!
      IF(STANEF_SUCCESS.EQ.0) THEN
         WRITE(OUTPUT,'(/A)') ' '
         STOP '    STANEF - Completed successfully'
      ELSE
         WRITE(OUTPUT,'(/A)') ' '
         STOP '    STANEF - Terminated abnormally!'
      END IF
!---MDC---
!
      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE RUN_STANEF
!
!     EXECUTES ENDF FILE STANDARDIZATION PROCESS
!
!-----------------------------------------------------------------------
!     INFIL = INPUT FILE SPECIFICATION
!     OUTFIL =OUTPUT FILE SPECIFICATION
!     LABEL =ID NUMBER OF THE OUTPUT FILE, TREATED AS FOLLOWS
!           =GREATER THAN ZERO---MODIFY TAPE ID OR CREATE A NEW ONE
!           =EQUAL TO ZERO---COPY EXISTING TAPE ID TO NEW FILE
!           =LESS THAN ZERO---NO LABEL ON OUTPUT FILE
!     LTEXT =REVISED LABEL TEXT
!           ="BLANK"---RETAIN TEXT FROM INPUT LABEL
!     MODE  =CONTROL FOR OUTPUT MODE
!           =EQUAL TO ZERO---CHARACTER
!           =EQUAL TO ONE---BINARY
!     INDX  =CONTROL FOR DIRECTORY REVISION OPTION
!           =N---COPY EXISTING DIRECTORY
!           =Y---CREATE OR REVISE THE DIRECTORY
!     IDLTES=CONTROL FOR STANDARDIZATION OF NUMERIC FIELDS
!           =N---DO NOT CONVERT NUMERIC FIELDS
!           =Y---CONVERT NUMERIC FIELDS TO STANDARD FORM
!     I151  =CONRTOL FOR STANDARDIZATION OF THE LIST DATA IN 151
!           =N---DO NOT CONVERT NUMERIC FIELDS
!           =Y---CONVERT NUMERIC FIELDS TO STANDARD FORM
!_______________________________________________________________________
!
      IMPLICIT NONE
!
!     Don't declare TRIM function (causes trouble for gfortran) cmattoon 10/2008
!     CHARACTER(LEN=*), INTRINSIC :: TRIM
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM, ICHAR
!
      LOGICAL(KIND=I4) IEXIST
      CHARACTER(LEN=1) :: IW
      CHARACTER(LEN=80) :: TREC
      CHARACTER(LEN=15) :: PARAMS
      CHARACTER(LEN=4) :: BUF
      CHARACTER(LEN=66) :: CTEXT
      INTEGER(KIND=I4) :: IPER,IC
      INTEGER(KIND=I4) :: J
!
      CHARACTER(LEN=11), PARAMETER :: FIELDF='----------+'
!
!     OUTPUT PROGRAM IDENTIFICATION
!
      STANEF_SUCCESS = 0
      IF(IMDC.LT.4) THEN
         WRITE(OUTPUT,'(/2A)') ' PROGRAM STANEF VERSION ',VERSION
      END IF
!
!     SEE IF INPUT PASSED FROM A DRIVER OR FROM THE COMMAND LINE
!
      IONEPASS = 0
      CALL GET_FROM_CLINE
!
!     INITIALIZE INPUT PARAMETERS
!
   10 IF(IMDC.LT.4) THEN
         STANEF_DATA%INFIL = '*'
         STANEF_DATA%OUTFIL = '*'
         STANEF_DATA%MODE = 0
         STANEF_DATA%LABEL = 0
         STANEF_DATA%LTEXT = ' '
         STANEF_DATA%INDX = 'Y'
         STANEF_DATA%IDLTES = 'Y'
         STANEF_DATA%I151 = 'Y'
      END IF
      SELECT CASE (IMDC)
         CASE (0)
            IW = 'N'
            IONEPASS = 0
         CASE(1,2,3)
            IF(ILENP.NE.0)  THEN
               CALL TOKEN(INPAR,'%',1,STANEF_DATA%INFIL)
               CALL TOKEN(INPAR,'%',2,STANEF_DATA%OUTFIL)
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
         IF(STANEF_DATA%INFIL.EQ.'*') THEN
            IF(IMDC.NE.0) THEN
               WRITE(OUTPUT,FMT=TFMT)                                   &       
     &             ' Input File Specification          - '
            END IF
            READ(INPUT,'(A)') STANEF_DATA%INFIL
         ELSE
            WRITE(OUTPUT,'(/2A)') ' Input file - ',                     &       
     &           TRIM(STANEF_DATA%INFIL)
         END IF
      END IF
!
!     SEE IF INPUT INDICATES FILE TERMINATION
!
      IF(STANEF_DATA%INFIL.EQ.' '.OR.STANEF_DATA%INFIL.EQ.'DONE') THEN
         IF(IONEPASS.EQ.1) STANEF_SUCCESS = 1
         GO TO 100
      END IF
!
!     MAKE SURE INPUT FILE EXISTS
!
      INQUIRE(FILE=STANEF_DATA%INFIL,EXIST=IEXIST)
      IF(.NOT.IEXIST)  THEN
         IF(IMDC.LT.4) THEN
            WRITE(OUTPUT,'(/A/)')  '       COULD NOT FIND INPUT FILE'
         END IF
         SELECT CASE (IMDC)
            CASE (0)
               CLOSE(UNIT=INPUT)
            CASE (1,2,3)
               IF(IONEPASS.EQ.0) GO TO 10
         END SELECT
         STANEF_SUCCESS = 1
         GO TO 90
      END IF
!
!     GET OUTPUT FILE SPECIFICATION
!
      IF(IMDC.LT.4) THEN
         IF(STANEF_DATA%OUTFIL.EQ.'*' ) THEN
            IF(IMDC.NE.0) THEN
               WRITE(OUTPUT,FMT=TFMT)                                   &       
     &           ' Output ENDF File Specification    - '
            END IF
            READ(INPUT,'(A)') STANEF_DATA%OUTFIL
         ELSE
            WRITE(OUTPUT,'(/2A)') ' Output file - ',                    &       
     &             TRIM(STANEF_DATA%OUTFIL)
         END IF
      END IF
      IF(STANEF_DATA%OUTFIL.EQ.' ')  THEN
         IF(IMDC.EQ.1) THEN
            STANEF_DATA%OUTFIL = STANEF_DATA%INFIL
         ELSE
            IPER = INDEX(STANEF_DATA%INFIL,'.')
            IF(IPER.NE.0) THEN
               STANEF_DATA%OUTFIL = STANEF_DATA%INFIL(1:IPER)//'STN'
            ELSE
               IPER = LEN_TRIM(STANEF_DATA%INFIL)
               STANEF_DATA%OUTFIL = STANEF_DATA%INFIL(1:IPER)//'.STN'
            END IF
        END IF
      END IF
!
!     GET TAPE LABEL
!
      IF(IW.NE.'Y') THEN
         IF(IMDC.LT.4) THEN
            IF(IMDC.NE.0) THEN
               WRITE(OUTPUT,FMT=TFMT)                                   &       
     &           '     Enter ENDF Tape Label Number - '
               READ(INPUT,'(BN,I4)',ERR=15) STANEF_DATA%LABEL
               IF(STANEF_DATA%LABEL.GT.0) THEN
                  WRITE(OUTPUT,'(/5X,A/7A)')                            &       
     &            'Enter ENDF Tape Label Text (66 characters maximum)', &       
     &            '        ',(FIELDF,J=1,6)
                  WRITE(OUTPUT,TFMT) '     ** '
                  READ(INPUT,'(A)') STANEF_DATA%LTEXT
              END IF
            ELSE
               READ(INPUT,'(A)') TREC
               CALL TOKEN(TREC,',',1,BUF)
               READ(BUF,'(BN,I4)',ERR=15) STANEF_DATA%LABEL
               CALL TOKEN(TREC,',',2,STANEF_DATA%LTEXT)
            END IF
         END IF
      END IF
!
!     CHECK IF STANDARD OPTIONS ARE WANTED
!
   15 IF(IW.EQ.'*') THEN
         IF(IMDC.GE.1.AND.IMDC.LE.3) THEN
   20       WRITE(OUTPUT,FMT=TFMT) ' Standard Options (Y(es),N(o),?) - '
            READ(INPUT,'(A)')  IW
            IC = ICHAR(IW)
            IF(IC.GT.96.AND.IC.LT.123)   IW = CHAR(IC-32)
            IF(IW.EQ.'?')  THEN
               IW = '*'
               WRITE(OUTPUT,25)
   25          FORMAT(10X,' STANDARD OPTIONS ARE'/                      &       
     &             10X,'    CHARACTER OUTPUT'/                          &       
     &             10X,'    CREATE NEW DIRECTORIES'/                    &       
     &             10X,'    CONVERT DATA FIELDS'/)
               GO TO 20
            END IF
         END IF
      ELSE IF(IMDC.EQ.0) THEN
         READ(INPUT,'(A)') PARAMS
         IF(PARAMS.EQ.' ')  IW = 'Y'
      END IF
!
!     GET OPTIONS
!
      IF(IW.EQ.'N')  THEN
         IF(IMDC.GE.1.AND.IMDC.LE.3) THEN
            WRITE(OUTPUT,TFMT)  '     Output Tape Mode '//              &       
     &                      'CHARACTER(0) or BINARY(1) - '
            READ(INPUT,'(BN,I1)',ERR=40) STANEF_DATA%MODE
         ELSE IF(IMDC.EQ.0) THEN
            CALL TOKEN(PARAMS,',',1,BUF)
            READ(BUF,'(BN,I1)',ERR=40)  STANEF_DATA%MODE
         END IF
  40     IF(STANEF_DATA%MODE.EQ.1) THEN
            STANEF_DATA%INDX = 'N'
            STANEF_DATA%IDLTES = 'N'
            STANEF_DATA%I151 = 'N'
         ELSE
            STANEF_DATA%MODE = 0
            IF(IMDC.GE.1.AND.IMDC.LE.3) THEN
               WRITE(OUTPUT,TFMT)   '     Create a New Directory '//    &       
     &                             '(Y(es),N(o))? - '
               READ(INPUT,'(A)') STANEF_DATA%INDX
            ELSE IF(IMDC.EQ.0) THEN
               CALL TOKEN(PARAMS,',',2,STANEF_DATA%INDX)
            END IF
            IC = ICHAR(STANEF_DATA%INDX)
            IF(IC.GT.96.AND.IC.LT.123)    STANEF_DATA%INDX = CHAR(IC-32)
            IF(STANEF_DATA%INDX.NE.'N') STANEF_DATA%INDX = 'Y'
            IF(IMDC.GE.1.AND.IMDC.LE.3) THEN
               WRITE(OUTPUT,TFMT)  '     Convert Data Fields to '//     &       
     &                            'Standard Form (Y(es),N(o))? - '
               READ(INPUT,'(A)') STANEF_DATA%IDLTES
            ELSE IF(IMDC.EQ.0) THEN
               CALL TOKEN(PARAMS,',',3,STANEF_DATA%IDLTES)
            END IF
            IC = ICHAR(STANEF_DATA%IDLTES)
            IF(IC.GT.96.AND.IC.LT.123) STANEF_DATA%IDLTES = CHAR(IC-32)
            IF(STANEF_DATA%IDLTES.NE.'N') STANEF_DATA%IDLTES = 'Y'
            IF(STANEF_DATA%IDLTES.EQ.'Y') THEN
               IF(IMDC.GE.1.AND.IMDC.LE.3) THEN
                  WRITE(OUTPUT,TFMT)   '     Convert Resonance '//      &       
     &                       'Parameter Fields to Standard '//          &       
     &                                'Form (Y(es),N(o))? - '
                  READ(INPUT,'(A)') STANEF_DATA%I151
               ELSE IF(IMDC.EQ.0) THEN
                  CALL TOKEN(PARAMS,',',4,STANEF_DATA%I151)
               END IF
               IC = ICHAR(STANEF_DATA%I151)
               IF(IC.GT.96.AND.IC.LT.123) STANEF_DATA%I151 = CHAR(IC-32)
               IF(STANEF_DATA%I151.NE.'N') STANEF_DATA%I151 = 'Y'
            END IF
         END IF
      END IF
!
!     OPEN INPUT AND OUTPUT FILES
!
   50 OPEN(UNIT=ITAPE,ACCESS='SEQUENTIAL',STATUS='OLD',                 &       
     &   FILE=STANEF_DATA%INFIL,ACTION='READ')
      IF(STANEF_DATA%MODE.EQ.0) THEN
!+++MDC+++
!...VMS
!/         OPEN(UNIT=OTAPE,ACCESS='SEQUENTIAL',STATUS=OSTATUS,          &       
!/     &       FILE=STANEF_DATA%OUTFIL,CARRIAGECONTROL='LIST')
!...WIN, DVF, UNX, LWI, ANS, MOD
         OPEN(UNIT=OTAPE,ACCESS='SEQUENTIAL',STATUS=OSTATUS,            &       
     &       FILE=STANEF_DATA%OUTFIL)
!---MDC---
      ELSE
         OPEN(UNIT=OTAPE,ACCESS='SEQUENTIAL',STATUS=OSTATUS,            &       
     &       FILE=STANEF_DATA%OUTFIL,FORM='UNFORMATTED')
      END IF
!
!     OUTPUT SELECTED OPTIONS
!
      IF(IMDC.LT.4) THEN
         WRITE(OUTPUT,'(/A)')' OPTION(S) SELECTED'
         IF(STANEF_DATA%LABEL.LT.0) THEN
            WRITE(OUTPUT,'(15X,A)') 'OUTPUT FILE WILL NOT BE LABELED'
         ELSE IF(STANEF_DATA%LABEL.EQ.0) THEN
            WRITE(OUTPUT,'(15X,A)') 'OUTPUT FILE WILL HAVE THE SAME'//  &       
     &                   ' LABEL AS THE INPUT'
         ELSE
            WRITE(OUTPUT,'(15X,A)') 'OUTPUT FILE WILL HAVE A NEW LABEL'
         END IF
         IF(STANEF_DATA%MODE.EQ.0) THEN
            WRITE(OUTPUT,'(15X,A)') 'MATERIALS WILL BE RESEQUENCED'
            IF(STANEF_DATA%INDX.EQ.'Y')  THEN
               WRITE(OUTPUT,'(15X,A)') 'A NEW INDEX WILL BE GENERATED'
            END IF
            IF(STANEF_DATA%IDLTES.EQ.'Y')   THEN
               IF(STANEF_DATA%I151.NE.'N') THEN
                  WRITE(OUTPUT,'(15X,A)')                               &       
     &                 'DATA FIELDS WILL BE CONVERTED TO STANDARD FORM'
               ELSE
                  WRITE(OUTPUT,'(15X,A/17X,A)')                         &       
     &                 'DATA FIELDS WILL BE CONVERTED TO STANDARD FORM',&       
     &                 'EXCEPT FOR MT=151'
               END IF
            END IF
         ELSE
            WRITE(OUTPUT,'(15X,A)')                                     &       
     &                'FILE WILL BE CONVERTED TO BINARY FORMAT'
         END IF
      END IF
!
!     IF NEW INDEX IS TO BE GENERATED CREATE CONTENTS IN A SCRATCH FILE
!      IN A PASS THROUGH THE TAPE
!
      IF(STANEF_DATA%INDX.EQ.'Y')   CALL STRIP
!
!     PROCESS TAPE LABEL
!
      READ(ITAPE,'(A,I4,I2,I3)') CTEXT,MAT,MF,MT
      IF(MF.NE.0.OR.MT.NE.0)   THEN
!        -- If MF>0 or MT>0 on first record, assume no header record
         REWIND (UNIT=ITAPE)
         CTEXT = ' '
         MAT = 0
      ELSE
!        -- If tape label is zero, change to 7777
         IF(MAT.EQ.0) MAT=7777
      END IF
!
!     NEW TAPE TO BE LABELED
!
      IF(STANEF_DATA%LABEL.GE.0)  THEN
!
!        COPY EXISTING LABEL?
!
         IF(STANEF_DATA%LABEL.EQ.0) THEN
            IF(MAT.EQ.0) GO TO 70
            STANEF_DATA%LABEL = MAT
            STANEF_DATA%LTEXT = CTEXT
         END IF
!
!        WRITE NEW LABEL IF NEEDED
!
         IF(STANEF_DATA%MODE.EQ.0)  THEN
            WRITE(OTAPE,'(A,I4,A)')  STANEF_DATA%LTEXT,                 &       
     &               STANEF_DATA%LABEL,' 0  0    0'
         ELSE
            WRITE(OTAPE) STANEF_DATA%LABEL,ZERO,ZERO,STANEF_DATA%LTEXT
         END IF
      END IF
!
!     CONSTRUCT NEW TAPE
!
   70 CALL COMPOZ
!
!     CLOSE DATA FILES
!
      IF(STANEF_DATA%INDX.EQ.'Y')    CLOSE(UNIT=IDIR,STATUS='DELETE')
      CLOSE(UNIT=ITAPE)
      CLOSE(UNIT=OTAPE)
!
   90 IF(IONEPASS.EQ.0) GO TO 10
!
!     TERMINATE JOB
!
  100 RETURN
!
      END SUBROUTINE RUN_STANEF
!
!***********************************************************************
!
      SUBROUTINE STRIP
!
!     ROUTINE TO PASS THROUGH AN ENDF TAPE TO EXTRACT DATA REQUIRED
!      FOR DIRECTORY GENERATION
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MFO,MTO
      INTEGER(KIND=I4) :: NC,MFT
!
!     OPEN SCRATCH FILE
!
      OPEN(UNIT=IDIR,STATUS='SCRATCH',FORM='FORMATTED')
!
!     CHECK FOR A TAPE LABEL
!
      NC = 0
      READ(ITAPE,'(66X,I4,I2,I3)')  MAT,MFO,MTO
      IF(MFO.NE.0.OR.MTO.NE.0)  NC = 1
!
!     PROCESS EACH RECORD ON TAPE
!
   20 READ(ITAPE,'(66X,I4,I2,I3)')  MAT,MF,MT
!
!     CHECK MATERIAL
!
      IF(MAT.LT.0) THEN     ! END OF TAPE
         GO TO 100
      ELSE IF(MAT.EQ.0) THEN     ! END OF MATERIAL
         WRITE(IDIR,'(A)') '    0      0     0'
         GO TO 20
      END IF
!
!     SKIP END OF FILE RECORD
!
      IF(MF.EQ.0)   GO TO 20
!
!     CHECK SECTION ID AND COUNT RECORDS IN A SECTION
!
      IF(MT.NE.0)   THEN
!
!     IF FIRST RECORD OF A NEW SECTION SO RESET MF AND MT
!
         IF(NC.EQ.0)   THEN
            MFO = MF
            MTO = MT
         END IF
         NC = NC + 1
      ELSE
!
!     END OF SECTION, SO PUT RECORD IN SCRATCH FILE
!
         MFT = 1000*MFO + MTO
         WRITE(IDIR,'(I5,I7,I6)')  MAT,MFT,NC
         NC = 0
      END IF
      GO TO 20
!
!     END OF TAPE
!
  100 ENDFILE (UNIT=IDIR)
      REWIND (UNIT=IDIR)
      REWIND (UNIT=ITAPE)
!
      RETURN
      END SUBROUTINE STRIP
!
!***********************************************************************
!
      SUBROUTINE COMPOZ
!
!     ROUTINE TO CONTROL CONSTRUCTION OF NEW ENDF TAPE
!
      IMPLICIT NONE
!
!     SIMPLE RESEQUENCING OF ASCII FORMAT FILE
!
      IF((STANEF_DATA%INDX.EQ.'N'.AND.STANEF_DATA%IDLTES.EQ.'N')        &       
     &        .AND.STANEF_DATA%MODE.EQ.0)      THEN
         NFOR = 6
         DO WHILE (MAT.GE.0)
            CALL COPYMA(1)
         END DO
         GO TO 100
      ELSE
         NSEQ = 0
         NFOR = 0
      END IF
!
!     READ FIRST RECORD OF NEXT FILE
!
   10 CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      IF(MAT.EQ.0) THEN
         NSEQ = 0
         NFOR = 0
         GO TO 10
      ELSE IF(MAT.LT.0) THEN
         GO TO 100
      ENDIF
!
!     BRANCH ON FILE TYPE
!
      SELECT CASE (MF)
!
         CASE (1)
            CALL FILE1
!*****IF FORMAT NOT CHANGED SIMPLY COPY AND RESEQUENCE REST OF MATERIAL
            IF(STANEF_DATA%MODE.EQ.0.AND.STANEF_DATA%IDLTES.EQ.'N') THEN
               CALL COPYMA(0)
               NSEQ = 0
               NFOR = 0
               GO TO 10
            END IF
         CASE (2)
            CALL FILE2
         CASE (3)
            CALL FILE3
         CASE (4)
            CALL FILE4
         CASE (5)
            CALL FILE5
         CASE (6)
            CALL FILE6
         CASE (7)
            CALL FILE7
         CASE (8)
            CALL FILE8
         CASE (9,10)
            CALL FILE9
         CASE (12,13)
            CALL FILE12
         CASE (14)
            CALL FILE14
         CASE (15)
            CALL FILE15
         CASE (23,27)
            CALL FILE23
         CASE (26)
            CALL FILE26
         CASE (28)
            CALL FILE28
         CASE (31,33)
            CALL FILE33
         CASE (32)
            CALL FILE32
         CASE (34)
            CALL FILE34
         CASE (35)
            CALL FILE35
         CASE (40)
            CALL FILE40
!
!       COPY FILE ONLY RESEQUENCING WHEN FILE TYPE NOT LEGAL
!
         CASE DEFAULT
            CALL COPYFL
!
      END SELECT
      GO TO 10
!
  100 RETURN
      END SUBROUTINE COMPOZ
!
!***********************************************************************
!
      SUBROUTINE COPYFL
!
!     ROUTINE TO COPY A FILE WITH NO CHANGE TO A NEW TAPE EXCEPT FOR
!       RESEQUENCING
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IOUTS
!
!     OUTPUT ERROR MESSAGE
!
      IOUTS = 1
      IF(IMDC.LT.4) THEN
         WRITE(OUTPUT,'(/5X,A,I2)') 'COPYING UNRECOGNIZED FILE ',MF
      END IF
!
!     READ EACH RECORD
!
   10 READ(ITAPE,'(A,I4,I2,I3,I5)') TEXT,MAT,MF,MT
      IF(IOUTS.EQ.1) THEN
         CALL OUT_STATUS
         IOUTS = 0
      END IF
!
!     WRITE EACH RECORD
!
      IF(STANEF_DATA%MODE.EQ.0)   THEN
         IF(NFOR.LT.6) THEN
           IF(MT.EQ.0.OR.MF.EQ.0) IOUTS = 1
           IF(MAT.GE.0)   THEN
              NSEQ = NSEQ+1
              IF(NSEQ.GT.99999) NSEQ=0
           ELSE
              NSEQ = 0
           END IF
         ELSE
            IF(MAT.LT.0) THEN
               NSEQ = 0
            ELSE IF(MF.EQ.0) THEN
               NSEQ = 0
               IOUTS = 1
            ELSE IF(MT.EQ.0) THEN
               NSEQ = 99999
               IOUTS = 1
            ELSE
               NSEQ = NSEQ + 1
            END IF
         END IF
         WRITE(OTAPE,'(A,I4,I2,I3,I5)') TEXT,MAT,MF,MT,NSEQ
         IF(NFOR.GT.5.AND.MT.EQ.0) NSEQ = 0
      ELSE
         IF(MF.NE.0)   GO TO 10
      END IF
!
      RETURN
      END SUBROUTINE COPYFL
!
!***********************************************************************
!
      SUBROUTINE COPYMA(ISEQ)
!
!     ROUTINE TO COPY A MATERIAL WITH NO CHANGE TO A NEW TAPE EXCEPT FOR
!       RESEQUENCING.
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: ISEQ
!
      INTEGER(KIND=I4), INTRINSIC :: MIN0
!
      INTEGER(KIND=I4) :: IOUTS
!
!     INITIALIZE SEQUENCE NUMBER FOR EACH NEW MATERIAL ENCOUNTERED
!
      IF(ISEQ.EQ.1)  THEN
         NSEQ = 0
         IF(IMDC.LT.4) WRITE(OUTPUT,'(/)')
         CALL OUT_STATUS
      END IF
      IOUTS = 0
!
!     READ EACH RECORD
!
   10 READ(ITAPE,'(A,I4,I2,I3,I5)') TEXT,MAT,MF,MT
      IF(IOUTS.EQ.1) THEN
         CALL OUT_STATUS
         IOUTS = 0
      END IF
!
!     WRITE EACH RECORD
!
      IF(NFOR.LT.6) THEN
        IF(MT.EQ.0.OR.MF.EQ.0) IOUTS = 1
        IF(MAT.GE.0) THEN
           NSEQ = MIN0(NSEQ+1,99999)
        ELSE
           NSEQ = 0
        END IF
      ELSE
         IF(MAT.LT.0) THEN
            NSEQ = 0
         ELSE IF(MF.EQ.0) THEN
            NSEQ = 0
            IOUTS = 1
         ELSE IF(MT.EQ.0) THEN
            NSEQ = 99999
            IOUTS = 1
         ELSE
            NSEQ = NSEQ + 1
         END IF
      END IF
      WRITE(OTAPE,'(A,I4,I2,I3,I5)') TEXT,MAT,MF,MT,NSEQ
      IF(NFOR.GT.5.AND.MT.EQ.0) NSEQ = 0
      IF(MAT.GT.0)   GO TO 10
!
      RETURN
      END SUBROUTINE COPYMA
!
!***********************************************************************
!
      SUBROUTINE FILE1
!
!     ROUTINE TO PROCESS FILE 1
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4), INTRINSIC :: IFIX, MOD, MIN0
!
      CHARACTER(LEN=1) :: IREV
      CHARACTER(LEN=4) :: ITEX
      CHARACTER(LEN=11) :: FLOAT1,FLOAT2
      CHARACTER(LEN=66), DIMENSION(3) :: TEM
      CHARACTER(LEN=66), DIMENSION(6) :: TEXTS
      CHARACTER(LEN=11) :: ZSA
      INTEGER(KIND=I4) :: NMOD,NLIB,NSUB,NVER
      integer(kind=i4) :: iza,iz,ia,imeta,iz1
      INTEGER(KIND=I4) :: NCDS,NXC,N2
      INTEGER(KIND=I4) :: L1,L2,N1
      INTEGER(KIND=I4) :: MATD,NTRD,NPROC,KL,NB,NDIF,NSEC
      INTEGER(KIND=I4) :: NOFF,NFREV,NCC,MODC
      INTEGER(KIND=I4) :: MFTC,MFC,MTC
      INTEGER(KIND=I4) :: K,N,NN
      INTEGER(KIND=I4) :: LFC, NFC
      REAL(KIND=R8) :: C1,C2,ZA
      REAL(KIND=4) ::ZZA
!
      INTEGER(KIND=I4), PARAMETER :: NSECMAX=1000
      INTEGER(KIND=I4), DIMENSION(NSECMAX) :: MFTD,NCD,SMODD
!
!     DEFINE ELEMENT SYMBOLS
!
      INTEGER(KIND=I4), PARAMETER :: IELM=113
      CHARACTER(LEN=2), DIMENSION(IELM), PARAMETER ::                   &       
     &   ELEMNT = (/'nn',                                               &       
     &       'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',         &       
     &       'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',         &       
     &       'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',         &       
     &       'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',         &       
     &       'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',         &       
     &       'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',         &       
     &       'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',         &       
     &       'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',         &       
     &       'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',         &       
     &       'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',         &       
     &       'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',         &       
     &       'Rg','XX'/)
!
      IF(IMDC.LT.4) WRITE(OUTPUT,'(/)')
      CALL OUT_STATUS
!
!     SAVE MATERIAL MOD NUMBER AND LIBRARY NUMBER
!
      ZA = C1H
      NMOD = N2H
      NLIB = N1H
!
!     READ SECOND CONT RECORD
!
      CALL CONT(C1,C2,L1,L2,N1,NFOR)
      IMETA = L2
      IF(NFOR.EQ.0)  NFOR = 5
!
!     BUILD Z-S-A FOR CARD 5
!
!
!     READ THIRD CONT RECORD IF ENDF-6 OR LATER
!
      IF(NFOR.GT.5)  THEN
         CALL CONT(C1,C2,L1,L2,NSUB,NVER)
      ELSE
         SELECT CASE (NLIB)
            CASE (2:4,35)
               NVER = 1
            CASE (5)
               NVER = 2
            CASE(6)
               NVER = 3
            CASE DEFAULT
               NVER = 5
         END SELECT
         NSUB = 10
      END IF
      ZSA = ' '
!     -- Special case: suppress redefining ZSA if ZA=MAT
      IF(NSUB.NE.12 .AND. NINT(ZA).NE.MAT) THEN
         ZZA = ZA
         IZA = IFIX(ZZA+.001)
         IA = MOD(IZA,1000)
         IZ = IZA/1000
         IZ1 = MIN0((IZ+1),IELM)
         WRITE(ZSA,'(I3,3A,I3)') IZ,'-',ELEMNT(IZ1),'-',IA
         IF(IMETA.GE.3) THEN
            ZSA(11:11) = 'O'
         ELSE IF(IMETA.GE.2) THEN
            ZSA(11:11) = 'N'
         ELSE IF(IMETA.GE.1) THEN
           ZSA(11:11) = 'M'
         END IF
      END IF
!
!     DIRECTORY NOT REVISED SO COPY NEXT CONT RECORD
!
      IF(STANEF_DATA%INDX.NE.'Y')  THEN
         CALL CONT(C1,C2,L1,L2,NCDS,NXC)
         NB = 1
         NDIF = 0
         GO TO 60
      ELSE
         NSEC = 0
      END IF
!
!     LOAD IN REVISED DIRECTORY CONTENTS
!
   10 READ(IDIR,'(I5,I7,I6)')  MATD,MFTD(NSEC+1),NCD(NSEC+1)
      SMODD(NSEC+1) = NMOD
      IF(MATD.NE.0)   THEN
         NSEC = NSEC + 1
         GO TO 10
      END IF
!
!     READ NEXT CONT CARD AND UP TO SIX HOLERITH COMMENTS
!
      READ(ITAPE,'(2E11.4,4I11,I4,I2,I3)') C1,C2,L1,L2,N1,N2,MAT,MF,MT
      NTRD = 0
      NDIF = 0
!
!     FOR ENDF-6 AND LATER FORMATS FIND SPECIAL FORMAT HOLERITH
!
      IF(NFOR.GE.6)   THEN
         NTRD = MIN0(N1,6)
         NPROC = 0
         DO N=1,NTRD
            NPROC = NPROC + 1
            READ(ITAPE,'(A)')   TEXTS(NPROC)
!           -- Special case: suppress redefining ZSA if ZA=MAT
            IF(NPROC.EQ.1.AND.NSUB.NE.12.AND.NINT(ZA).NE.MAT) THEN
               TEXTS(NPROC)(1:11) = ZSA
!           ELSE
!              CALL UPSTR(TEXTS(NPROC)(1:11))
            END IF
            ITEX = TEXTS(NPROC)(1:4)
            IF(ITEX.EQ.'----')   NPROC = NPROC - 1
            IF(NPROC.EQ.3)   THEN
               NTRD = N
               TEXTS(6) = TEXTS(3)
               GO TO 20
            END IF
         END DO
!
!        LEAVE SPACE FOR INSERTION OF SPECIAL HOLERITH RECORDS
!
   20    NDIF = NPROC - NTRD + 3
         IF(NDIF.NE.0)  THEN
            N1 = N1 + NDIF
            NCD(1) = NCD(1) + NDIF
         END IF
!
!        GENERATE SPECIAL HOLERITH RECORDS
!
         IF(NPROC.GE.2)  THEN
            IREV = TEXTS(2)(37:37)
         ELSE
            IREV = ' '
         END IF
         CALL TEXID(NFOR,NLIB,NVER,NSUB,MAT,IREV,TEM)
         NOFF = MIN0(NPROC,2)
         DO NN=1,3
            TEXTS(NN+NOFF) = TEM(NN)
         END DO
         NFREV = NPROC + 3
      END IF
!
!     SET INTEGERS FOR CONT RECORD
!
      NSEQ = MIN0(NSEQ+1,99999)
      NCDS = N1
      NXC = N2
!
!     FORMAT THE FLOATING POINT NUMBERS
!
      CALL NORMAL(C1,FLOAT1)
      CALL NORMAL(C2,FLOAT2)
!
!     OUTPUT THE CONT RECORD
!
      N2 = NSEC
      IF(STANEF_DATA%MODE.EQ.0)  THEN
         WRITE(OTAPE,'(2A,4I11,I4,I2,I3,I5)')  FLOAT1,FLOAT2,           &       
     &          L1,L2,N1,N2,MAT,MF,MT,NSEQ
      ELSE
         WRITE(OTAPE) MAT,MF,MT,C1,C2,L1,L2,N1,N2
      END IF
!
!     OUTPUT FIRST REVISED HOLERITH RECORDS
!
      NB = NTRD + 1
      IF(NB.GT.1)   THEN
         DO N=1,NFREV
            TEXT = TEXTS(N)
            CALL TEXTR(2)
         END DO
      END IF
!
!     PROCESS COMMENTS
!
   60 DO N=NB,NCDS-NDIF
         CALL TEXTR(1)
      END DO
!
!     PROCESS IF OLD DIRECTORY IS RETAINED
!
      IF(STANEF_DATA%INDX.EQ.'N')     THEN
         IF(NXC.GT.0)   THEN
            DO N=1,NXC
               READ(ITAPE,'(22X,3I11,I11,I4,I2,I3,I5)')  MFC,MTC,       &       
     &                NCC,MODC,MAT,MF,MT
               NSEQ = MIN0(NSEQ+1,99999)
               IF(STANEF_DATA%MODE.EQ.0)  THEN
                  WRITE(OTAPE,'(22X,3I11,I11,I4,I2,I3,I5)') MFC,MTC,    &       
     &                NCC,MODC,MAT,MF,MT,NSEQ
               ELSE
                  WRITE(OTAPE) MAT,MF,MT,ZERO,ZERO,MFC,MTC,NCC,MODC
               END IF
            END DO
         END IF
!
!     REVISED DIRECTORY
!
      ELSE
         NCD(1) = NCD(1) - NXC + NSEC
         IF(NXC.GT.0) THEN
            KL = 1
            DO N=1,NXC
               READ(ITAPE,'(22X,3I11,I11,I4,I2,I3,I5)')                 &       
     &                  MFC,MTC,NCC,MODC
               MFTC = 1000*MFC + MTC
               DO K=KL,NSEC
                  IF(MFTC.EQ.MFTD(K)) THEN
                     KL = K + 1
                     SMODD(K) = MODC
                     GO TO 50
                  ELSE IF(MFTC.GT.MFTD(K)) THEN
                     GO TO 50
                  END IF
               END DO
   50       END DO
         END IF
!
!        OUTPUT REVISED DIRECTORY
!
         DO N=1,NSEC
            NSEQ = MIN0(NSEQ+1,99999)
            MFC = MFTD(N)/1000
            MTC = MOD(MFTD(N),1000)
            WRITE(OTAPE,'(22X,3I11,I11,I4,I2,I3,I5)')  MFC,MTC,         &       
     &                NCD(N),SMODD(N),MAT,MF,MT,NSEQ
         END DO
      END IF
!
!     PROCESS SEND CARD
!
      CALL CONT(C1,C2,L1,L2,N1,N2)
!
!     IF DATA FIELDS ARE NOT TO BE CONVERTED, RETURN
!
      IF(STANEF_DATA%IDLTES.EQ.'Y'.OR.STANEF_DATA%MODE.EQ.1) THEN
!
!        PROCESS NEXT HEAD CARD
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
!
!        PROCESS ANY ADDITIONAL SECTIONS IN FILE 1
!
         IF(MF.NE.0)   CALL FILE1R
      END IF
!
      RETURN
      END SUBROUTINE FILE1
!
!***********************************************************************
!
      SUBROUTINE TEXID(NFOR,NLIB,NVER,NSUB,MAT,IREV,CARRAY)
!
!     ROUTINE TO BUILD FREE TEXT IDENTIIFIERS FOR INSERTION INTO 1/451
!
      IMPLICIT NONE
!
      CHARACTER(LEN=1) :: IREV
      CHARACTER(LEN=66), DIMENSION(3) :: CARRAY
      INTEGER(KIND=I4) :: NFOR,NLIB,NVER,NSUB,MAT
!
      INTEGER(KIND=I4), INTRINSIC :: MOD
!
      CHARACTER(LEN=24) :: BCD
      INTEGER(KIND=I4) :: NCBCD
      INTEGER(KIND=I4) :: IPART,ITYPE,IERR,N
!
      INTEGER(KIND=I4), PARAMETER :: NSUBS=21
      CHARACTER(LEN=40),DIMENSION(NSUBS) :: TSUBS
      DATA TSUBS/'PHOTONUCLEAR DATA                       ',            &       
     &           'PHOTON-INDUCED FISSION PRODUCT YIELDS   ',            &       
     &           'PHOTO-ATOMIC INTERACTION DATA           ',            &       
     &           'RADIOACTIVE DECAY DATA                  ',            &       
     &           'SPONTANEOUS FISSION PRODUCT YIELDS      ',            &       
     &           'ATOMIC RELAXATION DATA                  ',            &       
     &           'INCIDENT NEUTRON DATA                   ',            &       
     &           'NEUTRON-INDUCED FISSION PRODUCT YIELDS  ',            &       
     &           'THERMAL NEUTRON SCATTERING DATA         ',            &       
     &           'ELECTRO-ATOMIC INTERACTION DATA         ',            &       
     &           'INCIDENT PROTON DATA                    ',            &       
     &           'PROTON-INDUCED FISSION PRODUCT YIELDS   ',            &       
     &           'INCIDENT DEUTERON DATA                  ',            &       
     &           'DEUTERON-INDUCED FISSION PRODUCT YIELDS ',            &       
     &           'INCIDENT TRITON DATA                    ',            &       
     &           'TRITON-INDUCED FISSION PRODUCT YIELDS   ',            &       
     &           'INCIDENT HELIUM-3 DATA                  ',            &       
     &           'HELIUM-3-INDUCED FISSION PRODUCT YIELDS ',            &       
     &           'INCIDENT ALPHA DATA                     ',            &       
     &           'ALPHA-INDUCED FISSION PRODUCT YIELDS    ',            &       
     &           'UNKNOWN SUBLIBRARY                      '/
      INTEGER(KIND=I4), DIMENSION(NSUBS), PARAMETER ::                   &      
     &    ISUBS = (/0,1,3,4,5,6,10,11,12,113,10010,10011,10020,10021,   &       
     &          10030,10031,20030,20031,20040,20041,99999/)
!
!     INITIALIZE TEXT
!
      CARRAY(1) = '----'
      CARRAY(2) = '-----'
      CARRAY(3) = '------ENDF-  FORMAT'
!
!     BUILD LIB ID
!
      CALL LIBID(NLIB,NVER,CARRAY(1)(5:18))
      CARRAY(1)(23:30) = 'MATERIAL'
      WRITE(CARRAY(1)(32:35) ,'(I4)') MAT
      IF(IREV.NE.' ')  CARRAY(1)(45:54) = 'REVISION '//IREV
!
!     BUILD SUBLIB ID
!
      DO N=1,NSUBS-1
         IF(NSUB.EQ.ISUBS(N)) THEN
            CARRAY(2)(6:45) = TSUBS(N)
            GO TO 20
         END IF
      END DO
      IPART = NSUB/10
      ITYPE = MOD(NSUB,10)
      IF(ITYPE.GT.1)  THEN
         N = NSUBS
         CARRAY(2)(6:45) = TSUBS(N)
         GO TO 20
      END IF
      CALL ZAID(IPART,BCD,NCBCD,IERR)
      IF(IERR.EQ.1)   THEN
         N = NSUBS
         CARRAY(2)(6:45) = TSUBS(N)
         GO TO 20
      END IF
      IF(ITYPE.EQ.1)  THEN
         CARRAY(2)(6:61) = BCD(1:NCBCD)//                               &       
     &              '-INDUCED FISSION PRODUCT YIELDS'
      ELSE
         CARRAY(2)(6:61) = 'INCIDENT '//BCD(1:NCBCD)//' DATA'
      END IF
!
!     GENERATE FORMAT LINE
!
   20 WRITE(CARRAY(3)(12:12),'(I1)')  NFOR
!
      RETURN
      END SUBROUTINE TEXID
!
!***********************************************************************
!
      SUBROUTINE LIBID(NLIB,NVER,LIBBCD)
!
!     ROUTINE TO GENERATE LIBRARY ID FROM NLIB AND NVER
!
      IMPLICIT NONE
!
      CHARACTER(LEN=14) :: LIBBCD
      INTEGER(KIND=I4) :: NLIB,NVER
!
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM
!
      CHARACTER(LEN=1) :: TCHAR
      CHARACTER(LEN=14) :: LIBTMP
      INTEGER(KIND=I4) :: INVER,LNC,NCF
      INTEGER I,N
!
      INTEGER(KIND=I4), PARAMETER :: NLIBS=13
      INTEGER(KIND=I4), DIMENSION(NLIBS), PARAMETER ::                   &      
     &        LIBS = (/0,1,2,3,4,5,6,31,32,33,34,41,99/)
      CHARACTER(LEN=8), DIMENSION(NLIBS), PARAMETER ::                  &       
     &      TLIBS = (/'ENDF/B  ','ENDF/A  ','JEF     ','EFF     ',      &       
     &                'ENDF/HE ','CENDL   ','JENDL   ','INDL/V  ',      &       
     &                'INDL/A  ','INDL/F  ','IRDF    ','BROND   ',      &       
     &                'UNKNOWN '/)
!
      CHARACTER(LEN=5), DIMENSION(20), PARAMETER ::                     &       
     &     NUMS = (/'I    ','II   ','III  ','IV   ','V    ','VI   ',    &       
     &              'VII  ','VIII ','IX   ','X    ','XI   ','XII  ',    &       
     &              'XIII ','XIV  ','XV   ','XVI  ','XVII ','XVIII',    &       
     &              'XIX  ','?????'/)
!
!     FIND LIBRARY ID
!
      DO I=1,NLIBS-1
         IF(NLIB.EQ.LIBS(I))    GO TO 10
      END DO
      I = NLIBS
!
!     BUILD LIBRARY ID
!
   10 IF(NLIB.EQ.0.OR.NLIB.EQ.4) THEN
         IF(NVER.LT.1.OR.NVER.GT.19)  THEN
            INVER = 20
         ELSE
            INVER = NVER
         END IF
         LNC = LEN_TRIM(TLIBS(I))
         WRITE(LIBBCD,'(3A)')   TLIBS(I)(1:LNC),'-',NUMS(INVER)
      ELSE
         WRITE(LIBTMP,'(2A,I5)') TLIBS(I),'-',NVER
         LIBBCD = ' '
         NCF = 0
         DO N=1,14
            TCHAR = LIBTMP(N:N)
            IF(TCHAR.NE.' ')THEN
               NCF = NCF + 1
               LIBBCD(NCF:NCF) = TCHAR
            END IF
         END DO
      END IF
!
  100 RETURN
      END SUBROUTINE LIBID
!
!***********************************************************************
!
      SUBROUTINE ZAID(ZA,BCD,NCBCD,IERR)
!
!     CONSTRUCTION OF THE ALPHANUMERIC REPRESENTATION OF AN INCIDENT
!     PARTICLE HAVING A GIVEN ATOMIC NUMBER AND WEIGHT. IF THE ATOMIC
!     NUMBER IS EQUAL TO ZERO THE ISOTOPE IS ASSUMED TO BE A NATURAL
!     ELEMENT
!
!     SUBROUTINE ARGUMENTS ARE DEFINED AS FOLLOWS.....
!     ZA    = 1000*Z + A
!     BCD   =RESULTING ALPHANUMERIC EQUIVALENT.
!     NCBCD  =NUMBER OF CHARACTERS IN BCD STRING
!     IERR  =ERROR INDICATOR. SET EQUAL TO ONE IF IZ IS NOT IN THE RANGE
!            1 TO 111, OR IA IS NOT IN THE RANGE 1 TO 999. IT IS
!            SET EQUAL TO ZERO IF BOTH ARE IN ACCEPTABLE RANGE.
!
!     THE NAME OF THE ELEMENT FOLLOWED BY ITS ATOMIC NUMBER IS LEFT
!     ADJUSTED INTO ID. FOR NATURAL MATERIALS THE WORD 'NATURAL'
!     FOLLOWED BY THE ELEMENT NAME IS USED INSTEAD.
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: BCD
      INTEGER(KIND=I4) :: ZA,NCBCD,IERR
!
      CHARACTER(LEN=1), INTRINSIC :: CHAR
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM, ICHAR
!
      CHARACTER(LEN=24) :: BCDX
      CHARACTER(LEN=1) :: CCHAR
      INTEGER(KIND=I4) :: IZ,IA,NC,IC
      INTEGER(KIND=I4) :: IJK,I
!
!     DEFINE ELEMENT NAMES
!
      INTEGER(KIND=I4), PARAMETER :: IELM=111
      CHARACTER(LEN=13), DIMENSION(IELM) :: ELEMNT
      DATA (ELEMNT(IJK),IJK=1,10)/   'Hydrogen    ',                    &       
     &                               'Helium      ',                    &       
     &                               'Lithium     ',                    &       
     &                               'Beryllium   ',                    &       
     &                               'Boron       ',                    &       
     &                               'Carbon      ',                    &       
     &                               'Nitrogen    ',                    &       
     &                               'Oxygen      ',                    &       
     &                               'Fluorine    ',                    &       
     &                               'Neon        '/
      DATA (ELEMNT(IJK),IJK=11,20)/  'Sodium      ',                    &       
     &                               'Magnesium   ',                    &       
     &                               'Aluminum    ',                    &       
     &                               'Silicon     ',                    &       
     &                               'Phosphorus  ',                    &       
     &                               'Sulfur      ',                    &       
     &                               'Chlorine    ',                    &       
     &                               'Argon       ',                    &       
     &                               'Potassium   ',                    &       
     &                               'Calcium     '/
      DATA (ELEMNT(IJK),IJK=21,30)/  'Scandium    ',                    &       
     &                               'Titanium    ',                    &       
     &                               'Vanadium    ',                    &       
     &                               'Chromium    ',                    &       
     &                               'Manganese   ',                    &       
     &                               'Iron        ',                    &       
     &                               'Cobalt      ',                    &       
     &                               'Nickel      ',                    &       
     &                               'Copper      ',                    &       
     &                               'Zinc        '/
      DATA (ELEMNT(IJK),IJK=31,40)/  'Gallium     ',                    &       
     &                               'Germanium   ',                    &       
     &                               'Arsenic     ',                    &       
     &                               'Selenium    ',                    &       
     &                               'Bromine     ',                    &       
     &                               'Krypton     ',                    &       
     &                               'Rubidium    ',                    &       
     &                               'Strontium   ',                    &       
     &                               'Yttrium     ',                    &       
     &                               'Zirconium   '/
      DATA (ELEMNT(IJK),IJK=41,50)/  'Niobium     ',                    &       
     &                               'Molybdenum  ',                    &       
     &                               'Technetium  ',                    &       
     &                               'Ruthenium   ',                    &       
     &                               'Rhodium     ',                    &       
     &                               'Palladium   ',                    &       
     &                               'Silver      ',                    &       
     &                               'Cadmium     ',                    &       
     &                               'Indium      ',                    &       
     &                               'Tin         '/
      DATA (ELEMNT(IJK),IJK=51,60)/  'Antimony    ',                    &       
     &                               'Tellurium   ',                    &       
     &                               'Iodine      ',                    &       
     &                               'Xenon       ',                    &       
     &                               'Cesium      ',                    &       
     &                               'Barium      ',                    &       
     &                               'Lanthanum   ',                    &       
     &                               'Cerium      ',                    &       
     &                               'Praseodymium',                    &       
     &                               'Neodymium   '/
      DATA (ELEMNT(IJK),IJK=61,70)/  'Promethium  ',                    &       
     &                               'Samarium    ',                    &       
     &                               'Europium    ',                    &       
     &                               'Gadolinium  ',                    &       
     &                               'Terbium     ',                    &       
     &                               'Dysprosium  ',                    &       
     &                               'Holmium     ',                    &       
     &                               'Erbium      ',                    &       
     &                               'Thulium     ',                    &       
     &                               'Ytterbium   '/
      DATA (ELEMNT(IJK),IJK=71,80)/  'Lutetium    ',                    &       
     &                               'Hafnium     ',                    &       
     &                               'Tantalum    ',                    &       
     &                               'Tungsten    ',                    &       
     &                               'Rhenium     ',                    &       
     &                               'Osmium      ',                    &       
     &                               'Iridium     ',                    &       
     &                               'Platinum    ',                    &       
     &                               'Gold        ',                    &       
     &                               'Mercury     '/
      DATA (ELEMNT(IJK),IJK=81,90)/  'Thallium    ',                    &       
     &                               'Lead        ',                    &       
     &                               'Bismuth     ',                    &       
     &                               'Polonium    ',                    &       
     &                               'Astatine    ',                    &       
     &                               'Radon       ',                    &       
     &                               'Francium    ',                    &       
     &                               'Radium      ',                    &       
     &                               'Actinium    ',                    &       
     &                               'Thorium     '/
      DATA (ELEMNT(IJK),IJK=91,100)/ 'Protactinium',                    &       
     &                               'Uranium     ',                    &       
     &                               'Neptunium   ',                    &       
     &                               'Plutonium   ',                    &       
     &                               'Americium   ',                    &       
     &                               'Curium      ',                    &       
     &                               'Berkelium   ',                    &       
     &                               'Californium ',                    &       
     &                               'Einsteinium ',                    &       
     &                               'Fermium     '/
      DATA (ELEMNT(IJK),IJK=101,IELM)/'Mendelevium ',                   &       
     &                               'Nobelium    ',                    &       
     &                               'Lawrencium  ',                    &       
     &                               'Rutherfordium',                   &       
     &                               'Dubnium     ',                    &       
     &                               'Seaborgium  ',                    &       
     &                               'Bohrium     ',                    &       
     &                               'Hassium     ',                    &       
     &                               'Meitnerium  ',                    &       
     &                               'Darmstadtium',                    &       
     &                               'Roentgenium '/
!
!     INITIALIZE
!
      IERR=0
      BCD = ' '
!
!     DETERMINE IZ AND IA
!
      IZ = ZA/1000
      IA = ZA - 1000*IZ
!
!     DETERMINE IF IZ AND IA ARE IN THE ALLOWABLE RANGE.
!
      IF(IZ.LE.0.OR.IZ.GT.IELM) THEN
         IERR = 1
         GO TO 100
      ENDIF
      IF(IA.LT.0.OR.IA.GE.1000) THEN
         IERR = 1
         GO TO 100
      ENDIF
!
!     ZERO A CORRESPONDS TO NATURAL ELEMENTS
!
      IF(IA.EQ.0) THEN
!
!        CONSTRUCT IDENTIFICATION FOR NATURAL MATERIAL
!
         BCD = 'Natural '//ELEMNT(IZ)
         NC = LEN_TRIM(BCD)
      ELSE
!
!     BUILD ISOTOPE NAME
!
         WRITE(BCDX,'(2A,I3)')   ELEMNT(IZ),'-',IA
!
!     ELIMINATE IMBEDDED BLANKS
!
         NC = 0
         DO I=1,24
            CCHAR = BCDX(I:I)
            IF(CCHAR.NE.' ') THEN
               IC = ICHAR(CCHAR)
               IF(IC.GT.96.AND.IC.LT.123)   CCHAR = CHAR(IC-32)
               NC = NC + 1
               BCD(NC:NC) = CCHAR
            END IF
         END DO
      END IF
      NCBCD = NC
!
  100 RETURN
      END SUBROUTINE ZAID
!
!***********************************************************************
!
      SUBROUTINE FILE1R
!
!     ROUTINE TO PROCESS ALL BUT MF=451 OF FILE 1
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2,N,NE
      REAL(KIND=R8) :: C1,C2
      INTEGER(KIND=I4) :: LFC, NFC
!
   10 CALL OUT_STATUS
!
!     BRANCH TO CORRECT FORMAT
!
      SELECT CASE (MT)
!
!     TOTAL OR PROMPT NUBAR
!
         CASE (452,456)
            IF(L2H.EQ.1)    THEN
               CALL CANT(C1,C2,L1,L2,N1,N2)
            ELSE
               CALL CANT1(C1,C2,L1,L2,N1,N2)
            END IF
!
!     DELAYED NEUTRON YIELD DATA.
!
         CASE (455)
            IF(L1H.EQ.1) THEN
               CALL CANT2(C1,C2,L1,L2,N1,NE)
               DO N=1,NE
                  CALL CANT(C1,C2,L1,L2,N1,N2)
               END DO
            ELSE
               CALL CANT(C1,C2,L1,L2,N1,N2)
            END IF
            IF(L2H.EQ.1)    THEN
               CALL CANT(C1,C2,L1,L2,N1,N2)
            ELSE
               CALL CANT1(C1,C2,L1,L2,N1,N2)
            END IF
!
!     ENERGY RELEASE IN FISSION
!
         CASE (458)
            CALL CANT(C1,C2,L1,L2,N1,N2)
            NFC = N2H
            LFC = L2H
            IF (LFC.EQ.1) THEN ! Fission energy release as tables, just skipping over, not checking anything
              DO N=1,NFC
                CALL CANT1(C1,C2,L1,L2,N1,N2)
              END DO
            END IF
!
!     ENERGY RELEASE IN FISSION
!
         CASE (460)
            IF(L1H.EQ.1) THEN
              NE=N1H
              DO N=1,NE
                CALL CANT1(C1,C2,L1,L2,N1,N2)
              END DO
            ELSE
              CALL CANT(C1,C2,L1,L2,N1,N2)
            END IF
!
      END SELECT
!
!     PROCESS SEND CARD
!
      CALL CONT(C1,C2,L1,L2,N1,N2)
!
!     CHECK FOR END OF FILE 1
!
      CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      IF(MF.GT.0)   GO TO 10
!
  100 RETURN
      END SUBROUTINE FILE1R
!
!***********************************************************************
!
      SUBROUTINE FILE2
!
!     ROUTINE TO PROCESS FILE 2
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NIS,NLS,NJS,NSS,NCRE
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: LFW,NER,LRU,LRF,NRO,LBK,LPS
      INTEGER(KIND=I4) :: I,II,J,N,NCR,LIL,NL,NS,NJ
      REAL(KIND=R8) :: C1,C2
!
      CALL OUT_STATUS
!
!     SAVE ISOTOPE COUNT.
!
      NIS = N1H
!
!     LOOP OVER ISOTOPES
!
      DO I=1,NIS
         CALL CONT(C1,C2,L1,LFW,NER,N2)
!
!        LOOP OVER ENERGY RANGES.
!
         DO II=1,NER
            CALL CONT(C1,C2,LRU,LRF,NRO,N2)
!
!           BRANCH TO RESOLVED OR UNRESOLVED DATA REPRESENTATION
!
            SELECT CASE (LRU)
!
!              SCATTERING LENGTH ONLY
!
               CASE (0)
                  CALL CONT(C1,C2,L1,L2,N1,N2)
!
!              DATA IS FOR RESOLVED RANGE
!
               CASE (1)
                  IF(NRO.GT.0)  CALL CANT1(C1,C2,L1,L2,N1,N2)
                  CALL CONT(C1,C2,L1,L2,NLS,N2)
                  SELECT CASE (LRF)
                     CASE (1:3)   !  BREIT-WIGNER AND REICH MOORE
                        DO N=1,NLS
                           CALL CANT(C1,C2,L1,L2,N1,N2)
                        END DO
                     CASE (4)   !  ADLER-ADLER
                        CALL CANT(C1,C2,L1,L2,N1,N2)
                        DO N=1,NLS
                           CALL CONT(C1,C2,L1,L2,NJS,N2)
                           DO J=1,NJS
                              CALL CANT(C1,C2,L1,L2,N1,N2)
                           END DO
                        END DO
                     CASE (6)      ! HYBRID R-FUNCTION
                        CALL CONT(C1,C2,L1,L2,N1,NCRE)
                        CALL CONT(C1,C2,L1,L2,N1,N2)
                        CALL CANT(C1,C2,L1,L2,N1,N2)
                        IF(NCRE.GT.0)  THEN
                           DO NCR=1,NCRE
                              DO LIL=1,4
                                 CALL CANT1(C1,C2,L1,L2,N1,N2)
                              END DO
                           END DO
                        END IF
                        DO NL=1,NLS
                           CALL CONT(C1,C2,L1,L2,NSS,N2)
                           DO NS=1,NSS
                              CALL CONT(C1,C2,L1,L2,NJS,N2)
                              DO NJ=1,NJS
                                 CALL CANT(C1,C2,LBK,LPS,N1,N2)
                                 IF(LBK.NE.0)  THEN
                                    CALL CANT1(C1,C2,L1,L2,N1,N2)
                                    CALL CANT1(C1,C2,L1,L2,N1,N2)
                                 END IF
                                 IF(LPS.NE.0)  THEN
                                    CALL CANT1(C1,C2,L1,L2,N1,N2)
                                    CALL CANT1(C1,C2,L1,L2,N1,N2)
                                 END IF
                              END DO
                           END DO
                        END DO
                     CASE (7)      ! R-Matrix Limited
                        NJS=NLS
                        CALL CANT(C1,C2,L1,L2,N1,N2)
                        DO N=1,NJS
                            CALL CANT(C1,C2,L1,L2,N1,N2)
                            CALL CANT(C1,C2,L1,L2,N1,N2)
                        END DO
                  END SELECT
!
!              UNRESOLVED RESONANCE PARAMETERS
!
               CASE (2)
                  IF(NRO.GT.0)  CALL CANT1(C1,C2,L1,L2,N1,N2)
                  SELECT CASE (LRF)
!
!                  ONLY FISSION WIDTHS MAY BE ENERGY DEPENDENT
!
                     CASE (1)
                        IF(LFW.EQ.0)   THEN   ! NO FISSION WIDTHS
                           CALL CONT(C1,C2,L1,L2,NLS,N2)
                           DO N=1,NLS
                              CALL CANT(C1,C2,L1,L2,N1,N2)
                           END DO
                        ELSE   ! ENERGY DEPENDENT FISSION WIDTHS GIVEN
                           CALL CANT(C1,C2,L1,L2,N1,NLS)
                           DO N=1,NLS
                              CALL CONT(C1,C2,L1,L2,NJS,N2)
                              DO J=1,NJS
                                 CALL CANT(C1,C2,L1,L2,N1,N2)
                              END DO
                           END DO
                        END IF
!
!                   ALL PARAMETERS ENERGY DEPENDENT
!
                    CASE(2)
                       CALL CONT(C1,C2,L1,L2,NLS,N2)
                       DO N=1,NLS
                          CALL CONT(C1,C2,L1,L2,NJS,N2)
                          DO J=1,NJS
                             CALL CANT(C1,C2,L1,L2,N1,N2)
                          END DO
                       END DO
!
                  END SELECT
!
            END SELECT
         END DO
      END DO
!
!     PROCESS SEND AND FEND RECORDS
!
      CALL CONT(C1,C2,L1,L2,N1,N2)
      CALL CONT(C1,C2,L1,L2,N1,N2)
!
      RETURN
      END SUBROUTINE FILE2
!
!***********************************************************************
!
      SUBROUTINE FILE3
!
!     ROUTINE TO PROCESS FILE 3
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS IN FILE 3
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        PROCESS TAB1 RECORD
!
         CALL CANT1(C1,C2,L1,L2,N1,N2)
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE3
!
!***********************************************************************
!
      SUBROUTINE FILE4
!
!     ROUTINE TO PROCESS FILE 4 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LTT,LI,NE
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS IN FILE 4
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE DATA TYPE FLAG (LEGENDRE OR TABULATED).
!
         LTT = L2H
!
!        PROCESS THE TRANSFORMATION MATRIX
!
         IF(L1H.NE.0)  THEN
            CALL CANT(C1,C2,LI,L2,N1,N2)
         ELSE
            CALL CONT(C1,C2,LI,L2,N1,N2)
         END IF
!
!        TEST FOR ISOTROPIC DISTRIBUTION
!
         IF(LI.EQ.0) THEN
!
!           READ TAB2 FOR INTERPOLATION BETWEEN ENERGIES.
!
            CALL CANT2(C1,C2,L1,L2,N1,NE)
!
!           LEGENDRE REPRESENTATION
!
            IF(LTT.EQ.1)   THEN
               DO N=1,NE
                  CALL CANT(C1,C2,L1,L2,N1,N2)
               END DO
!
!           TABULAR REPRESENTATION
!
            ELSE IF(LTT.EQ.2) THEN
               DO N=1,NE
                  CALL CANT1(C1,C2,L1,L2,N1,N2)
               END DO
!
!           DUAL REPRESENTATION
!
            ELSE IF(LTT.EQ.3) THEN
               DO N=1,NE
                  CALL CANT(C1,C2,L1,L2,N1,N2)
               END DO
               CALL CANT2(C1,C2,L1,L2,N1,NE)
               DO N=1,NE
                  CALL CANT1(C1,C2,L1,L2,N1,N2)
               END DO
            END IF
         END IF
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE4
!
!***********************************************************************
!
      SUBROUTINE FILE5
!
!     ROUTINE TO PROCESS FILE FIVE DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NK,NE,LF
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: I,N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE LAW COUNT.
!
         NK = N1H
!
!        LOOP OVER LAWS.
!
         DO I=1,NK
            CALL CANT1(C1,C2,L1,LF,N1,N2)
!
!           BRANCH BASED ON TYPE OF LAW.
!
            SELECT CASE (LF)
!
!              TABULAR REPRESENTATION
!
               CASE (1)
                  CALL CANT2(C1,C2,L1,L2,N1,NE)
                  DO N=1,NE
                     CALL CANT1(C1,C2,L1,L2,N1,N2)
                  END DO
!
!             ALL OTHER VALID LAWS
!
               CASE(7,9,12)
                  CALL CANT1(C1,C2,L1,L2,N1,N2)
               CASE(5,11)
                  CALL CANT1(C1,C2,L1,L2,N1,N2)
                  CALL CANT1(C1,C2,L1,L2,N1,N2)
!
            END SELECT
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE5
!
!***********************************************************************
!
      SUBROUTINE FILE6
!
!     ROUTINE TO PROCESS FILE SIX DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NK,NE,LF,NMU
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: I,J,N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE LAW COUNT.
!
         NK = N1H
!
!        LOOP OVER LAWS.
!
         DO I=1,NK
            CALL CANT1(C1,C2,L1,LF,N1,N2)
!
!           BRANCH BASED ON TYPE OF LAW.
!
            SELECT CASE (LF)
!
!              TABULAR, DISCRETE TWO-BODY AND COULOMB ELASTIC LAWS
!
               CASE (1,2,5)
                  CALL CANT2(C1,C2,L1,L2,N1,NE)
                  DO N=1,NE
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
!
!              N-BODY PHASE SPACE
!
               CASE (6)
                  CALL CONT(C1,C2,L1,L2,N1,N2)
!
!              LABORATORY ENERGY ANGLE LAW
!
               CASE (7)
                  CALL CANT2(C1,C2,L1,L2,N1,NE)
                  DO N=1,NE
                     CALL CANT2(C1,C2,L1,L2,N1,NMU)
                     DO J=1,NMU
                        CALL CANT1(C1,C2,L1,L2,N1,N2)
                     END DO
                  END DO
!
            END SELECT
!
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE6
!
!***********************************************************************
!
      SUBROUTINE FILE7
!
!     ROUTINE TO PROCESS FILE SEVEN DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LTHR,LT,NB,NS
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: I,N,NN
      REAL(KIND=R4), DIMENSION(4) :: CFLAG
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE SCATTERING TYPE FLAG
!
         LTHR = L1H
         SELECT CASE (LTHR)
!
!        INCOHERENT INELASTIC SCATTERING
!
            CASE (0)
               CALL CANT(C1,C2,L1,L2,N1,N2)
               NS = N2
               DO NN=1,NS+1
                  CFLAG(NN) = YP(6*(NN-1)+1)
               END DO
               IF(CFLAG(1).NE.0.) THEN
                  CALL CANT2(C1,C2,L1,L2,N1,NB)
!
!                 LOOP OVER BETA VALUES.
!
                  DO I=1,NB
                     CALL CANT1(C1,C2,LT,L2,N1,N2)
                     IF(LT.GT.0)   THEN
                        DO N=1,LT
                           CALL CANT(C1,C2,L1,L2,N1,N2)
                        END DO
                     END IF
                  END DO
               END IF
!
!              PROCESS EFFECTIVE TEMPERATURE RECORD
!
               IF(NFOR.GE.6)   THEN
                  CALL CANT1(C1,C2,L1,L2,N1,N2)
                  IF(NS.GT.0)   THEN
                     DO NN=1,NS
                        IF(CFLAG(NN+1).EQ.0.)  THEN
                           CALL CANT1(C1,C2,L1,L2,N1,N2)
                        END IF
                     END DO
                  END IF
               END IF
!
!           COHERENT ELASTIC SCATTERING
!
            CASE (1)
               CALL CANT1(C1,C2,LT,L2,N1,N2)
               IF(LT.GT.0) THEN
                  DO N=1,LT
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
               END IF
!
!           INCOHERENT ELASTIC SCATTERING
!
            CASE (2)
               CALL CANT1(C1,C2,L1,L2,N1,N2)
         END SELECT
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE7
!
!***********************************************************************
!
      SUBROUTINE FILE8
!
!     ROUTINE TO PROCESS FILE EIGHT DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LEP1,NSP,LCON,LCOV,NER,NO,NS
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: N,NN
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        BRANCH ON SECTION FORMAT
!
         SELECT CASE (MT)
!
!           FISSION PRODUCT YIELDS
!
            CASE (454,459)
               LEP1 = L1H
               DO N=1,LEP1
                  CALL CANT(C1,C2,L1,L2,N1,N2)
               END DO
!
!           DECAY DATA
!
            CASE (457)
               CALL CANT(C1,C2,L1,L2,N1,N2)
               CALL CANT(C1,C2,L1,L2,N1,N2)
               NSP = N2H
               IF(NSP.GT.0)  THEN
                  DO N=1,NSP
                     CALL CANT(C1,C2,LCON,L2,N1,NER)
                     SELECT CASE (LCON)
                        CASE (0)
                           DO NN=1,NER
                              CALL CANT(C1,C2,L1,L2,N1,N2)
                           END DO
                        CASE (1)
                           CALL CANT1(C1,C2,L1,LCOV,N1,N2)
                           IF(LCOV.NE.0)   THEN
                              CALL CANT(C1,C2,L1,L2,N1,N2)
                           END IF
                        CASE (2)
                           DO NN=1,NER
                              CALL CANT(C1,C2,L1,L2,N1,N2)
                           END DO
                           CALL CANT1(C1,C2,L1,LCOV,N1,N2)
                           IF(LCOV.NE.0)   THEN
                              CALL CANT(C1,C2,L1,L2,N1,N2)
                           END IF
!
                     END SELECT
                  END DO
               END IF
!
!           REACTION PRODUCTS
!
            CASE DEFAULT
               NO = N2H
               NS = N1H
               DO N=1,NS
                  IF(NO.EQ.0)  THEN
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  ELSE
                     CALL CONT(C1,C2,L1,L2,N1,N2)
                  END IF
               END DO
!
         END SELECT
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE8
!
!***********************************************************************
!
      SUBROUTINE FILE9
!
!     ROUTINE TO PROCESS FILE NINE AND TEN DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NS
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        PROCESS ALL SUBSECTIONS
!
         NS = N1H
         DO N=1,NS
            CALL CANT1(C1,C2,L1,L2,N1,N2)
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE9
!
!***********************************************************************
!
      SUBROUTINE FILE12
!
!     ROUTINE TO PROCESS FILE TWELVE AND THIRTEEN DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LO,NK
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
         LO = L1H
!
!        MULTIPLICITIES OR CROSS SECTIONS
!
         IF(LO.EQ.1.OR.LO.EQ.0)   THEN
            NK = N1H
            IF(NK.GT.1)  CALL CANT1(C1,C2,L1,L2,N1,N2)
            DO N=1,NK
               CALL CANT1(C1,C2,L1,L2,N1,N2)
            END DO
!
!        TRANSITION PROBABILITIES
!
         ELSE
            CALL CANT(C1,C2,L1,L2,N1,N2)
         END IF
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE12
!
!***********************************************************************
!
      SUBROUTINE FILE14
!
!     ROUTINE TO PROCESS FILE FOURTEEN DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LI,LTT,NI,NK,NKL,NE
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: K,N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE CONTROL VARIABLES
!
         NI = N2H
         NK = N1H
         LTT = L2H
         LI = L1H
         IF(LI.EQ.0)    THEN
!
!           READ ISOTROPIC DISTRIBUTIONS
!
            IF(NI.GT.0)   THEN
               DO N=1,NI
                  CALL CONT(C1,C2,L1,L2,N1,N2)
               END DO
            END IF
!
!           ANY NON-ISOTROPIC RECORDS
!
            NKL = NK - NI
            IF(NKL.GT.0)  THEN
!
!              LEGENDRE COEFFICIENTS OR TABULATED
!
               IF(LTT.EQ.1)    THEN
!
!                 LEGENDRE COEFFICIENTS
!
                  DO N=1,NKL
                     CALL CANT2(C1,C2,L1,L2,N1,NE)
                     DO K=1,NE
                        CALL CANT(C1,C2,L1,L2,N1,N2)
                     END DO
                  END DO
!
!                 TABULATED
!
               ELSE
                  DO N=1,NKL
                     CALL CANT2(C1,C2,L1,L2,N1,NE)
                     DO K=1,NE
                        CALL CANT1(C1,C2,L1,L2,N1,N2)
                     END DO
                  END DO
               END IF
            END IF
         END IF
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE14
!
!***********************************************************************
!
      SUBROUTINE FILE15
!
!     ROUTINE TO PROCESS FILE FIFTEEN DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NC,NE
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: N,NN
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE SUBSECTION COUNT AND SET UP LOOP OVER SUBSECTIONS
!
         NC = N1H
         DO NN=1,NC
            CALL CANT1(C1,C2,L1,L2,N1,N2)
            CALL CANT2(C1,C2,L1,L2,N1,NE)
            DO N=1,NE
               CALL CANT1(C1,C2,L1,L2,N1,N2)
            END DO
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE15
!
!***********************************************************************
!
      SUBROUTINE FILE23
!
!     ROUTINE TO PROCESS FILE 23 AND 27 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
         CALL CANT1(C1,C2,L1,L2,N1,N2)
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE23
!
!***********************************************************************
!
      SUBROUTINE FILE26
!
!     ROUTINE TO PROCESS FILE 26
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NK,LAW,NE
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: I,N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE LAW COUNT.
!
         NK = N1H
!
!        LOOP OVER LAWS.
!
         DO I=1,NK
            CALL CANT1(C1,C2,L1,LAW,N1,N2)
            SELECT CASE (LAW)
!
!              TABULAR AND DISCRETE TWO-BODY ELASTIC LAWS
!
               CASE (1,2)
                  CALL CANT2(C1,C2,L1,L2,N1,NE)
                  DO N=1,NE
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
!
!              ENERGY TRANSFER FOR EXCITATION
!
               CASE (8)
                  CALL CANT1(C1,C2,L1,L2,N1,N2)
!
            END SELECT
         END DO
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE26
!
!***********************************************************************
!
      SUBROUTINE FILE28
!
!     ROUTINE TO PROCESS FILE 28
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NSS
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: I
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        SAVE SUBSHELL COUNT.
!
         NSS = N1H
!
!        LOOP OVER SUBSHELLS
!
         DO I=1,NSS
            CALL CANT(C1,C2,L1,L2,N1,N2)
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE28
!
!***********************************************************************
!
      SUBROUTINE FILE32
!
!     ROUTINE TO PROCESS FILE 32
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: LFW,NER,LRU,LRF,NRO,LCOMP,ISR
      INTEGER(KIND=I4) :: NIS,NIT,NLS,NJS,NLRS,NSRS
      INTEGER(KIND=I4) :: NRSA12,NRSA,NM,NNN,NJSX,I3
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: I,II,I2,N
      REAL(KIND=R8) :: C1,C2
!
      CALL OUT_STATUS
!
!     LOOP OVER ISOTOPES
!
      NIS = N1H
         DO I=1,NIS
         CALL CONT(C1,C2,L1,LFW,NER,N2)
!
!        LOOP OVER ENERGY RANGES.
!
         DO II=1,NER
            CALL CONT(C1,C2,LRU,LRF,NRO,N2)
!
!           BRANCH TO RESOLVED OR UNRESOLVED DATA REPRESENTATION
!
            IF(LRU.EQ.1) THEN
!
!              DATA IS FOR RESOLVED RANGE
!
!**************READ ENERGY DEPENDENT SCATTERING LENGTH
               IF(NRO.GT.0)  THEN
                  CALL CONT(C1,C2,L1,L2,N1,NIT)
                  DO I2=1,NIT
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
               END IF
               CALL CONT(C1,C2,L1,LCOMP,NLS,ISR)
!
!              SINGLE AND MULTILEVEL BREIT-WIGNER ENDF-5 STYLE
!
               IF(LCOMP.EQ.0)   THEN
!                 Process scattering radius uncertainty for LCOMP=0
                  IF(ISR.EQ.1) CALL CONT(C1,C2,L1,L2,N1,N2)
                  DO N=1,NLS
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
!
!              NEW STYLE FORMATS
!
               ELSE IF(LCOMP.EQ.1) THEN
!                 Process scattering radius uncertainty for LCOMP=1
                  IF(ISR.EQ.1) THEN
                    IF(LRF.EQ.1 .OR. LRF.EQ.2) THEN
                      CALL CONT(C1,C2,L1,L2,N1,N2)
                    ELSE IF(LRF.EQ.3 .OR. LRF.EQ.7) THEN
                      CALL CANT(C1,C2,L1,L2,N1,N2)
                    END IF
                  END IF
                  CALL CONT(C1,C2,L1,L2,NSRS,NLRS)
                  IF(NSRS.GT.0)  THEN
                     DO I2=1,NSRS
                        IF(LRF.EQ.7) THEN
                          CALL CONT(C1,C2,NJSX,L2,N1,N2)
                          DO I3=1,NJSX
                            CALL CANT(C1,C2,L1,L2,N1,N2)
                          END DO
                        END IF
                        CALL CANT(C1,C2,L1,L2,N1,N2)
                     END DO
                  END IF
                  IF(NLRS.GT.0)  THEN
                     DO I2=1,NLRS
                        CALL CANT(C1,C2,L1,L2,N1,N2)
                     END DO
                  END IF
               ELSE IF(LCOMP.EQ.2) THEN
!                 Process scattering radius uncertainty for LCOMP=2
                  IF(ISR.EQ.1) THEN
                    IF(LRF.EQ.1 .OR. LRF.EQ.2) THEN
                      CALL CONT(C1,C2,L1,L2,N1,N2)
                    ELSE IF(LRF.EQ.3 .OR. LRF.EQ.7) THEN
                      CALL CANT(C1,C2,L1,L2,N1,N2)
                    END IF
                  END IF
                  IF((LRF.GE.1 .AND. LRF.LE.3) .OR. LRF.EQ.7) THEN
                     CALL CANT(C1,C2,L1, L2,NRSA12,NRSA)
                     IF(LRF.EQ.7) THEN
                        DO I2=1,NLS
                           CALL CANT(C1,C2,L1,L2,N1,N2)
                           CALL CANT(C1,C2,L1,L2,N1,N2)
                        END DO
                     END IF
                     CALL CONT(C1,C2,L1,NNN,    NM,  N2)
                     DO I2=1,NM
                        CALL TEXTR(1)
                     END DO
                  ELSE
                     WRITE(*,*) 'Illegal LRF',LRF
                     STOP 'STANEF ERROR IN FILE32 - Unsupported LRF'
                  END IF
               ELSE
                  STOP 'STANEF ERROR IN FILE32 - Illegal LCOMP'
               END IF
!
!           UNRESOLVED RESONANCE PARAMETERS
!
            ELSE
               CALL CONT(C1,C2,L1,L2,NLS,N2)
               DO N=1,NLS
!original start
!                 CALL CONT(C1,C2,L1,L2,NJS,N2)
!                 DO NN=1,NJS
!                    CALL CANT(C1,C2,L1,L2,N1,N2)
!                 END DO
!original end
!trkov    start
                  CALL CANT(C1,C2,L1,L2,N1,NJS)
!trkov    end
               END DO
               CALL CANT(C1,C2,L1,L2,N1,N2)
            END IF
         END DO
      END DO
!
!     PROCESS SEND AND FEND RECORDS
!
      CALL CONT(C1,C2,L1,L2,N1,N2)
      CALL CONT(C1,C2,L1,L2,N1,N2)
!
      RETURN
      END SUBROUTINE FILE32
!
!***********************************************************************
!
      SUBROUTINE FILE33
!
!     ROUTINE TO PROCESS FILE 33 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: MTL,NL,NC,NI
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: I,J,N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        CHECK FOR LUMPED COVARIANCES
!
         MTL = L2H
         IF(MTL.EQ.0) THEN
!
!           LOOP OVER SUBSECTIONS
!
            NL = N2H
            DO N=1,NL
               CALL CONT(C1,C2,L1,L2,NC,NI)
!
!              LOOP OVER NC-TYPE SUB-SUBSECT.
!
               IF(NC.GT.0) THEN
                  DO I=1,NC
                     CALL CONT(C1,C2,L1,L2,N1,N2)
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
               END IF
!
!              LOOP OVER NI-TYPE SUB-SUBSECT.
!
               IF(NI.GT.0) THEN
                  DO J=1,NI
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
               END IF
            END DO
         END IF
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE33
!
!***********************************************************************
!
      SUBROUTINE FILE34
!
!     ROUTINE TO PROCESS FILE 34 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NI,NMT1,NL,NL1,NSS
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: J,N,NN
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        LOOP OVER SUBSECTIONS
!
         NMT1 = N2H
         DO N=1,NMT1
            CALL CONT(C1,C2,L1,L2,NL,NL1)
            IF(MT.EQ.L2) THEN
               NSS = NL*(NL+1)/2
            ELSE
               NSS = NL*NL1
            END IF
!
!           LOOP OVER ALL SUB-SUBSECTIONS
!
            DO NN=1,NSS
               CALL CONT(C1,C2,L1,L2,N1,NI)
!
!              LOOP OVER NI-TYPE SUB-SUB-SUBSECT.
!
               DO J=1,NI
                  CALL CANT(C1,C2,L1,L2,N1,N2)
               END DO
            END DO
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE34
!
!***********************************************************************
!
      SUBROUTINE FILE35
!
!     ROUTINE TO PROCESS FILE 35 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NK
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: N
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        LOOP OVER SUBSECTIONS
!
         NK = N1H
         DO N=1,NK
            CALL CANT(C1,C2,L1,L2,N1,N2)
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE35
!
!***********************************************************************
!
      SUBROUTINE FILE40
!
!     ROUTINE TO PROCESS FILE 40 DATA
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: NS,NL,NC,NI
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      INTEGER(KIND=I4) :: J,N,NN
      REAL(KIND=R8) :: C1,C2
!
!     PROCESS ALL SECTIONS
!
      DO WHILE (MF.GT.0)
!
         CALL OUT_STATUS
!
!        LOOP OVER SUBSECTIONS
!
         NS = N1H
         DO N=1,NS
            CALL CONT(C1,C2,L1,L2,N1,NL)
!
!           LOOP OVER ALL SUB-SUBSECTIONS
!
            DO NN=1,NL
               CALL CONT(C1,C2,L1,L2,NC,NI)
!
!              LOOP OVER NC-TYPE SUB-SUB-SUBSECT.
!
               IF(NC.GT.0)   THEN
                  DO J=1,NC
                     CALL CONT(C1,C2,L1,L2,N1,N2)
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
               END IF
!
!              LOOP OVER NI-TYPE SUB-SUB-SUBSECT.
!
               IF(NI.GT.0)  THEN
                  DO J=1,NI
                     CALL CANT(C1,C2,L1,L2,N1,N2)
                  END DO
               END IF
            END DO
         END DO
!
!        PROCESS SEND CARD
!
         CALL CONT(C1,C2,L1,L2,N1,N2)
!
!        CHECK FOR END OF FILE
!
         CALL CONT(C1H,C2H,L1H,L2H,N1H,N2H)
      END DO
!
      RETURN
      END SUBROUTINE FILE40
!
!***********************************************************************
!
      SUBROUTINE TEXTR(IPATH)
!
!     SUBROUTINE TO READ AND WRITE A TEXT RECORD IN CHARACTER FORMAT
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: IPATH
!
      INTEGER(KIND=I4), INTRINSIC :: MIN0
!
!     READ THE TEXT RECORD
!
      IF(IPATH.EQ.1)   READ(ITAPE,'(A,I4,I2,I3,I5)') TEXT
!
!     WRITE THE TEXT RECORD
!
      IF(STANEF_DATA%MODE.EQ.1)  THEN
         WRITE(OTAPE) MAT,MF,MT,TEXT
      ELSE
         NSEQ = MIN0(NSEQ+1,99999)
         WRITE(OTAPE,'(A,I4,I2,I3,I5)') TEXT,MAT,MF,MT,NSEQ
      END IF
!
      RETURN
      END SUBROUTINE TEXTR
!
!***********************************************************************
!
      SUBROUTINE CONT(C1,C2,L1,L2,N1,N2)
!
!     READ AND WRITE A CONT RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      REAL(KIND=R8) :: C1,C2
!
      INTEGER(KIND=I4), INTRINSIC :: MIN0
!
      CHARACTER(LEN=11), DIMENSION(2) ::  FLOATS
      CHARACTER(LEN=11), PARAMETER :: RZERO=' 0.000000+0',              &       
     &                                IZERO='          0'
!
!     READ CONT RECORD
!
      READ(ITAPE,'(2E11.4,4I11,I4,I2,I3)')                              &       
     &                    C1,C2,L1,L2,N1,N2,MAT,MF,MT
!
!     OUTPUT RECORD
!
      IF(STANEF_DATA%MODE.EQ.1)   THEN
!
!        BINARY MODE
!
         WRITE(OTAPE) MAT,MF,MT,C1,C2,L1,L2,N1,N2
      ELSE
!
!        FORMAT THE FLOATING POINT NUMBERS AND OUTPUT RECORD
!
         IF(MT.NE.0)    THEN
            CALL NORMAL(C1,FLOATS(1))
            CALL NORMAL(C2,FLOATS(2))
            NSEQ = MIN0(NSEQ+1,99999)
            WRITE(OTAPE,'(2A,4I11,I4,I2,I3,I5)')                        &       
     &               FLOATS(1),FLOATS(2),L1,L2,N1,N2,MAT,MF,MT,NSEQ
!
!        A SEND, FEND, MEND OR TEND RECORD
!
         ELSE
            IF(NFOR.LT.6) THEN
              NSEQ = MIN0(NSEQ+1,99999)
              IF(MAT.LT.0)   NSEQ = 0
            ELSE
               IF(MAT.LT.0) THEN
                  NSEQ = 0
               ELSE IF(MF.EQ.0) THEN
                  NSEQ = 0
               ELSE IF(MT.EQ.0) THEN
                  NSEQ = 99999
               ELSE
                  NSEQ = NSEQ + 1
               END IF
            END IF
            WRITE(OTAPE,'(6A,I4,I2,I3,I5)')                              &      
     &          RZERO,RZERO,IZERO,IZERO,IZERO,IZERO,MAT,MF,MT,NSEQ
            IF(NFOR.GT.5.AND.MT.EQ.0)  NSEQ = 0
         END IF
      END IF
!
      RETURN
      END SUBROUTINE CONT
!
!***********************************************************************
!
      SUBROUTINE CANT(C1,C2,L1,L2,N1,N2)
!
!     READ AND WRITE A LIST RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2,I1,I2
      REAL(KIND=R8) :: C1,C2
!
      INTEGER(KIND=I4), INTRINSIC :: MIN0
!
      CHARACTER(LEN=66) :: RECORD
      CHARACTER(LEN=66), DIMENSION(6) :: FLOATS
      INTEGER(KIND=I4) :: NCS,NU,MM
      INTEGER(KIND=I4) :: M,N
!
!     PROCESS THE CONT-LIKE PORTION OF RECORD
!
      IF(STANEF_DATA%MODE.EQ.0)  THEN
         CALL CONT(C1,C2,L1,L2,N1,N2)
      ELSE
         READ(ITAPE,'(2E11.4,4I11)') C1,C2,L1,L2,N1,N2
      END IF
!
!     SPECIAL PROCESSING FOR FILE 2 WHEN DATA NOT TO BE REFORMATED
!
      IF(STANEF_DATA%I151.NE.'Y'.AND.MT.EQ.151                          &       
     &             .AND.STANEF_DATA%MODE.EQ.0) THEN
         NCS = (N1+5)/6
         DO N=1,NCS
            READ(ITAPE,'(A)')  RECORD
            NSEQ = MIN0(NSEQ+1,99999)
            WRITE(OTAPE,'(A,I4,I2,I3,I5)') RECORD,MAT,MF,MT,NSEQ
         END DO
!
!     READ IN LIST FORMATTED ARRAY
!
      ELSE
         I1= 1
         I2=MIN(6*(POINTSMAX/6),N1)
!
         DO WHILE(I1.LE.N1)
           READ(ITAPE,'(6E11.4)') (YP(N-I1+1),N=I1,I2)
!
!          OUTPUT THE RECORD
!
           IF(STANEF_DATA%MODE.EQ.1)   THEN
!
!             OUTPUT IN BINARY MODE
!
              IF(I2.NE.N1) THEN
                STOP 'STANEF ERRORR - LIST too big for binary output'
              END IF
!
              WRITE(OTAPE) MAT,MF,MT,C1,C2,L1,L2,N1,N2,(YP(N),N=1,N1)
           ELSE
              MM = 1
              DO N=I1,I2,6
                 NU = MIN0(6,I2+1-N)
                 DO M=1,6
                    IF(M.LE.NU)   THEN
                       CALL NORMAL(YP(MM),FLOATS(M))
                    ELSE
                       FLOATS(M) = ' '
                    END IF
                    MM = MM + 1
                 END DO
!
!                OUTPUT NEXT SIX FIELD RECORD
!
                 NSEQ = NSEQ+1
                 IF(NSEQ.GT.99999) NSEQ=0
                 WRITE(OTAPE,'(6A11,I4,I2,I3,I5)') FLOATS,MAT,MF,MT,NSEQ
              END DO
           END IF
           I1=I2+1
           I2=MIN(I1-1+6*(POINTSMAX/6),N1)
         END DO
      END IF
!
      RETURN
      END SUBROUTINE CANT
!
!***********************************************************************
!
      SUBROUTINE CANT1(C1,C2,L1,L2,N1,N2)
!
!     READ AND WRITE A TAB1 RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      REAL(KIND=R8) :: C1,C2
!
      INTEGER(KIND=I4), INTRINSIC :: MIN0
!
      CHARACTER(LEN=66), DIMENSION(6) :: FLOATS
      INTEGER(KIND=I4)  :: NP,NU,MM,NN
      INTEGER(KIND=I4) :: M,N
!
!     PROCESS A TAB2-LIKE PORTION OF RECORD
!
      IF(STANEF_DATA%MODE.EQ.0) THEN
         CALL CANT2(C1,C2,L1,L2,N1,N2)
      ELSE
         READ(ITAPE,'(2E11.4,4I11)') C1,C2,L1,L2,N1,N2
         READ(ITAPE,'(6I11)') (NBT(N),INTRP(N),N=1,N1)
      END IF
!
!     READ IN DATA ARRAY
!
      NP = N2
      READ(ITAPE,'(6E11.4)')  (XP(N),YP(N),N=1,NP)
!
!     OUTPUT RECORD
!
      IF(STANEF_DATA%MODE.EQ.1)  THEN
!
!        BINARY MODE
!
         WRITE(OTAPE) MAT,MF,MT,C1,C2,L1,L2,N1,N2,                      &       
     &           (NBT(N),INTRP(N),N=1,N1),(XP(N),YP(N),N=1,NP)
!
!        FORMAT THE FLOATING POINT NUMBERS
!
      ELSE
         DO N=1,NP,3
            NU = MIN0(3,NP+1-N)
            MM = 1
            DO M=1,3
              NN = N + M -1
               IF(M.GT.NU)   THEN
                  FLOATS(MM) = ' '
                  FLOATS(MM+1) = ' '
               ELSE
                  CALL NORMAL(XP(NN),FLOATS(MM))
                  CALL NORMAL(YP(NN),FLOATS(MM+1))
               END IF
               MM = MM + 2
            END DO
!
!           OUTPUT NEXT SIX FIELD RECORD
!
            NSEQ = MIN0(NSEQ+1,99999)
            WRITE(OTAPE,'(6A11,I4,I2,I3,I5)') FLOATS,MAT,MF,MT,NSEQ
         END DO
      END IF
!
      RETURN
      END SUBROUTINE CANT1
!
!***********************************************************************
!
      SUBROUTINE CANT2(C1,C2,L1,L2,N1,N2)
!
!     READ AND WRITE A TAB2 RECORD
!
      IMPLICIT NONE
!
      INTEGER(KIND=I4) :: L1,L2,N1,N2
      REAL(KIND=R8) :: C1,C2
!
      INTEGER(KIND=I4), INTRINSIC :: MIN0
!
      CHARACTER(LEN=66), DIMENSION(6) :: FIXS
      INTEGER(KIND=I4)  :: N
      INTEGER(KIND=I4) :: I,J,K
!
!     PROCESS THE CONT-LIKE PORTION OF RECORD
!
      IF(STANEF_DATA%MODE.EQ.0)  THEN
         CALL CONT(C1,C2,L1,L2,N1,N2)
      ELSE
         READ(ITAPE,'(2E11.4,4I11)') C1,C2,L1,L2,N1,N2
      END IF
!
!     READ INTERPOLATION TABLE
!
      READ(ITAPE,'(6I11)')  (NBT(I),INTRP(I),I=1,N1)
!
!     WRITE INTERPOLATION TABLE
!
      IF(STANEF_DATA%MODE.EQ.1)   THEN
!
!        OUTPUT IN BINARY MODE
!
         WRITE(OTAPE) MAT,MF,MT,C1,C2,L1,L2,N1,N2,                      &       
     &                (NBT(N),INTRP(N),N=1,N1)
      ELSE
!
!        OUTPUT IN ASCII MODE
!
         DO N=1,N1,3
            NSEQ = MIN0(NSEQ+1,99999)
            IF(N+2.LE.N1)   THEN
               WRITE(OTAPE,'(6I11,I4,I2,I3,I5)')                        &       
     &                (NBT(I),INTRP(I),I=N,N+2),MAT,MF,MT,NSEQ
            ELSE
               J =  0
               DO K=N,N1
                  J = J + 1
                  WRITE(FIXS(J),'(I11)') NBT(K)
                  J = J + 1
                  WRITE(FIXS(J),'(I11)') INTRP(K)
               END DO
               DO K=J+1,6
                  FIXS(K) = ' '
               END DO
               WRITE(OTAPE,'(6A11,I4,I2,I3,I5)')  FIXS,MAT,MF,MT,NSEQ
            END IF
         END DO
      END IF
!
      RETURN
      END SUBROUTINE CANT2
!
!***********************************************************************
!
      SUBROUTINE NORMAL(X,CX)
!
!     ROUTINE TO CONVERT A FLOATING POINT NUMBER TO A CHARACTER STRING
!        IN ENDF FORMAT
!
      IMPLICIT NONE
!
      CHARACTER(LEN=11) :: CX
      REAL(KIND=R8) :: X
!
      INTEGER(KIND=I4), INTRINSIC :: IFIX, IABS
      REAL(KIND=R4), INTRINSIC :: ABS
!
      CHARACTER(LEN=1) :: CSIGN
      INTEGER(KIND=I4) :: IPOW
      REAL(KIND=R4) :: POWER
      REAL(KIND=R8) :: FNUM,ZERO
      REAL(KIND=4) :: SPOWER
!
!     INITIALIZE
!
      ZERO=0
      IF(X.NE.ZERO)   THEN
!
!        FIND POWER OF NUMBER
!
         POWER = LOG10(ABS(X)) + .000002
         SPOWER = POWER
         IPOW = IFIX(SPOWER)
         IF(POWER.LT.0.0)   THEN
            CSIGN = '-'
            IPOW = IPOW - 1
         ELSE
            CSIGN = '+'
         END IF
!
!        FIND NORMALIZED MANTISSA
!
         FNUM = X*(10.0D00**(-IPOW))
         IPOW = IABS(IPOW)
!
!        CHECK THE EXPONENT
!
         IF(IPOW.GE.1000)   THEN
            STOP 'STANEF ERROR - Exponent > 999'
!
!        THREE DIGIT EXPONENT
!
         ELSE IF(IPOW.GE.100)   THEN
            WRITE(CX,'(F7.4,A1,I3)')  FNUM,CSIGN,IPOW
!
!        TWO DIGIT EXPONENT
!
         ELSE IF(IPOW.GE.10)   THEN
            WRITE(CX,'(F8.5,A1,I2)')  FNUM,CSIGN,IPOW
!
!        ONE DIGIT EXPONENT
!
         ELSE
            WRITE(CX,'(F9.6,A1,I1)')  FNUM,CSIGN,IPOW
         END IF
      ELSE
         CX = ' 0.000000+0'
      END IF
!
  100 RETURN
      END SUBROUTINE NORMAL
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
      INTEGER(KIND=I4), INTRINSIC :: LEN_TRIM
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
      SUBROUTINE UPSTR(String)
!
!     ROUTINE TO CONVERT A STRING TO ALL UPPER CASE
!
!     STRING  -  CHARACTER STRING TO BE CONVERTED
!
      IMPLICIT NONE
!
!     Dummy arguments
!
      CHARACTER(LEN=*) :: String
!
!     Local variables
!
      CHARACTER(LEN=1),INTRINSIC :: CHAR
      INTEGER(KIND=4) :: i,ic,l
      INTEGER(KIND=4),INTRINSIC :: ICHAR
      INTEGER(KIND=4),INTRINSIC :: LEN_TRIM
!
      l = LEN_TRIM(String)
      DO i = 1,l
         ic = ICHAR(String(i:i))
         IF(ic.GT.96.AND.ic.LT.123)String(i:i) = CHAR(ic-32)
      END DO
!
      RETURN
      END SUBROUTINE UPSTR
!
!***********************************************************************
!
      SUBROUTINE GET_FROM_CLINE
!
!     GET CONTENTS OF COMMAND LINE OR NAME OF BATCH INPUT FILE
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
      INPUT = INPUT0
!
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
!/      WRITE(OUTPUT,'(A)')                                             &       
!/     &    ' Control File Specification        - '
!/      READ(INPUT,'(A)') CFILE
!/      INPUT = 19
!/      OPEN(UNIT=INPUT,FILE=CFILE,STATUS='OLD')
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
         WRITE(OUTPUT,'(5X,A,I5,A,I3,A,I4)')                            &       
     &         'PROCESSING MAT=',MAT,', MF=',MF,', MT=',MT
!---MDC---
      END IF
!
      RETURN
      END SUBROUTINE OUT_STATUS
!
!+++MDC+++
!...VMS, ANS, WIN, UNX
      END PROGRAM STANEF
!...LWI, DVF, MOD
!/      END MODULE STANEF
!---MDC---
