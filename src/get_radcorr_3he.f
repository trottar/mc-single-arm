      SUBROUTINE INIT_RADCORR_3HE(IMOD0,EBEAM_GEV,STAT)
C----------------------------------------------------------------------
C Load the selected 3He radiative-weight tables once before event
C generation.  The table ratio is XSrad_unp / XSborn_unp: it converts
C the MC's Born F1/F2 weight into a radiated weight.
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NRAD_ANG,NRAD_NU_MAX
      PARAMETER (NRAD_ANG=15,NRAD_NU_MAX=5200)
      INTEGER IMOD0,IA,IOERR,UNITNO,NVALID,K
      INTEGER RAD_NNU(NRAD_ANG),RAD_LOADED_IMOD
      REAL*8 EBEAM_GEV,RAD_ANGLE(NRAD_ANG)
      REAL*8 RAD_NU(NRAD_NU_MAX,NRAD_ANG)
      REAL*8 RAD_RATIO(NRAD_NU_MAX,NRAD_ANG)
      REAL*8 RAD_TABLE_EBEAM,COL(17),PREV_NU
      REAL*8 RAD_EBEAM_GEV,RAD_EBEAM_MEV,RAD_EBEAM_TOL_MEV
      CHARACTER*200 FILENAME
      CHARACTER*512 LINE
      LOGICAL STAT,RAD_INITIALIZED,VALID,SAW_VALID,SAW_INVALID
      COMMON /RADCORR_3HE_CACHE/ RAD_ANGLE,RAD_NU,RAD_RATIO,
     > RAD_TABLE_EBEAM,RAD_NNU,RAD_LOADED_IMOD,RAD_INITIALIZED
      PARAMETER (RAD_EBEAM_GEV=10.38D0,RAD_EBEAM_MEV=10380.D0,
     >           RAD_EBEAM_TOL_MEV=0.1D0)

      STAT=.FALSE.
      IF (IMOD0.LT.1 .OR. IMOD0.GT.5) THEN
         WRITE(6,*) 'ERROR: invalid radcorr SF model =',IMOD0
         RETURN
      ENDIF
      IF (DABS(EBEAM_GEV-RAD_EBEAM_GEV).GT.1.D-3) THEN
         WRITE(6,*) 'ERROR: 3He radcorr tables require Ebeam (GeV)=',
     >              RAD_EBEAM_GEV
         WRITE(6,*) '       requested Ebeam (GeV)=',EBEAM_GEV
         RETURN
      ENDIF

C A repeat initialization for the same model is a no-op.  A different
C model is allowed for standalone validation programs and reloads its
C selected 15 tables without involving the event lookup routine.
      IF (RAD_INITIALIZED .AND. RAD_LOADED_IMOD.EQ.IMOD0) THEN
         STAT=.TRUE.
         RETURN
      ENDIF
      RAD_INITIALIZED=.FALSE.
      RAD_LOADED_IMOD=0
      RAD_TABLE_EBEAM=0.D0
      DO IA=1,NRAD_ANG
         RAD_ANGLE(IA)=26.D0+0.5D0*DBLE(IA-1)
         RAD_NNU(IA)=0
      ENDDO

      DO IA=1,NRAD_ANG
         WRITE(FILENAME,100) IMOD0,RAD_ANGLE(IA)
  100    FORMAT('interp/radcorr_tables/Newfit_20260710_fullxquad_15',
     >   'angles/SF',I1,'_G1F1cmplt_QE95/radiated_data_',F4.1,
     >   'deg_short.csv')
         UNITNO=70+IA
         OPEN(UNITNO,FILE=FILENAME,STATUS='OLD',IOSTAT=IOERR)
         IF (IOERR.NE.0) THEN
            WRITE(6,*) 'ERROR: cannot open 3He radcorr table:'
            WRITE(6,*) FILENAME
            GOTO 900
         ENDIF

         READ(UNITNO,'(A)',IOSTAT=IOERR) LINE
         IF (IOERR.NE.0 .OR. INDEX(LINE,'XSborn_unp').EQ.0 .OR.
     >       INDEX(LINE,'XSrad_unp').EQ.0) THEN
            WRITE(6,*) 'ERROR: invalid header in 3He radcorr table:'
            WRITE(6,*) FILENAME
            CLOSE(UNITNO)
            GOTO 900
         ENDIF

         NVALID=0
         SAW_VALID=.FALSE.
         SAW_INVALID=.FALSE.
  200    CONTINUE
         READ(UNITNO,'(A)',IOSTAT=IOERR) LINE
         IF (IOERR.LT.0) GOTO 300
         IF (IOERR.NE.0) THEN
            WRITE(6,*) 'ERROR: failed reading 3He radcorr table:'
            WRITE(6,*) FILENAME
            CLOSE(UNITNO)
            GOTO 900
         ENDIF
         READ(LINE,*,IOSTAT=IOERR) (COL(K),K=1,17)
         IF (IOERR.NE.0) THEN
            WRITE(6,*) 'ERROR: invalid data row in 3He radcorr table:'
            WRITE(6,*) FILENAME
            CLOSE(UNITNO)
            GOTO 900
         ENDIF
         IF (DABS(COL(1)-RAD_EBEAM_MEV).GT.RAD_EBEAM_TOL_MEV) THEN
            WRITE(6,*) 'ERROR: unexpected table Ebeam (MeV)=',COL(1)
            WRITE(6,*) FILENAME
            CLOSE(UNITNO)
            GOTO 900
         ENDIF

         VALID=.TRUE.
         IF (COL(5).LE.0.D0 .OR. COL(10).LE.0.D0) VALID=.FALSE.
         IF (COL(5).NE.COL(5) .OR. COL(10).NE.COL(10)) VALID=.FALSE.
         IF (DABS(COL(5)).GT.1.D300 .OR.
     >       DABS(COL(10)).GT.1.D300) VALID=.FALSE.
         IF (VALID) THEN
            IF (SAW_INVALID) THEN
               WRITE(6,*) 'ERROR: non-contiguous finite radcorr grid:'
               WRITE(6,*) FILENAME
               CLOSE(UNITNO)
               GOTO 900
            ENDIF
            NVALID=NVALID+1
            IF (NVALID.GT.NRAD_NU_MAX) THEN
               WRITE(6,*) 'ERROR: radcorr table exceeds NRAD_NU_MAX='
     >                    ,NRAD_NU_MAX
               CLOSE(UNITNO)
               GOTO 900
            ENDIF
            IF (NVALID.GT.1 .AND. COL(2).LE.PREV_NU) THEN
               WRITE(6,*) 'ERROR: non-monotonic nu grid in:'
               WRITE(6,*) FILENAME
               CLOSE(UNITNO)
               GOTO 900
            ENDIF
            RAD_NU(NVALID,IA)=COL(2)
            RAD_RATIO(NVALID,IA)=COL(10)/COL(5)
            IF (RAD_RATIO(NVALID,IA).NE.RAD_RATIO(NVALID,IA) .OR.
     >          RAD_RATIO(NVALID,IA).LE.0.D0 .OR.
     >          DABS(RAD_RATIO(NVALID,IA)).GT.1.D300) THEN
               WRITE(6,*) 'ERROR: invalid rad-weight ratio in:'
               WRITE(6,*) FILENAME
               CLOSE(UNITNO)
               GOTO 900
            ENDIF
            PREV_NU=COL(2)
            SAW_VALID=.TRUE.
         ELSE IF (SAW_VALID) THEN
            SAW_INVALID=.TRUE.
         ENDIF
         GOTO 200

  300    CONTINUE
         CLOSE(UNITNO)
         IF (NVALID.LT.2) THEN
            WRITE(6,*) 'ERROR: too few finite rows in radcorr table:'
            WRITE(6,*) FILENAME
            GOTO 900
         ENDIF
         RAD_NNU(IA)=NVALID
      ENDDO

      RAD_TABLE_EBEAM=RAD_EBEAM_MEV
      RAD_LOADED_IMOD=IMOD0
      RAD_INITIALIZED=.TRUE.
      STAT=.TRUE.
      WRITE(6,*) 'Loaded 3He radiative-weight tables for SF',IMOD0
      WRITE(6,*) 'Table Ebeam (MeV)=',RAD_TABLE_EBEAM
      RETURN

  900 CONTINUE
      RAD_INITIALIZED=.FALSE.
      RAD_LOADED_IMOD=0
      RAD_TABLE_EBEAM=0.D0
      STAT=.FALSE.
      RETURN
      END


      SUBROUTINE GET_RADCORR_3HE(IMOD0,THETA_DEG,NU_GEV,
     >                           RAD_WEIGHT_FACTOR,STAT)
C----------------------------------------------------------------------
C Return XSrad_unp / XSborn_unp from the already initialized selected
C model.  This routine performs no file I/O.
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NRAD_ANG,NRAD_NU_MAX
      PARAMETER (NRAD_ANG=15,NRAD_NU_MAX=5200)
      INTEGER IMOD0,IA,IA1,IA2
      INTEGER RAD_NNU(NRAD_ANG),RAD_LOADED_IMOD
      REAL*8 THETA_DEG,NU_GEV,RAD_WEIGHT_FACTOR
      REAL*8 RAD_ANGLE(NRAD_ANG),RAD_NU(NRAD_NU_MAX,NRAD_ANG)
      REAL*8 RAD_RATIO(NRAD_NU_MAX,NRAD_ANG),RAD_TABLE_EBEAM
      REAL*8 THETA_LOOKUP,NU_TABLE,R1,R2,FRAC
      REAL*8 ANG_TOL,NU_TOL
      LOGICAL STAT,RAD_INITIALIZED,STAT1,STAT2
      COMMON /RADCORR_3HE_CACHE/ RAD_ANGLE,RAD_NU,RAD_RATIO,
     > RAD_TABLE_EBEAM,RAD_NNU,RAD_LOADED_IMOD,RAD_INITIALIZED
      PARAMETER (ANG_TOL=1.D-9,NU_TOL=1.D-6)

      STAT=.FALSE.
      RAD_WEIGHT_FACTOR=1.D0
      IF (.NOT.RAD_INITIALIZED .OR. IMOD0.NE.RAD_LOADED_IMOD) THEN
         RETURN
      ENDIF

      THETA_LOOKUP=THETA_DEG
      IF (THETA_LOOKUP.LT.RAD_ANGLE(1)-ANG_TOL .OR.
     >    THETA_LOOKUP.GT.RAD_ANGLE(NRAD_ANG)+ANG_TOL) RETURN
      IF (DABS(THETA_LOOKUP-RAD_ANGLE(1)).LE.ANG_TOL) THEN
         IA1=1
         IA2=1
         THETA_LOOKUP=RAD_ANGLE(1)
      ELSE IF (DABS(THETA_LOOKUP-RAD_ANGLE(NRAD_ANG)).LE.ANG_TOL) THEN
         IA1=NRAD_ANG
         IA2=NRAD_ANG
         THETA_LOOKUP=RAD_ANGLE(NRAD_ANG)
      ELSE
         IA1=0
         IA2=0
         DO IA=1,NRAD_ANG
            IF (DABS(THETA_LOOKUP-RAD_ANGLE(IA)).LE.ANG_TOL) THEN
               IA1=IA
               IA2=IA
               THETA_LOOKUP=RAD_ANGLE(IA)
               GOTO 400
            ENDIF
         ENDDO
         DO IA=1,NRAD_ANG-1
            IF (THETA_LOOKUP.GT.RAD_ANGLE(IA) .AND.
     >          THETA_LOOKUP.LT.RAD_ANGLE(IA+1)) THEN
               IA1=IA
               IA2=IA+1
               GOTO 400
            ENDIF
         ENDDO
      ENDIF

  400 CONTINUE
      IF (IA1.EQ.0 .OR. IA2.EQ.0) RETURN
      NU_TABLE=1000.D0*NU_GEV
      CALL RADCORR_NU_LOOKUP(IA1,NU_TABLE,R1,STAT1)
      IF (.NOT.STAT1) RETURN
      IF (IA2.EQ.IA1) THEN
         RAD_WEIGHT_FACTOR=R1
         STAT=.TRUE.
         RETURN
      ENDIF
      CALL RADCORR_NU_LOOKUP(IA2,NU_TABLE,R2,STAT2)
      IF (.NOT.STAT2) RETURN
      FRAC=(THETA_LOOKUP-RAD_ANGLE(IA1)) /
     >     (RAD_ANGLE(IA2)-RAD_ANGLE(IA1))
      RAD_WEIGHT_FACTOR=R1+FRAC*(R2-R1)
      IF (RAD_WEIGHT_FACTOR.NE.RAD_WEIGHT_FACTOR .OR.
     >    RAD_WEIGHT_FACTOR.LE.0.D0 .OR.
     >    DABS(RAD_WEIGHT_FACTOR).GT.1.D300) THEN
         RAD_WEIGHT_FACTOR=1.D0
         RETURN
      ENDIF
      STAT=.TRUE.
      RETURN
      END


      SUBROUTINE RADCORR_NU_LOOKUP(IA,NU_TABLE,RVALUE,STAT)
      IMPLICIT NONE
      INTEGER NRAD_ANG,NRAD_NU_MAX
      PARAMETER (NRAD_ANG=15,NRAD_NU_MAX=5200)
      INTEGER IA,J,NN
      INTEGER RAD_NNU(NRAD_ANG),RAD_LOADED_IMOD
      REAL*8 NU_TABLE,RVALUE,RAD_ANGLE(NRAD_ANG)
      REAL*8 RAD_NU(NRAD_NU_MAX,NRAD_ANG)
      REAL*8 RAD_RATIO(NRAD_NU_MAX,NRAD_ANG),RAD_TABLE_EBEAM
      REAL*8 NU_LOOKUP,NU_TOL,FRAC
      LOGICAL STAT,RAD_INITIALIZED
      COMMON /RADCORR_3HE_CACHE/ RAD_ANGLE,RAD_NU,RAD_RATIO,
     > RAD_TABLE_EBEAM,RAD_NNU,RAD_LOADED_IMOD,RAD_INITIALIZED
      PARAMETER (NU_TOL=1.D-6)

      STAT=.FALSE.
      RVALUE=1.D0
      IF (IA.LT.1 .OR. IA.GT.NRAD_ANG) RETURN
      NN=RAD_NNU(IA)
      IF (NN.LT.2) RETURN
      NU_LOOKUP=NU_TABLE
      IF (NU_LOOKUP.LT.RAD_NU(1,IA)-NU_TOL .OR.
     >    NU_LOOKUP.GT.RAD_NU(NN,IA)+NU_TOL) RETURN
      IF (DABS(NU_LOOKUP-RAD_NU(1,IA)).LE.NU_TOL) THEN
         RVALUE=RAD_RATIO(1,IA)
         STAT=.TRUE.
         RETURN
      ENDIF
      IF (DABS(NU_LOOKUP-RAD_NU(NN,IA)).LE.NU_TOL) THEN
         RVALUE=RAD_RATIO(NN,IA)
         STAT=.TRUE.
         RETURN
      ENDIF
      DO J=1,NN-1
         IF (NU_LOOKUP.GT.RAD_NU(J,IA) .AND.
     >       NU_LOOKUP.LT.RAD_NU(J+1,IA)) THEN
            FRAC=(NU_LOOKUP-RAD_NU(J,IA)) /
     >           (RAD_NU(J+1,IA)-RAD_NU(J,IA))
            RVALUE=RAD_RATIO(J,IA)+FRAC*
     >              (RAD_RATIO(J+1,IA)-RAD_RATIO(J,IA))
            STAT=.TRUE.
            RETURN
         ENDIF
      ENDDO
      RETURN
      END


      BLOCK DATA RADCORR_3HE_BLOCKDATA
      IMPLICIT NONE
      INTEGER NRAD_ANG,NRAD_NU_MAX
      PARAMETER (NRAD_ANG=15,NRAD_NU_MAX=5200)
      INTEGER RAD_NNU(NRAD_ANG),RAD_LOADED_IMOD
      REAL*8 RAD_ANGLE(NRAD_ANG),RAD_NU(NRAD_NU_MAX,NRAD_ANG)
      REAL*8 RAD_RATIO(NRAD_NU_MAX,NRAD_ANG),RAD_TABLE_EBEAM
      LOGICAL RAD_INITIALIZED
      COMMON /RADCORR_3HE_CACHE/ RAD_ANGLE,RAD_NU,RAD_RATIO,
     > RAD_TABLE_EBEAM,RAD_NNU,RAD_LOADED_IMOD,RAD_INITIALIZED
      DATA RAD_INITIALIZED /.FALSE./
      DATA RAD_LOADED_IMOD /0/
      DATA RAD_TABLE_EBEAM /0.D0/
      END
