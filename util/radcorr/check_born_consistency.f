      PROGRAM CHECK_BORN_CONSISTENCY
C Compare production GETSF_F1F2fit/Mott Born cross sections with exact
C XSborn_unp CSV rows.  This is QA only: no difference is a production cut.
      IMPLICIT NONE
      INTEGER NPOINT,IMOD(11),I,NFAIL
      PARAMETER (NPOINT=11)
      REAL*8 THETA(11),NU(11)
      REAL*8 EBEAM,MP,ALPHA,GEV2_TO_NB
      REAL*8 EPRIME,Q2_MC,X_MC,TH2,MOTT,F1,F2,FL,SIGMA_MC
      REAL*8 EBEAM_TAB,X_TAB,Q2_TAB,XS_BORN,XS_RAD,TABLE_RATIO
      REAL*8 READER_RATIO,DELTA
      REAL*8 SUM_ABS(5),MAX_ABS(5),MIN_DELTA(5),MAX_DELTA(5)
      INTEGER NMODEL(5)
      LOGICAL STAT,SF_STAT
      DATA IMOD /1,1,1,3,3,3,5,5,5,1,5/
      DATA THETA /30.D0,30.D0,30.D0,30.D0,30.D0,30.D0,
     >            30.D0,30.D0,30.D0,27.D0,32.D0/
      DATA NU /6.100D0,7.500D0,9.000D0,6.100D0,7.500D0,
     >         9.000D0,6.100D0,7.500D0,9.000D0,6.000D0,
     >         7.500D0/
      DATA EBEAM /10.38D0/
      DATA MP /0.93827208D0/
      DATA ALPHA /7.2973525693D-3/
      DATA GEV2_TO_NB /3.89379338D5/

      NFAIL=0
      DO I=1,5
         NMODEL(I)=0
         SUM_ABS(I)=0.D0
         MAX_ABS(I)=0.D0
         MIN_DELTA(I)=1.D99
         MAX_DELTA(I)=-1.D99
      ENDDO
      WRITE(6,*) 'Born consistency QA: sigma_mc is nb/sr; CSV XS units'
      WRITE(6,*) 'are compared numerically without an acceptance threshold.'

      DO I=1,NPOINT
         CALL INIT_RADCORR_3HE(IMOD(I),EBEAM,STAT)
         IF (.NOT.STAT) THEN
            WRITE(6,*) 'FAIL: radcorr initialization for SF',IMOD(I)
            NFAIL=NFAIL+1
            GOTO 200
         ENDIF
         CALL READ_RADCORR_EXACT(IMOD(I),THETA(I),1000.D0*NU(I),
     >        EBEAM_TAB,X_TAB,Q2_TAB,XS_BORN,XS_RAD,STAT)
         IF (.NOT.STAT) THEN
            WRITE(6,*) 'FAIL: exact CSV row for SF',IMOD(I)
            NFAIL=NFAIL+1
            GOTO 200
         ENDIF
         CALL GET_RADCORR_3HE(IMOD(I),THETA(I),NU(I),
     >        READER_RATIO,STAT)
         IF (.NOT.STAT) THEN
            WRITE(6,*) 'FAIL: reader lookup for SF',IMOD(I)
            NFAIL=NFAIL+1
            GOTO 200
         ENDIF
         TABLE_RATIO=XS_RAD/XS_BORN
         IF (DABS(READER_RATIO-TABLE_RATIO).GT.1.D-12) THEN
            WRITE(6,*) 'FAIL: reader/table ratio mismatch for SF',IMOD(I)
            WRITE(6,*) ' reader=',READER_RATIO,' table=',TABLE_RATIO
            NFAIL=NFAIL+1
         ENDIF

         EPRIME=EBEAM-NU(I)
         TH2=(THETA(I)/57.29577951308232D0)/2.D0
         Q2_MC=4.D0*EBEAM*EPRIME*SIN(TH2)**2
         X_MC=Q2_MC/(2.D0*MP*NU(I))
         CALL GETSF_F1F2fit(4,IMOD(I),X_MC,Q2_MC,F1,F2,FL,SF_STAT)
         IF (.NOT.SF_STAT) THEN
            WRITE(6,*) 'FAIL: GETSF_F1F2fit for SF',IMOD(I)
            NFAIL=NFAIL+1
            GOTO 200
         ENDIF
         MOTT=((ALPHA*COS(TH2)/(2.D0*EBEAM*SIN(TH2)*SIN(TH2)))**2)
     >        *GEV2_TO_NB
         SIGMA_MC=MOTT*(F2/NU(I)+2.D0*(F1/MP)*TAN(TH2)**2)
         DELTA=(SIGMA_MC-XS_BORN)/XS_BORN
         NMODEL(IMOD(I))=NMODEL(IMOD(I))+1
         SUM_ABS(IMOD(I))=SUM_ABS(IMOD(I))+DABS(DELTA)
         IF (DABS(DELTA).GT.MAX_ABS(IMOD(I))) THEN
            MAX_ABS(IMOD(I))=DABS(DELTA)
         ENDIF
         IF (DELTA.LT.MIN_DELTA(IMOD(I))) MIN_DELTA(IMOD(I))=DELTA
         IF (DELTA.GT.MAX_DELTA(IMOD(I))) MAX_DELTA(IMOD(I))=DELTA
         WRITE(6,*) 'SF=',IMOD(I),' theta=',THETA(I),' nu=',NU(I)
         WRITE(6,*) ' x_mc,x_csv=',X_MC,X_TAB
         WRITE(6,*) ' q2_mc,q2_csv=',Q2_MC,Q2_TAB
         WRITE(6,*) ' F1,F2=',F1,F2
         WRITE(6,*) ' sigma_born_mc,XSborn_unp=',SIGMA_MC,XS_BORN
         WRITE(6,*) ' XSrad_unp=',XS_RAD,' table ratio=',TABLE_RATIO
         WRITE(6,*) ' reader ratio=',READER_RATIO,' delta_born=',DELTA
  200    CONTINUE
      ENDDO

      DO I=1,5
         IF (NMODEL(I).GT.0) THEN
            WRITE(6,*) 'SF',I,' delta_born mean_abs=',
     >       SUM_ABS(I)/DBLE(NMODEL(I)),' max_abs=',MAX_ABS(I)
            WRITE(6,*) 'SF',I,' delta_born range=',MIN_DELTA(I),
     >       MAX_DELTA(I)
         ENDIF
      ENDDO
      IF (NFAIL.NE.0) THEN
         WRITE(6,*) 'FAIL: born consistency QA failures =',NFAIL
         STOP 1
      ENDIF
      WRITE(6,*) 'PASS: born consistency reader-ratio checks'
      END


      SUBROUTINE READ_RADCORR_EXACT(IMOD,THETA,NU_LOOKUP,
     > EBEAM_TAB,X_TAB,Q2_TAB,XS_BORN,XS_RAD,STAT)
      IMPLICIT NONE
      INTEGER IMOD,IOERR,UNITNO,K
      REAL*8 THETA,NU_LOOKUP,EBEAM_TAB,X_TAB,Q2_TAB,XS_BORN,XS_RAD
      REAL*8 COL(17)
      CHARACTER*200 FILENAME
      CHARACTER*512 LINE
      LOGICAL STAT

      STAT=.FALSE.
      EBEAM_TAB=0.D0
      X_TAB=0.D0
      Q2_TAB=0.D0
      XS_BORN=0.D0
      XS_RAD=0.D0
      WRITE(FILENAME,100) IMOD,THETA
  100 FORMAT('interp/radcorr_tables/Newfit_20260710_fullxquad_15',
     > 'angles/SF',I1,'_G1F1cmplt_QE95/radiated_data_',F4.1,
     > 'deg_short.csv')
      UNITNO=93
      OPEN(UNITNO,FILE=FILENAME,STATUS='OLD',IOSTAT=IOERR)
      IF (IOERR.NE.0) RETURN
      READ(UNITNO,'(A)',IOSTAT=IOERR) LINE
      IF (IOERR.NE.0) THEN
         CLOSE(UNITNO)
         RETURN
      ENDIF
   10 CONTINUE
      READ(UNITNO,'(A)',IOSTAT=IOERR) LINE
      IF (IOERR.NE.0) THEN
         CLOSE(UNITNO)
         RETURN
      ENDIF
      READ(LINE,*,IOSTAT=IOERR) (COL(K),K=1,17)
      IF (IOERR.NE.0) THEN
         CLOSE(UNITNO)
         RETURN
      ENDIF
      IF (DABS(COL(2)-NU_LOOKUP).GT.1.D-6) GOTO 10
      EBEAM_TAB=COL(1)
      X_TAB=COL(3)
      Q2_TAB=COL(4)/1.D6
      XS_BORN=COL(5)
      XS_RAD=COL(10)
      IF (XS_BORN.GT.0.D0 .AND. XS_RAD.GT.0.D0) STAT=.TRUE.
      CLOSE(UNITNO)
      RETURN
      END
