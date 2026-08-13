      PROGRAM TEST_GET_RADCORR_3HE
C Fixed-kinematics validation for the selected-model 3He table reader.
C Run this executable from src/ so its relative table paths resolve.
      IMPLICIT NONE
      LOGICAL STAT
      INTEGER NFAIL
      REAL*8 R,RLOW,RHIGH,RMID,EXPECTED
      REAL*8 SIGMA_BORN,SIGMA_USED

      NFAIL=0
      CALL INIT_RADCORR_3HE(1,10.38D0,STAT)
      CALL ASSERT_TRUE('SF1 initialization',STAT,NFAIL)

      CALL GET_RADCORR_3HE(1,30.0D0,6.007D0,R,STAT)
      CALL ASSERT_TRUE('SF1 exact table point',STAT,NFAIL)
      EXPECTED=1.1542D-5/3.8235D-4
      CALL ASSERT_CLOSE('SF1 exact ratio',R,EXPECTED,1.D-12,NFAIL)

      CALL GET_RADCORR_3HE(1,30.0D0,6.0075D0,RMID,STAT)
      CALL ASSERT_TRUE('SF1 nu interpolation',STAT,NFAIL)
      EXPECTED=0.5D0*(1.1542D-5/3.8235D-4+
     >                 2.0685D-5/6.8522D-4)
      CALL ASSERT_CLOSE('SF1 nu interpolation value',RMID,EXPECTED,
     >                  1.D-12,NFAIL)

      CALL GET_RADCORR_3HE(1,30.0D0,6.100D0,RLOW,STAT)
      CALL ASSERT_TRUE('SF1 lower angle point',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,30.5D0,6.100D0,RHIGH,STAT)
      CALL ASSERT_TRUE('SF1 upper angle point',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,30.25D0,6.100D0,R,STAT)
      CALL ASSERT_TRUE('SF1 angle interpolation',STAT,NFAIL)
      EXPECTED=0.5D0*(RLOW+RHIGH)
      CALL ASSERT_CLOSE('SF1 angle interpolation value',R,EXPECTED,
     >                  1.D-12,NFAIL)

      CALL GET_RADCORR_3HE(1,30.0D0,6.1005D0,RLOW,STAT)
      CALL ASSERT_TRUE('SF1 lower bilinear point',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,30.5D0,6.1005D0,RHIGH,STAT)
      CALL ASSERT_TRUE('SF1 upper bilinear point',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,30.25D0,6.1005D0,R,STAT)
      CALL ASSERT_TRUE('SF1 bilinear interpolation',STAT,NFAIL)
      EXPECTED=0.5D0*(RLOW+RHIGH)
      CALL ASSERT_CLOSE('SF1 bilinear interpolation value',R,EXPECTED,
     >                  1.D-12,NFAIL)

C Exact angle and nu endpoints must not require an out-of-range neighbor.
      CALL GET_RADCORR_3HE(1,26.0D0,5.260D0,R,STAT)
      CALL ASSERT_TRUE('26 degree lower nu endpoint',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,33.0D0,6.483D0,R,STAT)
      CALL ASSERT_TRUE('33 degree lower nu endpoint',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,30.0D0,9.910D0,R,STAT)
      CALL ASSERT_TRUE('9910 MeV upper nu endpoint',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,25.999D0,6.007D0,R,STAT)
      CALL ASSERT_FALSE('angle below coverage',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,33.001D0,6.483D0,R,STAT)
      CALL ASSERT_FALSE('angle above coverage',STAT,NFAIL)
      CALL GET_RADCORR_3HE(1,30.0D0,9.911D0,R,STAT)
      CALL ASSERT_FALSE('nu above finite coverage',STAT,NFAIL)

      CALL INIT_RADCORR_3HE(3,10.38D0,STAT)
      CALL ASSERT_TRUE('SF3 initialization',STAT,NFAIL)
      CALL GET_RADCORR_3HE(3,30.0D0,6.007D0,R,STAT)
      CALL ASSERT_TRUE('SF3 exact table point',STAT,NFAIL)
      EXPECTED=2.1172D-5/7.0138D-4
      CALL ASSERT_CLOSE('SF3 exact ratio',R,EXPECTED,1.D-12,NFAIL)

      CALL INIT_RADCORR_3HE(5,10.38D0,STAT)
      CALL ASSERT_TRUE('SF5 initialization',STAT,NFAIL)
      CALL GET_RADCORR_3HE(5,30.0D0,6.007D0,R,STAT)
      CALL ASSERT_TRUE('SF5 exact table point',STAT,NFAIL)
      EXPECTED=2.8060D-5/9.2955D-4
      CALL ASSERT_CLOSE('SF5 exact ratio',R,EXPECTED,1.D-12,NFAIL)

C The production weighting invariant is independent of random transport.
      SIGMA_BORN=7.25D0
      SIGMA_USED=SIGMA_BORN*R
      CALL ASSERT_CLOSE('RC-on weighting invariant',
     >                  SIGMA_USED/SIGMA_BORN,R,1.D-14,NFAIL)
      R=1.D0
      SIGMA_USED=SIGMA_BORN*R
      CALL ASSERT_CLOSE('RC-off weighting invariant',
     >                  SIGMA_USED/SIGMA_BORN,1.D0,1.D-14,NFAIL)

      IF (NFAIL.NE.0) THEN
         WRITE(6,*) 'FAIL: get_radcorr_3he tests =',NFAIL
         STOP 1
      ENDIF
      WRITE(6,*) 'PASS: get_radcorr_3he fixed-kinematics tests'
      END


      SUBROUTINE ASSERT_TRUE(NAME,VALUE,NFAIL)
      IMPLICIT NONE
      CHARACTER*(*) NAME
      LOGICAL VALUE
      INTEGER NFAIL
      IF (.NOT.VALUE) THEN
         WRITE(6,*) 'FAIL: ',NAME
         NFAIL=NFAIL+1
      ENDIF
      RETURN
      END


      SUBROUTINE ASSERT_FALSE(NAME,VALUE,NFAIL)
      IMPLICIT NONE
      CHARACTER*(*) NAME
      LOGICAL VALUE
      INTEGER NFAIL
      IF (VALUE) THEN
         WRITE(6,*) 'FAIL: ',NAME
         NFAIL=NFAIL+1
      ENDIF
      RETURN
      END


      SUBROUTINE ASSERT_CLOSE(NAME,ACTUAL,EXPECTED,TOL,NFAIL)
      IMPLICIT NONE
      CHARACTER*(*) NAME
      REAL*8 ACTUAL,EXPECTED,TOL
      INTEGER NFAIL
      IF (DABS(ACTUAL-EXPECTED).GT.TOL) THEN
         WRITE(6,*) 'FAIL: ',NAME
         WRITE(6,*) ' actual=',ACTUAL,' expected=',EXPECTED
         NFAIL=NFAIL+1
      ENDIF
      RETURN
      END
