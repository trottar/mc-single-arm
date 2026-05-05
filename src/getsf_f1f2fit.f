      SUBROUTINE GETSF_F1F2fit(isf0,imod0,xb,q2,F1,F2,FL,STAT)
      IMPLICIT NONE
      CHARACTER*100 filename
      real*8 PI/3.1416/
      integer imod0,isf0,imod,isf
      real*8 xb, q2, xb_r, q2_r, xb_r1, q2_r1, xb_r2, q2_r2
      real*8 data1(40)
      real*8 F1,F2,FL
      real*8 tmp1,tmp2
      integer ifile21
      integer i, ix, iq2, nxb_F1F2FIT/3000/, nq2_F1F2FIT/0/
      REAL*8 VQ2_F1F2FIT(500)
      REAL*8 VXB_F1F2FIT(3000)
      REAL*8 VF1_FIT(3000,500,4,0:5),VF2_FIT(3000,500,4,0:5),
     >     VFL_FIT(3000,500,4,0:5)
      logical STAT,SFDEBUG,SFDEBUG1
      logical FirstcallSF_F1F2FIT /.TRUE./
      integer IOERR
      save FirstcallSF_F1F2FIT,nxb_F1F2FIT,nq2_F1F2FIT,
     >     VXB_F1F2FIT,VQ2_F1F2FIT,VF1_FIT,VF2_FIT
      SFDEBUG=.false.
      SFDEBUG1=.false.
      if (FirstcallSF_F1F2FIT) then
         FirstcallSF_F1F2FIT = .FALSE.
         do imod=1,5
         WRITE(filename,'("src/interp/sf_tables/Table_3He_F1F2_SF",I1,".csv")')imod
         ifile21=20+imod
         OPEN(ifile21,FILE=filename,STATUS='OLD',err=101)
         ix=0
         iq2=0
         do ix=1,nxb_F1F2FIT
            vxb_F1F2FIT(ix)=ix/1000.
         enddo
         READ(ifile21,*)
 100     READ(ifile21,*,err=103,end=102,IOSTAT=IOERR) q2_r1,xb_r1,
     >        data1(1),data1(2),data1(3),data1(4),data1(5),data1(6),
     >        data1(7),data1(8),data1(9),data1(10),data1(11),data1(12)
         ix=int(xb_r1/0.001+0.5)
         if (imod.eq.1) then
            if (iq2.eq.0) then
               iq2=iq2+1
               if (vq2_F1F2FIT(iq2).eq.0) vq2_F1F2FIT(iq2)=q2_r1
            elseif (q2_r1.ne.vq2_F1F2FIT(iq2)) then
               iq2=iq2+1
               if (vq2_F1F2FIT(iq2).eq.0) vq2_F1F2FIT(iq2)=q2_r1
            endif
         else
            do iq2=1,nq2_f1f2fit
               if (q2_r1.eq.vq2_F1F2FIT(iq2)) goto 300
            enddo
 300        continue
         endif
         do isf=1,4
            vF2_FIT(ix,iq2,isf,imod) = data1((isf-1)*3+1)
            vFL_FIT(ix,iq2,isf,imod) = data1((isf-1)*3+2)
            vF1_FIT(ix,iq2,isf,imod) = data1((isf-1)*3+3)
         enddo
         goto 100
 102     continue
         nq2_F1F2FIT=iq2
         CLOSE(ifile21)
         enddo
      endif
      if (q2.lt.vq2_F1F2FIT(1)) q2=vq2_F1F2FIT(1)
      xb_r=xb
      q2_r=q2
      if((q2_r.lt.vq2_F1F2FIT(1)).or.(q2_r.gt.vq2_F1F2FIT(nq2_F1F2FIT))) then
         F1=0;F2=0;FL=0;STAT=.false.;return
      endif
      if((xb_r.lt.vxb_F1F2FIT(1)).or.(xb_r.gt.vxb_F1F2FIT(nxb_F1F2FIT))) then
         F1=0;F2=0;FL=0;STAT=.false.;return
      endif
      do ix=1,nxb_F1F2FIT-1
         if ((xb.ge.vxb_F1F2FIT(ix)).and.(xb.lt.vxb_F1F2FIT(ix+1))) then
            xb_r1=vxb_F1F2FIT(ix); xb_r2=vxb_F1F2FIT(ix+1); goto 104
         endif
      enddo
 104  do iq2=1,nq2_F1F2FIT-1
         if((q2.ge.vq2_F1F2FIT(iq2)).and.(q2.lt.vq2_F1F2FIT(iq2+1)))then
            q2_r1=vq2_F1F2FIT(iq2); q2_r2=vq2_F1F2FIT(iq2+1); goto 105
         endif
      enddo
 105  tmp1 = vf1_FIT(ix,iq2+1,isf0,imod0)*(q2-q2_r1)/(q2_r2-q2_r1)
     >     + vf1_FIT(ix,iq2,isf0,imod0)*(q2-q2_r2)/(q2_r1-q2_r2)
      tmp2 = vf1_FIT(ix+1,iq2+1,isf0,imod0)*(q2-q2_r1)/(q2_r2-q2_r1)
     >     + vf1_FIT(ix+1,iq2,isf0,imod0)*(q2-q2_r2)/(q2_r1-q2_r2)
      F1 = tmp2 * (xb-xb_r1)/(xb_r2-xb_r1)
     >     + tmp1 * (xb-xb_r2)/(xb_r1-xb_r2)
      tmp1 = vf2_FIT(ix,iq2+1,isf0,imod0)*(q2-q2_r1)/(q2_r2-q2_r1)
     >     + vf2_FIT(ix,iq2,isf0,imod0)*(q2-q2_r2)/(q2_r1-q2_r2)
      tmp2 = vf2_FIT(ix+1,iq2+1,isf0,imod0)*(q2-q2_r1)/(q2_r2-q2_r1)
     >     + vf2_FIT(ix+1,iq2,isf0,imod0)*(q2-q2_r2)/(q2_r1-q2_r2)
      F2 = tmp2 * (xb-xb_r1)/(xb_r2-xb_r1)
     >     + tmp1 * (xb-xb_r2)/(xb_r1-xb_r2)
      tmp1 = vfL_FIT(ix,iq2+1,isf0,imod0)*(q2-q2_r1)/(q2_r2-q2_r1)
     >     + vfL_FIT(ix,iq2,isf0,imod0)*(q2-q2_r2)/(q2_r1-q2_r2)
      tmp2 = vfL_FIT(ix+1,iq2+1,isf0,imod0)*(q2-q2_r1)/(q2_r2-q2_r1)
     >     + vfL_FIT(ix+1,iq2,isf0,imod0)*(q2-q2_r2)/(q2_r1-q2_r2)
      FL = tmp2 * (xb-xb_r1)/(xb_r2-xb_r1)
     >     + tmp1 * (xb-xb_r2)/(xb_r1-xb_r2)
      STAT=.true.
      RETURN
 101  STAT=.false.
      return
 103  STAT=.false.
      return
      END
