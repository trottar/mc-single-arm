      subroutine hms_hbook_init(filename,spec_ntuple)
C Initialize output names differently based on which spectrometer we
C are using
C
      implicit none
      include 'hbook.inc'
      character*80 filename
      logical spec_ntuple
      character*16      hut_nt_names(33)/
     >     'hsxfp', 'hsyfp', 'hsxpfp', 'hsypfp',
     >     'hsxtari','hsytari','hsxptari','hsyptari',
     >     'hsztari','hsdeltai','hsytar','hsxptar',
     >     'hsyptar','hsztar','hsdelta','fry',
     >     'xsnum','ysnum','xsieve','ysieve','stop_id',
     >     'hsvxi','hsvyi',
     >     'kin_ebeam','kin_eprime','kin_q2','kin_w2','kin_w',
     >     'kin_nu','kin_xbj','kin_theta','F1in21','F2in21' /
      integer*4 i
      NtupleIO=30
      NtupleSize=33

      open(NtupleIO,file=filename,form="unformatted",access="sequential")

      write(NtupleIO) NtupleSize
      do i=1,NtupleSize
         write(NtupleIO) hut_nt_names(i)
      enddo

      return
      end
      
