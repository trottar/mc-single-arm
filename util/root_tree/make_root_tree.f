      program make_root_tree
C     Program to convert simc .bin file to root tree
      implicit none

      character*80 rawname,filename,treefilename
      character*16 NtupleTag(80),varname
      

      integer i,j,nev,check
      integer*4 io
      integer*4 NtupleSize
      integer*4 weight_idx,rate_idx
      integer GetNumBranches
      external GetNumBranches
      logical*4 apply_norm
      real*8 generated_trials
      external read_trial_norm

      real*8 ntup(80)
      real*8 ntup_out(80)

      parameter(nev=10000000)
      io=99
      weight_idx = 0
      rate_idx = 0
      apply_norm = .false.
      generated_trials = 0.d0

c input filename
      write(6,*) 'Enter filename to convert (without .bin extension)'
      read(5,*) rawname
      call read_trial_norm(rawname,apply_norm,generated_trials)
      i=index(rawname,' ')
      filename='../../worksim/'//rawname(1:i-1)//'.bin'
      write(6,*) 'opening file: ',filename
      open(io,file=filename,form="unformatted",access="sequential")

c output filename
      treefilename='../../worksim/'//rawname(1:i-1)//'.root'
      
      read(io) NtupleSize
      call InitRootNT(treefilename,'RECREATE');

      write(6,*) 'Variables in output file:'
      do i=1,NtupleSize
         read(io) NtupleTag(i)
         if (NtupleTag(i).eq.'weight') weight_idx = i
         if (NtupleTag(i).eq.'rate_hz') rate_idx = i
         call AddNtBranch(ntup_out(i),NtupleTag(i))
         write(6,*) NtupleTag(i)
      enddo

      if (apply_norm) then
         write(6,*) 'Normalizing weight/rate_hz by generated trials = ',
     >              generated_trials
      else
         write(6,*) 'Using stored weight/rate_hz values',
     >              ' without renormalization.'
      endif

c now loop over events     
      do j=1,nev
         do i=1,NtupleSize
            read(io,iostat=check) ntup(i)
            if (check.lt.0) then
               call PrintNT()
               call RootNTOutp();
               stop
            endif
            ntup_out(i)=ntup(i)
         enddo ! loop over ntuple variables
         if (apply_norm.and.generated_trials.gt.0.d0) then
            if (weight_idx.gt.0) then
               ntup_out(weight_idx) =
     >              ntup_out(weight_idx)/generated_trials
            endif
            if (rate_idx.gt.0) then
               ntup_out(rate_idx) =
     >              ntup_out(rate_idx)/generated_trials
            endif
         endif
         call FillNTBranch('all')
      enddo ! loop over events

      call PrintNT()
      call RootNTOutp();

      end

      subroutine read_trial_norm(rawname,apply_norm,generated_trials)
      implicit none
      character*80 rawname
      logical*4 apply_norm
      real*8 generated_trials
      character*132 summaryfile,line
      integer*4 io2,check,i,itmp

      apply_norm = .false.
      generated_trials = 0.d0
      io2 = 98

      i=index(rawname,' ')
      if (i.eq.0) i=len(rawname)+1
      summaryfile='../../outfiles/'//rawname(1:i-1)//'.out'

      open(io2,file=summaryfile,status='old',iostat=check)
      if (check.ne.0) then
         write(6,*) 'WARNING: could not open summary file ',summaryfile
         return
      endif

 10   continue
      read(io2,'(A)',end=20,iostat=check) line
      if (check.ne.0) goto 20
      if (index(line,'Monte-Carlo trials generated').gt.0 .or.
     >    index(line,
     >    'Actual generated trials for normalization').gt.0 .or.
     >    index(line,
     >    'Generated-trial normalization denominator').gt.0 .or.
     >    index(line,
     >    'Event weight normalization denominator').gt.0) then
         read(line(1:11),'(i11)',iostat=check) itmp
         if (check.eq.0 .and. itmp.gt.0) then
            generated_trials = dble(itmp)
            apply_norm = .true.
         endif
         goto 20
      endif
      goto 10

 20   continue
      close(io2)
      return
      end
