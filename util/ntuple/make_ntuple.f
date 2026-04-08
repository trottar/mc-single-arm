      program make_ntuple
C     Program to convert simc .bin file to ntupl
      implicit none

*-sizes of CERNLIB working space

      integer HbookSize,HigzSize,KuipSize,PawSize
      parameter (HbookSize = 100000)
      parameter (HigzSize =  50000)
      parameter (KuipSize =  75000)
      parameter (PawSize = HigzSize+KuipSize+HbookSize+100000)

*-CERNLIB working space

      integer CernMemory(PawSize)
      common /PAWC/ CernMemory  !special nonstandard name!

C Ntuple ID stuff

      integer*4	defaultID
      parameter     	(defaultID = 666)
      character*132 	NtupleDirectory
      character*16    NtupleName
      integer*4       NtupleID,NtupleIO,NtupleSize
      integer*4 recl,status, cycle, bank

      character*16 NtupleTag(80)
      character*80 rawname,filename,ntfilename,directory
      character*16 title
      integer io,i,j,check
      integer chanout
      integer nev
      integer weight_idx,rate_idx
      real*8 ntup(80)
      real*4 ntup_out(80)
      logical*4 apply_norm
      real*8 generated_trials
      external read_weight_norm

      parameter(nev=10000000)

      io=10
      chanout=2
      NtupleID=defaultID
      recl=4096
      bank=8000
      weight_idx = 0
      rate_idx = 0
      apply_norm = .false.
      generated_trials = 0.d0

      NtupleName='MCNtuple'
      title='SIMTUPLE'

      call hlimit(PawSize)
c input filename
      write(6,*) 'Enter filename to convert (without .bin extension)'
      read(5,*) rawname
      call read_weight_norm(rawname,apply_norm,generated_trials)
      i=index(rawname,' ')
      filename='../../worksim/'//rawname(1:i-1)//'.bin'
      write(6,*) 'opening file: ',filename
      open(io,file=filename,form="unformatted",access="sequential")

c output filename
      ntfilename='../../worksim/'//rawname(1:i-1)//'.rzdat'
      call HCDIR(directory,'R')
      call HROPEN(chanout,NtupleName,ntfilename,'N',recl,status)

	if (status.ne.0)then
	  write(6,*) 'HROPEN error: istat=',status
	  stop
	endif


      read(io) NtupleSize
      write(6,*) 'Variables in output file: ',NtupleSize
      do i=1,NtupleSize
         read(io) NtupleTag(i)
         if (NtupleTag(i).eq.'weight') weight_idx = i
         if (NtupleTag(i).eq.'rate_hz') rate_idx = i
         write(6,*) NtupleTag(i)
      enddo

      if (apply_norm) then
         write(6,*) 'Normalizing weight/rate_hz by generated trials = ',
     >              generated_trials
      else
         write(6,*) 'Using stored weight/rate_hz values',
     >              ' without renormalization.'
      endif

      call HBOOKN(NtupleID,title,NtupleSize,NtupleName,bank,NtupleTag) !create Ntuple
      
      call HCDIR(NtupleDirectory,'R') !record Ntuple directory

      call HCDIR(directory,' ') !reset CERNLIB directory


c now loop over events     
      do j=1,nev
         do i=1,NtupleSize
            read(io,iostat=check) ntup(i)
            if (check.lt.0) then
               write(6,*) 'end of file'
               write(6,*) 'processed ',j, 'events'
               cycle=0
               call HCDIR(NtupleDirectory,' ')
               call HROUT(NtupleID,cycle,' ') !flush CERNLIB buffers
               call HREND(NtupleName) !CERNLIB close file
               write(6,*)'Closing file:',filename(1:60)
               close(chanout)
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
         call HFN(NtupleID,ntup_out)         
      enddo ! loop over events

      end

      subroutine read_weight_norm(rawname,apply_norm,generated_trials)
      implicit none
      character*80 rawname
      logical*4 apply_norm
      real*8 generated_trials
      character*132 summaryfile,line
      integer io2,check,i,itmp

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
      if (index(line,
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
