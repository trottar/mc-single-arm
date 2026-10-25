program sf_probe
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer :: ios, point_id, wfn, i
  real(8) :: ebeam,q2,w,w2,fp,flp,f2p,fn,fln,f2n,fd,f2d
  real(8) :: gamma,ypass,fy11,fy12,fy2,z,a
  logical :: firsty, exists
  character(40) :: tables(3)
  tables = [character(40) :: 'f1f2tables/f11-onshell-cdbonn.dat', &
    'f1f2tables/f12-onshell-cdbonn.dat', 'f1f2tables/f22-onshell-cdbonn.dat']
  do i=1,3
    inquire(file=trim(tables(i)),exist=exists)
    if (.not.exists) error stop 'required smearing table missing'
  enddo
  gamma=1.1d0
  ypass=1.d0
  wfn=2
  firsty=.true.
  call getfy(gamma,ypass,wfn,fy11,fy12,fy2,firsty)
  if (.not.(fy11>0.d0.and.fy2>0.d0)) error stop 'smearing initialization failed'
  do
    read(*,*,iostat=ios) point_id,ebeam,q2,w
    if(ios<0) exit
    if(ios/=0) error stop 'invalid input row'
    if(.not.all(ieee_is_finite([ebeam,q2,w]))) error stop 'nonfinite kinematics'
    if(ebeam<=0.d0.or.q2<=0.d0.or.w<=1.08d0.or.q2>=30.d0.or.w*w>=20.d0) &
      error stop 'outside benchmark domain'
    w2=w*w
    call sf(w2,q2,fp,flp,f2p,fn,fln,f2n)
    z=1.d0
    a=2.d0
    call f1f2in21(z,a,q2,w2,fd,f2d)
    if(.not.all(ieee_is_finite([fp,f2p,fn,f2n,fd,f2d]))) error stop 'nonfinite model'
    if(min(fp,fn,fd,f2p,f2n,f2d)<=0.d0) error stop 'nonpositive model'
    write(*,'(I6,1X,9(ES25.16E3,1X))') point_id,ebeam,q2,w,fp,f2p,fn,f2n,fd,f2d
  enddo
end program
