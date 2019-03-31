program chromibd
use pedmosaic
use unrelated
implicit none
integer ::i,j,k,l,io,narg
character(50):: args
character(50), allocatable ::infoline(:)

narg=command_argument_count()
if(narg<1)THEN
 write(*,*)"At least one argument should follow CHROMIBD invocation : Process will stop"
 stop
else
 allocate(infoline(narg))
 do i=1,narg
  call get_command_argument(i,args)
  infoline(i)=adjustl(args)
 end do
 select case(adjustl(infoline(1)))
  case("--pedigree")
   call pedigree(infoline,narg)
  case("--unrelated")
   call hmmref(infoline,narg)
  case default
   write(*,*)"Unrecognized method:",infoline(1)
   stop
 end select
endif

end program

