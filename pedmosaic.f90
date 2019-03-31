module pedmosaic
implicit none

! Author  : Tom DRUET
! Copyright (C) 2011, 2019

!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.

! works with up to 10 generations
! version: 30/11/2011

contains

subroutine pedigree(arguments,narg)
integer ::nmarq,nclust,nhap,i,j,k,l,ori,num,verif,maxall,checkparam,narg
integer ::all1,all2,io,MTRANS(3,2046,2046),recod_ancestors(2046),ancestors(2046),path(2046)
integer, allocatable ::hap(:,:),hapin(:),nall(:)
real*8 ::val,lik
integer ::ani,nani,pere,mere,add1,add2,parent,inbreeding
integer, allocatable ::sire(:),dam(:),sexe(:)
real*8,allocatable ::scaling(:),probrec(:),posi(:)
logical, allocatable ::haplotyped(:),founder(:),founder_above(:,:)
character*50 ::paramfile,infoline,arguments(50)
character*180 ::haplofile,pedfile,markfile,namefreq
integer ::ori0,ori1,ori2,meiosis(2046,0:10),common,opposite,nmeiosis
real*8,allocatable ::alpha(:,:),beta(:,:),gamma(:,:)
real*8,allocatable ::freq(:,:)
real*8 ::posit
real*8, allocatable ::genop(:,:)

!********************************************************************************************************
!************ Part I: reading data **********************************************************************
!********************************************************************************************************

 call read_data
 call find_founders
 call random_seed

print*,'###### RUNNING CHROMIBD vs 1.2 #########'

 open(101,file='HiddenAncestors')

 ! computation of a transition matrix
 if(inbreeding==0)then
  call transition

  do l=1,nani
   if(.not.haplotyped(l))cycle
    do ori0=1,2
      num=2*l-(2-ori0)
       if(.not.founder(num))then
        ancestors=0
        call find_ancestors(num,1,0)
        call recode_ancestors
        allocate(alpha(nclust,nmarq),beta(nclust,nmarq))
        allocate(gamma(nclust,nmarq))
        call forward(num)
        call backward(num)
        deallocate(alpha,beta,gamma)
       endif
    enddo
  enddo
 else if(inbreeding==1)then
  do l=1,nani
   if(.not.haplotyped(l))cycle
    do ori0=1,2
      num=2*l-(2-ori0)
       if(.not.founder(num))then
        ancestors=0;path=0
        call find_ancestors(num,1,0)
        call recode_ancestors
        call list_meioses
        call transition_F

        allocate(alpha(nclust,nmarq),beta(nclust,nmarq))
        allocate(gamma(nclust,nmarq))
        call forward(num)
        call backward(num)
        deallocate(alpha,beta,gamma)
       endif
    enddo
  enddo

 endif

 close(101)

contains

! **************************************************************
! ********* Reading data (observations) and parameters *********
! **************************************************************

subroutine read_data
implicit none

if(narg==1)then
 print*,'Name of file with known haplotypes?'
 read(*,*)haplofile
 print*,'Name of pedigree file?'
 read(*,*)pedfile
 print*,'Name of marker file ?'
 read(*,*)markfile
 print*,'Name of file with allele frequencies?'
 read(*,*)namefreq
 print*,'Model inbreeding or not (yes=1 / no=0)?'
 read(*,*)inbreeding
else
 checkparam=0
 do i=2,narg,2 !### take arguments by pairs
  infoline=arguments(i)
  select case(infoline)
  case("--haplotypes")
   haplofile=adjustl(arguments(i+1))
   checkparam=checkparam+1
  case("--ped")
   pedfile=adjustl(arguments(i+1))
   checkparam=checkparam+1
  case("--map")
   markfile=adjustl(arguments(i+1))
   checkparam=checkparam+1
  case("--freqs")
   namefreq=adjustl(arguments(i+1))
   checkparam=checkparam+1
  case("--inbreeding")
   if(adjustl(arguments(i+1)) .eq. "0")inbreeding=0
   if(adjustl(arguments(i+1)) .eq. "1")inbreeding=1
   checkparam=checkparam+1
  case default
   write(*,*)"Unrecognized option ::",infoline
   stop
  end select
 enddo
 if(checkparam /= 5)then
   print*,'Missing information in command line!'
   print*,'Found less than five additional information with pedigree method.'
   print*,'Arguments --haplotypes, --ped, --map, --freqs and --inbreeding are required.'
   stop
 endif
endif

open(50,file=markfile)

nmarq=0
do
read(50,*,iostat=io)num
if(io/=0)exit
nmarq=nmarq+1
enddo

if(num/=nmarq)then
 print*,'Error in marker file ::'
 print*,'Last marker is ::',num
 print*,'Number of markers red is ::',nmarq
 stop
endif

print*,'Number of markers ::',nmarq
rewind(50)

allocate(probrec(nmarq),posi(nmarq),nall(nmarq),freq(2,nmarq))
allocate(scaling(nmarq))
posi=0.0;freq=0.5

do
read(50,*,iostat=io)num,markfile,posit
if(io/=0)exit
if(num<=nmarq .and. num>0)posi(num)=posit/100.0
enddo

print*,'Length of marker map (in cM / Mb) ::',(posi(nmarq)-posi(1))*100.0

probrec=0.0

do i=1,nmarq-1
   probrec(i)=0.5d0*(1.d0-dexp(-2.d0*abs(posi(i+1)-posi(i))))
enddo

verif=0

! **************************************
! ***** Reading data (observations) ****
! **************************************
nhap=0;nani=0

print*,'Name of haplotypes file ::',haplofile
open(11,file=haplofile)

do 
read(11,*,iostat=io)num
if(io/=0)exit
nani=max(num,nani)
enddo

rewind(11)

print*,'Maximum number of animal red in haplotype file ::',nani

nhap=2*nani

allocate(hap(nhap,nmarq),hapin(nmarq))
allocate(haplotyped(nani),founder(2*nani),founder_above(nani,2))
allocate(sire(nani),dam(nani),sexe(nani))

allocate(genop(nani,4))

sire=0;dam=0;hap=0
founder=.false.
genop=0.00
sexe=0

open(9,file=pedfile)

do
read(9,*,iostat=io)ani,pere,mere
if(io/=0)exit
if(ani<=nani)then
 sire(ani)=pere
 dam(ani)=mere
 if(pere/=0)then
  if(sexe(pere)==0)then
    sexe(pere)=1
  else
    if(sexe(pere)==2)then

    endif
  endif
 endif 

 if(mere/=0)then
  if(sexe(mere)==0)then
    sexe(mere)=2
  else
    if(sexe(mere)==1)then

    endif
  endif
 endif 

endif
enddo
haplotyped=.false.;nall=0


do
!read(11,'(i6,i2,1x,<nmarq>i2)',iostat=io)num,ori,(hapin(j),j=1,nmarq)
read(11,*,iostat=io)num,ori,(hapin(j),j=1,nmarq)
if(io/=0)exit
if(ori/=1 .and. ori/=2)then
 print*,'Error in known haplotypes file: origin different from 1 or 2 ::',num
 stop
endif
haplotyped(num)=.true.
do i=1,nmarq
 nall(i)=max(nall(i),hapin(i))
enddo
hap(2*num-2+ori,:)=hapin(:)
enddo


maxall=0
do i=1,nmarq
 maxall=max(maxall,nall(i))
enddo

print*,'End of reading haplotype and pedigree files '

do i=1,nani
 if(haplotyped(i))verif=verif+1
enddo

print*,'Number of animals with haplotype ::',verif

open(99,file=namefreq)

do 
read(99,*,iostat=io)i,freq(1,i),freq(2,i)
if(io/=0)exit
enddo
if(i/=nmarq)then
 print*,'Incorrect number of markers in allelic frequencies file !'
 stop
endif


end subroutine

function emission(haplotype,ancestral,marker)
implicit none
integer ::haplotype,ancestral,marker,all1,all2
real*8 ::emission

  if(hap(haplotype,marker)/=0)then
    all1=hap(haplotype,marker)
    if(ancestral/=-1)then
       all2=hap(ancestral,marker)
       if(all1==all2)then
          emission=0.999
       else
         emission=0.001
       endif
       if(all2==0)emission=1.000
!       if(all2==0)emission=freq(all1,marker)
    else ! if ancestral == -1
        emission=freq(all1,marker)
    endif
    return
  else
    emission=1.000
    return
  endif

print*,'Emission  ::',haplotype,ancestral,marker,emission

end function

function jump(add1,add2,marker)
implicit none
integer ::marker,add1,add2
real*8 ::jump,v1,v2,v3

v1=(1.0-probrec(marker))**MTRANS(1,add1,add2)
v2=probrec(marker)**MTRANS(2,add1,add2)
v3=1.000/(2.000**MTRANS(3,add1,add2))
jump=v1*v2*v3

end function

!******* forward algorithm - computation of alpha ******************

subroutine forward(nh)
implicit none
integer ::nh

alpha=0.0;scaling=0.0

! initialisation
! pi(i,1) replaced by 1.00/(2**generation)

do i=1,nclust ! (2046 possible ancestors)
  alpha(i,1)=(1.00/(2**MTRANS(1,recod_ancestors(i),recod_ancestors(i))))*emission(nh,ancestors(recod_ancestors(i)),1)
  scaling(1)=scaling(1)+alpha(i,1)
enddo

 scaling(1)=1.0/scaling(1)
 alpha(:,1)=alpha(:,1)*scaling(1)

! induction: to get to i, completely defined by "jump"
 
 do k=2,nmarq
  do i=1,nclust
   do j=1,nclust
    alpha(i,k)=alpha(i,k)+alpha(j,k-1)*jump(recod_ancestors(j),recod_ancestors(i),k-1)*emission(nh,ancestors(recod_ancestors(i)),k)
   enddo
    scaling(k)=scaling(k)+alpha(i,k)
  enddo
 scaling(k)=1.0/scaling(k)
 alpha(:,k)=alpha(:,k)*scaling(k)
 enddo  

! termination

 val=-sum(log(scaling(:)))

lik=lik+val

end subroutine

!******* backward algorithm - estimation of beta ******************
!******* estimation of gamma = alpha*beta *************************
!******* gamma = probability of state i,j at marker k *************

subroutine backward(nh)
implicit none
integer ::nh

beta=0.0;gamma=0.0

! initialisation

do i=1,nclust ! nclust = 30 possible founders
 gamma(i,nmarq)=alpha(i,nmarq)*1.0  ! beta(i,j,nmarq)=1.0
 beta(i,nmarq)=1.0*scaling(nmarq)
enddo

! induction
! to arrive in k: with or without transition

do k=nmarq-1,1,-1
 do i=1,nclust
  do j=1,nclust
    beta(i,k)=beta(i,k)+jump(recod_ancestors(i),recod_ancestors(j),k)*emission(nh,ancestors(recod_ancestors(j)),k+1)*beta(j,k+1)
  enddo
  beta(i,k)=beta(i,k)*scaling(k)
  gamma(i,k)=alpha(i,k)*beta(i,k)/scaling(k)
 enddo
enddo

do k=1,nmarq
 do i=1,nclust
  write(101,*)nh,ancestors(recod_ancestors(i)),k,gamma(i,k)
 enddo
enddo

end subroutine

!****************************************************************************************************
!*********************** subroutines for ancestors contribution ************************************* 
!****************************************************************************************************

subroutine find_founders
implicit none

! animal are sorted
! first genotyped haplotypes are founders
! if unknown parents and genotyped=> founder
! if known parent but all ancestors ungenotyped => founder
! founder_above => true if genotyped or if parent genotyped

founder=.false.;founder_above=.false.

do i=1,nani
 if(haplotyped(i))then

   if(sire(i)==0)then
     founder(2*i-1)=.true.
     founder_above(i,1)=.true.
   else
     if(founder_above(sire(i),1) .or. founder_above(sire(i),2))then
        founder_above(i,1)=.true.
     else
        founder(2*i-1)=.true.
        founder_above(i,1)=.true.
     endif
   endif

   if(dam(i)==0)then
     founder(2*i)=.true.
     founder_above(i,2)=.true.
   else
     if(founder_above(dam(i),1) .or. founder_above(dam(i),2))then
        founder_above(i,2)=.true.
     else
        founder(2*i)=.true.
        founder_above(i,2)=.true.
     endif
   endif

 else

  if(sire(i)/=0)then
   if(founder_above(sire(i),1) .or. founder_above(sire(i),2))founder_above(i,1)=.true.
  endif

  if(dam(i)/=0)then
   if(founder_above(dam(i),1) .or. founder_above(dam(i),2))founder_above(i,2)=.true.
  endif
   
 endif
enddo

end subroutine

recursive subroutine find_ancestors(hapid,generation,add0)
integer ::i,add1,add2,addi,base
integer ::generation,id,add0,hapid,ori,parent

if(mod(hapid,2)==0)then
 id=hapid/2;ori=2
else
 id=(hapid+1)/2;ori=1
endif

add1=add0+1
add2=add0+2

if(ori==1)then
 parent=sire(id)
else
 parent=dam(id)
endif

if(parent/=0)then
 path(add1)=2*parent-1;path(add2)=2*parent
 if(haplotyped(parent))then
   ancestors(add1)=2*parent-1
   ancestors(add2)=2*parent
 endif
endif

base=0

do i=1,generation-1
 base=base+2**i
enddo
addi=base+2**generation

add1=addi+2*(add1-base-1)
add2=addi+2*(add2-base-1)

if(parent/=0 )then
 if(.not.haplotyped(parent))then
  if(founder_above(parent,1) .and. generation<11)then
   call find_ancestors(2*parent-1,generation+1,add1)
  else
   ancestors(add0+1)=-1
  endif
  if(founder_above(parent,2) .and. generation<11)then
   call find_ancestors(2*parent,generation+1,add2)
  else
   ancestors(add0+2)=-1
  endif
 endif
endif

end subroutine


subroutine transition
integer ::g1,g2,n1,n2,add1,add2,i,j,k,l
integer ::join_generation,num1,num2,n_join
real ::v1,v2

MTRANS=0
add1=0

do g1=1,10
 n1=2**g1
 do i=1,n1
  add1=add1+1
  add2=0
  do g2=1,10
   n2=2**g2
    do j=1,n2
      add2=add2+1
!   test if joined one generation above => add2 is ancestor of add1 or vice verse, no transition
      if(g2>g1)then
        if(int((j-1)/(n2/n1))+1==i)cycle ! add2 is ancestor of add1
      else if(g1>g2)then
        if(int((i-1)/(n1/n2))+1==j)cycle ! add1 is ancestor of add2
      endif

      if(add1==add2)then
        MTRANS(1,add1,add2)=g1
      else
        join_generation=min(g1,g2)
        n_join=2**join_generation
         num1=int((i-1)/(2*n1/n_join))+1
         num2=int((j-1)/(2*n2/n_join))+1
        do while(num1/=num2 .and. join_generation>-2)
          join_generation=join_generation-1
          n_join=2**join_generation
         num1=int((i-1)/(2*n1/n_join))+1
         num2=int((j-1)/(2*n2/n_join))+1
        enddo
        MTRANS(3,add1,add2)=g2-join_generation
        MTRANS(2,add1,add2)=1
        MTRANS(1,add1,add2)=join_generation-1
      endif
    enddo
   enddo
 enddo
enddo

end subroutine

subroutine recode_ancestors
integer ::i

recod_ancestors=0
nclust=0

do i=1,2046
 if(ancestors(i)/=0)then
  nclust=nclust+1
  recod_ancestors(nclust)=i
 endif
enddo

end subroutine

subroutine list_meioses
implicit none
integer ::add1,add2,g1,n1


meiosis=0;add1=0;add2=0
meiosis(:,0)=-1
do g1=1,10
 n1=2**(g1-1)
 do i=1,n1
! paternal
  add1=add1+1
  meiosis(add1,g1)=path(add1)
! maternal
  add1=add1+1
  meiosis(add1,g1)=path(add1)
  if(g1>1)then
    add2=add2+1
    meiosis(add1-1,1:g1-1)=meiosis(add2,1:g1-1)
    meiosis(add1,1:g1-1)=meiosis(add2,1:g1-1)
  endif
 enddo
enddo

end subroutine

subroutine count_events(pc1,pc2)
implicit none
integer ::pc1,pc2,ani1,ani2,i,j,k,l,io

common=0
opposite=0
nmeiosis=0

do i=0,9
 if(meiosis(pc1,i+1)==0)exit
 do j=0,9
  if(meiosis(pc2,j+1)==0)exit
  if(i==0)nmeiosis=nmeiosis+1
  if(meiosis(pc1,i)==meiosis(pc2,j) .and. meiosis(pc1,i+1)==meiosis(pc2,j+1))common=common+1
  if(meiosis(pc1,i)==meiosis(pc2,j) .and. meiosis(pc1,i+1)/=meiosis(pc2,j+1))opposite=opposite+1
 enddo
enddo

end subroutine

subroutine transition_F
integer ::g1,g2,n1,n2,add1,add2,i,j,k,l
integer ::join_generation,num1,num2,n_join
real ::v1,v2

MTRANS=0
add1=0

do g1=1,10
 n1=2**g1
 do i=1,n1
  add1=add1+1
  if(ancestors(add1)==0)cycle
  add2=0
  do g2=1,10
   n2=2**g2
    do j=1,n2
      add2=add2+1
      if(ancestors(add2)==0)cycle
!   test if joined one generation above => add2 is ancestor of add1 or vice verse, no transition
      if(g2>g1)then
        if(int((j-1)/(n2/n1))+1==i)cycle ! add2 is ancestor of add1
      else if(g1>g2)then
        if(int((i-1)/(n1/n2))+1==j)cycle ! add1 is ancestor of add2
      endif

      call count_events(add1,add2)

      MTRANS(3,add1,add2)=nmeiosis-opposite-common
      MTRANS(2,add1,add2)=opposite
      MTRANS(1,add1,add2)=common

    enddo
   enddo
 enddo
enddo

end subroutine

end subroutine

end module
