! PROGRAM : HiddenREFs
! Author  : Tom DRUET, Pierre Faux
! Copyright (C) 2015,2019
! version 27/03/2019


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

module unrelated
implicit none

contains

subroutine hmmref(arguments,narg)
implicit none
integer ::nmarq,nclust,nhap,i,j,k,l,t,ori,useold,maxall,io,nref,num,bestallele,narg
integer*1, allocatable ::hap(:,:),hapin(:),refhap(:,:)
integer*2, allocatable ::states(:,:), work1(:)
integer, allocatable ::nall(:), topfb(:)
real*8,allocatable ::pi(:,:),aij(:)
real*8 ::val,lik,lik2,bestlik,p1,p2
integer ::k1,k2,add1,add2,parent,seed,NbT,thread,algo,checkparam
real*8 ::verif,verif1,verif2,verif0,position,dosage,ngen,gerr,Ne,theta
real*8,allocatable ::probrec(:),posi(:)
character*20:: fmt1, fmt2
character*50 ::haplofile,haplofile2,genofile,pedfile,oldfile,statesfile,markfile
character*50 ::arguments(50),info1,info2

 
 call read_data
 call parameter_init

 write(fmt1,'(a,i7,a)') '(',nmarq,'i1)'
 write(fmt2,'(a,i7,a)') '(i6,', nmarq,'i4)'
 allocate(work1(nmarq))
if(algo==0 .or. algo==1)then
 open(51,file='imputed.dose')
 open(52,file='imputed.bestgeno')
 open(54,file='imputed.bestref')
 open(16,file='mosaic.fb')
  do l=1,nhap
    t=0
    if(mod(l,100)==0)print*,'Forward-Backward ::',l
    call forward_backward(l)
    forall(k=1:nmarq) work1(k)=refhap(topfb(k),k)-1
    write(54,fmt1) work1
    write(16,fmt2) l, topfb
  enddo
  close(54)
  close(16)
endif

if(algo==0 .or. algo==2)then
 open(14,file='mosaic.vtb')
 open(53,file='imputed.vtb')
  do l=1,nhap
    t=0
    if(mod(l,100)==0)print*,'Viterbi ::',l
    call viterbi(l)
  enddo

 do l=1, nhap
        forall(k=1:nmarq) work1(k)=refhap(states(k,l),k)-1
        write(53,fmt1) work1
        write(14,fmt2) l, states(:,l)
 enddo
 close(14)
 close(53)

endif

contains

! **************************************************************
! ********* Reading data (observations) and parameters *********
! **************************************************************

subroutine read_data
implicit none
integer ::nbr
real* 8 ::pld
character*9 ::option

algo=0;ngen=4.d0;gerr=0.001;Ne=0.d0
checkparam=0
do i=2,narg,2 !### take arguments by pairs
 info1=arguments(i)
 info2=adjustl(arguments(i+1))
 select case(info1)
 case("--refs")
  haplofile=adjustl(arguments(i+1))
  checkparam=checkparam+1
 case("--targets")
  haplofile2=adjustl(arguments(i+1))
  checkparam=checkparam+1
 case("--map")
  markfile=adjustl(arguments(i+1))
  checkparam=checkparam+1
 case("--ngen")
  read(info2,*)ngen
 case("--error")
  read(info2,*)gerr
 case("--algorithm")
  if(info2 .eq. 'FB')algo=1
  if(info2 .eq. 'VTB')algo=2
 case("--Ne")
  read(info2,*)Ne
 case default
  write(*,*)"Unrecognized option ::",info1
  stop
 end select
enddo

if(checkparam /= 3)then
 print*,'Missing information in command line!'
 print*,'Arguments --refs, --targets and --map are required with unrelated method.'
 stop
endif


if(algo==0)print*,"Running Forward-Backward and Viterbi algorithms"
if(algo==1)print*,"Running Forward-Backward algorithm"
if(algo==2)print*,"Running Viterbi algorithm"

open(50,file=markfile)

nmarq=0
do
        read(50,*,iostat=io)num
        if(io/=0)exit
        nmarq=nmarq+1
enddo

print*,'Number of markers ::',nmarq
rewind(50)


! **************************************
! ***** Reading data (observations) ****
! **************************************
nhap=0;nref=0

open(11,file=haplofile)

do 
read(11,*,iostat=io)num
if(io/=0)exit
nref=nref+1
enddo

rewind(11)

nclust=nref

print*,'Number of reference haplotypes ::',nref

open(12,file=haplofile2)

do 
read(12,*,iostat=io)num
if(io/=0)exit
nhap=nhap+1
enddo

rewind(12)

print*,'Number of target haplotypes ::',nhap

allocate(hap(nhap,nmarq),refhap(nref,nmarq),hapin(nmarq))
allocate(probrec(nmarq),posi(0:nmarq),nall(nmarq),states(nmarq,nhap))
allocate(topfb(nmarq))

hap=0;nall=0;refhap=0;nref=0

do
read(11,*,iostat=io)num,ori,(hapin(j),j=1,nmarq)
if(io/=0)exit
nref=nref+1
do i=1,nmarq
 nall(i)=max(nall(i),hapin(i))
enddo
refhap(nref,:)=hapin(:)
enddo

nhap=0
do
read(12,*,iostat=io)num,ori,(hapin(j),j=1,nmarq)
if(io/=0)exit
nhap=nhap+1
do i=1,nmarq
 nall(i)=max(nall(i),hapin(i))
enddo
hap(nhap,:)=hapin(:)
enddo

maxall=0
do i=1,nmarq
 maxall=max(maxall,nall(i))
enddo

print*,'End of reading haplotype files '

posi=0.d0

do
read(50,*,iostat=io)num,markfile,pld
if(io/=0)exit
if(num<=nmarq .and. num>0)posi(num)=pld/100.d0
enddo

print*,'Length of marker map (in cM / Mb) ::',(posi(nmarq)-posi(1))*100.0

probrec=0.0

do i=1,nmarq-1
  probrec(i)=0.5d0*(1.d0-dexp(-2.d0*abs(posi(i)-posi(i+1))))
enddo

end subroutine

!*********************************
!*** Parameters initialisation ***
!*********************************

subroutine parameter_init
implicit none
character*10 ::line
real*8 ::nG

allocate(pi(nclust,nmarq),aij(nmarq-1))

if(Ne > 0.d00)then
 ngen=4.d0*Ne/(1.d0*nref)
 theta=0.d00
 do i=1,(nref-1)
  theta=theta+1.d0/(1.d0*i)
 enddo
 theta=1.d0/theta
 gerr=theta/(2.d0*(theta+1.d0*nref))
endif

!nG=4.761905*1.d0
nG=ngen
do i=1,nmarq-1
!  aij(i)=1.00-(1.00-probrec(i))**nG
  aij(i)=1.d0-dexp(-abs(posi(i)-posi(i+1))*nG)
enddo

do i=1,nclust
 do j=1,nmarq
   pi(i,j)=1.d0/(1.d0*nclust)
 enddo
enddo

print*,'End initialisation of parameters, nG and eps :',nG,gerr


end subroutine

function emission(haplo,cluster1,marker)
implicit none
integer ::haplo,cluster1,marker,all1
real*8 ::emission,eps

!  eps=0.0005217585*1.d0
  eps=gerr
  emission=1.d0
  if(hap(haplo,marker)/=0 .and. refhap(cluster1,marker)/=0)then
    if(refhap(cluster1,marker)==hap(haplo,marker))emission=1.d0-eps*1.d0
    if(refhap(cluster1,marker)/=hap(haplo,marker))emission=eps*1.d0
    return
  endif

end function

function jump(marker)
implicit none
integer ::marker
real*8 ::jump(0:1)

! first indice 0 or 1 tells if the first haplotype jumped or not
! second indice for the second haplotype

jump(0)=(1.d0-1.d0*aij(marker))
jump(1)=1.d0*aij(marker)

end function


!******* forward ,algorithm - computation of alpha ******************
!******* backward algorithm - estimation of beta ******************
!******* estimation of gamma = alpha*beta *************************
!******* gamma = probability of state i,j at marker k *************

subroutine forward_backward(nh)
implicit none
integer ::nh,i,j,k,thr
real*8, allocatable ::scaling(:),alpha(:,:),beta(:,:),gamma(:,:)

allocate(alpha(nclust,nmarq),scaling(nmarq),beta(nclust,nmarq),gamma(nclust,nmarq))

!****************************
!***** forward algortihm ****
!****************************

alpha(:,:)=0.0;scaling(:)=0.0

! initialisation

do i=1,nclust
  alpha(i,1)=pi(i,1)*emission(nh,i,1)
  scaling(1)=scaling(1)+alpha(i,1)
enddo

 scaling(1)=1.0/scaling(1)
 alpha(:,1)=alpha(:,1)*scaling(1)

! induction: to get to i, two ways:
! 1) transition + jump into cluster i
! 2) no transition and in i at previous position
 
 do k=2,nmarq
  do i=1,nclust
   do j=1,nclust
    alpha(i,k)=alpha(i,k)+alpha(j,k-1)*aij(k-1)*pi(i,k)*emission(nh,i,k)
   enddo
    alpha(i,k)=alpha(i,k)+alpha(i,k-1)*(1.0-aij(k-1))*emission(nh,i,k)
    scaling(k)=scaling(k)+alpha(i,k)
  enddo
 scaling(k)=1.0/scaling(k)
 alpha(:,k)=alpha(:,k)*scaling(k)
 enddo  

!*****************************
!***** backward algortihm ****
!*****************************

beta(:,:)=0.0;gamma(:,:)=0.0

! initialisation

do i=1,nclust
 gamma(i,nmarq)=alpha(i,nmarq)*1.0  ! beta(i,j,nmarq)=1.0
 beta(i,nmarq)=1.0*scaling(nmarq)
enddo

! induction
! to arrive in k: with or without transition

do k=nmarq-1,1,-1
 do i=1,nclust
  do j=1,nclust
    beta(i,k)=beta(i,k)+aij(k)*pi(j,k+1)*emission(nh,j,k+1)*beta(j,k+1)
  enddo
  beta(i,k)=beta(i,k)+(1.0-aij(k))*emission(nh,i,k+1)*beta(i,k+1)
  beta(i,k)=beta(i,k)*scaling(k)
  gamma(i,k)=alpha(i,k)*beta(i,k)/scaling(k)
 enddo
enddo

!*******************************
!***** write results ***********
!*******************************

write(51,'(i6)',advance='no')nh
write(52,'(i6)',advance='no')nh

do k=1,nmarq
  topfb(k)=maxloc(gamma(:,k),1)
  dosage=0.d0
  do i=1,nclust
   dosage=dosage+gamma(i,k)*refhap(i,k) ! if no missing in reference
  enddo
  dosage=dosage-1.d00
  if(k<nmarq)write(51,'(1x,f9.6)',advance='no')dosage
  if(k==nmarq)write(51,'(1x,f9.6)')dosage
  bestallele=nint(dosage)
  if(k<nmarq)write(52,'(1x,i1)',advance='no')bestallele
  if(k==nmarq)write(52,'(1x,i1)')bestallele
enddo

deallocate(alpha,beta,gamma,scaling)

end subroutine

! ************************************************************************************************
! **** "optimal" sequence - choice: the single best state sequence (path) - Viterbi algorithm ****
! ************************************************************************************************

subroutine viterbi(nh)
implicit none
integer ::nh,i,j,k,k1,thr
integer, allocatable ::phi1(:,:)
real*8 ::prob_emission,pjump(0:1),pmax,pos(1)
real*8, allocatable ::delta(:,:), probs(:)

allocate(delta(nclust,nmarq),phi1(nclust,nmarq),probs(nclust))
delta=-1.d0;phi1=0

! initialisation
! phi1 = address i in t-1 which maximise path 

 do i=1,nclust
   prob_emission=emission(nh,i,1)
!   delta(i,1)=log(pi(i,1))+log(prob_emission)
   delta(i,1)=pi(i,1)*prob_emission
   phi1(i,1)=0;
 enddo

! recursion: search path which maximise probability to arrive in (i,j)
do k=2,nmarq
        pjump=jump(k-1)
        do i=1,nclust
                prob_emission=emission(nh,i,k)
                do k1=1,nclust
                        if (delta(k1,k-1)>0.d0) then
                                val=delta(k1,k-1)
                        else
                                val=1.d0*1E-20*1E-20*1E-20*1E-20*1E-20*1E-20
                        endif
                        if(k1/=i)then ! one jumps
                                val=val*pjump(1)*pi(i,k)
                        else if(k1==i)then ! 0 or 1 jumps
                                val=val*(pjump(0)+pjump(1)*pi(i,k))
                        endif
                        probs(k1)=val
                enddo
                !phi1(i,k)=findmax1(probs,nclust)
                phi1(i,k)=maxloc(probs,1)
                delta(i,k)=maxval(probs)*prob_emission
        enddo
enddo

! termination

pmax=maxval(delta(:,nmarq))
pos=maxloc(delta(:,nmarq))
states(nmarq,nh)=pos(1)

! path

do k=nmarq-1,1,-1
 states(k,nh)=phi1(states(k+1,nh),k+1)
enddo

deallocate(delta,phi1,probs)


end subroutine

function findmax1(x,l) result(p)
implicit none
integer:: i, l, p
real*8:: x(l), xm, mv
mv=maxval(x)
if (abs(mv)>0) then
        if (mv>0) then
                xm=mv-1E-12*mv*1.d0
                do i=1, l
                        if (x(i)>xm) exit
                enddo
        else
                xm=mv+1E-12*mv*1.d0
                do i=1, l
                        if (x(i)>xm) exit
                enddo
        endif
else
        xm=mv-1E-12
        do i=1, l
                if (x(i)>xm) exit
        enddo
endif
p=i
end function findmax1

end subroutine

end module

