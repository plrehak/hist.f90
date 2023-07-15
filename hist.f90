!       Pavel Rehak, October 2019
!       This program records the number of species (Ca2+, Citrate (if
!       applicable), Oxylate, and water)  in each cluster and
!       creates a histogram. In addition averages and stadndard
!       deviations are also calculated
        Program Main
        Implicit NONE

        integer :: i,j,Nf,ios,iu0,iu,iu2,ntot,dum,datnum,t
        integer :: clumax,iu3,ou2,upp
        real :: dum2,watnum,oxnum,upr
        integer,allocatable,dimension(:) :: clu,totnum
        real,allocatable,dimension(:) :: clurat,totrat
        character*18 :: title,typ,dum3
        character*18 ::output,input,input2,dumin,input3,out2
        
        datnum = 4
        dumin="Wat-cluster.xvg"
        iu0=10
        Nf=0
        ios=0
        open(file=dumin,unit=iu0,status='old')                
        read(iu0,*) title
        do while (ios==0)
            read(iu0,*,iostat=ios) dum!,dum2 ! count number of data points
!            upp = int(dum2)
!            if (upp > clumax) clumax=upp ! determine maximum cluster size
            Nf=Nf+1
        enddo
        Nf=Nf-1
!        clumax=clumax+1 ! pick up endpoints
        close(iu0)
         do t=1,datnum
            select case(t)
                case (1)
                   typ = "Wat" ! Calcium clusters
                case (2)
                   typ = "Ca" ! Water clusters
                case (3)
                   typ = "Ox" ! Oxylate clusters
                case (4)
                   typ = "Cit" ! Citrate Clusters
            end select
        
             input = trim(typ)//trim("-cluster.xvg")
!             output = trim(typ)//trim("-velz-r.txt")                
             iu = 11
!             ou = 20
        ! allocate arrays
          open(file=input,unit=iu,status='old')
          allocate(clu(Nf))
          allocate(clurat(Nf))
          clumax = 0
          upp = 0
          upr = 0.0
        read(iu,*) title ! read header
        ! read cluster sizes
        do i=1,Nf
         if(t==1) then
            read(iu,*) dum,dum2
         else
            read(iu,*) dum,dum2,dum3
         endif
          clu(i) = int(dum2)
          if (clu(i) > upp) upp=clu(i)
          if(t/=1)then
            if (dum3 == 'Inf') then
              clurat(i) = 0.0
            else
              read(dum3,'(f3.8)') clurat(i)
              if (clurat(i) > upr) upr= clurat(i)
            endif
          endif
        enddo
        clumax = max(upp,ceiling(upr))+1
        close(iu)

        call histograms(typ,clumax,Nf,clu,clurat)

        ! Do same calculations for combined Oxylate/Citrate species
        if (t==4) then
        ! Open data files
!          rewind(iu)
          iu2 = 21
          iu3 = 22
          ou2 = 32
          input2 = "Ox-cluster.xvg"
          input3 = "Wat-cluster.xvg"
          out2 = "Ox-Cit-cluster.xvg"
          open(file=input2,unit=iu2,status='old')
          open(file=input3,unit=iu3,status='old')
          read(iu2,*) title
          read(iu3,*) title
          open(file=out2,unit=ou2,status='unknown')
          write(ou2,*) "Time Cluster Ox+Cit  Wat:(Ox+Cit)"
! allocate necessary arrays
          allocate(totrat(Nf))
          allocate(totnum(Nf))
! Now determine ratio of Wat:(Citr+Oxyl)
          do i = 1, Nf
            read(iu2,*) dum2,oxnum
            read(iu3,*) dum2,watnum
            totnum(i) = clu(i)+int(oxnum)
            totrat(i)= watnum/real(totnum(i))
            write(ou2,*) dum2,totnum(i),totrat(i)
          enddo
        close(iu2)
        close(iu3)
        close(ou2)
        call histograms("Ox-Cit",clumax,Nf,totnum,totrat)
        endif
! Now we are finished deallocate arrays for next step
        deallocate(clu)
        deallocate(clurat)
!        deallocate(b)
!        deallocate(nclu)
!        deallocate(nclurat)
!        close(iu)
!        close(ou)
        enddo
!        close(iu0)
        END

      subroutine histograms(clutype,Nbins,Ndat,numarr,ratarr)
! This subroutine obtains arrays and places them in histograms and
! writes output files

        integer :: ou,Nbins,Ndat
        integer,dimension(Nbins+1):: b,nclu,nclurat
        integer,dimension(Ndat) :: numarr
        real,dimension(Ndat) :: ratarr
        character*18 :: clutype,output
        ! create bins
        do i=1,Nbins+1
          b(i) = i-1
          nclu(i)=0
          nclurat(i)=0
        enddo


! Check in which number bin cluster goes
        do i=1,Ndat
          do j=1,Nbins+1
             if (numarr(i)==b(j)) nclu(j) = nclu(j)+1
             if ((ratarr(i) >= real(b(j))) .and. &
                (ratarr(i) < real(b(j+1)))) nclurat(j) = nclurat(j)+1
          enddo
        enddo

! Write out histograms 
        output = trim(clutype)//trim("-hist.txt")
!        write(ou,*) clutype, "size Num-Clus Num-Wat:clus"
        ou = 2
        open(file=output,unit=ou,status='unknown')
        write(ou,*) clutype, "size Num-Clus Num-Wat:clus" 
       do i=1,Nbins+1
          if ((nclu(i) > 0.1) .or. (nclurat(i) > 0.1)) then
            write(ou,*) b(i),nclu(i),nclurat(i)
          else
            write(ou,*) b(i),0.0,0
          endif
        enddo
        end subroutine
