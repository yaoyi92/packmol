!  
!  Written by Leandro Martínez, 2009-2011.
!  Copyright (c) 2009-2011, Leandro Martínez, Jose Mario Martinez,
!  Ernesto G. Birgin.
!  
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  
!
! subroutine movebad: Move the worst molecules to new
!            positions
!

subroutine movebad(n,x,fx,movefrac,movebadrandom,precision,seed,hasbad,movebadprint)

  use sizes
  use molpa
  use usegencan
  implicit none

  ! For flashsort
  integer :: indflash(maxatom), lflash(maxatom), mflash

  ! Internal variables
  integer :: n, i, j, icart, itype, iatom, imol, ilubar, ilugan, &
             ilubar2, ilugan2, nbad, seed, igood, ibad, nmove
  double precision :: x(nn), fx, fmol(maxatom), movefrac, rnd, &
                      precision, frac, radiuswork(maxatom)
  logical :: movebadprint, hasbad, movebadrandom

  if(movebadprint) write(*,*) ' Moving worst molecules ... ' 

  icart = natfix
  do itype = 1, ntype
    if(.not.comptype(itype)) then
      icart = icart + nmols(itype)*natoms(itype)
    else
      do imol = 1, nmols(itype)
        do iatom = 1, natoms(itype)
          icart = icart + 1
          fatom(icart) = 0.d0
        end do
      end do
    end if
  end do

  move = .true.
  if(movebadprint) write(*,*) ' Function value before moving molecules:',fx
  do i = 1, ntotat
    radiuswork(i) = radius(i)
    radius(i) = radius_ini(i)
  end do
  call feasy(x,fx)
  move = .false.

  ! Moving the worst molecules

  hasbad = .false. 
  icart = natfix
  do itype = 1, ntype
    if(.not.comptype(itype)) then
      icart = icart + nmols(itype)*natoms(itype)
    else

      ! Checking the function value for each molecule

      nbad = 0                                             
      do imol = 1, nmols(itype)
        fmol(imol) = 0.d0
        do iatom = 1, natoms(itype)
          icart = icart + 1
          fmol(imol) = fmol(imol) + fatom(icart)
        end do
        if(fmol(imol).gt.precision) then 
          hasbad = .true.
          nbad = nbad + 1
        end if
      end do
      frac = dfloat(nbad)/dfloat(nmols(itype))
      if(movebadprint) write(*,"( a,i9,a,f8.2,a )") &
      '  Type ',itype,' molecules with non-zero contributions:', &
                     100.d0*frac,'%'

      if(nbad.gt.0) then

        frac = dmin1(movefrac,frac)

        ! Ordering molecules from best to worst

        mflash = 1 + nmols(itype)/10
        call flash1(fmol,nmols(itype),lflash,mflash,indflash)   

        ! Moving molecules

        nmove = max0(int(nmols(itype)*frac),1)
        if(movebadprint) then
          write(*,"( a,i9,a,i9 )") '  Moving ',nmove,' molecules of type ',itype
          if ( movebadrandom ) then
            write(*,*) ' New positions will be aleatory (movebadrandom is set) '
          else
            write(*,*) ' New positions will be based on good molecules (movebadrandom is not set) '
          end if
        end if
        imol = 0
        do i = 1, itype - 1
          if(comptype(i)) imol = imol + nmols(i) 
        end do
        write(*,"( '  Moving:|0 ',tr39,'100%|' )")
        write(*,"( '         |',$)") 
        j = 0
        do i = 1, nmove
          ibad = nmols(itype) - i + 1 
          igood = int(rnd(seed)*nmols(itype)*frac) + 1
          ilubar = 3*(indflash(ibad)+imol-1)
          ilugan = 3*(indflash(ibad)+imol-1)+3*ntotmol
          ilubar2 = 3*(indflash(igood)+imol-1)
          ilugan2 = 3*(indflash(igood)+imol-1)+3*ntotmol
          if ( movebadrandom ) then
            x(ilubar+1) = sizemin(1) + rnd(seed)*(sizemax(1)-sizemin(1))
            x(ilubar+2) = sizemin(2) + rnd(seed)*(sizemax(2)-sizemin(2)) 
            x(ilubar+3) = sizemin(3) + rnd(seed)*(sizemax(3)-sizemin(3)) 
          else
            x(ilubar+1) = x(ilubar2+1) - 0.3*dmax(itype)+0.6*rnd(seed)*dmax(itype) 
            x(ilubar+2) = x(ilubar2+2) - 0.3*dmax(itype)+0.6*rnd(seed)*dmax(itype)
            x(ilubar+3) = x(ilubar2+3) - 0.3*dmax(itype)+0.6*rnd(seed)*dmax(itype)
          end if
          x(ilugan+1) = x(ilugan2+1)
          x(ilugan+2) = x(ilugan2+2)
          x(ilugan+3) = x(ilugan2+3)
          call restmol(itype,ilubar,n,x,fx,.true.,precision,seed)
          do while( j < 45.d0*i/nmove ) 
            write(*,"('*',$)")
            j = j + 1
          end do
        end do             
        write(*,"('|')")
      end if
    end if
  end do

  call feasy(x,fx)
  if(movebadprint) write(*,*) ' Function value after moving molecules:', fx
  do i = 1, ntotat
    radius(i) = radiuswork(i)
  end do

return
end subroutine movebad

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                             c
!     Subroutine Flash1                                       c
!     SORTS ARRAY A WITH N ELEMENTS BY USE OF INDEX VECTOR L  c
!     OF DIMENSION M WITH M ABOUT 0.1 N.                      c
!     Karl-Dietrich Neubert, FlashSort1 Algorithm             c
!     in  Dr. Dobb's Journal Feb.1998,p.123                   c
!                                                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine flash1 (A, N, L, M, ind)

      implicit none
      double precision :: a(*), anmin, c1, hold, flash
      integer :: L(*), ind(*), i, n, nmax, m, k, ihold, nmove, j, iflash
!     ============================ CLASS FORMATION ===== 


      do i = 1, n
      ind(i) = i
      end do

      ANMIN=A(1)
      NMAX=1 
      DO I=1,N
         IF( A(I).LT.ANMIN) ANMIN=A(I)
         IF( A(I).GT.A(NMAX)) NMAX=I
      END DO

      IF (ANMIN.EQ.A(NMAX)) RETURN
      C1=(M - 1) / (A(NMAX) - ANMIN)
      DO K=1,M  
         L(K)=0
      END DO 
      DO I=1,N
         K=1 + INT(C1 * (A(I) - ANMIN))
         L(K)=L(K) + 1
      END DO
      DO K=2,M
         L(K)=L(K) + L(K - 1)
      END DO
      HOLD=A(NMAX)
      A(NMAX)=A(1) 
      A(1)=HOLD

      ihold = ind(nmax)
      ind(nmax) = ind(1)
      ind(1) = ihold


!     =============================== PERMUTATION ===== 
      NMOVE=0 
      J=1
      K=M
      DO WHILE (NMOVE.LT.N - 1)
         DO WHILE (J.GT.L(K)) 
            J=J + 1 
            K=1 + INT(C1 * (A(J) - ANMIN)) 
         END DO  
         FLASH=A(J)
         iflash=ind(j)

         DO WHILE (.NOT.(J.EQ.L(K) + 1)) 
            K=1 + INT(C1 * (FLASH - ANMIN))
            HOLD=A(L(K)) 
            ihold = ind(L(k))
            A(L(K))=FLASH
            ind(L(k)) = iflash
            iflash = ihold
            FLASH=HOLD
            L(K)=L(K) - 1
            NMOVE=NMOVE + 1 
         END DO
      END DO

!     ========================= STRAIGHT INSERTION =====
      DO I=N-2,1,-1
         IF  (A(I + 1).LT.A(I)) THEN
            HOLD=A(I)
            ihold = ind(i)
            J=I
            DO WHILE  (A(J + 1).LT.HOLD)
               A(J)=A(J + 1)
               ind(j) = ind(j+1) 
               J=J + 1 
            END DO
            A(J)=HOLD 
            ind(j) = ihold
         ENDIF
      END DO

!     =========================== RETURN,END FLASH1 =====
      RETURN
      END               




























