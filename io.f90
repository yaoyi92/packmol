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
! Subroutine getinp: subroutine that reads the input file
!

subroutine getinp(dism,precision,sidemax,&
                  ntype,nlines,nrest,&
                  natoms,idfirst,nconnect,maxcon,nmols,&
                  seed,&
                  discale,nloop,&
                  irestline,ityperest,linestrut,&
                  coor,amass,charge,restpars,&
                  pdbfile,name,ele,keyword,&
                  xyzout,writeout,writebad,&
                  tinker,pdb,xyz,moldy,check,chkgrad,&
                  randini,resnumbers,movefrac,movebadrandom,changechains,&
                  add_amber_ter,add_box_sides,add_sides_fix)

  use sizes
  use usegencan
  implicit none

  double precision :: dism, coor(maxatom,3), amass(maxatom),&
                      charge(maxatom), restpars(maxrest,9), clen,&
                      discale, movefrac, precision, add_sides_fix,&
                      sidemax

  integer :: i, j, k, ii, iarg, iline, nlines, idatom, iatom,&
             in, lixo, irest, itype, itest, seed, imark, ntype, nrest,&
             nmols(maxtype), natoms(maxtype), idfirst(maxtype),&
             ityperest(maxrest), irestline(maxrest),&
             linestrut(maxtype,2), maxcon(maxatom),&
             nconnect(maxatom,8), resnumbers(maxtype), charl, writeout,&
             nloop, ioerr
      
  character(len=200) :: inputfile(maxlines), keyword(maxlines,maxkeywords),&
                        record, blank, xyzout, name(maxtype)
  character(len=80) :: pdbfile(maxtype)
  character(len=3) :: ele(maxatom)

  logical :: tinker,pdb,xyz,moldy, ignore, randini, check, &
             chkgrad, writebad, add_amber_ter, add_box_sides, &
             changechains(maxtype), movebadrandom

  ! Clearing some character arrays

  do i = 1, 200
    do j = 1, 20
      keyword(i,j) = 'none'
    end do            
  end do
  do i = 1, 200
    blank(i:i) = ' '
  end do

  ! Instructions on how to run packmol

  write(*,*) ' Packmol must be run with: packmol < inputfile.inp '
  write(*,*)
  write(*,*) ' Userguide at: www.ime.unicamp.br/~martinez/packmol '
  write(*,*)
      
  ! Getting input lines from the input file

  write(*,*) ' Reading input file... (Control-C aborts)'
  nlines = 0
  do while (.true.)
    read(5,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit

    ! Checking if the line is empty or if it is a commentary

    i = 0
    ignore = .true.
    do while(ignore.and.i.lt.200)
      i = i + 1
      if(record(i:i).gt.' '.and.record(1:1).ne.'#') ignore = .false.
    end do

    ! If the line contains relevant data, save it in the inputfile array

    if(.not.ignore) then
      nlines = nlines + 1
      if ( nlines .gt. maxlines ) then
      write(*,*) ' ERROR: Input file too long.'
      write(*,*) ' Increase the maxlines parameter on sizes.f90 file.'
      stop
      end if
      inputfile(nlines) = record(1:200)
    end if
  
  end do

  ! Reading the keywords

  do iline = 1, nlines
    read(inputfile(iline),"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    call getkey(keyword,record,iline)
  end do

  ! Getting random seed and optional optimization parameters if set

  seed = 1234567
  randini = .false.
  check = .false.
  chkgrad = .false.
  iprint1 = 2
  iprint2 = 2
  discale = 1.1d0
  writeout = 10
  maxit = 20
  nloop = 0
  movefrac = 0.05
  movebadrandom = .false.
  precision = 1.d-2
  writebad = .false.
  add_amber_ter = .false.
  add_box_sides = .false.
  add_sides_fix = 0.d0
  sidemax = 1000.d0
  ioerr = 0
  do i = 1, nlines
    if(keyword(i,1).eq.'seed') then
      read(keyword(i,2),*,iostat=ioerr) seed
      if ( ioerr /= 0 ) exit
      if ( seed == -1 ) call seed_from_time(seed)
    else if(keyword(i,1).eq.'randominitialpoint') then
      randini = .true.
    else if(keyword(i,1).eq.'check') then
      check = .true.
    else if(keyword(i,1).eq.'writebad') then
      writebad = .true.
    else if(keyword(i,1).eq.'precision') then
      read(keyword(i,2),*,iostat=ioerr) precision
      if ( ioerr /= 0 ) exit
      write(*,*) ' Optional precision set: ', precision
    else if(keyword(i,1).eq.'movefrac') then
      read(keyword(i,2),*,iostat=ioerr) movefrac
      if ( ioerr /= 0 ) exit
      write(*,*) ' Optional movefrac set: ', movefrac
    else if(keyword(i,1).eq.'movebadrandom') then
      movebadrandom = .true.
      write(*,*) ' Will move randomly bad molecues (movebadrandom) '
    else if(keyword(i,1).eq.'chkgrad') then
      chkgrad = .true.
    else if(keyword(i,1).eq.'writeout') then
      read(keyword(i,2),*,iostat=ioerr) writeout
      if ( ioerr /= 0 ) exit
      write(*,*) ' Output frequency: ', writeout
    else if(keyword(i,1).eq.'maxit') then
      read(keyword(i,2),*,iostat=ioerr) maxit
      if ( ioerr /= 0 ) exit
      write(*,*) ' User defined GENCAN number of iterations: ', maxit
    else if(keyword(i,1).eq.'nloop') then
      read(keyword(i,2),*,iostat=ioerr) nloop
      if ( ioerr /= 0 ) exit
      write(*,*) ' User defined numer of GENCAN loops: ', nloop
    else if(keyword(i,1).eq.'discale') then
      read(keyword(i,2),*,iostat=ioerr) discale
      if ( ioerr /= 0 ) exit
      write(*,*) ' Optional initial tolerance scale: ', discale
    else if(keyword(i,1).eq.'sidemax') then
      read(keyword(i,2),*,iostat=ioerr) sidemax
      if ( ioerr /= 0 ) exit
      write(*,*) ' User set maximum system dimensions: ', sidemax
    else if(keyword(i,1).eq.'add_amber_ter') then
      add_amber_ter = .true.
      write(*,*) ' Will add the TER flag between molecules. '
    else if(keyword(i,1).eq.'add_box_sides') then
      add_box_sides = .true.
      write(*,*) ' Will print BOX SIDE informations. '
      read(keyword(i,2),*,iostat=ioerr) add_sides_fix
      if ( ioerr /= 0 ) cycle
      write(*,*) ' Will sum ', add_sides_fix,' to each side length on print'
    else if(keyword(i,1).eq.'iprint1') then 
      read(keyword(i,2),*,iostat=ioerr) iprint1
      if ( ioerr /= 0 ) exit
      write(*,*) ' Optional printvalue 1 set: ', iprint1
    else if(keyword(i,1).eq.'iprint2') then 
      read(keyword(i,2),*,iostat=ioerr) iprint2
      if ( ioerr /= 0 ) exit
      write(*,*) ' Optional printvalue 2 set: ', iprint2
    else if( keyword(i,1) /= 'tolerance' .and. &
             keyword(i,1) /= 'structure' .and. &
             keyword(i,1) /= 'end' .and. &
             keyword(i,1) /= 'atoms' .and. &
             keyword(i,1) /= 'output' .and. &
             keyword(i,1) /= 'filetype' .and. &
             keyword(i,1) /= 'number' .and. &
             keyword(i,1) /= 'inside' .and. &
             keyword(i,1) /= 'outside' .and. &
             keyword(i,1) /= 'fixed' .and. &
             keyword(i,1) /= 'center' .and. &
             keyword(i,1) /= 'centerofmass' .and. &
             keyword(i,1) /= 'over' .and. &
             keyword(i,1) /= 'below' .and. &
             keyword(i,1) /= 'constrain_rotation' .and. &
             keyword(i,1) /= 'radius' .and. &
             keyword(i,1) /= 'resnumbers' .and. &
             keyword(i,1) /= 'changechains' .and. &
             keyword(i,1) /= 'discale' .and. &
             keyword(i,1) /= 'maxit' .and. &
             keyword(i,1) /= 'movebadrandom' .and. &
             keyword(i,1) /= 'add_amber_ter' .and. &
             keyword(i,1) /= 'sidemax' .and. &
             keyword(i,1) /= 'seed' .and. &
             keyword(i,1) /= 'randominitialpoint' .and. &
             keyword(i,1) /= 'nloop' .and. &
             keyword(i,1) /= 'writeout' .and. &
             keyword(i,1) /= 'writebad' .and. &
             keyword(i,1) /= 'check' .and. &
             keyword(i,1) /= 'iprint1' .and. &
             keyword(i,1) /= 'iprint2' .and. &
             keyword(i,1) /= 'chkgrad' ) then
      write(*,*) ' ERROR: Keyword not recognized: ', trim(keyword(i,1))
      stop
    end if
  end do
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Some optional keyword was not used correctly: ', trim(keyword(i,1))
    stop
  end if
  write(*,*) ' Seed for random number generator: ', seed
  call init_random_number(seed)

  ! Searching for filetypes, default is pdb

  tinker = .false.
  pdb = .false.
  xyz = .false.
  moldy = .false.
  do i = 1, nlines
    if(keyword(i,1).eq.'filetype') then
      if(keyword(i,2).eq.'tinker') tinker = .true.
      if(keyword(i,2).eq.'pdb') pdb = .true.
      if(keyword(i,2).eq.'xyz') xyz = .true.
      if(keyword(i,2).eq.'moldy') moldy = .true.
    end if
  end do
  if(.not.pdb.and..not.tinker.and..not.xyz.and..not.moldy) then
    pdb = .true.
    write(*,*)' WARNING: File type not (correctly?) specified, using PDB'
  end if 

  ! Checking for the name of the output file to be created

  xyzout = '####'
  do iline = 1, nlines
    if(keyword(iline,1).eq.'output') then
      xyzout = keyword(iline,2)
      xyzout = xyzout(1:charl(xyzout))
    end if
  end do
  if(xyzout(1:4) == '####') then
    write(*,*)' ERROR: Output file not (correctly?) specified. '
    stop
  end if
  write(*,*)' Output file: ', xyzout(1:charl(xyzout))

  ! Reading structure files

  itype = 0
  do iline = 1, nlines
    if(keyword(iline,1).eq.'structure') then
      itype = itype + 1
     
      if(itype.gt.maxtype) then
        write(*,*) ' ERROR: Number of structures is greater than '
        write(*,*) '        maxtype. Please increase the maxtype '
        write(*,*) '        parameter on the sizes.i file and '
        write(*,*) '        recompile the program. '
        stop
      end if

      if(keyword(iline,2).eq.'none') then
        write(*,*) ' ERROR: structure without filename. '
        write(*,*) ' The syntax must be, for example: structure water.pdb '
        stop
      end if

      record = keyword(iline,2)
      write(*,*) ' Reading coordinate file: ', record(1:charl(record))

      ! Reading pdb input files

      if(pdb) then
        name(itype) = record(1:charl(record))
        record = keyword(iline,2)
        pdbfile(itype) = record(1:80)
        open(10,file=keyword(iline,2),status='old',iostat=ioerr)
        if ( ioerr /= 0 ) call failopen(keyword(iline,2))
        record(1:6) = '######'
        do while(record(1:4).ne.'ATOM'.and.record(1:6).ne.'HETATM')
          read(10,"( a200 )",iostat=ioerr) record
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Could not find any atom in PDB file. '
            stop
          end if
        end do
        natoms(itype) = 0
        do while(.true.)
          if(record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM') then
            natoms(itype)=natoms(itype) + 1
          end if
          read(10,"( a200 )",iostat=ioerr) record
          if ( ioerr /= 0 ) exit
        end do
        close(10)
        idfirst(itype) = 1
        do ii = itype - 1, 1, -1
          idfirst(itype) = idfirst(itype) + natoms(ii)
        end do
        open(10,file=keyword(iline,2),status='old',iostat=ioerr)
        if ( ioerr /= 0 ) call failopen(keyword(iline,2))
        record(1:6) = '######'
        do while(record(1:4).ne.'ATOM'.and.record(1:6).ne.'HETATM')
          read(10,"( a200 )") record
        end do
        idatom = idfirst(itype) - 1
        do while(idatom.lt.natoms(itype)+idfirst(itype)-1)
          if(record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM') then
            idatom = idatom + 1
            amass(idatom) = 1.d0
            read(record,"( t31,f8.3,t39,f8.3,t47,f8.3 )",iostat=ioerr) &
                 (coor(idatom,k),k=1,3)
            if( ioerr /= 0 ) then
              record = keyword(iline,2) 
              write(*,*) ' ERROR: Failed to read coordinates from', &
                         ' file: ', record(1:charl(record))
              write(*,*) ' Probably the coordinates are not in', &
                         ' standard PDB file format. '
              write(*,*) ' Standard PDB format specifications', &
                         ' can be found at: '
              write(*,*) ' www.rcsb.org/pdb '
              stop
            end if

            ! This only tests if residue numbers can be read, they are used 
            ! only for  output

            read(record(23:26),*,iostat=ioerr) itest
            if( ioerr /= 0 ) then
              record = pdbfile(itype)
              write(*,*) ' ERROR: Failed reading residue number',&
                         ' from PDB file: ',record(1:charl(record))
              write(*,*) ' Residue numbers are integers that',&
                         ' must be within columns 23 and 26. '
              write(*,*) ' Other characters within these columns',&
                         ' will cause input/output errors. '
              write(*,*) ' Standard PDB format specifications',&
                         ' can be found at: '
              write(*,*) ' www.rcsb.org/pdb '
              stop
            end if   
          end if
          read(10,"( a200 )",iostat=ioerr) record
          if ( ioerr /= 0 ) exit
        end do
        close(10)
      end if

      ! Reading tinker input files

      if(tinker) then
        open(10,file=keyword(iline,2),status='old',iostat=ioerr)
        if ( ioerr /= 0 ) call failopen(keyword(iline,2))
        idfirst(itype) = 1
        do ii = itype - 1, 1, -1
          idfirst(itype) = idfirst(itype) + natoms(ii)
        end do
        record = keyword(iline,2)
        call setcon(record(1:64),idfirst(itype),maxcon)
        open(10,file = keyword(iline,2), status = 'old')
        record = blank
        do while(record.le.blank)
          read(10,"( a200 )") record
        end do
        i = 1
        do while(record(i:i).le.' ')
          i = i + 1
        end do
        iarg = i
        do while(record(i:i).gt.' ')
          i = i + 1
        end do
        read(record(iarg:i-1),*) natoms(itype)
        do while(record(i:i).le.' ')
          i = i + 1
        end do
        iarg = i
        do while(record(i:i).gt.' ')
          i = i + 1
        end do
        read(record(iarg:i-1),"( a200 )") name(itype)
        record = name(itype)
        name(itype) = record(1:charl(record))
        if(name(itype).lt.' ') name(itype) = 'Without_title'
        idatom = idfirst(itype) - 1
        do iatom = 1, natoms(itype)
          idatom = idatom + 1
          record = blank
          do while(record.le.blank)
            read(10,"( a200 )") record
          end do
          i = 1
          do while(record(i:i).le.' ')
            i = i + 1
          end do
          iarg = i
          do while(record(i:i).gt.' ')
            i = i + 1
          end do
          read(record(iarg:i-1),*) in
          do while(record(i:i).le.' ')
            i = i + 1
          end do
          iarg = i
          do while(record(i:i).gt.' ')
            i = i + 1
          end do
          read(record(iarg:i-1),*) ele(idatom)    
          read(record(i:200),*) (coor(idatom,k), k = 1, 3),&
               (nconnect(idatom, k), k = 1, maxcon(idatom))
          amass(idatom) = 1.d0
        end do
        close(10)
      end if

      ! Reading xyz input files

      if(xyz) then
        open(10,file=keyword(iline,2),status='old',iostat=ioerr)
        if ( ioerr /= 0 ) call failopen(keyword(iline,2))
        read(10,*) natoms(itype)
        read(10,"( a200 )") name(itype)
        if(name(itype).lt.' ') name(itype) = 'Without_title'
        idfirst(itype) = 1
        do ii = itype - 1, 1, -1
          idfirst(itype) = idfirst(itype) + natoms(ii)
        end do
        idatom = idfirst(itype) - 1
        do iatom = 1, natoms(itype)
          idatom = idatom + 1
          record = blank
          read(10,"( a200 )") record
          read(record,*) ele(idatom), (coor(idatom,k),k=1,3)
          amass(idatom) = 1.d0
        end do
        close(10)
      end if
      
      ! Reading moldy input files
  
      if(moldy) then
        open(10,file=keyword(iline,2), status ='old',iostat=ioerr)
        if ( ioerr /= 0 ) call failopen(keyword(iline,2))
        read(10,*) name(itype), nmols(itype)
        natoms(itype) = 0
        do while(.true.)
          read(10,"( a200 )",iostat=ioerr) record
          if ( ioerr /= 0 ) exit
          if(record.gt.' '.and.record(1:3).ne.'end') & 
            natoms(itype) = natoms(itype) + 1
        end do
        close(10)
        idfirst(itype) = 1
        do ii = itype - 1, 1, -1
          idfirst(itype) = idfirst(itype) + natoms(ii)
        end do
        open(10,file=keyword(iline,2),status='old')
        read(10,"( a200 )") record
        idatom = idfirst(itype) - 1
        do iatom = 1, natoms(itype)
          idatom = idatom + 1
          read(10,"( a200 )") record
          read(record,*) lixo, (coor(idatom,k), k = 1, 3),&
                       amass(idatom), charge(idatom), ele(idatom)
        end do
        close(10)
      end if
    end if

  end do
  ntype = itype

  write(*,*) ' Number of independent structures: ', ntype
  write(*,*) ' The structures are: '

  do itype = 1, ntype
    record = name(itype)
    write(*,*) ' Structure ', itype, ':',record(1:charl(record)),&
               '(',natoms(itype),' atoms)'
  end do
  if(nloop.eq.0) nloop = 200*ntype
      
  ! Reading the restrictions that were set

  irest = 0
  ioerr = 0
  do iline = 1, nlines

    if(keyword(iline,1).eq.'fixed') then
      irest = irest + 1
      irestline(irest) = iline
      ityperest(irest) = 1
      read(keyword(iline,2),*,iostat=ioerr) restpars(irest,1)
      read(keyword(iline,3),*,iostat=ioerr) restpars(irest,2)
      read(keyword(iline,4),*,iostat=ioerr) restpars(irest,3)
      read(keyword(iline,5),*,iostat=ioerr) restpars(irest,4)
      read(keyword(iline,6),*,iostat=ioerr) restpars(irest,5)
      read(keyword(iline,7),*,iostat=ioerr) restpars(irest,6)
    end if
  
    if(keyword(iline,1).eq.'inside') then
      irest = irest + 1
      irestline(irest) = iline
      if(keyword(iline,2).eq.'cube') then
        ityperest(irest) = 2
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
      else if(keyword(iline,2).eq.'box') then
        ityperest(irest) = 3
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
        read(keyword(iline,7),*,iostat=ioerr) restpars(irest,5)
        read(keyword(iline,8),*,iostat=ioerr) restpars(irest,6)
      else if(keyword(iline,2).eq.'sphere') then
        ityperest(irest) = 4
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
      else if(keyword(iline,2).eq.'ellipsoid') then
        ityperest(irest) = 5
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
        read(keyword(iline,7),*,iostat=ioerr) restpars(irest,5)
        read(keyword(iline,8),*,iostat=ioerr) restpars(irest,6)
        read(keyword(iline,9),*,iostat=ioerr) restpars(irest,7)
      else if(keyword(iline,2).eq.'cylinder') then
        ityperest(irest) = 12
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1) 
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
        read(keyword(iline,7),*,iostat=ioerr) restpars(irest,5)
        read(keyword(iline,8),*,iostat=ioerr) restpars(irest,6)
        read(keyword(iline,9),*,iostat=ioerr) restpars(irest,7)  
        read(keyword(iline,10),*,iostat=ioerr) restpars(irest,9)  
        restpars(irest,8) = restpars(irest,4)**2 + &
                            restpars(irest,5)**2 + &
                            restpars(irest,6)**2
        if(restpars(irest,8).lt.1.d-10) then
           write(*,*) ' ERROR: The norm of the director vector', &
                      ' of the cylinder constraint cannot be zero.'
           ioerr = 1
        else
          clen = dsqrt(restpars(irest,8))
          restpars(irest,4) = restpars(irest,4) / clen
          restpars(irest,5) = restpars(irest,5) / clen
          restpars(irest,6) = restpars(irest,6) / clen
        end if                                     
      else
        ioerr = 1
      end if
    end if

    if(keyword(iline,1).eq.'outside') then
      irest = irest + 1
      irestline(irest) = iline
      if(keyword(iline,2).eq.'cube') then
        ityperest(irest) = 6
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
      else if(keyword(iline,2).eq.'box') then
        ityperest(irest) = 7
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
        read(keyword(iline,7),*,iostat=ioerr) restpars(irest,5)
        read(keyword(iline,8),*,iostat=ioerr) restpars(irest,6)
      else if(keyword(iline,2).eq.'sphere') then
        ityperest(irest) = 8
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
      else if(keyword(iline,2).eq.'ellipsoid') then
        ityperest(irest) = 9
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
        read(keyword(iline,7),*,iostat=ioerr) restpars(irest,5)
        read(keyword(iline,8),*,iostat=ioerr) restpars(irest,6)
        read(keyword(iline,9),*,iostat=ioerr) restpars(irest,7)
      else if(keyword(iline,2).eq.'cylinder') then
        ityperest(irest) = 13
        read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1) 
        read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
        read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
        read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
        read(keyword(iline,7),*,iostat=ioerr) restpars(irest,5)
        read(keyword(iline,8),*,iostat=ioerr) restpars(irest,6)
        read(keyword(iline,9),*,iostat=ioerr) restpars(irest,7)
        read(keyword(iline,10),*,iostat=ioerr) restpars(irest,9)
        restpars(irest,8) = restpars(irest,4)**2 + &
                            restpars(irest,5)**2 + &
                            restpars(irest,6)**2
        if(restpars(irest,8).lt.1.d-10) then
           write(*,*) ' ERROR: The norm of the director vector',&
                      ' of the cylinder constraint cannot be zero.'
           ioerr = 1
        else
          clen = dsqrt(restpars(irest,8))
          restpars(irest,4) = restpars(irest,4) / clen
          restpars(irest,5) = restpars(irest,5) / clen
          restpars(irest,6) = restpars(irest,6) / clen
        end if                                     
      else
        ioerr = 1
      end if
    end if

    if(keyword(iline,1).eq.'over') then
      irest = irest + 1
      irestline(irest) = iline
      ityperest(irest) = 10
      read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
      read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
      read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
      read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
      if(keyword(iline,2).ne.'plane') ioerr = 1
    end if

    if(keyword(iline,1).eq.'below') then
      irest = irest + 1
      irestline(irest) = iline
      ityperest(irest) = 11
      read(keyword(iline,3),*,iostat=ioerr) restpars(irest,1)
      read(keyword(iline,4),*,iostat=ioerr) restpars(irest,2)
      read(keyword(iline,5),*,iostat=ioerr) restpars(irest,3)
      read(keyword(iline,6),*,iostat=ioerr) restpars(irest,4)
      if(keyword(iline,2).ne.'plane') ioerr = 1 
    end if

    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Some restriction is not set correctly. '
      stop
    end if

  end do
  nrest = irest 
  write(*,*) ' Total number of restrictions: ', nrest

  ! Getting the tolerance

  ioerr = 1
  dism = -1.d0
  do iline = 1, nlines
    if(keyword(iline,1).eq.'tolerance') then
      read(keyword(iline,2),*,iostat=ioerr) dism
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Failed reading tolerance. '
        stop
      end if
      exit
    end if
  end do
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Overall tolerance not set. Use, for example: tolerance 2.0 '
    stop
  end if
  write(*,*) ' Distance tolerance: ', dism

  ! Assigning the input lines that correspond to each structure

  itype = 0
  iline = 0
  do while(iline.le.nlines)
    iline = iline + 1
    if(keyword(iline,1).eq.'structure') then
      itype = itype + 1
      linestrut(itype,1) = iline
      iline = iline + 1
      do while(keyword(iline,1).ne.'end'.or.&
               keyword(iline,2).ne.'structure')
        if(keyword(iline,1).eq.'structure'.or.&
           iline.eq.nlines) then
          write(*,*) ' Input file ERROR: structure specification',&
                     ' not ending with "end structure"'
          stop
        end if
        iline = iline + 1
      end do
      linestrut(itype,2) = iline 
    end if
  end do

  ! Getting the number of molecules of each type

  type: do itype = 1, ntype
    do iline = 1, nlines
      if(keyword(iline,1).eq.'number') then
        if(iline.gt.linestrut(itype,1).and.&
           iline.lt.linestrut(itype,2)) then
          read(keyword(iline,2),*,iostat=ioerr) nmols(itype)
          if ( ioerr /= 0 ) then
            write(*,*) ' ERROR: Error reading number of molecules of type ', itype
            stop
          else
            write(*,*) ' Number of molecules of type ', itype, ': ', nmols(itype)
            if(pdb.and.nmols(itype).gt.9999) then
              write(*,*) ' Warning: There will be more than 9999',&
                         ' molecules of type ',itype
              write(*,*) ' They will be divided into different chains in the PDB',&
                         ' output file. '
            end if
            cycle type
          end if 
        end if
      end if
    end do
    write(*,*) ' Warning: Number of molecules not set for type '&
               ,itype,': assuming 1 '
    nmols(itype) = 1
  end do type

  ! If pdb files, get the type of residue numbering output for each
  ! molecule

  if(pdb) then             
    do itype = 1, ntype
      resnumbers(itype) = -1
      changechains(itype) = .false.
      do iline = 1, nlines
        if(iline.gt.linestrut(itype,1).and.&
             iline.lt.linestrut(itype,2)) then
          if(keyword(iline,1).eq.'changechains') then
            changechains(itype) = .true.
          end if
          if(keyword(iline,1).eq.'resnumbers') then
            read(keyword(iline,2),*) resnumbers(itype)
          end if
        end if
      end do
      if ( resnumbers(itype) == -1 ) then
        write(*,*) ' Warning: Type of residue numbering not',&
                   ' set for structure ',itype
        call setrnum(pdbfile(itype),imark)
        if(imark.eq.1) resnumbers(itype) = 0
        if(imark.gt.1) resnumbers(itype) = 1 
      end if
      write(*,*) ' Residue numbering set for structure ',itype,':',&
                 resnumbers(itype)
      write(*,*) ' Swap chains of molecules of structure ',&
                 itype,':', changechains(itype) 
    end do
  end if

  return
end subroutine getinp

!
! Subroutine that stops if failed to open file
!

subroutine failopen(record)
  character(len=200) :: record
  write(*,*) 
  write(*,*) ' ERROR: Could not open file. '
  write(*,*) '        Could not find file: ',trim(record)
  write(*,*) '        Please check if all the input and structure ' 
  write(*,*) '        files are in the current directory or if the' 
  write(*,*) '        correct paths are provided.'
  write(*,*) 
  stop 
end subroutine failopen

!
! Subroutine that checks if a pdb structure has one or more than
! one residue
!

subroutine setrnum(file,nres)

  implicit none
  integer :: iread, ires, ireslast, nres, ioerr
  character(len=80) :: file
  character(len=200) :: record

  open(10,file=file,status='old')
  iread = 0
  nres = 1
  do while(nres.eq.1)
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if(record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM') then
      read(record(23:26),*,iostat=ioerr) ires
      if ( ioerr /= 0 ) cycle
      iread = iread + 1
      if(iread.gt.1) then
        if(ires.ne.ireslast) then
          nres = 2
          close(10)
          return
        end if
      end if
      ireslast = ires
    end if
  end do
  close(10)

  return
end subroutine setrnum

!
!  Subroutine that computes de number of connections of each atom
!  for tinker xyz files
! 

subroutine setcon(xyzfile, idfirst, maxcon)

  use sizes
  implicit none

  character(len=64) :: xyzfile
  character(len=120) :: record

  integer :: maxcon(maxatom)
  integer :: idfirst
  integer :: natoms, idatom, iatom, ic, i 

  open(10, file = xyzfile, status='old')
  read(10,*) natoms
  idatom = idfirst - 1
  do iatom = 1, natoms
    idatom = idatom + 1
    read(10,"( a120 )") record
    ic = 0
    do i = 1, 120
      if(record(i:i).gt.' '.and.record(i+1:i+1).le.' ') ic = ic + 1
    end do
    maxcon(idatom) = ic - 5
  end do
  close(10)

  return
end subroutine setcon
                
!
! Subroutine output: Subroutine that writes the output file
!
!
!

subroutine output(x,amass,&
                  irestline,linestrut,maxcon,ntcon,nconnect,&
                  nrest,ntfix,resnumbers,&
                  ele,pdbfile,xyzout,name,&
                  pdb,tinker,xyz,moldy,fix,&
                  add_amber_ter,add_box_sides,add_sides_fix,&
                  input_itype,thisisfixed,changechains)
 
  use sizes
  use molpa
  implicit none

  double precision :: x(nn)
  double precision :: tens(4,4), v(4,4), dv(4)
  double precision :: v1(3), v2(3), v3(3)
  double precision :: amass(maxatom)
  double precision :: xbar, ybar, zbar, beta, gama, teta, xcm, ycm, zcm
  double precision :: xlength, ylength, zlength
  double precision :: xxyx, xxyy, xxyz, xyyz, xyyy, xzyx,&
                   xzyy, xzyz, xyyx, xq, yq, zq, q0, q1, q2, q3
  double precision :: xtemp, ytemp, ztemp
  double precision :: sxmin, symin, szmin, sxmax, symax, szmax
  double precision :: add_sides_fix

  integer :: irestline(maxrest)
  integer :: linestrut(maxtype,2)
  integer :: nrest
  integer :: maxcon(maxatom)
  integer :: ntcon(9)
  integer :: nconnect(maxatom,8)
  integer :: k, ilugan, ilubar, itype, imol, idatom,&
          irest, iimol, ichain, iatom, irec, ilres, ifres,&
          iires, iconn, charl, irescount,&
          icart, i_ref_atom, ioerr
  integer :: ntfix, nr, nres, resnumbers(maxtype), imark  
  integer :: i_fixed, i_not_fixed, input_itype(maxtype)
  
  character :: chain, even_chain, odd_chain
  character(len=200) :: record
  character(len=64) :: title
  character(len=3) :: ele(maxatom)
  character(len=200) :: name(maxtype)
  character(len=80) :: pdbfile(maxtype), pdb_atom_line, pdb_hetatm_line,&
                       tinker_atom_line
  character(len=200) :: xyzout

  logical :: fix
  logical :: pdb, tinker, moldy, xyz
  logical :: add_amber_ter, add_box_sides
  logical :: thisisfixed(maxtype), changechains(maxtype)

  ! Job title

  title = ' Built with Packmol '

  ! Write the output (xyz file)

  if(xyz) then
    open(30,file=xyzout,status='unknown') 
    write(30,*) ntotat
    write(30,*) title 
    ilubar = 0 
    ilugan = ntotmol*3 
    icart = natfix
    i_not_fixed = 0
    i_fixed = ntype
    do itype = 1, ntfix
      if ( .not. thisisfixed(input_itype(itype)) ) then
        i_not_fixed = i_not_fixed + 1
        do imol = 1, nmols(i_not_fixed)
          xbar = x(ilubar+1) 
          ybar = x(ilubar+2) 
          zbar = x(ilubar+3)
          beta = x(ilugan+1)
          gama = x(ilugan+2)
          teta = x(ilugan+3)
          call eulerrmat(beta,gama,teta,v1,v2,v3)   
          idatom = idfirst(i_not_fixed) - 1      
          do iatom = 1, natoms(i_not_fixed) 
            icart = icart + 1
            idatom = idatom + 1
            call compcart(icart,xcart,&
                          xbar,ybar,zbar,&
                          coor(idatom,1),coor(idatom,2),&
                          coor(idatom,3),&
                          v1,v2,v3) 
            write(30,"( tr2,a3,tr2,3(tr2,f14.6) )") ele(idatom), (xcart(icart, k), k = 1, 3)
          end do 
          ilugan = ilugan + 3 
          ilubar = ilubar + 3 
        end do
      else
        i_fixed = i_fixed + 1
        idatom = idfirst(i_fixed) - 1
        do iatom = 1, natoms(i_fixed)
          idatom = idatom + 1
          write(30,"( tr2,a3,tr2,3(tr2,f14.6) )") ele(idatom), (coor(idatom,k),k=1,3)
        end do
      end if
    end do
    close(30)
  end if

  ! write the output as a MOLDY file

  if(moldy) then
    open(30,file=xyzout,status='unknown') 
    ! For square moldy boxes, this must be the side dimensions of the box 
    sxmin = 1.d30
    symin = 1.d30
    szmin = 1.d30
    sxmax = -1.d30
    symax = -1.d30
    szmax = -1.d30
    do irest = 1, nrest 
      if(ityperest(irest).eq.2) then
        sxmin = dmin1(restpars(irest,1),sxmin)
        symin = dmin1(restpars(irest,2),symin)
        szmin = dmin1(restpars(irest,3),szmin)
        sxmax = dmax1(restpars(irest,4)+restpars(irest,1),sxmax)
        symax = dmax1(restpars(irest,4)+restpars(irest,2),symax)
        szmax = dmax1(restpars(irest,4)+restpars(irest,3),szmax)
      else if(ityperest(irest).eq.3) then
        sxmin = dmin1(restpars(irest,1),sxmin)
        symin = dmin1(restpars(irest,2),symin)
        szmin = dmin1(restpars(irest,3),szmin)
        sxmax = dmax1(restpars(irest,4),sxmax)
        symax = dmax1(restpars(irest,5),symax)
        szmax = dmax1(restpars(irest,6),szmax)
      else
        write(*,*) ' WARNING: The first line of the moldy output'
        write(*,*) ' file contains the size of the sides of the'
        write(*,*) ' paralelogram that defines the system. '
        write(*,*) ' The numbers printed may not be correct in '
        write(*,*) ' this case because regions other than cubes '
        write(*,*) ' or boxes were used. '
        sxmin = dmin1(sxmin,sizemin(1))
        symin = dmin1(symin,sizemin(2))
        szmin = dmin1(szmin,sizemin(3))
        sxmax = dmax1(sxmax,sizemax(1))
        symax = dmax1(symax,sizemax(2))
        szmax = dmax1(szmax,sizemax(3))
      end if
    end do
    xlength = sxmax - sxmin
    ylength = symax - symin
    zlength = szmax - szmin
    write(30,"( 3(tr1,f12.6),' 90 90 90 1 1 1 ' )") xlength, ylength, zlength
    ilubar = 0 
    ilugan = ntotmol*3 
    i_not_fixed = 0
    i_fixed = ntype
    do itype = 1, ntfix
      if ( .not. thisisfixed(input_itype(itype)) ) then
        i_not_fixed = i_not_fixed + 1
        record = name(i_not_fixed)
        do imol = 1, nmols(i_not_fixed) 
          xbar = (x(ilubar+1) - sxmin) / xlength
          ybar = (x(ilubar+2) - symin) / ylength
          zbar = (x(ilubar+3) - szmin) / zlength
          beta = x(ilugan+1)
          gama = x(ilugan+2)
          teta = x(ilugan+3)
          call eulerrmat(beta,gama,teta,v1,v2,v3)  
 
          ! Computing cartesian coordinates and quaternions 
 
          xxyx = 0.d0
          xxyy = 0.d0
          xxyz = 0.d0
          xyyx = 0.d0
          xyyy = 0.d0
          xyyz = 0.d0
          xzyx = 0.d0
          xzyy = 0.d0
          xzyz = 0.d0 
          idatom = idfirst(i_not_fixed) - 1      
          do iatom = 1, natoms(i_not_fixed) 
            idatom = idatom + 1
            xq =   coor(idatom, 1)*v1(1) &
                 + coor(idatom, 2)*v2(1) &
                 + coor(idatom, 3)*v3(1)    
            yq =   coor(idatom, 1)*v1(2) &
                 + coor(idatom, 2)*v2(2) &
                 + coor(idatom, 3)*v3(2)    
            zq =   coor(idatom, 1)*v1(3) &
                 + coor(idatom, 2)*v2(3) &
                 + coor(idatom, 3)*v3(3)   

            ! Recovering quaternions for molecule imol

            xxyx = xxyx + xq * coor(idatom,1) * amass(idatom)
            xxyy = xxyy + xq * coor(idatom,2) * amass(idatom)
            xxyz = xxyz + xq * coor(idatom,3) * amass(idatom)
            xyyx = xyyx + yq * coor(idatom,1) * amass(idatom)
            xyyy = xyyy + yq * coor(idatom,2) * amass(idatom)
            xyyz = xyyz + yq * coor(idatom,3) * amass(idatom)
            xzyx = xzyx + zq * coor(idatom,1) * amass(idatom)
            xzyy = xzyy + zq * coor(idatom,2) * amass(idatom)
            xzyz = xzyz + zq * coor(idatom,3) * amass(idatom) 
          end do

          tens(1,1) = xxyx + xyyy + xzyz
          tens(1,2) = xzyy - xyyz
          tens(2,2) = xxyx - xyyy - xzyz
          tens(1,3) = xxyz - xzyx
          tens(2,3) = xxyy + xyyx
          tens(3,3) = xyyy - xzyz - xxyx
          tens(1,4) = xyyx - xxyy
          tens(2,4) = xzyx + xxyz
          tens(3,4) = xyyz + xzyy
          tens(4,4) = xzyz - xxyx - xyyy
          nr = 16
          call jacobi (tens, 4, 4, dv, v, nr)
          q0 = v(1,4)
          q1 = v(2,4)
          q2 = v(3,4)
          q3 = v(4,4)
          record = name(i_not_fixed) 
          xbar = dmin1(0.999999d0,xbar)     
          ybar = dmin1(0.999999d0,ybar)     
          zbar = dmin1(0.999999d0,zbar)     
          write(30,"( a10,tr1,7(f12.6) )") record(1:charl(record)), xbar, ybar, zbar, &
                         q0, q1, q2, q3
          ilugan = ilugan + 3 
          ilubar = ilubar + 3 
        end do
      else
        i_fixed = i_fixed + 1
        idatom = idfirst(i_fixed) - 1

        ! Getting the specified position of the molecule

        do irest = 1, nrest
          if(irestline(irest).gt.linestrut(i_fixed,1).and.&
             irestline(irest).lt.linestrut(i_fixed,2)) then
             xcm = restpars(irest,1) - sxmin
             ycm = restpars(irest,2) - symin
             zcm = restpars(irest,3) - szmin
             beta = -restpars(irest,4)
             gama = -restpars(irest,5)
             teta = -restpars(irest,6)
          end if
        end do
        call eulerrmat(beta,gama,teta,v1,v2,v3) 
  
        ! Computing cartesian coordinates and quaternions 
 
        xxyx = 0.d0
        xxyy = 0.d0
        xxyz = 0.d0
        xyyx = 0.d0
        xyyy = 0.d0
        xyyz = 0.d0
        xzyx = 0.d0
        xzyy = 0.d0
        xzyz = 0.d0 
        idatom = idfirst(i_fixed) - 1      
        do iatom = 1, natoms(i_fixed) 
          idatom = idatom + 1
          xtemp = coor(idatom,1) - xcm
          ytemp = coor(idatom,2) - ycm
          ztemp = coor(idatom,3) - zcm
          xq =   xtemp*v1(1) + ytemp*v2(1) + ztemp*v3(1)    
          yq =   xtemp*v1(2) + ytemp*v2(2) + ztemp*v3(2)    
          zq =   xtemp*v1(3) + ytemp*v2(3) + ztemp*v3(3)   
          xxyx = xxyx + xtemp * xq * amass(idatom)
          xxyy = xxyy + xtemp * yq * amass(idatom)
          xxyz = xxyz + xtemp * zq * amass(idatom)
          xyyx = xyyx + ytemp * xq * amass(idatom)
          xyyy = xyyy + ytemp * yq * amass(idatom)
          xyyz = xyyz + ytemp * zq * amass(idatom)
          xzyx = xzyx + ztemp * xq * amass(idatom)
          xzyy = xzyy + ztemp * yq * amass(idatom)
          xzyz = xzyz + ztemp * zq * amass(idatom)
        end do
        tens(1,1) = xxyx + xyyy + xzyz
        tens(1,2) = xzyy - xyyz
        tens(2,2) = xxyx - xyyy - xzyz
        tens(1,3) = xxyz - xzyx
        tens(2,3) = xxyy + xyyx
        tens(3,3) = xyyy - xzyz - xxyx
        tens(1,4) = xyyx - xxyy
        tens(2,4) = xzyx + xxyz
        tens(3,4) = xyyz + xzyy
        tens(4,4) = xzyz - xxyx - xyyy
        nr = 16
        call jacobi (tens, 4, 4, dv, v, nr)
        q0 = v(1,4)
        q1 = v(2,4)
        q2 = v(3,4)
        q3 = v(4,4)
        xcm = xcm / xlength
        ycm = ycm / ylength
        zcm = zcm / zlength
        record = name(itype)
        xcm = dmin1(0.999999d0,xcm)     
        ycm = dmin1(0.999999d0,ycm)     
        zcm = dmin1(0.999999d0,zcm)     
        write(30,"( a10,tr1,7(f12.6) )") record(1:charl(record)),&
                       xcm, ycm, zcm, q0, q1, q2, q3
      end if
    end do
    close(30) 
  end if 

  ! write the output as pdb file

  if(pdb) then

    pdb_atom_line = "( t1,a5,t7,i5,t12,a10,t22,a1,t23,&
                       i4,t27,a1,t31,f8.3,t39,f8.3,t47,&
                       f8.3,t55,a26 )"
    pdb_hetatm_line = "( t1,a6,t7,i5,t12,a10,t22,a1,&
                         t23,i4,t27,a1,t31,f8.3,t39,&
                         f8.3,t47,f8.3,t55,a26 )"

    open(30,file=xyzout,status='unknown') 
 
    write(30,"( & 
             'HEADER ',/&
             'TITLE    ', a64,/&
             'REMARK   Packmol generated pdb file ',/&
             'REMARK   Home-Page: ',&
             'http://www.ime.unicamp.br/~martinez/packmol',/,&
             'REMARK' )" ) title

    if(add_box_sides) then
      write(30,"( 'CRYST1',t7,f9.2,t16,f9.2,t25,f9.2,&
                      t34,f7.2,t41,f7.2,t48,f7.2,&
                      t56,'P 1           1' )") &
            sizemax(1)-sizemin(1) + add_sides_fix,&
            sizemax(2)-sizemin(2) + add_sides_fix,& 
            sizemax(3)-sizemin(3) + add_sides_fix,&
            90., 90., 90.
    end if
      
    ilubar = 0 
    ilugan = ntotmol*3 
    icart = natfix
    i_ref_atom = 0
    iimol = 0
    ichain = 0
    i_not_fixed = 0
    i_fixed = ntype
    irescount = 1
    do itype = 1, ntfix 
      if ( .not. thisisfixed(input_itype(itype)) ) then
        i_not_fixed = i_not_fixed + 1

        ! Counting the number of residues of this molecule

        open(15,file=pdbfile(i_not_fixed),status='old')
        ifres = 0
        do
          read(15,"( a80 )",iostat=ioerr) record
          if ( ioerr /= 0 ) exit
          if ( record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM' ) then
            read(record(23:26),*,iostat=ioerr) imark
            if ( ioerr /= 0 ) then
              record = pdbfile(i_not_fixed)
              write(*,*) ' ERROR: Failed reading residue number ',&
                         ' from PDB file: ', record(1:charl(record))
              write(*,*) ' Residue numbers are integers that must',&
                         ' be between columns 23 and 26. '
              write(*,*) ' Other characters within these columns',&
                         ' will cause input/output errors. '
              write(*,*) ' Standard PDB format specifications can',&
                         ' be found at: '
              write(*,*) ' www.rcsb.org/pdb '
              stop
            end if
            if ( ifres .eq. 0 ) ifres = imark
            ilres = imark
          end if
        end do
        nres = ilres - ifres + 1

        do irec = 1, 80
          record(irec:irec) = ' '
        end do

        mol: do imol = 1, nmols(i_not_fixed) 
          iimol = iimol + 1

          rewind(15)
          read(15,"( a80 )",iostat=ioerr) record
          do while(record(1:4).ne.'ATOM'.and.record(1:6).ne.'HETATM')
            read(15,"( a80 )",iostat=ioerr) record
          end do

          if(imol.eq.1.or.mod(imol,9999).eq.1) then
            ichain = ichain + 1
            if( changechains(i_not_fixed) ) then
              call chainc(ichain,odd_chain)
              ichain = ichain + 1
              call chainc(ichain,even_chain)
            else 
              call chainc(ichain,even_chain)
              odd_chain = even_chain
            end if
          end if
          if ( mod(imol,2) == 0 ) chain = even_chain
          if ( mod(imol,2) /= 0 ) chain = odd_chain

          xbar = x(ilubar+1) 
          ybar = x(ilubar+2) 
          zbar = x(ilubar+3) 
          beta = x(ilugan+1)
          gama = x(ilugan+2)
          teta = x(ilugan+3)

          call eulerrmat(beta,gama,teta,v1,v2,v3)

          idatom = idfirst(i_not_fixed) - 1     
          iatom = 0
          do while(iatom.lt.natoms(i_not_fixed))
            iatom = iatom + 1 
            icart = icart + 1
            idatom = idatom + 1
            i_ref_atom = i_ref_atom + 1

            call compcart(icart,xcart,&
                          xbar,ybar,zbar,&
                          coor(idatom,1),coor(idatom,2),&
                          coor(idatom,3),v1,v2,v3)

            if(iatom.gt.1) then
              read(15,"( a80 )",iostat=ioerr) record
              if ( ioerr /= 0 ) exit mol
            end if 
            if(record(1:4).ne.'ATOM'.and.record(1:6).ne.'HETATM') then
              iatom = iatom - 1
              icart = icart - 1
              idatom = idatom - 1
              i_ref_atom = i_ref_atom - 1
              cycle
            end if

            ! Setting residue numbers for this molecule

            imark = 0
            read(record(23:26),*,iostat=ioerr) imark
            if ( ioerr /= 0 ) imark = 1
            if(resnumbers(i_not_fixed).eq.0) then
              iires = mod(imol,9999)
            else if(resnumbers(i_not_fixed).eq.1) then
              iires = imark
            else if(resnumbers(i_not_fixed).eq.2) then
              iires = mod(imark-ifres+irescount,9999)
            else if(resnumbers(i_not_fixed).eq.3) then
              iires = mod(iimol,9999)
            end if
            if(iires.eq.0) iires = 9999

            ! Writing output line

            if(record(1:4).eq.'ATOM') then
              write(30, pdb_atom_line) record(1:5), i_ref_atom,&
                                       record(12:21), chain, iires,&
                                       record(27:27),&
                                       (xcart(icart,k), k = 1, 3),&
                                       record(55:80)
            end if

            if(record(1:6).eq.'HETATM') then
               write(30,pdb_hetatm_line) record(1:6), i_ref_atom,&
                                         record(12:21), chain, iires,&
                                         record(27:27),&
                                         (xcart(icart,k), k = 1, 3),&
                                         record(55:80)
            end if
          end do
          irescount = irescount + nres
          ilugan = ilugan + 3 
          ilubar = ilubar + 3 
 
          if(add_amber_ter) write(30,"('TER')") 
        end do mol
        close(15)

      else
        i_fixed = i_fixed + 1

        ! Counting the number of residues of this molecule

        open(15,file=pdbfile(i_fixed),status='old')
        ifres = 0
        do
          read(15,"( a80 )",iostat=ioerr) record
          if ( ioerr /= 0 ) exit
          if ( record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM' ) then
            read(record(23:26),*,iostat=ioerr) imark
            if ( ioerr /= 0 ) then
              record = pdbfile(i_not_fixed)
              write(*,*) ' ERROR: Failed reading residue number ',&
                         ' from PDB file: ', record(1:charl(record))
              write(*,*) ' Residue numbers are integers that must',&
                         ' be between columns 23 and 26. ' 
              write(*,*) ' Other characters within these columns',&
                         ' will cause input/output errors. ' 
              write(*,*) ' Standard PDB format specifications can',&
                         ' be found at: '
              write(*,*) ' www.rcsb.org/pdb '
              stop
            end if
            if ( ifres .eq. 0 ) ifres = imark
            ilres = imark
          end if
        end do
        nres = ilres - ifres + 1

        iimol = iimol + 1
        idatom = idfirst(i_fixed) - 1

        rewind(15)
        read(15,"( a80 )",iostat=ioerr) record
        do while(record(1:4).ne.'ATOM'.and.record(1:6).ne.'HETATM')
          read(15,"( a80 )",iostat=ioerr) record
        end do

        iatom = 0
        do while(iatom.lt.natoms(i_fixed))
          iatom = iatom + 1
          idatom = idatom + 1
          i_ref_atom = i_ref_atom + 1

          if(iatom.gt.1) then
            read(15,"( a80 )",iostat=ioerr) record
            if ( ioerr /= 0 ) exit
          end if

          if(record(1:4).ne.'ATOM'.and.record(1:6).ne.'HETATM') then
            write(30,"( a80 )") record(1:80)
            iatom = iatom - 1
            icart = icart - 1
            idatom = idatom - 1
            i_ref_atom = i_ref_atom - 1
            cycle
          end if    

          read(record(23:26),*) imark
          if(resnumbers(i_fixed).eq.0) then
            iires = 1
          else if(resnumbers(i_fixed).eq.1) then
            iires = imark
          else if(resnumbers(i_fixed).eq.2) then
            iires = mod(imark-ifres+irescount,9999) 
          else if(resnumbers(i_fixed).eq.3) then
            iires = mod(iimol,9999)
          end if

          if(record(1:4).eq.'ATOM') then
            write(30,pdb_atom_line) record(1:5), i_ref_atom,&
                                    record(12:21), record(22:22), iires,&
                                    record(27:27),&
                                    (coor(idatom,k), k = 1, 3),&
                                    record(55:80)
          end if

          if(record(1:6).eq.'HETATM') then
            write(30,pdb_hetatm_line) record(1:6), i_ref_atom,&
                                      record(12:21), record(22:22), iires,&
                                      record(27:27),&
                                      (coor(idatom,k), k = 1, 3),&
                                      record(55:80)
          end if

        end do
        irescount = irescount + nres
        close(15)
        if(add_amber_ter) write(30,"('TER')") 
      end if
    end do             
    write(30,"('END')")
    close(30) 
  end if 
 
  ! Write the output (tinker xyz file)

  if(tinker) then

    tinker_atom_line = "( i7,tr2,a3,3(tr2,f10.6),9(tr2,i7) )"

    open(30, file = xyzout,status='unknown') 
 
    write(30,"( i6,tr2,a64 )") ntotat, title 

    ilubar = 0 
    ilugan = ntotmol*3 
    icart = natfix
    i_ref_atom = 0
    i_not_fixed = 0
    i_fixed = ntype

    do itype = 1, ntfix

      if ( .not. thisisfixed(input_itype(itype)) ) then
        i_not_fixed = i_not_fixed + 1
  
        do imol = 1, nmols(i_not_fixed) 
 
          xbar = x(ilubar+1) 
          ybar = x(ilubar+2) 
          zbar = x(ilubar+3) 
          beta = x(ilugan+1)
          gama = x(ilugan+2)
          teta = x(ilugan+3)

          call eulerrmat(beta,gama,teta,v1,v2,v3) 

          iconn = icart
 
          idatom = idfirst(i_not_fixed) - 1      
          do iatom = 1, natoms(i_not_fixed) 

            icart = icart + 1
            idatom = idatom + 1
            i_ref_atom = i_ref_atom + 1

            call compcart(icart,xcart,&
                          xbar,ybar,zbar,&
                          coor(idatom,1),coor(idatom,2),&
                          coor(idatom,3),&
                          v1,v2,v3)    

            ntcon(1) = nconnect(idatom, 1)
            do k = 2, maxcon(idatom)
              ntcon(k) = nconnect(idatom, k) + iconn 
            end do
  
            write(30,tinker_atom_line) i_ref_atom,&
                                       ele(idatom), (xcart(icart, k), k = 1, 3),&
                                       ntcon(1),&
                                       (ntcon(k) - natfix, k = 2, maxcon(idatom))
          end do 
 
          ilugan = ilugan + 3 
          ilubar = ilubar + 3 
 
        end do 

      else 

        i_fixed = i_fixed + 1
        idatom = idfirst(i_fixed) - 1

        do iatom = 1, natoms(itype)
          idatom = idatom + 1
          i_ref_atom = i_ref_atom + 1
          ntcon(1) = nconnect(iatom,1)
          do k = 2, maxcon(idatom)
            ntcon(k) = nconnect(idatom,k) + i_ref_atom
          end do
        end do
        write(30,tinker_atom_line) i_ref_atom, ele(idatom),&
                                   (coor(idatom,k), k = 1, 3),&
                                   (ntcon(k), k = 1, maxcon(idatom))
      end if
    end do             
 
    close(30) 
  end if   

  return
end subroutine output

!
! Subroutine that writes the last point obtained when
! a solution was not found
!

subroutine checkpoint(n,x,amass,&
                      nrest,ntfix,nloop,&
                      irestline,linestrut,maxcon,ntcon,nconnect,&
                      ele,pdbfile,xyzout,name,&
                      pdb,tinker,xyz,moldy,fix,&
                      movefrac,movebadrandom,precision,seed,resnumbers,&
                      add_amber_ter,add_box_sides,add_sides_fix,&
                      input_itype,thisisfixed,changechains)

  use sizes
  use molpa
  use usegencan
  implicit none

  double precision :: x(*)
  double precision :: fx
  double precision :: amass(maxatom)
  double precision :: movefrac, precision, add_sides_fix

  integer :: irestline(maxrest)
  integer :: linestrut(maxtype,2)
  integer :: maxcon(maxatom)
  integer :: ntcon(9)
  integer :: nconnect(maxatom,8)
  integer :: i, charl
  integer :: n, nrest, ntfix
  integer :: resnumbers(maxtype)
  integer :: nloop, seed
  integer :: input_itype(maxtype)

  character(len=3) :: ele(maxatom)
  character(len=80) :: pdbfile(maxtype) 
  character(len=200) :: xyzout
  character(len=200) :: name(maxtype)
  character(len=80) :: dash1_line, dash2_line

  logical :: pdb, tinker, xyz, moldy, fix
  logical :: add_amber_ter, add_box_sides
  logical :: movebadprint, hasbad
  logical :: thisisfixed(maxtype), changechains(maxtype)
  logical :: movebadrandom

  dash1_line = "( 62('-'),/ )"
  dash2_line = "( /,62('#'),/ )"

  ! All molecules are important

   do i = 1, ntfix
     comptype(i) = .true.
   end do

  ! Call the subroutine that computes de function value

  call feasy(x,fx)

  write(*,"(/,62('-'),/,&
         ' Packmol was not able to find a solution to your',/,&
         ' packing problem with the desired distance tolerance.',/,/,&
         ' First of all, be sure if the molecules fit in the',/,&
         ' regions specified and if the constraints were set',/,&
         ' correctly. ',/,/,&
         ' Secondly, try simply running it again with a different ',/,&
         ' seed for the random number generator of the initial ',/,&
         ' point. This is done by adding the keyword seed to the',/,&
         ' input file, as in: ',/,/,&
         ' seed 192911 ',/,/,&
         ' The best configuration found has a function value of',/,&
         ' f = ', e14.7,/,/,&
         ' IMPORTANT: ',/,&
         ' If the number of molecules and the restraints are',/,&
         ' correct, it is still very likely that the current point',/,&
         ' fits your needs if your purpose is to run a MD',/,&
         ' simulation.',/,&
         ' Therefore, we recommend to minimize the energy of the',/,&
         ' solution found, equilibrate it and run with it as well.',/,/,&
         62('-'),/)") fx

  call output(x,amass,&
              irestline,linestrut,maxcon,ntcon,nconnect,&
              nrest,ntfix,resnumbers,&
              ele,pdbfile,xyzout,name,&
              pdb,tinker,xyz,moldy,fix,&
              add_amber_ter,add_box_sides,add_sides_fix,&
              input_itype,thisisfixed,changechains)

  write(*,*) ' The solution with the best function value was '
  write(*,*) ' written to the output file: ', xyzout(1:charl(xyzout))
  write(*,dash1_line)
  write(*,*) ' Forcing the solution to fit the constraints...'

  ! CALL GENCAN

  init1 = .true.
  do i = 1, nloop
    iprint1 = 0
    iprint2 = 0
    call pgencan(n,x,fx)
    movebadprint = .false.
    call movebad(n,x,fx,movefrac,movebadrandom,precision,seed,hasbad,movebadprint) 
  end do
  init1 = .false.

  write(*,*)
  write(*,dash1_line)
  xyzout = xyzout(1:charl(xyzout))//'_FORCED'
  call output(x,amass,&
              irestline,linestrut,maxcon,ntcon,nconnect,&
              nrest,ntfix,resnumbers,&
              ele,pdbfile,xyzout,name,&
              pdb,tinker,xyz,moldy,fix,&
              add_amber_ter,add_box_sides,add_sides_fix,&
              input_itype,thisisfixed,changechains)

  write(*,*) ' The forced point was writen to the '
  write(*,*) ' output file: ', xyzout(1:charl(xyzout)+7)
  write(*,*)
  write(*,*) ' If you want that the packing procedure continues'
  write(*,*) ' for a longer time, add the following keyword '
  write(*,*) ' to the input file: '
  write(*,*)
  write(*,*) ' nloop [integer]      (ex: nloop 200) '
  write(*,*)
  write(*,*) ' The default nloop value is 50 for each molecule.'
  write(*,*)

  write(*,dash2_line) 
  write(*,*) ' TERMINATION WITHOUT PERFECT PACKING: '
  write(*,*) ' The output file:'
  write(*,*)
  write(*,*) '   ',xyzout(1:charl(xyzout)-7) 
  write(*,*)
  write(*,*) ' contains the best solution found. '
  write(*,*)
  write(*,*) ' Very likely, if the input data was correct, '
  write(*,*) ' it is a reasonable starting configuration.'
  write(*,*) ' Check commentaries above for more details. '
  write(*,dash2_line) 
      
  return
end subroutine checkpoint

!
! function charl: set the length of a character string
! 

function charl(file)

  character(len=80) :: file
  integer :: charl

  charl = 1
  do while(file(charl:charl).le.' ')
    charl = charl + 1
  end do
  do while(file(charl:charl).gt.' ')
    charl = charl + 1
  end do
  charl = charl - 1

  return
end function charl

!
! Subroutine getkey: gets keywords and values from the input
!                    file in a more robust way
!

subroutine getkey(keyword,record,iline)

  use sizes
  implicit none
  character(len=200) :: keyword(maxlines,maxkeywords), record
  integer :: iline, i, ilast, ival

  i = 0
  ival = 0
  do while(i.le.200)
    i = i + 1
    ilast = i
    do while(record(i:i).gt.' '.and.i.lt.200)
      i = i + 1
    end do
    if(i.gt.ilast) then
      ival = ival + 1
      keyword(iline,ival) = record(ilast:i)
    end if
  end do

  return
end subroutine getkey

subroutine chainc(i,chain)

  implicit none
  integer :: i
  character :: chain

  if(i.eq.0) chain = ' '
  if(i.eq.1) chain = 'A'
  if(i.eq.2) chain = 'B'
  if(i.eq.3) chain = 'C'
  if(i.eq.4) chain = 'D'
  if(i.eq.5) chain = 'E'
  if(i.eq.6) chain = 'F'
  if(i.eq.7) chain = 'G'
  if(i.eq.8) chain = 'H'
  if(i.eq.9) chain = 'I'
  if(i.eq.10) chain = 'J'
  if(i.eq.11) chain = 'K'
  if(i.eq.12) chain = 'L'
  if(i.eq.13) chain = 'M'
  if(i.eq.14) chain = 'N'
  if(i.eq.15) chain = 'O'
  if(i.eq.16) chain = 'P'
  if(i.eq.17) chain = 'Q'
  if(i.eq.18) chain = 'R'
  if(i.eq.19) chain = 'S'
  if(i.eq.20) chain = 'T'
  if(i.eq.21) chain = 'U'
  if(i.eq.22) chain = 'V'
  if(i.eq.23) chain = 'W'
  if(i.eq.24) chain = 'X'
  if(i.eq.25) chain = 'Y'
  if(i.eq.26) chain = 'Z'
  if(i.eq.27) chain = '1'
  if(i.eq.28) chain = '2'
  if(i.eq.29) chain = '3'
  if(i.eq.30) chain = '4'
  if(i.eq.31) chain = '5'
  if(i.eq.32) chain = '6'
  if(i.eq.33) chain = '7'
  if(i.eq.34) chain = '8'
  if(i.eq.35) chain = '9'
  if(i.eq.36) chain = '0'
  if(i.gt.36) chain = '#'

  return
end subroutine chainc

subroutine clear(record)
      
  integer :: i
  character(len=80) :: record

  do i = 1, 80
    record(i:i) = ' '
  end do
      
  return
end subroutine clear

!
! JACOBI
! Jacobi diagonalizer with sorted output.  Same calling sequence as
! EISPACK routine, but must specify nrot!
!
      SUBROUTINE jacobi (a, n, np, d, v, nrot)
      IMPLICIT CHARACTER (A-Z)
!
      INTEGER n, np, nrot
      DOUBLEPRECISION a (np, n)
      DOUBLEPRECISION d (n)
      DOUBLEPRECISION v (np, n)
!
      DOUBLEPRECISION onorm, dnorm
      DOUBLEPRECISION b, dma, q, t, c, s
      DOUBLEPRECISION atemp, vtemp, dtemp
      INTEGER i, j, k, l
!
      DO 10000 j = 1, n
        DO 10010 i = 1, n
          v (i, j) = 0.0D0
10010     CONTINUE
        v (j, j) = 1.0D0
        d (j) = a (j, j)
10000   CONTINUE
!
      DO 20000 l = 1, nrot
        dnorm = 0.0D0
        onorm = 0.0D0
        DO 20100 j = 1, n
          dnorm = dnorm + ABS (d (j))
          DO 20110 i = 1, j - 1
            onorm = onorm + ABS (a (i, j))
20110       CONTINUE
20100     CONTINUE
        IF (onorm / dnorm .LE. 0.0D0) GOTO 19999
        DO 21000 j = 2, n
          DO 21010 i = 1, j - 1
            b = a (i, j)
            IF (ABS (b) .GT. 0.0D0) THEN
              dma = d (j) - d (i)
              IF (ABS (dma) + ABS (b) .LE. ABS (dma)) THEN
                t = b / dma
              ELSE
                q = 0.5D0 * dma / b
                t = SIGN (1.0D0 / (ABS (q) + SQRT (1.0D0 + q * q)), q)
              ENDIF
              c = 1.0D0 / SQRT (t * t + 1.0D0)
              s = t * c
              a (i, j) = 0.0D0
              DO 21110 k = 1, i - 1
                atemp    = c * a (k, i) - s * a (k, j)
                a (k, j) = s * a (k, i) + c * a (k, j)
                a (k, i) = atemp
21110           CONTINUE
              DO 21120 k = i + 1, j - 1
                atemp    = c * a (i, k) - s * a (k, j)
                a (k, j) = s * a (i, k) + c * a (k, j)
                a (i, k) = atemp
21120           CONTINUE
              DO 21130 k = j + 1, n
                atemp    = c * a (i, k) - s * a (j, k)
                a (j, k) = s * a (i, k) + c * a (j, k)
                a (i, k) = atemp
21130           CONTINUE
              DO 21140 k = 1, n
                vtemp    = c * v (k, i) - s * v (k, j)
                v (k, j) = s * v (k, i) + c * v (k, j)
                v (k, i) = vtemp
21140           CONTINUE
              dtemp = c * c * d (i) + s * s * d (j) -&
                      2.0D0 * c * s * b
              d (j) = s * s * d (i) + c * c * d (j) +&
                      2.0D0 * c * s * b
              d (i) = dtemp
              ENDIF
21010       CONTINUE
21000     CONTINUE
20000   CONTINUE
19999 CONTINUE
      nrot = l
!
      DO 30000 j = 1, n - 1
        k = j
        dtemp = d (k)
        DO 30100 i = j + 1, n
          IF (d (i) .LT. dtemp) THEN
            k = i
            dtemp = d (k)
            ENDIF
30100     CONTINUE
        IF (k .GT. j) THEN
          d (k) = d (j)
          d (j) = dtemp
          DO 30200 i = 1, n
            dtemp    = v (i, k)
            v (i, k) = v (i, j)
            v (i, j) = dtemp
30200       CONTINUE
          ENDIF
30000   CONTINUE
!
RETURN
END subroutine jacobi        

 
