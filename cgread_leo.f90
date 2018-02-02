! ------------------------------------------------------------------------------
!
!       This is a example program that shows how to read a cgns file and
!       exports it to BRU3D format arrays. It is designd to be
!       automatically integrated with BRU3D code. Howver, it shall be
!       further tested in the near future.
!
! ------------------------------------------------------------------------------
!
        program read_mesh
        USE CGNS
        implicit none

        ! The nodes coordinates shall be stored in an array:
        ! xyznode(i,1) = x
        ! xyznode(i,2) = y
        ! xyznode(i,3) = z

        ! The number of nodes per element shall be stored in an array
        ! called:
        ! ielenodes(i) = number of numdes for that volume.

        ! The connectivity shall be stored in an array called:
        ! itable(i,n_nodes) = index of the nodes per volume.

        ! The boundary conditions shall be stored in an array named:
        ! ifacetype(i) = type of the boundary condition.

        ! The nodes that form each one of the boundary conditions shall
        ! be stored in an extra data structure named:
        ! bc_nodes(i,n_index) = index of the boundary condition.

        end program read_mesh
!
! ------------------------------------------------------------------------------
!
!        subroutine BHDF_READ_RESTART: 
!
!       This subroutine receives a vector of properties and reads it in
!       the same way an in house code would.
!
        subroutine bhdf_read_restart
        use CGNS
        use HDF5
        implicit none

        end subroutine bhdf_read_restart
!
! ------------------------------------------------------------------------------
!
!        subroutine BHDF_WRITE_RESTART: 
!
!       This subroutine receives a vector of properties and writes it in
!       the same way an in-house code would.
!
        subroutine bhdf_write_restart
        use CGNS
        use HDF5
        implicit none

        end subroutine bhdf_write_restart
!
! ------------------------------------------------------------------------------
!
!        subroutine READ_MIXED_ELEMENTS: 
!
!        This subroutine reads an unstructured hdf cgns file and exports
!        the necessary datastructures for an in-house CFD code.
!
        subroutine read_mixed_elements
        use CGNS
        implicit none
!
        integer(kind=4) :: Cdim, Pdim, Idim, ier, count
        integer(kind=4) :: cg, base, zone, ZoneType, nelem, npe
        integer(kind=4) :: nbases, nzones, ncoords, nsections 
        integer(kind=4) :: i, n, sect, nbndry, type, parent_flag
!
        integer(cgsize_t) :: size_grid(1,3)
        integer(cgsize_t) :: range_min, range_max
        integer(cgsize_t) :: start, end
        integer(cgsize_t) :: elements(1000)
        integer(cgsize_t) :: ElementDataSize
        integer(cgsize_t) :: parent_data
!
        character(32) :: filename, nodename
!
        real(kind=8), allocatable,dimension(:) :: data_double
!
!       initialize
        ier = 0
!
!       Open cgns file.
        call cg_open_f('dummy_mix.cgns', MODE_READ, cg, ier)
        ! call cg_open_f('dummy.cgns', MODE_READ, cg, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
!
!       Dump...
        write(*,*)' -- File Opened and Read !'
!
!       Get the number of bases. Here, just one base is considered.
        call cg_nbases_f(cg, nbases, ier) 
        if (ier .eq. ERROR) call cg_error_exit_f
!
!       Check the number of bases.
        if (nbases.gt.1) write(*,*)'ERROR: This program reads only the first base'
!
!       Set base index. Here, we would have to loop through bases if
!       more than one was found.
        base = 1
        call cg_base_read_f(cg, base, nodename, Cdim, Pdim, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
!
!       Dump the informations about the base we found.
        write(*,*)'  --- CGNSBase_t node ---'
        write(*,*)'  Name = ', nodename
        write(*,*)'  CellDimension = ', Cdim
        write(*,*)'  PhysDimension = ', Pdim
        write(*,*)' '
!
!       Now, down in the base, lets read each zone. Note that here we
!       are still considering just one zone, but if we want to
!       partition.
        call cg_nzones_f(cg, base, nzones, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
!
!       Set one zone index.
        zone = 1
        write(*,*) ' -- The number of zones is: ',nzones
!
!       Read the zone.
        call cg_zone_read_f(cg,base,zone,nodename,size_grid,ier)
        if (ier .eq. ERROR) call cg_error_exit_f
        write(*,*) ' -- The name of zones is: ',nodename

        call cg_zone_type_f(cg, base, zone, ZoneType, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
        write(*,*) ' -- The type of zones is: ',nodename

        Idim=Cdim
        if (ZoneType.eq.Unstructured) Idim=1
!
        write(*,*)'  --- Zone_t node ---'
        write(*,*)'  Name     = ',nodename
        write(*,*)'  ZoneType = ',ZoneTypeName(ZoneType)
        write(*,*)'  Size     = ',size_grid(1,1)
!
!       Get the number of coordinates.
        call cg_ncoords_f(cg, base, zone, ncoords, ier)
        if (ier.eq.ERROR) call cg_error_exit_f
        write(*,*) ' -- ncoords = ',ncoords
!
        write(*,*) ' -- Idim    = ',Idim
        range_min = 1
        range_max = size_grid(1,1)
!
        allocate(data_double(range_max))
        data_double(:) = 0.0d0
!
!       Reads X.
        call cg_coord_read_f(cg, base, zone, 'CoordinateX', &
              RealDouble, range_min, range_max, data_double, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
        write(*,*) 'first point:',data_double(1)
        write(*,*) 'last point :',data_double(range_max)
!
!       Reads Y.
        call cg_coord_read_f(cg, base, zone, 'CoordinateY', &
              RealDouble, range_min, range_max, data_double, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
        write(*,*) 'first point:',data_double(1)
        write(*,*) 'last point :',data_double(range_max)
!
!       Reads Z.
        call cg_coord_read_f(cg, base, zone, 'CoordinateZ', &
              RealDouble, range_min, range_max, data_double, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
        write(*,*) 'first point:',data_double(1)
        write(*,*) 'last point :',data_double(range_max)
!
        write(*,*) "---------------------------------------------------------------"
!
!       Get eleent connectivity.
        write(*,*) ' --- Elements_t Nodes ---' 
        call cg_nsections_f(cg, base, zone, nsections, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
!
        do sect=1, nsections

            call cg_section_read_f(cg, base, zone, sect, &
                    nodename, type, start, end, nbndry, &
                    parent_flag, ier)
            if (ier .eq. ERROR) call cg_error_exit_f
            
            call cg_ElementDataSize_f(cg, base, zone, sect, &
                                      ElementDataSize, ier)
            if (ier .eq. ERROR) call cg_error_exit_f
            
            call cg_elements_read_f(cg, base, zone, sect, &
                       elements, parent_data, ier)
            if (ier .eq. ERROR) call cg_error_exit_f
!
            write(*,*) ' --- SECTION: ',sect
            write(*,*) 'Name= ',nodename
            write(*,*) 'Range= ',start,end
            if (nbndry .ne. 0) write(*,*) ' Sorted elements'
!
            call cg_npe_f(type, npe, ier)
            if (ier .eq. ERROR) call cg_error_exit_f
            write(*,*) "npe = ",npe
!
            write(*,*) 'ElementDataSize =', ElementDataSize
            write(*,*) 'Section Element Type= ', ElementTypeName(type)
!
            nelem = end-start+1
            write(*,*) 'nelem  = ', nelem 
            write(*,*) 'Element Connectivity:'
!
            if (type .lt. MIXED) then
                do i=1, nelem
                   write(*,*) (elements((i-1)*npe+n),n=1,npe)
                enddo
            elseif (type .eq. MIXED) then
!
                count = 0
!
                do i=1, nelem
                    count = count + 1
                    type = elements(count)
                    if (type .gt. NGON_n) then
                        npe = type - NGON_n
                        write(*,*) 'Element Type= NGON_n, npe =',npe
                    else
                        call cg_npe_f(type, npe, ier)
                        write(*,*) 'Element Type= ', ElementTypeName(type)
                    endif
!
                    write(*,*) (elements(count+n),n=1,npe)
                    count = count+npe
                enddo
            else
!
                count = 0
!
                do i=1, nelem
                    count = count + 1
                    npe = elements(count)
                    write(*,*) 'Element Number Points=',npe
                    write(*,*) (elements(count+n),n=1,npe)
                    count = count+npe
                enddo
             endif
         enddo
!
         call cg_close_f(cg, ier)
         if (ier.eq.ERROR) call cg_error_exit_f()
         write(*,*) 'CGNS File written & closed'
!
         end subroutine read_mixed_elements
