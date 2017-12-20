!
! --------------------------------------------------------------------
!
!     This program is the basic structure for outputing solution to a cgns
!     file based on the restart solution in HDF5 that the production code
!     exports.
!
!     This uses an already existed cgns file that is based on hdf5.
!     To compile: /home/leonardo/opt/hdf5/bin/h5pfc
!     -I/home/leonardo/opt/cgns/include/ cg_output.f90
!     -L/home/leonardo/opt/cgns/lib/ -lcgns -lhdf5
!
!     Remember that in order to deal with this type of datafile, it is needed
!     that the cgns compilation has been done with proper parallel and hdf5
!     support.
!     /home/leonardo/GitHub/CGNS/src/Test_UserGuideCode/Fortran_code/read_grid_unst.F90
! --------------------------------------------------------------------
!
      program cg_output                                                 

         use cgns
         implicit none

!        Error handling.
         integer(kind=4) :: ier

!        Identification of the file.
         integer(kind=4) :: index_file

!        Number of bases in the file.
         integer(kind=4) :: index_base = 1

!        Number of zones in the hdf file.
         integer(kind=4) :: index_zone = 1

!        Zone size and name.
         character(32) :: zonename
         integer(cgsize_t) :: isize(1,3)

!        Specific datatypes.
         integer(cgsize_t) irmin,irmax

!        Specific coordinates.
         real(kind=8), allocatable, dimension(:,:) :: corrd

         character(9) :: cg_infile = "dummy.hdf"

         write(*,*) " Openning file: ", cg_infile

!        Open CGNS file.
         call cg_open_f(cg_infile,CG_MODE_READ,index_file,ier)
         if (ier .ne. CG_OK) call cg_error_exit_f
!
         write(*,*) " The index of the file is: ", index_file
!
!        Get the number of zones.
         call cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,ier)
         if (ier .ne. CG_OK) call cg_error_exit_f
!
         write(*,*) " The name of the zone is: ", zonename
         write(*,*) " The index of the zone is: ", index_zone
         write(*,*) " isize(1,1) = ", isize(1,1)
         write(*,*) " isize(1,2) = ", isize(1,2)
         write(*,*) " isize(1,3) = ", isize(1,3)
!
!        Read the X-coordinates.
         call cg_coord_read_f(index_file,index_base,index_zone,'CoordinateX',RealSingle,irmin,irmax,x,ier)
         if (ier .ne. CG_OK) call cg_error_exit_f


      end program cg_output
