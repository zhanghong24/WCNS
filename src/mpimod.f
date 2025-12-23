! ===============================================================
!      <fantasy> - A Block Structured Compressible NS Solver
! ===============================================================
!  Copyright (c) 2005 Fox Liu
!  All Rights Reserved
!  Filename: mpimod.f90
!  Authors: Fox Liu
!  Date: 2006-2-26
! ===============================================================
!  Modification history:
!    Date             Programmer         Description
!    ------------------------------------------------
!
! ===============================================================



      module mpi
      implicit none
#ifdef PARALLEL
      include "mpif.h"
#endif      
      end module mpi