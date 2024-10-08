!!
!!  Copyright (C) 2019  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module sponge
!*******************************************************************************
! This module contains all of the subroutines associated with scalar transport
use types, only : rprec
implicit none

private
public :: sponge_init, sponge_force

! Sponge layer properties
logical, public :: use_sponge = .false.
real(rprec), public :: sponge_frequency = 3.9_rprec
real(rprec), public :: sponge_height = 0.75_rprec
real (rprec), dimension (:), allocatable, public :: sp
integer, public :: sp_loc ! GN
real(rprec), public :: z_sp_loc = 0.0_rprec ! GN
contains

!******************************************************************************
subroutine sponge_init()
!******************************************************************************
use param, only : nz, lbz, pi, L_z, coord
use types, only : rprec
use grid_m, only : grid
use functions, only : binary_search!GN
#ifdef PPMPI
use mpi ! GN
#endif
use param, only : ierr,comm,mpi_rprec,rank_of_coord ! GN
integer :: k
real(rprec) :: z_sp_loc_dummy = 0._rprec!GN
allocate (sp(lbz:nz)); sp = 0._rprec

!find position of sponge height !GN
if (sponge_height >= grid%z(1) .and. sponge_height <= grid%z(nz)) then
         sp_loc = binary_search(grid%z(1:), sponge_height)
         z_sp_loc = grid%z(sp_loc)
         !print *, "coord_sp_loc, sp_loc, sp_h, z_sp_loc",coord, sp_loc,  sponge_height, z_sp_loc
end if
#ifdef PPMPI
   call mpi_allreduce(z_sp_loc, z_sp_loc_dummy, 1, mpi_rprec, MPI_SUM, comm, ierr)
   z_sp_loc=z_sp_loc_dummy
#endif

do k = lbz, nz
    if (grid%z(k) > sponge_height) then
        sp(k) = 0.5_rprec*sponge_frequency*(1._rprec                           &
            - cos(pi*(grid%z(k) -sponge_height)/(L_z - sponge_height)))
    end if
end do

end subroutine sponge_init

!*******************************************************************************
subroutine sponge_force
!*******************************************************************************
! This subroutine calculates the sponge force term
use param, only : nx, ny, nz, jt_total, coord, nproc
use sim_param, only :  RHSx, RHSy, RHSz, u, v, w

!#ifdef PPSCALARS
!use scalars, only : theta, RHS_T
!#endif

integer :: k

if (use_sponge) then
    do k = 1, nz-1
        RHSx(1:nx,1:ny,k) = RHSx(1:nx,1:ny,k) - 0.5_rprec*(sp(k) + sp(k+1))    &
            * (u(1:nx,1:ny,k) - sum(u(1:nx,1:ny,k))/(nx*ny))
        RHSy(1:nx,1:ny,k) = RHSy(1:nx,1:ny,k) - 0.5_rprec*(sp(k) + sp(k+1))    &
            * (v(1:nx,1:ny,k) - sum(v(1:nx,1:ny,k))/(nx*ny))
        RHSz(1:nx,1:ny,k) = RHSz(1:nx,1:ny,k) - sp(k)                          &
            * (w(1:nx,1:ny,k) - sum(w(1:nx,1:ny,k))/(nx*ny))
!#ifdef PPSCALARS
!        RHS_T(1:nx,1:ny,k) = RHS_T(1:nx,1:ny,k) - 0.5_rprec*(sp(k) + sp(k+1))    &
!            * (theta(1:nx,1:ny,k) - sum(theta(1:nx,1:ny,k))/(nx*ny))
!#endif
    end do
end if

end subroutine sponge_force


end module sponge
