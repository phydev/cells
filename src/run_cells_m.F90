!! Copyright (C) 2016 M. Moreira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!!

module run_cells_m

  use global_m
  use derivatives_m
  use sim_init_m
  use misc_m

  implicit none

  private

  public :: run_cells


  contains

    subroutine run_cells()

      implicit none


      ! initializing parameters
      call  parameters_init(cell_radius, ntype, density, interface_width, tstep, dt, Lsize, dr, dir_name, iseed,&
           np_bndry, depletion_weight, output_period, periodic)

      ! number of points in the mesh
      np = 4*Lsize(1)*Lsize(2) ! number of points
      np_tt = 4*(Lsize(1)+2*np_bndry)*(Lsize(2)+2*np_bndry) ! number of points plus boundary points

      ALLOCATE(ncell(1:ntype))

      ncell(1:ntype) = (np*density(1:ntype))/(4*cell_radius**2)
      tcell = sum(ncell(1:ntype))

      ! partition mesh
      Lsize_part(1:2) = (/20, 20/)
      np_part = 4*Lsize_part(1)*Lsize_part(2)
      np_part_bndry = 2*Lsize(1)
      np_part_tt = 4*(Lsize_part(1) + 2*np_part_bndry)*(Lsize_part(2) + 2*np_part_bndry)
      ! allocating matrices and vectors

      ALLOCATE(lxyz(np_tt,1:2))
      ALLOCATE(lxyz_inv(-Lsize(1)-np_bndry:Lsize(1)+np_bndry, &
                        -Lsize(2)-np_bndry:Lsize(2)+np_bndry))
      ALLOCATE(lxyz_part(np_part_tt,1:2))
      ALLOCATE(lxyz_inv_part(-Lsize_part(1)-np_part_bndry:Lsize_part(1)+np_part_bndry, &
                        -Lsize_part(2)-np_part_bndry:Lsize_part(2)+np_part_bndry))

      ALLOCATE(cell(0:np_part,tcell))
      ALLOCATE(aux(np,ntype))

      ALLOCATE(gg(1:np,1:2))

      ALLOCATE(r(tcell))

      ALLOCATE(r_cm(1:tcell,1:2))

      call system('rm 001/*')

!      ALLOCATE(gammaw(ntype))
!      gammaw(1:2) = (/ 1.0, 2.0 /)
      cell(0,:)%phi = -1.d0
      cell(0,:)%mu = -2.d0
      cell(0,:)%lapl_mu = 0.d0
      cell(0,:)%lapl_phi = 0.d0
      ! initializing space matrices
      call space_init(Lsize, lxyz, lxyz_inv, np_bndry, np, .true.)  ! auxiliar fields
      call space_init(Lsize_part, lxyz_part, lxyz_inv_part, np_part_bndry, np_part, .false.) ! single cell fields periodic boundary

      ! generating sphere points
      call gen_cell_points(cell_radius,sphere,np_sphere)

      call single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first=.true.)

      call cahn_hilliard(cell(0:np_part,1), 20000, np_part, 1.0, lxyz_part, lxyz_inv_part)

      if(tcell.gt.1) then
         call single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first=.false.)
      end if


      ! initializing simulation box

      dr = (/ int(2*cell_radius), ntype*int(2*cell_radius) /)
      dri = (/ int(cell_radius), int(cell_radius) /)
      icell = 1
      do itype=1, ntype
         drf = (/ cell_radius, cell_radius + (ntype-itype)*(int(2*cell_radius)) /)
         dr = (/ int(2*cell_radius), ntype*int(2*cell_radius) /)
         nleap = tcell-ncell(itype)
         temp = ncell(itype) +(itype-1)*ncell(itype-1)
         call cell_pos_init(r, icell, nleap,  dr, dri, drf, temp, cell_radius, &
                                          lxyz, lxyz_inv, Lsize, np, sphere, np_sphere, iseed)
         dri(2) = dri(2) +   2*int(cell_radius)

      end do

      call hfield_calc(cell, aux, r, Lsize, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np )


      ! printing header
      call print_header(Lsize, tcell, ntype, ncell, dir_name, periodic)

      write(*,'(A)') "Initiating the core program..."

      ! saving the initial condition
      call output_aux(cell, 100, 'phi00',dir_name, 1, 1, np_part, lxyz_part)
      call output_aux(aux, 200, 'aux00',dir_name, 1, ntype, np, lxyz)
      cell(0,:)%phi = -1.d0
      cell(0,:)%mu = 0.d0

      nstep = 0

      do while(nstep<=tstep)
         nstep = nstep + 1

         call CPU_TIME(time_init)

         ! calculating gradient of vegf
         !call dderivatives_grad(cell, gg, np, lxyz, lxyz_inv, dr)

         call hfield_calc(cell, aux, r, Lsize, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np )

         ! calculating laplacian of phi

         do icell = 1, tcell
            call dderivatives_lapl(cell(0:np_part,icell)%phi, cell(1:np_part,icell)%lapl_phi ,&
                 np_part, dr, lxyz_part, lxyz_inv_part)
         end do

         ! Chemical potential:  phi**3 - phi - epsilon*laplacian(phi )

         do ip = 1, np_part
            do icell=1, tcell

              ip2 = lxyz_inv( lxyz_part(ip,1) + lxyz(r(icell),1), lxyz_part(ip,2) + lxyz(r(icell),2))
              fnu = 0.d0
              do itype=1, ntype
                 fnu = fnu + aux(ip2,itype)%phi - hfunc( cell(ip,icell)%phi )*deltak(cell(ip,icell)%itype,itype)
              end do

              cell(ip,icell)%mu = interface_width*cell(ip,icell)%lapl_phi +&
                   0.25*(cell(ip,icell)%phi+1)*(1.d0-cell(ip,icell)%phi)*( cell(ip,icell)%phi - 2.0*depletion_weight*fnu )

            end do
         end do


         ! Calculating laplacian of mu
         do icell = 1, tcell
            call dderivatives_lapl(cell(0:np_part,icell)%mu, cell(1:np_part,icell)%lapl_mu, np_part, dr, lxyz_part, lxyz_inv_part)
         end do


         cell(1:np_part, 1:tcell)%phi = cell(1:np_part,1:tcell)%phi - dt*(cell(1:np_part,1:tcell)%lapl_mu)

         call cm_calc(r_cm, cell, tcell, np, np_part, r, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part)




         call CPU_TIME(time_end)

         ctime = ctime + (time_end - time_init)

         if(nstep.eq.100) then
            ctime = (ctime*(tstep-nstep) )/6000.d0

            if( ctime>60.d0) then
               write(*,'(A,F10.2)') "Estimated time (hour): ",ctime/60.d0
            else
               write(*,'(A,F10.2)') "Estimated time (min): ",ctime
            end if
         end if



         ! output

         counter = counter + 1


         if(counter.eq.output_period) then
            write(*,*) nstep
            counter = 0

            write(file_name,'(I5)') nstep
            call output_aux(cell, nstep+1, trim(file_name),dir_name, 1, 1, np_part, lxyz_part)

            OPEN (UNIT=nstep,FILE=dir_name//'/phi'//trim(file_name)//'.xyz')

            do ip=1, np

              do itype = 1, ntype
                 if(aux(ip,itype)%phi>=0.9) then
                    write(nstep,'(I10,I10,F10.2,I10)') lxyz(ip,1:2),aux(ip,itype)%phi, itype
                 end if
              end do
            end do
            close(nstep)

         end if
         ! end of the output


      end do

      DEALLOCATE(ncell)
      DEALLOCATE(lxyz)
      DEALLOCATE(lxyz_inv)
      DEALLOCATE(lxyz_part)
      DEALLOCATE(lxyz_inv_part)
      DEALLOCATE(cell)
      DEALLOCATE(aux)
      DEALLOCATE(gg)
      DEALLOCATE(r)

    end subroutine run_cells



    subroutine cm_calc(r_cm, f, tcell, np, np_part, r, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part)

      implicit none

      integer, intent(in) :: np, tcell, np_part
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_part(:,:), lxyz_inv(:,:), lxyz_inv_part(:,:)
      type(mesh_t), allocatable, intent(inout) :: f(:,:)
      real, allocatable, intent(inout) :: r_cm(:,:)
      integer, allocatable, intent(inout) :: r(:)

      integer :: ip, icell, ncell, ip_part, ip_global_m
      real :: volume(1:tcell), delta_r(1:tcell,1:2), ftemp(1:np_part)

      volume = 0.d0
      r_cm(tcell,1:2) = 0.d0
      do ip=1, np

         do icell=1, tcell
            ip_global_m = r(icell)

            call vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

            if(f(ip_part,icell)%phi.gt.0.d0) then
               volume(icell) = volume(icell) + f(ip_part,icell)%phi
               r_cm(icell,1:2) = r_cm(icell,1:2) + f(ip_part,icell)%phi*lxyz(ip,1:2)
            end if
         end do
      end do

      r_cm(1:tcell,1) = r_cm(1:tcell,1)/volume(1:tcell)
      r_cm(1:tcell,2) = r_cm(1:tcell,2)/volume(1:tcell)

      do icell=1, tcell
        delta_r(icell,1:2) = int(r_cm(icell,1:2)-lxyz(r(icell),1:2))
        !write(*,*) "del__", delta_r(icell,1:2), icell
        !write(*,*) "pos_i", lxyz(r(icell),1:2)
        !write(*,*) "pos_f",r_cm(icell,1:2)
        !ftemp(:) = -1.d0
        do ip_part=1, np_part
          !if(f(lxyz_inv_part( lxyz_part(ip_part,1) - delta_r(icell,1), delta_r(icell,2) - delta_r(icell,2) ), icell )%phi.ge.0) then
            ftemp(ip_part) = f(lxyz_inv_part( lxyz_part(ip_part,1) - delta_r(icell,1), lxyz_part(ip_part,2) - delta_r(icell,2) ), icell )%phi

          !end if
        end do

        f(1:np_part,icell)%phi = ftemp(1:np_part)
        r(icell) = lxyz_inv(int(anint(r_cm(icell,1))) , anint(anint(r_cm(icell,2))))
      end do

      ! Volume = sum phi_i
      ! (sum phi_i r_i )/Volume
    end subroutine cm_calc


    subroutine vec_local2global(ip, ip_global_m, ip_local, lxyz, lxyz_inv, lxyz_part)

      implicit none

      integer, allocatable, intent(in) :: lxyz(:,:),  lxyz_inv(:,:), lxyz_part(:,:)
      integer, intent(in) :: ip_global_m, ip_local
      integer, intent(out) :: ip

      ip =  lxyz_inv( lxyz_part(ip_local,1) +lxyz(ip_global_m,1), lxyz_part(ip_local,2) + lxyz(ip_global_m,2) )

    end subroutine vec_local2global


    subroutine vec_global2local(ip_part, ip_global_m, ip, lxyz, lxyz_inv, lxyz_inv_part)

      implicit none

      integer, allocatable, intent(in) :: lxyz(:,:),  lxyz_inv(:,:),  lxyz_inv_part(:,:)
      integer, intent(in) :: ip_global_m, ip
      integer, intent(out) :: ip_part

      ip_part =  lxyz_inv_part( lxyz(ip,1) - lxyz(ip_global_m,1), lxyz(ip,2) - lxyz(ip_global_m,2))

    end subroutine  vec_global2local

end module run_cells_m
