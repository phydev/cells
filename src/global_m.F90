module global_m

  implicit none

  private

  type mesh_t
     real :: phi, mu, lapl_phi, lapl_mu, lapl_h
     integer :: itype
  end type mesh_t

  public :: mesh_t

  !> some mathematical constants
  REAL(8), public, parameter :: M_Pi        = 3.1415926535897932384626433832795029
  REAL(8), public, parameter :: M_E         = 2.7182818284590452353602874713526625
  REAL(8), public, parameter :: M_TWO       = 2.0
  REAL(8), public, parameter :: M_THREE     = 3.0
  REAL(8), public, parameter :: M_FOUR      = 4.0
  REAL(8), public, parameter :: M_FIVE      = 5.0
  REAL(8), public, parameter :: M_SIX       = 6.0
  REAL(8), public, parameter :: M_SEVEN     = 7.0
  REAL(8), public, parameter :: M_EIGHT     = 8.0
  REAL(8), public, parameter :: M_NINE      = 9.0
  REAL(8), public, parameter :: M_TEN       = 10.0
  REAL(8), public, parameter :: M_HALF      = 0.5

  !  begin parameters
  real, public :: cell_radius, interface_width,  dt
  integer, public :: tstep, iseed, temp, ip_part2, ip_part
  character(len=3), public :: dir_name
  character(len=6), public :: file_name, file_id
  real, public ::  depletion_weight, phi_total
  ! mesh variables
  integer, allocatable, public :: lxyz(:,:), lxyz_inv(:,:),  grid_cell_domain(:)
  integer, public :: Lsize(1:2), np_bndry
  integer, public :: np, np_tt, ip
  real, allocatable, public :: gg(:,:), density(:), r_cm(:,:), r_cm_part(:,:)
  type(mesh_t), allocatable, public :: cell(:,:), aux(:,:)

  ! new
  real, public :: self_int, sum_int, wgamma(10)
  integer, public :: jcell, icell, tcell
  integer, allocatable, public :: ncell(:), r(:)
  real, allocatable, public :: f(:), mu(:), lapl(:)
  ! misc
  integer, public :: nstep, counter = 0, ndim, output_period, extra_steps
  real, public :: hs, time_init, time_end, ctime, phi_max, hphi, fnu
  logical, public :: periodic

  integer, public :: np_sphere, sphere(1000,1:2)
  ! local variables
  integer,allocatable, public :: lxyz_part(:,:), lxyz_inv_part(:,:), gammaw(:)
  integer, public :: Lsize_part(2), np_part, np_part_tt, ntype, nleap, dri(2), drf(2), dr(2),&
  itype, np_part_bndry, ip2, ip_global

  ! coupling constants between the cell_radius
  real, allocatable, public :: gamma_l(:), eta(:,:), beta(:,:), volume(:)
  real, allocatable, public :: adhesion(:,:), hfield(:,:), hfield_lapl(:,:)
  real, public :: adh1, adh2, v1, v2, volume_target, vol_lagrangian
  integer, public :: cm_calc_counter


end module global_m
