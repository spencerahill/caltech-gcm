module grid_physics

  ! GRID_PHYSICS provides routines that can be combined to form a "physics"
  ! package in an idealized GCM.
  !
  ! The namelist GRID_PHYS.NML contains parameters for the parameterizations.
  !
  ! Radiative transfer is parameterized as Newtonian relaxation toward the
  ! radiative-equilibrium state of a semi-gray atmosphere with surface temperature
  ! specified in SURFACE_TEMPERATURE_FORCED.
  !
  ! The surface temperature of the radiative-equilibrium state depends on the
  ! parameters DELH and TSFC_SP, the pole-to-equator temperature contrast and the
  ! surface temperature at the south pole, respectively. The extent of the frictional
  ! boundary layer is determined by SIGMA_B. For the Newtonian relaxation
  ! parameterizations, SIGMA_B also determines the height of the low-latitudinal layer
  ! in which the relaxation coefficient is increased compared with the interior
  ! atmospheric relaxation coefficient.

  use constants_mod, only: KAPPA, CP_AIR, GRAV, RDGAS, RVGAS, STEFAN, HLV, &
       TFREEZE, PI, DEG_TO_RAD, RAD_TO_DEG, OMEGA, RADIUS, &
       DAY => SECONDS_PER_DAY

  use fms_mod, only: error_mesg, fatal, file_exist, open_namelist_file, &
       check_nml_error, mpp_pe, mpp_root_pe, close_file, &
       write_version_number, stdlog

  use time_manager_mod, only: time_type
  use diag_manager_mod, only: register_diag_field, send_data

  implicit none

  ! constants and defaults for parameters
  real, private :: &
       reference_sea_level_press = 1.e5, &
       delh = 120., &  ! pole-equator surface temperature contrast
       tsfc_sp = 320., &  ! surface temperature at south pole
       phi0 = 0., &   ! heating offset (degrees)
       lat_tropic = 32.092, &  ! latitude within/outside surface temps are flat
       t_strat = 200., &  ! "skin" temperature of atmosphere
       scale_height_ratio = 3.5, &  ! ratio of water vapor and pressure scale heights
       ka_days = 50., &   ! latitude-independent part of Newtonian
       ka, &  ! relaxation coefficient
       ks_days = 4., &  ! Newtonian damping coefficient in low
       ks, &  ! latitudes near surface
       Cdrag = 5.e-6, &  ! quadratic "drag coefficient" [1/m]

       ! Maximum surface temperature for Emanuel 95 critical profile
       e95_max_temp = 350., &
       !  Assumed tropopause-surface temperature difference
       e95_trop_sfc_temp_diff = 100., &
       ! Factor to make Emanuel 95 profiles sub- or super- critical.
       e95_crit_multiplier = 1., &
       ! Latitude (in degrees) at which to switch from the Emanuel 95 profile
       ! to its tangent, from here to the same latitude in the opposite
       ! hemisphere.
       e95_lat_switch_tangent = 10., &
       ! Extend the E95 profile this far in degrees latitude past the forcing
       ! maximum before switching to flat.
       e95_lat_past_max = 0.

  real, public :: &
       sigma_b = 0.85  ! extent of `PBL' (also used to compute relaxation rates)

  logical, private :: &
       flatxt = .false., &  ! Flatten temperatures in extratropics.
       flatt  = .false., &  ! Flatten temperatures in tropics.
       ! Use critical temperature profile from Emanuel 1995.
       do_e95_crit_temps = .false., &
       ! Flatten temperatures in winter hemisphere (which must be the southern
       ! hemisphere) and poleward of maximum to break equatorial symmetry.
       e95_break_eq_symm = .false., &
       ! Use dry adiabat as the vertical structure of the equilibrium
       ! temperature profile
       do_equil_adiabat = .false.

  namelist /grid_phys_list/    &
       reference_sea_level_press, delh, tsfc_sp, phi0, lat_tropic, t_strat, &
       scale_height_ratio, ka_days, ks_days, Cdrag, sigma_b, flatxt, flatt, &
       do_e95_crit_temps, e95_break_eq_symm, e95_max_temp, &
       e95_trop_sfc_temp_diff, e95_crit_multiplier, e95_lat_switch_tangent, &
       e95_lat_past_max, do_equil_adiabat

  private grid_phys_list

  character(len=128) :: version='$Id: grid_physics_simple.f90 $'
  character(len=128) :: tag='homemade'
  character(len=14)  :: mod_name = "grid_physics"

  ! diagnostic IDs
   integer id_tdt, id_teq, id_drag_x, id_drag_y
  real :: missing_value = -1.e10
contains

  !-------------------------------------------------------------------------------

  subroutine grid_phys_init(axes, Time)
    ! reads namelist file with parameters and echoes parameter values

    integer, intent(in) :: axes(4)
    type(time_type), intent(in) :: Time

    integer  unit, io, ierr

    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
          read  (unit, nml=grid_phys_list, iostat=io, end=100)
          ierr = check_nml_error (io, 'grid_phys_list')
       enddo
100     call close_file (unit)
    endif

!     ----- write version info and namelist to log file -----

    call write_version_number (version,tag)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=grid_phys_list)

!     ---------------write info to screen--------------------


    if(mpp_pe()==mpp_root_pe()) then
      if (do_e95_crit_temps) then
        write(*, *)
        write(*, *) 'Using critical temperature profile from Emanuel 1995'
        write(*,10) 'Offset of temperature maximum from equator:', phi0, 'deg'
        write(*,10) 'Emanuel 95 max surface temperature:', &
             e95_max_temp, 'K'
        write(*,10) 'Emanuel 95 tropopause-surface temperature difference:', &
             e95_trop_sfc_temp_diff, 'K'
        write(*,10) 'Emanuel 95 criticality multiplier:', &
             e95_crit_multiplier, '(dimensionless)'
        write(*,*) 'Emanuel 95 equatorial symmetry broken:', &
             e95_break_eq_symm
      else
        write(*, *) 'Parameters in radiative equilibrium surface temperature:'
        write(*,10) 'South pole surface temperature:             ', tsfc_sp, &
             'K'
        write(*,10) 'Equator-pole surface temperature contrast:  ', delh, 'K'
        write(*,10) 'Offset of temperature maximum from equator: ', phi0, 'deg'

        write(*, *) 'Relaxation toward radiative equilibrium state of a'
        write(*, *) 'semi-gray atmosphere with an isothermal stratosphere'
        write(*,10) 'at temperature t_strat =', t_strat, 'K.'
      endif

      write(*, *)
      write(*,10) 'Relaxation time scale in interior atmosphere:       ', &
           ka_days, 'days'
      write(*,10) 'Relaxation time scale in low latitudes near surface:', &
           ks_days, 'days'

      write(*, *)
      if (do_equil_adiabat) then
         write(*, *) 'Using dry adiabat for equilibrium temperature profile.'
      else
         write(*,10) 'Ratio of absorber scale height to pressure scale height:',&
              scale_height_ratio
      endif

      write(*, *)
      write(*, *) 'Quadratic friction in boundary layer.'
      write(*,10) 'Extent of frictional layer: sigma_b =', sigma_b
      write(*,20) 'Drag coeff. for quadratic PBL drag: Cdrag =', Cdrag, &
           'm**(-1)'
      write(*, *)
    endif

10  format(1x,a,1x,f8.2,1x,a)
20  format(1x,a,1x,e8.2,1x,a)

    ka = 1./(ka_days * day)
    ks = 1./(ks_days * day)


    ! register fields with diagnostic manager
    id_teq = register_diag_field( mod_name, 'teq',                         &
         axes(1:3), Time, 'equilibrium temperature', 'deg_K',              &
         missing_value=missing_value)

    id_tdt = register_diag_field( mod_name, 'tdt_total_grid_physics',      &
         axes(1:3), Time, 'total temperature tendency', 'deg_K/sec',       &
         missing_value=missing_value)

    id_drag_x = register_diag_field( mod_name, 'u_quad_drag_tend',         &
         axes(1:3), Time, 'zonal wind tendency from quadratic drag',       &
         'm/sec/sec', missing_value=missing_value)

    id_drag_y = register_diag_field( mod_name, 'v_quad_drag_tend',            &
         axes(1:3), Time, 'meridional wind tendency from quadratic drag',  &
         'm/sec/sec', missing_value=missing_value)

  end subroutine grid_phys_init

  !------------------------------------------------------------------

  subroutine compute_grid_physics( is,       ie,       js,       je,      &
                   Time, delta_t, lat,    dt_ug,    dt_vg,    dt_tg,      &
                      dt_qg,   dt_trg,       ug,       vg,       tg,      &
                         qg,      trg,   p_half,   p_full,                &
                ug_previous,        vg_previous,        tg_previous,      &
                qg_previous,       trg_previous,    p_half_previous,      &
            p_full_previous)

    integer, intent(in)                     :: is, ie, js, je
    type(time_type), intent(in)            ::  Time
    real, intent(in)                       ::  delta_t
    real, intent(in), dimension (:,:,:)     :: ug, vg, tg, qg, p_half, p_full
    real, intent(in), dimension (:,:,:,:)   :: trg
    real, intent(in), dimension (:,:)       :: lat

    real, optional, dimension (:,:,:)       :: ug_previous, vg_previous, &
         tg_previous, qg_previous, p_half_previous, p_full_previous
    real, optional, dimension (:,:,:,:)     :: trg_previous

    real, intent(inout), dimension(:,:,:)   :: dt_ug, dt_vg, dt_tg, dt_qg
    real, intent(inout), dimension(:,:,:,:) :: dt_trg


    if(present(ug_previous)) then
       call surface_drag(Time, is, js, ug_previous, vg_previous, &
            p_half_previous, p_full_previous, dt_ug, dt_vg)
       call diabatic_forcing(Time, is, js, tg_previous, p_half_previous, &
            p_full_previous, lat, dt_tg)
    else
       call error_mesg('compute_grid_physics',                             &
            'ug_previous is not present',fatal)
    end if

  end subroutine compute_grid_physics

  !---------------------------------------------------------------------

  subroutine compute_grid_adjustment(tg, qg, p_half, p_full)

    real, intent(in), dimension (:,:,:) :: tg, qg, p_full
    real, intent(in), dimension (:,:,:) :: p_half

    ! dummy routine

  end subroutine compute_grid_adjustment

  !-----------------------------------------------------------------

  subroutine surface_drag(Time, is, js, ug, vg, p_half, p_full, dt_ug, dt_vg)

    integer, intent(in)                   :: is, js
    type(time_type), intent(in)           :: Time
    real, intent(in), dimension (:,:,:)    :: ug, vg, p_full, p_half
    real, intent(inout), dimension(:,:,:)  :: dt_ug, dt_vg
    real, dimension(size(ug,1),size(ug,2)) :: sigma, sigma_norm, sigma_max
    real, dimension(size(ug, 1), size(ug, 2), size(ug, 3)) :: &
         vg_norm, drag_x, drag_y
    integer :: k, num_level
    logical :: used

    num_level = size(ug,3)    ! pseudo-parameter for number of levels

    vg_norm = sqrt(ug**2 + vg**2)
    do k = 1, num_level
       sigma(:, :)  = p_full(:,:,k) / p_half(:,:,num_level+1)
       sigma_norm   = (sigma - sigma_b) / (1.0 - sigma_b)
       sigma_max    = max(sigma_norm, 0.0)
       drag_x(:,:,k) = -1.0 * Cdrag * sigma_max * vg_norm(:,:,k) * ug(:,:,k)
       drag_y(:,:,k) = -1.0 * Cdrag * sigma_max * vg_norm(:,:,k) * vg(:,:,k)
       dt_ug(:,:,k) = dt_ug(:,:,k) + drag_x(:,:,k)
       dt_vg(:,:,k) = dt_vg(:,:,k) + drag_y(:,:,k)
    end do

    ! send data to diagnostic manager
    if (id_drag_x > 0) used = send_data(id_drag_x, drag_x, Time, is, js)
    if (id_drag_y > 0) used = send_data(id_drag_y, drag_y, Time, is, js)

  end subroutine surface_drag

  !-----------------------------------------------------------------

  subroutine diabatic_forcing(Time, is, js, tg, p_half, p_full, lat, dt_tg)

    implicit none

    integer, intent(in)                   :: is, js
    type(time_type), intent(in)           :: Time
    ! temperature and pressure fields
    real, intent(in), dimension(:, :, :)  :: tg, p_full
    real, intent(in), dimension(:, :, :)  :: p_half
    real, intent(in), dimension(:, :)     :: lat

    ! time tendency in temperature field
    real, intent(inout), dimension(:,:,:) :: dt_tg

    ! local variables
    real, dimension(size(tg, 1), size(tg, 2), size(tg, 3)) ::              &
         heating_rate

    integer :: num_level
    logical :: used


    ! pseudo-parameters
    num_level = size(tg, 3)

    if (size(lat,1) .ne. size(tg,1) &
         .or. size(lat,2) .ne. size(tg,2)) then
       call error_mesg('diabatic_forcing',                             &
            'mismatched argument list',fatal)
       ! terminate execution
    end if

    heating_rate = relax_to_gray_equilibrium(Time, tg, lat,      &
         p_full, p_half, is, js)

    dt_tg        = dt_tg + heating_rate

! send data to diagnostic manager
       if (id_tdt > 0) used = send_data(id_tdt, heating_rate, Time,      &
            is, js)

  end subroutine diabatic_forcing

  !-----------------------------------------------------------------------

  function relax_to_gray_equilibrium(Time, temp, lat, p_full, p_half, is, js)&
       result (heating_rate)

    ! computes heating rate from relaxation towards radiative
    ! equilibrium state of a gray atmosphere with an isothermal
    ! stratosphere

    implicit none

    type(time_type), intent(in)           :: Time

    integer, intent(in) :: is, js

    real, intent(in), dimension(:, :) ::                                   &
         lat
    real, intent(in)                  ::                                   &
         temp(:, :, :),                                                    &
         p_full(:, :, :),                                                  &
         p_half(:, :, :)

    ! local variables
    real, dimension(size(temp, 1), size(temp, 2)) ::                       &
         sin_lat_2,                                                        &
         cos_lat_2,                                                        &
         cos_lat_4,                                                        &
         cos_lat_8,                                                        &
         sigma,                                                            &
         sigma_norm,                                                       &
         sigma_max,                                                        &
         kt,                                                               &
         t_sfc,                                                            &
         t_strat_var,                                                      &
         optical_thickness, &
         t_adiabat

    real, dimension(size(temp, 1), size(temp, 2), size(temp, 3)) ::        &
         p_full_ref,                                                       &
         temp_eq,                                                          &
         heating_rate

    integer :: num_level, k
    logical used

    num_level = size(temp, 3)
    cos_lat_2 = cos(lat)*cos(lat)
    cos_lat_4 = cos_lat_2**2
    cos_lat_8 = cos_lat_2**4

    if (size(lat,1) .ne. size(p_full,1) &
         .or. size(lat,2) .ne. size(p_full,2)) then
       call error_mesg('relax_to_gray_equilibrium',                        &
            'mismatched argument list',fatal)
       ! terminate execution
    end if

    ! normalized full level pressure and functions thereof
    p_full_ref = p_full / reference_sea_level_press

    ! surface temperature in background state
    if (do_e95_crit_temps) then
       t_sfc = sfc_temp_emanuel95_crit(lat)
    else
       t_sfc = surface_temperature_forced(lat)
    end if

    ! latitude-dependent optical thickness
    optical_thickness = (t_sfc / t_strat)**4 - 1.

    do k = num_level, 1, -1
       ! compute current sigma-level to get height-dependence of cooling coeff
       sigma      = p_full(:, :, k) / p_half(:, :, num_level+1)
       sigma_norm = (sigma - sigma_b) / (1.0 - sigma_b)
       sigma_max  = max(sigma_norm, 0.0) ! sigma_max = 0 outside `PBL'

       ! Newtonian damping coefficient
       ! kt = ka + (ks - ka) * sigma_max * cos(lat - phi0*DEG_TO_RAD )**8
       ! Spencer Hill: removing latitudinal dependence.  Simulating somewhat
       ! exotic atmospheres for which the ascent location isn't well known a
       ! priori, even given the imposed temperature profile being relaxed to.
       ! So making the relaxation uniform in latitude seems the wisest choice.
       ! Also, the above formulations yields a damping coefficient that is
       ! non-monotonic latitude in the winter hemisphere for off-equatorial
       ! phi0.
       kt = ka + (ks - ka) * sigma_max


       ! Spencer Hill: relax towards an adiabat, rather than a radiative
       ! equilibrium profile
       if (do_equil_adiabat) then
          t_adiabat = t_sfc * &
               (p_full(:,:,k) / reference_sea_level_press)**KAPPA
          temp_eq(:,:,k) = max(t_adiabat, t_strat)
       else
         ! temperature in background state is radiative equilibrium
         ! of a gray atmosphere (cf. Goody and Yung, pp. 392 ff.)

         ! Spencer Hill: replacing p_full_ref array, which varies in space and
         ! time, with sigma, which does not.  For my simulations, conceptual
         ! simplicity is more important than physical realism, so it's better to
         ! use a fixed vertical structure of the relaxation field than one that
         ! more accurately approximates a true radiative-equilibrium temperature
         ! profile at each location, given it's surface pressure at that
         ! timestep.
         temp_eq(:,:,k) =                                                  &
              t_strat * (1. + optical_thickness *                          &
              sigma**scale_height_ratio)**(1./4.)
      endif

      ! relax towards background state
      heating_rate(:, :, k) = ( temp_eq(:,:,k)  - temp(:, :, k) ) * kt

    end do

    if (id_teq > 0) used = send_data ( id_teq, temp_eq,  Time, is, js)

  end function relax_to_gray_equilibrium

  !--------------------------------------------------------------------------

  function surface_temperature(temp, p_full, p_half)

    ! estimate the surface temperature from the temperature on lowest
    ! full level by assuming that the potential temperature on the surface
    ! is equal to the potential temperature on the lowest full level

    implicit none

    real, intent(in) ::                                                    &
         temp(:, :, :),    & ! temperature field
         p_full(:, :, :),  & ! full-level pressure
         p_half(:, :, :)     ! half-level pressure

    real, dimension(size(temp, 1), size(temp, 2)) ::                       &
         surface_temperature

    integer :: num_level
    num_level = size(temp, 3)

    surface_temperature = temp(:, :, num_level)                            &
         * ( p_half(:, :, num_level+1) / p_full(:, :, num_level) )**kappa

  end function surface_temperature

  !-----------------------------------------------------------------------

  function surface_temperature_forced(lat)

    ! returns a prescribed surface temperature

    implicit none

    real, intent(in) ::                                                    &
         lat(:, :)

    real, dimension(size(lat, 1), size(lat, 2)) ::                         &
         surface_temperature_forced,                                       &
         latitude_dependence

    integer :: i

    real match

    latitude_dependence = cos(lat)**2 + 2. * sin(phi0*DEG_TO_RAD) * &
         ( 1 + sin(lat) )
    surface_temperature_forced = tsfc_sp + delh * latitude_dependence

    ! make surface temperature flat inside/outside lat_tropic
    if(flatt .or. flatxt) then
       match = tsfc_sp + delh * cos(lat_tropic*DEG_TO_RAD)**2
       do i=1,size(lat,2)
          if(flatt .and. (abs(lat(1,i)*RAD_TO_DEG) < lat_tropic) ) &
               surface_temperature_forced(:, i) = match
          if(flatxt .and. (abs(lat(1,i)*RAD_TO_DEG) > lat_tropic) ) &
               surface_temperature_forced(:, i) = match
       enddo
    endif

  end function surface_temperature_forced

  !-----------------------------------------------------------------------

  function sfc_temp_emanuel95_crit(lat) &
    result(sfc_temp)

    ! Spencer Hill, 2017-09.  Replaces the standard relaxation temperature
    ! profile defined in `relax_to_gray_equilibrium` with a "critical" profile
    ! as defined in Eq. (11) of Emanuel (1995).  That expression is for
    ! subcloud equivalent potential temperature.  For a dry atmosphere, that
    ! effectively amounts to the boundary layer potential temperature.  We make
    ! the additional assumption that this approximately equals the actual
    ! surface temperature.

    implicit none

    real, intent(in), dimension(:, :) :: lat

    real, dimension(size(lat, 1), size(lat, 2)) :: &
         lat_deg, sfc_temp, cos_lat, cos_lat_2, lat_dependence, &
         sfc_temp_deriv, tangent_at_switch

    real :: chi, cos_phi0
    real :: &
         cos_lat_switch, tangent_multiplier, deriv_at_switch, &
         sfc_temp_at_switch, tangent_at_eq
    real :: lat_past_max, cos_lat_past_max, sfc_temp_past_max

    integer :: j

    lat_deg = lat*RAD_TO_DEG
    if (int(phi0) .eq. 0) then
       cos_phi0 = 1.
    else
       cos_phi0 = cos(phi0*DEG_TO_RAD)
    endif
    cos_lat = cos(lat)
    cos_lat_2 = cos_lat * cos_lat

    chi = e95_crit_multiplier * (OMEGA**2)*(RADIUS**2) / &
         (CP_AIR*e95_trop_sfc_temp_diff)
    lat_dependence = (cos_phi0**2 - cos_lat_2)**2 / cos_lat_2

    sfc_temp = e95_max_temp * exp(-0.5 * chi * lat_dependence)

    ! Spencer Hill 2017-10: break equatorial symmetry (if desired) by
    ! flattening temperatures poleward of the maximum and throughout the
    ! opposite hemisphere.
    if (e95_break_eq_symm) then

       ! Spencer Hill, 2017-12-07: to give rise to a nonzero gradient at the
       ! equator, switch from the E95 profile to its tangent within a given
       ! distance of the equator.  The intent of this modification is to
       ! prevent large equatorial jumps of the cross equatorial Hadley cell,
       ! c.f. the arguments by Pauluis (2004) JAS.
       cos_lat_switch = cos(real(e95_lat_switch_tangent*DEG_TO_RAD))

       ! Spencer Hill 2017-12-09: multiplying slope of tangent line by >1
       ! factor to get stronger cross equatorial gradient, to further suppress
       ! the equatorial jump.
       tangent_multiplier = 1.

       sfc_temp_at_switch = e95_max_temp * exp(-0.5 * chi * &
            (cos_phi0**2 - cos_lat_switch**2)**2 / cos_lat_switch**2)

       sfc_temp_deriv = chi * sfc_temp * sin(lat)*cos_lat * &
            (1. - (cos_phi0 / cos_lat)**4)

       deriv_at_switch = chi * sfc_temp_at_switch * &
            sin(e95_lat_switch_tangent*DEG_TO_RAD) * cos_lat_switch * &
            (1. - (cos_phi0 / cos_lat_switch)**4) * tangent_multiplier

       tangent_at_eq = sfc_temp_at_switch + &
            deriv_at_switch * -1 * e95_lat_switch_tangent*DEG_TO_RAD

       ! Spencer Hill 2018-07-20: allow Emanuel 1995 solution to extend
       ! poleward of the forcing maximum by a given latitudinal extent.
       lat_past_max = min(real(phi0 + e95_lat_past_max), 90.)
       cos_lat_past_max = cos(lat_past_max*DEG_TO_RAD)

       sfc_temp_past_max = e95_max_temp * exp(-0.5 * chi * &
            (cos_phi0**2 - cos_lat_past_max**2)**2 / cos_lat_past_max**2)

       do j = 1, size(lat, 2)
          tangent_at_switch(:, j) = sfc_temp_at_switch + deriv_at_switch * &
               (lat(:, j) - e95_lat_switch_tangent*DEG_TO_RAD)

          if (lat_deg(1, j) < -1*e95_lat_switch_tangent) then
             ! flat south of the switch latitude
             sfc_temp(:, j) = tangent_at_eq - &
                  (sfc_temp_at_switch - tangent_at_eq)
          elseif (lat_deg(1, j) < e95_lat_switch_tangent) then
             ! linear within the equatorial band"
             sfc_temp(:, j) = tangent_at_switch(:, j)
          elseif (lat_deg(1, j) > lat_past_max) then
             ! flat poleward of the (maximum + offset)
             sfc_temp(:, j) = sfc_temp_past_max
          endif
       enddo

    else
       ! Spencer Hill 2018-07-20: prevent surface forcing temperature from
       ! dropping below the specified stratospheric temperature, since
       ! otherwise these profiles will approach absolute zero at some latitude.
       sfc_temp = max(sfc_temp, t_strat)

    endif

  end function sfc_temp_emanuel95_crit

end module grid_physics
