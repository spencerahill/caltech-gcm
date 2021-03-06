module betts_miller_mod
  
  use		 constants_mod, only:    hlv,      hls,      pi             &
                                        ,cp,       rdgas,    dens_h2o       &
                                        ,tfreeze,  rvgas,    grav
  use		 utilities_mod, only:    open_file,          close_file     &
                                        ,check_nml_error
  use		       mpp_mod, only:    mpp_pe,   mpp_root_pe
  use                  fms_mod, only:    error_mesg,         stdlog         &
                                        ,write_version_number
  use             grid_physics, only:    surface_temperature
  use         time_manager_mod, only:    time_type
  use         diag_manager_mod, only:    register_diag_field, send_data
implicit none


!     namelist variables
!     -------- ---------
!
! typical values for timescales are (see comments in original code below)
!	      CBMTS =	14400 sec (all resolutions)
!		      { 14400 sec, n<=T21
!	      CBMTD = {	 7200 sec, T21<n<=T42
!		      {	 3600 sec, n>T42
!	      CBMEF =	 0.15	  (all resolutions).
!
  real :: cbmts = 14400.0  ! Adjustment timescale in seconds for shallow convection.
  real :: cbmtd = 14400.0  ! Adjustment timescale in seconds for deep convection.
  real :: cbmef = 0.15	  ! Downdraught efficiency for deep convection.
  real :: Ahlv  = 0.      ! Alternate latent heat of vaporization (can be 0)
  real :: Ahls  = 0.      ! Alternate latent heat of sublimation (can be 0)
  namelist /betts_miller_nml/ cbmts, cbmtd, cbmef, Ahlv, Ahls

  character(len=14) :: mod_name = 'betts_miller'
  character(len=128) :: version='$Id: betts_miller.f90 $'
  character(len=128) :: tag='homemade'

  real    :: dcvgr   !   global convective generation of rain. (not 
  real    :: dcvgs   !   global convective generation of snow.
  real    :: dcvmoi  !   global convective environmental moistening.

!     local constants 
!     ----- ---------
!
! constants for calculation of SVP
  real ::  vtmpc1,   vtmpc2,   c1es,     c2es,     c3ies,    c3les          &
          ,c4ies,    c4les,    c5les,    c5ies

  real, allocatable, dimension(:) ::     pk,       bk,       bk_full

  integer ::         num_levels,         lon_max

! diagnostic IDs
  integer ::         id_tdt,   id_qdt,   id_ztc,   id_dtd,   id_dqd,        &
           id_dts,   id_dqs,   id_dztref,          id_sztref,               &
           id_top,   id_base,  id_dztref0,         id_dztref1

  real :: missing_value = -1.e10

contains

!#######################################################################
  
  subroutine betts_miller_init(axes, Time, num_lats_in_subwindow,           &
       local_lon_max, bk_local)
    
    real, dimension(:), intent(in) :: bk_local
    integer, intent(in) :: axes(4)
    type(time_type), intent(in) :: Time

    integer, intent(in) :: local_lon_max, num_lats_in_subwindow

    integer :: ierr, namelist_unit, io
    
    num_levels = size(bk_local)-1
    lon_max = local_lon_max
    
    allocate (	      bk(num_levels+1) )
    allocate (	 bk_full(num_levels  ) )
    
    bk = bk_local
    bk_full=(bk(1:size(bk)-1)+bk(2:size(bk)))/2.0

    vtmpc1 = rvgas/rdgas-1.0
    c1es=610.78
    c2es=c1es*rdgas/rvgas
    c3les=17.269
    c3ies=21.875
    c4les=35.86
    c4ies=7.66
    c5les=c3les*(tfreeze-c4les)
    c5ies=c3ies*(tfreeze-c4ies)
    
    namelist_unit = open_file ('input.nml', action='read')
    ierr=1
    do while (ierr /= 0)
       read  (namelist_unit, nml=betts_miller_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'betts_miller_nml')
    enddo
10  call close_file (namelist_unit)
    
    call write_version_number (version,tag)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=betts_miller_nml)

! register fields with diagnostic manager
    id_tdt = register_diag_field ( mod_name, 'tdt_betts_miller', axes(1:3), &
         Time, 'convective temperature tendency', 'deg_K/sec',            &
         missing_value=missing_value)

    id_qdt = register_diag_field ( mod_name, 'qdt_betts_miller', axes(1:3), &
         Time, 'convective humidity tendency', 'deg_K/sec',            &
         missing_value=missing_value)

    id_ztc = register_diag_field ( mod_name, 'ztc_betts_miller', axes(1:3), &
         Time, 'parcel lift temperature', 'deg_K',            &
         missing_value=missing_value)

    id_dtd = register_diag_field ( mod_name, 'tdtd_betts_miller', axes(1:3), &
         Time, 'deep scheme temperature tendency', 'deg_K',            &
         missing_value=missing_value)

    id_dts = register_diag_field ( mod_name, 'tdts_betts_miller', axes(1:3), &
         Time, 'shallow scheme temperature tendency', 'deg_K',            &
         missing_value=missing_value)

    id_dztref0= register_diag_field ( mod_name, 'dztref_betts_miller0',      &
         axes(1:3), Time, 'deep scheme relaxation temperature', 'deg_K',     &
         missing_value=missing_value)

    id_dztref1= register_diag_field ( mod_name, 'dztref_betts_miller1',      &
         axes(1:3), Time, 'deep scheme relaxation temperature', 'deg_K',     &
         missing_value=missing_value)

    id_dztref= register_diag_field ( mod_name, 'dztref_betts_miller',        &
         axes(1:3), Time, 'deep scheme relaxation temperature', 'deg_K',     &
         missing_value=missing_value)

    id_sztref= register_diag_field ( mod_name, 'sztref_betts_miller',        &
         axes(1:3), Time, 'shallow scheme relaxation temperature', 'deg_K',  &
         missing_value=missing_value)

    id_top = register_diag_field ( mod_name, 'cloud_top_betts_miller',        &
         axes(1:2), Time, 'cloud_top', 'level_index',  &
         missing_value=missing_value)

    id_base = register_diag_field ( mod_name, 'cloud_base_betts_miller',      &
         axes(1:2), Time, 'cloud_base', 'level_index',  &
         missing_value=missing_value)


    dcvgr = 0.0
    dcvgs = 0.0
    dcvmoi = 0.0

  end subroutine betts_miller_init

!#######################################################################

  subroutine betts_miller(   is,       ie,       js,       je,              &
                 Time,  delta_t,        tg_previous,                        &
                    qg_previous,   p_half,   p_full,    dt_tg,              &
                          dt_qg )
    
    integer, intent(in)                    ::  is, ie, js, je
    type(time_type), intent(in)            ::  Time
    real, intent(in)                       ::  delta_t
    real, intent(in), dimension(:,:,:)     ::  tg_previous, qg_previous
    real, intent(in), dimension(:,:,:)     ::  p_half, p_full
    real, intent(inout), dimension(:,:,:)  ::  dt_tg, dt_qg
    
    ! local variables	(generally output from betts-miller scheme)
    real, dimension(size(dt_tg,1), size(dt_tg,2)) :: surf_temperature  ! surf. temp. (t-1)
    real, dimension(size(dt_tg,1), size(dt_tg,2)) :: rsfc      !	surface convective rainfall rate
    real, dimension(size(dt_tg,1), size(dt_tg,2)) :: ssfc      !	surface convective snowfall rate
    real, dimension(size(dt_tg,1), size(dt_tg,2)) :: aprc      !	accum. convective precip.
    real, dimension(size(dt_tg,1), size(dt_tg,2)) :: aprs      !	accum. convective snowfall
    real, dimension(size(dt_tg,1), size(dt_tg,2)) :: slmm      ! land/sea mask (1/0)
    logical, dimension(size(dt_tg,1), size(dt_tg,2)) :: lpland ! logical array for land points

    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: ztc  ! parcel lift temperature
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: dtd  ! deep temperature tendency
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: dqd  ! deep humidity tendency
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: dts  ! shallow temperature tendency
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: dqs  ! shallow humidity tendency
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: sztref  ! shallow reference temperature
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: dztref0 ! deep reference temperature
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: dztref1 ! deep reference temperature
    real, dimension(size(dt_tg,1), size(dt_tg,2), size(dt_tg,3)) :: dztref ! deep reference temperature
    integer, dimension(size(dt_tg,1), size(dt_tg,2)) :: ibase ! cloud base
    integer, dimension(size(dt_tg,1), size(dt_tg,2)) :: itop  ! cloud top

    logical :: lprad       !   switch for cloud-radiation diagnostics.
    integer :: neadjtop    !   imposed top level for betts-miller scheme    
    integer :: krow        !   latitude index
    logical :: used
    
    slmm = 0.
    lpland = .false.
    neadjtop = 1

    ssfc = 0.0
    rsfc = 0.0

    ztc = 0.0
    dtd = 0.0
    dqd = 0.0
    dts = 0.0
    dqs = 0.0
    dztref0 = 0.0
    dztref1 = 0.0
    dztref = 0.0
    sztref = 0.0
    itop = 0
    ibase = 0

    surf_temperature = surface_temperature(tg_previous, p_full, p_half)

    do krow = 1, size(tg_previous,2)
       call bmadj(                                                          &
            cbmts,	  cbmtd,         cbmef,	        neadjtop,           &
	    lon_max,	  num_levels,	 num_levels-1,  num_levels+1,       &
            bk_full,      bk,            delta_t,                           &
            krow,         Ahls,	         Ahlv,		pi,	            &
            cp,	          grav,	         rdgas,		dens_h2o,           &
            tfreeze,      vtmpc1,        vtmpc2,	c2es,	            &
            c3ies,	  c3les,         c4ies,		c4les,	            &
            c5les,        dcvgr,         dcvgs,         dcvmoi,	            &
            tg_previous(:, krow, :),     qg_previous(:, krow, :),           &
            dt_tg(:, krow, : ),		 dt_qg(:, krow, : ),	            &
            p_full(:, krow, : ),	 p_half(:, krow, : ),	            &
            slmm(:, krow),               surf_temperature(:, krow),         &
            rsfc(:, krow),	         ssfc(:, krow),                     &
	    lpland(:,krow),              aprc(:, krow),                     &
            aprs(:, krow),               ztc(:, krow, :),                   &
            dtd(:, krow, :),             dqd(:, krow, :),                   &
            dts(:, krow, :),             dqs(:, krow, :),                   &
            dztref0(:,krow,:),           dztref1(:, krow, :),               &
            dztref(:,krow,:),            sztref(:, krow, :),                &
            ibase(:,krow),               itop(:,krow)         )
    end do

! send data to diagnostic manager
    if(id_tdt > 0 ) used = send_data ( id_tdt, dt_tg, Time, is, js)
    if(id_qdt > 0 ) used = send_data ( id_qdt, dt_qg, Time, is, js)
    if(id_ztc > 0 ) used = send_data ( id_ztc, ztc,   Time, is, js)
    if(id_dtd > 0 ) used = send_data ( id_dtd, dtd,   Time, is, js)
    if(id_dts > 0 ) used = send_data ( id_dts, dts,   Time, is, js)
    if(id_dztref0 > 0 ) used = send_data ( id_dztref0, dztref0, Time, is, js)
    if(id_dztref1 > 0 ) used = send_data ( id_dztref1, dztref1, Time, is, js)
    if(id_dztref > 0 ) used = send_data ( id_dztref, dztref, Time, is, js)
    if(id_sztref > 0 ) used = send_data ( id_sztref, sztref, Time, is, js)
    if(id_base > 0 ) used = send_data ( id_base, float(ibase), Time, is, js)
    if(id_top  > 0 ) used = send_data ( id_top , float(itop), Time, is, js)

  end subroutine betts_miller

!*/ *********************************************************************
!*IDENT NDEEP4V
!*/    Update to UGAMP GCM version 2.0 : ugcm2.npl
!*/    Betts-Miller convective adjustment scheme.
!*/    N.B.  SUBROUTINES ONLY.	MUST BE USED WITH INTERFACE MODSET
!*/    MB940412, IN FILE "ugcm2_bmint2".
!*/ ---------------------------------------------------------------------
!*/    Existing comdecks modified:  - none -
!*/    Existing decks modified:	    - none -
!*/    New comdecks added:	    - none -
!*/    New decks added:		    BMDOC,BMINT,BMADJ,BMDEEP,BMSHAL
!*/ ---------------------------------------------------------------------
!*/
!*/ *************************************************************** BMDOC
!*/

  SUBROUTINE BMDOC

!
!**** *BMDOC* - DOCUMENTATION DUMMY ROUTINE FOR BETTS-MILLER CONVECTIVE
!		ADJUSTMENT SCHEME.
!
!     MIKE BLACKBURN		  U.G.A.M.P.		     31/05/94.
!
!  ---------------------------------------------------------------------
!
!     References.
!     -----------
!
!     The Betts-Miller convective adjustment scheme included in the
!     UGAMP GCM version 2.1 is based on code written by Alan Betts and
!     Martin Miller.  It contains essentially technical changes from
!     their supplied version.  Details are given later.
!
!     The basic method of the original scheme is described in
!
!     Betts, A.K.,  1986  Q.J.R. Meteorol. Soc., 112, 677-691.
!
!     Single column tests of the original scheme are described in
!
!     Betts, A.K. & M.J.Miller,	 1986  Q.J.R. Meteorol. Soc., 112,
!     693-709.
!
!     The current scheme, including an explicit parametrization of
!     downdraughts, is described by
!
!     Betts, A.K. & M.J.Miller,	 1994  American Met. Soc. Monograph
!     on Convective Parametrizations.
!
!     Implementation of the scheme in the UGCM is described by
!
!     Slingo, J.M. & M. Blackburn,  1992   UGAMP Technical Report 25.
!
!     The calculations for the saturation point (temperature of the
!     lifting condensation level) follow a functional fit derived by
!
!     Bolton, D.,  1980	 Mon. Wea. Rev., 108, 1046-1053.
!
!  ---------------------------------------------------------------------
!
!     Outline Documentation.
!     ----------------------
!
!     The Betts-Miller scheme contains the following subroutines:
!
!     BMINT:  Interface routine to the plug-compatible routine BMADJ.
!	      BMINT is set up for the UGAMP GCM.  It extracts common
!	      block variables and passes them in the argument list to
!	      BMADJ.  The memory manager is used to locate arrays in
!	      long term storage and these are similarly passed through
!	      the argument list.  By-pass processing has been extracted
!	      into this routine.
!
!     BMADJ:  Find conditionally unstable layers and select for deep or
!	      shallow moist convection.
!	      Create zonal diagnostics and parameters for diagnostic
!	      cloud scheme.
!	      Compute CAPE statistics, numbers of convecting points.
!
!     BMDEEP: Perform precipitating deep convective adjustment.
!	      Return temperature and moisture tendencies, and surface
!	      rainfall and snowfall rates.
!
!     BMSHAL: Perform non-precipitating shallow convective adjustment.
!	      Return temperature and moisture tendencies, and condens-
!	      ation rate for the cloud scheme.
!
!     The following tuneable parameters are provided:
!
!     CBMTS:  Adjustment timescale in seconds for shallow convection.
!	      This timescale may be weakly resolution dependent.  It
!	      must be chosen to maintain a balance between turbulent
!	      fluxes, grid-scale descent at inversion level and the
!	      convection.
!
!     CBMTD:  Adjustment timescale in seconds for deep convection.
!	      This timescale is resolution dependent and must be
!	      sufficiently short to prevent grid-scale saturation in
!	      the regions of strongest vertical motion.
!
!     CBMEF:  Downdraught efficiency for deep convection.  This is a
!	      fraction, in the range zero to unity.  It represents the
!	      proportion of precipitation generated by the adjustment
!	      above cloud base which is evaporated in the downdraught.
!	      It effectively defines the adjustment timescale in the
!	      planetary boundary layer for deep convection.
!
!     BMINT extracts values for these parameters from common, so they
!     must be defined in initialisation routines.  In the UGAMP GCM,
!     default values are defined in routine INIPHY and may be overridden
!     in namelist PHYSCTL.  The defaults are:
!
!	      CBMTS =	14400 sec (all resolutions)
!		      { 14400 sec, n<=T21
!	      CBMTD = {	 7200 sec, T21<n<=T42
!		      {	 3600 sec, n>T42
!	      CBMEF =	 0.15	  (all resolutions).
!
!     Note that these parameters have been tested only at T42L19 and
!     T21 resolutions in the UGAMP GCM.
!
!  ---------------------------------------------------------------------
!
!     Cautions.
!     ---------
!
!     1. This code has been tested and tuned only at T42L19 resolution
!	 in the UGAMP GCM, and has been tested at T21 horizontal resol-
!	 ution.
!
!     2. The scheme contains level dependent code which has only been
!	 tested with the standard 19 levels.
!	 a) The crossover between shallow and deep convection occurs at
!	    specific sigma (p/ps) values, so is essentially independent
!	    of model level placement.
!	 b) The downdraught part of the deep scheme assumes that the
!	    lowest three model levels constitute the planetary boundary
!	    layer.  This is checked but cannot be changed easily while
!	    retaining the downdraught component.
!	 c) The deep scheme contains a level-dependent definition to
!	    separate "low-lift" points from mid-level convection.  This
!	    is related to the downdraught parametrization in (b).
!	 d) All points up to a model level in the mid/upper troposphere,
!	    defined by a cut-off hybrid coordinate value, are tested for
!	    possible parcel ascents.  This is designed to limit deep
!	    convection to the troposphere and is essentially independent
!	    of model level placement.
!	 e) The scheme can be limited to operate below a prescibed model
!	    level, allowing it to be switched off above the tropopause.
!	    A warning is given if any parcel ascents are buoyant beyond
!	    this cut-off level.	 This option is largely redundent, given
!	    (d) above.
!
!     3. Note that routine BMINT here does not copy or update the
!	 diagnostic arrays containing zonally averaged tendencies from
!	 deep convection for wind components and KE (the NDU/V/ECUML/S
!	 arrays in AZDIA), since those particular arrays in UGCM 2.1
!	 are used by the gravity wave drag scheme after the convection
!	 call.
!
!     4. The Tetens formula used for the calculation of saturation
!	 specific humidity becomes invalid for pressures less than the
!	 saturated vapour pressure and can give large negative Qsat in
!	 these cases.  A correction has been applied to prevent this
!	 by modifying the ZCOR term from  ZCOR=1./(1.-VTMPC1*ZQSATC)
!	 to become  ZCOR=1./MAX(1.E-10,ABS(1.-VTMPC1*ZQSATC)) .
!	 The problem occurs only in the middle-atmosphere, generally at
!	 heights above 10hPa.  This should not be a problem with the
!	 current cut-off vertical coordinate value described in (3c)
!	 above.	 See UGAMP Technical Report No.5 for details.
!
!     5. Whilst an attempt has been made to extract all machine-specific
!	 code into the interface routine BMINT, the following Cray-
!	 specific code remains within the main subroutines:
!	 a) Calls to several Cray library routines.
!	 b) Work space local to each subroutine is allocated using the
!	    memory manager of the UGAMP/ECMWF models and Cray pointers.
!	    It will be possible to replace this with ANSI code in
!	    Fortran 90.
!
!  ---------------------------------------------------------------------
!
!     History.
!     --------
!
!     1. Code received from Martin Miller, ECMWF, May 1991.
!     2. Interface modset and routine interfaces modified for UGAMP
!	 GCM (UGCM version 1.2), May 1991.
!     3. Scheme tuned at T42L19 in UGCM, in collaboration with
!	 Martin Miller & Alan Betts, May - Dec 1991.  This involved
!	 the following changes (see UGAMP Technical Report No 25):
!	 a) The base temperature used to begin the deep first guess
!	    reference profile was modified from the environment
!	    temperature at the nominal cloud base to the average of
!	    the environment and cloud ascent temperatures at the
!	    same level.
!	 a) The deep scheme was modified to improve the reference
!	    profile structure above the freezing level, while retaining
!	    the near surface cooling and drying, to reduce the cloud
!	    top heating/cooling dipole.	 This involved changing from a
!	    linear to quadratic interpolation of the adiabatic deficit
!	    with pressure above freezing level.
!	 b) The relaxation timescales were modified to 2hours (deep)
!	    and 4hours (shallow).
!	 c) The downdraft efficiency parameter was reduced from 0.25
!	    to 0.15.
!     4. Used in the UGAMP AMIP integration, Jan 1992, as "NDEEP2F".
!
!  ---------------------------------------------------------------------
!
!     Further modification history:
!     -----------------------------
!
!     1.  Original NDEEP2F 12.91 consists of MJM/AKB code modified
!	  by JMS for UGCM.  Includes various switches and options
!	  tuned for T42 resolution, Autumn 1991.
!     2.  NDEEP3  : tidied version of code.		 JMS	04.92
!     3.  NDEEP3B : Logic modified to print CAPEs every timestep,
!		    (temporary).			  MB 14.05.92
!     4.  NDEEP4A : Modify modset only (interface), not decks. 18.5
!	  Note : Should be used with EPHYS3 modset, which includes
!	  code to set up NEADJTOP (not yet used here).
!     5.  NDEEP4B : Tidy deck CONVECT to end of sect 1.	  MB 21.05.92
!     6.  NDEEP4C : Reorder work array in CONVECT.	  MB 21.05.92
!     7.  NDEEP4D : Tidy CONVECT setion 2.		  MB 22.05.92
!     8.  NDEEP4E : Tidy CONVECT setion 3.		  MB 25.05.92
!     9.  NDEEP4F : Move CONVECT 4.1 to 3.7 and modify
!		    calc. of CAPE diagnostics (de-Cray).  MB 26.05.92
!     10. NDEEP4G : Reorganise remainder of CONVECT.	  MB 26.05.92
!     11. NDEEP4H : Skip CONVECT sects if no convection.  MB 03.06.92
!     12. NDEEP4I : Limit vertical loops in CONVECT.	  MB 03.06.92
!     13. NDEEP4J : Tidy up memory etc in DEEP to sect.1  MB 05.06.92
!     14. NDEEP4K : Fix ZCAPE bug in 371 loop (>= NDEEP4F).
!		    New multi-tasking code for CAPE diags MB 15.12.93
!     15. NDEEP4L : Full precision ALOG, EXP, but retain
!		    the XLG function.
!		    Fix bug for ZTREF in DEEP loop 222.
!		    Remove novector in DEEP loop 2201.	  MB 16.12.93
!     16. NDEEP4M : Remove Cray vector merges in CONVECT.
!		    Small mods to logic.		  MB 17.12.93
!     17. NDEEP4N : Remove Cray vector merges in DEEP.
!		    Tidy, resection, renumber.
!		    Separate loops for PBL integrals.	  MB 04.01.94
!     18. NDEEP4O : Remove Cray vector merges in SHALLOW.
!		    Tidy, resection, renumber.
!		    Remove unnecessary work arrays.
!		    Merge single-level loops.		  MB 06.01.94
!     19. NDEEP4P : Check top level of shallow tendency.
!		    Generalise level-specific IETLAB, NPBL.
!		    Fix bug in 5.5 for shallow AZDIA arrays
!		    which crept in in NDEEP4G(?).	  MB 12.01.94
!     20. NDEEP4Q : Include Qsat correction for p < SVP.  
!		    Modify by-pass processing in CONVECT.
!		    Improved vectorisation in CONVECT.	  MB 21.01.94
!     21. NDEEP4R : Revised CAPE diagnostics: normalisation
!		    corrected: count failed points.	  MB 18.03.94
!     22. NDEEP4S : Renumber/reorder DEEP sect.5. and
!		    optimise evaporation calcs.
!		    Arithmetic changes to optimise DEEP.  MB 30.03.94
!     23. NDEEP4T : Change names of routines, main switch.
!		    Timescales and downdraught efficiency
!		    become namelist variables.
!		    Checks for imposed top corrected.	  MB 18.04.94
!		    Further optimisation of BMADJ,BMDEEP. MB 20.04.94.
!		    By-pass processing corrected.	  MB 25.04.94.
!		    Computational constants changed/added.MB 25.04.94.
!	     4T12 : Correct IHITOP calculations in BMADJ. MB 27.04.94.
!	     4T13 : Avoid ILAB(i,0) use in BMADJ 3.5.	  MB 28.04.94.
!	     4T14 : Modify global diagnostics in BMADJ 6. MB 28.04.94.
!	     4T15 : Restrict BMDEEP 361 calcs to cloud.	  MB 29.04.94.
!	     4T16 : Replace redundent vector work arrays
!		    by scalars in BMSHAL.		  MB 29.04.94.
!     24. NDEEP4U : *CALL only required common blocks.	  MB 06.05.94.
!	     4U2  : No task info in DEEP/SHAL.	No DTIME. MB 09.05.94.
!	     4U3  : Plug-compatible routines BMDEEP/SHAL. MB 10.05.94.
!	     4U4  : Plug-compatible main routine BMADJ
!		    & new interface routine BMINT.	  MB 11.05.94.
!	     4U5  : Resolution checks extracted to BMINT,
!		    dummy argument name changes.	  MB 13.05.94.
!     25. NDEEP4V : By-pass processing extracted to BMINT MB 31.05.94.
!
!  ---------------------------------------------------------------------
!
!     Changes to results from NDEEP2F to current code:
!     ------------------------------------------------
!
!     1. The half precision exponential and log functions in NDEEP2F
!	 have been replaced by ANSI Fortran functions.
!     2. NDEEP2F contained a known bug, whereby the deep first guess
!	 reference temperature at the nominal cloud base was incorrectly
!	 modified at each pass of the loop over levels from the base to
!	 the freezing level.  This did not affect the reference profile
!	 above the base level.	The second line in (BM)DEEP loop 222:
!	    ZTMEAN=0.5*(ZTCD(JL,NUPBL)+ZTREF(JL,NUPBL))
!	 was replaced in NDEEP4L by:
!	    ZTMEAN=0.5*(ZTCD(JL,NUPBL)+ZTP1D(JL,NUPBL))
!	 Subsequently in NDEEP4N a separate loop (342) is executed for
!	 the reference profile at level NUPBL.
!     3. Changes to the arithmetic ordering in BMDEEP for optimisation
!	 purposes give rounding differences.  The exact changes are:
!	 a)  Regroup terms in the PBL evaporation integral in loop 411.
!	 b)  Multiply the ZFG factor by the latent heat in the single
!	     level loop 431 rather than nested loop 441.
!	 c)  Changes to the evaporation calculation in loop 531:
!	     remove redundent ELSE part of if-block for ZDQEV:
!	     reorder terms in expression for ZZDSP:
!	     introduce ZRINCM to avoid (ZRINC()-ZDR()) terms:
!	     use multiply in place of 2nd divide for ZDQEV.
!     4. The correction to the Qsat formula, for pressures less than the
!	 SVP, would alter the results if any cloud parcels encountered
!	 this condition.  It only occurs well above the tropopause.
!     5. The CAPE diagnostics, but not the meteorological results,
!	 differ following rewriting the macro-tasking code in BMADJ,
!	 in NDEEP4K and further, following the revisions in NDEEP4R to
!	 the normalisation (inclusion of failing non-zero CAPE points).
!
!  ---------------------------------------------------------------------
!
!     CAPE Diagnostics & logic of DEEP/SHALLOW switching.
!     ---------------------------------------------------
!
!     Subroutine BMADJ computes Convective Available Potential Energy
!     (CAPE) statistics at each latitude and globally, every 12 hours.
!     The logic of selecting deep and shallow moist convecting points
!     excludes two classes of convecting layer which are diagnosed as
!     having non-zero CAPE, causing the statistics in the original code
!     to be incorrectly normalised.  The errors were typically 1%
!     globally but were substantial at individual latitudes.  The
!     normalisation now includes the failing points, which are also
!     counted separately.
!
!     BMADJ section [3] tests all columns at the current latitude for
!     dry and moist conditional instability.  The CAPE is computed for
!     all moist unstable points and is set to zero for other points.
!     Points with CAPE below a pre-defined cut-off (which is negative)
!     and shallow unstable layers at the surface are excluded and have
!     CAPE reset to zero.
!
!     BMADJ section [4] selects points for deep convection.  These
!     require non-zero CAPE, a deep cloud-top and must include at least
!     3 levels from parcel origin to cloud top.	 Some points fail to
!     give positive precipitation, due to the convecting layer being
!     drier than the reference profile.	 These are excluded from the
!     count of deep convecting points and are swapped for possible
!     shallow convection.
!
!     BMADJ section [5] selects points for shallow convection.	These
!     require non-zero CAPE, zero or negative precipitation (i.e. swap
!     points can be included) and cloud-base must be below the imposed
!     shallow convective top.
!
!     Finally the CAPE statistics are computed in BMADJ section [6].
!     All points having non-zero CAPE are included in the zonal mean
!     and the land/sea averages.
!
!     Two classes of moist unstable point having non-zero diagnosed
!     CAPE fail to give successful convection, given the current logic
!     for deep and shallow selection.  They are
!
!     a)  2-level layers having a high cloud base,
!
!     b)  Swap points having a high cloud base.
!
!     All such points are mid-level convecting layers.	The former fail
!     the deep selection, having too few levels.  The latter are too
!     dry to give deep convection but cannot swap successfully.
!
!     In the original code, these failing points were included in the
!     CAPE sums but excluded from the normalisation and the count of
!     successful convecting points.  Such points typically account for
!     1% of the total globally, but could even dominate the statistics
!     at individual latitudes.	They are now included in the CAPE
!     normalisation and are printed in a separate count, together with
!     their average CAPE at each latitude.
!
!     An additional feature of the CAPE statistics in the original code
!     was that the zonal average could be larger than the maximum value.
!     This generally occurred when all non-zero CAPE points failed to
!     convect.	The zonal sums were then normalised by unity.  This can
!     no longer occur and rows with no convecting points correctly show
!     zero CAPE.
!
!     The deep-shallow selection logic also treats a further class of
!     points in a non-trivial way.  Points having a deep top (defined
!     as p/ps < 0.725 over land and p/ps < 0.810 over sea) but only 2
!     levels fail the deep test.  However those having cloud base below
!     the imposed shallow top (p/ps=0.725 at both land and sea points)
!     will be selected successfully for shallow convection.
!
!  ---------------------------------------------------------------------
!*ENDIF

    write(*,*) ' *** DUMMY ROUTINE BMDOC FOR DOCUMENTATION OF' &
	 ,' BETTS-MILLER CONVECTION SCHEME ***'
    
    RETURN
  END SUBROUTINE BMDOC
!*/
!*/ *************************************************************** BMADJ
!*/
   SUBROUTINE BMADJ(							 &
	  CBMTS,    CBMTD,    CBMEF,	NEADJTOP			 &
	 ,NLON,	    NLEV,     NLEVM1,	NLEVP1,	  CETA	                 &
	 ,CETAH,    TWODT,    KROW,	ALS,      ALV                    &
         ,API,      CPD,      G,	RD,	  RHOH2O	         &
	 ,TMELT,    VTMPC1,   VTMPC2,	C2ES,	  C3IES,    C3LES	 &
	 ,C4IES,    C4LES,    C5LES,	DCVGR	 &
	 ,DCVGS,    DCVMOI,   TM1,      QM1		 &
	 ,TE,	    QE,	      APP1,	APHP1,	  SLMM,	    TSM1M	 &
	 ,RSFC,	    SSFC,     LPLAND,	APRC,	  APRS,	    ztc          &
         ,dtd,      dqd,      dts,      dqs,      dztref0,  dztref1      &
         ,dztref,   sztref,   ibase,    itop)

     integer, intent(in) ::   neadjtop, nlon,     nlev,     nlevm1      &
         ,nlevp1,   krow 

     real, intent(in)    ::   cbmts,     cbmtd,   cbmef,    twodt        &
         ,als,      alv,      api,       cpd,     g,        rd           &
         ,rhoh2o,   tmelt,    vtmpc1,    vtmpc2,  c2es,     c3ies        &
         ,c3les,    c4ies,    c4les,     c5les

     real, intent(inout) ::   dcvgr,     dcvgs,   dcvmoi
!
!**** *BMADJ* - COMPUTES T AND Q TENDENCIES DUE TO DEEP AND SHALLOW
!		CONVECTION.
!
!     ORIGINAL VERSION		     M.J.MILLER	 E.C.M.W.F.  16/07/84.
!     VECTORIZED		       B.RITTER	 E.C.M.W.F.	12/84.
!     IMPLEMENTED & TUNED IN UGCM1   J.M.SLINGO	 U.G.A.M.P.	12/91.
!     IMPLEMENTED IN UGCM2	    M.BLACKBURN	 U.G.A.M.P.  31/05/94.
!
!     PURPOSE.
!     --------
!
!	   THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE TWO
!     PROGNOSTIC VARIABLES T AND Q DUE TO MOIST CONVECTIVE PROCESSES.
!     BOTH SHALLOW (NON-PRECIPITATING) CLOUDS AND DEEP (PRECIPITATING)
!     CLOUDS ARE INCLUDED VIA THE ADJUSTMENT TOWARDS REFERENCE PROFILES.
!     PRECIPITATION IS EITHER AS SNOW OR AS LIQUID WATER DEPENDING ON
!     TEMPERATURE.  MELTING (OR FREEZING) OF THE FALLING WATER IS ALSO
!     PERMITTED IN THE PARTITIONING OF SURFACE ACCUMULATIONS.
!
!**   INTERFACE.
!     ----------
!
!	   *BMADJ* IS CALLED FROM *BMINT*.
!	   THE ROUTINE IS ARGUMENT DRIVEN, TAKING ITS INPUT ENTIRELY
!     FROM THE VARIABLES AND ARRAYS SUPPLIED IN THE ARGUMENT LIST:
!     T,Q AND P AT FULL LEVELS, HALF LEVEL PRESSURE AND SURFACE FIELDS.
!     IT RETURNS ITS OUTPUT VIA THE ARGUMENT LIST:
!     MODIFIED TENDENCIES OF T AND Q, CONVECTIVE RATES OF SURFACE
!     PRECIPITATION FOR RAIN AND SNOW, VARIOUS DIAGNOSTICS.
!
!  ---------------------------------------------------------------------
!
!	   INPUT ARGUMENTS:
!
!     *CBMTS*	ADJUSTMENT TIMESCALE IN SECONDS FOR SHALLOW CONVECTION.
!     *CBMTD*	ADJUSTMENT TIMESCALE IN SECONDS FOR DEEP CONVECTION.
!     *CBMEF*	DOWNDRAUGHT EFFICIENCY PARAMETER FOR DEEP CONVECTION.
!     *NEADJTOP*IMOPSED TOP LEVEL FOR BETTS-MILLER SCHEME.
!     *NLON*	NUMBER OF LONGITUDES.
!     *NLEV*	NUMBER OF LEVELS.
!     *NLEVM1*	(NLEV-1).
!     *NLEVP1*	(NLEV+1).
!     *CETA*	HYBRID COORDINATE AT FULL LEVELS,      DIMENSION (NLEV).
!     *CETAH*	HYBRID COORDINATE AT HALF LEVELS,    DIMENSION (NLEVP1).
!     *TWODT*	2* TIMESTEP IN SECONDS.
!     *KROW*	INDEX FOR CURRENT LATITUDE ROW.
!     *ALS*	LATENT HEAT FOR SUBLIMATION.
!     *ALV*	LATENT HEAT FOR VAPORISATION.
!     *API*	PI.
!     *CPD*	SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR.
!     *G*	GRAVITATIONAL ACCELERATION.
!     *RD*	GAS CONSTANT FOR DRY AIR.
!     *RHOH2O*	DENSITY OF LIQUID WATER.
!     *TMELT*	TEMPERATURE OF FUSION OF ICE.
!     *VTMPC1*	CONSTANT FOR VIRTUAL EFFECTS, (RV/RD-1).
!     *VTMPC2*	CONSTANT FOR VIRTUAL EFFECTS, (CPV/CPD-1).
!     *C__ES*	CONSTANTS USED FOR COMPUTATION OF SATURATION SPECIFIC
!		HUMIDITY OVER LIQUID WATER (C_LES) OR ICE (C_IES).
!     *C2ES*	(RD/RV)*(SVP AT REFERENCE TEMPERATURE C4_ES).
!     *C3_ES*	CONSTANT FOR SVP.
!     *C4_ES*	REFERENCE TEMPERATURE FOR SVP.
!     *C5_ES*	(C3_ES*(TMELT-C4_ES)).
!     *PBUDW*	CURRENT LATITUDE WEIGHT FOR GAUSSIAN INTEGRAL.
!     *PTWOMU*	2* SINE(GAUSSIAN LATITUDE).
!     *DCVGR*	ACCUM. GLOBAL CONVECTIVE GENERATION OF RAIN.
!     *DCVGS*	ACCUM. GLOBAL CONVECTIVE GENERATION OF SNOW.
!     *DCVMOI*	ACCUM. GLOBAL CONVECTIVE ENVIRONMENTAL MOISTENING.
!     *ND_CUM*	INDEX FOR MASK DIAGNOSTICS OF SENSIBLE/LATENT HEAT
!		BY CONVECTION.
!     *ND_CUM_* INDEX OF ZONAL AVERAGE T/Q TENDENCY FOR LAND/SEA
!		FOR DEEP CONVECTION.
!     *ND_SCV_* INDEX OF ZONAL AVERAGE T/Q TENDENCY FOR LAND/SEA
!		FOR SHALLOW CONVECTION.
!     *TM1*	(T-1) TEMPERATURE,		  DIMENSION (NLON,NLEV).
!     *QM1*	(T-1) SPECIFIC HUMIDITY,	  DIMENSION (NLON,NLEV).
!     *TE*	TEMPERATURE TENDENCY,		  DIMENSION (NLON,NLEV).
!     *QE*	SPECIFIC HUMIDITY TENDENCY,	  DIMENSION (NLON,NLEV).
!     *APP1*	FULL LEVEL PRESSURE,		  DIMENSION (NLON,NLEV).
!     *APHP1*	HALF LEVEL PRESSURE,		DIMENSION (NLON,NLEVP1).
!     *SLMM*	LAND/SEA MASK (1/0),		       DIMENSION (NLON).
!     *TSM1M*	(T-1) SURFACE TEMPERATURE,	       DIMENSION (NLON).
!     *RSFC*	ZERO PRESET,			       DIMENSION (NLON).
!     *SSFC*	ZERO PRESET,			       DIMENSION (NLON).
!     *LPLAND*	LOGICAL ARRAY FOR LAND POINTS,	       DIMENSION (NLON).
!     *APRC*	ACCUM. CONVECTIVE PRECIP.,	       DIMENSION (NLON).
!     *APRS*	ACCUM. CONVECTIVE SNOWFALL,	       DIMENSION (NLON).
!
!	   OUTPUT ARGUMENTS:
!
!     *DCVGR*	MODIFIED GLOBAL CONVECTIVE GENERATION OF RAIN.
!     *DCVGS*	MODIFIED GLOBAL CONVECTIVE GENERATION OF SNOW.
!     *DCVMOI*	MODIFIED GLOBAL CONVECTIVE ENVIRONMENTAL MOISTENING.
!     *TE*	MODIFIED TEMPERATURE TENDENCY,	  DIMENSION (NLON,NLEV).
!     *QE*	MODIFIED MOISTURE TENDENCY,	  DIMENSION (NLON,NLEV).
!     *RSFC*	SURFACE CONVECTIVE RAINFALL RATE,      DIMENSION (NLON).
!     *SSFC*	SURFACE CONVECTIVE SNOWFALL RATE,      DIMENSION (NLON).
!     *APRC*	MODIFIED ACCUM. CONVECTIVE PRECIP.,    DIMENSION (NLON).
!     *APRS*	MODIFIED ACCUM. CONVECTIVE SNOWFALL,   DIMENSION (NLON).
!
!  ---------------------------------------------------------------------
!*ENDIF
!     METHOD.
!     -------
!
!	   AN INDEX OF CONVECTIVELY ACTIVE POINTS IS CREATED, DEFINED
!     BY THE EXISTENCE OF LOW LEVEL MOIST ADIABATIC INSTABILITY.
!     DEEP AND SHALLOW CONVECTIVE POINTS (DISTINGUISHED BY THE HEIGHT
!     OF THE UNSTABLE SLAB AND THE OCCURRENCE OF PRECIPITATION) ARE
!     TREATED SEPARATELY.  SUBROUTINES *BMDEEP* AND *BMSHAL* PERFORM
!     THE ADJUSTMENT PROCESS.
!
!     EXTERNALS.
!     ----------
!
!	   *ALLOCA*  ALLOCATE ARRAY SPACE (MEMORY MANAGER ROUTINE).
!	   *BMDEEP*  PERFORM THE ADJUSTMENT FOR DEEP CONVECTIVE
!		     SITUATIONS.
!	   *MASKDIA* COMPUTE MASK (LIMITED AREA) DIAGNOSTICS.
!	   *OFFLOCK* ALLOW SUBSEQUENT CODE TO BE EXECUTED MULTI-TASKED.
!		     (CRAY MACRO-TASKING LIBRARY ROUTINE).
!	   *ONLOCK*  FORCE SUBSEQUENT CODE TO BE EXECUTED SINGLE-TASKED.
!		     (CRAY MACRO-TASKING LIBRARY ROUTINE).
!	   *BMSHAL*  PERFORM THE ADJUSTMENT FOR SHALLOW CONVECTIVE
!		     SITUATIONS.
!	   *UNLOC*   FREE ARRAY SPACE (MEMORY MANAGER ROUTINE).
!	   *WHENEQ*  SEARCH A VECTOR FOR ELEMENTS EQUAL TO A TARGET,
!		     TO GATHER CONVECTING POINTS (CRAY LIBRARY ROUTINE).
!
!     REFERENCES.
!     -----------
!
!	   A BASIC DESCRIPTION OF THE METHOD CAN BE FOUND IN BETTS
!     (1986 : QJRMS 112, 677-691) AND BETTS & MILLER (1986 : QJRMS,
!     112, 693-709).
!	   A MORE DETAILED DESCRIPTION OF THE SCHEME IS CONTAINED IN
!     BETTS & MILLER (1994 : AMERICAN METEOROLOGICAL SOCIETY MONOGRAPH
!     ON CONVECTIVE PARAMETRIZATION).
!	   IMPLEMENTATION OF THE SCHEME IN THE UGAMP GCM IS DESCRIBED
!     IN UGAMP TECHNICAL REPORT NO. 25 (SLINGO & BLACKBURN, 1992).
!
!	   SATURATION POINT CALCULATIONS FOLLOW BOLTON (1980, MON. WEA.
!     REV., 108, 1046-1053).  SATURATION SPECIFIC HUMIDITY CALCULATIONS
!     USE THE TETENS FORMULA (LOWE, 1977, J.APPL.MET., 16, 100-103).
!
!     ------------------------------------
      LOGICAL LO,LOA,LOB,LO2,LO3,LOIS,LOVS
      LOGICAL LOPRT
!     --------------------------
!
! global integrals won't work since code is running in parallel--fix some day
!!$	 SAVE NROWSUM
!!$	 SAVE NGCONV,NGSHAL,NGDEEP,NGSWAP,NGFAIL,NGTOT,GCAPE
!
      integer :: ndeep, nshal

      REAL								    &
	 TM1   (NLON,NLEV)						    &
	,QM1   (NLON,NLEV)						    &
	,TE    (NLON,NLEV)						    &
	,QE    (NLON,NLEV)						    &
	,APP1  (NLON,NLEV)						    &
	,APHP1 (NLON,NLEVP1)
!
      REAL								    &
	 SLMM  (NLON)							    &
	,TSM1M (NLON)							    &
	,RSFC  (NLON)							    &
	,SSFC  (NLON)							    &
	,APRC  (NLON)							    &
	,APRS  (NLON)							    &
	,ARPRC (NLON)
!
      INTEGER								    &
	 NTOPC (NLON)							    &
	,NBASEC(NLON)							    
!
      LOGICAL								    &
	 LPLAND(NLON)
!
      REAL								    &
	 CETA  (NLEV)							    &
	,CETAH (NLEVP1)
!
      integer, dimension(nlon,nlev) ::					    &
	   ilab
!
      real, dimension(nlon,nlev) ::					    &
	   ztp1								    &
	  ,zqp1								    &
	  ,zdpp1							    &
	  ,zdpkpk							    &
	  ,ztc								    &
	  ,zqc								    &
	  ,zdt								    &
	  ,zdq                                                              &
          ,dtd                                                              &
          ,dqd                                                              &
          ,dts                                                              &
          ,dqs                                                              &
          ,dztref0                                                          &
          ,dztref1                                                          &
          ,dztref                                                           &
          ,sztref
!
      integer, dimension(nlon) ::					    &
	   iqcd								    &
	  ,itest							    &
	  ,ilift							    &
	  ,ibase							    &
	  ,itop								    &
	  ,ishtop							    &
	  ,idx								   
!
      logical, dimension(nlon) ::					    &
	   lloc								    &
	  ,loswap                                                           &
	  ,llloc
!
      real, dimension(nlon) ::						    &
	   ztcb								    &
	  ,zsgtp							    &
	  ,zcape							    &
	  ,zsc
!
      real xlg, arg1, arg2
      XLG(ARG1,ARG2)=EXP(ARG2*ALOG(ARG1))
!
!*    DATA STATEMENTS.
!     ---- -----------
!
!     *NROWSUM*	  A LATITUDE-ROW COUNTER FOR THE CURRENT TIMESTEP.
!		  IT ENSURES REPRODUCIBLE DIAGNOSTICS IF MULTI-TASKING
!		  CAUSES ROWS TO BE PROCESSED IN AN ARBITRARY ORDER.
!
!!$      DATA NROWSUM / 0 /
!
!* local variables 
!  ----- ---------
!
      real ::        zbmts,    zbmtd,    zbmef,    zmincp,   zctops         &
          ,zctopss,  zgamma,   ztcrit,   zmu,      zetlab,   zepcor         &
          ,zepq,     zinicp,   zpwarn,   zdiagt,   zdiagw,   zctopsl        &
          ,zcons1,   zcons2,   zc1,      zc2,      zc3,      zc4            &
          ,zc5,      zc6,      zc7,      ztmst,    zpmax,    zqsatc         &
          ,zcor,     zqcd,     zzq,      zp1,      ztsp,     zzsp           &
          ,zratdp,   zapp1,    ztcld   

      integer ::     ietlab,   jk,       jl,       ihitop,   isum           &
          ,iit,      nldeep,   nsdeep,   iswap,    nlshal,   inpsec         &
          ,ihstop,   nsshal,   ihctop,   ilobas,   zsig1,    zsig2          

!*    PHYSICAL CONSTANTS.
!     -------- ----------
!
!     *ZBMTS*	  ADJUSTMENT TIMESCALE IN SECONDS FOR SHALLOW CONVECTION
!     *ZBMTD*	  ADJUSTMENT TIMESCALE IN SECONDS FOR DEEP CONVECTION.
!     *ZBMEF*	  DOWNDRAUGHT EFFICIENCY FOR DEEP CONVECTION.
!     *ZMINCP*	  MINIMUM CONVECTIVE AVAILABLE POTENTIAL ENERGY (CAPE)
!		  FOR WHICH AN ADJUSTMENT PROCESS WILL BE PERFORMED.
!     *ZCTOPSL/S* THRESHOLD VALUE OF CLOUD TOP SIGMA FOR SHALLOW
!		  CONVECTION FOR LAND/SEA.
!     *ZTCRIT*	  MARGIN OF STABILITY, WHICH DETERMINES WHETHER
!		  SLIGHTLY STABLE LAYERS MAY OR MAY NOT SUPRESS
!		  FURTHER ASCENT.
!     *ZGAMMA*	  CLOUD-TOP MIXING FRACTION USED FOR PARCEL ASCENT
!		  TEMPERATURE PROFILE (ZGAMMA=0.0 FOR NO MIXING).
!     *ZMU*	  CLOUD-TOP MIXING FRACTION USED FOR PARCEL ASCENT
!		  MOISTURE PROFILE (ZMU=0.0 FOR NO MIXING).
!     *ZETLAB*	  HYBRID COORD VALUE BELOW WHICH NEW PARCEL ASCENTS ARE
!		  ATTEMPTED (SHOULD BE IN MID/UPPER TROPOSPHERE).
!     *IETLAB*	  MODEL LEVEL UP TO WHICH NEW PARCEL ASCENTS ARE BEGUN.
!
      llloc=.false.
      ZBMTS=CBMTS
      ZBMTD=CBMTD
      ZBMEF=CBMEF
      ZMINCP=-10./RD ! original value
      ZCTOPSL=1.0 ! set to 1 to use only deep scheme
      ZCTOPSS=1.0 ! set to 1 to use only deep scheme
!!$      ZCTOPSL=0.725
!!$      ZCTOPSS=0.810
      ZTCRIT=0.0
!!$      ZTCRIT=1.0
!      ZGAMMA=0.2
      ZGAMMA=0.0
      ZMU=0.0 ! the default
      ZETLAB=0.3 ! the default
      ZETLAB=0.2
      IETLAB=1
      DO 11 JK=1,NLEV
      IF (CETA(JK).LT.ZETLAB) IETLAB=JK+1
   11 CONTINUE
!
!*    SECURITY PARAMETERS.
!     --------------------
!
!     *ZEPCOR*	  MINIMUM VALUE OF DENOMINATOR IN QSAT CALCULATION.
!     *ZEPQ*	  MINIMUM SPECIFIC HUMIDITY TO AVOID DIVERGENCE OF THE
!		  SATURATION POINT CALCULATIONS.
!     *ZINICP*	  PRESET FOR MAX. CAPE, TO CAPTURE NEGATIVE VALUES.
!     *ZPWARN*	  PRESSURE (IN PA) BELOW WHICH IMPOSED TOP OF
!		  CONVECTION IS TOO NEAR TROPOPAUSE.
!
      ZEPCOR=1.E-10
      ZEPQ=0.000002
      ZINICP=-1.E10
      ZPWARN=5000.
!
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------
!
!     *ZC1-ZC5*	  CONSTANTS FOR SATURATION POINT CALCULATIONS.
!     *ZC6*	  (RD/RV) USED FOR VAPOUR PRESSURE.
!     *ZC7*	  (1-RD/RV) USED FOR VAPOUR PRESSURE.
!     *INPSEC*	  PRINT FREQUENCY FOR CAPE DIAGNOSTICS.
!     *LOPRT*	  SWITCH FOR CAPE DIAGNOSTICS AT CURRENT STEP.
!
      ZTMST=TWODT
      ZDIAGT=0.5*TWODT
      ZDIAGW=ZDIAGT/RHOH2O
!
      ZCONS1=RD/CPD
      ZCONS2=ALV/CPD
!
      ZC1=1./ZCONS1
      ZC2=55.
      ZC3=2840.
      ZC4=3.5
      ZC5=0.2
      ZC6=0.622
      ZC7=0.378
!
      INPSEC=43200
!
!     ------------------------------------------------------------------
!
!*	   1.	  ALLOCATE SPACE AND POSITION VARIABLES.
!		  -------- ----- --- -------- ----------
!		  (removed in favor of automatic arrays)
!
  100 CONTINUE

!
!     ------------------------------------------------------------------
!
!*	   2.	  PRELIMINARY COMPUTATIONS.
!		  ----------- -------------
!
  200 CONTINUE
!
!
!*	   2.1	  T, Q, DELTA-P AT TIME LEVEL T+1.
!		  PRESET TEST FLAG AND TENDENCIES.
!
  210 CONTINUE
!
      DO 212 JK=1,NLEV
      DO 211 JL=1,NLON
      ZTP1(JL,JK)=TM1(JL,JK)+ZTMST*TE(JL,JK)
      ZQP1(JL,JK)=QM1(JL,JK)+ZTMST*QE(JL,JK)
      ZDPP1(JL,JK)=APHP1(JL,JK+1)-APHP1(JL,JK)
      ILAB(JL,JK)=0
      ZDT(JL,JK)=0.
      ZDQ(JL,JK)=0.
! look at tendencies for each scheme separately
      dts(jl,jk)=0.
      dqs(jl,jk)=0.
      dtd(jl,jk)=0.
      dqd(jl,jk)=0.
  211 CONTINUE
  212 CONTINUE
!
!	  2.2	 CHECK THAT ANY IMPOSED TOP LEVEL DOES NOT EXTEND
!		  DOWN TOO NEAR THE TROPOPAUSE.
!
  220 CONTINUE
!
      IF (NEADJTOP.GT.1) THEN
	 ZPMAX=0.
	 DO 221 JL=1,NLON
	 ZPMAX=MAX(ZPMAX,APP1(JL,NEADJTOP))
  221	 CONTINUE
	 IF (ZPMAX.GT.ZPWARN) THEN
	    write(*,*) ' *** WARNING : ATTEMPT TO SWITCH OFF BMADJ'	    &
		      ,' DOWN TO ',ZPMAX/100.,'HPA AT ROW ',KROW
	 ENDIF
      ENDIF
!
!     ----------------------------------------------------------------
!
!	  3.	 CLOUD ASCENT AND FLAGGING.
!		  ----- ------ --- ---------
!
  300 CONTINUE
!		  PRESET HIGHEST LEVEL IN ROW REACHED BY CONVECTION.
!
      IHITOP=NLEVM1
!
!*	  3.1	 INITIALIZE PARCEL PROPERTIES AND CLOUD INDICES.
!
  310 CONTINUE
!
      DO 311 JL=1,NLON
      ZTC(JL,NLEV)=ZTP1(JL,NLEV)
      ZQC(JL,NLEV)=ZQP1(JL,NLEV)
      IBASE(JL)=NLEV
      ITOP(JL)=NLEV
      ISHTOP(JL)=1
      ILIFT(JL)=1
  311 CONTINUE
!
!		  *ILAB*=1 INDICATES THAT A NEW PARCEL CAN BE STARTED.
!		  *ISHTOP* INDICATES THE IMPOSED SHALLOW CONVECTIVE TOP.
!
      DO 313 JK=IETLAB,NLEV
      DO 312 JL=1,NLON
      ILAB(JL,JK)=1
      IF (APP1(JL,JK).LT.(ZCTOPSL*APHP1(JL,NLEVP1))) ISHTOP(JL)=JK+1
  312 CONTINUE
  313 CONTINUE
!
!	  3.2	 START LOOP OVER LEVELS FOR PARCEL ASCENTS.
!
  320 CONTINUE
!
      DO 341 JK=NLEVM1,1,-1
!
!		  EXTRACT PARCEL-BASE TEMPERATURE FOR VECTORISATION.
!
      DO 321 JL=1,NLON
      ZTCB(JL)=ZTC(JL,IBASE(JL))
      ZDPKPK(JL,JK)=XLG((APP1(JL,JK)/APP1(JL,JK+1)),ZCONS1)
  321 CONTINUE
!
!	  3.3	 PARCEL CALCULATIONS AT CURRENT LEVEL, SET FLAGS.
!
  330 CONTINUE
!
      ISUM=0
      DO 331 JL=1,NLON
!
!		  LIFT PARCEL UP A DRY ADIABAT.
!		  *ZMU* IS A MOISTURE MIXING FRACTION FROM CLOUD-TOP
!		  ENTRAINMENT.
!
      ZTC(JL,JK)=ZTC(JL,JK+1)*ZDPKPK(JL,JK)
      ZQC(JL,JK)=ZQC(JL,JK+1)*(1.-ZMU)+ZQP1(JL,JK)*ZMU
!
!		  SUPERSATURATION AND CORRECTION FOR MOIST ADIABAT.
!
      ZQSATC=C2ES*EXP(C3LES*(ZTC(JL,JK)-TMELT)*				    &
			(1./(ZTC(JL,JK)-C4LES)))/APP1(JL,JK)
      ZCOR=1./MAX(ZEPCOR,ABS(1.-VTMPC1*ZQSATC))
      ZQSATC=ZQSATC*ZCOR
      ZQCD=MAX(0.,(ZQC(JL,JK)-ZQSATC)/(1.+C5LES*ZCONS2*ZQSATC*ZCOR*	    &
				      (1./(ZTC(JL,JK)-C4LES))**2))
      IF (ZQCD.EQ.0.) THEN
	 IQCD(JL)=0
      ELSE
	 IQCD(JL)=1
	 ZTC(JL,JK)=ZTC(JL,JK)+ZCONS2*ZQCD
	 ZQC(JL,JK)=ZQC(JL,JK)-ZQCD
!
!		  SECOND ITERATION FOR MOIST ADIABAT.
!
	 ZQSATC=C2ES*EXP(C3LES*(ZTC(JL,JK)-TMELT)*			    &
			  (1./(ZTC(JL,JK)-C4LES)))/APP1(JL,JK)
	 ZCOR=1./MAX(ZEPCOR,(1.-VTMPC1*ZQSATC))
	 ZQSATC=ZQSATC*ZCOR 
	 ZQCD=(ZQC(JL,JK)-ZQSATC)/(1.+C5LES*ZCONS2*ZQSATC*ZCOR*		    &
				  (1./(ZTC(JL,JK)-C4LES))**2)
	 ZTC(JL,JK)=ZTC(JL,JK)+ZCONS2*ZQCD
	 ZQC(JL,JK)=ZQC(JL,JK)-ZQCD
      ENDIF
!
!		  PARCEL TEMPERATURE FOR BUOYANCY.
!		  *ZGAMMA* IS A MIXING FRACTION FROM CLOUD-TOP
!		  ENTRAINMENT.
!
      ZZQ=MAX(ZEPQ,ZQP1(JL,JK))
      ZP1=ZZQ*APP1(JL,JK)/(ZC6+ZC7*ZZQ)
      ZTSP=ZC2+ZC3/(ZC4*ALOG(ZTP1(JL,JK))-ALOG(ZP1)-ZC5)
      ZZSP=APP1(JL,JK)*XLG((ZTSP/ZTP1(JL,JK)),ZC1)
      ZRATDP=(APP1(JL,NLEV)-ZZSP)/(APP1(JL,NLEV)-APP1(JL,JK))
      ZAPP1=APP1(JL,IBASE(JL))
      ZTCLD=ZTC(JL,JK)*(1.0-ZGAMMA*ZRATDP)+ZGAMMA*(ZTP1(JL,JK)+		    &
	    ZTCB(JL)*XLG((APP1(JL,JK)/ZAPP1),ZCONS1)*			    &
	    (ZRATDP-1.0))
!
!		  BUOYANCY CONSIDERATION - SET FLAGS.
!		  *LOVS* ACCOUNTS FOR MARGINAL STABILITY.
!		  *ILAB* = 2   DRY: UNSTABLE OR MARGINALLY STABLE.
!		  *ILAB* = 3   MOIST: UNSTABLE, OR MARGINALLY STABLE
!			       NEAR CLOUD-BASE.
!
      LOVS=(ZTP1(JL,JK)-ZTCLD).GT.ZTCRIT
      LOIS=(ZTP1(JL,JK)-ZTCLD).LT.0.
      LLOC(JL)=IQCD(JL).NE.0 ! always false for dry runs (subsaturated)
      LO2=.NOT.(LLOC(JL).OR.LOVS)
      LO3=LLOC(JL).AND.(LOIS.OR.(.NOT.(LOVS.OR.(JK.LT.(IBASE(JL)-2)))))
!
      IF (LO2) THEN    
	 ILAB(JL,JK)=2
         if(ilab(jl,jk+1).eq.1.and.itop(jl).lt.nlev) &
              llloc(jl) = .true.
	 itop(jl)=jk ! treat dry unstable parcels like moist unstable parcels:
      ENDIF

      IF (LO3) THEN
	 ILAB(JL,JK)=3
	 ITOP(JL)=JK
      ELSE
	 IBASE(JL)=IBASE(JL)-1
      ENDIF

!
!		  IF LIFT FROM BELOW IS NOT BUOYANT, RESET PARCEL TO
!		  ENVIRONMENT FOR POSSIBLE NEW ASCENT.
!
! if lift from below is not buoyant AND it's not above a lower cloud
      IF (ILAB(JL,JK).EQ.1.AND..not.(ITOP(JL).LT.NLEV)) THEN
	 ZTC(JL,JK)=ZTP1(JL,JK)
	 ZQC(JL,JK)=ZQP1(JL,JK)
      ENDIF
! original code:
!!$      IF (ILAB(JL,JK).EQ.1) THEN
!!$	 ZTC(JL,JK)=ZTP1(JL,JK)
!!$	 ZQC(JL,JK)=ZQP1(JL,JK)
!!$      ENDIF
!
!		  AVOID ATTEMPT FOR NEW LIFT IF ABOVE MIDTROPOSPHERE
!		  OR ABOVE LOWER CLOUD.
!
      LO= (ILAB(JL,JK).EQ.0)
! original code:
!!$      LO= (ILAB(JL,JK).EQ.0).OR.					    & 
!!$	 ((ILAB(JL,JK).EQ.1).AND.(ITOP(JL).LT.NLEV))
      IF (LO) THEN
	 ZTC(JL,JK)=150.
	 ZQC(JL,JK)=0.
      ENDIF
!
      ISUM=ISUM+ILAB(JL,JK)
  331 CONTINUE
!
!		  END ASCENT IF NO CONVECTING OR NEW START POINTS.
!		  *IHITOP* FLAGS DEEPEST ACTUAL MOIST CONVECTION
!			   OR POSSIBLE LIFT POINT.
!
      IF (ISUM.EQ.0) GO TO 342 
      IHITOP=JK
!
!*	   3.4	  END OF VERTICAL LOOP.
!		  CHECK THAT ALL ASCENTS ARE COMPLETED BY TOP LEVEL.
!		  RESTRICT HIGHEST PARCEL TOP IF ABOVE IMPOSED TOP.
!
  340 CONTINUE
!
  341 CONTINUE
!
  342 CONTINUE
!
      IF (IHITOP.LT.NEADJTOP) THEN
	 write(*,*) ' *** WARNING : CONVECTION EXTENDS ABOVE IMPOSED TOP'   &
	       ,' TO LEVEL ',IHITOP,' AT ROW ',KROW
	 write(*,*) '		  : BUT IS LIMITED TO IMPOSED TOP LEVEL '   &
	       ,NEADJTOP
	 DO 343 JL=1,NLON
	 ITOP(JL)=MAX(ITOP(JL),NEADJTOP)
  343	 CONTINUE
      ELSE IF (IHITOP.EQ.1) THEN
	 write(*,*) ' *** WARNING : CONVECTION REACHES TOP MODEL'	    &
	       ,' LEVEL ',JK,' AT ROW ',KROW
      ENDIF
!
!*	   3.5	  SET INDECES FLAGGING CLOUD LAYERS AND DIAGNOSE CAPES.
!		  *ITOP* ALREADY INDICATES THE TOP IN-CLOUD LEVEL.
!		  *IHITOP* IS RESET TO THE HIGHEST CLOUD TOP IN ROW.
!		  *IBASE* INDICATES THE TOP OF THE SUB-CLOUD LAYER.
!		  *ILIFT* INDICATES THE START OF THE PARCEL ASCENT.
!
  350 CONTINUE
!
      IHITOP=NLEV
      DO 351 JL=1,NLON
      ZCAPE(JL)=0.
      IBASE(JL)=1
      IIT=ITOP(JL)
      ZSGTP(JL)=APP1(JL,IIT)/APHP1(JL,NLEVP1)
      IHITOP=MIN(IHITOP,IIT)
  351 CONTINUE
!
      DO 354 JK=IHITOP,NLEV
      DO 352 JL=1,NLON
!!$      LO3=ILAB(JL,JK).EQ.3 ! original code
      LO3=(ILAB(JL,JK).EQ.3).or.(ilab(jl,jk).eq.2) ! calculate CAPE for dry
      IF (LO3) IBASE(JL)=JK+1                      ! unstable regions too
      IF (LO3.OR.(JK.EQ.IBASE(JL))) THEN	      
	 ZCAPE(JL)=ZCAPE(JL)+						    &
		  (ZTC(JL,JK)-ZTP1(JL,JK))*ZDPP1(JL,JK)/APP1(JL,JK)
      ENDIF
!      if(llloc(jl)) write(*,*) jl,zcape(jl),zcape(jl).lt.zmincp
      
  352 CONTINUE
      IF (JK.GE.MAX(2,IETLAB)) THEN
	 DO 353 JL=1,NLON
	 LO=(ILAB(JL,JK-1).GT.1).AND.(ILAB(JL,JK).EQ.1)			    &
	    .AND.(ILIFT(JL).EQ.1)
	 IF (LO) ILIFT(JL)=JK
  353	 CONTINUE
      ENDIF
  354 CONTINUE
!
!		  AVOID COUNTING SHALLOW LAYERS IN BOUNDARY LAYER.
!		  ALLOW FOR SMALL NEGATIVE CAPE FROM OVERSHOOTS.
!
      DO 355 JL=1,NLON
! be mindful of this change:
      IF (ZCAPE(JL).LT.ZMINCP) ZCAPE(JL)=0.
!!$      IF (ZCAPE(JL).LT.ZMINCP.OR.ITOP(JL).GE.NLEVM1) ZCAPE(JL)=0. ! orig.
  355 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   4.	  DEEP ADJUSTMENT PROCESS.
!		  ---- ---------- --------
!
  400 CONTINUE
!
!*	   4.1	  SAMPLE DEEP CONVECTIVE POINTS (LLOC, ITEST=1).
!
  410 CONTINUE
!
!
      DO 411 JL=1,NLON
      IF (LPLAND(JL)) THEN
	 ZCTOPS=ZCTOPSL
      ELSE
	 ZCTOPS=ZCTOPSS
      ENDIF
! be mindful of this change:
      LLOC(JL)=.NOT.( (ZSGTP(JL).GE.ZCTOPS).OR.(ZCAPE(JL).EQ.0.))
! original code:
!!$      LLOC(JL)=.NOT.( (ZSGTP(JL).GE.ZCTOPS).OR.(ZCAPE(JL).EQ.0.)             &
!!$		     .OR.((ILIFT(JL)-1).EQ.ITOP(JL)) )
      IF (LLOC(JL)) THEN
	 ITEST(JL)=1
      ELSE
	 ITEST(JL)=0
      ENDIF
  411 CONTINUE

!
!		  FIND NUMBER OF DEEP CONVECTIVE POINTS AND
!		  FORM ARRAY OF INDECES (IDX) FOR GATHERS IN BMDEEP.
!		  PRESET OTHER COUNTERS.
!
      CALL WHENEQ(NLON,ITEST,1,1,IDX,NDEEP)
      NLDEEP=0
      NSDEEP=0
      ISWAP=0
!
      IF (NDEEP.EQ.0) GOTO 500
!
!*	   4.2	  PREFORM DEEP CONVECTIVE ADJUSTMENT.
!
  420 CONTINUE
! be mindful of this change -- relaxes the whole column
!!$      itop = ietlab
!
!!$      write(*,*) 'got a deep convecting point!'
      CALL BMDEEP(							    &
	  NDEEP,    ZBMTD,    ZBMEF,	IDX				    &
	 ,NLON,	    NLEV,     NLEVM1,	CETAH				    &
	 ,ALS,	    ALV,      CPD,	G,	  RD,	    TMELT	    &
	 ,VTMPC1,   VTMPC2						    &
	 ,C2ES,	    C3IES,    C3LES,	C4IES,	  C4LES,    C5LES	    &
	 ,APP1,	    TSM1M,    RSFC,	SSFC				    &
	 ,ZDPP1,    ZDPKPK,   ZTP1,	ZQP1				    &
	 ,ZTC,	    ZQC,      ZDT,	ZDQ				    &
	 ,IBASE,    ITOP,     ILIFT,    dztref0,  dztref1,  dztref        &
	 )
!
!*	   4.3	  RESET NEGATIVE PRECIPITATION AND CORRESPONDING
!		  TENDENCIES.  FLAG SUCH POINTS AS "SWAPPED" FOR
!		  POSSIBLE SHALLOW CONVECTION.
!
  430 CONTINUE
!
      DO 431 JL=1,NLON
      LOSWAP(JL)=(SSFC(JL)+RSFC(JL)).LT.0.
      IF (LOSWAP(JL).and.int(alv).ne.0) THEN !points aren't swapped in dry runs
	 ISWAP=ISWAP+1
	 LLOC(JL)=.FALSE.
	 ITEST(JL)=-1
	 SSFC(JL)=0.
	 RSFC(JL)=0.
      ELSE
	 NLDEEP=NLDEEP+ITEST(JL)*NINT(SLMM(JL))
      ENDIF
  431 CONTINUE
      NDEEP=NDEEP-ISWAP
      NSDEEP=NDEEP-NLDEEP
!
      DO 433 JK=IHITOP,NLEV
      DO 432 JL=1,NLON
      IF (LOSWAP(JL).and.int(alv).ne.0) THEN !points aren't swapped in dry runs
	 ZDT(JL,JK)=0.
	 ZDQ(JL,JK)=0.
      ENDIF
  432 CONTINUE
  433 CONTINUE
!
!*	   4.4	  STORE TENDENCIES FOR DEEP CONVECTION POINTS ONLY.
!
  440 CONTINUE
!
! be mindful of this change
!!$      DO 442 JK=ietlab,NLEV
      DO 442 JK=IHITOP,NLEV ! original code
      DO 441 JL=1,NLON
      IF (LLOC(JL)) THEN
         dtd(jl,jk)=zdt(jl,jk)
         dqd(jl,jk)=zdq(jl,jk)
	 TE(JL,JK)=TE(JL,JK)+ZDT(JL,JK)
	 QE(JL,JK)=QE(JL,JK)+ZDQ(JL,JK)
      ENDIF
  441 CONTINUE
  442 CONTINUE
!
!*	   4.5	  DEEP CONVECTION PARAMETERS FOR CLOUD SCHEME.
!                              (removed by CW)
  450 CONTINUE
!
!
!*	   4.6	  ZONAL MEAN DIAGNOSTICS FOR DEEP CONVECTION.
!                              (removed by CW)
  460 CONTINUE
!
!
!     ------------------------------------------------------------------
!
!*	   5.	  SHALLOW CONVECTION PROCESS.
!		  ------- ---------- --------
!
  500 CONTINUE
!
!*	   5.1	  SAMPLE SHALLOW CONVECTIVE POINTS (LLOC, ITEST=2).
!		  SWAPPED POINTS HAVE TOP SET TO IMPOSED MAXIMUM.
!		  PRESET CONDENSATION RATE.
!
  510 CONTINUE
!
      IHSTOP=NLEV
      DO 511 JL=1,NLON
      ZSC(JL)=0.
      LLOC(JL)=.NOT.( (ZCAPE(JL).EQ.0.).OR.((RSFC(JL)+SSFC(JL)).GT.0.)	    &
		     .OR.(IBASE(JL).LE.ISHTOP(JL)) )
      IF (LLOC(JL)) THEN
	 ITEST(JL)=2
	 ITOP(JL)=MAX(ITOP(JL),ISHTOP(JL))
	 IHSTOP=MIN(IHSTOP,ITOP(JL))
      ENDIF
  511 CONTINUE
!		  CHECK HIGHEST LEVEL OF SHALLOW CONVECTIVE TENDENCIES
!		  (WHICH INCLUDES A LEVEL ABOVE CLOUD TOP) AGAINST THE
!		  IMPOSED TOP.	IN THE UNLIKELY EVENT THAT THE IMPOSED
!		  TOP IS PASSED, RESTRICT TENDENCIES TO IMPOSED TOP.
!
      IF (IHSTOP.LT.NEADJTOP+1) THEN
	 write(*,*) ' *** WARNING : SHALLOW CONVECTIVE TENDENCIES EXTEND'   &
		,' ABOVE IMPOSED TOP TO LEVEL ',IHSTOP-1		    &
		,' AT STEP ROW ',KROW
	 write(*,*) '		  : BUT ARE LIMITED TO IMPOSED TOP LEVEL'   &
		,NEADJTOP
	 DO 512 JL=1,NLON
	 IF (LLOC(JL)) ITOP(JL)=MAX(ITOP(JL),NEADJTOP+1)
  512	 CONTINUE
      ENDIF
!
!		  FIND NUMBER OF SHALLOW CONVECTIVE POINTS.
!		  FORM ARRAY OF INDECES (IDX) FOR GATHERS IN BMSHAL.
!		  PRESET OTHER COUNTERS.
!
      CALL WHENEQ(NLON,ITEST,1,2,IDX,NSHAL)
      NLSHAL=0
      NSSHAL=0
!
      IF (NSHAL.EQ.0) GOTO 600
!
!*	   5.2	  PERFORM SHALLOW CONVECTIVE ADJUSTMENT.
!
  520 CONTINUE
!
      write(*,*) 'shallow convecting point'
      CALL BMSHAL(							    &
	   NSHAL,    ZBMTS,    IDX					    &
	  ,NLON,     NLEV,     NLEVM1,	 NLEVP1				    &
	  ,CPD,	     G,	       RD,	 TMELT				    &
	  ,C2ES,     C3LES,    C4LES					    &
	  ,APP1,     APHP1						    &
	  ,ZTP1,     ZQP1						    &
	  ,ZDT,	     ZDQ,      ZSC					    &
	  ,IBASE,    ITOP						    &
	  ,IHCTOP,   ILOBAS,   sztref                                       &
	  ) 
!
!*	   5.3	  DIAGNOSE NUMBERS OF CONVECTING POINTS IN ROW.
!
  530 CONTINUE
!
      DO 531 JL=1,NLON
      IF (LLOC(JL)) NLSHAL=NLSHAL+NINT(SLMM(JL))
  531 CONTINUE
      NSSHAL=NSHAL-NLSHAL
!
!*	   5.4	  STORE TENDENCIES FOR SHALLOW CONVECTION POINTS ONLY.
!
  540 CONTINUE
!
      DO 542 JK=IHCTOP,ILOBAS
      DO 541 JL=1,NLON
      IF (LLOC(JL)) THEN
         dts(jl,jk)=zdt(jl,jk)
         dqs(jl,jk)=zdq(jl,jk)
	 TE(JL,JK)=TE(JL,JK)+ZDT(JL,JK)
	 QE(JL,JK)=QE(JL,JK)+ZDQ(JL,JK)
      ENDIF
  541 CONTINUE
  542 CONTINUE
!
!*	   5.5	  SHALLOW CONVECTION PARAMETERS FOR CLOUD SCHEME.
!                               (removed by CW)
  550 CONTINUE
!
!*	   5.6	  ZONAL MEAN DIAGNOSTICS FOR SHALLOW CONVECTION.
!                               (removed by CW)

  560 CONTINUE
!
!
!     ------------------------------------------------------------------
!
!*	   6.	  STORE PRECIPITATION AND COMPUTE DIAGNOSTICS.
!		  ----- ------------- --- ------- ------------
!
  600 CONTINUE
!
      ZSIG1=0.
      ZSIG2=0.
      IF (NDEEP.GT.0) THEN
	 DO 601 JL=1,NLON
	 APRC(JL)=APRC(JL)+ZDIAGW*(RSFC(JL)+SSFC(JL))
	 APRS(JL)=APRS(JL)+ZDIAGW*SSFC(JL)
	 ZSIG1=ZSIG1+RSFC(JL)
	 ZSIG2=ZSIG2+SSFC(JL)
  601	 CONTINUE
      ENDIF
!
!		  GLOBAL DIAGNOSTICS.
!		  *DCVGR*  IS THE RAINFALL GENERATION.
!		  *DCVGS*  IS THE SNOWFALL GENERATION.
!		  *DCVMOI* IS THE ENVIRONMENTAL MOISTENING, JUST THE
!			   NEGATIVE OF THE PRECIPITATION.
!		  EVAPORATION / MELTING RATES ARE NOT RETURNED BY
!		  *BMDEEP*, SO THE GENERATION TERMS ARE OBTAINED FROM
!		  THE SURFACE PRECIPITATION TERMS.
!
!
!!$      DCVGR =DCVGR +ZDIAGW*PBUDW*ZSIG1
!!$      DCVGS =DCVGS +ZDIAGW*PBUDW*ZSIG2
!!$      DCVMOI=DCVMOI+ZDIAGW*PBUDW*(-ZSIG1-ZSIG2)
!
!
!		  MASK DIAGNOSTICS (TOTAL CONVECTIVE TERMS).
!			       (removed by cw)
!
!*	   6.1	  COMPUTE CAPE STATISTICS FOR CURRENT LATITUDE ROW.
!
  610 CONTINUE
!
!!$      IF (LOPRT) THEN
!!$!
!!$	 ZMCAPE=0.
!!$	 ZMLCAPE=0.
!!$	 ZMFCAPE=0.
!!$	 INTOT=0
!!$	 INLAND=0
!!$	 INFAIL=0
!!$	 ZMXCAPE=ZINICP
!!$	 IMXCAPE=-1
!!$!
!!$	 DO 611 JL=1,NLON
!!$	 IF (ZCAPE(JL).NE.0.) THEN
!!$	    ZMCAPE=ZMCAPE+ZCAPE(JL)
!!$	    INTOT=INTOT+1
!!$	    ZMLCAPE=ZMLCAPE+ZCAPE(JL)*SLMM(JL)
!!$	    INLAND=INLAND+NINT(SLMM(JL))
!!$	    IF (ITEST(JL).LE.0) THEN
!!$	       ZMFCAPE=ZMFCAPE+ZCAPE(JL)
!!$	       INFAIL=INFAIL+1
!!$	    ENDIF
!!$	 ENDIF
!!$	 IF (ZCAPE(JL).GT.ZMXCAPE) THEN
!!$	    ZMXCAPE=ZCAPE(JL)
!!$	    IMXCAPE=JL
!!$	 ENDIF
!!$  611	 CONTINUE
!!$	 ZMSCAPE=ZMCAPE-ZMLCAPE
!!$	 INSEA=INTOT-INLAND
!!$!
!!$	 IF (INTOT.GT.0) THEN
!!$	    ZNCAPE=ZMCAPE*RD/FLOAT(INTOT)
!!$	 ELSE
!!$	    ZNCAPE=0.
!!$	 ENDIF
!!$	 IF (INLAND.GT.0) THEN
!!$	    ZNLCAPE=ZMLCAPE*RD/FLOAT(INLAND)
!!$	 ELSE
!!$	    ZNLCAPE=0.
!!$	 ENDIF
!!$	 IF (INSEA.GT.0) THEN
!!$	    ZNSCAPE=ZMSCAPE*RD/FLOAT(INSEA)
!!$	 ELSE
!!$	    ZNSCAPE=0.
!!$	 ENDIF
!!$	 IF (INFAIL.GT.0) THEN
!!$	    ZNFCAPE=ZMFCAPE*RD/FLOAT(INFAIL)
!!$	 ELSE
!!$	    ZNFCAPE=0.
!!$	 ENDIF
!!$	 ZNXCAPE=ZMXCAPE*RD
!!$	 ZLAT=180.*ASIN(0.5*PTWOMU)/API
!!$!
!!$	 NCONV=NSHAL+NDEEP
!!$	 NLCONV=NLDEEP+NLSHAL
!!$	 NSCONV=NSDEEP+NSSHAL
!
!*	   6.2	  COMPUTE GLOBAL CONVECTIVE STATISTICS.
!		  USE LOCKS AND INDEPENDENT LATITUDE COUNTER
!		  FOR REPRODUCIBLE RESULTS WHEN MULTI-TASKING.
!
  620	 CONTINUE
!
!
!!$	    NROWSUM=NROWSUM+1
!!$
!!$	    IF (NROWSUM.EQ.1) THEN
!!$	       PRINT 9001,NSTEP
!!$	       GCAPE=0.
!!$	       NGTOT=0
!!$	       NGCONV=0
!!$	       NGDEEP=0
!!$	       NGSHAL=0
!!$	       NGSWAP=0
!!$	       NGFAIL=0
!!$	    ENDIF
!!$
!!$	    NGTOT =NGTOT +INTOT
!!$	    NGCONV=NGCONV+NCONV
!!$	    NGDEEP=NGDEEP+NDEEP
!!$	    NGSHAL=NGSHAL+NSHAL
!!$	    NGSWAP=NGSWAP+ISWAP
!!$	    NGFAIL=NGFAIL+INFAIL
!!$	    GCAPE=GCAPE+ZMCAPE*RD
!!$!
!!$	    write(*,9002) ZLAT,ZNCAPE,ZNLCAPE,ZNSCAPE,ZNXCAPE,IMXCAPE	      &
!!$		     ,NCONV,NLCONV,NSCONV,NDEEP,NSHAL,ISWAP		      &
!!$		     ,INFAIL,ZNFCAPE
!!$!
!!$	    IF (NROWSUM.EQ.NGL) THEN
!!$	       write(*,9003)
!!$	       write(*,9004) NGCONV,NGDEEP,NGSHAL,NGSWAP,NGFAIL,GCAPE/NGTOT
!!$	       NROWSUM=0
!!$	    ENDIF
!!$!
!!$!
!!$      ENDIF
!     ------------------------------------------------------------------
!
!*	   7.	  RETURN WORKSPACE.
!		  ------ ----------
!                 
  700 CONTINUE
!
!
!     ------------------------------------------------------------------
!
!*	   8.	  FORMAT STATEMENTS.
!		  ------ -----------
!
  800 CONTINUE
!
 9001 FORMAT (/,1X,'CONVECTIVE ACTIVITY STATISTICS AT STEP ',I5,/,	    &
	4X,'LATITUDE	    ZONAL MEAN CAPES',12X,'MAX CAPE',		    &
	4X,'NUMBER OF CONVECTIVE POINTS',/,				    &
       60X,'-------- SUCCESSFUL --------	---- FAIL -----',//,	    &
       16X,'TOTAL     LAND	 SEA	  VALUE	 INDEX',		    &
	2X,'TOTAL LAND	SEA   DEEP SHALL  SWAP	TOTAL MEAN CAPE',/)
 9002 FORMAT (1X,5F10.3,8I6,F10.3)
 9003 FORMAT (/,8X,17('*'),' GLOBAL CONVECTION STATISTICS ',17('*'),/,	    &
	8X,' NUMBER OF CONVECTIVE POINTS      FAIL',			    &
       10X,'GLOBAL MEAN CAPE',/,					    &
	8X,'TOTAL    DEEP  SHALLOW   SWAP    POINTS',/)
 9004 FORMAT (8X,5(I5,3X),8X,F12.3,/)
      
      RETURN

    END SUBROUTINE BMADJ

!*/
!*/ ************************************************************** BMDEEP
!*/
   SUBROUTINE BMDEEP(							    &
	   NDEEP,    CBMTD,    CBMEF,	 KDX				    &
	  ,NLON,     NLEV,     NLEVM1,	 CETAH				    &
	  ,ALS,	     ALV,      CPD,	 G,	   RD,	     TMELT	    &
	  ,VTMPC1,   VTMPC2						    &
	  ,C2ES,     C3IES,    C3LES,	 C4IES,	   C4LES,    C5LES	    &
	  ,APP1,     TSM1M,    RSFC,	 SSFC				    &
	  ,PDPP1,    PDPKPK,   PTP1,	 PQP1				    &
	  ,PTC,	     PQC,      PDT,	 PDQ				    &
	  ,KBASE,    KTOP,     KLIFT,    dztref0,  dztref1,  dztref         &
	  )

     integer, intent(in) ::    ndeep,    nlon,     nlev,     nlevm1
     real, intent(in)    ::    cbmtd,    cbmef,    als,	     alv            &
          ,cpd,	     g,        rd,       tmelt,    vtmpc1,   vtmpc2         &
          ,c2es,     c3ies,    c3les,    c4ies,    c4les,    c5les 


!**** *BMDEEP* - PERFORMS SATURATION POINT ADJUSTMENT
!		 FOR DEEP CONVECTION.
!
!     ORIGINAL VERSION		       B.RITTER	 E.C.M.W.F.  31/10/84.
!     MODIFIED			     M.J.MILLER	 E.C.M.W.F.  22/07/91.
!     IMPLEMENTED & TUNED IN UGCM1   J.M.SLINGO	 U.G.A.M.P.	12/91.
!     IMPLEMENTED IN UGCM2	    M.BLACKBURN	 U.G.A.M.P.  13/05/94.
!
!     PURPOSE.
!     --------
!
!	   THIS ROUTINE ADJUSTS THE TEMPERATURE AND SPECIFIC HUMIDITY
!     PROFILES OF DEEP CONVECTIVE POINTS TOWARDS AN INTERNALLY SPECIFIED
!     REFERENCE PROFILE.
!
!**   INTERFACE.
!     ----------
!
!	   *BMDEEP* IS CALLED FROM *BMADJ*.
!	   THE ROUTINE IS ARGUMENT DRIVEN, TAKING ITS INPUT ENTIRELY
!     FROM THE VARIABLES AND ARRAYS SUPPLIED IN THE ARGUMENT LIST.
!     FOR VECTORISATION, THE INPUT IS GATHERED IN TEMPORARY ARRAYS OF
!     LENGTH *NDEEP* (THE NUMBER OF GRIDPOINTS WITH DEEP CONVECTION
!     ON A LATITUDE ROW) USING INDEX ARRAY *KDX*.  THE ROUTINE RETURNS
!     ITS OUTPUT AS ARGUMENTS, HAVING SCATTERED TO THE FULL GRID.
!
!  ---------------------------------------------------------------------
!
!	   INPUT ARGUMENTS:
!
!     *NDEEP*	NUMBER OF GRIDPOINTS WITH DEEP CONVECTION IN ROW.
!     *CBMTD*	ADJUSTMENT TIMESCALE IN SECONDS FOR DEEP CONVECTION.
!     *CBMEF*	DOWNDRAUGHT EFFICIENCY FOR DEEP CONVECTION.
!     *KDX*	INDEX OF CONVECT. POINTS ON FULL GRID, DIMENSION (NLON).
!     *NLON*	NUMBER OF LONGITUDES.
!     *NLEV*	NUMBER OF LEVELS.
!     *NLEVM1*	(NLEV-1).
!     *CETAH*	HYBRID COORDINATE AT HALF LEVELS,    DIMENSION (NLEV+1).
!     *ALS*	LATENT HEAT FOR SUBLIMATION.
!     *ALV*	LATENT HEAT FOR VAPORISATION.
!     *CPD*	SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR.
!     *G*	GRAVITATIONAL ACCELERATION.
!     *RD*	GAS CONSTANT FOR DRY AIR.
!     *TMELT*	TEMPERATURE OF FUSION OF ICE.
!     *VTMPC1*	CONSTANT FOR VIRTUAL EFFECTS, (RV/RD-1).
!     *VTMPC2*	CONSTANT FOR VIRTUAL EFFECTS, (CPV/CPD-1).
!     *C__ES*	CONSTANTS USED FOR COMPUTATION OF SATURATION SPECIFIC
!		HUMIDITY OVER LIQUID WATER (*C_LES*) OR ICE (*C_IES*).
!     *C2ES*	(RD/RV)*(SVP AT REFERENCE TEMPERATURE C4_ES).
!     *C3_ES*	CONSTANT FOR SVP.
!     *C4_ES*	REFERENCE TEMPERATURE FOR SVP.
!     *C5_ES*	(C3_ES*(TMELT-C4_ES)).
!     *APP1*	FULL LEVEL PRESSURE,		  DIMENSION (NLON,NLEV).
!     *TSM1M*	(T-1) SURFACE TEMPERATURE,	       DIMENSION (NLON).
!     *PDPP1*	LAYER THICKNESS (DELTA-P),	  DIMENSION (NLON,NLEV).
!     *PDPKPK*	(P(K)/P(K+1))**(RD/CPD),	  DIMENSION (NLON,NLEV).
!     *PTP1*	PRELIMINARY (T+1) TEMPERATURE,	  DIMENSION (NLON,NLEV).
!     *PQP1*	PRELIMINARY (T+1) SPECIFIC HUMID, DIMENSION (NLON,NLEV).
!     *PTC*	PARCEL ASCENT TEMPERATURE,	  DIMENSION (NLON,NLEV).
!     *PQC*	PARCEL ASCENT SPECIFIC HUMIDITY,  DIMENSION (NLON,NLEV).
!     *KBASE*	HIGHEST LEVEL IN SUB-CLOUD LAYER,      DIMENSION (NLON).
!     *KTOP*	HIGHEST LEVEL IN CLOUD LAYER,	       DIMENSION (NLON).
!     *KLIFT*	LEVEL OF PARCEL ASCENT ORIGIN,	       DIMENSION (NLON).
!
!	   OUTPUT ARGUMENTS:
!
!     *RSFC*	SURFACE RAINFALL RATE,		       DIMENSION (NLON).
!     *SSFC*	SURFACE SNOWFALL RATE,		       DIMENSION (NLON).
!     *PDT*	TEMPERATURE TENDENCY,		  DIMENSION (NLON,NLEV).
!     *PDQ*	SPECIFIC HUMIDITY TENDENCY,	  DIMENSION (NLON,NLEV).
!
!  ---------------------------------------------------------------------
!*ENDIF
!	   OUTPUT ARRAYS ARE ONLY ASSIGNED AT DEEP CONVECTING POINTS.
!     ALL ARRAY ELEMENTS SHOULD BE PRESET TO ZERO ON INPUT.
!
!     METHOD.
!     -------
!
!	   THE BASIC REFERENCE PROFILE IS CONSTRUCTED, BASED ON THE
!     MEAN OF THE ENVIRONMENT AND PARCEL TEMPERATURES AT THE FIRST
!     LEVEL ABOVE THE PLANETARY BOUNDARY LAYER.	 THE PBL IS ASSUMED
!     TO CONSIST OF THE LOWEST *NPBL* LEVELS, SPECIFIED IN THE CODE.
!     THE SLOPE OF THE REFERENCE TEMPERATURE PROFILE IS A PRESCRIBED
!     FRACTION OF THE MOIST ADIABAT (TO APPROXIMATE A VIRTUAL MOIST
!     ADIABAT) UP TO THE FREEZING LEVEL.  ABOVE THIS IT RETURNS, AT
!     CLOUD TOP, TO THE MOIST ADIABAT FROM CLOUD BASE.	THE REFERENCE
!     MOISTURE PROFILE IS LINEARLY INTERPOLATED BETWEEN PRESCRIBED
!     SUBSATURATION VALUES AT CLOUD BASE, FREEZING LEVEL AND CLOUD TOP.
!	   THE BASIC REFERENCE PROFILE IS MODIFIED TO GUARANTEE ENERGY
!     CONSERVATION THROUGH THE ADJUSTMENT PROCESS.  FINALLY, HEATING
!     RATES AND MOISTURE CHANGES ARE CALCULATED AS A RELAXATION TOWARDS
!     THE REFERENCE PROFILE.
!	   MODIFICATION OF THE PBL BY DOWNDRAUGHTS HAS BEEN INCLUDED.
!     ******************************************************************
!     * N.B.  THE RELAXATION TIME *ZTAU* AND DOWNDRAUGHT EFFICIENCY    *
!     * PARAMETER *ZALPHA* ARE RESOLUTION DEPENDENT.  SEE COMMENTS     *
!     * BELOW WHERE THESE PHYSICAL CONSTANTS ARE DEFINED.	       *
!     ******************************************************************
!     * N.B.  THE DOWNDRAUGHT MODIFICATION OF THE BOUNDARY LAYER USES  *
!     * LEVEL-SPECIFIC CODE.  SEE COMMENTS BELOW WHERE THE PHYSICAL    *
!     * CONSTANT *NPBL* IS DEFINED.				       *
!     ******************************************************************
!
!     EXTERNALS.
!     ----------
!
!	   *ABORT*   CAUSE A FATAL ERROR CONDITION.
!	   *ALLOCA*  ALLOCATE ARRAY SPACE (MEMORY MANAGER ROUTINE).
!	   *GATHER*  GATHER SELECTED ARRAY ELEMENTS INTO A VECTOR.
!		     (CRAY LIBRARY ROUTINE).
!	   *SCATTER* SCATTER A VECTOR INTO SELECTED ARRAY ELEMENTS.
!		     (CRAY LIBRARY ROUTINE).
!	   *UNLOC*   FREE ARRAY SPACE (MEMORY MANAGER ROUTINE).
!
!     REFERENCES.
!     -----------
!
!	   A BASIC DESCRIPTION OF THE METHOD CAN BE FOUND IN BETTS
!     (1986 : QJRMS 112, 677-691) AND BETTS & MILLER (1986 : QJRMS,
!     112, 693-709).
!	   A MORE DETAILED DESCRIPTION OF THE SCHEME IS CONTAINED IN
!     BETTS & MILLER (1994 : AMERICAN METEOROLOGICAL SOCIETY MONOGRAPH
!     ON CONVECTIVE PARAMETRIZATION).
!	   IMPLEMENTATION OF THE SCHEME IN THE UGAMP GCM IS DESCRIBED
!     IN UGAMP TECHNICAL REPORT NO. 25 (SLINGO & BLACKBURN, 1992).
!
!	   SATURATION POINT CALCULATIONS FOLLOW BOLTON (1980, MON. WEA.
!     REV., 108, 1046-1053).  SATURATION SPECIFIC HUMIDITY CALCULATIONS
!     USE THE TETENS FORMULA (LOWE, 1977, J.APPL.MET., 16, 100-103).
!
!     -------------------------
      LOGICAL LO,LOA,LOFR,LOLIM
!!$	 LOGICAL LOLIFT,LOMELT,LLOC
!     -------------------------
!
      REAL CETAH (NLEV+1)
!
      REAL								    &
	   APP1	 (NLON,NLEV)						    &
	  ,TSM1M (NLON)							    &
	  ,RSFC	 (NLON)							    &
	  ,SSFC	 (NLON)							    &
	  ,PDPP1 (NLON,NLEV)						    &
	  ,PDPKPK(NLON,NLEV)						    &
	  ,PTP1	 (NLON,NLEV)						    &
	  ,PQP1	 (NLON,NLEV)						    &
	  ,PTC	 (NLON,NLEV)						    &
	  ,PQC	 (NLON,NLEV)						    &
	  ,PDT	 (NLON,NLEV)						    &
	  ,PDQ	 (NLON,NLEV)
!
      real dztref0(nlon, nlev)
      real dztref1(nlon, nlev)
      real dztref(nlon, nlev)
!
      INTEGER								    &
	   KDX	 (NLON)							    &
	  ,KBASE (NLON)							    &
	  ,KTOP	 (NLON)							    &
	  ,KLIFT (NLON)
!
      logical, dimension(ndeep, nlev) ::				    &
	   lloc								  
!
      real, dimension(ndeep,nlev) ::					    &
	   zpp1d							    &
	  ,zdpd								    &
	  ,zdpkd							    &
	  ,ztp1d							    &
	  ,zqp1d							    &
	  ,ztcd								    &
	  ,zqcd								    &
	  ,ztref0							    &
	  ,ztref1							    &
	  ,ztref							    &
	  ,zqref							    &
	  ,zsp
!
      integer, dimension(ndeep) ::					    &
	   itopd							    &
	  ,ibased							    &
	  ,iliftd							    &
	  ,icount							    &
	  ,ifreez							   
!      
      logical, dimension(ndeep) ::					    &
	   lolift							    &
	  ,lomelt
!
      real,  dimension(ndeep) ::					    &
	   ztsd								    &
	  ,zrain							    &
	  ,zsnow							    &
	  ,zdtd								    &
	  ,zdqd								    &
	  ,zpfr								    &
	  ,zptop							    &
	  ,zdtfr							    &
	  ,zdttop							    &
	  ,zdhtot							    &
	  ,zdptot							    &
	  ,zrinc							    &
	  ,zrtpbl							   
!      

      real, dimension(ndeep) ::						    &
	   zalph							    &
	  ,zepbl							    &
	  ,zfpbl							    &
	  ,zlfgbl 							    &
	  ,zlv								    &
	  ,zc3es							    &
	  ,zc4es
!
      real :: xlg, arg1, arg2
      XLG(ARG1,ARG2)=EXP(ARG2*ALOG(ARG1))
!
!*    local variables
!     ----- ---------
!
      real ::        ztau,     zalpha,   zalphs,   zdsplow,  zdspfr         &
          ,zdsptop,  zdsphi,   zdspel,   zstab,    zflim,    zetpbl         &
          ,zepcor,   zepq,     zepsec,   zepsrn,   zepvap                   &
          ,zcons1,   zcons2,   zcons3,   zc1,      zc2,      zc3            &
          ,zc4,      zc5,      zc6,      zc7,      zqsatc,   zcor           &
          ,zdqcd,    zndp,     zzdspfr,  zdsp,     ztsp,     zdevap         &
          ,zfi,      zf,       zgi,      zg,       zzdp,     zzdsptop       &
          ,ztsp1,    zqsp,     zdhdt,    zzdq,     zzq,      zzt            &
          ,zp1,      zzsp,     zzq1,     zzt1,     zp11,     zzsp1          &
          ,zzdsp,    zdqev,    zdr,      zrincm

      integer ::     npbl,     jk,       nlevm2,   nlevm3,   nupbl          &
          ,nddlev,   nhitop,   nhibas,   jl,       init,     jit            &
          ,isum,     ilfr,     iit,      iif,      ilfrm1,   nlevm

!*    PHYSICAL CONSTANTS.
!     -------- ----------
!
!     *ZTAU*	  ADJUSTMENT TIMESCALE IN SECONDS.
!		  ****************************************************
!		  * THE ADJUSTMENT TIMESCALE IS RESOLUTION DEPENDENT.*
!		  * *ZTAU* MUST BE SUFFICIENTLY SMALL TO MAINTAIN    *
!		  * SUBSATURATION IN THE REGIONS OF STRONGEST ASCENT.*
!		  * BOTH *ZTAU* AND DOWNDRAUGHT EFFICIENCY FACTOR    *
!		  * *ZALPHA* HAVE BEEN TESTED AND TUNED IN THE UGAMP *
!		  * GCM ONLY AT HORIZONTAL RESOLUTIONS UP TO T42.    *
!		  *		 T21	   T42			     *
!		  * *ZTAU*	14400	   7200			     *
!		  * *ZALPHA*	 0.15	   0.15			     *
!		  ****************************************************
!     *ZALPHA*	  DOWNDRAUGHT EFFICIENCY PARAMETER.
!     *ZALPHS*	  FACTOR TO REDUCE DOWNDRAUGHT EVAPORATION EFFICIENCY
!		  FOR SNOW.
!     *ZDSPLOW*	  SATURATION POINT DIFFERENCE FOR THE LOWEST LEVEL.
!     *ZDSPFR*	  SATURATION POINT DIFFERENCE FOR THE FREEZING LEVEL.
!     *ZDSPTOP*	  SATURATION POINT DIFFERENCE AT CLOUD TOP.
!     *ZDSPHI*	  SATURATION POINT DIFFERENCE FOR UPPER LEVEL LIFT.
!     *ZDSPEL*	  LIMIT FOR SATURATION POINT DIFFERENCE IN SUB-CLOUD
!		  LAYER DUE TO EVAPORATION (HIGH-LIFT POINTS).
!     *ZSTAB*	  FRACTION OF MOIST ADIABAT SLOPE USED FOR REFERENCE
!		  PROFILE.
!     *ZFLIM*	  LIMIT OF *ZF* IN DOWNDRAUGHT EVAPORATION CALCULATION.
!     *ZETPBL*	  HYBRID COORD VALUE AT THE TOP OF THE SPECIFIED
!		  PLANETARY BOUNDARY LAYER (PBL).  ALL LAYERS WITH
!		  HALF-LEVEL ETA GREATER THAN THIS ARE DESIGNATED
!		  PART OF THE PBL.
!     *NPBL*	  NUMBER OF LEVELS IN PLANETARY BOUNDARY LAYER,
!		  USED FOR DOWNDRAUGHT CALCULATIONS.
!		  **************************************************
!		  * THE DOWNDRAUGHT CALCS CONTAIN LEVEL-SPECIFIC   *
!		  * CODE WHICH REQUIRES *NPBL*=3.  ALL OTHER VALUES*
!		  * SWITCH OFF THE DOWNDRAUGHT PART OF THE SCHEME. *
!		  * *NPBL* IS COMPUTED FROM *ZETPBL* TO CHECK THAT *
!		  * THE IMPLIED CLOUD-BASE IS REASONABLE.	   *
!		  * THE SCHEME HAS ONLY BEEN TUNED FOR *NPBL*=3.   *
!		  **************************************************
!     *NUPBL*	  FIRST CLOUD LEVEL ABOVE SPECIFIED PBL.
!     *NDDLEV*	  DOWNDRAUGHT INFLOW LEVEL.

      ZTAU = CBMTD
      ZALPHA=CBMEF
      ZALPHS=0.1
      ZDSPLOW=-2500.
      ZDSPFR =-4000.
      ZDSPTOP=-2000.
      ZDSPHI =-4000.
      ZDSPEL =-2000.
      ZSTAB=0.85
      ZFLIM=-0.5
!      ZETPBL=0.825 ! what I've been using for L20 runs with equal sigma levels
      ZETPBL=0.9  ! this is the default
      NPBL=0
      DO 10 JK=NLEV,1,-1
      IF (CETAH(JK).GT.ZETPBL) NPBL=NPBL+1
   10 CONTINUE
!
      if(int(alv).eq.0) npbl = 0 ! turn off boundary layer for dry runs
!
      NLEVM2=NLEV-2
      NLEVM3=NLEV-3
      NUPBL=NLEV-NPBL
      NDDLEV=NUPBL-1
!
!!$      IF (NPBL.NE.3) THEN
!!$	 write(*,*) ' **************************************************'
!!$	 write(*,*) ' * ABORT IN BETTS-MILLER DEEP CONVECTION ROUTINE. *'
!!$	 write(*,*) ' * THE DOWNDRAUGHT CALCS CONTAIN LEVEL-SPECIFIC   *'
!!$	 write(*,*) ' * CODE WHICH REQUIRES *NPBL*=3.  ALL OTHER VALUES*'
!!$	 write(*,*) ' * SWITCH OFF THE DOWNDRAUGHT PART OF THE SCHEME. *'
!!$	 write(*,*) ' * *NPBL* IS COMPUTED FROM *ZETPBL* TO CHECK THAT *'
!!$	 write(*,*) ' * THE IMPLIED CLOUD-BASE IS REASONABLE.	       *'
!!$	 write(*,*) ' * THE SCHEME HAS ONLY BEEN TUNED FOR *NPBL*=3.   *'
!!$	 write(*,*) ' **************************************************'
!!$	 write(*,*) ' COMPUTED VALUE OF NPBL = ',NPBL
!!$	 CALL ABORT
!!$      ENDIF
!
!*    SECURITY PARAMETERS.
!     -------- -----------
!
!     *ZEPCOR*	  MINIMUM VALUE OF DENOMINATOR IN QSAT CALCULATION.
!     *ZEPQ*	  MINIMUM SPECIFIC HUMIDITY TO AVOID DIVERGENCE OF THE
!		  SATURATION POINT CALCULATIONS.
!     *ZEPSEC*	  IS A SECURITY FOR THE FREEZING LEVEL PRESSURE TO BE
!		  NOT EXACTLY IDENTICAL TO THE CLOUD TOP PRESSURE.
!     *ZEPSRN*	  ENSURES THAT ALL RAIN DOES NOT EVAPORATE UNDER HIGH
!		  BASES.
!     *ZEPVAP*	  MINIMUM VALUE OF DOWNDRAUGHT EVAPORATION INTEGRAL
!		  OVER THE PBL (PA.KG/KG).
!
      ZEPCOR=1.E-10
      ZEPQ=0.000002
      ZEPSEC=0.01
      ZEPSRN=1.E-9
      ZEPVAP=0.01
!
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------
!
!     *ZC1-ZC5*	  CONSTANTS FOR SATURATION POINT CALCULATIONS.
!     *ZC6*	  (RD/RV) USED FOR VAPOUR PRESSURE.
!     *ZC7*	  (1-RD/RV) USED FOR VAPOUR PRESSURE.
!
      ZCONS1=RD/CPD
      ZCONS2=ALV/CPD
      ZCONS3=1./ZTAU
!
      ZC1=CPD/RD
      ZC2=55.
      ZC3=2840.
      ZC4=3.5
      ZC5=0.2
      ZC6=0.622
      ZC7=0.378
!
!     ------------------------------------------------------------------
!
!*	   1.	  ALLOCATE SPACE AND POSITION VARIABLES.
!		  -------- ----- --- -------- ----------
!		  (removed in favor of automatic arrays)
!
  100 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   2.	  COLLECT DEEP CONVECTIVE POINTS.
!		  ------- ---- ---------- -------
!
  200 CONTINUE
!
      CALL GATHER_ir(NDEEP,ZTSD  ,TSM1M,KDX)
      CALL GATHER_ii(NDEEP,ITOPD ,KTOP ,KDX)
      CALL GATHER_ii(NDEEP,IBASED,KBASE,KDX)
      CALL GATHER_ii(NDEEP,ILIFTD,KLIFT,KDX)
!
      NHITOP=NLEV
      NHIBAS=NLEV
      DO 201 JL=1,NDEEP
      NHITOP=MIN(NHITOP,ITOPD(JL))
      NHIBAS=MIN(NHIBAS,IBASED(JL))
  201 CONTINUE

      DO 202 JK=NHITOP,NLEV
      CALL GATHER_ir(NDEEP,ZPP1D(:,JK),APP1  (:,JK),KDX)
      CALL GATHER_ir(NDEEP,ZDPD (:,JK),PDPP1 (:,JK),KDX)
      CALL GATHER_ir(NDEEP,ZDPKD(:,JK),PDPKPK(:,JK),KDX)
      CALL GATHER_ir(NDEEP,ZTP1D(:,JK),PTP1  (:,JK),KDX)
      CALL GATHER_ir(NDEEP,ZQP1D(:,JK),PQP1  (:,JK),KDX)
      CALL GATHER_ir(NDEEP,ZTCD (:,JK),PTC   (:,JK),KDX)
      CALL GATHER_ir(NDEEP,ZQCD (:,JK),PQC   (:,JK),KDX)

  202 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   3.	  PRELIMINARY REFERENCE PROFILE.
!		  ----------- --------- --------
!
  300 CONTINUE
!
!*	   3.1	  COMPLETE MOIST ADIABAT BELOW CLOUD BASE,
!		  WITH SECOND ITERATION AT EACH LEVEL.
!
  310 CONTINUE
!
      IF (NHIBAS.LE.NLEVM1) THEN
	 INIT=2
	 DO 313 JK=NHIBAS,NLEV
	 DO 312 JIT=1,INIT
	 DO 311 JL=1,NDEEP
	 LO=JK.GE.IBASED(JL)
	 IF (LO) THEN
	    IF (JIT.EQ.1) THEN
	       ZTCD(JL,JK)=ZTCD(JL,JK-1)/ZDPKD(JL,JK-1)
	       ZQCD(JL,JK)=ZQCD(JL,JK-1)
	    ENDIF
	    ZQSATC=C2ES*EXP(C3LES*(ZTCD(JL,JK)-TMELT)*			    &
			   (1./(ZTCD(JL,JK)-C4LES)))/(ZPP1D(JL,JK))
	    ZCOR=1./MAX(ZEPCOR,(1.-VTMPC1*ZQSATC))
	    ZQSATC=ZQSATC*ZCOR
	    ZDQCD=(ZQCD(JL,JK)-ZQSATC)/(1.+C5LES*ZCONS2*ZQSATC*ZCOR*	    &
				       (1./(ZTCD(JL,JK)-C4LES))**2)
	    ZTCD(JL,JK)=ZTCD(JL,JK)+ZCONS2*ZDQCD
	    ZQCD(JL,JK)=ZQCD(JL,JK)-ZDQCD
	 ENDIF
  311	 CONTINUE
  312	 CONTINUE
  313	 CONTINUE
      ENDIF
!
!*	   3.2	  BOUNDARY LAYER REFERENCE PROFILE, BASED ON A
!		  DOWNDRAUGHT FED BY ENVIRONMENTAL AIR AT *NDDLEV*.
!
  320 CONTINUE
!
      IF (NUPBL.EQ.NLEVM3)  THEN
	 DO 322 JK=NUPBL+1,NLEV
	 DO 321 JL=1,NDEEP
	 ZTREF(JL,JK)=ZTP1D(JL,NDDLEV)+(ZTCD(JL,JK)-ZTCD(JL,NDDLEV))
	 ZQREF(JL,JK)=ZQP1D(JL,NDDLEV)+(ZQCD(JL,JK)-ZQCD(JL,NDDLEV))
  321	 CONTINUE
  322	 CONTINUE
      ENDIF
!
!*	   3.3	  PRESET REFERENCE PROFILE ABOVE BOUNDARY LAYER.
!
  330 CONTINUE
!
      DO 332 JK=NHITOP,NUPBL
      DO 331 JL=1,NDEEP
      ZTREF(JL,JK)=ZTP1D(JL,JK)
      ZQREF(JL,JK)=ZQP1D(JL,JK)
  331 CONTINUE
  332 CONTINUE
!
!	   3.4	  FOR PARCELS LIFTED FROM LOW LEVELS, CONSTRUCT THE
!		  TEMPERATURE PROFILE FROM CLOUD BASE UP TO FREEZING
!		  LEVEL.  APPROXIMATE A VIRTUAL MOIST ADIABAT BY USING
!		  A FRACTION OF THE MOIST ADIABAT STATIC STABILITY.
!		  OMIT PARCELS LIFTED FROM HIGHER LEVELS, BY SETTING 
!		  NOMINAL FREEZING LEVEL TO THE FIRST IN-CLOUD LEVEL.
!
!		  *ICOUNT=0* BELOW NOMINAL FREEZING LEVEL.

!		  *IFREEZ*   INDICATES NOMINAL FREEZING LEVEL.
!
  340 CONTINUE
!
      DO 341 JL=1,NDEEP
      LOLIFT(JL)=(ILIFTD(JL)-1).GE.(NLEV-4)
      IF (LOLIFT(JL)) THEN
	 ICOUNT(JL)=0
	 IFREEZ(JL)=1
      ELSE
	 ICOUNT(JL)=1
	 IFREEZ(JL)=IBASED(JL)-1
      ENDIF
!
!		  ***********************************************
!		  START OF REFERENCE PROFILE CAN BE VARIED.
!		  TUNED AT T42L19 TO USE MEAN OF ENVIRONMENT AND
!		  CLOUD TEMPERATURES AT THE FIRST IN-CLOUD LEVEL.
!		  ***********************************************
!                  Without a boundary level (nbpl=0, nupbl=nlev),
!                  ztref(jl,nupbl) defaults to ztd1d(jl,nupbl), or
!                  the lowest-level temperature

      LO=(ICOUNT(JL).EQ.0).AND.(NUPBL.EQ.NLEVM3)
      IF (LO) ZTREF(JL,NUPBL)=0.5*(ZTCD(JL,NUPBL)+ZTP1D(JL,NUPBL))
  341 CONTINUE
!
      NLEVM=NUPBL-1
      DO 343 JK=NLEVM,NHITOP,-1
!
      ISUM=0
      DO 342 JL=1,NDEEP
      IF (ICOUNT(JL).EQ.0) THEN
! BE MINDFUL OF THIS CHANGE!
         ztref(jl,jk)= ZTREF(JL,JK+1) +                           &
              0.95*( ZTREF(JL,JK+1)*ZDPKD(JL,JK) - ZTREF(JL,JK+1) )
! this is the original code:
!!$	 ZTREF(JL,JK)=ZTREF(JL,JK+1)*ZDPKD(JL,JK)+			    &
!!$		     ZSTAB*(ZTCD(JL,JK)-ZTCD(JL,JK+1)*ZDPKD(JL,JK))
      ENDIF
! BE MINDFUL OF THIS CHANGE
      LOFR=(ZTREF(JL,JK).LE.0.0).OR.(JK.LE.ITOPD(JL))
!!$      LOFR=(ZTREF(JL,JK).LE.TMELT).OR.(JK.LE.ITOPD(JL)) ! original code
      IF (LOFR.AND.IFREEZ(JL).EQ.1) THEN
	 ICOUNT(JL)=1
	 IFREEZ(JL)=JK
      ENDIF
      ISUM=ISUM+ICOUNT(JL)
  342 CONTINUE
!
!		  SKIP OUT OF LEVEL-LOOP IF ALL FREEZING LEVELS
!		  HAVE BEEN FOUND.
!
      IF (ISUM.EQ.NDEEP) GO TO 344
!
  343 CONTINUE
!
  344 CONTINUE
!
!*	   3.5	  CONSTRUCT TEMPERATURE PROFILE ABOVE FREEZING LEVEL,
!		  BY INTERPOLATING (QUADRATICALLY) BETWEEN TEMPERATURE
!		  DEFICITS RELATIVE TO THE MOIST ADIABAT AT FREEZING
!		  LEVEL AND CLOUD TOP.
!
!		  SINCE POINT OF NEUTRAL STABILITY DOES NOT NECESSARILY
!		  FALL ON A MODEL LEVEL AN ASSUMPTION ON ITS POSITION
!		  HAS TO BE MADE.
!	      *** CLOUD TOP ADIABATIC DEFICIT SET TO ZERO: NEVER USED.
!
  350 CONTINUE
!
      ILFR=NHITOP
      DO 351 JL=1,NDEEP
      IIT=ITOPD(JL)
      ZPTOP(JL)=ZPP1D(JL,IIT)
!***  ZDTTOP(JL)=(ZTCD(JL,IIT)-ZTP1D(JL,IIT))
      ZDTTOP(JL)=0.0
      IIF=IFREEZ(JL)
      ZPFR(JL)=ZPP1D(JL,IIF)+ZEPSEC
      ZDTFR(JL)=(ZTCD(JL,IIF)-ZTREF(JL,IIF))
      ILFR=MAX(ILFR,IIF)
  351 CONTINUE
!
!		  *ILFR* IS LOWEST FREEZING LEVEL AROUND LATITUDE.
!
      IF (ILFR.GT.NHITOP) THEN
	 ILFRM1=ILFR-1
	 DO 353 JK=NHITOP,ILFRM1
	 DO 352 JL=1,NDEEP
	 LO=(JK.LE.IFREEZ(JL)).AND.(JK.GE.ITOPD(JL)) ! original code
! BE MINDFUL OF THIS CHANGE--turns off interpolation to DALR at cloud top
         lo =.false.
	 IF (LO) THEN
	    ZNDP=(ZPFR(JL)-ZPP1D(JL,JK))/(ZPFR(JL)-ZPTOP(JL))
	    ZTREF(JL,JK)=ZTCD(JL,JK)-ZDTFR(JL)*(1.0-ZNDP*ZNDP)
	 ENDIF
  352	 CONTINUE
  353	 CONTINUE
      ENDIF
!
!*	   3.6	  SATURATION POINTS AND REFERENCE HUMIDITY
!		  FOR ALL IN-CLOUD LEVELS ABOVE BOUNDARY LAYER.
!
!		  RETAIN PRESET ENVIRONMENT PROFILE ABOVE CLOUD
!		  AND BELOW NOMINAL FREEZING LEVEL (CLOUD-BASE)
!		  FOR CLOUDS LIFTED FROM HIGHER LEVELS.
!
  360 CONTINUE
!
      DO 362 JK=NHITOP,NUPBL
      DO 361 JL=1,NDEEP
      IF (LOLIFT(JL)) THEN
	 ZZDSPFR=ZDSPFR
	 ZZDSPTOP=ZDSPTOP
      ELSE
	 ZZDSPFR=ZDSPHI
	 ZZDSPTOP=ZDSPHI
      ENDIF
      LOFR=JK.LE.IFREEZ(JL)
      IF (LOFR) THEN
	 ZDSP=ZZDSPFR-(ZZDSPFR-ZZDSPTOP)*(ZPFR(JL)-ZPP1D(JL,JK))	    &
					/(ZPFR(JL)-ZPTOP(JL))
      ELSE
	 ZDSP=ZDSPLOW-(ZDSPLOW-ZDSPFR)*(ZPP1D(JL,NUPBL)-ZPP1D(JL,JK))	    &
				      /(ZPP1D(JL,NUPBL)-ZPFR(JL))
      ENDIF
      LLOC(JL,JK)=(JK.GE.ITOPD(JL)).AND.(LOFR.OR.LOLIFT(JL))
      IF (LLOC(JL,JK)) THEN
	 ZSP(JL,JK)=ZPP1D(JL,JK)+ZDSP
	 ZTSP=ZTREF(JL,JK)*XLG((ZSP(JL,JK)/ZPP1D(JL,JK)),ZCONS1)
	 ZQREF(JL,JK)=C2ES*EXP(C3LES*(ZTSP-TMELT)/(ZTSP-C4LES))		    &
		      /ZSP(JL,JK)
      ENDIF
  361 CONTINUE
  362 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   4.	  CONSIDERATION OF ENERGY CONSERVATION.
!		  ------------- -- ------ -------------
!
  400 CONTINUE
!
!*	   4.1	  COMPUTE THE DOWNDRAUGHT EVAPORATION INTEGRAL
!		  OVER THE BOUNDARY LAYER.
!		  SET LOGICAL ARRAY FOR SURFACE TEMPERATURE.
!
  410 CONTINUE
!
      DO 411 JL=1,NDEEP
      ZDEVAP=(ZQCD(JL,NLEV  )-ZQCD(JL,NDDLEV))*ZDPD(JL,NLEV  )+		    &
	     (ZQCD(JL,NLEVM1)-ZQCD(JL,NDDLEV))*ZDPD(JL,NLEVM1)+		    &
	     (ZQCD(JL,NLEVM2)-ZQCD(JL,NDDLEV))*ZDPD(JL,NLEVM2)
      ZEPBL(JL)=MAX(ZEPVAP,ZDEVAP)
      LOMELT(JL)=ZTSD(JL).GT.TMELT
      IF (LOMELT(JL)) THEN
	 ZLV(JL)=ALV
	 ZC3ES(JL)=C3LES
	 ZC4ES(JL)=C4LES
      ELSE
	 ZLV(JL)=ALS
	 ZC3ES(JL)=C3IES
	 ZC4ES(JL)=C4IES
      ENDIF
  411 CONTINUE
!
!*	   4.2	  BEGIN LOOP FOR TWO ITERATIONS OF ENERGY CORRECTION.
!
  420 CONTINUE
!
      DO 461 JIT=1,2
!
!*	   4.3	  INITIALISE ENERGY AND PRESSURE INTEGRALS.
!		  COMPUTE BOUNDARY LAYER INTEGRALS.
!
  430 CONTINUE
!
      DO 431 JL=1,NDEEP
      ZDHTOT(JL)=0.
      ZDPTOT(JL)=0.
      IF (LOMELT(JL)) THEN
	 ZALPH(JL)=ZALPHA
      ELSE
	 ZALPH(JL)=ZALPHS*ZALPHA
      ENDIF
      ZFI=(ZQREF(JL,NLEV  )-ZQP1D(JL,NLEV  ))*ZDPD(JL,NLEV  )+		    &
	  (ZQREF(JL,NLEVM1)-ZQP1D(JL,NLEVM1))*ZDPD(JL,NLEVM1)+		    &
	  (ZQREF(JL,NLEVM2)-ZQP1D(JL,NLEVM2))*ZDPD(JL,NLEVM2)
      ZF=ZALPH(JL)*ZFI/ZEPBL(JL)
      LOLIM=ZF.LT.ZFLIM
      IF (LOLIM) THEN
	 ZALPH(JL)=ZALPH(JL)*(ZFLIM/ZF)
	 ZF=ZFLIM
      ENDIF
      ZGI=(ZTREF(JL,NLEV  )-ZTP1D(JL,NLEV  ))*ZDPD(JL,NLEV  )+		    &
	  (ZTREF(JL,NLEVM1)-ZTP1D(JL,NLEVM1))*ZDPD(JL,NLEVM1)+		    &
	  (ZTREF(JL,NLEVM2)-ZTP1D(JL,NLEVM2))*ZDPD(JL,NLEVM2)
      if(int(zlv(jl)).ne.0) then ! prevent NaN (zcons2=0 when als=0)
	 ZG=ZALPH(JL)*ZGI/(ZEPBL(JL)*ZCONS2)
	 ZLFGBL(JL)=ZLV(JL)*(1.0-ZG)/(1.0+ZF)
      else
	 zlfgbl(jl)=0.0
      endif
  431 CONTINUE
!
!*	   4.4	  COMPUTE ENERGY INTEGRAL.
!
  440 CONTINUE
!
      DO 442 JK=NHITOP,NUPBL
      DO 441 JL=1,NDEEP
      IF (LLOC(JL,JK)) THEN
	 ZDPTOT(JL)=ZDPTOT(JL)+ZDPD(JL,JK)
	 ZDHTOT(JL)=ZDHTOT(JL)+ZDPD(JL,JK)*				    &
	      (CPD*(1.+VTMPC2*ZQP1D(JL,JK))*(ZTP1D(JL,JK)-ZTREF(JL,JK))+    &
	       ZLFGBL(JL)*(ZQP1D(JL,JK)-ZQREF(JL,JK)))
      ENDIF
  441 CONTINUE
  442 CONTINUE
!
!		  NORMALISE INTEGRAL, AVOIDING ZERO PRESSURE INTERVAL.
!
      DO 443 JL=1,NDEEP
      ZZDP=MAX(ZEPSEC,ZDPTOT(JL))
      ZDHTOT(JL)=ZDHTOT(JL)/ZZDP
  443 CONTINUE
!
!*	   4.5	  MODIFY ORIGINAL PROFILE FOR ENERGY CONSERVATION.
!		  INCLUDE CLOUDTOP TO ALLOW CORRECTION WHEN ZDTTOP=0.
!
  450 CONTINUE
!
      DO 452 JK=NHITOP,NUPBL
      DO 451 JL=1,NDEEP
      IF (LLOC(JL,JK)) THEN
	 ZTSP1=(ZTREF(JL,JK)+1.)*(XLG((ZSP(JL,JK)/ZPP1D(JL,JK)),ZCONS1))
	 ZQSP=C2ES*EXP(ZC3ES(JL)*(ZTSP1-TMELT)/(ZTSP1-ZC4ES(JL)))	    &
	      /ZSP(JL,JK)
	 ZDHDT=CPD*(1.+VTMPC2*ZQP1D(JL,JK))+ZLV(JL)*(ZQSP-ZQREF(JL,JK))
         if(jit.eq.2) then
            ztref1(jl,jk)=ztref(jl,jk)
         elseif(jit.eq.1) then
            ztref0(jl,jk)=ztref(jl,jk)
         endif
	 ZTREF(JL,JK)=ZTREF(JL,JK)+ZDHTOT(JL)/ZDHDT
	 ZTSP=ZTREF(JL,JK)*(XLG((ZSP(JL,JK)/ZPP1D(JL,JK)),ZCONS1))
	 ZQREF(JL,JK)=C2ES*EXP(ZC3ES(JL)*(ZTSP-TMELT)/(ZTSP-ZC4ES(JL)))	    &
		      /ZSP(JL,JK)
      ENDIF
  451 CONTINUE
  452 CONTINUE
!
!*	   4.6	  END OF LOOP OVER ITERATIONS.
!
  461 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   5.	  PRECIPITATION AND EVAPORATION, FINAL TENDENCIES.
!		  ------------- --- ------------ ----- -----------
!
  500 CONTINUE
!
!*	   5.1	  PRESET ARRAYS.
!		  COMPUTE A BOUNDARY LAYER MOISTURE INTEGRAL.
!
  510 CONTINUE
!
      DO 511 JL=1,NDEEP
      ZRAIN(JL)=0.
      ZRINC(JL)=0.
      IF (LOMELT(JL)) THEN
	 ZALPH(JL)=ZALPHA
      ELSE
	 ZALPH(JL)=0.1*ZALPHA
      ENDIF
      ZFI=(ZQREF(JL,NLEV  )-ZQP1D(JL,NLEV  ))*ZDPD(JL,NLEV  )+		    &
	  (ZQREF(JL,NLEVM1)-ZQP1D(JL,NLEVM1))*ZDPD(JL,NLEVM1)+		    &
	  (ZQREF(JL,NLEVM2)-ZQP1D(JL,NLEVM2))*ZDPD(JL,NLEVM2)
      ZF=ZALPH(JL)*ZFI/ZEPBL(JL)
      LOLIM=ZF.LT.ZFLIM
      IF (LOLIM) THEN
	 ZALPH(JL)=ZALPH(JL)*(ZFLIM/ZF)
	 ZF=ZFLIM
      ENDIF
      ZFPBL(JL)=ZF
  511 CONTINUE
!
!*	   5.2	  COMPUTE PRECIPITATION FROM INTEGRAL OF MOISTURE
!		  TENDENCY OVER CLOUD LEVELS.
!		  *ZRAIN* CONTAINS TOTAL PRECIPITATION IN SECTIONS
!		  5.2, 5.3, IRRESPECTIVE OF SURFACE TEMPERATURE.
!
  520 CONTINUE
!
      DO 522 JK=NHITOP,NUPBL
      DO 521 JL=1,NDEEP
      ZZDQ=(ZQREF(JL,JK)-ZQP1D(JL,JK))*ZCONS3/(1.0+ZFPBL(JL))
      ZRAIN(JL)=ZRAIN(JL)-ZZDQ*ZDPD(JL,JK)/G
  521 CONTINUE
  522 CONTINUE
!
!*	   5.3	  EVAPORATION BELOW CLOUD-BASE FOR PARCELS LIFTED
!		  FROM HIGH LEVELS.
!
  530 CONTINUE
!
      if(int(als).ne.0) then ! don't execute this code for dry runs--doesn't
                             ! matter for T, since there is no evaporation, but
                             ! CPD/ALV and ZZDSP/((ZZSP1-ZZSP)*ZCONS2) terms 
                             ! cause problems in dry case.

      DO 532 JK=NHIBAS,NLEV
      DO 531 JL=1,NDEEP
      LO=     (JK.GE.IBASED(JL))					    &
	 .AND.(.NOT.LOLIFT(JL))						    &
	 .AND.(ZRINC(JL).LE.ZRAIN(JL))
      IF (LO) THEN
	 ZZQ=MAX(ZEPQ,ZQP1D(JL,JK))
	 ZZT=ZTP1D(JL,JK)
	 ZP1=ZZQ*ZPP1D(JL,JK)/(ZC6+ZC7*ZZQ)
	 ZTSP=ZC2+ZC3/(ZC4*ALOG(ZZT)-ALOG(ZP1)-ZC5)
	 ZZSP=ZPP1D(JL,JK)*XLG((ZTSP/ZZT),ZC1)
	 ZZQ1=ZZQ+CPD/ALV ! here
	 ZZT1=ZZT-1.
	 ZP11=ZZQ1*ZPP1D(JL,JK)/(ZC6+ZC7*ZZQ1)
	 ZTSP1=ZC2+ZC3/(ZC4*ALOG(ZZT1)-ALOG(ZP11)-ZC5)
	 ZZSP1=ZPP1D(JL,JK)*XLG((ZTSP1/ZZT1),ZC1)
	 ZZDSP=ZPP1D(JL,JK)-ZZSP+ZDSPEL
	 IF (ZZDSP.GT.0.) THEN
            ZDQEV=ZZDSP/((ZZSP1-ZZSP)*ZCONS2) ! here
	    ZDR=ZDPD(JL,JK)*ZDQEV*ZCONS3/G
	    ZRINCM=ZRINC(JL)
	    ZRINC(JL)=ZRINC(JL)+ZDR
	    LOA=(ZRAIN(JL).LT.ZRINC(JL)).AND.(ZRAIN(JL).GT.ZRINCM)
	    IF (LOA) THEN
	       ZDQEV=(ZRAIN(JL)-ZRINCM-ZEPSRN)*G/(ZCONS3*ZDPD(JL,JK))
	    ENDIF
	    ZTREF(JL,JK)=ZTREF(JL,JK)-ZDQEV*ZCONS2
	    ZQREF(JL,JK)=ZQREF(JL,JK)+ZDQEV
	 ENDIF
      ENDIF
  531 CONTINUE
  532 CONTINUE

      endif

!
!*	   5.4	  COMPUTE FINAL TENDENCIES.
!		  SEPARATE ADJUSTMENT TIMESCALE FOR LEVELS IN THE
!		  BOUNDARY LAYER AFFECTED BY THE DOWNDRAUGHT.
!		  RECALCULATE PRECIPITATION AFTER EVAPORATON.
!		  SCATTER RESULTS BACK TO ORIGINAL GRID.
!
  540 CONTINUE
!
      DO 541 JL=1,NDEEP
      ZRTPBL(JL)=(ZALPH(JL)*G*ZRAIN(JL))/ZEPBL(JL)
      ZRAIN(JL)=0.
      ZSNOW(JL)=0.
  541 CONTINUE
!
      DO 543 JK=NHITOP,NLEV
!
      LO=(JK.LE.NUPBL).OR.(NUPBL.NE.NLEVM3)
!
      DO 542 JL=1,NDEEP
      IF (LO) THEN
	 ZDTD(JL)=(ZTREF(JL,JK)-ZTP1D(JL,JK))*ZCONS3
	 ZDQD(JL)=(ZQREF(JL,JK)-ZQP1D(JL,JK))*ZCONS3
      ELSE
	 ZDTD(JL)=(ZTREF(JL,JK)-ZTP1D(JL,JK))*ZRTPBL(JL)
	 ZDQD(JL)=(ZQREF(JL,JK)-ZQP1D(JL,JK))*ZRTPBL(JL)
      ENDIF
      IF (LOMELT(JL)) THEN
	 ZRAIN(JL)=ZRAIN(JL)-ZDQD(JL)*ZDPD(JL,JK)/G
      ELSE
	 ZSNOW(JL)=ZSNOW(JL)-ZDQD(JL)*ZDPD(JL,JK)/G
      ENDIF
  542 CONTINUE
!
      CALL SCATTER(NDEEP,PDT(:,JK),KDX,ZDTD)
      CALL SCATTER(NDEEP,PDQ(:,JK),KDX,ZDQD)
      call scatter(ndeep,dztref0(:,jk), kdx, ztref0(:,jk))
      call scatter(ndeep,dztref1(:,jk), kdx, ztref1(:,jk))
      call scatter(ndeep,dztref(:,jk), kdx, ztref(:,jk))
  543 CONTINUE
!
!		  SCATTER PRECIPITATION BACK TO ORIGINAL GRID.
!		  CHECK FOR NEGATIVE PRECIP. OCCURS IN THE MAIN ROUTINE,
!		  WHERE SUCH POINTS ARE TREATED AS SHALLOW CONVECTION.
!		 
      CALL SCATTER(NDEEP,RSFC,KDX,ZRAIN)
      CALL SCATTER(NDEEP,SSFC,KDX,ZSNOW)
!
!     ------------------------------------------------------------------
!
!*	   6.	  RETURN WORKSPACE.
!		  ------ ----------
!
  600 CONTINUE
!
!

      RETURN
    END SUBROUTINE BMDEEP
!*/
!*/ ************************************************************** BMSHAL
!*/
      SUBROUTINE BMSHAL(						    &
	   NSHAL,    CBMTS,    KDX					    &
	  ,NLON,     NLEV,     NLEVM1,	 NLEVP1				    &
	  ,CPD,	     G,	       RD,	 TMELT				    &
	  ,C2ES,     C3LES,    C4LES					    &
	  ,APP1,     APHP1						    &
	  ,PTP1,     PQP1						    &
	  ,PDT,	     PDQ,      PSC					    &
	  ,KBASE,    KTOP						    &
	  ,KHCTOP,   KLOBAS,   sztref                                       &
	  )

        real, intent(in) ::    cbmts,    cpd,      g,        rd             &
          ,tmelt,    c2es,     c3les,    c4les    

        integer, intent(in) :: nshal,    nlon,     nlev,     nlevm1         &
          ,nlevp1

        integer, intent(out) :: khctop,   klobas

!
!**** *BMSHAL* - PERFORMS SATURATION POINT ADJUSTMENT
!		 FOR SHALLOW CONVECTION.
!
!     ORIGINAL VERSION		       B.RITTER	 E.C.M.W.F.  31/10/84.
!     IMPLEMENTED & TUNED IN UGCM1   J.M.SLINGO	 U.G.A.M.P.	12/91.
!     IMPLEMENTED IN UGCM2	    M.BLACKBURN	 U.G.A.M.P.  13/05/94.
!
!     PURPOSE.
!     --------
!
!	   THIS ROUTINE ADJUSTS THE TEMPERATURE AND SPECIFIC HUMIDITY
!     PROFILES OF SHALLOW CONVECTIVE POINTS TO AN INTERNALLY SPECIFIED
!     REFERENCE PROFILE.
!
!**   INTERFACE.
!     ----------
!
!	   *BMSHAL* IS CALLED FROM *BMADJ*.
!	   THE ROUTINE IS ARGUMENT DRIVEN, TAKING ITS INPUT ENTIRELY
!     FROM THE VARIABLES AND ARRAYS SUPPLIED IN THE ARGUMENT LIST.
!     FOR VECTORISATION, THE INPUT IS GATHERED IN TEMPORARY ARRAYS OF
!     LENGTH *NSHAL* (THE NUMBER OF GRIDPOINTS WITH SHALLOW CONVECTION
!     ON A LATITUDE ROW) USING INDEX ARRAY *KDX*.  THE ROUTINE RETURNS
!     ITS OUTPUT AS ARGUMENTS, HAVING SCATTERED TO THE FULL GRID.
!
!  ---------------------------------------------------------------------
!
!	   INPUT ARGUMENTS:
!
!     *NSHAL*	NUMBER OF GRIDPOINTS WITH SHALLOW CONVECTION IN ROW.
!     *CBMTS*	ADJUSTMENT TIMESCALE IN SECONDS FOR SHALLOW CONVECTION.
!     *KDX*	INDEX OF CONVECT. POINTS ON FULL GRID, DIMENSION (NLON).
!     *NLON*	NUMBER OF LONGITUDES.
!     *NLEV*	NUMBER OF LEVELS.
!     *NLEVM1*	(NLEV-1).
!     *NLEVP1*	(NLEV+1).
!     *CPD*	SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR.
!     *G*	GRAVITATIONAL ACCELERATION.
!     *RD*	GAS CONSTANT FOR DRY AIR.
!     *TMELT*	TEMPERATURE OF FUSION OF ICE.
!     *C__ES*	CONSTANTS USED FOR COMPUTATION OF SATURATION SPECIFIC
!		HUMIDITY OVER LIQUID WATER (*C_LES*) OR ICE (*C_IES*).
!     *C2ES*	(RD/RV)*(SVP AT REFERENCE TEMPERATURE C4_ES).
!     *C3_ES*	CONSTANT FOR SVP.
!     *C4_ES*	REFERENCE TEMPERATURE FOR SVP.
!     *APP1*	FULL LEVEL PRESSURE,		  DIMENSION (NLON,NLEV).
!     *APHP1*	HALF LEVEL PRESSURE,		DIMENSION (NLON,NLEVP1).
!     *PTP1*	PRELIMINARY (T+1) TEMPERATURE,	  DIMENSION (NLON,NLEV).
!     *PQP1*	PRELIMINARY (T+1) SPECIFIC HUMID, DIMENSION (NLON,NLEV).
!     *KBASE*	HIGHEST LEVEL IN SUB-CLOUD LAYER,      DIMENSION (NLON).
!     *KTOP*	HIGHEST LEVEL IN CLOUD LAYER,	       DIMENSION (NLON).
!
!	   OUTPUT ARGUMENTS:
!
!     *PDT*	TEMPERATURE TENDENCY,		  DIMENSION (NLON,NLEV).
!     *PDQ*	SPECIFIC HUMIDITY TENDENCY,	  DIMENSION (NLON,NLEV).
!     *PSC*	CONDENSATION RATE,		       DIMENSION (NLON).
!     *KHCTOP*	UPPERMOST LEVEL AT WHICH TENDENCIES ARE COMPUTED.
!     *KLOBAS*	LOWEST LEVEL AT WHICH TENDENCIES ARE COMPUTED.
!
!  ---------------------------------------------------------------------
!*ENDIF
!	   OUTPUT ARRAYS ARE ONLY ASSIGNED AT SHALLOW CONVECTING POINTS.
!     ALL ARRAY ELEMENTS SHOULD BE PRESET TO ZERO ON INPUT.
!
!     METHOD.
!     -------
!
!	   THE BASIC REFERENCE PROFILE IS DEFINED BY THE SATURATION
!     POINTS OF THE LOWEST LEVEL AND A LEVEL ABOVE THE INVERSION.
!     THESE DETERMINE THE MIXING LINE.	THE STRENGTH OF THE INVERSION
!     IS CONSIDERED AS WELL. THE REFERENCE PROFILES ARE CONSTRAINED
!     BY THE REQUIREMENT OF CONSERVATION OF SENSIBLE AND LATENT HEAT
!     SEPERATELY.
!	   FINALLY HEATING RATES AND MOISTURE CHANGES ARE CALCULATED
!     AS A RELAXATION TO THE ADJUSTMENT PROFILE.  NOTE THAT TENDENCIES
!     ARE APPLIED UP TO THE INVERSION LEVEL ABOVE THE INPUT CLOUD TOP.
!     ******************************************************************
!     * N.B.  THE RELAXATION TIME *ZTAU* IS RESOLUTION DEPENDENT.      *
!     * SEE COMMENTS BELOW WHERE THIS PHYSICAL CONSTANT IS DEFINED.    *
!     ******************************************************************
!
!     EXTERNALS.
!     ----------
!
!	   *ALLOCA*  ALLOCATE ARRAY SPACE (MEMORY MANAGER ROUTINE).
!	   *GATHER*  GATHER SELECTED ARRAY ELEMENTS INTO A VECTOR.
!		     (CRAY LIBRARY ROUTINE).
!	   *SCATTER* SCATTER A VECTOR INTO SELECTED ARRAY ELEMENTS.
!		     (CRAY LIBRARY ROUTINE).
!	   *UNLOC*   FREE ARRAY SPACE (MEMORY MANAGER ROUTINE).
!
!     REFERENCES.
!     -----------
!
!	   A BASIC DESCRIPTION OF THE METHOD CAN BE FOUND IN BETTS
!     (1986 : QJRMS 112, 677-691) AND BETTS & MILLER (1986 : QJRMS,
!     112, 693-709).
!	   A MORE DETAILED DESCRIPTION OF THE SCHEME IS CONTAINED IN
!     BETTS & MILLER (1994 : AMERICAN METEOROLOGICAL SOCIETY MONOGRAPH
!     ON CONVECTIVE PARAMETRIZATION).
!	   IMPLEMENTATION OF THE SCHEME IN THE UGAMP GCM IS DESCRIBED
!     IN UGAMP TECHNICAL REPORT NO. 25 (SLINGO & BLACKBURN, 1992).
!
!	   SATURATION POINT CALCULATIONS FOLLOW BOLTON (1980, MON. WEA.
!     REV., 108, 1046-1053).  SATURATION SPECIFIC HUMIDITY CALCULATIONS
!     USE THE TETENS FORMULA (LOWE, 1977, J.APPL.MET., 16, 100-103).
!
!     ------------------
      LOGICAL LO,LOA,LOB
!     ------------------
!
      REAL								    &
	   APP1	 (NLON,NLEV)						    &
	  ,APHP1 (NLON,NLEVP1)						    &
	  ,PTP1	 (NLON,NLEV)						    &
	  ,PQP1	 (NLON,NLEV)						    &
	  ,PDT	 (NLON,NLEV)						    &
	  ,PDQ	 (NLON,NLEV)						    &
	  ,PSC	 (NLON)
!
      real sztref(nlon,nlev)
!
      INTEGER								    &
	   KDX	 (NLON)							    &
	  ,KTOP	 (NLON)							    &
	  ,KBASE (NLON)
!
      real, dimension(nshal, nlev) ::				            &
	   zpp1s							    &
	  ,zphp1s							    &
	  ,ztp1s							    &
	  ,zqp1s							    &
	  ,ztref							    &
	  ,zqref							    &
	  ,zdp	 
!      
      integer, dimension(nshal) ::					    &
	   itops							    &
	  ,ibases							    
!
      real, dimension(nshal) ::						    &
	   zdts								    &
	  ,zdqs								    &
	  ,zmix								    &
	  ,zbinv							    &
	  ,ztsum							    &
	  ,zqsum							    &
	  ,zqint							  
!
      real :: xlg, arg1, arg2
      XLG(ARG1,ARG2)=EXP(ARG2*ALOG(ARG1))

!     local variables
!     ----- ---------
!
      real ::        ztau,     zspdiff,  zbitop,   zbmin,    zbmax          &
          ,zbshal,   zstabm,   zpzero,   zepdsp,   zepmix,   zepq           &
          ,zcons1,   zcons2,   zcons3,   zc1,      zc2,      zc3            &
          ,zc4,      zc5,      zc6,      zc7,      zptopm1                  &
          ,zttopm1,  zqtopm1,  zptopp1,  zqtopp1,  zdptop,   zqnlm1         &
          ,zqmx,     zthmx,    zzt,      zp1,      ztsp,     zsp1           &
          ,zspt,     zthbot,   zspti,    zttopp1,  ztmx,     zdpmix         &
          ,zspbi,    zbeta,    zsp,      zts,      zcond,    zspl 

      integer ::     jk,       jl,       itopm1,   itopp,    itopp1         &
           ,ibot
!
!*    PHYSICAL CONSTANTS.
!     -------- ----------
!
!     *ZTAU*	  RELAXATION TIMESCALE IN SECONDS.
!		  ****************************************************
!		  * THE ADJUSTMENT TIMESCALE IS RESOLUTION DEPENDENT *
!		  * *ZTAU* DETERMINES THE RATE AT WHICH MOISTURE IS  *
!		  * TRANSFERRED THROUGH THE PBL AND INTO THE FREE    *
!		  * ATMOSPHERE.	 IT HAS BEEN TESTED AND TUNED IN THE *
!		  * UGAMP GCM, TOGETHER WITH *ZTAU* AND *ZALPHA* IN  *
!		  * *BMDEEP*, ONLY AT HORIZONTAL RESOLUTIONS UP TO   *
!		  * T42.					     *
!		  *		 T21	   T42			     *
!		  * *ZTAU*	14400	  14400			     *
!		  ****************************************************
!     *ZSPDIFF*	  SATURATION POINT DIFFERENCE FOR THE LOWEST LEVEL.
!     *ZBITOP*	  MIXING PARAMETER FOR THE INVERSION.
!     *ZBMIN*	  MINIMUM LIMIT OF INVERSION MIXING PARAMETER.
!     *ZBMAX*	  MAXIMUM LIMIT OF INVERSION MIXING PARAMETER.
!     *ZBSHAL*	  MIXING PARAMETER.
!     *ZSTABM*	  MIXING LINE WEIGHT.
!     *ZPZERO*	  REFERENCE PRESSURE FOR DRY ADIABATS.
!
      ZTAU=CBMTS
      ZSPDIFF=-5000.
      ZBITOP=1.2
      ZBMIN=1.0
      ZBMAX=2.5
      ZBSHAL=1.0
      ZSTABM=0.85
      ZPZERO=100000.

!
!*    SECURITY CONSTANTS.
!     -------- ----------
!
!     *NLEVM1*	  (DUMMY ARGUMENT)  LOWEST CLOUD BASE ALLOWED.
!     *ZEPDSP*	  MINIMUM SATURATION POINT DIFFERENCE OVER LAYER,
!		  TO AVOID DIVIDE BY ZERO.
!     *ZEPMIX*	  MINIMUM VALUE OF MIXING LINE GRADIENT.
!     *ZEPQ*	  MINIMUM SPECIFIC HUMIDITY TO AVOID DIVERGENCE OF THE
!		  SATURATION POINT CALCULATIONS.
!
      ZEPDSP=-1.E-2
      ZEPMIX=-0.0007
      ZEPQ=0.000002
!
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------
!
!     *ZC1-ZC5*	  CONSTANTS FOR SATURATION POINT CALCULATIONS.
!     *ZC6*	  (RD/RV) USED FOR VAPOUR PRESSURE.
!     *ZC7*	  (1-RD/RV) USED FOR VAPOUR PRESSURE.
!
      ZCONS1=RD/CPD
      ZCONS3=1./ZTAU
!
      ZC1=CPD/RD
      ZC2=55.
      ZC3=2840.
      ZC4=3.5
      ZC5=0.2
      ZC6=0.622
      ZC7=0.378
!
!     ------------------------------------------------------------------
!
!*	   1.	  ALLOCATE SPACE AND POSITION VARIABLES.
!		  -------- ----- --- -------- ----------
!
  100 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   2.	  COLLECT SHALLOW CONVECTIVE POINTS.
!		  ------- ------- ---------- -------
!
  200 CONTINUE
!
      CALL GATHER_ii(NSHAL, ITOPS,  KTOP,  KDX)
      CALL GATHER_ii(NSHAL, IBASES, KBASE, KDX)
!
!*	   2.1	  FIND LOWEST CLOUD-BASE AND HIGHEST TOP.
!		  MODIFY ALL BASES TO FIRST IN-CLOUD LEVEL.
!		  MODIFY ALL TOPS TO INCLUDE INVERSION LEVEL.
!
  210 CONTINUE
!
      KLOBAS=1
      KHCTOP=NLEV
      DO 211 JL=1,NSHAL
      ITOPS(JL)=ITOPS(JL)-1
      IBASES(JL)=IBASES(JL)-1
      KLOBAS=MAX(KLOBAS,IBASES(JL))
      KHCTOP=MIN(KHCTOP,ITOPS(JL))
  211 CONTINUE
!
!		  GATHER MULTI-LEVEL ARRAYS UP TO THE SECOND LEVEL
!		  ABOVE THE ORIGINAL HIGHEST CLOUD-TOP, WHICH IS
!		  USED FOR THE MIXING LINE CALCULATIONS.
!
      DO 212 JK=KHCTOP-1,NLEV
      CALL GATHER_ir(NSHAL,ZTP1S (:,JK),PTP1 (:,JK),KDX)
      CALL GATHER_ir(NSHAL,ZQP1S (:,JK),PQP1 (:,JK),KDX)
      CALL GATHER_ir(NSHAL,ZPP1S (:,JK),APP1 (:,JK),KDX)
      CALL GATHER_ir(NSHAL,ZPHP1S(:,JK),APHP1(:,JK),KDX)
  212 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   3.	  PRELIMINARY REFERENCE PROFILE.
!		  ----------- --------- --------
!
  300 CONTINUE
!
!	   3.1	  PRESET REFERENCE PROFILE.
!
      DO 312 JK=KHCTOP,NLEV
      DO 311 JL=1,NSHAL
      ZTREF(JL,JK)=ZTP1S(JL,JK)
      ZQREF(JL,JK)=ZQP1S(JL,JK)
  311 CONTINUE
  312 CONTINUE
!
!	   3.2	  SINGLE LEVEL CALCULATIONS FOR MIXING LINE.
!
  320 CONTINUE
!
      DO 321 JL=1,NSHAL
!
!		  UPPER LEVEL T AND Q FOR MIXING LINE.
!
      ITOPM1=ITOPS(JL)-1
      ITOPP=ITOPS(JL)
      ITOPP1=ITOPS(JL)+1
      ZPTOPM1=ZPP1S(JL,ITOPM1)
      ZTTOPM1=ZTP1S(JL,ITOPM1)
      ZQTOPM1=MAX(ZEPQ,ZQP1S(JL,ITOPM1))
      ZPTOPP1=ZPP1S(JL,ITOPP1)
      ZTTOPP1=ZTP1S(JL,ITOPP1)
      ZQTOPP1=MAX(ZEPQ,ZQP1S(JL,ITOPP1))
      ZDPTOP=ZPP1S(JL,ITOPP)-ZPP1S(JL,ITOPM1)
      ZQNLM1=MAX(ZEPQ,ZQP1S(JL,NLEVM1))
!
!		  MEAN THETA AND Q FOR MIXING LINE.
!
      ZQMX=0.5*(ZQTOPM1+ZQNLM1)
      ZTHMX=0.5*(ZTTOPM1*XLG((ZPZERO/ZPTOPM1),ZCONS1)+			    &
	    ZTP1S(JL,NLEVM1)*XLG((ZPZERO/ZPP1S(JL,NLEVM1)),ZCONS1))
      ZTMX=ZTHMX*XLG((ZPTOPM1/ZPZERO),ZCONS1)
!
!		  SATURATION POINT FOR NLEVM1.
!
      ZZT=ZTP1S(JL,NLEVM1)
      ZP1=ZQNLM1*ZPP1S(JL,NLEVM1)/(ZC6+ZC7*ZQNLM1)
      ZTSP=ZC2+ZC3/(ZC4*ALOG(ZZT)-ALOG(ZP1)-ZC5)
      ZSPL=ZPP1S(JL,NLEVM1)*XLG((ZTSP/ZZT),ZC1)
!
!		  SATURATION POINT FOR UPPER MIXING LEVEL.
!
      ZP1=ZQMX*ZPTOPM1/(ZC6+ZC7*ZQMX)
      ZTSP=ZC2+ZC3/(ZC4*ALOG(ZTMX)-ALOG(ZP1)-ZC5)
      ZSPT=ZPTOPM1*XLG((ZTSP/ZTMX),ZC1)
!
!		  SLOPE OF MIXING LINE.
!
      ZDPMIX=MIN(ZSPT-ZSPL,ZEPDSP)
      ZTHBOT=ZTP1S(JL,NLEVM1)*XLG((ZPZERO/ZPP1S(JL,NLEVM1)),ZCONS1)
      ZMIX(JL)=(ZTHMX-ZTHBOT)/ZDPMIX
      ZMIX(JL)=MAX(ZEPMIX,MIN(ZMIX(JL),0.))
!
!		  SATURATION POINT FOR LEVEL ABOVE INVERSION.
!
      ZP1=ZQTOPM1*ZPTOPM1/(ZC6+ZC7*ZQTOPM1)
      ZTSP=ZC2+ZC3/(ZC4*ALOG(ZTTOPM1)-ALOG(ZP1)-ZC5)
      ZSPTI=ZPTOPM1*XLG((ZTSP/ZTTOPM1),ZC1)
!
!		  SATURATION POINT FOR LEVEL BELOW INVERSION.
!
      ZP1=ZQTOPP1*ZPTOPP1/(ZC6+ZC7*ZQTOPP1)
      ZTSP=ZC2+ZC3/(ZC4*ALOG(ZTTOPP1)-ALOG(ZP1)-ZC5)
      ZSPBI=ZPTOPP1*XLG((ZTSP/ZTTOPP1),ZC1)
!
!		  MIXING PARAMETER FOR INVERSION.
!
      ZBINV(JL)=ZBITOP*(ZSPBI-ZSPTI-ZDPTOP)				    &
	       /(ZPTOPP1-ZPTOPM1-ZDPTOP)
      ZBINV(JL)=MAX(ZBMIN,MIN(ZBMAX,ZBINV(JL)))
!
  321 CONTINUE
!
!*	   3.3	  REFERENCE TEMPERATURE AND HUMIDITY.
!
  330 CONTINUE
!
      DO 332 JK=KLOBAS,KHCTOP,-1
      DO 331 JL=1,NSHAL
      IF (JK.EQ.ITOPS(JL)) THEN
	 ZBETA=ZBINV(JL)
      ELSE
	 ZBETA=ZBSHAL
      ENDIF

      if(zepmix.gt.min(zmix(jl),0.)) zbeta=0.0 ! beta=0 gives DALR

      LOA=(JK.GE.ITOPS(JL)).AND.(JK.LT.IBASES(JL))
      IF (LOA) THEN
            ZTREF(JL,JK)=                                                   &
                 (ZTREF(JL,JK+1)*XLG((ZPZERO/ZPP1S(JL,JK+1)),ZCONS1)+	    &
                 ZBETA*ZSTABM*ZMIX(JL)*(ZPP1S(JL,JK)-ZPP1S(JL,JK+1)))*	    &
                 XLG((ZPP1S(JL,JK)/ZPZERO),ZCONS1)
      ENDIF
      LOB=(JK.GE.ITOPS(JL)).AND.(JK.LE.IBASES(JL))
      IF (LOB) THEN 
	 ZSP=ZSPDIFF-(ZBETA-1.)*(ZPP1S(JL,JK+1)-ZPP1S(JL,JK))		    &
	     +ZPP1S(JL,JK)
	 ZTS=ZTREF(JL,JK)*XLG((ZSP/ZPP1S(JL,JK)),ZCONS1)
	 ZQREF(JL,JK)=C2ES*EXP(C3LES*(ZTS-TMELT)/(ZTS-C4LES))/ZSP
      ENDIF
  331 CONTINUE
  332 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   4.	  CONSIDERATION OF ENERGY CONSERVATION.
!		  ------------- -- ------ -------------
!
  400 CONTINUE
!
!*	   4.1	  SEPERATE ENERGY INTEGRALS FOR T AND Q.
!
  410 CONTINUE
!
      DO 411 JL=1,NSHAL
      ZTSUM(JL)=0.
      ZQSUM(JL)=0.
      ZQINT(JL)=0.
  411 CONTINUE
!
      DO 413 JK=KHCTOP,KLOBAS
      DO 412 JL=1,NSHAL
      ZDP(JL,JK)=ZPHP1S(JL,JK+1)-ZPHP1S(JL,JK)
      ZTSUM(JL)=ZTSUM(JL)+(ZTP1S(JL,JK)-ZTREF(JL,JK))*ZDP(JL,JK)
      ZQSUM(JL)=ZQSUM(JL)+(ZQP1S(JL,JK)-ZQREF(JL,JK))*ZDP(JL,JK)
  412 CONTINUE
  413 CONTINUE
!
      DO 414 JL=1,NSHAL
      IBOT=IBASES(JL)+1
      ITOPP=ITOPS(JL)
      ZTSUM(JL)=ZTSUM(JL)/(ZPHP1S(JL,IBOT)-ZPHP1S(JL,ITOPP))
      ZQSUM(JL)=ZQSUM(JL)/(ZPHP1S(JL,IBOT)-ZPHP1S(JL,ITOPP))
  414 CONTINUE
!
!*	   4.2	  MODIFICATION OF PRELIMINARY REFERENCE PROFILE.
!		  COMPUTE INTEGRATED CONDENSATION RATE.
!
  420 CONTINUE
!
      DO 422 JK=KHCTOP,KLOBAS
      DO 421 JL=1,NSHAL
      LO=(JK.GE.ITOPS(JL)).AND.(JK.LE.IBASES(JL))
      IF (LO) THEN
	 ZTREF(JL,JK)=ZTREF(JL,JK)+ZTSUM(JL)
	 ZQREF(JL,JK)=ZQREF(JL,JK)+ZQSUM(JL)
	 ZCOND=(ZQREF(JL,JK)-ZQP1S(JL,JK))*ZCONS3
	 IF (ZCOND.GT.0.) THEN
	    ZQINT(JL)=ZQINT(JL)+ZCOND*ZDP(JL,JK)/G
	 ENDIF
      ENDIF
  421 CONTINUE
  422 CONTINUE
!
!     ------------------------------------------------------------------
!
!*	   5.	  FINAL TENDENCIES, SCATTER TO FULL GRID.
!		  ----- ----------- ------- -- ---- -----
!
  500 CONTINUE
!

      DO 502 JK=KHCTOP,KLOBAS
!
      DO 501 JL=1,NSHAL
      ZDQS(JL)=(ZQREF(JL,JK)-ZQP1S(JL,JK))*ZCONS3
      ZDTS(JL)=(ZTREF(JL,JK)-ZTP1S(JL,JK))*ZCONS3
  501 CONTINUE
!
      CALL SCATTER(NSHAL,PDT(:,JK),KDX,ZDTS)
      CALL SCATTER(NSHAL,PDQ(:,JK),KDX,ZDQS)
      call scatter(nshal,sztref(:,jk),kdx,ztref(:,jk))
!
  502 CONTINUE
!
      CALL SCATTER(NSHAL,PSC,KDX,ZQINT)
!
!     ------------------------------------------------------------------
!
!*	   6.	  RETURN WORKSPACE.
!		  ------ ----------
!
  600 CONTINUE
!
!
      RETURN
    END SUBROUTINE BMSHAL

!----------------------------------------------------------------------------!
! Following were added by cw so code can be run with portland group compiler !
!----------------------------------------------------------------------------!


! Search an integer vector and return the number of values equal to
! the target and return an integer vector containing their locations.
     subroutine wheneq( n, x, incx, ftarget, indx, nn )
	 integer,		intent(in)  :: n, incx
	 integer,		intent(in)  :: ftarget
	 integer, dimension(:), intent(in)  :: x
	 integer, dimension(:), intent(out) :: indx
	 integer,		intent(out) :: nn
	 integer			    :: i, first, last, lloc

	 nn = 0

!	  if( n .gt. 0 .and. incx .ne. 0 ) then
         last = size(x)
         
         ! Begin search
         lloc = 0
         do i = 1, last, incx
            lloc = lloc + 1
            if( x(i) == ftarget ) then
               nn = nn + 1
               indx(nn) = lloc
            endif
         enddo
         !	  endif

       end subroutine wheneq
   
!-----------------------------------------------------------------------------
! rewrite these as an interface!
       subroutine gather_ir(nsmall, small, large, ind_large)
         
         real small(:), large(:)
         integer nsmall, ind_large(:), i
         
         do i=1, nsmall
            small(i) = large(ind_large(i))
         enddo
         
         return
       end subroutine gather_ir

!-------------------------------

       subroutine gather_ii(nsmall, small, large, ind_large)
         
         integer small(:), large(:)
         integer nsmall, ind_large(:), i
         
         do i=1, nsmall
            small(i) = large(ind_large(i))
         enddo
         
         return
       end subroutine gather_ii

!-----------------------------------------------------------------------------

      subroutine scatter(nsmall, large, ind_large, small)

      real small(:), large(:)
      integer nsmall, ind_large(:), i

      do i=1, nsmall
	 large(ind_large(i)) = small(i)
      enddo

      return
    end subroutine scatter

  end module betts_miller_mod
