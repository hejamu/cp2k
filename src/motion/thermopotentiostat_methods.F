!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for the Thermopotentiostat ensemble
!>        The thermopotentiostat varies the core correction (zeff_correction) of a selected
!>        atom type based on the system dipole moment, allowing control of electrochemical
!>        potential in surface simulations.
!>
!>        This follows the approach described in:
!>        SS, MT, MWF, JN PRL, 2018, 120, 246801
!>        where the electron density is adjusted according to the applied core correction.
!>
!> \par History
!>      02.2026 created [thermopotentiostat implementation]
!> \author CP2K developers
! **************************************************************************************************
MODULE thermopotentiostat_methods

   USE atomic_kind_types,               ONLY: atomic_kind_type, &
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output, &
                                              cp_print_key_unit_nr
   USE external_potential_types,        ONLY: get_potential, &
                                              gth_potential_type, &
                                              set_potential
   USE input_section_types,             ONLY: section_vals_get_subs_vals, &
                                              section_vals_type, &
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length, &
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi
   USE parallel_rng_types,              ONLY: GAUSSIAN, &
                                              next_rng_seed, &
                                              rng_stream_type
   USE physcon,                         ONLY: bohr, &
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env, &
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set, &
                                              qs_kind_type
   USE simpar_types,                    ONLY: simpar_type
   USE thermopotentiostat_types,        ONLY: thermopotentiostat_type, &
                                              thermopotentiostat_create

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermopotentiostat_methods'

   PUBLIC :: thermopotentiostat_init, &
             thermopotentiostat_update, &
             thermopotentiostat_compute_potential, &
             update_zeff_correction_for_dipole, &
             print_thermopotentiostat_status

CONTAINS

! **************************************************************************************************
!> \brief Initialize the thermopotentiostat from input
!> \param thermo_pot the thermopotentiostat structure to initialize
!> \param md_section the MD input section
!> \param qs_env the QS environment (to find atom kind)
! **************************************************************************************************
   SUBROUTINE thermopotentiostat_init(thermo_pot, md_section, qs_env)

      TYPE(thermopotentiostat_type), POINTER             :: thermo_pot
      TYPE(section_vals_type), POINTER                   :: md_section
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=40)                                  :: rng_name
      CHARACTER(LEN=default_string_length)               :: element_symbol
      INTEGER                                            :: direction, ikind, nkind
      LOGICAL                                            :: element_found
      REAL(dp)                                           :: electrode_charge, geometric_capacitance, &
                                                            target_potential, timecon, zeff_correction
      REAL(dp), DIMENSION(3, 2)                          :: initial_seed
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: thermo_section

      NULLIFY (thermo_section, qs_kind_set, gth_potential, dft_control)

      thermo_section => section_vals_get_subs_vals(md_section, "THERMOPOTENTIOSTAT")

      ! Read input parameters
      CALL section_vals_val_get(thermo_section, "TARGET_POTENTIAL", r_val=target_potential)
      CALL section_vals_val_get(thermo_section, "DIPOLE_DIRECTION", i_val=direction)
      CALL section_vals_val_get(thermo_section, "GEOMETRIC_CAPACITANCE", r_val=geometric_capacitance)
      CALL section_vals_val_get(thermo_section, "TIMECON", r_val=timecon)
      CALL section_vals_val_get(thermo_section, "ELEMENT", c_val=element_symbol)
      CALL section_vals_val_get(thermo_section, "ELECTRODE_CHARGE", r_val=electrode_charge)

      ! Validate input parameters
      IF (LEN_TRIM(element_symbol) == 0) THEN
         CPABORT("THERMOPOTENTIOSTAT requires ELEMENT to be specified")
      END IF
      IF (geometric_capacitance <= 0.0_dp) THEN
         CPABORT("THERMOPOTENTIOSTAT: GEOMETRIC_CAPACITANCE must be positive")
      END IF
      IF (timecon < 0.0_dp) THEN
         CPABORT("THERMOPOTENTIOSTAT: TIMECON must be non-negative")
      END IF

      ! Create the thermopotentiostat structure
      CALL thermopotentiostat_create(thermo_pot, target_potential, geometric_capacitance, &
                                     timecon, direction, element_symbol)

      ! Initialize random number generator for stochastic dynamics
      initial_seed = next_rng_seed()
      WRITE (UNIT=rng_name, FMT="(A)") "Thermopotentiostat RNG"
      thermo_pot%rng_stream = rng_stream_type( &
                              name=rng_name, distribution_type=GAUSSIAN, &
                              extended_precision=.TRUE., seed=initial_seed)

      ! Find the atom kind index for the specified element
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      nkind = SIZE(qs_kind_set)

      ! Verify required DFT settings are active
      IF (.NOT. dft_control%correct_surf_dip) &
         CPABORT("THERMOPOTENTIOSTAT requires SURFACE_DIPOLE_CORRECTION")
      IF (.NOT. dft_control%correct_el_density_dip) &
         CPABORT("THERMOPOTENTIOSTAT requires CORE_CORR_DIP")

      element_found = .FALSE.
      DO ikind = 1, nkind
         IF (TRIM(qs_kind_set(ikind)%element_symbol) == TRIM(element_symbol)) THEN
            element_found = .TRUE.
            thermo_pot%atom_kind_index = ikind
            ! Check that GTH potential is available
            IF (.NOT. ASSOCIATED(qs_kind_set(ikind)%gth_potential)) THEN
               CALL cp_abort(__LOCATION__, &
                             "THERMOPOTENTIOSTAT: Element "//TRIM(element_symbol)// &
                             " must have a GTH potential for zeff_correction to work")
            END IF
            ! Get initial zeff_correction
            CALL get_potential(potential=qs_kind_set(ikind)%gth_potential, &
                               zeff_correction=zeff_correction)
            thermo_pot%initial_zeff_corr = zeff_correction
            thermo_pot%current_zeff_corr = zeff_correction
            EXIT
         END IF
      END DO

      IF (.NOT. element_found) THEN
         CALL cp_abort(__LOCATION__, &
                       "THERMOPOTENTIOSTAT: Element "//TRIM(element_symbol)//" not found in system")
      END IF

      ! On restart, override current_zeff_corr with saved electrode charge
      ! and apply the correction to the GTH potential before the first force calculation
      IF (ABS(electrode_charge) > 0.0_dp) THEN
         thermo_pot%current_zeff_corr = electrode_charge
         CPWARN("THERMOPOTENTIOSTAT: Restoring electrode charge from restart file")
         CALL update_zeff_correction_for_dipole(thermo_pot, qs_env)
      END IF

   END SUBROUTINE thermopotentiostat_init

! **************************************************************************************************
!> \brief Update the thermopotentiostat based on current dipole moment
!>        This is called at each MD step to adjust zeff_correction
!>        The potential is calculated from the dipole: V = -4*pi*dipole/area
!>        Dissipative term: C * (V - V_target) * (1 - exp(-dt/tau))
!>        Fluctuation term: N * sqrt(k_B * T * C * (1 - exp(-2*dt/tau)))
!> \param thermo_pot the thermopotentiostat structure
!> \param qs_env the QS environment
!> \param simpar the simulation parameters (timestep, temperature, etc.)
! **************************************************************************************************
   SUBROUTINE thermopotentiostat_update(thermo_pot, qs_env, simpar)

      TYPE(thermopotentiostat_type), INTENT(INOUT)       :: thermo_pot
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(simpar_type), INTENT(IN)                      :: simpar

      CHARACTER(LEN=*), PARAMETER :: routineN = 'thermopotentiostat_update'

      INTEGER                                            :: handle
      REAL(dp)                                           :: current_potential, dt, &
                                                            dissipative_term, exp_factor, &
                                                            exp_factor_2, fluctuation_term, &
                                                            new_zeff_corr

      CALL timeset(routineN, handle)

      dt = simpar%dt

      ! Compute current potential from surface dipole and store in thermo_pot
      CALL thermopotentiostat_compute_potential(thermo_pot, qs_env)
      current_potential = thermo_pot%last_current_potential

      ! Calculate exponential factors for Langevin dynamics
      ! exp_factor = 1 - exp(-dt/tau)
      ! exp_factor_2 = 1 - exp(-2*dt/tau)
      IF (thermo_pot%timecon > 0.0_dp) THEN
         exp_factor = 1.0_dp - EXP(-dt / thermo_pot%timecon)
         exp_factor_2 = 1.0_dp - EXP(-2.0_dp * dt / thermo_pot%timecon)
      ELSE
         ! Instantaneous response if timecon is zero
         exp_factor = 1.0_dp
         exp_factor_2 = 1.0_dp
      END IF

      ! Dissipative term: C * (V_current - V_target) * (1 - exp(-dt/tau))
      dissipative_term = - 1.0_dp * thermo_pot%geometric_capacitance * &
                           (current_potential - thermo_pot%target_potential) * exp_factor

      ! Fluctuation term: N * sqrt(k_B * T * C * (1 - exp(-2*dt/tau)))
      ! simpar%temp_ext is already k_B*T in atomic units (Hartree)
      fluctuation_term = thermo_pot%rng_stream%next() * &
                           SQRT(simpar%temp_ext * thermo_pot%geometric_capacitance * exp_factor_2)

      ! Store Langevin terms for output (potential already stored by thermopotentiostat_compute_potential)
      thermo_pot%last_dissipative_term = dissipative_term
      thermo_pot%last_fluctuation_term = fluctuation_term

      ! Update zeff_correction with dissipative and fluctuation contributions
      new_zeff_corr = thermo_pot%current_zeff_corr + dissipative_term + fluctuation_term

      thermo_pot%current_zeff_corr = new_zeff_corr

      CALL timestop(handle)

   END SUBROUTINE thermopotentiostat_update

! **************************************************************************************************
!> \brief Apply the updated zeff_correction to the specified atom kind
!>        This modifies the GTH potential and updates qs_env%total_zeff_corr
!>        The electron density will be adjusted according to the applied core correction
!>        following the approach in PRL 120, 246801 (2018)
!> \param thermo_pot the thermopotentiostat structure
!> \param qs_env the QS environment to modify
! **************************************************************************************************
   SUBROUTINE update_zeff_correction_for_dipole(thermo_pot, qs_env)

      TYPE(thermopotentiostat_type), INTENT(IN)          :: thermo_pot
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'update_zeff_correction_for_dipole'

      INTEGER                                            :: handle, ikind, natom_kind
      INTEGER, DIMENSION(:), POINTER                     :: elec_conf
      REAL(dp)                                           :: alpha, ccore_charge, cerf_ppl, &
                                                            new_zeff, total_zeff_corr, zeff_corr, &
                                                            zeff_corr_per_atom
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, gth_potential, elec_conf, atomic_kind_set)

      IF (.NOT. thermo_pot%active) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ikind = thermo_pot%atom_kind_index
      zeff_corr = thermo_pot%current_zeff_corr

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      ! Validate atom_kind_index is within bounds
      IF (ikind < 1 .OR. ikind > SIZE(qs_kind_set)) THEN
         CPABORT("THERMOPOTENTIOSTAT: atom_kind_index out of bounds")
      END IF

      ! Get the number of atoms of this kind
      CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_kind)

      ! Distribute the total correction across all atoms of this kind
      ! The total system charge correction must be divided by the number of electrode atoms
      IF (natom_kind > 0) THEN
         zeff_corr_per_atom = zeff_corr / REAL(natom_kind, dp)
      ELSE
         CPABORT("THERMOPOTENTIOSTAT: No atoms found for the specified element")
      END IF

      IF (.NOT. ASSOCIATED(qs_kind_set(ikind)%gth_potential)) THEN
         CALL cp_abort(__LOCATION__, &
                       "THERMOPOTENTIOSTAT: GTH potential for element "// &
                       TRIM(thermo_pot%element_symbol)//" is no longer associated")
      END IF

      gth_potential => qs_kind_set(ikind)%gth_potential

      ! Get current parameters
      CALL get_potential(potential=gth_potential, &
                         elec_conf=elec_conf, &
                         alpha_core_charge=alpha)

      ! Calculate new zeff and derived quantities using per-atom correction
      new_zeff = REAL(SUM(elec_conf), dp) + zeff_corr_per_atom

      ! Update the core charge coefficient: ccore_charge = zeff * sqrt((alpha/pi)^3)
      ccore_charge = new_zeff * SQRT((alpha/pi)**3)
      cerf_ppl = ccore_charge  ! Same formula for error function coefficient

      ! Set the updated values on the potential
      CALL set_potential(gth_potential, &
                         zeff=new_zeff, &
                         zeff_correction=zeff_corr_per_atom, &
                         ccore_charge=ccore_charge, &
                         cerf_ppl=cerf_ppl)

      ! Update total_zeff_corr in qs_env for electron density adjustment
      ! This is used by qs_mo_occupation to adjust the electron count
      ! See qs_mo_occupation.F lines ~465-520
      CALL get_qs_kind_set(qs_kind_set, total_zeff_corr=total_zeff_corr)
      qs_env%total_zeff_corr = total_zeff_corr

      CALL timestop(handle)

   END SUBROUTINE update_zeff_correction_for_dipole

! **************************************************************************************************
!> \brief Compute and store the current potential from the surface dipole moment
!>        This is used to populate the thermopotentiostat output values at step 0,
!>        before the first Langevin update. Only computes and stores the potential;
!>        does NOT perform the Langevin dynamics (dissipative/fluctuation terms remain zero).
!> \param thermo_pot the thermopotentiostat structure
!> \param qs_env the QS environment
! **************************************************************************************************
   SUBROUTINE thermopotentiostat_compute_potential(thermo_pot, qs_env)

      TYPE(thermopotentiostat_type), INTENT(INOUT)       :: thermo_pot
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: dir
      REAL(dp)                                           :: area, dipole_moment, current_potential
      TYPE(cell_type), POINTER                           :: cell

      IF (.NOT. thermo_pot%active) RETURN

      NULLIFY (cell)
      CALL get_qs_env(qs_env, cell=cell)

      ! Calculate surface area perpendicular to dipole direction
      dir = thermo_pot%direction
      SELECT CASE (dir)
      CASE (1)  ! X direction: area = Y * Z
         area = cell%hmat(2, 2) * cell%hmat(3, 3)
      CASE (2)  ! Y direction: area = X * Z
         area = cell%hmat(1, 1) * cell%hmat(3, 3)
      CASE (3)  ! Z direction: area = X * Y
         area = cell%hmat(1, 1) * cell%hmat(2, 2)
      CASE DEFAULT
         CPABORT("THERMOPOTENTIOSTAT: Invalid DIPOLE_DIRECTION")
      END SELECT

      ! Get surface dipole moment
      dipole_moment = qs_env%surface_dipole_moment * bohr
      thermo_pot%last_dipole = dipole_moment

      ! Calculate current potential from dipole: V = -4*pi*dipole / area
      current_potential = - 4.0_dp * pi * dipole_moment / area
      thermo_pot%last_current_potential = current_potential

   END SUBROUTINE thermopotentiostat_compute_potential

! **************************************************************************************************
!> \brief Print thermopotentiostat status to output file and dedicated file
!> \param thermopotentiostat the thermopotentiostat structure
!> \param itimes current timestep
!> \param time current simulation time in fs
!> \param motion_section the MOTION input section for print keys
!> \param logger the CP2K logger
!> \param my_pos file position (APPEND/REWIND)
!> \param my_act file action (WRITE)
! **************************************************************************************************
   SUBROUTINE print_thermopotentiostat_status(thermopotentiostat, itimes, time, &
                                              motion_section, logger, my_pos, my_act)

      TYPE(thermopotentiostat_type), POINTER             :: thermopotentiostat
      INTEGER, POINTER                                   :: itimes
      REAL(dp), INTENT(IN)                               :: time
      TYPE(section_vals_type), POINTER                   :: motion_section
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(LEN=*), INTENT(IN)                       :: my_pos, my_act

      INTEGER                                            :: iw
      LOGICAL                                            :: new_file

      IF (.NOT. ASSOCIATED(thermopotentiostat)) RETURN
      IF (.NOT. thermopotentiostat%active) RETURN

      ! Print to dedicated thermopotentiostat file
      iw = cp_print_key_unit_nr(logger, motion_section, "MD%PRINT%THERMOPOTENTIOSTAT", &
                                extension=".thermopotentiostat", file_position=my_pos, &
                                file_action=my_act, is_new_file=new_file)
      IF (iw > 0) THEN
         IF (new_file) THEN
            ! Write header
            WRITE (iw, '("#",5X,A,8X,A,8X,A,7X,A,8X,A,12X,A)') &
               "Step Nr.", "Time[fs]", "Potential [V]", "Charge [a.u.]", &
               "Disp [a.u.]", "Fluc [a.u.]"
         END IF
         WRITE (iw, "(I10,F20.6,F20.9,F20.9,F20.9,F20.9)") &
            itimes, time, thermopotentiostat%last_current_potential * evolt, &
            thermopotentiostat%current_zeff_corr, &
            thermopotentiostat%last_dissipative_term, &
            thermopotentiostat%last_fluctuation_term
         CALL m_flush(iw)
      END IF
      CALL cp_print_key_finished_output(iw, logger, motion_section, "MD%PRINT%THERMOPOTENTIOSTAT")

      ! Print summary to main output (PROGRAM_RUN_INFO)
      iw = cp_print_key_unit_nr(logger, motion_section, "MD%PRINT%PROGRAM_RUN_INFO", &
                                extension=".mdLog", file_position=my_pos, file_action=my_act)
      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') "THERMOPOTENTIOSTAT| Status:"
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Current potential                [V]:", &
            thermopotentiostat%last_current_potential * evolt
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Target potential                 [V]:", &
            thermopotentiostat%target_potential * evolt
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Potential error                  [V]:", &
            (thermopotentiostat%last_current_potential - thermopotentiostat%target_potential) * evolt
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Current zeff correction       [a.u.]:", &
            thermopotentiostat%current_zeff_corr
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Dissipative dN                [a.u.]:", &
            thermopotentiostat%last_dissipative_term
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Fluctuation dN                [a.u.]:", &
            thermopotentiostat%last_fluctuation_term
      END IF
      CALL cp_print_key_finished_output(iw, logger, motion_section, "MD%PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE print_thermopotentiostat_status

END MODULE thermopotentiostat_methods
