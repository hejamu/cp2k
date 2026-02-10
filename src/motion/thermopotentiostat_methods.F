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

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger, &
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file, &
                                              cp_print_key_finished_output, &
                                              cp_print_key_should_output, &
                                              cp_print_key_unit_nr
   USE cp_result_methods,               ONLY: get_results, &
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get, &
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE force_env_types,                 ONLY: force_env_get, &
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals, &
                                              section_vals_type, &
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length, &
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi
   USE md_environment_types,            ONLY: get_md_env, &
                                              md_environment_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_type
   USE molecule_types,                  ONLY: global_constraint_type, &
                                              molecule_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: debye
   USE qs_environment_types,            ONLY: get_qs_env, &
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set, &
                                              qs_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE simpar_types,                    ONLY: simpar_type
   USE external_potential_types,        ONLY: get_potential, &
                                              gth_potential_type, &
                                              set_potential
   USE thermopotentiostat_types,        ONLY: thermopotentiostat_type, &
                                              thermopotentiostat_create, &
                                              thermopotentiostat_release

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermopotentiostat_methods'

   PUBLIC :: thermopotentiostat_init, &
             thermopotentiostat_update, &
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

      CHARACTER(LEN=default_string_length)               :: element_symbol
      INTEGER                                            :: direction, ikind, nkind
      LOGICAL                                            :: use_integral_control, element_found
      REAL(dp)                                           :: coupling_strength, integral_gain, &
                                                            max_zeff_correction, min_zeff_correction, &
                                                            target_dipole, zeff_correction
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: thermo_section

      NULLIFY (thermo_section, qs_kind_set, gth_potential, dft_control)

      thermo_section => section_vals_get_subs_vals(md_section, "THERMOPOTENTIOSTAT")

      ! Read input parameters
      CALL section_vals_val_get(thermo_section, "TARGET_DIPOLE", r_val=target_dipole)
      CALL section_vals_val_get(thermo_section, "DIPOLE_DIRECTION", i_val=direction)
      CALL section_vals_val_get(thermo_section, "COUPLING_STRENGTH", r_val=coupling_strength)
      CALL section_vals_val_get(thermo_section, "ELEMENT", c_val=element_symbol)
      CALL section_vals_val_get(thermo_section, "MAX_ZEFF_CORRECTION", r_val=max_zeff_correction)
      CALL section_vals_val_get(thermo_section, "MIN_ZEFF_CORRECTION", r_val=min_zeff_correction)
      CALL section_vals_val_get(thermo_section, "USE_INTEGRAL_CONTROL", l_val=use_integral_control)
      CALL section_vals_val_get(thermo_section, "INTEGRAL_GAIN", r_val=integral_gain)

      ! Check that element was specified
      IF (LEN_TRIM(element_symbol) == 0) THEN
         CPABORT("THERMOPOTENTIOSTAT requires ELEMENT to be specified")
      END IF

      ! Create the thermopotentiostat structure
      CALL thermopotentiostat_create(thermo_pot, target_dipole, coupling_strength, &
                                     direction, element_symbol, max_zeff_correction, &
                                     min_zeff_correction, use_integral_control, integral_gain)

      ! Find the atom kind index for the specified element
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      nkind = SIZE(qs_kind_set)

      ! Enable CORE_CORR_DIP automatically for thermopotentiostat
      ! This is required for electron density to respond to zeff_correction changes
      ! Note: SURFACE_DIPOLE_CORRECTION is already validated in md_run.F
      IF (.NOT. dft_control%correct_el_density_dip) THEN
         dft_control%correct_el_density_dip = .TRUE.
         CPWARN("THERMOPOTENTIOSTAT: Auto-enabling CORE_CORR_DIP")
      END IF

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

   END SUBROUTINE thermopotentiostat_init

! **************************************************************************************************
!> \brief Update the thermopotentiostat based on current dipole moment
!>        This is called at each MD step to adjust zeff_correction
!> \param thermo_pot the thermopotentiostat structure
!> \param qs_env the QS environment
!> \param dt the timestep
!> \param dipole_moment output: the current dipole moment
!> \param dipole_available output: whether dipole was available
! **************************************************************************************************
   SUBROUTINE thermopotentiostat_update(thermo_pot, qs_env, dt, dipole_moment, dipole_available)

      TYPE(thermopotentiostat_type), INTENT(INOUT)       :: thermo_pot
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(IN)                               :: dt
      REAL(dp), DIMENSION(3), INTENT(OUT)                :: dipole_moment
      LOGICAL, INTENT(OUT)                               :: dipole_available

      CHARACTER(LEN=*), PARAMETER :: routineN = 'thermopotentiostat_update'
      CHARACTER(LEN=default_string_length)               :: description

      INTEGER                                            :: handle
      REAL(dp)                                           :: delta_zeff, dipole_error, &
                                                            new_zeff_corr, proportional_term, &
                                                            integral_term
      TYPE(cp_result_type), POINTER                      :: results

      CALL timeset(routineN, handle)

      NULLIFY (results)
      dipole_moment = 0.0_dp
      dipole_available = .FALSE.

      ! Get the dipole moment from results
      CALL get_qs_env(qs_env, results=results)
      description = "[DIPOLE]"

      IF (test_for_result(results, description=description)) THEN
         CALL get_results(results, description=description, values=dipole_moment)
         dipole_available = .TRUE.

         ! Store last dipole for output
         thermo_pot%last_dipole = dipole_moment

         ! Calculate the error (deviation from target)
         dipole_error = dipole_moment(thermo_pot%direction) - thermo_pot%target_dipole

         ! Proportional control term
         proportional_term = thermo_pot%coupling_strength * dipole_error

         ! Integral control term (if enabled)
         integral_term = 0.0_dp
         IF (thermo_pot%use_integral_control) THEN
            thermo_pot%integral_error = thermo_pot%integral_error + dipole_error * dt
            integral_term = thermo_pot%integral_gain * thermo_pot%integral_error
         END IF

         ! Total change in zeff_correction
         delta_zeff = (proportional_term + integral_term) * dt

         ! Update zeff_correction with limits
         new_zeff_corr = thermo_pot%current_zeff_corr + delta_zeff
         new_zeff_corr = MIN(new_zeff_corr, thermo_pot%max_zeff_correction)
         new_zeff_corr = MAX(new_zeff_corr, thermo_pot%min_zeff_correction)

         thermo_pot%current_zeff_corr = new_zeff_corr

      ELSE
         CPWARN("THERMOPOTENTIOSTAT: Dipole not available. Enable MOMENTS")
      END IF

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
                                                            new_zeff, total_zeff_corr, zeff_corr
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, gth_potential, elec_conf)

      IF (.NOT. thermo_pot%active) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ikind = thermo_pot%atom_kind_index
      zeff_corr = thermo_pot%current_zeff_corr

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)

      IF (ASSOCIATED(qs_kind_set(ikind)%gth_potential)) THEN
         gth_potential => qs_kind_set(ikind)%gth_potential

         ! Get current parameters
         CALL get_potential(potential=gth_potential, &
                            elec_conf=elec_conf, &
                            alpha_core_charge=alpha)

         ! Calculate new zeff and derived quantities
         new_zeff = REAL(SUM(elec_conf), dp) + zeff_corr

         ! Update the core charge coefficient: ccore_charge = zeff * sqrt((alpha/pi)^3)
         ccore_charge = new_zeff * SQRT((alpha/pi)**3)
         cerf_ppl = ccore_charge  ! Same formula for error function coefficient

         ! Set the updated values on the potential
         CALL set_potential(gth_potential, &
                            zeff=new_zeff, &
                            zeff_correction=zeff_corr, &
                            ccore_charge=ccore_charge, &
                            cerf_ppl=cerf_ppl)

      END IF

      ! Update total_zeff_corr in qs_env for electron density adjustment
      ! This is used by qs_mo_occupation to adjust the electron count
      ! See qs_mo_occupation.F lines ~465-520
      CALL get_qs_kind_set(qs_kind_set, total_zeff_corr=total_zeff_corr)
      qs_env%total_zeff_corr = total_zeff_corr

      CALL timestop(handle)

   END SUBROUTINE update_zeff_correction_for_dipole

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
            WRITE (iw, '("#",5X,A,10X,A,8X,A,8X,A,8X,A,8X,A)') &
               "Step Nr.", "Time[fs]", "Dipole_X[a.u.]", "Dipole_Y[a.u.]", &
               "Dipole_Z[a.u.]", "Zeff_Corr[a.u.]"
         END IF
         WRITE (iw, "(I10,F20.6,F20.9,F20.9,F20.9,F20.9)") &
            itimes, time, thermopotentiostat%last_dipole(1), &
            thermopotentiostat%last_dipole(2), thermopotentiostat%last_dipole(3), &
            thermopotentiostat%current_zeff_corr
         CALL m_flush(iw)
      END IF
      CALL cp_print_key_finished_output(iw, logger, motion_section, "MD%PRINT%THERMOPOTENTIOSTAT")

      ! Print summary to main output (PROGRAM_RUN_INFO)
      iw = cp_print_key_unit_nr(logger, motion_section, "MD%PRINT%PROGRAM_RUN_INFO", &
                                extension=".mdLog", file_position=my_pos, file_action=my_act)
      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') "THERMOPOTENTIOSTAT| Status:"
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Current dipole (target dir.)  [a.u.]:", &
            thermopotentiostat%last_dipole(thermopotentiostat%direction)
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Target dipole                 [a.u.]:", &
            thermopotentiostat%target_dipole
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Dipole error                  [a.u.]:", &
            thermopotentiostat%last_dipole(thermopotentiostat%direction) - thermopotentiostat%target_dipole
         WRITE (iw, '(T2,A,T61,F20.9)') "THERMOPOTENTIOSTAT| Current zeff correction       [a.u.]:", &
            thermopotentiostat%current_zeff_corr
      END IF
      CALL cp_print_key_finished_output(iw, logger, motion_section, "MD%PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE print_thermopotentiostat_status

END MODULE thermopotentiostat_methods
