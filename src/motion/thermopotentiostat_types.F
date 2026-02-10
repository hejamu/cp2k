!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for the Thermopotentiostat ensemble
!>        The thermopotentiostat varies the core correction (zeff_correction) of a selected
!>        atom type based on the system dipole moment, allowing control of electrochemical
!>        potential in surface simulations.
!> \par History
!>      02.2026 created [thermopotentiostat implementation]
!> \author CP2K developers
! **************************************************************************************************
MODULE thermopotentiostat_types

   USE kinds,                           ONLY: dp, &
                                              default_string_length

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermopotentiostat_types'

   PUBLIC :: thermopotentiostat_type, &
             thermopotentiostat_create, &
             thermopotentiostat_release

! **************************************************************************************************
!> \brief Type to store thermopotentiostat parameters and state
! **************************************************************************************************
   TYPE thermopotentiostat_type
      ! Target dipole moment (atomic units)
      REAL(dp)                              :: target_dipole = 0.0_dp
      ! Coupling strength / feedback constant (controls response rate)
      REAL(dp)                              :: coupling_strength = 0.1_dp
      ! Current zeff correction value being applied
      REAL(dp)                              :: current_zeff_corr = 0.0_dp
      ! Maximum allowed zeff correction magnitude
      REAL(dp)                              :: max_zeff_correction = 1.0_dp
      ! Minimum allowed zeff correction magnitude
      REAL(dp)                              :: min_zeff_correction = -1.0_dp
      ! Direction for dipole control: 1=x, 2=y, 3=z
      INTEGER                               :: direction = 3
      ! Index of the atom kind to modify
      INTEGER                               :: atom_kind_index = -1
      ! Element symbol of the atom kind to modify
      CHARACTER(LEN=default_string_length)  :: element_symbol = ""
      ! Whether to use integral (PI) control in addition to proportional
      LOGICAL                               :: use_integral_control = .FALSE.
      ! Integral term accumulator
      REAL(dp)                              :: integral_error = 0.0_dp
      ! Integral gain (for PI control)
      REAL(dp)                              :: integral_gain = 0.01_dp
      ! Store the initial zeff_correction from input
      REAL(dp)                              :: initial_zeff_corr = 0.0_dp
      ! Last measured dipole moment for output
      REAL(dp), DIMENSION(3)                :: last_dipole = [0.0_dp, 0.0_dp, 0.0_dp]
      ! Flag indicating if thermopotentiostat is active
      LOGICAL                               :: active = .FALSE.
   END TYPE thermopotentiostat_type

CONTAINS

! **************************************************************************************************
!> \brief Create and initialize a thermopotentiostat structure
!> \param thermo_pot the thermopotentiostat to create
!> \param target_dipole target dipole moment
!> \param coupling_strength feedback coupling constant
!> \param direction dipole direction (1=x, 2=y, 3=z)
!> \param element_symbol element symbol of atom kind to modify
!> \param max_zeff_correction maximum allowed zeff correction
!> \param min_zeff_correction minimum allowed zeff correction
!> \param use_integral_control whether to use PI control
!> \param integral_gain integral term gain for PI control
! **************************************************************************************************
   SUBROUTINE thermopotentiostat_create(thermo_pot, target_dipole, coupling_strength, &
                                        direction, element_symbol, max_zeff_correction, &
                                        min_zeff_correction, use_integral_control, integral_gain)

      TYPE(thermopotentiostat_type), POINTER             :: thermo_pot
      REAL(dp), INTENT(IN)                               :: target_dipole
      REAL(dp), INTENT(IN)                               :: coupling_strength
      INTEGER, INTENT(IN)                                :: direction
      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol
      REAL(dp), INTENT(IN)                               :: max_zeff_correction
      REAL(dp), INTENT(IN)                               :: min_zeff_correction
      LOGICAL, INTENT(IN)                                :: use_integral_control
      REAL(dp), INTENT(IN)                               :: integral_gain

      CPASSERT(.NOT. ASSOCIATED(thermo_pot))
      ALLOCATE (thermo_pot)

      thermo_pot%target_dipole = target_dipole
      thermo_pot%coupling_strength = coupling_strength
      thermo_pot%direction = direction
      thermo_pot%element_symbol = TRIM(element_symbol)
      thermo_pot%max_zeff_correction = max_zeff_correction
      thermo_pot%min_zeff_correction = min_zeff_correction
      thermo_pot%use_integral_control = use_integral_control
      thermo_pot%integral_gain = integral_gain
      thermo_pot%current_zeff_corr = 0.0_dp
      thermo_pot%integral_error = 0.0_dp
      thermo_pot%active = .TRUE.

   END SUBROUTINE thermopotentiostat_create

! **************************************************************************************************
!> \brief Release a thermopotentiostat structure
!> \param thermo_pot the thermopotentiostat to release
!> \note The caller is responsible for deallocating the structure after this call
! **************************************************************************************************
   SUBROUTINE thermopotentiostat_release(thermo_pot)

      TYPE(thermopotentiostat_type), POINTER             :: thermo_pot

      ! No internal pointers to clean up
      ! Deallocation is handled by the caller (consistent with other release routines)
      ! Do not NULLIFY here - that would prevent the caller from deallocating
      MARK_USED(thermo_pot)

   END SUBROUTINE thermopotentiostat_release

END MODULE thermopotentiostat_types
