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
   USE parallel_rng_types,              ONLY: rng_stream_type

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
      ! Target potential in atomic units (converted from Volts)
      REAL(dp)                              :: target_potential = 0.0_dp
      ! Geometric capacitance in atomic units (converted from attofarads)
      REAL(dp)                              :: geometric_capacitance = 0.1_dp
      ! Time constant for thermopotentiostat response in atomic units (converted from fs)
      REAL(dp)                              :: timecon = 0.0_dp
      ! Current zeff correction value being applied
      REAL(dp)                              :: current_zeff_corr = 0.0_dp
      ! Direction for dipole control: 1=x, 2=y, 3=z
      INTEGER                               :: direction = 3
      ! Index of the atom kind to modify
      INTEGER                               :: atom_kind_index = -1
      ! Element symbol of the atom kind to modify
      CHARACTER(LEN=default_string_length)  :: element_symbol = ""
      ! Store the initial zeff_correction from input
      REAL(dp)                              :: initial_zeff_corr = 0.0_dp
      ! Last measured dipole moment for output
      REAL(dp), DIMENSION(3)                :: last_dipole = [0.0_dp, 0.0_dp, 0.0_dp]
      ! Last calculated current potential for output
      REAL(dp)                              :: last_current_potential = 0.0_dp
      ! Last dissipative term (change in electron count) for output
      REAL(dp)                              :: last_dissipative_term = 0.0_dp
      ! Last fluctuation term (change in electron count) for output
      REAL(dp)                              :: last_fluctuation_term = 0.0_dp
      ! Flag indicating if thermopotentiostat is active
      LOGICAL                               :: active = .FALSE.
      ! Random number generator for stochastic dynamics
      TYPE(rng_stream_type)                 :: rng_stream = rng_stream_type()
   END TYPE thermopotentiostat_type

CONTAINS

! **************************************************************************************************
!> \brief Create and initialize a thermopotentiostat structure
!> \param thermo_pot the thermopotentiostat to create
!> \param target_potential target potential in atomic units
!> \param geometric_capacitance geometric capacitance in atomic units
!> \param timecon time constant in atomic units
!> \param direction dipole direction (1=x, 2=y, 3=z)
!> \param element_symbol element symbol of atom kind to modify
! **************************************************************************************************
   SUBROUTINE thermopotentiostat_create(thermo_pot, target_potential, geometric_capacitance, &
                                        timecon, direction, element_symbol)

      TYPE(thermopotentiostat_type), POINTER             :: thermo_pot
      REAL(dp), INTENT(IN)                               :: target_potential
      REAL(dp), INTENT(IN)                               :: geometric_capacitance
      REAL(dp), INTENT(IN)                               :: timecon
      INTEGER, INTENT(IN)                                :: direction
      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol

      CPASSERT(.NOT. ASSOCIATED(thermo_pot))
      ALLOCATE (thermo_pot)

      thermo_pot%target_potential = target_potential
      thermo_pot%geometric_capacitance = geometric_capacitance
      thermo_pot%timecon = timecon
      thermo_pot%direction = direction
      thermo_pot%element_symbol = TRIM(element_symbol)
      thermo_pot%current_zeff_corr = 0.0_dp
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
