# Thermopotentiostat Implementation Summary

## Overview

This document describes the implementation of the **Thermopotentiostat ensemble** in CP2K. The thermopotentiostat varies the `zeff_correction` of an atom type based on the system dipole moment, allowing control of electrochemical potential in surface simulations.

The implementation follows the approach from:
> SS, MT, MWF, JN PRL, 2018, 120, 246801

where electron density is adjusted according to the applied core correction.

## Files Created

### 1. `src/motion/thermopotentiostat_types.F`
Type definitions for the thermopotentiostat structure including:
- Target dipole moment
- Coupling strength (feedback gain)
- Current zeff correction value
- Min/max limits for zeff correction
- PI (proportional-integral) control parameters
- `thermopotentiostat_create()` and `thermopotentiostat_release()` subroutines

### 2. `src/motion/thermopotentiostat_methods.F`
Core logic for the thermopotentiostat:
- `thermopotentiostat_init()` - Initialize from MD input section
- `thermopotentiostat_update()` - Update zeff correction based on current dipole moment
- `update_zeff_correction_for_dipole()` - Apply changes to GTH potential and update `qs_env%total_zeff_corr` for electron density adjustment

## Files Modified

### 3. `src/input_constants.F`
Added new ensemble constant:
```fortran
thermopotentiostat_ensemble = 15
```

### 4. `src/motion/input_cp2k_md.F`
- Added `THERMOPOTENTIOSTAT` to the ensemble enum values
- Created `create_thermopotentiostat_section()` subroutine with keywords:
  - `TARGET_DIPOLE` - Target dipole moment (in Debye)
  - `DIPOLE_DIRECTION` - X, Y, or Z component to control
  - `COUPLING_STRENGTH` - Proportional feedback gain
  - `ELEMENT` - Element symbol of atom kind to modify
  - `MAX_ZEFF_CORRECTION` - Maximum allowed correction
  - `MIN_ZEFF_CORRECTION` - Minimum allowed correction
  - `USE_INTEGRAL_CONTROL` - Enable PI control
  - `INTEGRAL_GAIN` - Integral term gain

### 5. `src/motion/integrator.F`
Added `thermopotentiostat()` subroutine:
- Based on NVE integrator structure
- Adds feedback control of zeff_correction after force calculation
- Updates GTH potential parameters at each timestep

### 6. `src/motion/velocity_verlet_control.F`
Registered the new ensemble in the velocity Verlet dispatch:
```fortran
CASE (thermopotentiostat_ensemble)
   CALL thermopotentiostat(md_env, globenv)
```

### 7. `src/motion/md_environment_types.F`
Added thermopotentiostat pointer to `md_environment_type`:
- Added `thermopotentiostat` pointer member
- Extended `get_md_env()` and `set_md_env()` interfaces
- Added cleanup in `md_env_release()`

### 8. `src/motion/md_energies.F`
Added output support:
- Imports `thermopotentiostat_type` and `print_thermopotentiostat_status`
- Calls `print_thermopotentiostat_status()` in `md_write_output`

## Output Files

### Main Output (`.mdLog` or stdout)
At each timestep, the main output includes:
```
THERMOPOTENTIOSTAT| Status:
THERMOPOTENTIOSTAT| Current dipole (target dir.)  [a.u.]:          0.123456789
THERMOPOTENTIOSTAT| Target dipole                 [a.u.]:          0.000000000
THERMOPOTENTIOSTAT| Dipole error                  [a.u.]:          0.123456789
THERMOPOTENTIOSTAT| Current zeff correction       [a.u.]:         -0.012345678
```

### Dedicated File (`.thermopotentiostat`)
A dedicated output file with columns:
| Column | Description |
|--------|-------------|
| Step Nr. | MD timestep number |
| Time[fs] | Simulation time in femtoseconds |
| Dipole_X[a.u.] | X-component of dipole moment |
| Dipole_Y[a.u.] | Y-component of dipole moment |
| Dipole_Z[a.u.] | Z-component of dipole moment |
| Zeff_Corr[a.u.] | Current zeff correction value |

Control output with print key:
```
&MOTION
  &MD
    &PRINT
      &THERMOPOTENTIOSTAT
        &EACH
          MD 1
        &END EACH
      &END THERMOPOTENTIOSTAT
    &END PRINT
  &END MD
&END MOTION
```

## Usage Example

```
&FORCE_EVAL
  &DFT
    &PRINT
      &MOMENTS ON
      &END MOMENTS
    &END PRINT
  &END DFT
&END FORCE_EVAL

&MOTION
  &MD
    ENSEMBLE THERMOPOTENTIOSTAT
    TIMESTEP 0.5
    STEPS 1000
    TEMPERATURE 300
    &THERMOPOTENTIOSTAT
      TARGET_DIPOLE 0.0        ! Target dipole in Debye
      DIPOLE_DIRECTION Z       ! Control Z-component
      COUPLING_STRENGTH 0.1    ! Feedback gain
      ELEMENT Li               ! Modify Li atoms
      MAX_ZEFF_CORRECTION 1.0
      MIN_ZEFF_CORRECTION -1.0
      USE_INTEGRAL_CONTROL .FALSE.
      INTEGRAL_GAIN 0.01
    &END THERMOPOTENTIOSTAT
  &END MD
&END MOTION
```

## Key Integration Points

1. **Electron Density Adjustment**: The `update_zeff_correction_for_dipole()` routine modifies:
   - GTH potential's `zeff`
   - `ccore_charge` coefficient
   - `cerf_ppl` coefficient
   - `qs_env%total_zeff_corr`

2. **Automatic MO Occupation**: The existing code in `src/qs_mo_occupation.F` (lines ~465-520) automatically adjusts electron occupation based on `total_zeff_corr`.

3. **Dipole Calculation**: Requires `DFT%PRINT%MOMENTS` to be enabled so the dipole moment is computed and stored in results.

## CMake Build System Integration

The new source files must be added to `src/CMakeLists.txt`. The motion sources are listed 
in the `CP2K_SRCS_F` list, with motion files starting around line 1100.

### Location in CMakeLists.txt

Add the two new files in alphabetical order around **line 1178** (after `thermal_region_utils.F`):

```cmake
  motion/thermal_region_types.F
  motion/thermal_region_utils.F
  motion/thermopotentiostat_methods.F    # <-- ADD THIS
  motion/thermopotentiostat_types.F      # <-- ADD THIS
  motion/thermostat/al_system_dynamics.F
```

### Build Commands

After adding the files to CMakeLists.txt:

```bash
# Configure (from build directory)
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build -j $(nproc)
```

---

## Restart File Support

To enable continuation of thermopotentiostat simulations across restarts, the following changes
are needed to preserve the `current_zeff_corr` state.

### 1. Add Restart Section to Input Definition

In `src/motion/input_cp2k_md.F`, add a restart subsection to the THERMOPOTENTIOSTAT section:

```fortran
! Inside create_thermopotentiostat_section(), after existing keywords:

CALL section_create(subsection, __LOCATION__, name="RESTART", &
     description="Restart information for thermopotentiostat.", &
     n_keywords=1, n_subsections=0, repeats=.FALSE.)

CALL keyword_create(keyword, __LOCATION__, name="CURRENT_ZEFF_CORR", &
     description="Current zeff correction value from previous run.", &
     usage="CURRENT_ZEFF_CORR 0.05", &
     default_r_val=0.0_dp, unit_str="au_e")
CALL section_add_keyword(subsection, keyword)
CALL keyword_release(keyword)

CALL section_add_subsection(section, subsection)
CALL section_release(subsection)
```

### 2. Extend md_environment_types.F

The `md_environment_type` needs a pointer to the thermopotentiostat. In 
`src/motion/md_environment_types.F`:

```fortran
! Add to USE statements:
USE thermopotentiostat_types, ONLY: thermopotentiostat_type

! Add to md_environment_type:
TYPE(thermopotentiostat_type), POINTER :: thermopotentiostat => NULL()

! Add to get_md_env interface:
TYPE(thermopotentiostat_type), OPTIONAL, POINTER, INTENT(OUT) :: thermopotentiostat

! Add to set_md_env interface:
TYPE(thermopotentiostat_type), OPTIONAL, POINTER, INTENT(IN) :: thermopotentiostat
```

### 3. Add Write Logic to input_cp2k_restarts.F

In `src/motion/input_cp2k_restarts.F`, inside the `write_restart_md` subroutine (around line 500),
add code to save the current zeff correction:

```fortran
! Add USE statement at top of module:
USE thermopotentiostat_types, ONLY: thermopotentiostat_type

! In write_restart_md, declare variable:
TYPE(thermopotentiostat_type), POINTER :: thermopotentiostat

! After getting md_env data:
CALL get_md_env(md_env=md_env, thermopotentiostat=thermopotentiostat)

! Write the value if thermopotentiostat is active:
IF (ASSOCIATED(thermopotentiostat)) THEN
   work_section => section_vals_get_subs_vals(motion_section, &
                   "MD%THERMOPOTENTIOSTAT%RESTART")
   CALL section_vals_val_set(work_section, "CURRENT_ZEFF_CORR", &
                             r_val=thermopotentiostat%current_zeff_corr)
END IF
```

### 4. Add Read Logic in thermopotentiostat_init

In `src/motion/thermopotentiostat_methods.F`, modify `thermopotentiostat_init` to check for
restart values:

```fortran
! After reading normal input, check for restart value:
restart_section => section_vals_get_subs_vals(thermopotentiostat_section, "RESTART")
CALL section_vals_val_get(restart_section, "CURRENT_ZEFF_CORR", r_val=restart_zeff)

IF (ABS(restart_zeff) > 1.0E-12_dp) THEN
   thermopotentiostat%current_zeff_corr = restart_zeff
   ! Also initialize integral term for PI controller stability
   IF (thermopotentiostat%use_integral_control .AND. &
       ABS(thermopotentiostat%integral_gain) > 1.0E-12_dp) THEN
      thermopotentiostat%integral_error = restart_zeff / thermopotentiostat%integral_gain
   END IF
END IF
```

---

## Testing

Regression tests should be added in `tests/QS/regtest-*/`. A minimal test input would verify:

1. The ensemble runs without crashing
2. Dipole moment changes in response to zeff correction
3. Restart correctly preserves thermopotentiostat state

---

## Future Enhancements

1. Binary restart support for large-scale simulations
2. Detailed output logging with history of zeff corrections
3. Multi-region support for different electrode surfaces
4. Multiple element types with independent corrections
5. Adaptive timestep coupling strength

The thermopotentiostat uses a feedback control loop:

```
dipole_error = current_dipole(direction) - target_dipole

# Proportional term
delta_zeff = coupling_strength * dipole_error * dt

# Optional integral term (PI control)
if use_integral_control:
    integral_error += dipole_error * dt
    delta_zeff += integral_gain * integral_error

# Apply with limits
new_zeff_corr = clip(current_zeff_corr + delta_zeff, min, max)
```

This drives the system dipole towards the target value by adjusting the effective nuclear charge of the specified element, which in turn modifies the electron density distribution.
