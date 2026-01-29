# LaMEM Development Notes

## Project Overview

LaMEM (Lithosphere and Mantle Evolution Model) is a parallel 3D thermo-mechanical geodynamics code.

- **Language**: C++ with Fortran bindings
- **Build System**: Make-based (`src/Makefile`)
- **Parallelization**: PETSc-based distributed computing
- **Current Branch**: `yl/fastscape-coupling` - adds coupling with FastScape landscape evolution model

## Build Commands

```bash
make mode=opt all           # Standard optimized build
make mode=optFS all         # With FastScape coupling
make mode=deb all           # Debug build
```

## Key Directories

| Directory | Purpose |
|-----------|---------|
| `src/` | Core C++ implementation |
| `test/` | 38+ test cases |
| `input_models/` | Pre-built model examples |
| `doc/` | Documentation |

## FastScape Coupling Module

New files added on this branch:
- `src/fastscape.cpp` (3,094 lines)
- `src/fastscape.h` (382 lines)

---

## Fixed Issues: 64-bit Compilation Problems in fastscape.cpp

The following issues have been fixed (2026-01-28):

### 1. Fortran Interface Type Mismatch (FIXED)

The Fortran functions in `fastscape.h` are declared with `int` (32-bit), but were being called with `PetscInt` (64-bit on 64-bit builds).

**Solution applied**: Added local `int` variables (`nx_f`, `ny_f`, `ibc_f`, `istep_f`) and copy values before/after Fortran calls:
```cpp
int nx_f = (int)FSLib->nx_solve;
int ny_f = (int)FSLib->ny_solve;
fastscape_set_nx_ny_(&nx_f, &ny_f);
```

### 2. Printf Format Specifiers (FIXED)

All `%d` format specifiers for `PetscInt` values have been replaced with `%" PetscInt_FMT "`:
```cpp
PetscPrintf(PETSC_COMM_WORLD, "[nodeX, nodeY]: [%" PetscInt_FMT ", %" PetscInt_FMT "]\n", nx, ny);
```

### 3. Lambda Parameter Types (FIXED)

Changed `int` to `PetscInt` in lambda parameters:
- `update_coords` lambda: `int count` -> `PetscInt count`
- `manage_output` lambda: `int flag` -> `PetscInt flag`

### 4. Loop Variable Types (FIXED)

Changed loop variables comparing with `size_t` to use `size_t`:
```cpp
for (size_t ii = 0; ii < output_arrays.size(); ii++)
```

### 5. Deprecated PETSC_NULL (FIXED)

Replaced all `PETSC_NULL` with `PETSC_NULLPTR` to silence PETSc 3.19+ deprecation warnings.

---

## Remaining Minor Warnings

The following warnings remain but don't affect functionality:
- Float-to-int conversion warnings in grid calculations (lines 503, 884, 909, 1048)
- Unused variable warnings in `InterpolationFor2DNonUniformGrid` function

These are code style issues that can be addressed in future cleanup.
