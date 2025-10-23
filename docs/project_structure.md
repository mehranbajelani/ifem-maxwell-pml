# iFEM-Maxwell-PML Project Structure

## Overview
This project extends the iFEM Maxwell solver to implement realistic electromagnetic scattering simulations with:
- Spherical air domain
- Perfectly Matched Layer (PML) absorption
- Silver-Müller Scattering Boundary Condition (SBC)
- Scattering cross-section post-processing

## Directory Structure

```
ifem-maxwell-pml/
├── external/ifem/              # Original iFEM (vendored)
├── src/                        # Main source code
│   ├── geometry/              # Meshing and geometry utilities
│   ├── physics/               # SBC, PML, material functions
│   ├── post/                  # Post-processing utilities
│   └── examples/              # Driver scripts and examples
├── tests/                     # Unit tests and validation
├── docs/                      # Documentation
├── LICENSE                    # Project license
├── README.md                  # Project overview
└── PML_plan_10_21_2025.pdf   # Detailed implementation plan
```

## Implementation Phases

### Phase 1: Geometry and Domain Definition ✅
- **Goal**: Build three-region geometry (SiC cube ⊂ air sphere ⊂ PML shell)
- **Files**: `src/geometry/makeSphericalAirAndPML.m`
- **Output**: Tagged mesh with region IDs, SBC boundary faces, PML grading
- **Status**: Implemented with hybrid region tagging and PML mesh grading

### Phase 2: Scattering Boundary Condition (SBC)
- **Goal**: Implement Silver-Müller boundary condition
- **Files**: `src/physics/addSBC.m`, `src/physics/incidentPlaneWave.m`
- **Output**: SBC matrix and vector assembly

### Phase 3: Perfectly Matched Layer (PML)
- **Goal**: Add spherical PML absorption
- **Files**: `src/physics/pmlSphericalTensors.m`
- **Output**: Anisotropic tensor assembly

### Phase 4: Post-Processing
- **Goal**: Compute scattering cross-section
- **Files**: `src/post/sigma_s_post.m`
- **Output**: Normalized scattering cross-section

### Phase 5: Convergence Testing
- **Goal**: Verify numerical stability
- **Files**: `tests/convergence_tests.m`
- **Output**: Convergence validation

### Phase 6: Real Material Data
- **Goal**: Use real SiC permittivity data
- **Files**: `src/physics/eps_SiC_6H_400C.m`
- **Output**: Wavelength-dependent results

## Key Parameters

- **SiC Cube**: L = 1 μm (side length)
- **Air Sphere**: r_SBC = 6 μm
- **PML Shell**: t_PML = 2 μm, r_PML = 8 μm
- **Wavelength**: λ = 5 μm
- **Frequency**: f = c/λ ≈ 60 THz
- **Incident Field**: E_0 = 1 V/m, x-polarized, -z propagation

## Validation Strategy

1. **Phase 1**: Visual inspection of geometry and mesh quality
2. **Phase 2**: Plane wave injection test (no scatterer)
3. **Phase 3**: PML reflection < 1%
4. **Phase 4**: Scattering cross-section convergence
5. **Phase 5**: Mesh and PML parameter convergence
6. **Phase 6**: Comparison with COMSOL and Pouria et al. (2025)
